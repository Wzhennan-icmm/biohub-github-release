use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::error::Error;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

const VERSION: &str = "0.2.1";

fn print_help() {
    println!(
        "biohub {version}

Usage:
  biohub <command> [subcommand] [options]

Commands:
  rename hjjn-genes      -i <input> [-o <output>]
  rename scaffolds       -i <input> -l <map> [-o <output>]
  rename fasta-scaffolds -i <input> -l <map> [-o <output>]
  blast reciprocal       -i <blast> -r <reverse> [-o <output>]
  gff filter-ncbi       --input/--gff <gff> [-o <output>]
  gff filter-gemoma     --input/--gff <gff> [-o <output>]
  gff convert-ty1-hjjn  --gff <gff> --bed <bed> [-o <output>]
  fasta longest-transcript -f <fasta> [-o <output>]
  stats coverage-ratio   -i <input> -r <reference> [-o <output>]
  stats hic-matrix-reindex -b <bed> -m <matrix> -p <group> [-o <output>]
  stats wgcna-weight     -i <weight-file> [-o <output>]
  scripts catalog
  scripts run <script-id> [options]
  psmc merge             -d <dir> [-p <pattern>] [-o <output>]

  --help         Show this message
  -h             Show this message
  --version      Show version
",
        version = VERSION
    );
}

fn needs_help(args: &[String]) -> bool {
    args.iter()
        .any(|a| a == "--help" || a == "-h" || a == "help")
}

fn expand_path(path: &str) -> String {
    if path == "~" {
        return env::var("HOME").unwrap_or_else(|_| path.to_string());
    }
    if let Some(rest) = path.strip_prefix("~/") {
        if let Ok(home) = env::var("HOME") {
            return Path::new(&home).join(rest).to_string_lossy().to_string();
        }
    }
    path.to_string()
}

fn open_reader(path: &str) -> Result<BufReader<File>> {
    let p = expand_path(path);
    let fh = File::open(&p)?;
    Ok(BufReader::new(fh))
}

fn open_writer(path: Option<&str>) -> Result<Box<dyn Write>> {
    match path {
        None => Ok(Box::new(BufWriter::new(io::stdout()))),
        Some("-") => Ok(Box::new(BufWriter::new(io::stdout()))),
        Some(out) => {
            let p = expand_path(out);
            let fh = File::create(p)?;
            Ok(Box::new(BufWriter::new(fh)))
        }
    }
}

fn parse_required(args: &[String], idx: &mut usize, flag: &str) -> Result<String> {
    if *idx + 1 >= args.len() {
        return Err(format!("missing value for {flag}").into());
    }
    *idx += 1;
    Ok(args[*idx].clone())
}

fn parse_two_col_mapping(path: &str, reverse: bool) -> Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    let mut reader = open_reader(path)?;
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let row = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if row.is_empty() || row.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = row.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        if reverse {
            map.insert(cols[1].to_string(), cols[0].to_string());
        } else {
            map.insert(cols[0].to_string(), cols[1].to_string());
        }
    }
    Ok(map)
}

#[derive(Clone)]
struct ScriptSpec {
    id: String,
    source: String,
    description: String,
    status: String,
    note: String,
}

#[derive(Clone)]
struct ScatterPoint {
    x: f64,
    y: f64,
    group: String,
}

const SVGPLOT_COLORS: [&str; 12] = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#17becf",
    "#bcbd22",
    "#e6550d",
    "#6b6bd6",
];

fn svg_escape_xml(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    for c in input.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '\'' => out.push_str("&apos;"),
            '"' => out.push_str("&quot;"),
            _ => out.push(c),
        }
    }
    out
}

fn svg_color(idx: usize) -> &'static str {
    SVGPLOT_COLORS[idx % SVGPLOT_COLORS.len()]
}

fn downsample<T: Clone>(data: &[T], max_points: usize) -> Vec<T> {
    if data.len() <= max_points {
        return data.to_vec();
    }
    let step = (data.len() as f64 / max_points as f64).ceil() as usize;
    data.iter().step_by(step.max(1)).cloned().collect()
}

fn write_scatter_svg(
    path: &Path,
    title: &str,
    x_label: &str,
    y_label: &str,
    points: &[ScatterPoint],
    point_radius: f64,
    width: usize,
    height: usize,
    custom_x_ticks: Option<&[(f64, String)]>,
) -> Result<i32> {
    if points.is_empty() {
        return Err("no data points for plot".into());
    }

    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;

    let mut groups: Vec<String> = Vec::new();
    let mut group_to_idx: HashMap<String, usize> = HashMap::new();
    for p in points {
        x_min = x_min.min(p.x);
        x_max = x_max.max(p.x);
        y_min = y_min.min(p.y);
        y_max = y_max.max(p.y);
        if !group_to_idx.contains_key(&p.group) {
            let idx = groups.len();
            groups.push(p.group.clone());
            group_to_idx.insert(p.group.clone(), idx);
        }
    }

    if (x_max - x_min).abs() < 1e-9 {
        x_max += 1.0;
    }
    if (y_max - y_min).abs() < 1e-9 {
        y_max += 1.0;
    }

    let margin_left = 80.0;
    let margin_right = 24.0;
    let margin_top = 52.0;
    let margin_bottom = 64.0;

    let w = width as f64;
    let h = height as f64;
    let plot_w = w - margin_left - margin_right;
    let plot_h = h - margin_top - margin_bottom;

    let map_x = |x: f64| margin_left + ((x - x_min) / (x_max - x_min)) * plot_w;
    let map_y = |y: f64| h - margin_bottom - ((y - y_min) / (y_max - y_min)) * plot_h;

    let mut out = BufWriter::new(File::create(path)?);
    writeln!(
        out,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">",
        w = w,
        h = h
    )?;
    writeln!(out, "<rect x=\"0\" y=\"0\" width=\"{w}\" height=\"{h}\" fill=\"#ffffff\"/>")?;
    writeln!(
        out,
        "<text x=\"{title_x}\" y=\"{title_y}\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"22\" text-anchor=\"middle\">{}</text>",
        svg_escape_xml(title),
        title_x = margin_left + plot_w / 2.0,
        title_y = 28.0
    )?;

    writeln!(
        out,
        "<rect x=\"{l}\" y=\"{t}\" width=\"{pw}\" height=\"{ph}\" fill=\"none\" stroke=\"#333\" stroke-width=\"1\"/>",
        l = margin_left,
        t = margin_top,
        pw = plot_w,
        ph = plot_h
    )?;

    let x_ticks = 8usize;
    for i in 0..=x_ticks {
        let ratio = i as f64 / x_ticks as f64;
        let x = margin_left + ratio * plot_w;
        let value = x_min + ratio * (x_max - x_min);
        writeln!(
            out,
            "<line x1=\"{x}\" y1=\"{t}\" x2=\"{x}\" y2=\"{y2}\" stroke=\"#ddd\" stroke-width=\"1\"/>",
            t = margin_top + plot_h,
            y2 = margin_top
        )?;
        writeln!(
            out,
            "<text x=\"{x}\" y=\"{yt}\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"11\" text-anchor=\"middle\" fill=\"#333\">{value:.3}</text>",
            x = x,
            yt = margin_top + plot_h + 22.0,
            value = value
        )?;
    }

    if let Some(ticks) = custom_x_ticks {
        for (value, label) in ticks.iter() {
            if *value < x_min || *value > x_max {
                continue;
            }
            let x = map_x(*value);
            writeln!(
                out,
                "<line x1=\"{x}\" y1=\"{t}\" x2=\"{x}\" y2=\"{b}\" stroke=\"#555\" stroke-width=\"1\" stroke-dasharray=\"3 3\"/>",
                t = margin_top + plot_h,
                b = margin_top
            )?;
            writeln!(
                out,
                "<text x=\"{x}\" y=\"{yt}\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"10\" text-anchor=\"middle\" fill=\"#555\" transform=\"rotate(45 {x} {yt})\">{}</text>",
                svg_escape_xml(label),
                x = x,
                yt = margin_top + plot_h + 36.0
            )?;
        }
    }

    let y_ticks = 8usize;
    for i in 0..=y_ticks {
        let ratio = i as f64 / y_ticks as f64;
        let y = margin_top + plot_h - ratio * plot_h;
        let value = y_min + ratio * (y_max - y_min);
        writeln!(
            out,
            "<line x1=\"{l}\" y1=\"{y}\" x2=\"{r}\" y2=\"{y}\" stroke=\"#ddd\" stroke-width=\"1\"/>",
            l = margin_left,
            r = margin_left + plot_w
        )?;
        writeln!(
            out,
            "<text x=\"{x}\" y=\"{y}\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"11\" text-anchor=\"end\" dy=\"4\" fill=\"#333\">{value:.3}</text>",
            x = margin_left - 8.0,
            y = y,
            value = value
        )?;
    }

    let x_axis_y = margin_top + plot_h;
    let y_axis_x = margin_left;
    writeln!(
        out,
        "<line x1=\"{y_axis_x}\" y1=\"{x_axis_y}\" x2=\"{x_axis_x2}\" y2=\"{x_axis_y}\" stroke=\"#333\" stroke-width=\"1.2\"/>",
        x_axis_x2 = margin_left + plot_w
    )?;
    writeln!(
        out,
        "<line x1=\"{y_axis_x}\" y1=\"{y_axis_y}\" x2=\"{y_axis_x}\" y2=\"{y_axis_end}\" stroke=\"#333\" stroke-width=\"1.2\"/>",
        y_axis_y = margin_top,
        y_axis_end = margin_top + plot_h
    )?;

    writeln!(
        out,
        "<text x=\"{}\" y=\"{}\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"14\" text-anchor=\"middle\">{}</text>",
        margin_left + plot_w / 2.0,
        h - 20.0,
        svg_escape_xml(x_label)
    )?;
    writeln!(
        out,
        "<text x=\"{x}\" y=\"{y}\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"14\" text-anchor=\"middle\" transform=\"rotate(-90 {x} {y})\">{label}</text>",
        x = margin_left - 46.0,
        y = margin_top + plot_h / 2.0,
        label = svg_escape_xml(y_label)
    )?;

    for p in points {
        let cx = map_x(p.x);
        let cy = map_y(p.y);
        let gi = *group_to_idx.get(&p.group).unwrap_or(&0);
        let color = svg_color(gi);
        writeln!(
            out,
            "<circle cx=\"{:.3}\" cy=\"{:.3}\" r=\"{r}\" fill=\"{color}\" fill-opacity=\"0.45\" stroke=\"#333\" stroke-width=\"0.2\"/>",
            cx,
            cy,
            r = point_radius,
            color = color
        )?;
    }

    let legend_x = margin_left + plot_w - 140.0;
    let mut legend_y = margin_top + 12.0;
    for (idx, group) in groups.iter().enumerate() {
        let color = svg_color(idx);
        writeln!(
            out,
            "<circle cx=\"{}\" cy=\"{}\" r=\"5\" fill=\"{}\" />",
            legend_x,
            legend_y,
            color
        )?;
        writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"11\" fill=\"#333\">{}</text>",
            legend_x + 8.0,
            legend_y + 4.0,
            svg_escape_xml(group)
        )?;
        legend_y += 15.0;
    }

    writeln!(out, "</svg>")?;
    Ok(0)
}
const SCRIPT_CATALOG_TEXT: &str = include_str!("script_catalog.tsv");

fn load_script_catalog() -> Vec<ScriptSpec> {
    let mut out = Vec::new();
    for line in SCRIPT_CATALOG_TEXT.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut cells = line.split('\t');
        let (Some(id), Some(source), Some(description), Some(status), Some(note)) = (
            cells.next(),
            cells.next(),
            cells.next(),
            cells.next(),
            cells.next(),
        ) else {
            continue;
        };

        out.push(ScriptSpec {
            id: id.to_string(),
            source: source.to_string(),
            description: description.to_string(),
            status: status.to_string(),
            note: note.to_string(),
        });
    }
    out
}

fn print_script_catalog() {
    println!(
        "{:<32} {:<45} {:<12} DESCRIPTION",
        "ID", "SOURCE", "STATUS"
    );
    println!("{:-<105}", "");
    for spec in load_script_catalog() {
        println!(
            "{:<32} {:<45} {:<12} {}",
            spec.id, spec.source, spec.status, spec.description
        );
        if !spec.note.is_empty() {
            println!("{:<44} NOTE: {}", "", spec.note);
        }
    }
}

fn get_opt(args: &[String], keys: &[&str]) -> Option<String> {
    let mut i = 0usize;
    while i + 1 < args.len() {
        let arg = args[i].as_str();
        if keys.iter().any(|k| arg == *k) {
            return Some(args[i + 1].clone());
        }
        i += 1;
    }
    None
}

fn get_required_opt(args: &[String], keys: &[&str], label: &str) -> Result<String> {
    get_opt(args, keys).ok_or_else(|| format!("missing required option for {label}").into())
}

fn parse_usize_arg(args: &[String], keys: &[&str], label: &str) -> Result<usize> {
    let raw = get_required_opt(args, keys, label)?;
    raw.parse::<usize>()
        .map_err(|_| format!("invalid {label}: {raw}").into())
}

fn run_hjjn_genes(input: &str, output: Option<&str>) -> Result<i32> {
    fn normalize_gene(name: &str) -> String {
        let lower = name.to_lowercase();
        let Some(pos) = lower.rfind("gene") else {
            return name.to_string();
        };
        let num = &name[pos + 4..];
        if let Some(dot) = num.find('.') {
            let num_part = &num[..dot];
            let ext = &num[dot..];
            if num_part.chars().all(|c| c.is_ascii_digit()) {
                return format!("{}gene{:0>5}{}", &name[..pos], num_part, ext);
            }
        } else if num.chars().all(|c| c.is_ascii_digit()) {
            return format!("{}gene{:0>5}", &name[..pos], num.parse::<usize>().unwrap_or(0));
        }
        name.to_string()
    }

    let mut reader = open_reader(input)?;
    let mut out = open_writer(output)?;
    let mut row = String::new();

    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        writeln!(out, "{}\t{}", cols[0], normalize_gene(cols[1]))?;
    }
    Ok(0)
}

fn run_scaffold_rename(input: &str, map_file: &str, output: Option<&str>) -> Result<i32> {
    let mut map: HashMap<String, String> = HashMap::new();
    let mut mreader = open_reader(map_file)?;
    let mut row = String::new();
    while mreader.read_line(&mut row)? > 0 {
        let t = row.trim_end_matches(['\n', '\r']).trim().to_string();
        row.clear();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = t.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        map.insert(fields[1].to_string(), fields[0].to_string());
    }

    let mut reader = open_reader(input)?;
    let mut out = open_writer(output)?;
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let fields: Vec<&str> = raw.split('\t').collect();
        if fields.is_empty() {
            continue;
        }
        let replaced = map
            .get(fields[0])
            .cloned()
            .unwrap_or_else(|| fields[0].to_string());
        let mut cols: Vec<String> = fields.iter().map(|f| f.to_string()).collect();
        cols[0] = replaced;
        writeln!(out, "{}", cols.join("\t"))?;
    }
    Ok(0)
}

fn run_fasta_scaffold_rename(input: &str, map_file: &str, output: Option<&str>) -> Result<i32> {
    let map = parse_two_col_mapping(map_file, false)?;
    let mut reader = open_reader(input)?;
    let mut out = open_writer(output)?;
    let mut line = String::new();

    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.starts_with('>') {
            let h = raw.trim_start_matches('>');
            let id = h.split_whitespace().next().unwrap_or("");
            let mapped = map.get(id).cloned().unwrap_or_else(|| id.to_string());
            writeln!(out, ">{mapped}")?;
        } else if !raw.is_empty() {
            writeln!(out, "{}", raw)?;
        }
    }
    Ok(0)
}

fn run_reciprocal(blast: &str, reverse: &str, output: Option<&str>) -> Result<i32> {
    fn read_pairs(path: &str) -> Result<HashMap<String, String>> {
        let mut map = HashMap::new();
        let mut reader = open_reader(path)?;
        let mut row = String::new();
        while reader.read_line(&mut row)? > 0 {
            let raw = row.trim_end_matches(['\n', '\r']).to_string();
            row.clear();
            if raw.is_empty() {
                continue;
            }
            let fields: Vec<&str> = raw.split_whitespace().collect();
            if fields.len() < 2 {
                continue;
            }
            map.entry(fields[0].to_string()).or_insert_with(|| fields[1].to_string());
        }
        Ok(map)
    }

    let a = read_pairs(blast)?;
    let b = read_pairs(reverse)?;
    let mut out = open_writer(output)?;

    for (q, t) in a {
        if b.get(&t).is_some_and(|v| v == &q) {
            writeln!(out, "{q}\t{t}")?;
        }
    }
    Ok(0)
}

fn flush_ncbi_gene(
    gene_line: &Option<String>,
    records: &[String],
    gene_num: u32,
    out: &mut dyn Write,
) -> Result<()> {
    let Some(gene_line) = gene_line else {
        return Ok(());
    };
    if records.is_empty() {
        return Ok(());
    }
    let gcols: Vec<&str> = gene_line.split('\t').collect();
    if gcols.len() < 9 {
        return Ok(());
    }

    let gene_id = format!("TA-gene{gene_num:05}");
    let gene_name = format!("TA{gene_num:05}");
    writeln!(out, "{}\tID={gene_id};Name={gene_name}", gcols[..8].join("\t"))?;

    let mut mrna_num = 1u32;
    let mut exon_num = 1u32;
    let mut current_mrna = String::new();

    for raw in records {
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        match cols[2] {
            "mRNA" => {
                current_mrna = format!("{gene_id}.{mrna_num}");
                let name = format!("{gene_name}.{mrna_num}");
                mrna_num += 1;
                writeln!(
                    out,
                    "{}\tID={};Parent={};Name={}",
                    cols[..8].join("\t"),
                    current_mrna,
                    gene_name,
                    name
                )?;
            }
            "CDS" => {
                if current_mrna.is_empty() {
                    continue;
                }
                let mut exon_cols = cols.clone();
                exon_cols[2] = "exon";
                writeln!(
                    out,
                    "{}\tID={}.exon{exon_num};Parent={}",
                    exon_cols[..8].join("\t"),
                    current_mrna,
                    current_mrna
                )?;
                exon_num += 1;
                writeln!(
                    out,
                    "{}\tID=cds.{};Parent={}",
                    cols[..8].join("\t"),
                    current_mrna,
                    current_mrna
                )?;
            }
            _ => {}
        }
    }
    Ok(())
}

fn run_filter_ncbi(input: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(input)?;
    let mut out = open_writer(output)?;
    let mut row = String::new();
    let mut current_gene: Option<String> = None;
    let mut records: Vec<String> = Vec::new();
    let mut gene_num: u32 = 1;

    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        if raw.starts_with('#') {
            writeln!(out, "{}", raw)?;
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        if cols[2] == "gene" {
            flush_ncbi_gene(&current_gene, &records, gene_num, &mut *out)?;
            gene_num += 1;
            current_gene = Some(raw);
            records.clear();
        } else if current_gene.is_some() && (cols[2] == "mRNA" || cols[2] == "CDS") {
            records.push(raw);
        }
    }
    flush_ncbi_gene(&current_gene, &records, gene_num, &mut *out)?;
    Ok(0)
}

fn run_filter_gemoma(input: &str, output: Option<&str>) -> Result<i32> {
    fn flush_tx(current: &mut Vec<String>, len: &mut i64, bag: &mut Vec<(i64, Vec<String>)>) {
        if !current.is_empty() {
            bag.push((*len, std::mem::take(current)));
            *len = 0;
        }
    }

    fn flush_gene(
        gene_line: &Option<String>,
        transcripts: &mut Vec<(i64, Vec<String>)>,
        out: &mut dyn Write,
) -> Result<()> {
        let Some(gene_line) = gene_line else {
            return Ok(());
        };
        writeln!(out, "{}", gene_line)?;
        if transcripts.is_empty() {
            return Ok(());
        }
        transcripts.sort_by(|a, b| b.0.cmp(&a.0));
        let longest = &transcripts[0].1;
        for line in longest {
            writeln!(out, "{}", line)?;
        }
        transcripts.clear();
        Ok(())
    }

    let mut reader = open_reader(input)?;
    let mut out = open_writer(output)?;

    let mut current_gene: Option<String> = None;
    let mut transcripts: Vec<(i64, Vec<String>)> = Vec::new();
    let mut current_tx: Vec<String> = Vec::new();
    let mut current_len: i64 = 0;
    let mut tx_started = false;
    let mut line = String::new();

    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        match cols[2] {
            "gene" => {
                if tx_started {
                    flush_tx(&mut current_tx, &mut current_len, &mut transcripts);
                    tx_started = false;
                }
                flush_gene(&current_gene, &mut transcripts, &mut *out)?;
                current_gene = Some(raw);
            }
            _ if current_gene.is_none() => {}
            "mRNA" => {
                if tx_started {
                    flush_tx(&mut current_tx, &mut current_len, &mut transcripts);
                }
                current_tx.push(raw);
                current_len = 0;
                tx_started = true;
            }
            _ => {
                if tx_started {
                    current_tx.push(raw.clone());
                    if cols[2] == "CDS" {
                        if cols.len() > 4 {
                            if let (Ok(st), Ok(ed)) = (cols[3].parse::<i64>(), cols[4].parse::<i64>()) {
                                current_len += if ed > st { ed - st } else { 0 };
                            }
                        }
                    }
                }
            }
        }
    }
    if tx_started {
        flush_tx(&mut current_tx, &mut current_len, &mut transcripts);
    }
    flush_gene(&current_gene, &mut transcripts, &mut *out)?;
    Ok(0)
}

fn run_convert_ty1_hjjn(gff: &str, bed: &str, output: &str) -> Result<i32> {
    let mut bed_reader = open_reader(bed)?;
    let mut row = String::new();

    let mut mapping: HashMap<String, (String, i64, i64, i64, String)> = HashMap::new();
    let mut err_scaffolds: Vec<String> = Vec::new();

    while bed_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let fields: Vec<&str> = raw.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }
        let Ok(chr_start) = fields[1].parse::<i64>() else { continue };
        let Ok(sca_start) = fields[4].parse::<i64>() else { continue };
        let Ok(sca_end) = fields[5].parse::<i64>() else { continue };

        let direction = if fields.len() >= 8 && fields[6] == "1" {
            if fields[7] == "+" || fields[7] == "-" {
                fields[7].to_string()
            } else {
                "+".to_string()
            }
        } else if fields.len() >= 7 && (fields[6] == "+" || fields[6] == "-") {
            fields[6].to_string()
        } else {
            err_scaffolds.push(raw.to_string());
            "+".to_string()
        };

        mapping.insert(
            fields[3].to_string(),
            (fields[0].to_string(), chr_start, sca_start, sca_end, direction),
        );
    }

    let out_path = Path::new(output);
    let stem = out_path
        .file_stem()
        .unwrap_or_else(|| OsStr::new("Results"))
        .to_string_lossy();
    let parent = out_path.parent().unwrap_or_else(|| Path::new("."));

    let not_in = parent.join(format!("{stem}-not-in-bed.txt"));
    let err_file = parent.join(format!("{stem}-errors-scaffolds.txt"));
    let annot_log = parent.join("change-annot.log");
    let split1 = parent.join("change-gene-on-splitSca1.txt");
    let split2 = parent.join("change-gene-on-splitSca.txt");

    let mut out = BufWriter::new(File::create(expand_path(output))?);
    let mut not_in_writer = BufWriter::new(File::create(not_in)?);
    let mut err_writer = BufWriter::new(File::create(err_file)?);
    let mut log_writer = BufWriter::new(File::create(annot_log)?);
    let mut split_writer = BufWriter::new(File::create(split1)?);
    let _ = BufWriter::new(File::create(split2)?);

    for line in err_scaffolds.iter() {
        writeln!(err_writer, "{line}")?;
    }

    let mut gff_reader = open_reader(gff)?;
    let mut gline = String::new();

    while gff_reader.read_line(&mut gline)? > 0 {
        let raw = gline.trim_end_matches(['\n', '\r']).to_string();
        gline.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 8 {
            continue;
        }
        let ctg = cols[0];
        let Some((chr_name, chr_start, sca_start, sca_end, direction)) = mapping.get(ctg).cloned() else {
            writeln!(not_in_writer, "{raw}")?;
            continue;
        };

        let Ok(start) = cols[3].parse::<i64>() else { continue };
        let Ok(end) = cols[4].parse::<i64>() else { continue };
        if start < sca_start || end > sca_end {
            writeln!(split_writer, "{raw}")?;
            continue;
        }

        let mut out_cols: Vec<String> = cols.iter().map(|c| c.to_string()).collect();
        if direction == "+" {
            out_cols[3] = (start - sca_start + chr_start).to_string();
            out_cols[4] = (end - sca_start + chr_start).to_string();
        } else if direction == "-" {
            out_cols[3] = (sca_end - end + chr_start).to_string();
            out_cols[4] = (sca_end - start + chr_start).to_string();
            if out_cols.len() >= 7 {
                out_cols[6] = if out_cols[6] == "+" { "-".to_string() } else { "+".to_string() };
            }
        } else {
            writeln!(log_writer, "{raw}")?;
            continue;
        }
        out_cols[0] = chr_name;
        writeln!(out, "{}", out_cols.join("\t"))?;
    }

    for line in err_scaffolds {
        writeln!(log_writer, "{line}")?;
    }
    Ok(0)
}

fn extract_gene_name(header: &str) -> String {
    let find_key = |name: &str| -> Option<String> {
        for sep in ['=', ':'] {
            let token = format!("{name}{sep}");
            if let Some(mut pos) = header.find(&token) {
                pos += token.len();
                let bytes = header.as_bytes();
                let mut end = pos;
                while end < bytes.len() {
                    let ch = bytes[end] as char;
                    if ch.is_whitespace() || ch == ';' || ch == '|' {
                        break;
                    }
                    end += 1;
                }
                return Some(header[pos..end].to_string());
            }
        }
        None
    };

    if let Some(v) = find_key("gene") {
        return v;
    }
    if let Some(v) = find_key("gene_name") {
        return v;
    }

    let fields: Vec<&str> = header.split_whitespace().collect();
    if fields.len() > 3 {
        return fields[3].to_string();
    }
    if !fields.is_empty() {
        return fields[0].to_string();
    }
    header.to_string()
}

fn run_longest_transcript(fasta: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(fasta)?;
    let mut out = open_writer(output)?;
    let mut line = String::new();

    let mut selected: HashMap<String, (usize, String)> = HashMap::new();
    let mut seq_store: HashMap<String, Vec<String>> = HashMap::new();
    let mut order: Vec<String> = Vec::new();

    let mut current_header: Option<String> = None;
    let mut seq_lines: Vec<String> = Vec::new();

    let finalize_record = |header: String, seq: Vec<String>,
                           selected: &mut HashMap<String, (usize, String)>,
                           store: &mut HashMap<String, Vec<String>>,
                           order: &mut Vec<String>| {
        if header.trim().is_empty() {
            return;
        }
        let gene = extract_gene_name(&header);
        let len: usize = seq.iter().map(|s| s.len()).sum();
        let mut need_insert = false;
        match selected.get(&gene) {
            Some((best, _)) => {
                if len > *best {
                    need_insert = true;
                }
            }
            None => {
                need_insert = true;
                order.push(gene.clone());
            }
        }
        if need_insert {
            selected.insert(gene.clone(), (len, header.clone()));
            store.insert(header, seq);
        }
    };

    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.trim().is_empty() {
            continue;
        }
        if let Some(rest) = raw.strip_prefix('>') {
            if let Some(h) = current_header.take() {
                finalize_record(h, std::mem::take(&mut seq_lines), &mut selected, &mut seq_store, &mut order);
            }
            current_header = Some(rest.trim().to_string());
        } else {
            seq_lines.push(raw.trim().to_string());
        }
    }
    if let Some(h) = current_header.take() {
        finalize_record(h, std::mem::take(&mut seq_lines), &mut selected, &mut seq_store, &mut order);
    }

    for gene in &order {
        if let Some((_, header)) = selected.get(gene) {
            writeln!(out, ">{header}")?;
            if let Some(seqs) = seq_store.get(header) {
                for seq in seqs {
                    writeln!(out, "{seq}")?;
                }
            }
        }
    }
    Ok(0)
}

fn run_change_scaffolds_name(input: &str, name_list: &str, output: Option<&str>) -> Result<i32> {
    let mut mapping: HashMap<String, String> = HashMap::new();
    let mut mapping_reader = open_reader(name_list)?;
    let mut line = String::new();
    while mapping_reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        mapping.insert(cols[0].to_string(), cols[1].to_string());
    }

    let mut reader = open_reader(input)?;
    let mut out = open_writer(output)?;
    let mut row = String::new();
    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.is_empty() {
            continue;
        }
        let mut out_cols: Vec<&str> = cols;
        if let Some(new_name) = mapping.get(out_cols[0]) {
            out_cols[0] = new_name;
        }
        writeln!(out, "{}", out_cols.join("\t"))?;
    }
    Ok(0)
}

fn run_change_scaffolds_name_fasta(input: &str, name_list: &str, output: Option<&str>) -> Result<i32> {
    let mut name_map: HashMap<String, String> = HashMap::new();
    let mut map_order: Vec<String> = Vec::new();
    let mut name_reader = open_reader(name_list)?;
    let mut row = String::new();
    while name_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        let old_name = cols[1].to_string();
        let new_name = cols[0].to_string();
        if !name_map.contains_key(&old_name) {
            map_order.push(old_name.clone());
        }
        name_map.insert(old_name, new_name);
    }

    let mut seqs: HashMap<String, String> = HashMap::new();
    let mut reader = open_reader(input)?;
    let mut row = String::new();
    let mut current = String::new();
    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        if let Some(rest) = raw.strip_prefix('>') {
            current = rest.split_whitespace().next().unwrap_or("").to_string();
            seqs.entry(current.clone()).or_default();
        } else if !current.is_empty() {
            seqs.entry(current.clone()).or_default().push_str(&raw);
        }
    }

    let mut out = open_writer(output)?;
    for old_name in map_order {
        if let Some(seq) = seqs.get(&old_name) {
            if let Some(new_name) = name_map.get(&old_name) {
                writeln!(out, ">{}", new_name)?;
                writeln!(out, "{}", seq)?;
            }
        }
    }
    Ok(0)
}

fn run_change_seqname_for_fasta(input_dir: &str, name_list: &str, output_dir: &str) -> Result<i32> {
    let mut mapping: HashMap<String, String> = HashMap::new();
    let mut name_reader = open_reader(name_list)?;
    let mut line = String::new();
    while name_reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        let seq_name = format!(">{}", cols[1]);
        mapping.insert(seq_name, format!(">{}", cols[0]));
    }

    fs::create_dir_all(output_dir)?;
    for entry in fs::read_dir(expand_path(input_dir))? {
        let entry = entry?;
        let path = entry.path();
        let filename = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
        if !filename.ends_with(".aln") {
            continue;
        }
        let stem = Path::new(filename)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output")
            .to_string();
        let out_file = Path::new(output_dir).join(format!("{stem}.newName.fa"));
        let mut out = BufWriter::new(File::create(out_file)?);
        let mut reader = open_reader(&path.to_string_lossy())?;
        let mut seq_name = String::new();
        let mut row = String::new();
        while reader.read_line(&mut row)? > 0 {
            let raw = row.trim_end_matches(['\n', '\r']).to_string();
            row.clear();
            if raw.is_empty() {
                continue;
            }
            if let Some(rest) = raw.strip_prefix('>') {
                let key = format!(">{}", rest.split_whitespace().next().unwrap_or(""));
                let replaced = mapping.get(&key).cloned().unwrap_or(key);
                seq_name = replaced;
                writeln!(out, "{}", seq_name)?;
            } else if !seq_name.is_empty() {
                writeln!(out, "{}", raw)?;
            }
        }
    }
    Ok(0)
}

fn run_convert_3line2one(input: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(input)?;
    let mut out = open_writer(output)?;
    let mut line = String::new();
    let mut group = 0usize;
    let mut first_depth = String::new();
    let mut prefix: Vec<String> = Vec::new();
    let mut second_depth = String::new();

    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 4 {
            continue;
        }
        group = (group + 1) % 3;
        match group {
            1 => {
                first_depth = cols[3].to_string();
                prefix = cols[0..3].iter().map(|s| s.to_string()).collect();
            }
            2 => {
                second_depth = cols[3].to_string();
            }
            0 => {
                let out_line = format!(
                    "{}\t{}\t{}\t{}",
                    prefix.join("\t"),
                    first_depth,
                    second_depth,
                    cols[3]
                );
                writeln!(out, "{}", out_line)?;
            }
            _ => {}
        }
    }
    Ok(0)
}

fn run_filter_seq_by_length(input: &str, min_len: usize, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(input)?;
    let mut out = open_writer(output)?;
    let mut row = String::new();
    let mut header = String::new();
    let mut seq = String::new();

    let flush_record = |h: &mut String, s: &mut String, out: &mut dyn Write| -> Result<()> {
        if !h.is_empty() && s.len() >= min_len {
            writeln!(out, "{}", h)?;
            writeln!(out, "{}", s)?;
        }
        Ok(())
    };

    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        if raw.starts_with('>') {
            flush_record(&mut header, &mut seq, &mut *out)?;
            header = raw.to_string();
            seq.clear();
        } else {
            seq.push_str(&raw);
        }
    }
    flush_record(&mut header, &mut seq, &mut *out)?;
    Ok(0)
}

fn run_filter_gff_by_id(gff: &str, ids: &str, output: Option<&str>) -> Result<i32> {
    let mut id_set: HashSet<String> = HashSet::new();
    let mut id_reader = open_reader(ids)?;
    let mut row = String::new();
    while id_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).trim().to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        id_set.insert(raw.to_string());
    }

    let mut reader = open_reader(gff)?;
    let mut out = open_writer(output)?;
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        if raw.starts_with('#') {
            writeln!(out, "{}", raw)?;
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 9 {
            continue;
        }
        let first_attr = cols[8].split(';').next().unwrap_or("");
        if id_set.contains(first_attr) {
            writeln!(out, "{}", raw)?;
        }
    }
    Ok(0)
}

fn run_filter_gtf_ctg(input: &str, id_list: &str, output: &str) -> Result<i32> {
    let mut filtered: HashSet<String> = HashSet::new();
    let mut id_reader = open_reader(id_list)?;
    let mut line = String::new();
    while id_reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if !raw.is_empty() {
            filtered.insert(raw.to_string());
        }
    }

    let mut reader = open_reader(input)?;
    let mut out = open_writer(Some(output))?;
    let mut gline = String::new();
    while reader.read_line(&mut gline)? > 0 {
        let raw = gline.trim_end_matches(['\n', '\r']).to_string();
        gline.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.is_empty() {
            continue;
        }
        if filtered.contains(cols[0]) {
            continue;
        }
        writeln!(out, "{}", raw)?;
    }
    Ok(0)
}

fn run_merge_two_txt(input: &str, output: Option<&str>) -> Result<i32> {
    let mut best_map: HashMap<String, Vec<String>> = HashMap::new();
    let mut reader = open_reader(input)?;
    let mut row = String::new();
    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 3 {
            continue;
        }
        let key = cols[0].to_string();
        let score = cols[2].parse::<f64>().unwrap_or(0.0);
        if let Some(cur) = best_map.get(&key) {
            if let Ok(cur_score) = cur[2].parse::<f64>() {
                if score > cur_score {
                    best_map.insert(key, cols.iter().map(|v| v.to_string()).collect());
                }
            }
        } else {
            best_map.insert(key, cols.iter().map(|v| v.to_string()).collect());
        }
    }

    let mut out = open_writer(output)?;
    for fields in best_map.values() {
        writeln!(out, "{}", fields.join("\t"))?;
    }
    Ok(0)
}

fn run_compare_two_blast(primary: &str, reverse: &str, output: Option<&str>) -> Result<i32> {
    let mut primary_map: HashMap<String, String> = HashMap::new();
    let mut reverse_map: HashMap<String, String> = HashMap::new();
    let mut reader = open_reader(primary)?;
    let mut row = String::new();
    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        primary_map.insert(cols[0].to_string(), cols[1].to_string());
    }

    let mut rev_reader = open_reader(reverse)?;
    let mut rr = String::new();
    while rev_reader.read_line(&mut rr)? > 0 {
        let raw = rr.trim_end_matches(['\n', '\r']).to_string();
        rr.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        reverse_map.insert(cols[0].to_string(), cols[1].to_string());
    }

    let mut out = open_writer(output)?;
    for (query, target) in primary_map {
        if let Some(rev_target) = reverse_map.get(&target) {
            if rev_target == &query {
                writeln!(out, "{}\t{}", query, target)?;
            }
        }
    }
    Ok(0)
}

fn run_get_best_idy(input: &str, output: Option<&str>) -> Result<i32> {
    let mut best_map: HashMap<String, String> = HashMap::new();
    let mut reader = open_reader(input)?;
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 3 {
            continue;
        }
        let query = cols[0];
        let idy = cols[2];
        match best_map.get(query) {
            Some(best) => {
                if idy.parse::<f64>().unwrap_or(0.0) > best.parse::<f64>().unwrap_or(0.0) {
                    best_map.insert(query.to_string(), idy.to_string());
                }
            }
            None => {
                best_map.insert(query.to_string(), idy.to_string());
            }
        }
    }

    let mut out = open_writer(output)?;
    for (query, idy) in best_map {
        writeln!(out, "{query}\t{idy}")?;
    }
    Ok(0)
}

fn run_get_best_hit_based_on_idy(input: &str, output: Option<&str>) -> Result<i32> {
    let mut best_score: HashMap<String, String> = HashMap::new();
    let mut best_line: HashMap<String, String> = HashMap::new();
    let mut reader = open_reader(input)?;
    let mut row = String::new();
    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 12 {
            continue;
        }
        let query = cols[0];
        let score = cols[11];
        match best_score.get(query) {
            Some(best) => {
                if score.parse::<f64>().unwrap_or(0.0) > best.parse::<f64>().unwrap_or(0.0) {
                    best_score.insert(query.to_string(), score.to_string());
                    best_line.insert(query.to_string(), raw.to_string());
                }
            }
            None => {
                best_score.insert(query.to_string(), score.to_string());
                best_line.insert(query.to_string(), raw.to_string());
            }
        }
    }

    let mut out = open_writer(output)?;
    for line in best_line.values() {
        writeln!(out, "{line}")?;
    }
    Ok(0)
}

#[derive(Clone)]
struct BestHitState {
    best_score: f64,
    best_raw: String,
    best_target: String,
    second_score: f64,
    second_target: String,
}

fn run_get_best_hit_genes(input: &str, output: Option<&str>) -> Result<i32> {
    let mut e_dict: HashMap<String, BestHitState> = HashMap::new();
    let mut query_file = open_reader(input)?;
    let mut line = String::new();

    while query_file.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 11 {
            continue;
        }
        let gene_a = cols[0];
        let gene_b = cols[1];
        let evalue = cols[10].parse::<f64>().unwrap_or(f64::INFINITY);

        for (query, target) in &[(gene_a.to_string(), gene_b.to_string()), (gene_b.to_string(), gene_a.to_string())] {
            let state = e_dict.entry(query.clone()).or_insert(BestHitState {
                best_score: f64::INFINITY,
                best_raw: "100".to_string(),
                best_target: "Noop".to_string(),
                second_score: 100.0,
                second_target: "Noop".to_string(),
            });
            if evalue < state.best_score {
                state.second_score = state.best_score;
                state.second_target = state.best_target.clone();
                state.best_score = evalue;
                state.best_raw = cols[10].to_string();
                state.best_target = target.clone();
            } else if evalue < state.second_score {
                state.second_score = evalue;
                state.second_target = target.clone();
            }
        }
    }

    let mut out = open_writer(output)?;
    for (gene, state) in &e_dict {
        if let Some(best_target_state) = e_dict.get(&state.best_target) {
            if best_target_state.best_target == *gene {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    gene,
                    state.best_target,
                    state.best_raw,
                    state.second_score,
                    best_target_state.best_raw,
                    best_target_state.second_score,
                    state.second_target,
                    best_target_state.second_target
                )?;
            }
        }
    }
    Ok(0)
}

fn run_compare_as_and_no_as(no_as_file: &str, as_file: &str, output: Option<&str>) -> Result<i32> {
    let mut no_as: HashMap<String, String> = HashMap::new();
    let mut no_reader = open_reader(no_as_file)?;
    let mut row = String::new();
    while no_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 3 {
            continue;
        }
        no_as.insert(cols[0].to_string(), cols[2].to_string());
    }

    let mut as_reader = open_reader(as_file)?;
    let mut out = open_writer(output)?;
    row.clear();
    while as_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 3 {
            continue;
        }
        if cols[2] != "" && no_as.get(cols[0]).is_some_and(|v| v != "") {
            writeln!(out, "{}", raw)?;
        }
    }
    Ok(0)
}

fn run_compare_busco(asc_busco: &str, offspring_busco: &str) -> Result<i32> {
    let mut asc: Vec<String> = Vec::new();
    let mut off: Vec<String> = Vec::new();
    let mut asc_reader = open_reader(asc_busco)?;
    let mut row = String::new();
    while asc_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        if !asc.contains(&raw.to_string()) {
            asc.push(raw.to_string());
        }
    }

    let mut off_reader = open_reader(offspring_busco)?;
    while off_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        if !off.contains(&raw.to_string()) {
            off.push(raw.to_string());
        }
    }

    let mut just_have = BufWriter::new(File::create("justHave.txt")?);
    for item in &off {
        if !asc.contains(item) {
            writeln!(just_have, "{item}")?;
        }
    }

    let mut just_lost = BufWriter::new(File::create("justLost.txt")?);
    for item in &asc {
        if !off.contains(item) {
            writeln!(just_lost, "{item}")?;
        }
    }
    Ok(0)
}

fn run_merge_fpkm_file(input_dir: &str, output_fpkm: &str, output_profile: &str) -> Result<i32> {
    let mut records: BTreeMap<String, String> = BTreeMap::new();
    let mut exp_profile: HashMap<String, f64> = HashMap::new();
    let mut file_names: Vec<String> = Vec::new();

    let mut files: Vec<PathBuf> = fs::read_dir(expand_path(input_dir))?
        .map(|e| e.map(|f| f.path()))
        .collect::<std::result::Result<Vec<_>, io::Error>>()?;
    files.sort_by(|a, b| a.file_name().cmp(&b.file_name()));

    for path in &files {
        if !path.is_file() {
            continue;
        }
        let file_name = path.to_string_lossy().to_string();
        file_names.push(file_name.clone());
        let mut file_reader = open_reader(&file_name)?;
        let mut row = String::new();
        while file_reader.read_line(&mut row)? > 0 {
            let raw = row.trim_end_matches(['\n', '\r']).to_string();
            row.clear();
            if raw.is_empty() {
                continue;
            }
            let cols: Vec<&str> = raw.split_whitespace().collect();
            if cols.len() < 2 {
                continue;
            }
            let gene = cols[0];
            let val = cols[1];
            let entry = records.entry(gene.to_string()).or_default();
            if !entry.is_empty() {
                entry.push('\t');
            }
            entry.push_str(val);
            let expr = val.parse::<f64>().unwrap_or(0.0);
            *exp_profile.entry(gene.to_string()).or_insert(0.0) += expr;
        }
    }

    let mut out_fpkm = BufWriter::new(File::create(output_fpkm)?);
    let mut out_profile = BufWriter::new(File::create(output_profile)?);
    writeln!(out_fpkm, "gene_id\t{}", file_names.join("\t"))?;
    for (gene, seqs) in records {
        writeln!(out_fpkm, "{gene}\t{seqs}")?;
        let profile = exp_profile.remove(&gene).unwrap_or(0.0);
        writeln!(out_profile, "{gene}\t{profile}")?;
    }
    Ok(0)
}

fn run_save_go(input_file: &str, output: Option<&str>) -> Result<i32> {
    let mut out = open_writer(output)?;
    let mut reader = open_reader(input_file)?;
    let mut row = String::new();
    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.is_empty() {
            continue;
        }
        let gene = cols[0];
        let gos: Vec<String> = cols.iter().skip(1).filter_map(|x| {
            if x.contains("GO") {
                Some((*x).to_string())
            } else {
                None
            }
        }).collect();
        writeln!(out, "{}\t{}", gene, gos.join("\t"))?;
    }
    Ok(0)
}

fn run_merge_gos(swissprot: &str, nr: &str, trembl: &str, output: Option<&str>) -> Result<i32> {
    fn load_go(source: &str, map: &mut HashMap<String, Vec<String>>) {
        if let Ok(mut reader) = open_reader(source) {
            let mut row = String::new();
            let mut seen: HashSet<String> = HashSet::new();
            while let Ok(size) = reader.read_line(&mut row) {
                if size == 0 {
                    break;
                }
                let raw = row.trim_end_matches(['\n', '\r']).to_string();
                row.clear();
                if raw.is_empty() || raw.starts_with('#') {
                    continue;
                }
                let cols: Vec<&str> = raw.split_whitespace().collect();
                if cols.len() < 2 {
                    continue;
                }
                let gene = cols[0];
                let go = cols[1];
                if go.starts_with("EC") {
                    continue;
                }
                let row = map.entry(gene.to_string()).or_default();
                let key = format!("{gene}\t{go}");
                if !seen.contains(&key) {
                    seen.insert(key);
                    row.push(go.to_string());
                }
            }
        }
    }

    let mut go_map: HashMap<String, Vec<String>> = HashMap::new();
    load_go(swissprot, &mut go_map);
    load_go(nr, &mut go_map);
    load_go(trembl, &mut go_map);

    let mut out = open_writer(output)?;
    for (gene, gos) in go_map {
        for go_name in gos {
            writeln!(out, "{gene}\t{go_name}")?;
        }
    }
    Ok(0)
}

fn run_merge_blastp_best_jcvi(path_dir: &str, output: Option<&str>) -> Result<i32> {
    let mut files: Vec<PathBuf> = fs::read_dir(expand_path(path_dir))?
        .map(|e| e.map(|f| f.path()))
        .collect::<std::result::Result<Vec<_>, io::Error>>()?;
    files.sort_by(|a, b| a.file_name().cmp(&b.file_name()));

    let mut gene_order: Vec<String> = Vec::new();
    let mut gene_rows: HashMap<String, Vec<String>> = HashMap::new();
    let mut file_count = 0usize;

    for path in &files {
        if !path.is_file() {
            continue;
        }
        let mut reader = open_reader(&path.to_string_lossy())?;
        let mut row = String::new();
        while reader.read_line(&mut row)? > 0 {
            let raw = row.trim_end_matches(['\n', '\r']).to_string();
            row.clear();
            if raw.is_empty() {
                continue;
            }
            let cols: Vec<&str> = raw.split_whitespace().collect();
            if cols.len() < 2 {
                continue;
            }
            let orthogene = cols[0].to_string();
            let gene = cols[1].to_string();
            let item = gene_rows.entry(gene.clone()).or_insert_with(|| {
                gene_order.push(gene.clone());
                Vec::new()
            });
            if item.len() < file_count {
                item.resize(file_count, ".".to_string());
            }
            item.push(orthogene);
        }

        file_count += 1;
        for gene in &gene_order {
            if let Some(item) = gene_rows.get_mut(gene) {
                if item.len() < file_count {
                    item.resize(file_count, ".".to_string());
                }
            }
        }
    }

    let mut out = open_writer(output)?;
    for gene in gene_order {
        if let Some(items) = gene_rows.get(&gene) {
            writeln!(out, "{}\t{}", gene, items.join("\t"))?;
        }
    }
    Ok(0)
}

fn run_orthogenes(input: &str, output: &str) -> Result<i32> {
    let mut reader = open_reader(input)?;
    let mut out = open_writer(Some(output))?;
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let mut split = raw.splitn(2, ':');
        let group_name = split.next().unwrap_or("").trim();
        let members = split.next().unwrap_or("");
        let mut uniq: Vec<String> = Vec::new();
        let mut seen: HashSet<String> = HashSet::new();
        for member in members.split(',') {
            let m = member.trim();
            if m.is_empty() {
                continue;
            }
            if seen.insert(m.to_string()) {
                uniq.push(m.to_string());
            }
        }
        if uniq.len() > 1 && !group_name.is_empty() {
            writeln!(out, "{group_name}:{}", uniq.join(","))?;
        }
    }
    Ok(0)
}

fn run_genome_gc(input: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(input)?;
    let mut row = String::new();
    let mut gc_count = 0usize;
    let mut atcg_count = 0usize;
    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        for ch in raw.chars() {
            if ch == '>' {
                continue;
            }
            match ch {
                'G' | 'g' | 'C' | 'c' => {
                    gc_count += 1;
                    atcg_count += 1;
                }
                'A' | 'a' | 'T' | 't' => atcg_count += 1,
                _ => {}
            }
        }
    }
    let ratio = if atcg_count == 0 {
        0.0
    } else {
        gc_count as f64 / atcg_count as f64
    };
    let mut out = open_writer(output)?;
    writeln!(out, "{ratio}")?;
    Ok(0)
}

fn run_get_the_longest_seq(input: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(input)?;
    let mut out = open_writer(output)?;
    let mut row = String::new();
    let mut lengths: HashMap<String, usize> = HashMap::new();
    let flush = |name: &str, seq_len: usize, lengths: &mut HashMap<String, usize>| -> Option<(String, String)> {
        if name.is_empty() {
            return None;
        }
        let len_parts: Vec<&str> = name.split('.').collect();
        if len_parts.len() < 2 {
            return None;
        }
        let gene_name = len_parts[1].to_string();
        let trans_name = len_parts[0].to_string();
        if let Some(last_len) = lengths.get(&gene_name) {
            if seq_len > *last_len {
                lengths.insert(gene_name.clone(), seq_len);
                return Some((gene_name, trans_name));
            }
            None
        } else {
            lengths.insert(gene_name.clone(), seq_len);
            Some((gene_name, trans_name))
        }
    };

    let mut selected: HashMap<String, String> = HashMap::new();
    let mut current_name = String::new();
    let mut current_len = 0usize;
    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        if let Some(rest) = raw.strip_prefix('>') {
            if !current_name.is_empty() {
                if let Some((gene, trans)) = flush(&current_name, current_len, &mut lengths) {
                    selected.insert(gene, trans);
                }
            }
            current_name = rest.to_string();
            current_len = 0;
        } else {
            current_len += raw.len();
        }
    }
    if !current_name.is_empty() {
        if let Some((gene, trans)) = flush(&current_name, current_len, &mut lengths) {
            selected.insert(gene, trans);
        }
    }

    for (gene, trans) in selected {
        writeln!(out, "{trans}\t{gene}")?;
    }
    Ok(0)
}

fn run_extract_longest_pep_from_ensembl_download(input: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(input)?;
    let mut row = String::new();
    let mut current_name = String::new();
    let mut records: HashMap<String, Vec<String>> = HashMap::new();
    let mut dup_count = 0usize;
    let mut current_seq = String::new();

    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        if let Some(rest) = raw.strip_prefix('>') {
            if !current_name.is_empty() {
                let key = if records.contains_key(rest) {
                    let v = format!("{rest}-{dup_count}");
                    dup_count += 1;
                    v
                } else {
                    rest.to_string()
                };
                records.insert(key, vec![current_seq.clone()]);
            }
            current_name = rest.to_string();
            current_seq.clear();
        } else {
            current_seq.push_str(&raw);
        }
    }
    if !current_name.is_empty() {
        let key = if records.contains_key(current_name.as_str()) {
            let v = format!("{current_name}-{dup_count}");
            v
        } else {
            current_name.clone()
        };
        records.insert(key, vec![current_seq]);
    }

    let mut longest: HashMap<String, (String, usize)> = HashMap::new();
    for (name, seqs) in records {
        let seq = seqs.join("");
        let seq_len = seq.len();
        let true_name = name.split('-').next().unwrap_or(&name).to_string();
        if let Some((best_seq, best_len)) = longest.get(&true_name) {
            if seq_len > *best_len {
                longest.insert(true_name, (seq, seq_len));
            } else if best_seq.is_empty() {
                longest.insert(true_name, (seq, seq_len));
            }
        } else {
            longest.insert(true_name, (seq, seq_len));
        }
    }

    let mut out = open_writer(output)?;
    for (name, (seq, _)) in longest {
        writeln!(out, ">{name}")?;
        writeln!(out, "{seq}")?;
    }
    Ok(0)
}

fn run_convert_lastz2_jcvi(
    bed_file: &str,
    _ref_len_file: &str,
    _query_len_file: &str,
    ref_name: &str,
    query_name: &str,
) -> Result<i32> {
    let mut bed_reader = open_reader(bed_file)?;
    let mut row = String::new();
    let mut query_out = BufWriter::new(File::create("query.bed")?);
    let mut ref_out = BufWriter::new(File::create("ref.bed")?);
    let mut simple_out = BufWriter::new(File::create("ref_query.simple")?);

    let mut query_num = 0usize;
    let mut ref_num = 0usize;
    let gene_len = 100i64;
    let gap = 500i64;

    let emit_segments = |
        chrom: &str,
        mut start: i64,
        end: i64,
        strand: &str,
        name_prefix: &str,
        counter: &mut usize,
        out: &mut BufWriter<File>,
    | -> Result<(String, String)> {
        let mut first_name = String::new();
        let mut names: Vec<String> = Vec::new();

        if start >= end {
            let name = format!("{chrom}{name_prefix}.gene.psuedgene.{counter}");
            *counter += 1;
            writeln!(out, "{chrom}\t{start}\t{end}\t{name}\t0\t{strand}")?;
            return Ok((name.clone(), name));
        }

        loop {
            let this_end = start + gene_len;
            let name = format!("{chrom}{name_prefix}.gene.psuedgene.{counter}");
            if first_name.is_empty() {
                first_name = name.clone();
            }
            names.push(name.clone());
            *counter += 1;

            if this_end < end && start < end {
                writeln!(out, "{chrom}\t{start}\t{this_end}\t{name}\t0\t{strand}")?;
                start = this_end + gap;
                continue;
            }

            if this_end > end && start < end {
                writeln!(out, "{chrom}\t{start}\t{end}\t{name}\t0\t{strand}")?;
                break;
            }

            let adj_start = start - gap;
            writeln!(out, "{chrom}\t{adj_start}\t{end}\t{name}\t0\t{strand}")?;
            break;
        }

        let last_name = names.last().cloned().unwrap_or_else(|| first_name.clone());
        Ok((first_name, last_name))
    };

    while bed_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 7 {
            continue;
        }

        let q_chrom = cols[0];
        let q_start = cols[1].parse::<i64>()?;
        let q_end = cols[2].parse::<i64>()?;
        let r_chrom = cols[3];
        let mut r_start = cols[4].parse::<i64>()?;
        let mut r_end = cols[5].parse::<i64>()?;
        let strand = cols[6];

        if r_start > r_end {
            std::mem::swap(&mut r_start, &mut r_end);
        }

        let (q_first, q_last) = emit_segments(
            q_chrom,
            q_start,
            q_end,
            strand,
            ref_name,
            &mut query_num,
            &mut query_out,
        )?;
        let (r_first, r_last) = emit_segments(
            r_chrom,
            r_start,
            r_end,
            strand,
            query_name,
            &mut ref_num,
            &mut ref_out,
        )?;
        writeln!(simple_out, "{q_first}\t{q_last}\t{r_first}\t{r_last}\t322\t{strand}")?;
    }

    Ok(0)
}

fn run_extract_pasa_results(input: &str, out_seq: &str, out_gff: &str) -> Result<i32> {
    let mut reader = open_reader(input)?;
    let mut out_gff_writer = BufWriter::new(File::create(out_gff)?);
    let mut out_seq_writer = BufWriter::new(File::create(out_seq)?);

    let mut row = String::new();
    let mut in_gene = false;
    let mut current_mrna: Vec<String> = Vec::new();
    let mut best_mrna: Vec<String> = Vec::new();
    let mut prot_lines: Vec<String> = Vec::new();

    let flush_gene = |
        out: &mut dyn Write,
        current_mrna: &mut Vec<String>,
        best_mrna: &mut Vec<String>,
    | -> Result<()> {
        if current_mrna.len() > best_mrna.len() {
            *best_mrna = std::mem::take(current_mrna);
        }
        for line in best_mrna.iter() {
            writeln!(out, "{line}")?;
        }
        best_mrna.clear();
        Ok(())
    };

    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() || raw.starts_with('#') {
            prot_lines.push(raw);
            continue;
        }

        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.is_empty() {
            continue;
        }

        let feature = cols[2];
        if feature == "gene" {
            if in_gene {
                flush_gene(
                    &mut out_gff_writer,
                    &mut current_mrna,
                    &mut best_mrna,
                )?;
            }
            writeln!(out_gff_writer, "{raw}")?;
            in_gene = true;
            current_mrna.clear();
            continue;
        }

        if !in_gene {
            continue;
        }

        if feature == "mRNA" {
            if current_mrna.len() > best_mrna.len() {
                best_mrna = std::mem::take(&mut current_mrna);
            } else {
                current_mrna.clear();
            }
            current_mrna.push(raw);
            continue;
        }

        current_mrna.push(raw);
    }

    if in_gene {
        flush_gene(
            &mut out_gff_writer,
            &mut current_mrna,
            &mut best_mrna,
        )?;
    }

    for prot in prot_lines {
        if prot.starts_with("#PROT") {
            let fields: Vec<&str> = prot.split_whitespace().collect();
            if fields.len() >= 4 {
                writeln!(out_seq_writer, ">{}-{}\n{}", fields[1], fields[2], fields[3])?;
            }
        }
    }

    Ok(0)
}

fn run_convert_gemoma_gff3(input_gff: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(input_gff)?;
    let mut out = open_writer(output)?;
    let mut line = String::new();

    let mut gene_id = 1usize;
    let mut mrna_id = 1usize;
    let mut exon_id = 1usize;
    let mut current_mrna = String::new();

    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = raw.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let feature = fields[2];
        match feature {
            "gene" => {
                let g_id = format!("Plants{}gene{gene_id:05}", fields[0].to_uppercase());
                let g_name = format!("Plant{gene_id:05}");
                writeln!(out, "{}\tID={g_id};Name={g_name}", fields[..8].join("\t"))?;
                current_mrna.clear();
                gene_id += 1;
                mrna_id = 1;
                exon_id = 1;
            }
            "mRNA" => {
                let m_id = format!("Plants{}gene{:05}.{}", fields[0].to_uppercase(), gene_id - 1, mrna_id);
                let m_name = format!("Plant{:05}.{}", gene_id - 1, mrna_id);
                writeln!(
                    out,
                    "{}\tID={m_id};Parent=Plant{:05};Name={m_name}",
                    fields[..8].join("\t"),
                    gene_id - 1
                )?;
                current_mrna = m_id;
                mrna_id += 1;
            }
            _ => {
                let exon = format!("{}.exon{exon_id}", current_mrna);
                let cds = format!("cds.{current_mrna}");
                let mut exon_cols: Vec<&str> = Vec::with_capacity(8);
                exon_cols.extend_from_slice(&fields[..2]);
                exon_cols.push("exon");
                exon_cols.extend_from_slice(&fields[3..8]);
                writeln!(out, "{}\tID={exon};Parent={}", exon_cols.join("\t"), current_mrna)?;
                writeln!(out, "{}\tID={cds};Parent={}", fields[..8].join("\t"), current_mrna)?;
                exon_id += 1;
            }
        }
    }
    Ok(0)
}

fn run_convert_gene_annotation_contigs2chr_pasa(
    gff_file: &str,
    background_file: &str,
    output: Option<&str>,
) -> Result<i32> {
    let mut bg_reader = open_reader(background_file)?;
    let mut row = String::new();
    let mut map: HashMap<String, Vec<String>> = HashMap::new();

    while bg_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let mut fields: Vec<String> = raw.split_whitespace().map(|v| v.to_string()).collect();
        if fields.len() < 5 {
            continue;
        }
        let ctg = fields[3].clone();
        fields.remove(3);
        map.insert(ctg, fields);
    }

    let mut gff_reader = open_reader(gff_file)?;
    let mut out = open_writer(output)?;
    let mut gline = String::new();

    while gff_reader.read_line(&mut gline)? > 0 {
        let raw = gline.trim_end_matches(['\n', '\r']).to_string();
        gline.clear();
        if raw.is_empty() {
            continue;
        }
        if raw.starts_with('#') {
            writeln!(out, "{raw}")?;
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }

        let ctg = cols[0];
        let mut out_cols: Vec<String> = cols.iter().map(|x| x.to_string()).collect();
        let info = match map.get(ctg) {
            Some(v) => v,
            None => {
                writeln!(out, "{raw}")?;
                continue;
            }
        };
        if info.len() < 5 {
            writeln!(out, "{raw}")?;
            continue;
        }

        let chr = &info[0];
        let chr_start = info[1].parse::<i64>().unwrap_or(0);
        let orientation = info[3].as_str();
        let ctg_end = info[4].parse::<i64>().unwrap_or(0);
        let s = cols[3].parse::<i64>().unwrap_or(0);
        let e = cols[4].parse::<i64>().unwrap_or(0);

        out_cols[0] = chr.clone();
        match orientation {
            "0" => {
                out_cols[3] = (s + chr_start - 1).to_string();
                out_cols[4] = (e + chr_start - 1).to_string();
            }
            "1" => {
                out_cols[3] = (ctg_end - e + 1 + chr_start - 1).to_string();
                out_cols[4] = (ctg_end - s + 1 + chr_start - 1).to_string();
                if out_cols.len() >= 7 {
                    out_cols[6] = if out_cols[6] == "+" {
                        "-".to_string()
                    } else {
                        "+".to_string()
                    };
                }
            }
            _ => {
                writeln!(out, "{raw}")?;
                continue;
            }
        }
        writeln!(out, "{}", out_cols.join("\t"))?;
    }
    Ok(0)
}

fn run_convert_gene_annotation_scaffold2chr_nextgenomics(
    gff_file: &str,
    background_file: &str,
    output: Option<&str>,
) -> Result<i32> {
    let mut bg_reader = open_reader(background_file)?;
    let mut row = String::new();
    let mut map: HashMap<String, Vec<String>> = HashMap::new();
    let mut split_scaf: Vec<String> = Vec::new();

    while bg_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 7 {
            continue;
        }
        if cols[6] != "1" {
            split_scaf.push(raw);
            continue;
        }
        let mut kept: Vec<String> = Vec::new();
        for (idx, c) in cols.iter().enumerate() {
            if idx == 3 || idx == 4 {
                continue;
            }
            kept.push(c.to_string());
        }
        if cols.len() > 3 {
            map.insert(cols[3].to_string(), kept);
        }
    }

    let mut out = open_writer(output.or(Some("Change-annot.gff3")))?;
    let mut no_change: Vec<String> = Vec::new();
    let mut no_change_log = BufWriter::new(File::create("change-gene-on-splitSca.txt")?);
    let mut no_change_log1 = BufWriter::new(File::create("change-gene-on-splitSca1.txt")?);
    let mut log = BufWriter::new(File::create("change-annot.log")?);
    let mut err_scaf = BufWriter::new(File::create("change-errors-scaffolds.txt")?);

    for line in &split_scaf {
        writeln!(err_scaf, "{line}")?;
    }

    let mut gff_reader = open_reader(gff_file)?;
    let mut gline = String::new();
    while gff_reader.read_line(&mut gline)? > 0 {
        let raw = gline.trim_end_matches(['\n', '\r']).to_string();
        gline.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let ctg = cols[0];
        let info = match map.get(ctg) {
            Some(v) => v,
            None => {
                writeln!(log, "{raw}")?;
                continue;
            }
        };
        if info.len() < 6 {
            no_change.push(raw);
            continue;
        }
        let chr = &info[0];
        let chr_start = info[1].parse::<i64>().unwrap_or(0);
        let chr_dir = info[3].as_str();
        let sca_start = info[4].parse::<i64>().unwrap_or(0);
        let sca_end = info[5].parse::<i64>().unwrap_or(0);
        let start = cols[3].parse::<i64>().unwrap_or(0);
        let end = cols[4].parse::<i64>().unwrap_or(0);

        let mut out_cols: Vec<String> = cols.iter().map(|x| x.to_string()).collect();
        if (chr_dir == "+" && start >= sca_start && end <= sca_end)
            || (chr_dir == "-" && start >= sca_start && end <= sca_end)
        {
            out_cols[0] = chr.clone();
            if chr_dir == "+" {
                out_cols[3] = (start - sca_start + chr_start).to_string();
                out_cols[4] = (end - sca_start + chr_start).to_string();
            } else {
                out_cols[3] = (sca_end - end + chr_start).to_string();
                out_cols[4] = (sca_end - start + chr_start).to_string();
                if out_cols.len() >= 7 {
                    out_cols[6] = if out_cols[6] == "+" {
                        "-".to_string()
                    } else {
                        "+".to_string()
                    };
                }
            }
            writeln!(out, "{}", out_cols.join("\t"))?;
        } else {
            no_change.push(raw);
        }
    }

    for n in &no_change {
        writeln!(no_change_log1, "{n}")?;
        writeln!(no_change_log1, "q")?;
    }

    // Try one extra-pass for records on split scaffolds.
    for err_line in split_scaf {
        let e = err_line.split('\t').collect::<Vec<_>>();
        if e.len() < 8 {
            continue;
        }
        let err_ctg = e[3];
        let chr = e[0];
        let chr_dir = e[5];
        let sca_s = e[6].parse::<i64>().unwrap_or(0);
        let sca_e = e[7].parse::<i64>().unwrap_or(0);
        let chr_s = e[1].parse::<i64>().unwrap_or(0);

        let mut remaining = Vec::new();
        for n in no_change {
            let n_cols: Vec<&str> = n.split('\t').collect();
            if n_cols.len() < 5 || n_cols[0] != err_ctg {
                remaining.push(n);
                continue;
            }
            let ns = n_cols[3].parse::<i64>().unwrap_or(0);
            let ne = n_cols[4].parse::<i64>().unwrap_or(0);
            if ns < sca_s || ne > sca_e {
                remaining.push(n.to_string());
                continue;
            }
            if chr_dir == "+" {
                writeln!(
                    out,
                    "{chr}\t{}\t{}\t{}\t{}\t{}",
                    n_cols[1],
                    n_cols[2],
                    ns - sca_s + chr_s,
                    ne - sca_s + chr_s,
                    n_cols[5..].join("\t")
                )?;
            } else {
                let mut ncols: Vec<String> = n_cols.iter().map(|v| (*v).to_string()).collect();
                if ncols.len() >= 7 {
                    ncols[3] = (sca_e - ne + chr_s).to_string();
                    ncols[4] = (sca_e - ns + chr_s).to_string();
                    if ncols[6] == "-" {
                        ncols[6] = "+".to_string();
                    } else {
                        ncols[6] = "-".to_string();
                    }
                }
                ncols[0] = chr.to_string();
                writeln!(out, "{}", ncols.join("\t"))?;
            }
        }
        remaining.retain(|x| !x.is_empty());
        no_change = remaining;
    }

    for n in no_change {
        writeln!(no_change_log, "{n}")?;
    }
    Ok(0)
}

fn run_filter_gemoma_as(input_gff: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(input_gff)?;
    let mut out = open_writer(output)?;
    let mut line = String::new();

    let mut in_gene = false;
    let mut best_mrna: Vec<String> = Vec::new();
    let mut current_mrna: Vec<String> = Vec::new();
    let mut best_len = 0i64;
    let mut cur_len = 0i64;

    let flush_gene = |
        out: &mut dyn Write,
        in_gene: &mut bool,
        best_mrna: &mut Vec<String>,
        best_len: &mut i64,
        cur_mrna: &mut Vec<String>,
        cur_len: &mut i64,
    | -> Result<()> {
        if !*in_gene {
            return Ok(());
        }
        if *cur_len > *best_len {
            *best_len = *cur_len;
            *best_mrna = std::mem::take(cur_mrna);
        }
        for l in best_mrna.iter() {
            writeln!(out, "{l}")?;
        }
        best_mrna.clear();
        *best_len = 0;
        *cur_len = 0;
        *in_gene = false;
        Ok(())
    };

    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        if raw.starts_with('#') {
            writeln!(out, "{raw}")?;
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let feature = cols[2];
        if feature == "gene" {
            flush_gene(
                &mut out,
                &mut in_gene,
                &mut best_mrna,
                &mut best_len,
                &mut current_mrna,
                &mut cur_len,
            )?;
            writeln!(out, "{raw}")?;
            in_gene = true;
            continue;
        }

        if !in_gene {
            continue;
        }

        if feature == "mRNA" {
            if cur_len > best_len {
                best_len = cur_len;
                best_mrna = std::mem::take(&mut current_mrna);
            } else {
                current_mrna.clear();
                cur_len = 0;
            }
            current_mrna.push(raw);
            continue;
        }

        if feature == "CDS" {
            if cols.len() > 4 {
                let s = cols[3].parse::<i64>().unwrap_or(0);
                let e = cols[4].parse::<i64>().unwrap_or(0);
                cur_len += e - s + 1;
            }
        }
        current_mrna.push(raw);
    }

    flush_gene(
        &mut out,
        &mut in_gene,
        &mut best_mrna,
        &mut best_len,
        &mut current_mrna,
        &mut cur_len,
    )?;
    Ok(0)
}

fn run_filter_gemoma_as2(input_gff: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(input_gff)?;
    let mut out = open_writer(output)?;
    let mut line = String::new();

    let mut current_gene: Option<String> = None;
    let mut transcripts: Vec<(Vec<String>, i64)> = Vec::new();
    let mut current_lines: Vec<String> = Vec::new();
    let mut current_len = 0i64;

    let flush = |
        out: &mut dyn Write,
        gene_line: &mut Option<String>,
        transcripts: &mut Vec<(Vec<String>, i64)>,
    | -> Result<()> {
        let Some(gene_line) = gene_line.take() else {
            return Ok(());
        };
        writeln!(out, "{gene_line}")?;
        if transcripts.is_empty() {
            return Ok(());
        }
        let mut best_idx = 0usize;
        let mut best_len = -1i64;
        for (idx, (_, len)) in transcripts.iter().enumerate() {
            if *len > best_len {
                best_len = *len;
                best_idx = idx;
            }
        }
        for line in transcripts[best_idx].0.iter() {
            writeln!(out, "{line}")?;
        }
        transcripts.clear();
        Ok(())
    };

    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        if raw.starts_with('#') {
            writeln!(out, "{raw}")?;
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let feature = cols[2];

        if feature == "gene" {
            flush(&mut out, &mut current_gene, &mut transcripts)?;
            current_gene = Some(raw);
            transcripts.clear();
            current_lines.clear();
            current_len = 0;
            continue;
        }

        let Some(_) = current_gene else {
            continue;
        };

        if feature == "mRNA" {
            if !current_lines.is_empty() {
                transcripts.push((std::mem::take(&mut current_lines), current_len));
            }
            current_lines.push(raw);
            current_len = 0;
            continue;
        }

        if feature == "CDS" && cols.len() > 4 {
            let s = cols[3].parse::<i64>().unwrap_or(0);
            let e = cols[4].parse::<i64>().unwrap_or(0);
            current_len += (e - s).abs();
        }
        current_lines.push(raw);
    }

    if current_gene.is_some() && !current_lines.is_empty() {
        transcripts.push((std::mem::take(&mut current_lines), current_len));
    }
    flush(&mut out, &mut current_gene, &mut transcripts)?;
    Ok(0)
}

fn run_get_best_hit_by_score(
    query_file: &str,
    ref_file: &str,
    output: Option<&str>,
) -> Result<i32> {
    let mut q_reader = open_reader(query_file)?;
    let mut row = String::new();
    let mut qbest: HashMap<String, (f64, String)> = HashMap::new();
    while q_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 12 {
            continue;
        }
        let q = cols[0];
        let t = cols[1];
        let score = cols[11].parse::<f64>().unwrap_or(f64::NEG_INFINITY);
        qbest
            .entry(q.to_string())
            .and_modify(|(s, tgt)| {
                if score > *s {
                    *s = score;
                    *tgt = t.to_string();
                }
            })
            .or_insert((score, t.to_string()));
    }

    let mut r_reader = open_reader(ref_file)?;
    let mut rrow = String::new();
    let mut rbest: HashMap<String, (f64, String)> = HashMap::new();
    while r_reader.read_line(&mut rrow)? > 0 {
        let raw = rrow.trim_end_matches(['\n', '\r']).to_string();
        rrow.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 12 {
            continue;
        }
        let q = cols[0];
        let t = cols[1];
        let score = cols[11].parse::<f64>().unwrap_or(f64::NEG_INFINITY);
        rbest
            .entry(q.to_string())
            .and_modify(|(s, tgt)| {
                if score > *s {
                    *s = score;
                    *tgt = t.to_string();
                }
            })
            .or_insert((score, t.to_string()));
    }

    let mut out = open_writer(output)?;
    for (q, (_score, t)) in qbest {
        if let Some((_r, rev_q)) = rbest.get(&t) {
            if rev_q == &q {
                writeln!(out, "{q}\t{t}")?;
            }
        }
    }
    Ok(0)
}

fn run_get_best_hit_by_score_one_file(
    blast_file: &str,
    out_prefix: &str,
    output: Option<&str>,
) -> Result<i32> {
    let mut reader = open_reader(blast_file)?;
    let mut row = String::new();

    let mut qbest: HashMap<String, (f64, String, String)> = HashMap::new();
    let mut rbest: HashMap<String, (f64, String)> = HashMap::new();
    let mut lines: Vec<String> = Vec::new();

    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        lines.push(raw);
    }

    for raw in &lines {
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 12 {
            continue;
        }
        let q = cols[0];
        let t = cols[1];
        let score = cols[11].parse::<f64>().unwrap_or(f64::NEG_INFINITY);
        let idy = cols[2].to_string();
        qbest
            .entry(q.to_string())
            .and_modify(|(s, target, id)| {
                if score > *s {
                    *s = score;
                    *target = t.to_string();
                    *id = idy.to_string();
                }
            })
            .or_insert((score, t.to_string(), idy));
    }
    for raw in &lines {
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 12 {
            continue;
        }
        let q = cols[0];
        let t = cols[1];
        let score = cols[11].parse::<f64>().unwrap_or(f64::NEG_INFINITY);
        rbest
            .entry(t.to_string())
            .and_modify(|(s, qv)| {
                if score > *s {
                    *s = score;
                    *qv = q.to_string();
                }
            })
            .or_insert((score, q.to_string()));
    }

    let out_file = output.unwrap_or(&format!("{out_prefix}.idy.txt")).to_string();
    let mut out = BufWriter::new(File::create(out_file)?);
    for (q, (_s, t, idy)) in qbest {
        if let Some((_score, rev_q)) = rbest.get(&t) {
            if rev_q == &q {
                writeln!(out, "{q}\t{t}\t{idy}")?;
            }
        }
    }
    Ok(0)
}

fn run_get_best_hit_from_blast(
    path_dir: &str,
    num_species: usize,
    output: Option<&str>,
) -> Result<i32> {
    let mut all_map: HashMap<String, Vec<String>> = HashMap::new();
    for entry in fs::read_dir(expand_path(path_dir))? {
        let entry = entry?;
        if !entry.path().is_file() {
            continue;
        }
        let file = entry.path().to_string_lossy().to_string();
        let mut reader = open_reader(&file)?;
        let mut row = String::new();
        while reader.read_line(&mut row)? > 0 {
            let raw = row.trim_end_matches(['\n', '\r']).to_string();
            row.clear();
            if raw.is_empty() {
                continue;
            }
            let cols: Vec<&str> = raw.split_whitespace().collect();
            if cols.len() < 2 {
                continue;
            }
            all_map.entry(cols[0].to_string()).or_default().push(cols[1].to_string());
        }
    }

    let need = num_species.saturating_sub(1);
    let mut temp: HashMap<String, Vec<String>> = HashMap::new();
    for (q, ts) in all_map {
        if ts.len() == need {
            temp.insert(q, ts);
        }
    }

    let mut filtered: HashMap<String, Vec<String>> = HashMap::new();
    for (q, ts) in temp.iter() {
        let mut keep = true;
        for t in ts {
            match temp.get(t) {
                Some(v) if v.contains(q) => {}
                _ => {
                    keep = false;
                    break;
                }
            }
        }
        if keep {
            filtered.insert(q.clone(), ts.clone());
        }
    }

    let mut to_remove: HashSet<String> = HashSet::new();
    let keys: Vec<String> = filtered.keys().cloned().collect();
    for k in keys {
        if to_remove.contains(&k) {
            continue;
        }
        if let Some(ts) = filtered.get(&k) {
            for t in ts {
                if filtered.contains_key(t) {
                    to_remove.insert(t.clone());
                }
            }
        }
    }
    for k in to_remove {
        filtered.remove(&k);
    }

    let mut out = open_writer(output)?;
    for (q, ts) in filtered {
        writeln!(out, "{q}\t{}", ts.join("\t"))?;
    }
    Ok(0)
}

fn run_extract_gene_family_info(
    gene_len_file: &str,
    expr_file: &str,
    parent_file: &str,
    coverage: f64,
    output: Option<&str>,
) -> Result<i32> {
    let mut gene_len: HashMap<String, f64> = HashMap::new();
    let mut parent_map: HashMap<String, String> = HashMap::new();
    let mut expr_map: HashMap<String, Vec<String>> = HashMap::new();
    let mut expr_row = String::new();

    let out_dir = output.unwrap_or(".");
    fs::create_dir_all(expand_path(out_dir))?;

    let mut len_reader = open_reader(gene_len_file)?;
    let mut row = String::new();
    while len_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() || raw.starts_with("gene_id") {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 3 {
            continue;
        }
        gene_len.insert(cols[0].to_string(), cols[2].parse::<f64>().unwrap_or(0.0));
    }

    let mut pn_reader = open_reader(parent_file)?;
    let mut p_row = String::new();
    while pn_reader.read_line(&mut p_row)? > 0 {
        let raw = p_row.trim_end_matches(['\n', '\r']).to_string();
        p_row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 3 {
            continue;
        }
        parent_map.insert(cols[0].to_string(), cols[2].to_string());
    }

    let mut ex_reader = open_reader(expr_file)?;
    let _ = ex_reader.read_line(&mut expr_row)?;
    while ex_reader.read_line(&mut expr_row)? > 0 {
        let raw = expr_row.trim_end_matches(['\n', '\r']).to_string();
        expr_row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.is_empty() {
            continue;
        }
        let gene = cols[0];
        if gene.starts_with("evm") || gene.starts_with("gene") || gene.starts_with("DMRT") {
            continue;
        }
        expr_map.insert(gene.to_string(), cols[1..].iter().map(|v| v.to_string()).collect());
    }

    let mut part: HashMap<String, Vec<String>> = HashMap::new();
    let mut full: HashMap<String, Vec<String>> = HashMap::new();

    for (gene, expr_fields) in expr_map {
        let family = gene.split('.').next().unwrap_or("");
        let Some(parent_gene) = parent_map.get(family) else {
            continue;
        };
        let Some(g_len) = gene_len.get(&gene) else {
            continue;
        };
        let Some(p_len) = gene_len.get(parent_gene) else {
            continue;
        };
        if p_len == &0.0 {
            continue;
        }

        let line = format!("{}\t{}", gene, expr_fields.join("\t"));
        if g_len / p_len >= coverage {
            full.entry(family.to_string())
                .or_default()
                .push(line);
        } else {
            part.entry(family.to_string()).or_default().push(line);
        }
    }

    for (family, lines) in part {
        let out_path = Path::new(out_dir).join(format!("{family}.partial.geneExpression.txt"));
        let mut out = BufWriter::new(File::create(out_path)?);
        for l in lines {
            writeln!(out, "{l}")?;
        }
    }

    for (family, lines) in full {
        let out_path = Path::new(out_dir).join(format!("{family}.fullenth.geneExpression.txt"));
        let mut out = BufWriter::new(File::create(out_path)?);
        for l in lines {
            writeln!(out, "{l}")?;
        }
    }

    Ok(0)
}

fn run_extract_gene_family_matrix(
    gene_len_file: &str,
    expr_file: &str,
    family_name_file: &str,
    gene_family_file: &str,
    coverage: f64,
    output: Option<&str>,
) -> Result<i32> {
    let mut gene_len: HashMap<String, f64> = HashMap::new();
    let mut family_name_set: HashSet<String> = HashSet::new();
    let mut gene_to_family: HashMap<String, String> = HashMap::new();
    let mut expr_avg: HashMap<String, f64> = HashMap::new();

    let mut r1 = String::new();

    let mut len_reader = open_reader(gene_len_file)?;
    let mut row = String::new();
    while len_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() || raw.starts_with("gene_id") {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        gene_len.insert(cols[0].to_string(), cols[1].parse::<f64>().unwrap_or(0.0));
    }

    let mut fn_reader = open_reader(family_name_file)?;
    let mut fn_row = String::new();
    while fn_reader.read_line(&mut fn_row)? > 0 {
        let raw = fn_row.trim_end_matches(['\n', '\r']).to_string();
        fn_row.clear();
        if raw.is_empty() {
            continue;
        }
        family_name_set.insert(raw);
    }

    let mut gf_reader = open_reader(gene_family_file)?;
    let mut gf_row = String::new();
    while gf_reader.read_line(&mut gf_row)? > 0 {
        let raw = gf_row.trim_end_matches(['\n', '\r']).to_string();
        gf_row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        gene_to_family.insert(cols[1].to_string(), cols[0].to_string());
    }

    let mut expr_reader = open_reader(expr_file)?;
    let _ = expr_reader.read_line(&mut r1)?;
    while expr_reader.read_line(&mut r1)? > 0 {
        let raw = r1.trim_end_matches(['\n', '\r']).to_string();
        r1.clear();
        if raw.is_empty() {
            continue;
        }
        let mut it = raw.split_whitespace();
        let Some(gene) = it.next() else { continue };
        let nums: Vec<f64> = it.filter_map(|v| v.parse::<f64>().ok()).collect();
        if nums.is_empty() {
            continue;
        }
        let avg = nums.iter().sum::<f64>() / nums.len() as f64;
        expr_avg.insert(gene.to_string(), avg);
    }

    let mut family_genes: HashMap<String, Vec<String>> = HashMap::new();
    for (gene, family) in gene_to_family {
        if family_name_set.contains(&family) {
            family_genes.entry(family).or_default().push(gene);
        }
    }

    let mut out = open_writer(output.or(Some("final.matrix.txt")))?;
    for (family, genes) in family_genes {
        if genes.is_empty() {
            continue;
        }
        let mut available: Vec<(String, f64)> = Vec::new();
        for g in &genes {
            if let Some(l) = gene_len.get(g) {
                available.push((g.clone(), *l));
            }
        }
        if available.is_empty() {
            continue;
        }
        let mut sorted_idx = 0usize;
        let mut max_len = available[0].1;
        for (i, (_, l)) in available.iter().enumerate() {
            if *l > max_len {
                max_len = *l;
                sorted_idx = i;
            }
        }
        let parent = available[sorted_idx].0.clone();
        let parent_len = gene_len.get(&parent).copied().unwrap_or(0.0);

        let mut partial = Vec::new();
        let mut full = Vec::new();
        full.push(parent.clone());

        for g in genes {
            if g == parent {
                continue;
            }
            let l = gene_len.get(&g).copied().unwrap_or(0.0);
            if parent_len > 0.0 && l / parent_len >= coverage {
                full.push(g);
            } else {
                partial.push(g);
            }
        }

        let p_avg = if partial.is_empty() {
            0.0
        } else {
            partial.iter().map(|g| expr_avg.get(g).copied().unwrap_or(0.0)).sum::<f64>() / partial.len() as f64
        };
        let f_avg = full.iter().map(|g| expr_avg.get(g).copied().unwrap_or(0.0)).sum::<f64>() / full.len() as f64;
        writeln!(out, "{family}\t{p_avg}\t{f_avg}")?;
    }

    Ok(0)
}

fn normalize_bam_name(name: &str, trim_suffix: usize) -> String {
    if trim_suffix == 0 || trim_suffix >= name.len() {
        name.to_string()
    } else {
        name[..name.len() - trim_suffix].to_string()
    }
}

fn read_sam_header(path: &str) -> Result<Vec<String>> {
    let mut cmd = Command::new("samtools")
        .args(["view", "-H", path])
        .stdout(Stdio::piped())
        .spawn()?;
    let mut out = Vec::new();
    let stdout = cmd.stdout.take().ok_or("missing samtools stdout")?;
    let mut r = BufReader::new(stdout);
    let mut line = String::new();
    while r.read_line(&mut line)? > 0 {
        out.push(line.trim_end_matches(['\n', '\r']).to_string());
        line.clear();
    }
    let status = cmd.wait()?;
    if !status.success() {
        return Err("samtools view -H failed".into());
    }
    Ok(out)
}

fn read_sam_records(path: &str, trim_suffix: usize) -> Result<HashMap<String, String>> {
    let mut cmd = Command::new("samtools")
        .arg("view")
        .arg(path)
        .stdout(Stdio::piped())
        .spawn()?;
    let stdout = cmd.stdout.take().ok_or("missing samtools stdout")?;
    let mut r = BufReader::new(stdout);
    let mut line = String::new();
    let mut data = HashMap::new();

    while r.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let mut it = raw.split('\t');
        let name = it.next().unwrap_or("");
        data.insert(normalize_bam_name(name, trim_suffix), raw);
    }

    let status = cmd.wait()?;
    if !status.success() {
        return Err("samtools view failed".into());
    }
    Ok(data)
}

fn run_merge_two_end_bam_internal(
    r1: &str,
    r2: &str,
    out1: &str,
    out2: &str,
    trim_suffix: usize,
) -> Result<i32> {
    let header = read_sam_header(r1)?;
    let r2_map = read_sam_records(r2, trim_suffix)?;

    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|_| "clock error")?
        .as_nanos();
    let temp_dir = std::env::temp_dir();
    let r1_tmp = temp_dir.join(format!("biohub_merge_{}_r1.sam", nanos));
    let r2_tmp = temp_dir.join(format!("biohub_merge_{}_r2.sam", nanos));
    let mut out1_s = BufWriter::new(File::create(&r1_tmp)?);
    let mut out2_s = BufWriter::new(File::create(&r2_tmp)?);

    for h in &header {
        writeln!(out1_s, "{h}")?;
        writeln!(out2_s, "{h}")?;
    }

    let mut cmd = Command::new("samtools")
        .arg("view")
        .arg(r1)
        .stdout(Stdio::piped())
        .spawn()?;
    let stdout = cmd.stdout.take().ok_or("missing r1 stdout")?;
    let mut r1_reader = BufReader::new(stdout);
    let mut line = String::new();
    while r1_reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let key = normalize_bam_name(raw.split('\t').next().unwrap_or(""), trim_suffix);
        if let Some(r2_line) = r2_map.get(&key) {
            writeln!(out1_s, "{raw}")?;
            writeln!(out2_s, "{r2_line}")?;
        }
    }
    let status = cmd.wait()?;
    if !status.success() {
        return Err("samtools view failed".into());
    }

    let st1 = Command::new("samtools")
        .args(["view", "-bS", r1_tmp.to_string_lossy().as_ref(), "-o", out1])
        .status()?;
    if !st1.success() {
        return Err("samtools convert r1 failed".into());
    }
    let st2 = Command::new("samtools")
        .args(["view", "-bS", r2_tmp.to_string_lossy().as_ref(), "-o", out2])
        .status()?;
    if !st2.success() {
        return Err("samtools convert r2 failed".into());
    }

    let _ = fs::remove_file(&r1_tmp);
    let _ = fs::remove_file(&r2_tmp);
    Ok(0)
}

fn run_merge_two_end_bam(r1: &str, r2: &str, out1: Option<&str>, out2: Option<&str>) -> Result<i32> {
    let out_r1 = out1.unwrap_or("test13.h1.R1.outReads.bam");
    let out_r2 = out2.unwrap_or("test13.h1.R2.outReads.bam");
    run_merge_two_end_bam_internal(r1, r2, out_r1, out_r2, 0)
}

fn run_merge_two_end_bam1(r1: &str, r2: &str, out1: Option<&str>, out2: Option<&str>) -> Result<i32> {
    let out_r1 = out1.unwrap_or("R1.outReads.bam");
    let out_r2 = out2.unwrap_or("R2.outReads.bam");
    run_merge_two_end_bam_internal(r1, r2, out_r1, out_r2, 0)
}

fn run_merge_two_end_bam_for_mgi(r1: &str, r2: &str, out1: Option<&str>, out2: Option<&str>) -> Result<i32> {
    let out_r1 = out1.unwrap_or("R1.outReads.bam");
    let out_r2 = out2.unwrap_or("R2.outReads.bam");
    run_merge_two_end_bam_internal(r1, r2, out_r1, out_r2, 2)
}

fn run_zhouxiaoxuan_merge_xls(
    first_file: &str,
    second_file: &str,
    output: Option<&str>,
) -> Result<i32> {
    let mut map: HashMap<String, (String, Vec<String>)> = HashMap::new();
    let mut base_name_map: HashMap<String, String> = HashMap::new();

    let mut first_reader = open_reader(first_file)?;
    let mut line = String::new();
    while first_reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.is_empty() {
            continue;
        }
        let gene = cols[0].to_string();
        let extras: Vec<String> = cols.iter().skip(1).map(|v| v.to_string()).collect();
        map.insert(gene.clone(), (gene.clone(), extras));
        if let Some(base_prefix) = gene.split('.').next() {
            if let Some((_lhs, rhs)) = base_prefix.split_once(':') {
                base_name_map.insert(rhs.to_string(), gene.clone());
            }
        }
    }

    let mut out = open_writer(output)?;
    let mut second_reader = open_reader(second_file)?;
    let mut unmatched: Vec<String> = Vec::new();
    let mut row = String::new();
    while second_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let key = raw.split('.').next().unwrap_or("");
        if let Some(name) = base_name_map.get(key) {
            if let Some((_, extras)) = map.get(name) {
                writeln!(out, "{}\t{}\t{}", raw, name, extras.join("\t"))?;
            }
        } else {
            unmatched.push(raw);
        }
    }
    for l in unmatched {
        writeln!(out, "{l}")?;
    }
    Ok(0)
}

fn run_check_duplication_gene_pairs(input: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(input)?;
    let mut out = open_writer(output)?;
    let mut seen: HashSet<String> = HashSet::new();
    let mut row = String::new();

    while reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 2 {
            continue;
        }
        let a = cols[0];
        let b = cols[1];
        let k1 = format!("{a}\t{b}");
        let k2 = format!("{b}\t{a}");
        if !seen.contains(&k1) && !seen.contains(&k2) {
            seen.insert(k1.clone());
            writeln!(out, "{k1}")?;
        }
    }
    Ok(0)
}

fn read_fasta_records(path: &str) -> Result<Vec<(String, String)>> {
    let mut reader = open_reader(path)?;
    let mut line = String::new();
    let mut records: Vec<(String, String)> = Vec::new();
    let mut name: Option<String> = None;
    let mut seq = String::new();

    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        if let Some(rest) = raw.strip_prefix('>') {
            if let Some(h) = name.take() {
                records.push((h, seq.clone()));
                seq.clear();
            }
            name = Some(rest.trim().to_string());
        } else {
            seq.push_str(raw.trim());
        }
    }
    if let Some(h) = name {
        records.push((h, seq));
    }

    Ok(records)
}

fn reverse_complement(seq: &str) -> String {
    let mut out = String::with_capacity(seq.len());
    for b in seq.as_bytes().iter().rev() {
        out.push(match b {
            b'A' | b'a' => 'T',
            b'T' | b't' => 'A',
            b'G' | b'g' => 'C',
            b'C' | b'c' => 'G',
            b'N' | b'n' => 'N',
            _ => 'N',
        });
    }
    out
}

fn generate_rotations(seq: &str) -> Vec<String> {
    let mut out = Vec::new();
    let seq = seq.as_bytes();
    if seq.is_empty() {
        return out;
    }
    for i in 0..seq.len() {
        let mut rotated = String::with_capacity(seq.len());
        for j in 0..seq.len() {
            rotated.push(seq[(i + j) % seq.len()] as char);
        }
        out.push(rotated);
    }
    out
}

fn four_fold_codon_set() -> HashSet<String> {
    let aa_codons: Vec<Vec<&str>> = vec![
        vec!["CTT", "CTA", "CTC", "CTG"],
        vec!["GTT", "GTA", "GTC", "GTG"],
        vec!["TCT", "TCA", "TCC", "TCG"],
        vec!["CCT", "CCA", "CCC", "CCG"],
        vec!["ACT", "ACA", "ACC", "ACG"],
        vec!["GCT", "GCA", "GCC", "GCG"],
        vec!["CGT", "CGA", "CGC", "CGG"],
        vec!["GGT", "GGA", "GGC", "GGG"],
        vec!["GTT", "GTA", "GTC", "GTG"],
    ];
    let mut out: HashSet<String> = HashSet::new();
    for group in aa_codons {
        for c in group {
            out.insert(c.to_string());
        }
    }
    out
}

fn run_get_longest_transcript(input: &str, output: Option<&str>) -> Result<i32> {
    let records = read_fasta_records(input)?;
    let mut seq_by_header: HashMap<String, String> = HashMap::new();
    let mut by_gene: HashMap<String, (String, usize)> = HashMap::new();
    for (header, seq) in &records {
        seq_by_header.insert(header.clone(), seq.clone());

        if header.is_empty() {
            continue;
        }
        let gene = header
            .split('.')
            .collect::<Vec<_>>()
            .split_last()
            .map(|(_, rest)| rest.join("."))
            .unwrap_or_else(String::new);
        if gene.is_empty() {
            continue;
        }
        let len = seq.len();
        match by_gene.get_mut(&gene) {
            Some((best_header, best_len)) => {
                if len > *best_len {
                    *best_len = len;
                    *best_header = header.to_string();
                }
            }
            None => {
                by_gene.insert(gene, (header.clone(), len));
            }
        }
    }

    let mut out = open_writer(output)?;
    for (_, (header, _)) in by_gene {
        if let Some(seq) = seq_by_header.get(&header) {
            writeln!(out, ">{header}")?;
            writeln!(out, "{seq}")?;
        }
    }
    Ok(0)
}

fn run_orthofiner_to_pal2nal(path_of_prot: &str, out_dir: &str, nucl_cds: &str) -> Result<i32> {
    fs::create_dir_all(out_dir)?;

    let nucl_records = read_fasta_records(nucl_cds)?;
    let mut nucl_map: HashMap<String, String> = HashMap::new();
    for (h, s) in nucl_records {
        nucl_map.insert(h, s);
    }

    let mut files: Vec<PathBuf> = Vec::new();
    for e in fs::read_dir(path_of_prot)? {
        let e = e?;
        let p = e.path();
        if p.is_file() {
            if let Some(name) = p.file_name().and_then(|x| x.to_str()) {
                if name.starts_with('.') {
                    continue;
                }
            }
            files.push(p);
        }
    }
    files.sort();

    for prot_path in files {
        let recs = read_fasta_records(prot_path.to_string_lossy().as_ref())?;
        if recs.is_empty() {
            continue;
        }

        let stem = prot_path
            .file_name()
            .and_then(|s| s.to_str())
            .ok_or("invalid protein file name")?;
        let out_aln = Path::new(out_dir).join(stem);
        let out_nucl = Path::new(out_dir).join(format!("{stem}.cds.fasta"));
        let out_nucl_aln = out_nucl.with_extension("aln.fasta");

        let mut nucl_out = BufWriter::new(File::create(&out_nucl)?);
        for (h, _) in &recs {
            let key = if h.starts_with("Migut") {
                h.split('.').take(3).collect::<Vec<_>>().join(".")
            } else {
                h.to_string()
            };
            if let Some(seq) = nucl_map.get(&key) {
                writeln!(nucl_out, ">{key}")?;
                writeln!(nucl_out, "{seq}")?;
            }
        }
        drop(nucl_out);

        let mafft_status = Command::new("mafft")
            .args(["--maxiterate", "1000", "--localpair", &prot_path.to_string_lossy()])
            .stdout(Stdio::from(File::create(&out_aln)?))
            .status()?;
        if !mafft_status.success() {
            return Err(format!("mafft failed for {}", prot_path.display()).into());
        }

        let pal2nal_status = Command::new("pal2nal.pl")
            .args([
                out_aln.to_string_lossy().as_ref(),
                out_nucl.to_string_lossy().as_ref(),
                "-output",
                "paml",
                "-o",
                out_nucl_aln.to_string_lossy().as_ref(),
            ])
            .status()?;
        if !pal2nal_status.success() {
            return Err(format!("pal2nal failed for {}", prot_path.display()).into());
        }
    }

    Ok(0)
}

fn run_get_diff_sites_from_orthology(pal2nal_dir: &str, out_dir: &str) -> Result<i32> {
    fs::create_dir_all(out_dir)?;
    let four_fold = four_fold_codon_set();

    let mut sample_seq: HashMap<String, String> = HashMap::new();
    let mut sample_set: Vec<String> = Vec::new();

    let mut files: Vec<PathBuf> = Vec::new();
    for e in fs::read_dir(pal2nal_dir)? {
        let e = e?;
        if e.path().is_file() {
            files.push(e.path());
        }
    }
    files.sort();

    let mut lengths: Vec<usize> = Vec::new();
    for p in &files {
        for (h, seq) in read_fasta_records(&p.to_string_lossy())? {
            let sample = h.chars().take(9).collect::<String>();
            let rec = sample_seq.entry(sample.clone()).or_insert_with(String::new);
            if rec.is_empty() {
                sample_set.push(sample);
            }
            rec.push_str(&seq);
            lengths.push(seq.len());
        }
    }

    if sample_seq.is_empty() {
        return Ok(0);
    }
    let step_len = if lengths.is_empty() { 0 } else { lengths[0] };

    let mut first = HashMap::new();
    let mut second = HashMap::new();
    let mut third = HashMap::new();
    let mut four_site = HashMap::new();
    let mut four_code = HashMap::new();

    for name in sample_seq.keys() {
        first.insert(name.clone(), String::new());
        second.insert(name.clone(), String::new());
        third.insert(name.clone(), String::new());
        four_site.insert(name.clone(), String::new());
        four_code.insert(name.clone(), String::new());
    }

    let mut four_pos: Vec<usize> = Vec::new();
    let mut i = 0usize;
    while i + 3 <= step_len {
        let mut hit = 0usize;
        for seq in sample_seq.values() {
            if i + 3 <= seq.len() {
                let codon = &seq[i..i + 3];
                if four_fold.contains(&codon.to_string()) {
                    hit += 1;
                }
            }
        }
        let threshold = sample_seq.len();
        if hit >= threshold {
            four_pos.push(i + 2);
            for (name, seq) in &sample_seq {
                if let Some(f) = first.get_mut(name) {
                    f.push(seq.chars().nth(i).unwrap_or('N'));
                }
                if let Some(s) = second.get_mut(name) {
                    s.push(seq.chars().nth(i + 1).unwrap_or('N'));
                }
                if let Some(t) = third.get_mut(name) {
                    t.push(seq.chars().nth(i + 2).unwrap_or('N'));
                }
                if let Some(fs) = four_site.get_mut(name) {
                    fs.push(seq.chars().nth(i + 2).unwrap_or('N'));
                }
                if let Some(fc) = four_code.get_mut(name) {
                    fc.push_str(&seq[i..i + 3.min(seq.len() - i)]);
                }
            }
        }
        i += 3;
    }

    let mut all_site = File::create(Path::new(out_dir).join("allSingleCopyGeneSeq.fasta"))?;
    for name in &sample_set {
        if let Some(seq) = sample_seq.get(name) {
            writeln!(all_site, ">{name}")?;
            writeln!(all_site, "{seq}")?;
        }
    }

    let write_sites = |path: &str, map: &HashMap<String, String>| -> Result<()> {
        let mut out = BufWriter::new(File::create(Path::new(out_dir).join(path))?);
        for name in &sample_set {
            if let Some(seq) = map.get(name) {
                writeln!(out, ">{name}")?;
                writeln!(out, "{seq}")?;
            }
        }
        Ok(())
    };

    write_sites("TheFirstSite.fasta", &first)?;
    write_sites("TheSecondSite.fasta", &second)?;
    write_sites("TheThirdSite.fasta", &third)?;
    write_sites("fourDegenerateSite.fasta", &four_site)?;
    write_sites("fourDegenerateCodenFile.fasta", &four_code)?;

    let mut fp = File::create(Path::new(out_dir).join("filesImportScripts.txt"))?;
    for p in files {
        writeln!(fp, "{}", p.to_string_lossy())?;
    }

    eprintln!("fourfold sites found: {}", four_pos.len());
    Ok(0)
}

fn run_get_four_degenerate_sites(pal2nal_dir: &str, out_dir: &str) -> Result<i32> {
    fs::create_dir_all(out_dir)?;
    let four_fold = four_fold_codon_set();

    let mut files: Vec<PathBuf> = Vec::new();
    for e in fs::read_dir(pal2nal_dir)? {
        let e = e?;
        if e.path().is_file() {
            files.push(e.path());
        }
    }
    files.sort();

    let mut sample_site: HashMap<String, String> = HashMap::new();
    let mut total_site_count = 0usize;
    let mut total_len = 0usize;

    for p in files {
        let recs = read_fasta_records(&p.to_string_lossy())?;
        if recs.is_empty() {
            continue;
        }
        if let Some((_, ref_seq)) = recs.first() {
            total_len += ref_seq.len() / 3;
        }

        let mut all_sites = Vec::new();
        let align_len = recs.iter().map(|(_, s)| s.len()).min().unwrap_or(0);
        let mut i = 0usize;
        while i + 3 <= align_len {
            let codon = recs[0].1[i..i + 3].to_string();
            if four_fold.contains(&codon) {
                let mut unique: HashSet<String> = HashSet::new();
                for (_, s) in &recs {
                    unique.insert(s[i..i + 3].to_string());
                }
                if unique.len() == 1 {
                    all_sites.push(i + 2);
                    total_site_count += 1;
                }
            }
            i += 3;
        }

        let mut sample_records: HashMap<String, String> = HashMap::new();
        for (h, s) in &recs {
            let sample = h.chars().take(9).collect::<String>();
            sample_records.entry(sample.clone()).or_insert_with(|| s.clone());
            sample_site.entry(sample).or_insert_with(String::new);
        }

        for (sample, entries) in sample_records {
            let seq = entries.as_str();
            let acc = sample_site.get_mut(&sample).unwrap();
            for pos in &all_sites {
                if *pos < seq.len() {
                    acc.push(seq.chars().nth(*pos).unwrap_or('N'));
                }
            }
        }
    }

    let mut out = BufWriter::new(File::create(Path::new(out_dir).join("fourDegenerateSite.fasta"))?);
    for (sample, seq) in sample_site {
        writeln!(out, ">{sample}")?;
        writeln!(out, "{seq}")?;
    }

    let mut stat = BufWriter::new(File::create(Path::new(out_dir).join("fourDegenerateSite.stat"))?);
    writeln!(stat, "#totalSites\t{}", total_site_count)?;
    writeln!(stat, "#allSites\t{}", total_len)?;

    eprintln!("fourdegenerate sites collected: {}", total_site_count);
    Ok(0)
}

fn run_plot_depth_pandepth(input: &str, output: &str, min_chr_length_mb: f64) -> Result<i32> {
    run_plot_depth_pandepth_impl(input, output, min_chr_length_mb, false)
}

fn run_plot_depth_pandepth_impl(
    input: &str,
    output: &str,
    min_chr_length_mb: f64,
    styled: bool,
) -> Result<i32> {
    let min_len_bp = (min_chr_length_mb * 1_000_000.0) as usize;
    fs::create_dir_all(output)?;
    let mut reader = open_reader(input)?;
    let mut line = String::new();
    let mut header: Vec<String> = Vec::new();
    let mut chr_idx = 0usize;
    let mut start_idx = 1usize;
    let mut end_idx = 2usize;
    let mut depth_idx = 3usize;
    let mut gc_idx = 4usize;

    let mut rows: Vec<(String, usize, usize, f64, f64)> = Vec::new();
    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        if raw.starts_with('#') {
            if header.is_empty() {
                let cols = raw.trim_start_matches('#').trim();
                if !cols.is_empty() && !cols.starts_with('#') {
                    header = cols.split('\t').map(|v| v.to_string()).collect();
                    chr_idx = header.iter().position(|v| v == "Chr").unwrap_or(0);
                    start_idx = header.iter().position(|v| v == "Start").unwrap_or(1);
                    end_idx = header.iter().position(|v| v == "End").unwrap_or(2);
                    depth_idx = header.iter().position(|v| v == "MeanDepth").unwrap_or(3);
                    gc_idx = header.iter().position(|v| v == "GC(%)").unwrap_or(4);
                }
            }
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 5 || cols.is_empty() {
            continue;
        }
        let chr = cols.get(chr_idx).unwrap_or(&cols[0]).to_string();
        let start: usize = cols
            .get(start_idx)
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(0);
        let end: usize = cols
            .get(end_idx)
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(0);
        let depth: f64 = cols
            .get(depth_idx)
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.0);
        let gc: f64 = cols
            .get(gc_idx)
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.0);
        if start > end {
            continue;
        }
        rows.push((chr, start, end, depth, gc));
    }

    if rows.is_empty() {
        return Err("no valid rows found for plot-depth-pandepth".into());
    }

    let mut len_map: HashMap<String, usize> = HashMap::new();
    let mut depth_sum: HashMap<String, (f64, usize)> = HashMap::new();
    let mut gc_sum: HashMap<String, (f64, usize)> = HashMap::new();
    for (chr, _start, end, depth, gc) in &rows {
        len_map
            .entry(chr.clone())
            .and_modify(|v| {
                if *end > *v {
                    *v = *end;
                }
            })
            .or_insert(*end);
        depth_sum
            .entry(chr.clone())
            .and_modify(|(s, n)| {
                *s += *depth;
                *n += 1;
            })
            .or_insert((*depth, 1));
        gc_sum
            .entry(chr.clone())
            .and_modify(|(s, n)| {
                *s += *gc;
                *n += 1;
            })
            .or_insert((*gc, 1));
    }

    let mut filtered: Vec<(String, usize, usize, f64, f64)> = Vec::new();
    let mut plot_points: Vec<ScatterPoint> = Vec::new();
        for row in rows {
        let chr = row.0.clone();
        let len = *len_map.get(&chr).unwrap_or(&0);
        if len >= min_len_bp {
            filtered.push(row.clone());
            plot_points.push(ScatterPoint {
                x: row.4,
                y: row.3,
                group: chr,
            });
        }
    }

    if filtered.is_empty() {
        return Err(format!(
            "no chromosome longer than {min_len_bp} bp in input"
        )
        .into());
    }

    let mut stats_rows: Vec<String> = len_map.keys().cloned().collect();
    stats_rows.sort_unstable();
    let mut out_stats = BufWriter::new(File::create(Path::new(output).join("chromosome_stats.tsv"))?);
    writeln!(out_stats, "Chr\tlength\tmean_depth\tmean_gc")?;
    for chr in stats_rows {
        let len = *len_map.get(&chr).unwrap_or(&0);
        let (dsum, dcnt) = depth_sum.get(&chr).copied().unwrap_or((0.0, 0));
        let (gsum, gcnt) = gc_sum.get(&chr).copied().unwrap_or((0.0, 0));
        let dmean = if dcnt == 0 { 0.0 } else { dsum / dcnt as f64 };
        let gmean = if gcnt == 0 { 0.0 } else { gsum / gcnt as f64 };
        writeln!(out_stats, "{chr}\t{len}\t{dmean}\t{gmean}")?;
    }

    let mut out_data = BufWriter::new(File::create(Path::new(output).join("filtered_depth.tsv"))?);
    writeln!(out_data, "Chr\tStart\tEnd\tMeanDepth\tGC(%)")?;
    for (chr, start, end, depth, gc) in filtered {
        writeln!(out_data, "{chr}\t{start}\t{end}\t{depth}\t{gc}")?;
    }

    let svg_points = downsample(&plot_points, if styled { 40000 } else { 25000 });
    let svg_file = if styled {
        Path::new(output).join("depth_gc_styled.svg")
    } else {
        Path::new(output).join("depth_gc_scatter.svg")
    };
    write_scatter_svg(
        &svg_file,
        if styled {
            "Pandepth depth vs GC (styled)"
        } else {
            "Pandepth depth vs GC"
        },
        "GC(%)",
        "MeanDepth",
        &svg_points,
        if styled { 1.0 } else { 1.6 },
        1400,
        900,
        None,
    )?;

    Ok(0)
}

fn run_plot_depth_pandepth2(input: &str, output: &str, min_chr_length_mb: f64) -> Result<i32> {
    run_plot_depth_pandepth_impl(input, output, min_chr_length_mb, true)
}

fn run_plot_mosdepth_point(input: &str, output: &str, min_length: usize) -> Result<i32> {
    fs::create_dir_all(output)?;
    let mut reader = open_reader(input)?;
    let mut line = String::new();
    let mut data: Vec<(String, usize, usize, f64)> = Vec::new();
    let mut len_map: HashMap<String, usize> = HashMap::new();

    while reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 4 {
            continue;
        }
        let chrom = cols[0].to_string();
        let start = cols[1].parse::<usize>().unwrap_or(0);
        let end = cols[2].parse::<usize>().unwrap_or(0);
        let cov = cols[3].parse::<f64>().unwrap_or(0.0);
        data.push((chrom.clone(), start, end, cov));
        len_map.entry(chrom).and_modify(|x| if end > *x { *x = end }).or_insert(end);
    }

    if len_map.values().all(|l| *l < min_length) {
        return Err(format!("No chromosome longer than {min_length} bp").into());
    }

    let mut out = BufWriter::new(File::create(Path::new(output).join("mosdepth_points.tsv"))?);
    writeln!(out, "chrom\tstart\tend\tcoverage\tstatus")?;
    for (c, s, e, cov) in data.iter() {
        if let Some(len) = len_map.get(c.as_str()) {
            if *len >= min_length {
                writeln!(out, "{c}\t{s}\t{e}\t{cov}\tOK")?;
            }
        }
    }

    let mut chr_rows: Vec<(String, usize)> = len_map.iter().map(|(c, l)| (c.clone(), *l)).collect();
    chr_rows.sort_by(|a, b| a.0.cmp(&b.0));
    chr_rows = chr_rows
        .into_iter()
        .filter(|(_, l)| *l >= min_length)
        .collect();
    if chr_rows.is_empty() {
        return Err(format!("No chromosome longer than {min_length} bp").into());
    }

    let mut chr_order: Vec<String> = Vec::new();
    let mut chr_offsets: HashMap<String, f64> = HashMap::new();
    let mut cursor: f64 = 0.0;
    for (chr, len) in &chr_rows {
        chr_order.push(chr.clone());
        chr_offsets.insert(chr.clone(), cursor);
        cursor += *len as f64 + 1_000_000.0_f64.max(*len as f64 * 0.01);
    }

    let mut points: Vec<ScatterPoint> = Vec::new();
    let mut custom_ticks: Vec<(f64, String)> = Vec::new();
    for (c, start, end, cov) in &data {
        if let Some(len) = len_map.get(c) {
            if *len < min_length {
                continue;
            }
            let offset = chr_offsets.get(c).cloned().unwrap_or(0.0);
            let center = (*start + *end) as f64 / 2.0 + offset;
            points.push(ScatterPoint {
                x: center,
                y: *cov,
                group: c.clone(),
            });
        }
    }
    if points.is_empty() {
        return Err("No data points after filtering by chromosome length".into());
    }

    let mut cumulative = 0.0;
    for (idx, chr) in chr_order.iter().enumerate() {
        if idx == 0 {
            cumulative = *chr_offsets.get(chr).unwrap_or(&0.0);
        }
        if let Some(len) = len_map.get(chr) {
            let center = cumulative + (*len as f64) / 2.0;
            custom_ticks.push((center, chr.clone()));
            cumulative += (*len as f64) + 1_000_000.0_f64.max(*len as f64 * 0.01);
        }
    }

    let sampled = downsample(&points, 60000);
    let svg_path = Path::new(output).join("mosdepth_scatter.svg");
    write_scatter_svg(
        &svg_path,
        "Mosdepth point coverage",
        "Cumulative coordinate (bp)",
        "Coverage",
        &sampled,
        1.2,
        1600,
        900,
        Some(&custom_ticks),
    )?;
    Ok(0)
}

fn run_trim_ttaggg_fastq(input: &str, output: &str, motif: &str) -> Result<i32> {
    if input.ends_with(".gz") || output.ends_with(".gz") {
        return Err("gzip FASTQ is not supported in this Rust build".into());
    }
    let mut reader = open_reader(input)?;
    let mut writer = open_writer(Some(output))?;

    let rc = reverse_complement(motif);
    let mut keep: HashSet<String> = HashSet::new();
    for m in generate_rotations(motif) {
        keep.insert(m);
    }
    for m in generate_rotations(&rc) {
        keep.insert(m);
    }

    loop {
        let mut h = String::new();
        let mut s = String::new();
        let mut p = String::new();
        let mut q = String::new();

        if reader.read_line(&mut h)? == 0 {
            break;
        }
        if reader.read_line(&mut s)? == 0 || reader.read_line(&mut p)? == 0 || reader.read_line(&mut q)? == 0 {
            break;
        }
        let seq = s.trim_end_matches(['\n', '\r']).to_uppercase();
        let head = seq.chars().take(6).collect::<String>();
        let tail = if seq.len() > 6 {
            seq.chars().rev().take(6).collect::<String>().chars().rev().collect()
        } else {
            String::new()
        };
        if !keep.contains(&head) && !keep.contains(&tail) {
            writer.write_all(h.as_bytes())?;
            writer.write_all(s.as_bytes())?;
            writer.write_all(p.as_bytes())?;
            writer.write_all(q.as_bytes())?;
        }

        h.clear();
        s.clear();
        p.clear();
        q.clear();
    }

    Ok(0)
}

#[derive(Clone)]
struct IntervalFeature {
    start: usize,
    end: usize,
    name: String,
    transcript: String,
}

fn parse_gff_attrs(attrs: &str) -> HashMap<String, String> {
    let mut out = HashMap::new();
    for raw in attrs.split(';') {
        let raw = raw.trim();
        if raw.is_empty() {
            continue;
        }
        if let Some(eq_pos) = raw.find('=') {
            let key = raw[..eq_pos].trim().to_string();
            let mut value = raw[eq_pos + 1..].trim().to_string();
            if value.starts_with('"') && value.ends_with('"') && value.len() >= 2 {
                value = value[1..value.len() - 1].to_string();
            }
            if !key.is_empty() {
                out.insert(key, value);
            }
            continue;
        }

        if let Some(space_pos) = raw.find(' ') {
            let key = raw[..space_pos].trim().to_string();
            let value = raw[space_pos + 1..].trim().trim_matches('"').to_string();
            if !key.is_empty() {
                out.insert(key, value);
            }
        }
    }
    out
}

fn choose_alt_first(alt: &str) -> &str {
    if let Some((first, _)) = alt.split_once(',') {
        first
    } else {
        alt
    }
}

fn var_kind(ref_allele: &str, alt_allele: &str) -> &'static str {
    if ref_allele.len() == 1 && alt_allele.len() == 1 {
        "SNV"
    } else if ref_allele.len() < alt_allele.len() {
        "INS"
    } else if ref_allele.len() > alt_allele.len() {
        "DEL"
    } else {
        "MNV"
    }
}

fn run_annotation_vcf(reference: &str, gff: &str, vcf: &str, out_path: &str, format: &str) -> Result<i32> {
    let mut fasta_reader = open_reader(reference)?;
    let mut _line = String::new();
    let mut reference_chromosomes: HashMap<String, usize> = HashMap::new();
    while fasta_reader.read_line(&mut _line)? > 0 {
        let raw = _line.trim_end_matches(['\n', '\r']).to_string();
        _line.clear();
        if raw.starts_with('>') {
            let header = raw[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            reference_chromosomes.insert(header, 0);
        }
    }

    let mut gene_intervals: HashMap<String, Vec<(usize, usize, String, char)>> = HashMap::new();
    let mut transcript_name_to_gene: HashMap<String, String> = HashMap::new();
    let mut transcript_intervals: HashMap<String, Vec<(String, usize, usize, String, char)>> = HashMap::new();
    let mut cds_intervals: HashMap<String, Vec<IntervalFeature>> = HashMap::new();
    let mut exon_intervals: HashMap<String, Vec<IntervalFeature>> = HashMap::new();

    let mut gff_reader = open_reader(gff)?;
    let mut line = String::new();
    while gff_reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let chr = cols[0];
        let feature = cols[2];
        let start = cols[3].parse::<usize>().unwrap_or(0);
        let end = cols[4].parse::<usize>().unwrap_or(0);
        let strand = cols[6].chars().next().unwrap_or('+');
        let attrs = parse_gff_attrs(cols[8]);
        let gene_id = attrs
            .get("gene_id")
            .cloned()
            .or_else(|| attrs.get("gene").cloned())
            .or_else(|| attrs.get("Name").cloned())
            .or_else(|| attrs.get("ID").cloned())
            .unwrap_or_else(|| "NA".to_string());
        let tx_id = attrs
            .get("transcript_id")
            .or_else(|| attrs.get("Parent"))
            .or_else(|| attrs.get("ID"))
            .map(|v| v.to_string())
            .unwrap_or_else(|| "NA".to_string());

        match feature {
            "gene" => {
                gene_intervals
                    .entry(chr.to_string())
                    .or_default()
                    .push((start, end, gene_id, strand));
            }
            "mRNA" | "transcript" => {
                transcript_name_to_gene.insert(tx_id.clone(), gene_id.clone());
                transcript_intervals
                    .entry(chr.to_string())
                    .or_default()
                    .push((tx_id.clone(), start, end, gene_id, strand));
            }
            "CDS" => {
                let parent = attrs.get("Parent").cloned().unwrap_or_else(|| tx_id.clone());
                let tx = transcript_name_to_gene.get(&parent).cloned().unwrap_or(parent.clone());
                cds_intervals.entry(chr.to_string()).or_default().push(IntervalFeature {
                    start,
                    end,
                    name: tx_id.clone(),
                    transcript: tx,
                });
            }
            "exon" => {
                let parent = attrs.get("Parent").cloned().unwrap_or_else(|| tx_id.clone());
                let tx = transcript_name_to_gene.get(&parent).cloned().unwrap_or(parent.clone());
                exon_intervals.entry(chr.to_string()).or_default().push(IntervalFeature {
                    start,
                    end,
                    name: tx_id.clone(),
                    transcript: tx,
                });
            }
            _ => {}
        }
    }

    for chr in exon_intervals.values_mut() {
        chr.sort_by_key(|x| x.start);
    }
    for chr in cds_intervals.values_mut() {
        chr.sort_by_key(|x| x.start);
    }
    for chr in transcript_intervals.values_mut() {
        chr.sort_by_key(|x| x.1);
    }
    for chr in gene_intervals.values_mut() {
        chr.sort_by_key(|x| x.0);
    }

    let format_lower = format.to_lowercase();
    let emit_json = format_lower == "json" || format_lower == "pickle" || out_path.ends_with(".json") || out_path.ends_with(".pickle");
    let out_tsv = !emit_json || format_lower == "txt" || format_lower == "tsv" || out_path.ends_with(".txt") || out_path.ends_with(".tsv");
    let mut out = BufWriter::new(File::create(out_path)?);
    let mut out_json = if emit_json {
        Some(BufWriter::new(File::create(format!(
            "{}.json",
            out_path.trim_end_matches(".pickle").trim_end_matches(".json")
        ))?))
    } else {
        None
    };

    let overlaps_gene = |intervals: &Vec<(usize, usize, String, char)>, pos: usize| -> Vec<String> {
        let mut out = Vec::new();
        for (s, e, id, _) in intervals {
            if pos >= *s && pos <= *e {
                if !out.iter().any(|v| v == id) {
                    out.push(id.clone());
                }
            }
        }
        out
    };

    let overlaps_transcript = |intervals: &Vec<(String, usize, usize, String, char)>, pos: usize| -> Vec<(String, String)> {
        let mut out = Vec::new();
        for (tx, s, e, gene, _) in intervals {
            if pos >= *s && pos <= *e {
                out.push((tx.clone(), gene.clone()));
            }
        }
        out
    };

    let overlaps_interval = |intervals: &Vec<IntervalFeature>, pos: usize| -> Vec<(String, String)> {
        let mut out = Vec::new();
        for it in intervals {
            if pos >= it.start && pos <= it.end {
                out.push((it.name.clone(), it.transcript.clone()));
            }
        }
        out
    };

    if out_tsv {
        writeln!(out, "chrom\tpos\tref\talt\ttype\tstatus\tgene_ids\ttranscript_ids\tfeature")?;
    }
    if let Some(json_out) = out_json.as_mut() {
        writeln!(json_out, "[")?;
    }

    let mut vcf_reader = open_reader(vcf)?;
    let mut rline = String::new();
    let mut first_json = true;
    while vcf_reader.read_line(&mut rline)? > 0 {
        let raw = rline.trim_end_matches(['\n', '\r']).to_string();
        rline.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split('\t').collect();
        if cols.len() < 5 {
            continue;
        }
        let chr = cols[0].trim().to_string();
        if !reference_chromosomes.is_empty() && !reference_chromosomes.contains_key(&chr) {
            continue;
        }
        let pos = cols[1].parse::<usize>().unwrap_or(0);
        if pos == 0 {
            continue;
        }
        let reference = cols[3];
        let alt_raw = cols[4];
        let alt = choose_alt_first(alt_raw);
        let vtype = var_kind(reference, alt);

        let gene_hits = gene_intervals.get(&chr).map_or_else(Vec::new, |v| overlaps_gene(v, pos));
        let transcript_hits = transcript_intervals
            .get(&chr)
            .map_or_else(Vec::new, |v| overlaps_transcript(v, pos));
        let cds_hits = cds_intervals.get(&chr).map_or_else(Vec::new, |v| overlaps_interval(v, pos));
        let exon_hits = exon_intervals.get(&chr).map_or_else(Vec::new, |v| overlaps_interval(v, pos));

        let status = if !cds_hits.is_empty() {
            "CDS"
        } else if !exon_hits.is_empty() {
            "Exon"
        } else if !transcript_hits.is_empty() {
            "Intron"
        } else if !gene_hits.is_empty() {
            "Gene-body"
        } else {
            "Intergenic"
        };

        let mut tx_set = Vec::new();
        for (name, _) in transcript_hits.iter() {
            if !tx_set.contains(name) {
                tx_set.push(name.clone());
            }
        }
        for (name, _) in cds_hits.iter() {
            if !tx_set.contains(name) {
                tx_set.push(name.clone());
            }
        }
        for (name, _) in exon_hits.iter() {
            if !tx_set.contains(name) {
                tx_set.push(name.clone());
            }
        }
        let mut genes = gene_hits.clone();
        if genes.is_empty() {
            genes = transcript_hits
                .iter()
                .map(|(tx, _)| {
                    transcript_name_to_gene
                        .get(tx)
                        .cloned()
                        .unwrap_or_else(|| tx.clone())
                })
                .collect();
            if genes.is_empty() {
                genes = cds_hits
                    .iter()
                    .map(|(_, tx)| tx.clone())
                    .collect();
            }
            if genes.is_empty() {
                genes = exon_hits
                    .iter()
                    .map(|(_, tx)| tx.clone())
                    .collect();
            }
        }

        let gene_ids = if genes.is_empty() {
            "NA".to_string()
        } else {
            genes.join(";")
        };
        let tx_ids = if tx_set.is_empty() {
            "NA".to_string()
        } else {
            tx_set.join(";")
        };

        let feature = if !cds_hits.is_empty() {
            "CDS"
        } else if !exon_hits.is_empty() {
            "exon"
        } else if !transcript_hits.is_empty() {
            "transcript"
        } else if !gene_hits.is_empty() {
            "gene"
        } else {
            "intergenic"
        };

        if out_tsv {
            writeln!(
                out,
                "{chr}\t{pos}\t{reference}\t{alt}\t{vtype}\t{status}\t{gene_ids}\t{tx_ids}\t{feature}"
            )?;
        }

        if let Some(json_out) = out_json.as_mut() {
            if !first_json {
                writeln!(json_out, ",")?;
            }
            first_json = false;
            let gene_json = gene_ids.replace('\"', "\\\"");
            let tx_json = tx_ids.replace('\"', "\\\"");
            writeln!(
                json_out,
                "{{\"chrom\":\"{}\",\"pos\":{},\"ref\":\"{}\",\"alt\":\"{}\",\"type\":\"{}\",\"status\":\"{}\",\"gene_ids\":\"{}\",\"transcript_ids\":\"{}\",\"feature\":\"{}\"}}",
                chr,
                pos,
                reference,
                alt,
                vtype,
                status,
                gene_json,
                tx_json,
                feature
            )?;
        }
    }
    if let Some(json_out) = out_json.as_mut() {
        writeln!(json_out, "]")?;
    }

    if !out_tsv {
        out.flush()?;
    }
    Ok(0)
}

fn run_scripts(args: &[String]) -> i32 {
    if args.is_empty() || args[0] == "catalog" {
        print_script_catalog();
        return 0;
    }
    if args[0] != "run" {
        eprintln!("unknown scripts command: {}", args[0]);
        return 1;
    }
    if args.len() < 2 {
        eprintln!("missing script-id");
        return 1;
    }

    let script = args[1].as_str();
    let script_args = &args[2..];
    let run = |res: Result<i32>| match res {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{e}");
            1
        }
    };

    match script {
        "change-scaffolds-name" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "-f"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let names = match get_required_opt(script_args, &["-l", "--nameList", "--list"], "name list") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_change_scaffolds_name(&input, &names, get_opt(script_args, &["-o", "--output"]).as_deref()))
        }
        "change-scaffolds-name-fasta" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--inputFile"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let names = match get_required_opt(script_args, &["-l", "--nameList", "--name-list"], "name list") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_change_scaffolds_name_fasta(
                &input,
                &names,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "change-seqname-for-fasta" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--input-dir"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let names = match get_required_opt(script_args, &["-l", "--nameList", "--name-list"], "name list") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out_dir = match get_required_opt(script_args, &["-o", "--outDir", "--output-dir", "--output"], "output directory") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_change_seqname_for_fasta(&input, &names, &out_dir))
        }
        "convert-3line2one" => {
            let input = match get_required_opt(script_args, &["-i", "--input"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_convert_3line2one(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "filter-seq-by-length" => {
            let input = match get_required_opt(script_args, &["-i", "--input"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let len = match parse_usize_arg(script_args, &["-l", "--length"], "length") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_filter_seq_by_length(
                &input,
                len,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "filter-gff-by-id" => {
            let gff = match get_required_opt(script_args, &["-gff", "--inputGff", "--input"], "input gff") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let ids = match get_required_opt(script_args, &["-id", "--IDlist"], "id list") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_filter_gff_by_id(
                &gff,
                &ids,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "filter-gtf-ctg" => {
            let input = match get_required_opt(script_args, &["-i", "--input"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let ids = match get_required_opt(script_args, &["-id", "--idlist", "--id-list"], "id list") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let output = match get_required_opt(script_args, &["-o", "--output"], "output") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_filter_gtf_ctg(&input, &ids, &output))
        }
        "merge-two-txt" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--inputFile"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_merge_two_txt(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "compare-two-blast" => {
            let blast = match get_required_opt(script_args, &["-i", "--blastRes", "--input"], "blast results") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let reverse = match get_required_opt(script_args, &["-r", "--resBlastRes"], "reverse blast results") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_compare_two_blast(
                &blast,
                &reverse,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "get-best-idy" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--inputFile"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_get_best_idy(&input, get_opt(script_args, &["-o", "--output"]).as_deref()))
        }
        "get-best-hit-based-on-idy" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--inputFile"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_get_best_hit_based_on_idy(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "get-best-hit-genes" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--inputFile"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_get_best_hit_genes(&input, get_opt(script_args, &["-o", "--output"]).as_deref()))
        }
        "compare-as-and-noAS" => {
            let no_as = match get_required_opt(script_args, &["-nA", "--noASFile", "--noas"], "noAS file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let as_file = match get_required_opt(script_args, &["-AS", "--ASFile", "--as"], "AS file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_compare_as_and_no_as(
                &no_as,
                &as_file,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "compare-busco-results" => {
            let asc = match get_required_opt(script_args, &["-a", "--ascBusco", "--asc"], "asc busco") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let off = match get_required_opt(script_args, &["-o", "--offSpringBusco", "--offspring"], "offspring busco") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_compare_busco(&asc, &off))
        }
        "merge-fpkm-file" => {
            let dir = match get_required_opt(script_args, &["-i", "--inDir", "--dir", "--input"], "input directory") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out_fpkm = match get_required_opt(script_args, &["-oF", "--outFpkm", "--out"], "outFpkm") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out_profile = match get_required_opt(script_args, &["-oP", "--outProfile", "--out-profile"], "outProfile") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_merge_fpkm_file(&dir, &out_fpkm, &out_profile))
        }
        "save-go" => {
            let input = match get_required_opt(script_args, &["-i", "--input"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_save_go(&input, get_opt(script_args, &["-o", "--output"]).as_deref()))
        }
        "merge-gos" => {
            let swiss = match get_required_opt(script_args, &["-s", "--swissprot", "--swiss"], "swissprot") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let nr = match get_required_opt(script_args, &["-n", "--nr", "--nr"], "nr") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let trembl = match get_required_opt(script_args, &["-T", "--Trembl", "--trembl"], "trembl") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_merge_gos(
                &swiss,
                &nr,
                &trembl,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "merge-blastp-best-jcvi" => {
            let path = match get_required_opt(script_args, &["-p", "--pathDir", "--path", "--dir"], "path") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_merge_blastp_best_jcvi(
                &path,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "orthogenes" => {
            let input = match get_required_opt(script_args, &["-i", "--infile", "--input"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let output = match get_required_opt(script_args, &["-o", "--output"], "output") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_orthogenes(&input, &output))
        }
        "genome-gc" => {
            let input = match get_required_opt(script_args, &["-f", "--fasta", "--input", "--inputFile"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_genome_gc(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "get-the-longest-seq" => {
            let input = match get_required_opt(script_args, &["-i", "--protFasta", "--input"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_get_the_longest_seq(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "get-longest-transcript" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--fasta"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_get_longest_transcript(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "extract-longest-pep" | "extract-gene-family-info-alt" => {
            let input = match get_required_opt(script_args, &["-f", "--fasta", "--fastaFile", "--input"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_extract_longest_pep_from_ensembl_download(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "convert-lastz2-jcvi" => {
            let bed = match get_required_opt(script_args, &["-i", "--bed", "--input"], "bed file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let ref_len = match get_required_opt(script_args, &["-r", "--refLen", "--ref-len"], "reference len file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let query_len = match get_required_opt(script_args, &["-q", "--queryLen", "--query-len"], "query len file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let ref_name = get_opt(script_args, &["--refName"]).unwrap_or_else(|| "Ref".to_string());
            let query_name = get_opt(script_args, &["--queryName"]).unwrap_or_else(|| "Query".to_string());
            run(run_convert_lastz2_jcvi(&bed, &ref_len, &query_len, &ref_name, &query_name))
        }
        "extract-pasa-results" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--pasa"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out_seq = match get_required_opt(
                script_args,
                &["-s", "--outSeq", "--out-seq", "--seqOutput"],
                "output seq",
            ) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out_gff = match get_required_opt(
                script_args,
                &["-g", "--outGff", "--out-gff", "--gffOutput"],
                "output gff",
            ) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_extract_pasa_results(&input, &out_seq, &out_gff))
        }
        "convert-gemoma-gff3" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--gff"], "input gff") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_convert_gemoma_gff3(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "convert-gene-annotation-contigs2chr-PASA" | "convert-gene-annotation-contigs2chr" => {
            let gff = match get_required_opt(script_args, &["-gff", "--input"], "gff") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let background = match get_required_opt(
                script_args,
                &["-b", "--background", "--backgroundFile"],
                "background",
            ) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_convert_gene_annotation_contigs2chr_pasa(
                &gff,
                &background,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "convert-gene-annotation-scaffold2chr-nextgenomics" | "convert-gene-annotation-legacy-alias" => {
            let gff = match get_required_opt(script_args, &["-gff", "--input"], "gff") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let background = match get_required_opt(
                script_args,
                &["-b", "--background", "--backgroundFile"],
                "background",
            ) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_convert_gene_annotation_scaffold2chr_nextgenomics(
                &gff,
                &background,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "filter-gemoma-as" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--gff"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_filter_gemoma_as(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "filter-gemoma-as2" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--gff"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_filter_gemoma_as2(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "get-best-hit-by-score" => {
            let query = match get_required_opt(script_args, &["-i", "--query", "--input"], "query file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let refs = match get_required_opt(script_args, &["-r", "--refs", "--ref"], "reference file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_get_best_hit_by_score(
                &query,
                &refs,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "get-best-hit-by-score-one-file" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--blast"], "blast file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out_prefix = match get_required_opt(
                script_args,
                &["-p", "--outPrefix", "--prefix", "--out-prefix"],
                "output prefix",
            ) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_get_best_hit_by_score_one_file(
                &input,
                &out_prefix,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "get-best-hit-from-blast" | "get-best-hit-for-scan" => {
            let dir = match get_required_opt(script_args, &["-i", "--input", "--dir", "-p", "--path"], "directory") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let num = match parse_usize_arg(script_args, &["-n", "--num", "--numSpecies"], "num species") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_get_best_hit_from_blast(
                &dir,
                num,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "extract-gene-family-info" => {
            let gene_len = match get_required_opt(script_args, &["-l", "--len", "--len-file"], "gene length file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let expr = match get_required_opt(script_args, &["-e", "--expr", "--expr-file"], "expression file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let parent = match get_required_opt(script_args, &["-p", "--parent", "--parent-file"], "parent file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let coverage = match get_required_opt(script_args, &["-c", "--coverage", "--ratio"], "coverage") {
                Ok(v) => {
                    match v.parse::<f64>() {
                        Ok(vv) => vv,
                        Err(_) => {
                            eprintln!("invalid coverage: {v}");
                            return 1;
                        }
                    }
                }
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_extract_gene_family_info(
                &gene_len,
                &expr,
                &parent,
                coverage,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "extract-gene-family-matrix" => {
            let gene_len = match get_required_opt(script_args, &["-l", "--len", "--len-file"], "gene length file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let expr = match get_required_opt(script_args, &["-e", "--expr", "--expr-file"], "expression file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let family_names = match get_required_opt(
                script_args,
                &["-f", "--family-names", "--family-list"],
                "family name file",
            ) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let gene_family = match get_required_opt(
                script_args,
                &["-g", "--gene-family", "--gene-family-file", "--geneFamily"],
                "gene family file",
            ) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let coverage = match get_required_opt(script_args, &["-c", "--coverage", "--ratio"], "coverage") {
                Ok(v) => match v.parse::<f64>() {
                    Ok(vv) => vv,
                    Err(_) => {
                        eprintln!("invalid coverage: {v}");
                        return 1;
                    }
                },
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_extract_gene_family_matrix(
                &gene_len,
                &expr,
                &family_names,
                &gene_family,
                coverage,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "merge-two-end-bam" => {
            let r1 = match get_required_opt(script_args, &["-i", "--r1", "--read1", "--input1"], "r1") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let r2 = match get_required_opt(script_args, &["-j", "--r2", "--read2", "--input2"], "r2") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out1 = get_opt(script_args, &["-o1", "--output1", "--outR1"]);
            let out2 = get_opt(script_args, &["-o2", "--output2", "--outR2"]);
            run(run_merge_two_end_bam(&r1, &r2, out1.as_deref(), out2.as_deref()))
        }
        "merge-two-end-bam1" => {
            let r1 = match get_required_opt(script_args, &["-i", "--r1", "--read1", "--input1"], "r1") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let r2 = match get_required_opt(script_args, &["-j", "--r2", "--read2", "--input2"], "r2") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out1 = get_opt(script_args, &["-o1", "--output1", "--outR1"]);
            let out2 = get_opt(script_args, &["-o2", "--output2", "--outR2"]);
            run(run_merge_two_end_bam1(&r1, &r2, out1.as_deref(), out2.as_deref()))
        }
        "merge-two-end-bam-forMGI" => {
            let r1 = match get_required_opt(script_args, &["-i", "--r1", "--read1", "--input1"], "r1") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let r2 = match get_required_opt(script_args, &["-j", "--r2", "--read2", "--input2"], "r2") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out1 = get_opt(script_args, &["-o1", "--output1", "--outR1"]);
            let out2 = get_opt(script_args, &["-o2", "--output2", "--outR2"]);
            run(run_merge_two_end_bam_for_mgi(&r1, &r2, out1.as_deref(), out2.as_deref()))
        }
        "zhouxiaoxuan-mergexls" => {
            let first = match get_required_opt(script_args, &["-a", "--first", "--file1"], "first file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let second = match get_required_opt(script_args, &["-b", "--second", "--file2"], "second file") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_zhouxiaoxuan_merge_xls(
                &first,
                &second,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "check-duplication-gene-pairs" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--file"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_check_duplication_gene_pairs(
                &input,
                get_opt(script_args, &["-o", "--output"]).as_deref(),
            ))
        }
        "orthofiner-to-pal2nal" => {
            let prot_dir = match get_required_opt(script_args, &["-p", "--pathOfprot", "--path", "--input"], "protein path") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out_dir = match get_required_opt(script_args, &["-o", "--outPutPath", "--output", "--outputDir", "--out"], "output path") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let nucl = match get_required_opt(script_args, &["-n", "--nuclOfcds", "--nucl"], "cds fasta") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            run(run_orthofiner_to_pal2nal(&prot_dir, &out_dir, &nucl))
        }
        "get-diff-sites-from-orthology" | "compare-orthology" => {
            let input = match get_required_opt(script_args, &["-i", "--pal2nalResPath", "--input", "--dir"], "input directory") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out_dir = get_opt(script_args, &["-o", "--output", "--outDir", "--output-dir"]).unwrap_or_else(|| input.clone());
            run(run_get_diff_sites_from_orthology(&input, &out_dir))
        }
        "get-four-degenerate-sites" => {
            let input = match get_required_opt(script_args, &["-i", "--pal2nalResPath", "--input", "--dir"], "input directory") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let out_dir = get_opt(script_args, &["-o", "--output", "--outDir", "--output-dir"]).unwrap_or_else(|| input.clone());
            run(run_get_four_degenerate_sites(&input, &out_dir))
        }
        "plot-depth-pandepth" | "plot-depth-pandepth2" => {
            let input = match get_required_opt(script_args, &["-i", "--input"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let output = match get_required_opt(script_args, &["-o", "--output", "--outDir", "--output-dir"], "output") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let min_len = get_opt(script_args, &["-l", "--min_length", "--min-length"]).and_then(|v| v.parse::<f64>().ok()).unwrap_or(10.0);
            if script == "plot-depth-pandepth" {
                run(run_plot_depth_pandepth(&input, &output, min_len))
            } else {
                run(run_plot_depth_pandepth2(&input, &output, min_len))
            }
        }
        "plot-mosdepth-point" => {
            let input = match get_required_opt(script_args, &["-i", "--input"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let output = match get_required_opt(script_args, &["-o", "--output", "--outDir", "--output-dir"], "output") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let min_len = get_opt(script_args, &["-l", "--min_length", "--min-length"]).and_then(|v| v.parse::<usize>().ok()).unwrap_or(0);
            run(run_plot_mosdepth_point(&input, &output, min_len))
        }
        "trim-ttaggg-fastq" => {
            let input = match get_required_opt(script_args, &["-i", "--input", "--in"], "input") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let output = match get_required_opt(script_args, &["-o", "--output", "--out"], "output") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let motif = get_opt(script_args, &["-s", "--sequence", "--seq"]).unwrap_or_else(|| "TTAGGG".to_string());
            run(run_trim_ttaggg_fastq(&input, &output, &motif))
        }
        "annotation-vcf" => {
            let reference = match get_required_opt(script_args, &["-r", "--reference", "--ref"], "reference fasta") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let gff = match get_required_opt(script_args, &["-g", "--gff"], "gff") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let vcf = match get_required_opt(script_args, &["-v", "--vcf", "--input"] , "vcf") {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            let format = get_opt(script_args, &["-f", "--format"]).unwrap_or_else(|| "tsv".to_string());
            run(run_annotation_vcf(
                &reference,
                &gff,
                &vcf,
                get_opt(script_args, &["-o", "--output", "--pickle"]).as_deref().unwrap_or("annotation_vcf.txt"),
                &format,
            ))
        }
        _ => {
            eprintln!("script not implemented: {script}");
            1
        }
    }
}

fn run_coverage_ratio(input: &str, reference: &str, output: Option<&str>) -> Result<i32> {
    let mut ref_map: HashMap<String, f64> = HashMap::new();
    let mut ref_reader = open_reader(reference)?;
    let mut row = String::new();
    while ref_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']).to_string();
        row.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        if let Ok(v) = cols[1].parse::<f64>() {
            ref_map.insert(cols[0].to_string(), v);
        }
    }

    let mut out = open_writer(output)?;
    let mut in_reader = open_reader(input)?;
    let mut line = String::new();
    while in_reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        if raw.is_empty() {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        if let Some(total) = ref_map.get(cols[0]) {
            if let Ok(v) = cols[1].parse::<f64>() {
                let pct = if *total == 0.0 { 0.0 } else { v / *total };
                writeln!(out, "{}\t{:.6}", cols[0], pct)?;
            }
        }
    }
    Ok(0)
}

fn run_wgcna_weight(weight_file: &str, output: Option<&str>) -> Result<i32> {
    let mut reader = open_reader(weight_file)?;
    let mut header = String::new();
    reader.read_line(&mut header)?;
    let names: Vec<&str> = header.split_whitespace().collect();
    if names.is_empty() {
        return Ok(1);
    }
    let mut out = open_writer(output)?;
    for i in 0..names.len() {
        let mut row = Vec::with_capacity(names.len() + 1);
        row.push(names[i].to_string());
        for j in 0..names.len() {
            row.push(if i == j { "1".into() } else { "0".into() });
        }
        writeln!(out, "{}", row.join("\t"))?;
    }
    Ok(0)
}

fn sort_positions(mut values: Vec<String>) -> Vec<String> {
    let numeric = values.iter().all(|v| v.parse::<i64>().is_ok());
    if numeric {
        values.sort_by_key(|v| v.parse::<i64>().unwrap_or(0));
    } else {
        values.sort();
    }
    values
}

fn run_hic_matrix_reindex(
    bed: &str,
    matrix: &str,
    group: &str,
    output: Option<&str>,
) -> Result<i32> {
    let mut bed_reader = open_reader(bed)?;
    let mut rows: HashMap<String, Vec<String>> = HashMap::new();
    let mut line = String::new();
    while bed_reader.read_line(&mut line)? > 0 {
        let raw = line.trim_end_matches(['\n', '\r']).to_string();
        line.clear();
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 4 {
            continue;
        }
        rows.entry(cols[0].to_string())
            .or_default()
            .push(cols[3].to_string());
    }

    let group_path = expand_path(group);
    let gp = Path::new(&group_path);
    if !gp.exists() {
        return Err(format!("group path not found: {group}").into());
    }

    let mut groups: Vec<(String, String)> = Vec::new();
    if gp.is_file() {
        let mut gr = open_reader(gp.to_string_lossy().as_ref())?;
        let mut gline = String::new();
        while gr.read_line(&mut gline)? > 0 {
            let raw = gline.trim_end_matches(['\n', '\r']).to_string();
            gline.clear();
            if raw.is_empty() || raw.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = raw.split_whitespace().collect();
            if cols.len() < 3 {
                continue;
            }
            groups.push((cols[1].to_string(), cols[2].to_string()));
        }
    } else if gp.is_dir() {
        let mut files: Vec<PathBuf> = Vec::new();
        for e in fs::read_dir(gp)? {
            let e = e?;
            let name = e.file_name().to_string_lossy().to_string();
            if name.starts_with("group") {
                files.push(e.path());
            }
        }
        files.sort();
        for f in files {
            let mut gr = open_reader(&f.to_string_lossy())?;
            let mut gline = String::new();
            while gr.read_line(&mut gline)? > 0 {
                let raw = gline.trim_end_matches(['\n', '\r']).to_string();
                gline.clear();
                if raw.is_empty() || raw.starts_with('#') {
                    continue;
                }
                let cols: Vec<&str> = raw.split_whitespace().collect();
                if cols.len() < 3 {
                    continue;
                }
                groups.push((cols[1].to_string(), cols[2].to_string()));
            }
        }
    }

    let mut remap: HashMap<String, usize> = HashMap::new();
    let mut idx = 1usize;
    for (ctg, dir) in groups {
        if let Some(locations) = rows.get(&ctg) {
            let mut ordered = sort_positions(locations.clone());
            if dir == "1" {
                ordered.reverse();
            }
            for loc in ordered {
                remap.insert(loc, idx);
                idx += 1;
            }
        }
    }

    let mut out = open_writer(output)?;
    let mut matrix_reader = open_reader(matrix)?;
    let mut mline = String::new();
    while matrix_reader.read_line(&mut mline)? > 0 {
        let raw = mline.trim_end_matches(['\n', '\r']).to_string();
        mline.clear();
        if raw.trim_start().starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw.split_whitespace().collect();
        if cols.len() < 2 {
            continue;
        }
        let mut row = Vec::new();
        let left = remap.get(cols[0]).map_or_else(|| cols[0].to_string(), |v| v.to_string());
        let right = remap.get(cols[1]).map_or_else(|| cols[1].to_string(), |v| v.to_string());
        row.push(left);
        row.push(right);
        row.extend(cols.iter().skip(2).map(|s| s.to_string()));
        writeln!(out, "{}", row.join("\t"))?;
    }
    Ok(0)
}

fn run_psmc_merge(dir: &str, pattern: &str, output: &str) -> Result<i32> {
    let mut files: Vec<PathBuf> = Vec::new();
    for entry in fs::read_dir(expand_path(dir))? {
        let entry = entry?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        if path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .ends_with(pattern)
        {
            files.push(path);
        }
    }
    files.sort();

    let mut out = BufWriter::new(File::create(expand_path(output))?);
    writeln!(out, "Sample\tTime\tNe")?;
    for p in files {
        let sample = p
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .split('.')
            .next()
            .unwrap_or("");
        let mut r = open_reader(&p.to_string_lossy())?;
        let mut row = String::new();
        while r.read_line(&mut row)? > 0 {
            let raw = row.trim_end_matches(['\n', '\r']).to_string();
            row.clear();
            let cols: Vec<&str> = raw.split_whitespace().collect();
            if cols.len() < 2 {
                continue;
            }
            writeln!(out, "{sample}\t{}\t{}", cols[0], cols[1])?;
        }
    }
    Ok(0)
}

fn run_rename(args: &[String]) -> i32 {
    if args.is_empty() {
        return 1;
    }
    let sub = args[0].as_str();
    let mut input: Option<String> = None;
    let mut map: Option<String> = None;
    let mut output: Option<String> = None;

    let mut i = 1usize;
    while i < args.len() {
        let flag = args[i].as_str();
        match flag {
            "-i" | "--input" | "-f" => match parse_required(args, &mut i, flag) {
                Ok(v) => input = Some(v),
                Err(_) => return 1,
            },
            "-l" | "--map" => match parse_required(args, &mut i, flag) {
                Ok(v) => map = Some(v),
                Err(_) => return 1,
            },
            "-o" | "--output" => match parse_required(args, &mut i, flag) {
                Ok(v) => output = Some(v),
                Err(_) => return 1,
            },
            _ => {}
        }
        i += 1;
    }

    match sub {
        "hjjn-genes" => {
            let Some(i) = input else { return 1 };
            match run_hjjn_genes(&i, output.as_deref()) {
                Ok(c) => c,
                Err(_) => 1,
            }
        }
        "scaffolds" => {
            let (Some(i), Some(m)) = (input, map) else {
                return 1;
            };
            match run_scaffold_rename(&i, &m, output.as_deref()) {
                Ok(c) => c,
                Err(_) => 1,
            }
        }
        "fasta-scaffolds" => {
            let (Some(i), Some(m)) = (input, map) else {
                return 1;
            };
            match run_fasta_scaffold_rename(&i, &m, output.as_deref()) {
                Ok(c) => c,
                Err(_) => 1,
            }
        }
        _ => 1,
    }
}

fn run_blast(args: &[String]) -> i32 {
    if args.is_empty() || args[0] != "reciprocal" {
        return 1;
    }
    let mut blast: Option<String> = None;
    let mut reverse: Option<String> = None;
    let mut output: Option<String> = None;

    let mut i = 1usize;
    while i < args.len() {
        let flag = args[i].as_str();
        match flag {
            "-i" | "--blast" => match parse_required(args, &mut i, flag) {
                Ok(v) => blast = Some(v),
                Err(_) => return 1,
            },
            "-r" | "--reverse" => match parse_required(args, &mut i, flag) {
                Ok(v) => reverse = Some(v),
                Err(_) => return 1,
            },
            "-o" | "--output" => match parse_required(args, &mut i, flag) {
                Ok(v) => output = Some(v),
                Err(_) => return 1,
            },
            _ => {}
        }
        i += 1;
    }

    let (Some(a), Some(b)) = (blast, reverse) else {
        return 1;
    };
    match run_reciprocal(&a, &b, output.as_deref()) {
        Ok(c) => c,
        Err(_) => 1,
    }
}

fn run_gff(args: &[String]) -> i32 {
    if args.is_empty() {
        return 1;
    }
    let sub = args[0].as_str();

    let mut input: Option<String> = None;
    let mut gff_file: Option<String> = None;
    let mut bed_file: Option<String> = None;
    let mut output: Option<String> = None;

    let mut i = 1usize;
    while i < args.len() {
        let flag = args[i].as_str();
        match flag {
            "-i" | "--input" => match parse_required(args, &mut i, flag) {
                Ok(v) => input = Some(v),
                Err(_) => return 1,
            },
            "-gff" => match parse_required(args, &mut i, "-gff") {
                Ok(v) => gff_file = Some(v),
                Err(_) => return 1,
            },
            "--gff" => match parse_required(args, &mut i, "--gff") {
                Ok(v) => gff_file = Some(v),
                Err(_) => return 1,
            },
            "-b" | "--bed" => match parse_required(args, &mut i, flag) {
                Ok(v) => bed_file = Some(v),
                Err(_) => return 1,
            },
            "-o" | "--output" => match parse_required(args, &mut i, flag) {
                Ok(v) => output = Some(v),
                Err(_) => return 1,
            },
            _ => {}
        }
        i += 1;
    }

    match sub {
        "filter-ncbi" => {
            let output = output.or(Some("TA-filtered.gff3".to_string()));
            let Some(file) = gff_file.or(input) else { return 1 };
            match run_filter_ncbi(&file, output.as_deref()) {
                Ok(c) => c,
                Err(_) => 1,
            }
        }
        "filter-gemoma" => {
            let output = output.or(Some("gemoma-longest.gff3".to_string()));
            let Some(file) = gff_file.or(input) else { return 1 };
            match run_filter_gemoma(&file, output.as_deref()) {
                Ok(c) => c,
                Err(_) => 1,
            }
        }
        "convert-ty1-hjjn" => {
            let Some(g) = gff_file else { return 1 };
            let Some(b) = bed_file else { return 1 };
            let out = output.unwrap_or_else(|| "Results.gff3".to_string());
            match run_convert_ty1_hjjn(&g, &b, &out) {
                Ok(c) => c,
                Err(_) => 1,
            }
        }
        _ => 1,
    }
}

fn run_fasta(args: &[String]) -> i32 {
    if args.is_empty() || args[0] != "longest-transcript" {
        return 1;
    }
    let mut input: Option<String> = None;
    let mut output: Option<String> = None;

    let mut i = 1usize;
    while i < args.len() {
        let flag = args[i].as_str();
        match flag {
            "-f" | "--fasta" => match parse_required(args, &mut i, flag) {
                Ok(v) => input = Some(v),
                Err(_) => return 1,
            },
            "-o" | "--output" => match parse_required(args, &mut i, flag) {
                Ok(v) => output = Some(v),
                Err(_) => return 1,
            },
            _ => {}
        }
        i += 1;
    }
    let Some(fa) = input else { return 1 };
    match run_longest_transcript(&fa, output.as_deref()) {
        Ok(c) => c,
        Err(_) => 1,
    }
}

fn run_stats(args: &[String]) -> i32 {
    if args.is_empty() {
        return 1;
    }
    let sub = args[0].as_str();

    match sub {
        "coverage-ratio" => {
            let mut input: Option<String> = None;
            let mut reference: Option<String> = None;
            let mut output: Option<String> = None;
            let mut i = 1usize;
            while i < args.len() {
                let flag = args[i].as_str();
        match flag {
                    "-i" | "--input" => match parse_required(args, &mut i, flag) {
                        Ok(v) => input = Some(v),
                        Err(_) => return 1,
                    },
                    "-r" | "--reference" => match parse_required(args, &mut i, flag) {
                        Ok(v) => reference = Some(v),
                        Err(_) => return 1,
                    },
                    "-o" | "--output" => match parse_required(args, &mut i, flag) {
                        Ok(v) => output = Some(v),
                        Err(_) => return 1,
                    },
                    _ => {}
                }
                i += 1;
            }
            let Some(i) = input else { return 1 };
            let Some(r) = reference else { return 1 };
            match run_coverage_ratio(&i, &r, output.as_deref()) {
                Ok(c) => c,
                Err(_) => 1,
            }
        }
        "hic-matrix-reindex" => {
            let mut bed: Option<String> = None;
            let mut matrix: Option<String> = None;
            let mut group: Option<String> = None;
            let mut output: Option<String> = None;
            let mut i = 1usize;
            while i < args.len() {
                let flag = args[i].as_str();
        match flag {
                    "-b" | "--bed" => match parse_required(args, &mut i, flag) {
                        Ok(v) => bed = Some(v),
                        Err(_) => return 1,
                    },
                    "-m" | "--matrix" => match parse_required(args, &mut i, flag) {
                        Ok(v) => matrix = Some(v),
                        Err(_) => return 1,
                    },
                    "-p" | "--group" => match parse_required(args, &mut i, flag) {
                        Ok(v) => group = Some(v),
                        Err(_) => return 1,
                    },
                    "-o" | "--output" => match parse_required(args, &mut i, flag) {
                        Ok(v) => output = Some(v),
                        Err(_) => return 1,
                    },
                    _ => {}
                }
                i += 1;
            }
            let (Some(b), Some(m), Some(g)) = (bed, matrix, group) else {
                return 1;
            };
            match run_hic_matrix_reindex(&b, &m, &g, output.as_deref()) {
                Ok(c) => c,
                Err(_) => 1,
            }
        }
        "wgcna-weight" => {
            let mut weight_file: Option<String> = None;
            let mut output: Option<String> = None;
            let mut i = 1usize;
            while i < args.len() {
                let flag = args[i].as_str();
        match flag {
                    "-i" | "--weight-file" => match parse_required(args, &mut i, flag) {
                        Ok(v) => weight_file = Some(v),
                        Err(_) => return 1,
                    },
                    "-o" | "--output" => match parse_required(args, &mut i, flag) {
                        Ok(v) => output = Some(v),
                        Err(_) => return 1,
                    },
                    _ => {}
                }
                i += 1;
            }
            let Some(w) = weight_file else { return 1 };
            match run_wgcna_weight(&w, output.as_deref()) {
                Ok(c) => c,
                Err(_) => 1,
            }
        }
        _ => 1,
    }
}

fn run_psmc(args: &[String]) -> i32 {
    if args.is_empty() || args[0] != "merge" {
        return 1;
    }
    let mut dir: Option<String> = None;
    let mut pattern = ".0.txt".to_string();
    let mut output = "merge.psmc.0.txt".to_string();

    let mut i = 1usize;
    while i < args.len() {
        let flag = args[i].as_str();
        match flag {
            "-d" | "--dir" => match parse_required(args, &mut i, flag) {
                Ok(v) => dir = Some(v),
                Err(_) => return 1,
            },
            "-p" | "--pattern" => match parse_required(args, &mut i, flag) {
                Ok(v) => pattern = v,
                Err(_) => return 1,
            },
            "-o" | "--output" => match parse_required(args, &mut i, flag) {
                Ok(v) => output = v,
                Err(_) => return 1,
            },
            _ => {}
        }
        i += 1;
    }

    let Some(d) = dir else {
        return 1;
    };
    match run_psmc_merge(&d, &pattern, &output) {
        Ok(c) => c,
        Err(_) => 1,
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() <= 1 || needs_help(&args[1..]) {
        print_help();
        if args.len() <= 1 {
            std::process::exit(1);
        }
        return;
    }
    if args[1] == "--version" || args[1] == "-V" {
        println!("{VERSION}");
        return;
    }

    let code = match args[1].as_str() {
        "rename" => run_rename(&args[2..]),
        "blast" => run_blast(&args[2..]),
        "gff" => run_gff(&args[2..]),
        "fasta" => run_fasta(&args[2..]),
        "stats" => run_stats(&args[2..]),
        "scripts" => run_scripts(&args[2..]),
        "psmc" => run_psmc(&args[2..]),
        _ => {
            print_help();
            1
        }
    };

    std::process::exit(code);
}
