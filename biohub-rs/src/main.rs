use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

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
        let row = line.trim_end_matches(['\n', '\r']);
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
        let raw = row.trim_end_matches(['\n', '\r']);
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
        let t = row.trim_end_matches(['\n', '\r']).trim();
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
        let raw = line.trim_end_matches(['\n', '\r']);
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
        let raw = line.trim_end_matches(['\n', '\r']);
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
            let raw = row.trim_end_matches(['\n', '\r']);
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
        let raw = row.trim_end_matches(['\n', '\r']);
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
        let raw = gline.trim_end_matches(['\n', '\r']);
        gline.clear();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let mut cols: Vec<&str> = raw.split('\t').collect();
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
            let head = rest.trim().to_string();
            let to_emit = current_header.replace(head);
            if let Some(h) = to_emit {
                finalize_record(h, std::mem::take(&mut seq_lines), &mut selected, &mut seq_store, &mut order);
            }
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

fn run_coverage_ratio(input: &str, reference: &str, output: Option<&str>) -> Result<i32> {
    let mut ref_map: HashMap<String, f64> = HashMap::new();
    let mut ref_reader = open_reader(reference)?;
    let mut row = String::new();
    while ref_reader.read_line(&mut row)? > 0 {
        let raw = row.trim_end_matches(['\n', '\r']);
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
        let raw = line.trim_end_matches(['\n', '\r']);
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
                writeln!(out, "{}\t{pct:.6f}", cols[0])?;
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

    let gp = Path::new(&expand_path(group));
    if !gp.exists() {
        return Err(format!("group path not found: {group}").into());
    }

    let mut groups: Vec<(String, String)> = Vec::new();
    if gp.is_file() {
        let mut gr = open_reader(gp.to_string_lossy().as_ref())?;
        let mut gline = String::new();
        while gr.read_line(&mut gline)? > 0 {
            let raw = gline.trim_end_matches(['\n', '\r']);
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
                let raw = gline.trim_end_matches(['\n', '\r']);
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
        let raw = mline.trim_end_matches(['\n', '\r']);
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
            let raw = row.trim_end_matches(['\n', '\r']);
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
        match args[i].as_str() {
            "-i" | "--blast" => match parse_required(args, &mut i, args[i].as_str()) {
                Ok(v) => blast = Some(v),
                Err(_) => return 1,
            },
            "-r" | "--reverse" => match parse_required(args, &mut i, args[i].as_str()) {
                Ok(v) => reverse = Some(v),
                Err(_) => return 1,
            },
            "-o" | "--output" => match parse_required(args, &mut i, args[i].as_str()) {
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
        match args[i].as_str() {
            "-i" | "--input" => match parse_required(args, &mut i, args[i].as_str()) {
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
            "-b" | "--bed" => match parse_required(args, &mut i, args[i].as_str()) {
                Ok(v) => bed_file = Some(v),
                Err(_) => return 1,
            },
            "-o" | "--output" => match parse_required(args, &mut i, args[i].as_str()) {
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
        match args[i].as_str() {
            "-f" | "--fasta" => match parse_required(args, &mut i, args[i].as_str()) {
                Ok(v) => input = Some(v),
                Err(_) => return 1,
            },
            "-o" | "--output" => match parse_required(args, &mut i, args[i].as_str()) {
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
                match args[i].as_str() {
                    "-i" | "--input" => match parse_required(args, &mut i, args[i].as_str()) {
                        Ok(v) => input = Some(v),
                        Err(_) => return 1,
                    },
                    "-r" | "--reference" => match parse_required(args, &mut i, args[i].as_str()) {
                        Ok(v) => reference = Some(v),
                        Err(_) => return 1,
                    },
                    "-o" | "--output" => match parse_required(args, &mut i, args[i].as_str()) {
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
                match args[i].as_str() {
                    "-b" | "--bed" => match parse_required(args, &mut i, args[i].as_str()) {
                        Ok(v) => bed = Some(v),
                        Err(_) => return 1,
                    },
                    "-m" | "--matrix" => match parse_required(args, &mut i, args[i].as_str()) {
                        Ok(v) => matrix = Some(v),
                        Err(_) => return 1,
                    },
                    "-p" | "--group" => match parse_required(args, &mut i, args[i].as_str()) {
                        Ok(v) => group = Some(v),
                        Err(_) => return 1,
                    },
                    "-o" | "--output" => match parse_required(args, &mut i, args[i].as_str()) {
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
                match args[i].as_str() {
                    "-i" | "--weight-file" => match parse_required(args, &mut i, args[i].as_str()) {
                        Ok(v) => weight_file = Some(v),
                        Err(_) => return 1,
                    },
                    "-o" | "--output" => match parse_required(args, &mut i, args[i].as_str()) {
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
        match args[i].as_str() {
            "-d" | "--dir" => match parse_required(args, &mut i, args[i].as_str()) {
                Ok(v) => dir = Some(v),
                Err(_) => return 1,
            },
            "-p" | "--pattern" => match parse_required(args, &mut i, args[i].as_str()) {
                Ok(v) => pattern = v,
                Err(_) => return 1,
            },
            "-o" | "--output" => match parse_required(args, &mut i, args[i].as_str()) {
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
        "psmc" => run_psmc(&args[2..]),
        _ => {
            print_help();
            1
        }
    };

    std::process::exit(code);
}
