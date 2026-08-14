use std::process::Command;
use std::{env, fs};

const SCRIPT_CATALOG: &str = include_str!("../src/script_catalog.tsv");
const CHINESE_USER_GUIDE: &str = include_str!("../../docs/USER_GUIDE.zh-CN.md");
const ROOT_README: &str = include_str!("../../README.md");

fn biohub(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_biohub"))
        .args(args)
        .output()
        .expect("run biohub test binary")
}

#[test]
fn help_surfaces_exit_successfully() {
    for args in [
        vec!["--help"],
        vec!["catalog", "--help"],
        vec!["doctor", "--help"],
        vec!["recipe", "--help"],
        vec!["r", "--help"],
        vec!["run", "--help"],
        vec!["scripts", "--help"],
        vec!["run", "dotplot", "--help"],
    ] {
        let output = biohub(&args);
        assert!(output.status.success(), "help failed for {args:?}");
        assert!(!output.stdout.is_empty(), "help was empty for {args:?}");
    }
}

#[test]
fn catalog_json_exposes_backend_metadata() {
    let output = biohub(&["catalog", "--format", "json"]);
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("catalog is UTF-8");
    assert!(stdout.starts_with("[\n"));
    assert!(stdout.contains("\"id\":\"dotplot\""));
    assert!(stdout.contains("\"backend\":\"r\""));
    assert!(stdout.contains("\"dependencies\":[\"Rscript\"]"));
    assert!(stdout.contains("\"id\":\"selection-branch-site\""));
    assert_eq!(stdout.matches("\"kind\":\"recipe\"").count(), 13);
}

#[test]
fn recipe_init_copies_template_schema_and_readme() {
    let directory = env::temp_dir().join(format!(
        "biohub_recipe_init_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ));
    let recipe_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("repository root")
        .join("recipes");
    let output = Command::new(env!("CARGO_BIN_EXE_biohub"))
        .args([
            "recipe",
            "init",
            "family-denovo-rate",
            "--workdir",
            directory.to_str().expect("workdir path"),
        ])
        .env("BIOHUB_RECIPE_DIR", recipe_root)
        .output()
        .expect("run recipe init");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(directory.join("config.yaml").is_file());
    assert!(directory.join("config.schema.yaml").is_file());
    assert!(directory.join("README.txt").is_file());
    assert_eq!(
        fs::read_dir(&directory)
            .expect("list init directory")
            .count(),
        3
    );
    fs::remove_dir_all(directory).expect("remove test directory");
}

#[test]
fn unknown_command_returns_failure() {
    let output = biohub(&["not-a-command"]);
    assert!(!output.status.success());
}

#[test]
fn chinese_user_guide_documents_every_catalog_command_once() {
    let ids: Vec<&str> = SCRIPT_CATALOG
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .filter_map(|line| line.split('\t').next())
        .collect();

    assert_eq!(ids.len(), 57, "unexpected catalog size");
    for id in ids {
        let heading = format!("### `{id}`");
        let occurrences = CHINESE_USER_GUIDE
            .lines()
            .filter(|line| *line == heading)
            .count();
        assert_eq!(
            occurrences, 1,
            "guide must contain exactly one reference heading for {id}"
        );
    }
}

#[test]
fn readme_links_chinese_user_guide_without_local_paths() {
    assert!(ROOT_README.contains("docs/USER_GUIDE.zh-CN.md"));
    assert!(CHINESE_USER_GUIDE.contains("文档版本：0.4.0"));
    assert!(!CHINESE_USER_GUIDE.contains("/Users/"));
}

#[test]
fn filter_gff_by_fasta_drops_invalid_models_and_descendants() {
    let directory = env::temp_dir().join(format!(
        "biohub_filter_gff_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ));
    fs::create_dir(&directory).expect("create test directory");
    let lengths = directory.join("reference.fai");
    let gff = directory.join("input.gff3");
    let output = directory.join("filtered.gff3");
    fs::write(&lengths, "chr1\t100\t0\t0\t0\n").expect("write lengths");
    fs::write(
        &gff,
        concat!(
            "##gff-version 3\n",
            "chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=g1\n",
            "chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n",
            "chr1\ttest\texon\t80\t110\t.\t+\t.\tParent=t1\n",
            "chr1\ttest\tgene\t1\t120\t.\t+\t.\tID=g2\n",
            "chr1\ttest\tmRNA\t1\t80\t.\t+\t.\tID=t2;Parent=g2\n",
            "chr1\ttest\tgene\t1\t80\t.\t+\t.\tID=g3\n",
            "chr1\ttest\tmRNA\t1\t80\t.\t+\t.\tID=t3;Parent=g3\n",
            "malformed row retained\n",
        ),
    )
    .expect("write GFF");

    let result = biohub(&[
        "run",
        "filter-gff-by-fasta",
        "--gff",
        gff.to_str().expect("GFF path"),
        "--fai",
        lengths.to_str().expect("length path"),
        "--output",
        output.to_str().expect("output path"),
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let filtered = fs::read_to_string(&output).expect("read filtered GFF");
    assert!(filtered.contains("ID=g1"));
    assert!(!filtered.contains("ID=t1"));
    assert!(!filtered.contains("ID=g2"));
    assert!(!filtered.contains("ID=t2"));
    assert!(filtered.contains("ID=g3"));
    assert!(filtered.contains("ID=t3"));
    assert!(filtered.contains("malformed row retained"));
    assert!(String::from_utf8_lossy(&result.stderr).contains("skipped_records=4"));

    fs::remove_dir_all(directory).expect("remove test directory");
}
