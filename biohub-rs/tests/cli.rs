use std::process::Command;

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

    assert_eq!(ids.len(), 54, "unexpected catalog size");
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
    assert!(CHINESE_USER_GUIDE.contains("文档版本：0.3.0"));
    assert!(!CHINESE_USER_GUIDE.contains("/Users/"));
}
