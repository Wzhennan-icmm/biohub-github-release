use std::process::Command;

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
