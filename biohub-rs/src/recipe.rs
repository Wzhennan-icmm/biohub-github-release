use std::env;
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Component, Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

type RecipeResult<T> = std::result::Result<T, String>;

const RECIPE_CATALOG_TEXT: &str = include_str!("recipe_catalog.tsv");

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RecipeSpec {
    pub id: String,
    pub domain: String,
    pub description: String,
    pub status: String,
    pub workflow: String,
    pub config_template: String,
    pub config_schema: String,
    pub dependencies: Vec<String>,
    pub container: String,
    pub license: String,
}

#[derive(Debug)]
struct RunOptions {
    recipe_id: String,
    config: PathBuf,
    workdir: PathBuf,
    profile: String,
    cores: usize,
    dry_run: bool,
    resume: bool,
}

#[derive(Debug)]
struct PreparedRun {
    recipe: RecipeSpec,
    workflow: PathBuf,
    config: PathBuf,
    workdir: PathBuf,
    profile: Option<PathBuf>,
    config_sha256: String,
    started_unix: u64,
}

#[derive(Debug)]
struct RunState<'a> {
    recipe_id: &'a str,
    status: &'a str,
    config_sha256: &'a str,
    started_unix: u64,
    updated_unix: u64,
    exit_code: Option<i32>,
    message: &'a str,
}

pub(crate) fn load_recipe_catalog() -> RecipeResult<Vec<RecipeSpec>> {
    let mut recipes = Vec::new();
    let mut ids = std::collections::HashSet::new();
    for (line_number, line) in RECIPE_CATALOG_TEXT.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() != 10 {
            return Err(format!(
                "invalid recipe catalog row {}: expected 10 tab-separated fields, found {}",
                line_number + 1,
                cells.len()
            ));
        }
        let id = cells[0].trim();
        if !valid_id(id) {
            return Err(format!(
                "invalid recipe id at row {}: {id}",
                line_number + 1
            ));
        }
        if !ids.insert(id.to_string()) {
            return Err(format!("duplicate recipe id: {id}"));
        }
        for value in [cells[4], cells[5], cells[6]] {
            validate_relative_path(value)?;
        }
        recipes.push(RecipeSpec {
            id: id.to_string(),
            domain: cells[1].to_string(),
            description: cells[2].to_string(),
            status: cells[3].to_string(),
            workflow: cells[4].to_string(),
            config_template: cells[5].to_string(),
            config_schema: cells[6].to_string(),
            dependencies: split_csv(cells[7]),
            container: cells[8].to_string(),
            license: cells[9].to_string(),
        });
    }
    Ok(recipes)
}

fn split_csv(value: &str) -> Vec<String> {
    value
        .split(',')
        .map(str::trim)
        .filter(|part| !part.is_empty() && *part != "-")
        .map(ToString::to_string)
        .collect()
}

fn valid_id(value: &str) -> bool {
    !value.is_empty()
        && value
            .chars()
            .all(|ch| ch.is_ascii_lowercase() || ch.is_ascii_digit() || ch == '-')
}

fn validate_relative_path(value: &str) -> RecipeResult<()> {
    let path = Path::new(value);
    if value.is_empty() || path.is_absolute() {
        return Err(format!(
            "recipe path must be non-empty and relative: {value}"
        ));
    }
    if path.components().any(|component| {
        matches!(
            component,
            Component::ParentDir | Component::RootDir | Component::Prefix(_)
        )
    }) {
        return Err(format!("recipe path cannot escape recipe root: {value}"));
    }
    Ok(())
}

fn recipe_by_id(id: &str) -> RecipeResult<RecipeSpec> {
    load_recipe_catalog()?
        .into_iter()
        .find(|recipe| recipe.id == id)
        .ok_or_else(|| format!("unknown recipe-id: {id}"))
}

fn recipe_root() -> RecipeResult<PathBuf> {
    let mut candidates = Vec::new();
    if let Ok(value) = env::var("BIOHUB_RECIPE_DIR") {
        candidates.push(PathBuf::from(value));
    }
    if let Ok(cwd) = env::current_dir() {
        candidates.push(cwd.join("recipes"));
        if let Some(parent) = cwd.parent() {
            candidates.push(parent.join("recipes"));
        }
    }
    if let Ok(executable) = env::current_exe() {
        if let Some(bin_dir) = executable.parent() {
            candidates.push(bin_dir.join("../share/biohub/recipes"));
        }
    }

    for candidate in candidates {
        if candidate.is_dir() {
            return candidate
                .canonicalize()
                .map_err(|error| format!("cannot resolve recipe directory: {error}"));
        }
    }
    Err("recipe directory not found; set BIOHUB_RECIPE_DIR".to_string())
}

fn recipe_file(root: &Path, relative: &str) -> RecipeResult<PathBuf> {
    validate_relative_path(relative)?;
    let path = root.join(relative);
    if !path.is_file() {
        return Err(format!(
            "packaged recipe file not found: {}",
            path.display()
        ));
    }
    Ok(path)
}

fn json_escape(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    for ch in input.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

fn json_string_array(values: &[String]) -> String {
    values
        .iter()
        .map(|value| format!("\"{}\"", json_escape(value)))
        .collect::<Vec<_>>()
        .join(",")
}

pub(crate) fn recipe_catalog_json_objects() -> RecipeResult<Vec<String>> {
    Ok(load_recipe_catalog()?
        .iter()
        .map(|recipe| {
            format!(
                "{{\"id\":\"{}\",\"kind\":\"recipe\",\"domain\":\"{}\",\"description\":\"{}\",\"status\":\"{}\",\"backend\":\"snakemake\",\"workflow\":\"{}\",\"config_schema\":\"{}\",\"dependencies\":[{}],\"container\":\"{}\",\"version\":\"{}\",\"license\":\"{}\"}}",
                json_escape(&recipe.id),
                json_escape(&recipe.domain),
                json_escape(&recipe.description),
                json_escape(&recipe.status),
                json_escape(&recipe.workflow),
                json_escape(&recipe.config_schema),
                json_string_array(&recipe.dependencies),
                json_escape(&recipe.container),
                env!("CARGO_PKG_VERSION"),
                json_escape(&recipe.license)
            )
        })
        .collect())
}

fn print_list_table() -> RecipeResult<()> {
    println!("{:<32} {:<24} {:<13} DESCRIPTION", "ID", "DOMAIN", "STATUS");
    println!("{:-<112}", "");
    for recipe in load_recipe_catalog()? {
        println!(
            "{:<32} {:<24} {:<13} {}",
            recipe.id, recipe.domain, recipe.status, recipe.description
        );
    }
    Ok(())
}

pub(crate) fn print_recipe_catalog_table() -> RecipeResult<()> {
    print_list_table()
}

fn print_list_json() -> RecipeResult<()> {
    println!("[");
    let objects = recipe_catalog_json_objects()?;
    for (index, object) in objects.iter().enumerate() {
        println!(
            "  {}{}",
            object,
            if index + 1 == objects.len() { "" } else { "," }
        );
    }
    println!("]");
    Ok(())
}

fn print_description(recipe: &RecipeSpec, json: bool) {
    if json {
        let dependencies = json_string_array(&recipe.dependencies);
        println!(
            "{{\"id\":\"{}\",\"kind\":\"recipe\",\"domain\":\"{}\",\"description\":\"{}\",\"status\":\"{}\",\"backend\":\"snakemake\",\"workflow\":\"{}\",\"config_template\":\"{}\",\"config_schema\":\"{}\",\"dependencies\":[{}],\"container\":\"{}\",\"version\":\"{}\",\"license\":\"{}\"}}",
            json_escape(&recipe.id),
            json_escape(&recipe.domain),
            json_escape(&recipe.description),
            json_escape(&recipe.status),
            json_escape(&recipe.workflow),
            json_escape(&recipe.config_template),
            json_escape(&recipe.config_schema),
            dependencies,
            json_escape(&recipe.container),
            env!("CARGO_PKG_VERSION"),
            json_escape(&recipe.license)
        );
        return;
    }
    println!("Recipe: {}", recipe.id);
    println!("Domain: {}", recipe.domain);
    println!("Status: {}", recipe.status);
    println!("Backend: Snakemake");
    println!("Description: {}", recipe.description);
    println!("Dependencies: {}", recipe.dependencies.join(", "));
    println!("Container: {}", recipe.container);
    println!("BioHub version: {}", env!("CARGO_PKG_VERSION"));
    println!("License: {}", recipe.license);
    println!("Config template: {}", recipe.config_template);
    println!("Config schema: {}", recipe.config_schema);
}

fn option_value(args: &[String], flag: &str) -> RecipeResult<Option<String>> {
    let mut found = None;
    let mut index = 0usize;
    while index < args.len() {
        if args[index] == flag {
            if index + 1 >= args.len() {
                return Err(format!("missing value for {flag}"));
            }
            if found.is_some() {
                return Err(format!("option specified more than once: {flag}"));
            }
            found = Some(args[index + 1].clone());
            index += 2;
        } else {
            index += 1;
        }
    }
    Ok(found)
}

fn required_option(args: &[String], flag: &str) -> RecipeResult<String> {
    option_value(args, flag)?.ok_or_else(|| format!("missing required option: {flag}"))
}

fn contains_flag(args: &[String], flag: &str) -> bool {
    args.iter().any(|value| value == flag)
}

fn reject_unknown_options(
    args: &[String],
    value_flags: &[&str],
    bool_flags: &[&str],
) -> RecipeResult<()> {
    let mut index = 0usize;
    while index < args.len() {
        let value = args[index].as_str();
        if value_flags.contains(&value) {
            if index + 1 >= args.len() {
                return Err(format!("missing value for {value}"));
            }
            index += 2;
        } else if bool_flags.contains(&value) {
            index += 1;
        } else {
            return Err(format!("unknown option: {value}"));
        }
    }
    Ok(())
}

fn init_recipe(recipe: &RecipeSpec, workdir: &Path) -> RecipeResult<()> {
    let root = recipe_root()?;
    let source = recipe_file(&root, &recipe.config_template)?;
    let schema_source = recipe_file(&root, &recipe.config_schema)?;
    let created = if workdir.exists() {
        if !workdir.is_dir() {
            return Err(format!(
                "workdir exists and is not a directory: {}",
                workdir.display()
            ));
        }
        if workdir
            .read_dir()
            .map_err(|error| format!("cannot inspect workdir: {error}"))?
            .next()
            .is_some()
        {
            return Err(format!("workdir must be empty: {}", workdir.display()));
        }
        false
    } else {
        fs::create_dir_all(workdir)
            .map_err(|error| format!("cannot create workdir {}: {error}", workdir.display()))?;
        true
    };

    let config_path = workdir.join("config.yaml");
    let schema_path = workdir.join("config.schema.yaml");
    let readme_path = workdir.join("README.txt");
    let result = (|| {
        copy_create_new(&source, &config_path)?;
        copy_create_new(&schema_source, &schema_path)?;
        let mut readme = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&readme_path)
            .map_err(|error| format!("cannot create README.txt: {error}"))?;
        writeln!(readme, "BioHub recipe: {}", recipe.id).map_err(|e| e.to_string())?;
        writeln!(readme, "{}", recipe.description).map_err(|e| e.to_string())?;
        writeln!(readme).map_err(|e| e.to_string())?;
        writeln!(
            readme,
            "1. Edit config.yaml; config.schema.yaml defines accepted fields."
        )
        .map_err(|e| e.to_string())?;
        writeln!(
            readme,
            "2. Validate: biohub recipe validate {} --config config.yaml",
            recipe.id
        )
        .map_err(|e| e.to_string())?;
        writeln!(
            readme,
            "3. Run: biohub recipe run {} --config config.yaml --workdir . --cores 1",
            recipe.id
        )
        .map_err(|e| e.to_string())?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&readme_path);
        let _ = fs::remove_file(&schema_path);
        let _ = fs::remove_file(&config_path);
        if created {
            let _ = fs::remove_dir(workdir);
        }
    }
    result
}

fn copy_create_new(source: &Path, destination: &Path) -> RecipeResult<()> {
    let mut input =
        File::open(source).map_err(|error| format!("cannot open {}: {error}", source.display()))?;
    let mut output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(destination)
        .map_err(|error| format!("cannot create {}: {error}", destination.display()))?;
    match std::io::copy(&mut input, &mut output) {
        Ok(_) => Ok(()),
        Err(error) => {
            drop(output);
            let _ = fs::remove_file(destination);
            Err(format!("cannot copy {}: {error}", destination.display()))
        }
    }
}

fn now_unix() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs()
}

fn unique_validation_directory(recipe_id: &str) -> RecipeResult<PathBuf> {
    let base = env::temp_dir();
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    for attempt in 0..100usize {
        let path = base.join(format!(
            "biohub_validate_{}_{}_{}_{}",
            recipe_id,
            std::process::id(),
            nanos,
            attempt
        ));
        match fs::create_dir(&path) {
            Ok(()) => return Ok(path),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(format!(
                    "cannot create validation directory {}: {error}",
                    path.display()
                ));
            }
        }
    }
    Err("cannot allocate unique validation directory".to_string())
}

fn sha256_file(path: &Path) -> RecipeResult<String> {
    let attempts: [(&str, &[&str]); 2] = [("sha256sum", &[]), ("shasum", &["-a", "256"])];
    for (program, prefix_args) in attempts {
        let output = Command::new(program).args(prefix_args).arg(path).output();
        let Ok(output) = output else {
            continue;
        };
        if !output.status.success() {
            continue;
        }
        let stdout = String::from_utf8_lossy(&output.stdout);
        let Some(value) = stdout.split_whitespace().next() else {
            continue;
        };
        if value.len() == 64 && value.chars().all(|ch| ch.is_ascii_hexdigit()) {
            return Ok(value.to_ascii_lowercase());
        }
    }
    Err("cannot calculate SHA256: install sha256sum or shasum".to_string())
}

fn shell_quote(value: &str) -> String {
    if value.is_empty() {
        return "''".to_string();
    }
    format!("'{}'", value.replace('\'', "'\"'\"'"))
}

fn write_atomic(path: &Path, content: &str) -> RecipeResult<()> {
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| format!("invalid output path: {}", path.display()))?;
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    let temporary =
        path.with_file_name(format!(".{file_name}.{}.{}.tmp", std::process::id(), nonce));
    fs::write(&temporary, content)
        .map_err(|error| format!("cannot write {}: {error}", temporary.display()))?;
    fs::rename(&temporary, path)
        .map_err(|error| format!("cannot finalize {}: {error}", path.display()))?;
    Ok(())
}

fn write_run_state(workdir: &Path, state: &RunState<'_>) -> RecipeResult<()> {
    let exit_code = state
        .exit_code
        .map(|code| code.to_string())
        .unwrap_or_else(|| "null".to_string());
    let content = format!(
        "{{\n  \"schema_version\": 1,\n  \"recipe_id\": \"{}\",\n  \"status\": \"{}\",\n  \"config_sha256\": \"{}\",\n  \"started_unix\": {},\n  \"updated_unix\": {},\n  \"exit_code\": {},\n  \"message\": \"{}\"\n}}\n",
        json_escape(state.recipe_id),
        json_escape(state.status),
        json_escape(state.config_sha256),
        state.started_unix,
        state.updated_unix,
        exit_code,
        json_escape(state.message)
    );
    write_atomic(&workdir.join("run.json"), &content)
}

fn read_trimmed(path: &Path) -> RecipeResult<String> {
    fs::read_to_string(path)
        .map(|value| value.trim().to_string())
        .map_err(|error| format!("cannot read {}: {error}", path.display()))
}

fn read_started_unix(path: &Path) -> RecipeResult<u64> {
    let content = fs::read_to_string(path)
        .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
    let prefix = "\"started_unix\":";
    for line in content.lines() {
        let value = line.trim().strip_prefix(prefix).map(str::trim);
        if let Some(value) = value {
            return value
                .trim_end_matches(',')
                .parse::<u64>()
                .map_err(|_| format!("invalid started_unix in {}", path.display()));
        }
    }
    Err(format!("missing started_unix in {}", path.display()))
}

fn resolve_profile(root: &Path, profile: &str) -> RecipeResult<Option<PathBuf>> {
    if profile == "local" {
        return Ok(None);
    }
    let candidate = if profile == "slurm" {
        root.join("profiles/slurm")
    } else {
        PathBuf::from(profile)
    };
    if !candidate.is_dir() {
        return Err(format!(
            "Snakemake profile not found: {}",
            candidate.display()
        ));
    }
    candidate
        .canonicalize()
        .map(Some)
        .map_err(|error| format!("cannot resolve profile: {error}"))
}

fn snakemake_version() -> String {
    match Command::new("snakemake").arg("--version").output() {
        Ok(output) if output.status.success() => {
            String::from_utf8_lossy(&output.stdout).trim().to_string()
        }
        _ => "unavailable".to_string(),
    }
}

fn normalized_version_output(output: std::process::Output) -> String {
    if !output.status.success() {
        return "version probe failed".to_string();
    }
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    let value = if stdout.trim().is_empty() {
        stderr.trim()
    } else {
        stdout.trim()
    };
    if value.is_empty() {
        return "available; version not reported".to_string();
    }
    value
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .collect::<Vec<_>>()
        .join(" | ")
        .replace('\t', " ")
}

fn dependency_version(dependency: &str) -> String {
    if let Some(package) = dependency.strip_prefix("R-package:") {
        if !dependency_available(dependency) {
            return "unavailable".to_string();
        }
        return match Command::new("Rscript")
            .arg("--vanilla")
            .arg("-e")
            .arg("cat(as.character(utils::packageVersion(commandArgs(TRUE)[1])))")
            .arg(package)
            .output()
        {
            Ok(output) => normalized_version_output(output),
            Err(_) => "version probe failed".to_string(),
        };
    }
    if !command_available(dependency) {
        return "unavailable".to_string();
    }
    let supports_version_flag = matches!(
        dependency,
        "snakemake"
            | "mafft"
            | "cafe5"
            | "minimap2"
            | "syri"
            | "plink2"
            | "vcftools"
            | "Rscript"
            | "python2"
            | "python3"
    );
    if !supports_version_flag {
        return "available; version probe not defined".to_string();
    }
    match Command::new(dependency).arg("--version").output() {
        Ok(output) => normalized_version_output(output),
        Err(_) => "version probe failed".to_string(),
    }
}

fn snakemake_command(
    workflow: &Path,
    config: &Path,
    workdir: &Path,
    profile: Option<&Path>,
    cores: usize,
    dry_run: bool,
) -> Command {
    let mut command = Command::new("snakemake");
    command
        .arg("--snakefile")
        .arg(workflow)
        .arg("--configfile")
        .arg(config)
        .arg("--directory")
        .arg(workdir)
        .arg("--cores")
        .arg(cores.to_string())
        .arg("--printshellcmds")
        .arg("--rerun-incomplete");
    if let Some(profile) = profile {
        command.arg("--profile").arg(profile);
    }
    if dry_run {
        command.arg("--dry-run");
    }
    command
}

fn command_as_shell(
    workflow: &Path,
    config: &Path,
    workdir: &Path,
    profile: Option<&Path>,
    cores: usize,
    dry_run: bool,
) -> String {
    let mut parts = vec![
        "snakemake".to_string(),
        "--snakefile".to_string(),
        shell_quote(&workflow.to_string_lossy()),
        "--configfile".to_string(),
        shell_quote(&config.to_string_lossy()),
        "--directory".to_string(),
        shell_quote(&workdir.to_string_lossy()),
        "--cores".to_string(),
        cores.to_string(),
        "--printshellcmds".to_string(),
        "--rerun-incomplete".to_string(),
    ];
    if let Some(profile) = profile {
        parts.push("--profile".to_string());
        parts.push(shell_quote(&profile.to_string_lossy()));
    }
    if dry_run {
        parts.push("--dry-run".to_string());
    }
    format!(
        "#!/usr/bin/env bash\nset -euo pipefail\n{}\n",
        parts.join(" ")
    )
}

fn prepare_run(options: &RunOptions) -> RecipeResult<PreparedRun> {
    if options.cores == 0 {
        return Err("--cores must be greater than zero".to_string());
    }
    let recipe = recipe_by_id(&options.recipe_id)?;
    let root = recipe_root()?;
    let workflow = recipe_file(&root, &recipe.workflow)?
        .canonicalize()
        .map_err(|error| format!("cannot resolve workflow: {error}"))?;
    let _schema = recipe_file(&root, &recipe.config_schema)?;
    let config = options.config.canonicalize().map_err(|error| {
        format!(
            "cannot resolve config {}: {error}",
            options.config.display()
        )
    })?;
    if !options.workdir.exists() {
        fs::create_dir_all(&options.workdir).map_err(|error| {
            format!(
                "cannot create workdir {}: {error}",
                options.workdir.display()
            )
        })?;
    }
    if !options.workdir.is_dir() {
        return Err(format!(
            "workdir is not a directory: {}",
            options.workdir.display()
        ));
    }
    let workdir = options
        .workdir
        .canonicalize()
        .map_err(|error| format!("cannot resolve workdir: {error}"))?;
    let profile = resolve_profile(&root, &options.profile)?;
    let config_sha256 = sha256_file(&config)?;
    let mut started = now_unix();

    let state_path = workdir.join("run.json");
    if state_path.exists() {
        if !options.resume {
            return Err(format!(
                "run already initialized in {}; use --resume with unchanged config or choose a new workdir",
                workdir.display()
            ));
        }
        let prior_recipe = read_trimmed(&workdir.join("recipe.id"))?;
        let prior_hash = read_trimmed(&workdir.join("config.sha256"))?;
        if prior_recipe != recipe.id {
            return Err(format!(
                "resume recipe mismatch: existing {prior_recipe}, requested {}",
                recipe.id
            ));
        }
        if prior_hash != config_sha256 {
            return Err(
                "resume config mismatch; create a new workdir for changed configuration"
                    .to_string(),
            );
        }
        started = read_started_unix(&state_path)?;
    } else if options.resume {
        return Err("--resume requested but workdir has no run.json".to_string());
    } else {
        let managed_paths = [
            "config.resolved.yaml",
            "recipe.id",
            "config.sha256",
            "command.sh",
            "versions.tsv",
            "provenance.json",
            "inputs.manifest.tsv",
            "checksums.sha256",
            "logs",
            "results",
            "report",
            ".snakemake",
        ];
        let collisions: Vec<_> = managed_paths
            .iter()
            .filter(|relative| workdir.join(relative).exists())
            .copied()
            .collect();
        if !collisions.is_empty() {
            return Err(format!(
                "workdir contains BioHub-managed paths without run.json: {}; choose a new workdir",
                collisions.join(", ")
            ));
        }
    }

    Ok(PreparedRun {
        recipe,
        workflow,
        config,
        workdir,
        profile,
        config_sha256,
        started_unix: started,
    })
}

fn initialize_run_files(
    options: &RunOptions,
    recipe: &RecipeSpec,
    workflow: &Path,
    config: &Path,
    profile: Option<&Path>,
    config_sha256: &str,
) -> RecipeResult<()> {
    for directory in ["logs", "results", "report"] {
        fs::create_dir_all(options.workdir.join(directory))
            .map_err(|error| format!("cannot create {directory} directory: {error}"))?;
    }
    let resolved = options.workdir.join("config.resolved.yaml");
    if !resolved.exists() {
        copy_create_new(config, &resolved)?;
    }
    write_atomic(
        &options.workdir.join("recipe.id"),
        &format!("{}\n", recipe.id),
    )?;
    write_atomic(
        &options.workdir.join("config.sha256"),
        &format!("{config_sha256}\n"),
    )?;
    let command = command_as_shell(
        workflow,
        &resolved,
        &options.workdir,
        profile,
        options.cores,
        options.dry_run,
    );
    write_atomic(&options.workdir.join("command.sh"), &command)?;
    let mut dependencies = recipe.dependencies.clone();
    dependencies.push("snakemake".to_string());
    dependencies.sort();
    dependencies.dedup();
    let mut versions = format!(
        "software\tversion\nbiohub\t{}\nsnakemake\t{}\n",
        env!("CARGO_PKG_VERSION"),
        snakemake_version()
    );
    for dependency in dependencies {
        if dependency != "snakemake" {
            versions.push_str(&format!(
                "{}\t{}\n",
                dependency.replace('\t', " "),
                dependency_version(&dependency)
            ));
        }
    }
    write_atomic(&options.workdir.join("versions.tsv"), &versions)?;
    write_recipe_source_checksums(workflow, &options.workdir.join("recipe.sources.sha256"))?;
    let workflow_sha256 = sha256_file(workflow)?;
    let runtime_container_digest = env::var("BIOHUB_CONTAINER_DIGEST")
        .ok()
        .filter(|value| !value.trim().is_empty())
        .map(|value| format!("\"{}\"", json_escape(value.trim())))
        .unwrap_or_else(|| "null".to_string());
    let provenance = format!(
        "{{\n  \"schema_version\": 1,\n  \"recipe_id\": \"{}\",\n  \"biohub_version\": \"{}\",\n  \"domain\": \"{}\",\n  \"status\": \"{}\",\n  \"config_sha256\": \"{}\",\n  \"workflow\": \"{}\",\n  \"workflow_sha256\": \"{}\",\n  \"recipe_sources_manifest\": \"recipe.sources.sha256\",\n  \"profile\": \"{}\",\n  \"container_hint\": \"{}\",\n  \"runtime_container_digest\": {}\n}}\n",
        json_escape(&recipe.id),
        env!("CARGO_PKG_VERSION"),
        json_escape(&recipe.domain),
        json_escape(&recipe.status),
        json_escape(config_sha256),
        json_escape(&workflow.to_string_lossy()),
        workflow_sha256,
        json_escape(&options.profile),
        json_escape(&recipe.container),
        runtime_container_digest
    );
    write_atomic(&options.workdir.join("provenance.json"), &provenance)?;
    Ok(())
}

fn collect_files(root: &Path, current: &Path, output: &mut Vec<PathBuf>) -> RecipeResult<()> {
    let mut entries: Vec<_> = fs::read_dir(current)
        .map_err(|error| format!("cannot list {}: {error}", current.display()))?
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|error| format!("cannot list {}: {error}", current.display()))?;
    entries.sort_by_key(|entry| entry.file_name());
    for entry in entries {
        let path = entry.path();
        let file_type = entry
            .file_type()
            .map_err(|error| format!("cannot inspect {}: {error}", path.display()))?;
        if file_type.is_symlink() {
            continue;
        }
        if file_type.is_dir() {
            collect_files(root, &path, output)?;
        } else if file_type.is_file() {
            let relative = path
                .strip_prefix(root)
                .map_err(|error| format!("cannot relativize {}: {error}", path.display()))?;
            if relative != Path::new("checksums.sha256")
                && relative != Path::new("run.json")
                && !relative.starts_with(".snakemake")
            {
                output.push(relative.to_path_buf());
            }
        }
    }
    Ok(())
}

fn write_checksums(workdir: &Path) -> RecipeResult<()> {
    let mut files = Vec::new();
    collect_files(workdir, workdir, &mut files)?;
    files.sort();
    let mut lines = String::new();
    for relative in files {
        let hash = sha256_file(&workdir.join(&relative))?;
        lines.push_str(&format!("{hash}  {}\n", relative.to_string_lossy()));
    }
    write_atomic(&workdir.join("checksums.sha256"), &lines)
}

fn write_recipe_source_checksums(workflow: &Path, destination: &Path) -> RecipeResult<()> {
    let recipe_dir = workflow
        .parent()
        .ok_or_else(|| format!("workflow has no parent directory: {}", workflow.display()))?;
    let recipes_root = recipe_dir
        .parent()
        .ok_or_else(|| format!("recipe has no package root: {}", recipe_dir.display()))?;
    let mut files = Vec::new();
    collect_files(recipes_root, recipe_dir, &mut files)?;
    let shared = recipes_root.join("_lib/provenance.py");
    if !shared.is_file() {
        return Err(format!(
            "missing shared recipe source: {}",
            shared.display()
        ));
    }
    files.push(PathBuf::from("_lib/provenance.py"));
    files.retain(|relative| {
        !relative
            .components()
            .any(|component| component.as_os_str() == "__pycache__")
            && relative.extension().and_then(|value| value.to_str()) != Some("pyc")
    });
    files.sort();
    files.dedup();
    let mut lines = String::new();
    for relative in files {
        let hash = sha256_file(&recipes_root.join(&relative))?;
        lines.push_str(&format!("{hash}  {}\n", relative.to_string_lossy()));
    }
    write_atomic(destination, &lines)
}

fn run_workflow(mut options: RunOptions) -> RecipeResult<i32> {
    if !command_available("snakemake") {
        return Err("missing dependency: snakemake".to_string());
    }
    let PreparedRun {
        recipe,
        workflow,
        config,
        workdir,
        profile,
        config_sha256,
        started_unix: started,
    } = prepare_run(&options)?;
    options.workdir = workdir;
    initialize_run_files(
        &options,
        &recipe,
        &workflow,
        &config,
        profile.as_deref(),
        &config_sha256,
    )?;
    let running_status = if options.dry_run {
        "validating"
    } else {
        "running"
    };
    write_run_state(
        &options.workdir,
        &RunState {
            recipe_id: &recipe.id,
            status: running_status,
            config_sha256: &config_sha256,
            started_unix: started,
            updated_unix: now_unix(),
            exit_code: None,
            message: "Snakemake started",
        },
    )?;

    let resolved = options.workdir.join("config.resolved.yaml");
    let status = match snakemake_command(
        &workflow,
        &resolved,
        &options.workdir,
        profile.as_deref(),
        options.cores,
        options.dry_run,
    )
    .env("BIOHUB_RUN_DIR", &options.workdir)
    .status()
    {
        Ok(status) => status,
        Err(error) => {
            write_run_state(
                &options.workdir,
                &RunState {
                    recipe_id: &recipe.id,
                    status: "failed",
                    config_sha256: &config_sha256,
                    started_unix: started,
                    updated_unix: now_unix(),
                    exit_code: None,
                    message: "failed to start Snakemake",
                },
            )?;
            return Err(format!("failed to start snakemake: {error}"));
        }
    };
    let exit_code = status.code().unwrap_or(1);
    if !status.success() {
        write_run_state(
            &options.workdir,
            &RunState {
                recipe_id: &recipe.id,
                status: "failed",
                config_sha256: &config_sha256,
                started_unix: started,
                updated_unix: now_unix(),
                exit_code: Some(exit_code),
                message: "Snakemake failed; inspect logs and resume after correction",
            },
        )?;
        return Ok(exit_code);
    }

    if !options.dry_run {
        let finalization = (|| {
            if !options.workdir.join("inputs.manifest.tsv").exists() {
                write_atomic(
                    &options.workdir.join("inputs.manifest.tsv"),
                    "logical_name\tpath\tsize_bytes\tsha256\n",
                )?;
            }
            write_checksums(&options.workdir)
        })();
        if let Err(error) = finalization {
            write_run_state(
                &options.workdir,
                &RunState {
                    recipe_id: &recipe.id,
                    status: "finalization_failed",
                    config_sha256: &config_sha256,
                    started_unix: started,
                    updated_unix: now_unix(),
                    exit_code: Some(exit_code),
                    message: "workflow finished but provenance finalization failed",
                },
            )?;
            return Err(error);
        }
    }
    write_run_state(
        &options.workdir,
        &RunState {
            recipe_id: &recipe.id,
            status: if options.dry_run {
                "validated"
            } else {
                "complete"
            },
            config_sha256: &config_sha256,
            started_unix: started,
            updated_unix: now_unix(),
            exit_code: Some(exit_code),
            message: "Snakemake finished",
        },
    )?;
    Ok(exit_code)
}

fn validate_recipe(recipe_id: &str, config: &Path, workdir: Option<&Path>) -> RecipeResult<i32> {
    let temporary = if workdir.is_none() {
        Some(unique_validation_directory(recipe_id)?)
    } else {
        None
    };
    let target = workdir
        .map(Path::to_path_buf)
        .or_else(|| temporary.clone())
        .ok_or_else(|| "cannot select validation directory".to_string())?;
    let options = RunOptions {
        recipe_id: recipe_id.to_string(),
        config: config.to_path_buf(),
        workdir: target.clone(),
        profile: "local".to_string(),
        cores: 1,
        dry_run: true,
        resume: false,
    };
    let result = run_workflow(options);
    if let Some(temporary) = temporary {
        let _ = fs::remove_dir_all(temporary);
    }
    result
}

fn parse_run_options(args: &[String]) -> RecipeResult<RunOptions> {
    if args.is_empty() {
        return Err("missing recipe-id".to_string());
    }
    let recipe_id = args[0].clone();
    reject_unknown_options(
        &args[1..],
        &["--config", "--workdir", "--profile", "--cores"],
        &["--dry-run", "--resume"],
    )?;
    let config = PathBuf::from(required_option(&args[1..], "--config")?);
    let workdir = PathBuf::from(required_option(&args[1..], "--workdir")?);
    let profile = option_value(&args[1..], "--profile")?.unwrap_or_else(|| "local".to_string());
    let raw_cores = option_value(&args[1..], "--cores")?.unwrap_or_else(|| "1".to_string());
    let cores = raw_cores
        .parse::<usize>()
        .map_err(|_| format!("invalid --cores value: {raw_cores}"))?;
    Ok(RunOptions {
        recipe_id,
        config,
        workdir,
        profile,
        cores,
        dry_run: contains_flag(&args[1..], "--dry-run"),
        resume: contains_flag(&args[1..], "--resume"),
    })
}

fn generate_report(workdir: &Path, force: bool) -> RecipeResult<i32> {
    if !command_available("snakemake") {
        return Err("missing dependency: snakemake".to_string());
    }
    let workdir = workdir
        .canonicalize()
        .map_err(|error| format!("cannot resolve workdir: {error}"))?;
    let recipe_id = read_trimmed(&workdir.join("recipe.id"))?;
    let recipe = recipe_by_id(&recipe_id)?;
    let root = recipe_root()?;
    let workflow = recipe_file(&root, &recipe.workflow)?;
    let config = workdir.join("config.resolved.yaml");
    if !config.is_file() {
        return Err(format!("missing resolved config: {}", config.display()));
    }
    fs::create_dir_all(workdir.join("report"))
        .map_err(|error| format!("cannot create report directory: {error}"))?;
    let report = workdir.join("report/report.html");
    if report.exists() && !force {
        return Err(format!(
            "report exists: {}; pass --force to replace",
            report.display()
        ));
    }
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    let temporary_report = workdir.join(format!(
        "report/.report.{}.{}.html",
        std::process::id(),
        nonce
    ));
    let status = Command::new("snakemake")
        .arg("--snakefile")
        .arg(workflow)
        .arg("--configfile")
        .arg(config)
        .arg("--directory")
        .arg(&workdir)
        .arg("--report")
        .arg(&temporary_report)
        .status()
        .map_err(|error| format!("failed to start snakemake report: {error}"))?;
    if status.success() {
        fs::rename(&temporary_report, &report)
            .map_err(|error| format!("cannot finalize {}: {error}", report.display()))?;
        write_checksums(&workdir)?;
    } else {
        let _ = fs::remove_file(&temporary_report);
    }
    Ok(status.code().unwrap_or(1))
}

fn path_is_executable(path: &Path) -> bool {
    let Ok(metadata) = path.metadata() else {
        return false;
    };
    if !metadata.is_file() {
        return false;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        metadata.permissions().mode() & 0o111 != 0
    }
    #[cfg(not(unix))]
    {
        true
    }
}

pub(crate) fn command_available(command: &str) -> bool {
    let Some(path_var) = env::var_os("PATH") else {
        return false;
    };
    env::split_paths(&path_var).any(|directory| {
        let candidate = directory.join(command);
        path_is_executable(&candidate)
            || cfg!(windows) && path_is_executable(&candidate.with_extension("exe"))
    })
}

pub(crate) fn dependency_available(dependency: &str) -> bool {
    let Some(package) = dependency.strip_prefix("R-package:") else {
        return command_available(dependency);
    };
    if package.is_empty() || !command_available("Rscript") {
        return false;
    }
    Command::new("Rscript")
        .arg("--vanilla")
        .arg("-e")
        .arg("quit(status=if (requireNamespace(commandArgs(TRUE)[1], quietly=TRUE)) 0 else 1)")
        .arg(package)
        .status()
        .is_ok_and(|status| status.success())
}

pub(crate) fn recipe_dependencies(id: &str) -> RecipeResult<Vec<String>> {
    let mut dependencies = vec!["snakemake".to_string()];
    dependencies.extend(recipe_by_id(id)?.dependencies);
    dependencies.sort();
    dependencies.dedup();
    Ok(dependencies)
}

fn print_help() {
    println!(
        "Usage:\n  biohub recipe list [--format table|json]\n  biohub recipe describe <recipe-id> [--format table|json]\n  biohub recipe init <recipe-id> --workdir DIR\n  biohub recipe validate <recipe-id> --config FILE [--workdir DIR]\n  biohub recipe run <recipe-id> --config FILE --workdir DIR [--profile local|slurm|PATH] [--cores N] [--dry-run] [--resume]\n  biohub recipe report --workdir DIR [--force]"
    );
}

pub(crate) fn run_recipe_cli(args: &[String]) -> i32 {
    let result = (|| -> RecipeResult<i32> {
        if args.is_empty() || matches!(args[0].as_str(), "--help" | "-h" | "help") {
            print_help();
            return Ok(0);
        }
        match args[0].as_str() {
            "list" => {
                reject_unknown_options(&args[1..], &["--format"], &[])?;
                match option_value(&args[1..], "--format")?
                    .unwrap_or_else(|| "table".to_string())
                    .as_str()
                {
                    "table" => print_list_table()?,
                    "json" => print_list_json()?,
                    other => return Err(format!("unsupported format: {other}")),
                }
                Ok(0)
            }
            "describe" => {
                if args.len() < 2 {
                    return Err("missing recipe-id".to_string());
                }
                reject_unknown_options(&args[2..], &["--format"], &[])?;
                let recipe = recipe_by_id(&args[1])?;
                let format =
                    option_value(&args[2..], "--format")?.unwrap_or_else(|| "table".to_string());
                if format != "table" && format != "json" {
                    return Err(format!("unsupported format: {format}"));
                }
                print_description(&recipe, format == "json");
                Ok(0)
            }
            "init" => {
                if args.len() < 2 {
                    return Err("missing recipe-id".to_string());
                }
                reject_unknown_options(&args[2..], &["--workdir"], &[])?;
                let recipe = recipe_by_id(&args[1])?;
                let workdir = PathBuf::from(required_option(&args[2..], "--workdir")?);
                init_recipe(&recipe, &workdir)?;
                println!("initialized {} in {}", recipe.id, workdir.display());
                Ok(0)
            }
            "validate" => {
                if args.len() < 2 {
                    return Err("missing recipe-id".to_string());
                }
                reject_unknown_options(&args[2..], &["--config", "--workdir"], &[])?;
                let config = PathBuf::from(required_option(&args[2..], "--config")?);
                let workdir = option_value(&args[2..], "--workdir")?.map(PathBuf::from);
                validate_recipe(&args[1], &config, workdir.as_deref())
            }
            "run" => run_workflow(parse_run_options(&args[1..])?),
            "report" => {
                reject_unknown_options(&args[1..], &["--workdir"], &["--force"])?;
                let workdir = PathBuf::from(required_option(&args[1..], "--workdir")?);
                generate_report(&workdir, contains_flag(&args[1..], "--force"))
            }
            other => Err(format!("unknown recipe command: {other}")),
        }
    })();
    match result {
        Ok(code) => code,
        Err(error) => {
            eprintln!("{error}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalog_is_well_formed_and_unique() {
        let recipes = load_recipe_catalog().expect("load recipe catalog");
        assert!(!recipes.is_empty());
        let ids: std::collections::HashSet<_> = recipes.iter().map(|item| &item.id).collect();
        assert_eq!(ids.len(), recipes.len());
        assert!(recipes
            .iter()
            .all(|item| item.dependencies.contains(&"snakemake".to_string())));
    }

    #[test]
    fn relative_paths_reject_parent_escape() {
        assert!(validate_relative_path("comparative/Snakefile").is_ok());
        assert!(validate_relative_path("../Snakefile").is_err());
        assert!(validate_relative_path("/tmp/Snakefile").is_err());
    }

    #[test]
    fn shell_quoting_handles_spaces_and_quotes() {
        assert_eq!(shell_quote("a b"), "'a b'");
        assert_eq!(shell_quote("a'b"), "'a'\"'\"'b'");
    }
}
