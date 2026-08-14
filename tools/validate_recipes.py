#!/usr/bin/env python3
"""Static release contract checks for packaged BioHub recipes."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from pathlib import Path
import re
import shutil
import subprocess


REQUIRED_FILES = {"Snakefile", "config.example.yaml", "config.schema.yaml"}
RECIPE_GUIDE_SECTIONS = (
    "### 适用场景与边界",
    "### 输入准备",
    "### 配置字段",
    "### 运行步骤",
    "### 产物与解读",
    "### 失败、恢复与发表核对",
)
DISALLOWED_R = (
    "install.packages(",
    "install_github(",
    "file.choose(",
    "setwd(",
    'source("http',
    "source('http",
    "/Users/",
)


def parse_options() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=Path(__file__).resolve().parents[1])
    parser.add_argument("--require-r", action="store_true")
    return parser.parse_args()


def read_catalog(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("#"):
        raise ValueError("recipe catalog needs commented header")
    rows = list(csv.DictReader([lines[0][1:], *lines[1:]], delimiter="\t"))
    expected = {
        "id",
        "domain",
        "description",
        "status",
        "workflow",
        "config_template",
        "config_schema",
        "dependencies",
        "container",
        "license",
    }
    if not rows or set(rows[0]) != expected:
        raise ValueError("recipe catalog header/rows do not match contract")
    ids = [row["id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate recipe IDs")
    return rows


def example_config_keys(path: Path) -> set[str]:
    return {
        match.group(1)
        for line in path.read_text(encoding="utf-8").splitlines()
        if (match := re.match(r"^([A-Za-z0-9_]+):", line))
    }


def yaml_contract(path: Path, example: Path) -> None:
    schema_lines = path.read_text(encoding="utf-8").splitlines()
    properties = set()
    required = set()
    in_required = False
    in_properties = False
    for line in schema_lines:
        if line == "required:":
            in_required = True
            in_properties = False
            continue
        if line == "properties:":
            in_required = False
            in_properties = True
            continue
        if line and not line.startswith(" "):
            in_required = False
            in_properties = False
        if in_required:
            match = re.fullmatch(r"  - ([A-Za-z0-9_]+)", line)
            if match:
                required.add(match.group(1))
        if in_properties:
            match = re.match(r"^  ([A-Za-z0-9_]+):", line)
            if match:
                properties.add(match.group(1))
    example_keys = example_config_keys(example)
    if not properties or required != properties:
        raise ValueError(f"all schema properties must be required: {path}")
    if example_keys != properties:
        raise ValueError(f"example/schema key mismatch: {example}")


def validate_local_markdown_links(paths: list[Path]) -> None:
    link_pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
    for source in paths:
        text = source.read_text(encoding="utf-8")
        for raw_target in link_pattern.findall(text):
            target = raw_target.strip().split(" ", 1)[0].strip("<>")
            if not target or target.startswith(("#", "http://", "https://", "mailto:")):
                continue
            relative = target.split("#", 1)[0]
            if relative and not (source.parent / relative).resolve().exists():
                raise ValueError(f"broken local Markdown link in {source}: {target}")


def validate_recipe_guide(root: Path, catalog: list[dict[str, str]]) -> None:
    guide_path = root / "docs/RECIPES.zh-CN.md"
    guide = guide_path.read_text(encoding="utf-8")
    if "文档版本：0.4.0" not in guide:
        raise ValueError("recipe guide version does not match release")
    for private_marker in ("/Users/", "/mnt/kobe/"):
        if private_marker in guide:
            raise ValueError(f"private path in recipe guide: {private_marker}")

    for row in catalog:
        marker = f"## `{row['id']}`"
        if guide.splitlines().count(marker) != 1:
            raise ValueError(f"recipe guide needs exactly one section: {row['id']}")
        start = guide.index(marker) + len(marker)
        next_heading = guide.find("\n## ", start)
        section = guide[start:] if next_heading == -1 else guide[start:next_heading]
        missing_sections = [name for name in RECIPE_GUIDE_SECTIONS if name not in section]
        if missing_sections:
            raise ValueError(
                f"recipe guide section {row['id']} lacks {missing_sections}"
            )
        config = root / "recipes" / row["config_template"]
        missing_keys = sorted(
            key for key in example_config_keys(config) if f"`{key}`" not in section
        )
        if missing_keys:
            raise ValueError(
                f"recipe guide section {row['id']} lacks config keys {missing_keys}"
            )

    readme = (root / "README.md").read_text(encoding="utf-8")
    command_guide = (root / "docs/USER_GUIDE.zh-CN.md").read_text(encoding="utf-8")
    if "docs/RECIPES.zh-CN.md" not in readme:
        raise ValueError("README does not link recipe guide")
    if "RECIPES.zh-CN.md" not in command_guide:
        raise ValueError("command guide does not link recipe guide")
    validate_local_markdown_links(
        [root / "README.md", root / "docs/USER_GUIDE.zh-CN.md", guide_path]
    )


def validate_migration_matrix(path: Path) -> int:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("empty migration matrix")
    expected_columns = {
        "inventory_id",
        "source_artifact",
        "language",
        "decision",
        "target",
        "verification",
        "note",
    }
    if set(rows[0]) != expected_columns:
        raise ValueError("migration matrix columns do not match contract")
    expected_ids = [f"{index:03d}" for index in range(1, len(rows) + 1)]
    if [row["inventory_id"] for row in rows] != expected_ids:
        raise ValueError("migration matrix IDs must be sequential")
    if len({row["source_artifact"] for row in rows}) != len(rows):
        raise ValueError("migration matrix source_artifact values must be unique")
    allowed = {"integrated", "recipe", "superseded", "retained", "deferred", "excluded"}
    for row in rows:
        if row["decision"] not in allowed:
            raise ValueError(f"invalid migration decision: {row['decision']}")
        if not row["source_artifact"] or Path(row["source_artifact"]).is_absolute():
            raise ValueError("migration source paths must be non-empty and sanitized")
        if row["decision"] in {"deferred", "excluded"}:
            if row["target"] != "-":
                raise ValueError("deferred/excluded migration rows must use target '-'")
        elif row["target"] == "-":
            raise ValueError("resolved migration rows need a target")
        if not row["verification"] or not row["note"]:
            raise ValueError("migration rows need verification and note")
    expected_decisions = {
        "integrated": 63,
        "recipe": 5,
        "superseded": 11,
        "retained": 2,
        "deferred": 5,
        "excluded": 5,
    }
    if Counter(row["decision"] for row in rows) != expected_decisions:
        raise ValueError("migration decision counts differ from dated snapshot")
    return len(rows)


def main() -> int:
    options = parse_options()
    root = Path(options.root).resolve()
    recipes = root / "recipes"
    catalog = read_catalog(root / "biohub-rs/src/recipe_catalog.tsv")
    catalog_ids = {row["id"] for row in catalog}
    directory_ids = {
        path.name
        for path in recipes.iterdir()
        if path.is_dir() and not path.name.startswith("_") and path.name != "profiles"
    }
    if directory_ids != catalog_ids:
        raise ValueError(
            f"recipe directories differ from catalog: {directory_ids ^ catalog_ids}"
        )

    python_count = 0
    r_count = 0
    for row in catalog:
        dependencies = row["dependencies"].split(",")
        if any(not value for value in dependencies) or len(dependencies) != len(set(dependencies)):
            raise ValueError(f"invalid or duplicate dependencies: {row['id']}")
        if not {"snakemake", "python3"}.issubset(dependencies):
            raise ValueError(f"recipe missing runtime dependency: {row['id']}")
        if row["status"] != "experimental" or row["license"] != "MIT":
            raise ValueError(f"unexpected recipe status/license: {row['id']}")
        expected_paths = {
            "workflow": f"{row['id']}/Snakefile",
            "config_template": f"{row['id']}/config.example.yaml",
            "config_schema": f"{row['id']}/config.schema.yaml",
        }
        if any(row[field] != value for field, value in expected_paths.items()):
            raise ValueError(f"recipe paths do not follow package layout: {row['id']}")
        if row["id"] == "kmer-gwas":
            if row["container"] != "not-provided-python2-eol":
                raise ValueError("kmer-gwas must not advertise a Python 2 image")
        elif not row["container"].startswith("ghcr.io/wzhennan-icmm/biohub-"):
            raise ValueError(f"invalid domain container: {row['id']}")
        directory = recipes / row["id"]
        missing = REQUIRED_FILES - {
            path.name for path in directory.iterdir() if path.is_file()
        }
        if missing:
            raise ValueError(
                f"missing recipe files for {row['id']}: {sorted(missing)}"
            )
        for field in ("workflow", "config_template", "config_schema"):
            candidate = recipes / row[field]
            if not candidate.is_file() or recipes not in candidate.resolve().parents:
                raise ValueError(f"invalid catalog path: {row[field]}")
        snakefile = (directory / "Snakefile").read_text(encoding="utf-8")
        for marker in (
            "validate(config",
            "write_input_manifest",
            "write_readme",
            "rule archive",
        ):
            if marker not in snakefile:
                raise ValueError(f"{row['id']} Snakefile missing {marker}")
        yaml_contract(
            directory / "config.schema.yaml", directory / "config.example.yaml"
        )

        for script in sorted(directory.glob("*.py")):
            compile(script.read_text(encoding="utf-8"), str(script), "exec")
            python_count += 1
        for script in sorted(directory.glob("*.R")):
            source = script.read_text(encoding="utf-8")
            for marker in DISALLOWED_R:
                if marker in source:
                    raise ValueError(
                        f"disallowed R behavior {marker!r}: {script}"
                    )
            r_count += 1

    provenance = recipes / "_lib/provenance.py"
    compile(provenance.read_text(encoding="utf-8"), str(provenance), "exec")
    python_count += 1

    rscript = shutil.which("Rscript")
    if options.require_r and not rscript:
        raise ValueError("Rscript is required")
    if rscript:
        scripts = sorted(root.glob("r/*.R")) + sorted(recipes.glob("*/*.R"))
        for script in scripts:
            subprocess.run(
                [
                    rscript,
                    "-e",
                    "parse(file=commandArgs(TRUE)[1])",
                    str(script),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
            )

    profile = recipes / "profiles/slurm/config.yaml"
    if (
        not profile.is_file()
        or "executor: slurm" not in profile.read_text(encoding="utf-8")
    ):
        raise ValueError("missing Slurm executor profile")
    migration_count = validate_migration_matrix(
        root / "docs/SCRIPT_MIGRATION_MATRIX.tsv"
    )
    validate_recipe_guide(root, catalog)
    print(
        f"recipes={len(catalog)} python={python_count} "
        f"r={r_count} migration_rows={migration_count} recipe_guide=ok"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
