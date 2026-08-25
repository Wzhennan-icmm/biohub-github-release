#!/usr/bin/env python3
"""Build and verify BioHub domain-review evidence without granting approval."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Any


PACK_IDS = (
    "annotation-coordinates",
    "orthology-codon",
    "visualization",
    "statistics",
)
PACK_SCHEMA_VERSION = 1
GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


class ValidationError(ValueError):
    """Raised when evidence or pack contract is invalid."""


def options() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build reproducible review evidence; never approve reviews",
    )
    parser.add_argument("--root", default=Path(__file__).resolve().parents[1])
    subparsers = parser.add_subparsers(dest="action", required=True)

    build = subparsers.add_parser("build", help="execute packs and create evidence")
    build.add_argument("--pack", choices=(*PACK_IDS, "all"), default="all")
    build.add_argument("--output", help="evidence root; defaults to validation/evidence")
    build.add_argument("--biohub", help="BioHub executable")
    build.add_argument("--snakemake", help="Snakemake executable used by recipe packs")
    build.add_argument(
        "--force",
        action="store_true",
        help="replace only selected existing evidence directories",
    )
    build.add_argument(
        "--keep-failed",
        action="store_true",
        help="retain selected failed pack output for diagnosis; status remains failed",
    )

    verify = subparsers.add_parser("verify", help="verify hashes and comparisons")
    verify.add_argument("--pack", choices=(*PACK_IDS, "all"), default="all")
    verify.add_argument("--output", help="evidence root; defaults to validation/evidence")

    summary = subparsers.add_parser("summary", help="show pack and review status")
    summary.add_argument("--format", choices=("text", "json"), default="text")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def files_under(root: Path) -> list[Path]:
    return sorted(path for path in root.rglob("*") if path.is_file())


def write_hash_manifest(base: Path, paths: list[Path], output: Path) -> None:
    rows = []
    for path in sorted(paths):
        resolved = path.resolve()
        try:
            relative = resolved.relative_to(base.resolve())
        except ValueError as error:
            raise ValidationError(f"manifest path escapes evidence: {path}") from error
        rows.append(f"{sha256_file(resolved)}  {relative.as_posix()}")
    output.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")


def verify_hash_manifest(base: Path, manifest: Path) -> int:
    if not manifest.is_file():
        raise ValidationError(f"missing hash manifest: {manifest}")
    count = 0
    for line_number, row in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        expected, separator, relative = row.partition("  ")
        if not separator or not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise ValidationError(f"invalid hash row {manifest}:{line_number}")
        path = (base / relative).resolve()
        try:
            path.relative_to(base.resolve())
        except ValueError as error:
            raise ValidationError(f"hash path escapes evidence: {relative}") from error
        if not path.is_file():
            raise ValidationError(f"hashed file missing: {relative}")
        actual = sha256_file(path)
        if actual != expected:
            raise ValidationError(f"SHA256 mismatch: {relative}")
        count += 1
    return count


def load_pack(pack_dir: Path, expected_pack_id: str | None = None) -> dict[str, Any]:
    manifest_path = pack_dir / "pack.json"
    try:
        pack = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValidationError(f"cannot read pack manifest: {manifest_path}: {error}") from error
    required = {
        "schema_version",
        "pack_id",
        "title",
        "inventory_ids",
        "fixture_license",
        "acceptance",
        "manual_checks",
        "reference_commands",
        "steps",
    }
    missing = required - set(pack)
    if missing:
        raise ValidationError(f"pack {pack_dir.name} missing fields: {sorted(missing)}")
    if pack["schema_version"] != PACK_SCHEMA_VERSION:
        raise ValidationError(f"unsupported pack schema: {pack['schema_version']}")
    expected = expected_pack_id or pack_dir.name
    if pack["pack_id"] != expected or pack["pack_id"] not in PACK_IDS:
        raise ValidationError(f"pack ID/path mismatch: {pack['pack_id']} != {expected}")
    if not pack["inventory_ids"] or not all(
        re.fullmatch(r"\d{3}", item) for item in pack["inventory_ids"]
    ):
        raise ValidationError(f"invalid inventory IDs in {pack['pack_id']}")
    step_ids = [step.get("id") for step in pack["steps"]]
    if any(not value or not re.fullmatch(r"[A-Za-z0-9_.-]+", value) for value in step_ids):
        raise ValidationError(f"invalid step ID in {pack['pack_id']}")
    if len(step_ids) != len(set(step_ids)):
        raise ValidationError(f"duplicate step ID in {pack['pack_id']}")
    return pack


def render(value: Any, context: dict[str, str]) -> Any:
    if isinstance(value, str):
        try:
            return value.format_map(context)
        except KeyError as error:
            raise ValidationError(f"unknown pack placeholder: {error.args[0]}") from error
    if isinstance(value, list):
        return [render(item, context) for item in value]
    if isinstance(value, dict):
        return {key: render(item, context) for key, item in value.items()}
    return value


def resolve_executable(root: Path, requested: str | None) -> Path:
    candidates: list[Path] = []
    if requested:
        requested_path = Path(requested).expanduser()
        found = shutil.which(requested) if requested_path.parent == Path(".") else None
        candidates.append(Path(found) if found else requested_path)
    else:
        candidates.extend(
            [root / "biohub-rs/target/release/biohub", root / "biohub-rs/target/debug/biohub"]
        )
        found = shutil.which("biohub")
        if found:
            candidates.append(Path(found))
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise ValidationError("BioHub executable not found; pass --biohub")


def resolve_snakemake(root: Path, requested: str | None) -> Path | None:
    values = [requested] if requested else []
    values.extend([str(root / ".venv-ci/bin/snakemake"), "snakemake"])
    for value in values:
        if not value:
            continue
        path = Path(value).expanduser()
        found = shutil.which(value) if path.parent == Path(".") else None
        candidate = Path(found) if found else path
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    return None


def command_version(argv: list[str], environment: dict[str, str]) -> str:
    try:
        result = subprocess.run(
            argv,
            text=True,
            capture_output=True,
            check=False,
            env=environment,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return f"unavailable: {error}"
    combined = (result.stdout + "\n" + result.stderr).strip()
    first = combined.splitlines()[0] if combined else f"exit={result.returncode}"
    return first.replace("\t", " ")


def git_state(root: Path) -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "-c", f"safe.directory={root}", "rev-parse", "HEAD"],
        cwd=root,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    dirty = bool(
        subprocess.run(
            ["git", "-c", f"safe.directory={root}", "status", "--porcelain"],
            cwd=root,
            text=True,
            capture_output=True,
            check=True,
        ).stdout.strip()
    )
    return commit, dirty


def run_step(
    root: Path,
    evidence: Path,
    source: Path,
    step: dict[str, Any],
    biohub: Path,
    snakemake: Path | None,
    environment: dict[str, str],
) -> tuple[list[str], Path]:
    step_base = evidence / "steps" / step["id"]
    step_base.mkdir(parents=True)
    runner = step.get("runner")
    output = step_base
    context = {
        "root": str(root),
        "source": str(source),
        "fixtures": str(source / "fixtures"),
        "expected": str(source / "expected"),
        "evidence": str(evidence),
        "steps": str(evidence / "steps"),
        "output": str(output),
        "biohub": str(biohub),
    }
    if runner == "biohub":
        command = [str(biohub), *map(str, render(step.get("args", []), context))]
        cwd = step_base
    elif runner == "recipe":
        if snakemake is None:
            raise ValidationError(f"step {step['id']} needs Snakemake; pass --snakemake")
        run_dir = step_base / "run"
        output = run_dir
        context["output"] = str(output)
        config_source = source / step["config"]
        config = json.loads(config_source.read_text(encoding="utf-8"))
        config = render(config, context)
        config_path = step_base / "config.resolved.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        command = [
            str(biohub),
            "recipe",
            "run",
            step["recipe"],
            "--config",
            str(config_path),
            "--workdir",
            str(run_dir),
            "--cores",
            "1",
        ]
        cwd = root
    else:
        raise ValidationError(f"unsupported runner {runner!r} in step {step['id']}")

    step_environment = environment.copy()
    step_environment["BIOHUB_RECIPE_DIR"] = str(root / "recipes")
    step_environment["XDG_CACHE_HOME"] = str(evidence / ".cache")
    if snakemake:
        step_environment["PATH"] = str(snakemake.parent) + os.pathsep + step_environment["PATH"]
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
        env=step_environment,
    )
    (step_base / "stdout.log").write_text(result.stdout, encoding="utf-8")
    (step_base / "stderr.log").write_text(result.stderr, encoding="utf-8")
    (step_base / "command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
    if result.returncode != 0:
        diagnostic = (result.stderr.strip() or result.stdout.strip() or "no subprocess output")[-4000:]
        raise ValidationError(
            f"step {step['id']} failed with exit {result.returncode}:\n{diagnostic}"
        )
    declared = [output / path for path in step.get("outputs", [])]
    missing = [path for path in declared if not path.is_file()]
    if missing:
        raise ValidationError(f"step {step['id']} missing outputs: {missing}")
    return [str(value) for value in command], output


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def parse_number(value: str, label: str) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        return math.nan
    try:
        return float(value)
    except ValueError as error:
        raise ValidationError(f"non-numeric value for {label}: {value!r}") from error


def numeric_close(actual: float, expected: float, absolute: float, relative: float) -> bool:
    if math.isnan(actual) or math.isnan(expected):
        return math.isnan(actual) and math.isnan(expected)
    return abs(actual - expected) <= max(absolute, relative * abs(expected))


def keyed_rows(
    path: Path, key_columns: list[str]
) -> tuple[list[str], dict[tuple[str, ...], dict[str, str]]]:
    header, rows = read_tsv(path)
    missing = set(key_columns) - set(header)
    if missing:
        raise ValidationError(f"{path} lacks key columns: {sorted(missing)}")
    keyed: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = tuple(row[column] for column in key_columns)
        if key in keyed:
            raise ValidationError(f"duplicate TSV key {key} in {path}")
        keyed[key] = row
    return header, keyed


def fnv1a64(data: bytes) -> int:
    value = 0xCBF29CE484222325
    for byte in data:
        value ^= byte
        value = (value * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return value


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current: str | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            if current in records:
                raise ValidationError(f"duplicate FASTA ID {current} in {path}")
            records[current] = ""
        elif current is None:
            raise ValidationError(f"sequence before FASTA header in {path}")
        else:
            records[current] += line.upper()
    return records


def read_paml(path: Path) -> dict[str, str]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines or len(lines[0].split()) < 2:
        raise ValidationError(f"invalid PAML header: {path}")
    try:
        expected_records, expected_length = map(int, lines[0].split()[:2])
    except ValueError as error:
        raise ValidationError(f"invalid PAML dimensions: {path}") from error
    records: dict[str, str] = {}
    index = 1
    while index < len(lines):
        fields = lines[index].split()
        if len(fields) >= 2:
            name, sequence = fields[0], "".join(fields[1:])
            index += 1
        else:
            name = fields[0]
            index += 1
            if index >= len(lines):
                raise ValidationError(f"missing PAML sequence for {name}")
            sequence = "".join(lines[index].split())
            index += 1
        records[name] = sequence.upper()
    if len(records) != expected_records:
        raise ValidationError(f"PAML record count differs: {len(records)} != {expected_records}")
    if any(len(sequence) != expected_length for sequence in records.values()):
        raise ValidationError(f"PAML sequence length differs from header in {path}")
    return records


def translate(sequence: str) -> str:
    if len(sequence) % 3:
        raise ValidationError(f"codon sequence length is not divisible by 3: {len(sequence)}")
    amino = []
    for index in range(0, len(sequence), 3):
        codon = sequence[index : index + 3]
        if "-" in codon:
            amino.append("-")
        else:
            amino.append(GENETIC_CODE.get(codon, "X"))
    return "".join(amino)


def compare_check(source: Path, output: Path, check: dict[str, Any]) -> str:
    kind = check.get("type")
    actual = output / check["path"]
    if not actual.is_file():
        raise ValidationError(f"comparison output missing: {actual}")
    if kind == "exact":
        expected = source / check["expected"]
        if actual.read_bytes() != expected.read_bytes():
            raise ValidationError(f"exact comparison failed: {actual} != {expected}")
        return f"exact bytes match {expected.relative_to(source)}"

    if kind == "tsv_numeric":
        expected = source / check["expected"]
        keys = check["key_columns"]
        actual_header, actual_rows = keyed_rows(actual, keys)
        expected_header, expected_rows = keyed_rows(expected, keys)
        if actual_header != expected_header:
            raise ValidationError(f"TSV headers differ: {actual} != {expected}")
        if set(actual_rows) != set(expected_rows):
            raise ValidationError(f"TSV key sets differ: {actual} != {expected}")
        exact_columns = check.get("exact_columns", [])
        numeric_columns = check.get("numeric_columns", [])
        absolute = float(check.get("absolute_tolerance", 0.0))
        relative = float(check.get("relative_tolerance", 0.0))
        tolerances = check.get("tolerances", {})
        for key in sorted(actual_rows):
            for column in exact_columns:
                if actual_rows[key][column] != expected_rows[key][column]:
                    raise ValidationError(f"TSV exact value differs at {key}/{column}")
            for column in numeric_columns:
                observed = parse_number(actual_rows[key][column], f"{key}/{column}")
                reference = parse_number(expected_rows[key][column], f"{key}/{column}")
                column_tolerance = tolerances.get(column, {})
                column_absolute = float(column_tolerance.get("absolute", absolute))
                column_relative = float(column_tolerance.get("relative", relative))
                if not numeric_close(observed, reference, column_absolute, column_relative):
                    raise ValidationError(
                        f"TSV numeric value differs at {key}/{column}: "
                        f"{observed} != {reference}"
                    )
        detail = (
            f"{len(actual_rows)} keyed rows match; tolerance="
            f"max({absolute}, {relative}*abs(reference))"
        )
        if tolerances:
            detail += "; column_overrides=" + json.dumps(tolerances, sort_keys=True)
        return detail

    if kind == "tsv_sign_invariant":
        expected = source / check["expected"]
        keys = check["key_columns"]
        _, actual_rows = keyed_rows(actual, keys)
        _, expected_rows = keyed_rows(expected, keys)
        if set(actual_rows) != set(expected_rows):
            raise ValidationError(f"score-table key sets differ: {actual} != {expected}")
        threshold = float(check.get("minimum_absolute_correlation", 0.999999))
        for column in check["numeric_columns"]:
            observed_scores = [
                parse_number(actual_rows[key][column], column) for key in sorted(actual_rows)
            ]
            reference_scores = [
                parse_number(expected_rows[key][column], column) for key in sorted(expected_rows)
            ]
            observed_mean = sum(observed_scores) / len(observed_scores)
            reference_mean = sum(reference_scores) / len(reference_scores)
            numerator = sum(
                (left - observed_mean) * (right - reference_mean)
                for left, right in zip(observed_scores, reference_scores)
            )
            denominator = math.sqrt(
                sum((value - observed_mean) ** 2 for value in observed_scores)
                * sum((value - reference_mean) ** 2 for value in reference_scores)
            )
            if denominator == 0:
                correlation = 1.0 if observed_scores == reference_scores else 0.0
            else:
                correlation = numerator / denominator
            if abs(correlation) < threshold:
                raise ValidationError(
                    f"sign-invariant correlation failed for {actual}/{column}: {correlation}"
                )
        return f"axis scores have abs(correlation) >= {threshold}"

    if kind == "svg":
        data = actual.read_bytes()
        text = data.decode("utf-8")
        for fragment in check.get("contains", []):
            if fragment not in text:
                raise ValidationError(f"SVG lacks required fragment {fragment!r}: {actual}")
        if "NaN" in text:
            raise ValidationError(f"SVG contains NaN: {actual}")
        expected_circles = check.get("circle_count")
        if expected_circles is not None and text.count("<circle ") != expected_circles:
            raise ValidationError(f"SVG circle count differs: {actual}")
        expected_fingerprint = check.get("fnv1a64")
        if expected_fingerprint is not None and fnv1a64(data) != int(expected_fingerprint):
            raise ValidationError(f"SVG FNV-1a fingerprint differs: {actual}")
        return f"SVG semantics and fingerprint pass; bytes={len(data)}"

    if kind == "codon_translation":
        proteins = read_fasta(source / check["protein_fasta"])
        codons = read_paml(actual)
        if set(proteins) != set(codons):
            raise ValidationError("protein/codon ID sets differ")
        for identifier in sorted(proteins):
            sequence = codons[identifier]
            amino = translate(sequence)
            if "*" in amino[:-1]:
                raise ValidationError(f"unexpected internal stop for {identifier}")
            if amino.rstrip("*").replace("-", "") != proteins[identifier].rstrip("*"):
                raise ValidationError(f"translation differs for {identifier}")
        return f"{len(codons)} codon sequences translate to source proteins"

    raise ValidationError(f"unsupported comparison type: {kind!r}")


def compare_pack(evidence: Path, pack: dict[str, Any]) -> list[tuple[str, str, str]]:
    source = evidence / "source"
    rows: list[tuple[str, str, str]] = []
    for step in pack["steps"]:
        output = evidence / "steps" / step["id"]
        if step["runner"] == "recipe":
            output /= "run"
        for index, check in enumerate(step.get("checks", []), 1):
            detail = compare_check(source, output, check)
            rows.append((step["id"], f"{check['type']}:{index}", detail))
    return rows


def review_document(
    pack: dict[str, Any], commit: str, dirty: bool, generated_at: str, evidence: Path
) -> str:
    input_manifest = sha256_file(evidence / "inputs.sha256")
    output_manifest = sha256_file(evidence / "outputs.sha256")
    acceptance = "\n".join(f"- {item}" for item in pack["acceptance"])
    manual = "\n".join(f"- [ ] {item}" for item in pack["manual_checks"])
    inventory = ", ".join(pack["inventory_ids"])
    return f"""# BioHub validation review: {pack['pack_id']}

- Title: {pack['title']}
- Inventory IDs: {inventory}
- BioHub commit: `{commit}`
- Worktree dirty during build: `{str(dirty).lower()}`
- Evidence generated at: `{generated_at}`
- Fixture license: `{pack['fixture_license']}`
- Input-manifest SHA256: `{input_manifest}`
- Output-manifest SHA256: `{output_manifest}`
- Executed commands: `commands.tsv`
- Reference commands: `reference-commands.txt`
- Software versions: `versions.tsv`
- Automated differences: `comparisons.tsv`
- Automated comparison: `passed`
- Human decision: `pending`

## Acceptance contract

{acceptance}

## Manual review

{manual}

## Human sign-off

- Reviewer name:
- Reviewer affiliation:
- Review date (`YYYY-MM-DD`):
- Decision (`approved` or `rejected`):
- Differences found and resolution:
- Data/license confirmation:

This file is evidence template, not approval. Only explicit human decision may update
`validation/reviews.tsv` and migration matrix.
"""


def write_gallery(pack: dict[str, Any], evidence: Path) -> None:
    items = pack.get("gallery", [])
    if not items:
        return
    figures = []
    for item in items:
        relative = Path("steps") / item["step"] / item["path"]
        target = evidence / relative
        if not target.is_file():
            raise ValidationError(f"gallery image missing: {target}")
        figures.append(
            "<figure><img src=\"{}\" alt=\"{}\"><figcaption>{}</figcaption></figure>".format(
                html.escape(relative.as_posix(), quote=True),
                html.escape(item["caption"], quote=True),
                html.escape(item["caption"]),
            )
        )
    document = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>BioHub visual review</title>
<style>body{{font:16px system-ui;margin:2rem;max-width:1100px}}figure{{margin:0 0 3rem}}
img{{border:1px solid #bbb;max-width:100%;height:auto}}figcaption{{font-weight:600;margin-top:.5rem}}</style>
</head><body><h1>BioHub visual review</h1><p>Human decision: pending.</p>
{}\n</body></html>\n""".format("\n".join(figures))
    (evidence / "gallery.html").write_text(document, encoding="utf-8")


def build_pack(
    root: Path,
    pack_id: str,
    output_root: Path,
    biohub: Path,
    snakemake: Path | None,
    force: bool,
    keep_failed: bool,
) -> Path:
    source_pack = root / "validation/packs" / pack_id
    load_pack(source_pack)
    output_root.mkdir(parents=True, exist_ok=True)
    final = output_root / pack_id
    if final.exists() and not force:
        raise ValidationError(f"evidence exists: {final}; use --force to replace selected pack")
    temporary = Path(tempfile.mkdtemp(prefix=f".{pack_id}.", dir=output_root))
    try:
        source = temporary / "source"
        shutil.copytree(source_pack, source)
        copied_pack = load_pack(source, pack_id)
        commit, dirty = git_state(root)
        generated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
        environment = os.environ.copy()
        probe_environment = environment.copy()
        if snakemake:
            probe_environment["PATH"] = (
                str(snakemake.parent) + os.pathsep + probe_environment["PATH"]
            )
        command_rows: list[tuple[str, str]] = []
        for step in copied_pack["steps"]:
            command, _ = run_step(
                root, temporary, source, step, biohub, snakemake, environment
            )
            command_rows.append((step["id"], shlex.join(command)))

        comparison_rows = compare_pack(temporary, copied_pack)
        write_gallery(copied_pack, temporary)
        with (temporary / "comparisons.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["step", "check", "status", "detail"])
            for step_id, check_id, detail in comparison_rows:
                writer.writerow([step_id, check_id, "passed", detail])
        with (temporary / "commands.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["step", "command"])
            writer.writerows(command_rows)
        (temporary / "reference-commands.txt").write_text(
            "\n".join(copied_pack["reference_commands"]) + "\n",
            encoding="utf-8",
        )

        version_rows = [
            ("biohub", str(biohub), command_version([str(biohub), "--version"], probe_environment)),
            ("python", sys.executable, sys.version.splitlines()[0]),
        ]
        if snakemake:
            version_rows.append(
                ("snakemake", str(snakemake), command_version([str(snakemake), "--version"], probe_environment))
            )
        for dependency in copied_pack.get("dependencies", []):
            found = shutil.which(dependency, path=probe_environment.get("PATH"))
            if found:
                version_rows.append(
                    (dependency, found, command_version([found, "--version"], probe_environment))
                )
            else:
                version_rows.append((dependency, "-", "not found"))
        rscript = shutil.which("Rscript", path=probe_environment.get("PATH"))
        if copied_pack.get("r_packages") and rscript:
            version_rows.append(
                ("R", rscript, command_version([rscript, "--version"], probe_environment))
            )
            for package in copied_pack["r_packages"]:
                expression = f"cat(as.character(packageVersion({json.dumps(package)})))"
                version_rows.append(
                    (
                        f"R-package:{package}",
                        rscript,
                        command_version([rscript, "-e", expression], probe_environment),
                    )
                )
        with (temporary / "versions.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["component", "path", "version"])
            writer.writerows(version_rows)

        input_files = files_under(source / "fixtures") if (source / "fixtures").is_dir() else []
        write_hash_manifest(temporary, input_files, temporary / "inputs.sha256")
        output_files = files_under(temporary / "steps")
        write_hash_manifest(temporary, output_files, temporary / "outputs.sha256")
        write_hash_manifest(temporary, files_under(source), temporary / "source.sha256")
        (temporary / "review.md").write_text(
            review_document(copied_pack, commit, dirty, generated_at, temporary),
            encoding="utf-8",
        )
        metadata = {
            "schema_version": 1,
            "pack_id": pack_id,
            "inventory_ids": copied_pack["inventory_ids"],
            "biohub_commit": commit,
            "worktree_dirty": dirty,
            "generated_at": generated_at,
            "automated_status": "passed",
            "human_status": "pending",
        }
        (temporary / "evidence.json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )
        support_files = [
            temporary / name
            for name in (
                "commands.tsv",
                "comparisons.tsv",
                "evidence.json",
                "inputs.sha256",
                "outputs.sha256",
                "reference-commands.txt",
                "review.md",
                "source.sha256",
                "versions.tsv",
            )
        ]
        if (temporary / "gallery.html").is_file():
            support_files.append(temporary / "gallery.html")
        write_hash_manifest(temporary, support_files, temporary / "evidence.sha256")
        if final.exists():
            shutil.rmtree(final)
        temporary.rename(final)
        return final
    except Exception as error:
        if keep_failed:
            failed = output_root / f"{pack_id}.failed"
            if failed.exists():
                shutil.rmtree(failed)
            temporary.rename(failed)
            raise ValidationError(f"{error}; failed output retained at {failed}") from error
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def verify_pack(output_root: Path, pack_id: str) -> tuple[int, int, int]:
    evidence = output_root / pack_id
    if not evidence.is_dir():
        raise ValidationError(f"evidence missing: {evidence}")
    metadata = json.loads((evidence / "evidence.json").read_text(encoding="utf-8"))
    if metadata.get("pack_id") != pack_id:
        raise ValidationError(f"evidence pack ID mismatch: {pack_id}")
    if metadata.get("human_status") != "pending":
        raise ValidationError("generated evidence must retain human_status=pending")
    source_count = verify_hash_manifest(evidence, evidence / "source.sha256")
    input_count = verify_hash_manifest(evidence, evidence / "inputs.sha256")
    output_count = verify_hash_manifest(evidence, evidence / "outputs.sha256")
    support_count = verify_hash_manifest(evidence, evidence / "evidence.sha256")
    pack = load_pack(evidence / "source", pack_id)
    comparison_rows = compare_pack(evidence, pack)
    return source_count, input_count, output_count + support_count + len(comparison_rows)


def selected_pack_ids(value: str) -> tuple[str, ...]:
    return PACK_IDS if value == "all" else (value,)


def summary(root: Path, output_format: str) -> None:
    with (root / "validation/reviews.tsv").open(encoding="utf-8", newline="") as handle:
        reviews = {row["inventory_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    result = []
    for pack_id in PACK_IDS:
        pack = load_pack(root / "validation/packs" / pack_id)
        statuses = {reviews[item]["status"] for item in pack["inventory_ids"]}
        result.append(
            {
                "pack_id": pack_id,
                "inventory_ids": pack["inventory_ids"],
                "review_status": "approved" if statuses == {"approved"} else "pending",
                "recipe_status": "experimental" if any(
                    step.get("runner") == "recipe" for step in pack["steps"]
                ) else "not-applicable",
            }
        )
    if output_format == "json":
        print(json.dumps(result, indent=2))
        return
    for row in result:
        print(
            f"{row['pack_id']}\tinventory={','.join(row['inventory_ids'])}"
            f"\treview={row['review_status']}\trecipe={row['recipe_status']}"
        )


def main() -> int:
    args = options()
    root = Path(args.root).resolve()
    if args.action == "summary":
        summary(root, args.format)
        return 0
    output_root = Path(args.output).expanduser().resolve() if args.output else root / "validation/evidence"
    pack_ids = selected_pack_ids(args.pack)
    if args.action == "build":
        biohub = resolve_executable(root, args.biohub)
        snakemake = resolve_snakemake(root, args.snakemake)
        for pack_id in pack_ids:
            built = build_pack(
                root,
                pack_id,
                output_root,
                biohub,
                snakemake,
                args.force,
                args.keep_failed,
            )
            print(f"built\t{pack_id}\t{built}\thuman_status=pending")
        return 0
    for pack_id in pack_ids:
        source_count, input_count, checks = verify_pack(output_root, pack_id)
        print(
            f"verified\t{pack_id}\tsource_files={source_count}"
            f"\tinput_files={input_count}\toutputs_and_checks={checks}"
            "\thuman_status=pending"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValidationError, OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"validation review failed: {error}") from error
