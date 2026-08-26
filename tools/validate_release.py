#!/usr/bin/env python3
"""Validate BioHub release-readiness evidence without inventing approvals."""

from __future__ import annotations

import argparse
import csv
from datetime import date
from pathlib import Path
import re
import subprocess


REVIEW_COLUMNS = {
    "inventory_id",
    "target",
    "review_class",
    "status",
    "evidence",
    "reviewer",
    "reviewed_on",
    "note",
}
REVIEW_CLASSES = {"golden", "domain", "external", "visual", "recipe"}
REVIEW_STATUSES = {"automated", "approved", "pending"}
RELEASABLE_STATUSES = {"automated", "approved"}
RELEASED_DECISIONS = {"integrated", "recipe", "retained"}
GENERIC_AUTHOR_MARKERS = {"biohub contributors", "contributors", "anonymous"}


def options() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=Path(__file__).resolve().parents[1])
    parser.add_argument("--release", action="store_true")
    parser.add_argument("--tag", help="Expected release tag, for example v0.4.0")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_date(value: str, label: str) -> None:
    try:
        date.fromisoformat(value)
    except ValueError as error:
        raise ValueError(f"invalid ISO date for {label}: {value}") from error


def approval_record_field(text: str, label: str, evidence: Path) -> str:
    prefix = f"- {label}:"
    values = [
        line.removeprefix(prefix).strip()
        for line in text.splitlines()
        if line.startswith(prefix)
    ]
    if len(values) != 1 or not values[0]:
        raise ValueError(f"approval record {evidence} needs exactly one non-empty {label} field")
    return values[0]


def validate_approval_record(evidence: Path, review: dict[str, str]) -> None:
    text = evidence.read_text(encoding="utf-8")
    inventory_ids = set(
        re.findall(
            r"(?<!\d)\d{3}(?!\d)",
            approval_record_field(text, "Inventory IDs", evidence),
        )
    )
    if review["inventory_id"] not in inventory_ids:
        raise ValueError(
            f"approval record {evidence} does not cover inventory {review['inventory_id']}"
        )
    reviewer = approval_record_field(text, "Reviewer", evidence)
    if reviewer != review["reviewer"]:
        raise ValueError(f"approval record reviewer differs for {review['inventory_id']}")
    affiliation = approval_record_field(text, "Reviewer affiliation", evidence)
    if affiliation in {"-", "unknown"}:
        raise ValueError(f"approval record affiliation missing for {review['inventory_id']}")
    reviewed_on = approval_record_field(text, "Review date", evidence)
    validate_date(reviewed_on, f"approval record {review['inventory_id']}")
    if reviewed_on != review["reviewed_on"]:
        raise ValueError(f"approval record date differs for {review['inventory_id']}")
    if approval_record_field(text, "Decision", evidence).lower() != "approved":
        raise ValueError(f"approval record decision differs for {review['inventory_id']}")
    commit = approval_record_field(text, "BioHub head commit", evidence).strip("`")
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError(f"approval record commit is invalid for {review['inventory_id']}")
    for label in ("Input-manifest SHA256", "Output-manifest SHA256"):
        digest = approval_record_field(text, label, evidence).strip("`")
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError(f"approval record {label} is invalid for {review['inventory_id']}")
    data_license = approval_record_field(text, "Data/license confirmation", evidence)
    if data_license in {"-", "unknown"}:
        raise ValueError(f"approval record data/license confirmation missing for {review['inventory_id']}")


def package_version(root: Path) -> str:
    cargo = (root / "biohub-rs/Cargo.toml").read_text(encoding="utf-8")
    match = re.search(r'^version = "([^"]+)"$', cargo, flags=re.MULTILINE)
    if not match:
        raise ValueError("Cargo package version missing")
    return match.group(1)


def valid_orcid(value: str) -> bool:
    identifier = value.removeprefix("https://orcid.org/").replace("-", "")
    if not re.fullmatch(r"\d{15}[\dX]", identifier):
        return False
    total = 0
    for character in identifier[:15]:
        total = (total + int(character)) * 2
    result = (12 - total % 11) % 11
    expected = "X" if result == 10 else str(result)
    return identifier[-1] == expected


def validate_citation_metadata(
    citation: str,
    release: bool,
    changelog_state: str,
) -> None:
    lowered = citation.lower()
    for marker in GENERIC_AUTHOR_MARKERS:
        if re.search(
            rf"^\s*-?\s*name:\s*[\"']?{re.escape(marker)}[\"']?\s*$",
            lowered,
            re.MULTILINE,
        ):
            raise ValueError("CITATION.cff uses a generic author")
    if not re.search(r"^\s+-?\s*family-names:\s*\S+", citation, re.MULTILINE):
        raise ValueError("CITATION.cff needs structured family-names")
    if not re.search(r"^\s+-?\s*given-names:\s*\S+", citation, re.MULTILINE):
        raise ValueError("CITATION.cff needs structured given-names")
    if not re.search(r"^\s+affiliation:\s*\S+", citation, re.MULTILINE):
        raise ValueError("CITATION.cff needs author affiliation")
    if not re.search(r"^abstract:\s*(?:>|\S)", citation, re.MULTILINE):
        raise ValueError("CITATION.cff abstract missing")
    if not re.search(r"^keywords:\s*$", citation, re.MULTILINE):
        raise ValueError("CITATION.cff keywords missing")

    orcids = re.findall(r"https://orcid\.org/(\d{4}-\d{4}-\d{4}-[\dX]{4})", citation)
    if not orcids:
        raise ValueError("CITATION.cff needs an ORCID URL")
    for identifier in orcids:
        if not valid_orcid(identifier):
            raise ValueError(f"CITATION.cff has invalid ORCID checksum: {identifier}")

    released = re.search(
        r"^date-released:\s*[\"']?(\d{4}-\d{2}-\d{2})[\"']?\s*$",
        citation,
        re.MULTILINE,
    )
    if changelog_state == "Unreleased" and released:
        raise ValueError("CITATION.cff date-released must be absent while CHANGELOG is Unreleased")
    if release:
        if not released:
            raise ValueError("formal release requires CITATION.cff date-released")
        validate_date(released.group(1), "CITATION.cff date-released")
        if released.group(1) != changelog_state:
            raise ValueError("CITATION.cff date-released differs from CHANGELOG release date")


def validate_versions(root: Path, release: bool, tag: str | None) -> str:
    version = package_version(root)
    if release and not tag:
        raise ValueError("formal release validation requires --tag")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
    command_guide = (root / "docs/USER_GUIDE.zh-CN.md").read_text(encoding="utf-8")
    recipe_guide = (root / "docs/RECIPES.zh-CN.md").read_text(encoding="utf-8")
    if f"version: {version}" not in citation:
        raise ValueError("CITATION.cff version differs from Cargo version")
    if f"文档版本：{version}" not in command_guide or f"文档版本：{version}" not in recipe_guide:
        raise ValueError("Chinese guide version differs from Cargo version")
    heading = re.search(rf"^## {re.escape(version)} - (.+)$", changelog, flags=re.MULTILINE)
    if not heading:
        raise ValueError("CHANGELOG release heading missing")
    changelog_state = heading.group(1).strip()
    if release and changelog_state == "Unreleased":
        raise ValueError("CHANGELOG still marks release as Unreleased")
    if changelog_state != "Unreleased":
        validate_date(changelog_state, "CHANGELOG release heading")
    validate_citation_metadata(citation, release, changelog_state)
    if tag and tag != f"v{version}":
        raise ValueError(f"tag {tag} differs from package version v{version}")
    return version


def validate_reviews(root: Path, release: bool) -> tuple[int, list[str]]:
    matrix_path = root / "docs/SCRIPT_MIGRATION_MATRIX.tsv"
    reviews_path = root / "validation/reviews.tsv"
    matrix = read_tsv(matrix_path)
    reviews = read_tsv(reviews_path)
    if not matrix or not reviews:
        raise ValueError("migration matrix and review register must be non-empty")
    if set(reviews[0]) != REVIEW_COLUMNS:
        raise ValueError("review register columns differ from contract")

    matrix_by_id = {row["inventory_id"]: row for row in matrix}
    if len(matrix_by_id) != len(matrix):
        raise ValueError("duplicate migration inventory ID")
    review_by_id = {row["inventory_id"]: row for row in reviews}
    if len(review_by_id) != len(reviews):
        raise ValueError("duplicate review inventory ID")

    pending_matrix = {
        row["inventory_id"]
        for row in matrix
        if "pending" in row["verification"] and row["decision"] in RELEASED_DECISIONS
    }
    missing_reviews = pending_matrix - set(review_by_id)
    if missing_reviews:
        raise ValueError(f"pending migration rows lack review entries: {sorted(missing_reviews)}")

    pending: list[str] = []
    for inventory_id, review in review_by_id.items():
        migration = matrix_by_id.get(inventory_id)
        if migration is None:
            raise ValueError(f"review references unknown inventory ID: {inventory_id}")
        if review["target"] != migration["target"]:
            raise ValueError(f"review target mismatch for {inventory_id}")
        if review["review_class"] not in REVIEW_CLASSES:
            raise ValueError(f"invalid review class for {inventory_id}")
        if review["status"] not in REVIEW_STATUSES:
            raise ValueError(f"invalid review status for {inventory_id}")
        evidence = Path(review["evidence"])
        if evidence.is_absolute() or not review["evidence"]:
            raise ValueError(f"evidence path must be relative for {inventory_id}")
        evidence_path = root / evidence
        if not evidence_path.is_file():
            raise ValueError(f"evidence file missing for {inventory_id}: {evidence}")

        matrix_is_pending = "pending" in migration["verification"]
        if review["status"] == "pending":
            pending.append(f"{inventory_id}:{review['target']}")
            if not matrix_is_pending:
                raise ValueError(f"pending review has non-pending matrix status: {inventory_id}")
            continue

        if matrix_is_pending:
            raise ValueError(f"completed review still marked pending in matrix: {inventory_id}")
        if review["status"] == "automated" and review["reviewer"] != "ci":
            raise ValueError(f"automated review must use reviewer=ci: {inventory_id}")
        if review["status"] == "approved" and review["reviewer"] in {"", "-", "ci"}:
            raise ValueError(f"approved review needs named human reviewer: {inventory_id}")
        validate_date(review["reviewed_on"], inventory_id)
        if review["status"] == "approved":
            validate_approval_record(evidence_path, review)

    if release:
        matrix_pending = [
            f"{row['inventory_id']}:{row['target']}"
            for row in matrix
            if row["decision"] in RELEASED_DECISIONS and "pending" in row["verification"]
        ]
        blocked = sorted(set(pending + matrix_pending))
        if blocked:
            raise ValueError("release blocked by pending validation: " + ", ".join(blocked))
        invalid = [
            f"{row['inventory_id']}:{row['status']}"
            for row in reviews
            if row["status"] not in RELEASABLE_STATUSES
        ]
        if invalid:
            raise ValueError("release review register incomplete: " + ", ".join(invalid))

    return len(reviews), sorted(pending)


def validate_private_markers(root: Path) -> None:
    result = subprocess.run(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
        cwd=root,
        check=True,
        capture_output=True,
    )
    private_markers = ("/" + "Users/", "/mnt/" + "kobe/")
    secret_pattern = re.compile(r"(?<![A-Za-z0-9])sk-[A-Za-z0-9_-]{20,}")
    for relative in result.stdout.decode("utf-8").split("\0"):
        if not relative:
            continue
        path = root / relative
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for marker in private_markers:
            if marker in text:
                raise ValueError(f"private marker {marker!r} in {path.relative_to(root)}")
        if secret_pattern.search(text):
            raise ValueError(f"possible API credential in {path.relative_to(root)}")


def main() -> int:
    args = options()
    root = Path(args.root).resolve()
    version = validate_versions(root, args.release, args.tag)
    review_count, pending = validate_reviews(root, args.release)
    validate_private_markers(root)
    print(
        f"version={version} reviews={review_count} pending={len(pending)} "
        f"release_gate={'pass' if args.release else 'not-requested'}"
    )
    if pending and not args.release:
        print("pending_reviews=" + ",".join(pending))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        raise SystemExit(f"release validation failed: {error}") from error
