#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix.lower() in {".gz", ".bgz"} else path.open(encoding="utf-8")


def unique_samples(values: list[str], source: Path) -> list[str]:
    if not values or any(not value for value in values) or len(values) != len(set(values)):
        raise ValueError(f"empty or duplicate IID values: {source}")
    return values


def vcf_samples(path: Path) -> list[str]:
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                return unique_samples(fields[9:], path)
    raise ValueError("VCF has no #CHROM header")


def psam_samples(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as handle:
        rows = [line.split() for line in handle if line.strip() and not line.startswith("##")]
    if len(rows) < 2:
        raise ValueError(f"PSAM has no sample rows: {path}")
    header = [value.lstrip("#") for value in rows[0]]
    if "IID" not in header:
        raise ValueError(f"PSAM needs IID column: {path}")
    index = header.index("IID")
    if any(len(row) <= index for row in rows[1:]):
        raise ValueError(f"short PSAM row: {path}")
    return unique_samples([row[index] for row in rows[1:]], path)


def fam_samples(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as handle:
        rows = [line.split() for line in handle if line.strip()]
    if any(len(row) < 2 for row in rows):
        raise ValueError(f"FAM row needs FID and IID: {path}")
    return unique_samples([row[1] for row in rows], path)


def table_samples(path: Path, required_column: str | None = None):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"empty table: {path}")
        normalized = {name.lstrip("#"): name for name in reader.fieldnames}
        if "IID" not in normalized:
            raise ValueError(f"table needs IID column: {path}")
        if required_column and required_column not in reader.fieldnames:
            raise ValueError(f"missing phenotype column {required_column}: {path}")
        rows = list(reader)
    samples = [row[normalized["IID"]] for row in rows]
    return unique_samples(samples, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genotype-kind", required=True, choices=["vcf", "pfile", "bfile"])
    parser.add_argument("--genotype", required=True, type=Path)
    parser.add_argument("--traits", required=True, type=Path)
    parser.add_argument("--covariates", type=Path)
    parser.add_argument("--require-complete", required=True, choices=["true", "false"])
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    with args.traits.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"trait_id", "phenotype_file", "phenotype_column"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"traits manifest needs {sorted(required)}")
        traits = list(reader)
    if not traits or len({row["trait_id"] for row in traits}) != len(traits):
        raise ValueError("trait_id values must be non-empty and unique")

    genotype_samples = None
    if args.genotype_kind == "vcf":
        if not args.genotype.is_file():
            raise FileNotFoundError(args.genotype)
        genotype_samples = vcf_samples(args.genotype)
    else:
        required_suffixes = [".pgen", ".pvar", ".psam"] if args.genotype_kind == "pfile" else [".bed", ".bim", ".fam"]
        for suffix in required_suffixes:
            path = Path(str(args.genotype) + suffix)
            if not path.is_file():
                raise FileNotFoundError(path)
        genotype_samples = (
            psam_samples(Path(str(args.genotype) + ".psam"))
            if args.genotype_kind == "pfile"
            else fam_samples(Path(str(args.genotype) + ".fam"))
        )

    normalized = []
    all_trait_samples = set()
    for row in traits:
        phenotype = Path(row["phenotype_file"]).expanduser().resolve()
        samples = table_samples(phenotype, row["phenotype_column"])
        all_trait_samples.update(samples)
        normalized.append({**row, "phenotype_file": str(phenotype), "samples": len(samples)})
        missing = sorted(set(samples) - set(genotype_samples))
        absent_pheno = sorted(set(genotype_samples) - set(samples))
        if missing or (args.require_complete == "true" and absent_pheno):
            raise ValueError(
                f"sample mismatch for {row['trait_id']}: phenotype-only={missing[:10]} genotype-only={absent_pheno[:10]}"
            )
    if args.covariates:
        covariate_samples = set(table_samples(args.covariates))
        missing = sorted(all_trait_samples - covariate_samples)
        if missing:
            raise ValueError(f"covariates missing phenotype samples: {missing[:10]}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["trait_id", "phenotype_file", "phenotype_column", "samples"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(normalized)
    args.summary.write_text(
        "metric\tvalue\n"
        f"traits\t{len(normalized)}\n"
        f"unique_phenotype_samples\t{len(all_trait_samples)}\n"
        f"genotype_samples\t{len(genotype_samples)}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
