#!/usr/bin/env python3
"""Build every packaged recipe DAG from synthetic, non-biological inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=Path(__file__).resolve().parents[1])
    parser.add_argument("--snakemake", default="snakemake")
    parser.add_argument(
        "--biohub",
        help="Execute safe recipes through this BioHub binary and verify run bundles",
    )
    parser.add_argument(
        "--execute-safe",
        action="store_true",
        help="Also execute recipes needing only Python, R, DESeq2, and vegan",
    )
    parser.add_argument(
        "--recipe",
        action="append",
        dest="recipes",
        help="Restrict checks to one recipe ID; repeat for multiple recipes",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def verify_run_bundle(run_dir: Path, recipe_id: str) -> None:
    state = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    if state.get("recipe_id") != recipe_id or state.get("status") != "complete":
        raise ValueError(f"invalid completed run state for {recipe_id}: {state}")
    checksum_path = run_dir / "checksums.sha256"
    rows = checksum_path.read_text(encoding="utf-8").splitlines()
    if not rows:
        raise ValueError(f"empty checksum manifest for {recipe_id}")
    paths = set()
    for row in rows:
        expected, separator, relative = row.partition("  ")
        if not separator or not relative or len(expected) != 64:
            raise ValueError(f"invalid checksum row for {recipe_id}: {row}")
        candidate = run_dir / relative
        if not candidate.is_file() or sha256_file(candidate) != expected:
            raise ValueError(f"checksum mismatch for {recipe_id}: {relative}")
        paths.add(relative)
    if "run.json" in paths:
        raise ValueError("mutable run.json must not appear in checksums.sha256")
    required = {
        "config.resolved.yaml",
        "command.sh",
        "versions.tsv",
        "provenance.json",
        "recipe.sources.sha256",
        "inputs.manifest.tsv",
    }
    if not required.issubset(paths):
        raise ValueError(f"run bundle lacks checksummed files: {sorted(required - paths)}")


class Fixture:
    def __init__(self, root: Path):
        self.root = root

    def write(self, relative: str, content: str = "fixture\n") -> str:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return str(path.resolve())

    def directory(self, relative: str) -> str:
        path = self.root / relative
        path.mkdir(parents=True, exist_ok=True)
        return str(path.resolve())


def configs(fixture: Fixture) -> dict[str, dict[str, object]]:
    reference = fixture.write("inputs/reference.fa", ">chr1\nATGATGATGATG\n")
    query = fixture.write("inputs/query.fa", ">chr1\nATGATGATGATG\n")
    protein_dir = fixture.directory("inputs/protein_groups")
    fixture.write("inputs/protein_groups/OG0001.fa", ">tax1\nM\n>tax2\nM\n")
    cds = fixture.write("inputs/cds.fa", ">tax1\nATG\n>tax2\nATG\n")

    assembly_manifest = fixture.write(
        "inputs/assemblies.tsv", f"assembly_id\tfasta\nassembly1\t{query}\n"
    )
    synteny_manifest = fixture.write(
        "inputs/synteny.tsv",
        f"pair_id\treference_fasta\tquery_fasta\npair1\t{reference}\t{query}\n",
    )

    counts = fixture.write(
        "inputs/gene_counts.tsv",
        "Orthogroup\ttax1\ttax2\nOG0001\t1\t1\n",
    )
    tree = fixture.write("inputs/tree.nwk", "(tax1:1,tax2:1);\n")

    alignment = fixture.write("inputs/alignment.phy", " 2 3\ntax1  ATG\ntax2  ATG\n")
    marked_tree = fixture.write("inputs/marked_tree.nwk", "(tax1 #1,tax2);\n")
    selection_manifest = fixture.write(
        "inputs/selection.tsv",
        f"test_id\talignment\tmarked_tree\tforeground\ntest1\t{alignment}\t{marked_tree}\ttax1\n",
    )

    stage1 = fixture.write("inputs/stage1.ctl", "ndata = 1\nusedata = 3\n")
    stage2 = fixture.write("inputs/stage2.ctl", "ndata = 1\nusedata = 2\n")
    dating_manifest = fixture.write(
        "inputs/dating.tsv",
        f"run_id\tstage1_ctl\tstage2_ctl\nrun1\t{stage1}\t{stage2}\n",
    )

    phenotype = fixture.write("inputs/phenotype.tsv", "FID\tIID\ttrait\n0\ts1\t1\n")
    traits = fixture.write(
        "inputs/traits.tsv",
        f"trait_id\tphenotype_file\tphenotype_column\ntrait1\t{phenotype}\ttrait\n",
    )
    vcf = fixture.write(
        "inputs/genotypes.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n",
    )
    pop1 = fixture.write("inputs/pop1.txt", "s1\n")
    pop2 = fixture.write("inputs/pop2.txt", "s2\n")
    comparisons = fixture.write(
        "inputs/comparisons.tsv",
        f"comparison_id\tpopulation1_samples\tpopulation2_samples\ncmp1\t{pop1}\t{pop2}\n",
    )

    kmer_traits = fixture.write(
        "inputs/kmer_traits.tsv", f"trait_id\tphenotype_file\ntrait1\t{phenotype}\n"
    )
    kmer_prefix = str((fixture.root / "inputs/kmer_table").resolve())
    fixture.write("inputs/kmer_table.gz")
    kmer_script = fixture.write("inputs/kmers_gwas.py", "#!/usr/bin/env python2\n")

    candidates = fixture.write(
        "inputs/candidates.tsv",
        "chrom\tposition\tref\talt\ttier\tevidence_class\nchr1\t2\tT\tC\tprimary\tcandidate\n",
    )
    callable_bed = fixture.write("inputs/callable.bed", "chr1\t0\t12\n")
    denovo_manifest = fixture.write(
        "inputs/denovo.tsv",
        f"family_id\tpair_id\tcandidates_tsv\tcallable_bed\nfamily1\tpair1\t{candidates}\t{callable_bed}\n",
    )

    foreground = fixture.write("inputs/foreground.tsv", "set_id\tgene_id\nset1\tgene1\n")
    background = fixture.write("inputs/background.txt", "gene1\ngene2\n")
    associations = fixture.write(
        "inputs/associations.tsv",
        "gene_id\tterm_id\tsource\tterm_name\ngene1\tGO:1\tGO\tterm one\ngene2\tGO:1\tGO\tterm one\n",
    )

    sample_ids = [f"s{index}" for index in range(1, 9)]
    random_generator = random.Random(20260814)
    expression_lines = ["gene_id\t" + "\t".join(sample_ids)]
    for gene in range(1, 101):
        mean = 20.0 + gene * 1.5
        shape = (1.5, 4.0, 12.0)[gene % 3]
        controls = [
            max(0, round(random_generator.gammavariate(shape, mean / shape)))
            for _ in range(4)
        ]
        treatment_mean = mean * (2.0 if gene % 4 == 0 else 0.55 if gene % 5 == 0 else 1.0)
        treated = [
            max(
                0,
                round(
                    random_generator.gammavariate(
                        shape, treatment_mean / shape
                    )
                ),
            )
            for _ in range(4)
        ]
        expression_lines.append(
            f"gene{gene}\t" + "\t".join(map(str, controls + treated))
        )
    expression_counts = fixture.write(
        "inputs/expression_counts.tsv", "\n".join(expression_lines) + "\n"
    )
    samples = fixture.write(
        "inputs/samples.tsv",
        "sample_id\tcondition\n"
        + "\n".join(
            f"{sample}\t{'control' if index < 4 else 'treated'}"
            for index, sample in enumerate(sample_ids)
        )
        + "\n",
    )
    contrasts = fixture.write(
        "inputs/contrasts.tsv",
        "contrast_id\tfactor\tnumerator\tdenominator\ntreated_vs_control\tcondition\ttreated\tcontrol\n",
    )

    features = fixture.write(
        "inputs/features.tsv",
        "feature_id\ts1\ts2\ts3\ts4\ts5\ts6\n"
        "f1\t2\t1\t3\t8\t9\t7\n"
        "f2\t8\t7\t9\t2\t1\t3\n"
        "f3\t1\t2\t1\t3\t2\t4\n"
        "f4\t4\t3\t5\t1\t2\t1\n",
    )
    metadata = fixture.write(
        "inputs/metadata.tsv",
        "sample_id\tph\ns1\t5.0\ns2\t5.2\ns3\t5.5\ns4\t6.5\ns5\t6.8\ns6\t7.0\n",
    )

    return {
        "comparative-orthology-codon": {
            "protein_groups_dir": protein_dir,
            "cds_fasta": cds,
            "expected_taxa": 2,
            "biohub_executable": "biohub",
        },
        "gene-family-cafe": {
            "gene_counts": counts,
            "tree": tree,
            "input_format": "orthofinder",
            "max_family_size": 100,
            "ultrametric_tolerance": 0.001,
            "cafe_executable": "cafe5",
            "cores": 1,
            "replicates": 2,
            "model_k": 1,
            "root_distribution": 1.0,
            "error_model": None,
            "lambda_value": None,
            "alpha_value": None,
            "family_pvalue": 0.05,
            "minimum_converged_replicates": 2,
            "likelihood_tolerance": 0.01,
            "parameter_cv_tolerance": 0.1,
        },
        "selection-branch-site": {
            "tests_manifest": selection_manifest,
            "expected_taxa": 2,
            "codeml_executable": "codeml",
            "codon_frequency": 2,
            "initial_kappa": 2.0,
            "alternative_initial_omega": 1.5,
            "clean_data": 1,
            "optimizer_method": 0,
            "beb_thresholds": [0.95, 0.99],
            "bh_family": "all_manifest_tests",
        },
        "dating-mcmctree": {
            "runs_manifest": dating_manifest,
            "mcmctree_executable": "mcmctree",
            "expected_loci": 1,
            "expected_internal_nodes": 1,
            "age_column_regex": "^t_n[0-9]+$",
            "burnin_samples": 0,
            "minimum_ess": 10.0,
        },
        "assembly-t2t-evaluate": {
            "assemblies_manifest": assembly_manifest,
            "reference_fasta": reference,
            "expected_chromosomes": 1,
            "minimum_contig_length": 1,
            "telomere_motif": "TTAGGG",
            "telomere_window_bp": 12,
            "minimum_telomere_hits": 1,
            "minimap2_executable": "minimap2",
            "minimap_preset": "asm5",
            "threads": 1,
        },
        "synteny-sv": {
            "pairs_manifest": synteny_manifest,
            "require_matching_sequence_ids": True,
            "minimap2_executable": "minimap2",
            "minimap_preset": "asm5",
            "threads": 1,
            "minimum_alignment_length": 1,
            "minimum_mapq": 0,
            "syri_executable": "syri",
            "syri_cores": 1,
            "syri_use_filtered_paf": False,
            "syri_include_cigar": True,
            "syri_include_snps": False,
        },
        "population-gwas": {
            "genotype_kind": "vcf",
            "genotype": vcf,
            "traits_manifest": traits,
            "covariates_file": None,
            "covariate_columns": [],
            "require_complete_samples": True,
            "plink2_executable": "plink2",
            "minor_allele_frequency": 0.05,
            "maximum_variant_missingness": 0.1,
            "maximum_sample_missingness": 0.1,
            "hardy_weinberg_pvalue": 0.000001,
            "significance_threshold": 0.000001,
        },
        "population-selection": {
            "vcf": vcf,
            "comparisons_manifest": comparisons,
            "vcftools_executable": "vcftools",
            "window_size_bp": 1000,
            "window_step_bp": 500,
            "minor_allele_frequency": 0.05,
            "minimum_site_call_rate": 0.9,
            "candidate_fst_threshold": 0.5,
        },
        "kmer-gwas": {
            "traits_manifest": kmer_traits,
            "kmers_table_prefix": kmer_prefix,
            "kmers_gwas_script": kmer_script,
            "python2_executable": "python2",
            "kmer_length": 31,
            "threads": 1,
            "minimum_data_points": 2,
            "minor_allele_count": 1,
            "minor_allele_frequency": 0.05,
            "kmers_number": 1,
            "permutations": 10,
            "permutation_tail_percent": 5,
        },
        "family-denovo-rate": {
            "pairs_manifest": denovo_manifest,
            "included_tiers": ["primary"],
            "included_evidence_classes": ["candidate"],
            "ploidy": 2,
            "confidence_level": 0.95,
            "require_all_candidates_callable": True,
        },
        "rnaseq-deseq2": {
            "counts_matrix": expression_counts,
            "samples_tsv": samples,
            "contrasts_manifest": contrasts,
            "design": "~ condition",
            "minimum_count": 1,
            "minimum_samples": 2,
            "alpha": 0.05,
            "minimum_absolute_log2_fold_change": 1.0,
            "independent_filtering": True,
        },
        "functional-enrichment": {
            "foreground_tsv": foreground,
            "background_genes": background,
            "associations_tsv": associations,
            "sources": ["GO"],
            "minimum_term_size": 1,
            "maximum_term_size": 10,
            "minimum_overlap": 1,
            "fdr": 0.05,
            "adjustment_scope": "set_source",
            "require_foreground_in_background": True,
            "plot_top_terms": 5,
        },
        "microbiome-rda": {
            "feature_table": features,
            "metadata_tsv": metadata,
            "constraints": ["ph"],
            "condition_covariates": [],
            "transform": "hellinger",
            "minimum_prevalence": 0.0,
            "minimum_total_abundance": 0.0,
            "drop_incomplete_samples": False,
            "permutations": 9,
            "random_seed": 1,
            "scaling": 1,
        },
    }


def main() -> int:
    options = parse_args()
    root = Path(options.root).resolve()
    with tempfile.TemporaryDirectory(prefix="biohub-snakemake-smoke-") as temporary:
        fixture = Fixture(Path(temporary))
        recipe_configs = configs(fixture)
        catalog_ids = {
            line.split("\t", 1)[0]
            for line in (root / "biohub-rs/src/recipe_catalog.tsv")
            .read_text(encoding="utf-8")
            .splitlines()
            if line and not line.startswith("#")
        }
        if set(recipe_configs) != catalog_ids:
            raise SystemExit(
                f"smoke configs differ from catalog: {set(recipe_configs) ^ catalog_ids}"
            )
        selected_recipes = set(options.recipes or recipe_configs)
        unknown = selected_recipes - catalog_ids
        if unknown:
            raise SystemExit(f"unknown --recipe values: {sorted(unknown)}")
        safe_recipes = {
            "family-denovo-rate",
            "functional-enrichment",
            "microbiome-rda",
            "rnaseq-deseq2",
        }
        for recipe_id in sorted(recipe_configs):
            if recipe_id not in selected_recipes:
                continue
            dry_run_dir = fixture.root / "dry-runs" / recipe_id
            dry_run_dir.mkdir(parents=True)
            config_path = fixture.root / "configs" / f"{recipe_id}.json"
            config_path.parent.mkdir(parents=True, exist_ok=True)
            config_path.write_text(
                json.dumps(recipe_configs[recipe_id], indent=2) + "\n",
                encoding="utf-8",
            )
            command = [
                options.snakemake,
                "--snakefile",
                str(root / "recipes" / recipe_id / "Snakefile"),
                "--configfile",
                str(config_path),
                "--directory",
                str(dry_run_dir),
                "--cores",
                "1",
                "--dry-run",
                "--printshellcmds",
            ]
            environment = os.environ.copy()
            environment["XDG_CACHE_HOME"] = str(fixture.root / ".cache")
            result = subprocess.run(
                command, text=True, capture_output=True, env=environment
            )
            if result.returncode != 0:
                raise SystemExit(
                    f"Snakemake dry-run failed: {recipe_id}\n{result.stdout}\n{result.stderr}"
                )
            print(f"dry-run OK: {recipe_id}")
            if options.execute_safe and recipe_id in safe_recipes:
                run_dir = fixture.root / "runs" / recipe_id
                if options.biohub:
                    execute = [
                        options.biohub,
                        "recipe",
                        "run",
                        recipe_id,
                        "--config",
                        str(config_path),
                        "--workdir",
                        str(run_dir),
                        "--cores",
                        "1",
                    ]
                    environment["BIOHUB_RECIPE_DIR"] = str(root / "recipes")
                else:
                    run_dir.mkdir(parents=True)
                    execute = [value for value in command if value != "--dry-run"]
                    execute[execute.index("--directory") + 1] = str(run_dir)
                result = subprocess.run(
                    execute, text=True, capture_output=True, env=environment
                )
                if result.returncode != 0:
                    raise SystemExit(
                        f"Snakemake execution failed: {recipe_id}\n"
                        f"{result.stdout}\n{result.stderr}"
                    )
                if options.biohub:
                    verify_run_bundle(run_dir, recipe_id)
                print(f"execution OK: {recipe_id}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
