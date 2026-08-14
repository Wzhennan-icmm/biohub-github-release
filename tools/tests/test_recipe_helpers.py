from __future__ import annotations

import csv
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]


def run_script(relative: str, *arguments: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(ROOT / relative), *map(str, arguments)],
        text=True,
        capture_output=True,
        check=False,
    )


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class RecipeHelperTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is not installed")
    def test_psmc_plot_writes_pdf_and_rejects_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "psmc.pdf"
            command = [
                shutil.which("Rscript"),
                str(ROOT / "r/psmc_plot.R"),
                "--input", str(ROOT / "examples/psmc/merged.tsv"),
                "--stages", str(ROOT / "examples/psmc/stages.tsv"),
                "--output", str(output),
            ]
            first = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertGreater(output.stat().st_size, 100)
            second = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertNotEqual(second.returncode, 0)

    def test_private_inventory_requires_force_and_is_idempotent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "script.py").write_text("print('fixture')\n", encoding="utf-8")
            output = root / "inventory.tsv"
            first = run_script(
                "tools/audit_script_inventory.py",
                "--root", root,
                "--output", output,
            )
            self.assertEqual(first.returncode, 0, first.stderr)
            expected = output.read_bytes()
            second = run_script(
                "tools/audit_script_inventory.py",
                "--root", root,
                "--output", output,
            )
            self.assertNotEqual(second.returncode, 0)
            forced = run_script(
                "tools/audit_script_inventory.py",
                "--root", root,
                "--output", output,
                "--force",
            )
            self.assertEqual(forced.returncode, 0, forced.stderr)
            self.assertEqual(output.read_bytes(), expected)

    def test_assembly_coverage_uses_all_fasta_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            query = root / "query.fa"
            target = root / "target.fa"
            paf = root / "alignments.paf"
            output = root / "summary.tsv"
            query.write_text(">q1\n" + "A" * 100 + "\n>q2\n" + "C" * 100 + "\n", encoding="utf-8")
            target.write_text(">t1\n" + "A" * 100 + "\n>t2\n" + "G" * 100 + "\n", encoding="utf-8")
            paf.write_text("q1\t100\t0\t50\t+\tt1\t100\t0\t50\t50\t50\t60\n", encoding="utf-8")
            result = run_script(
                "recipes/assembly-t2t-evaluate/summarize_paf.py",
                "--assembly-id", "assembly1",
                "--paf", paf,
                "--query-fasta", query,
                "--target-fasta", target,
                "--output", output,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            row = read_rows(output)[0]
            self.assertEqual(row["query_sequences"], "2")
            self.assertEqual(row["target_sequences"], "2")
            self.assertAlmostEqual(float(row["query_union_coverage"]), 0.25)
            self.assertAlmostEqual(float(row["target_union_coverage"]), 0.25)

    def test_syri_summary_reads_documented_annotation_column(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pair = root / "runs" / "pair1"
            pair.mkdir(parents=True)
            (pair / "pair.syri.out").write_text(
                "Chr1\t1\t10\t-\t-\tChr1\t1\t10\tSYN1\t-\tSYN\t-\n",
                encoding="utf-8",
            )
            output = root / "summary.tsv"
            log = root / "validation.log"
            result = run_script(
                "recipes/synteny-sv/summarize.py",
                "--pair-id", "pair1",
                "--root", root / "runs",
                "--output", output,
                "--log", log,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(read_rows(output)[0]["syri_type_counts"], "SYN:1")

    def test_mcmctree_summary_accepts_whitespace_chain(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run = root / "runs" / "run1"
            run.mkdir(parents=True)
            lines = ["Gen t_n1 lnL"]
            lines.extend(f"{index} {1 + index / 100:.3f} {-10 - index / 10:.3f}" for index in range(1, 31))
            (run / "mcmc.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
            output = root / "nodes.tsv"
            summary = root / "summary.tsv"
            log = root / "validation.log"
            result = run_script(
                "recipes/dating-mcmctree/summarize.py",
                "--runs", root / "runs",
                "--run-id", "run1",
                "--age-column-regex", "^t_n[0-9]+$",
                "--burnin-samples", 0,
                "--expected-nodes", 1,
                "--minimum-ess", 1,
                "--output", output,
                "--summary", summary,
                "--log", log,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(read_rows(output)[0]["node"], "t_n1")

    def test_plink_hybrid_result_is_not_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trait = root / "gwas" / "trait1"
            trait.mkdir(parents=True)
            (trait / "trait1.PHENO1.glm.logistic.hybrid").write_text(
                "#CHROM\tPOS\tID\tTEST\tP\n1\t42\trs1\tADD\t0.001\n",
                encoding="utf-8",
            )
            status = root / "status.tsv"
            lead = root / "lead.tsv"
            log = root / "validation.log"
            result = run_script(
                "recipes/population-gwas/summarize.py",
                "--trait-id", "trait1",
                "--root", root / "gwas",
                "--threshold", 0.01,
                "--output", status,
                "--lead", lead,
                "--log", log,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(read_rows(status)[0]["tested_rows"], "1")
            self.assertEqual(read_rows(lead)[0]["variant_id"], "rs1")

    def test_population_selection_reports_pi_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            comparison = root / "comparisons" / "cmp1"
            comparison.mkdir(parents=True)
            (comparison / "fst.windowed.weir.fst").write_text(
                "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST\n"
                "1\t1\t100\t5\t0.3\t0.4\n",
                encoding="utf-8",
            )
            for population, values in [("population1", [0.1, 0.2]), ("population2", [0.3, 0.5])]:
                (comparison / f"{population}.windowed.pi").write_text(
                    "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\n"
                    + "\n".join(
                        f"1\t{index * 100 + 1}\t{(index + 1) * 100}\t5\t{value}"
                        for index, value in enumerate(values)
                    )
                    + "\n",
                    encoding="utf-8",
                )
            candidates = root / "candidates.tsv"
            summary = root / "summary.tsv"
            pi = root / "pi.tsv"
            log = root / "validation.log"
            result = run_script(
                "recipes/population-selection/summarize.py",
                "--comparison-id", "cmp1",
                "--root", root / "comparisons",
                "--threshold", 0.2,
                "--output", candidates,
                "--summary", summary,
                "--pi-summary", pi,
                "--log", log,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = read_rows(pi)
            self.assertEqual(len(rows), 2)
            self.assertAlmostEqual(float(rows[0]["mean_window_pi"]), 0.15)
            self.assertAlmostEqual(float(rows[1]["mean_window_pi"]), 0.4)


if __name__ == "__main__":
    unittest.main()
