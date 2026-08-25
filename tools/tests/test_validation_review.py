from __future__ import annotations

import csv
import json
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import validation_review  # noqa: E402


class ValidationReviewTests(unittest.TestCase):
    def test_packs_cover_every_pending_review_once(self) -> None:
        with (ROOT / "validation/reviews.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            pending = {
                row["inventory_id"]
                for row in csv.DictReader(handle, delimiter="\t")
                if row["status"] == "pending"
            }
        covered: list[str] = []
        for pack_id in validation_review.PACK_IDS:
            pack = validation_review.load_pack(ROOT / "validation/packs" / pack_id)
            covered.extend(pack["inventory_ids"])
            self.assertEqual(pack["fixture_license"], "CC0-1.0 synthetic fixture")
            self.assertTrue(pack["acceptance"])
            self.assertTrue(pack["manual_checks"])
        self.assertEqual(len(covered), len(set(covered)))
        self.assertEqual(set(covered), pending)

    def test_reviewed_recipe_packs_remain_experimental(self) -> None:
        statuses = {}
        for row in (ROOT / "biohub-rs/src/recipe_catalog.tsv").read_text(
            encoding="utf-8"
        ).splitlines():
            if not row or row.startswith("#"):
                continue
            fields = row.split("\t")
            statuses[fields[0]] = fields[3]
        self.assertEqual(statuses["comparative-orthology-codon"], "experimental")
        self.assertEqual(statuses["microbiome-rda"], "experimental")
        self.assertEqual(statuses["functional-enrichment"], "experimental")
        self.assertEqual(statuses["rnaseq-deseq2"], "experimental")

    def test_tool_exposes_no_approval_command(self) -> None:
        result = subprocess.run(
            [sys.executable, str(ROOT / "tools/validation_review.py"), "approve"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("invalid choice", result.stderr)

    def test_summary_is_machine_readable_and_stays_pending(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools/validation_review.py"),
                "summary",
                "--format",
                "json",
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = json.loads(result.stdout)
        self.assertEqual([row["pack_id"] for row in rows], list(validation_review.PACK_IDS))
        self.assertEqual({row["review_status"] for row in rows}, {"pending"})

    def test_numeric_tolerance_is_max_not_sum(self) -> None:
        self.assertTrue(validation_review.numeric_close(100.00009, 100.0, 1e-8, 1e-6))
        self.assertFalse(validation_review.numeric_close(100.00011, 100.0, 1e-8, 1e-6))
        self.assertTrue(validation_review.numeric_close(0.0, 0.0, 1e-8, 1e-6))

    def test_explicit_commit_describes_clean_snapshot_without_git(self) -> None:
        commit, dirty, state = validation_review.source_state(ROOT, "A" * 40)
        self.assertEqual(commit, "a" * 40)
        self.assertFalse(dirty)
        self.assertEqual(state, "explicit-clean-snapshot")
        with self.assertRaisesRegex(validation_review.ValidationError, "--commit"):
            validation_review.source_state(ROOT, "not-a-commit")

    def test_evidence_permissions_allow_external_artifact_reader(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "evidence"
            nested = root / "steps"
            nested.mkdir(parents=True, mode=0o700)
            output = nested / "result.tsv"
            output.write_text("value\n", encoding="utf-8")
            root.chmod(0o700)
            nested.chmod(0o700)
            output.chmod(0o600)
            validation_review.make_evidence_readable(root)
            self.assertTrue(root.stat().st_mode & stat.S_IXOTH)
            self.assertTrue(nested.stat().st_mode & stat.S_IXOTH)
            self.assertTrue(output.stat().st_mode & stat.S_IROTH)

    def test_sign_invariant_score_comparison(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            output = root / "output"
            source.mkdir()
            output.mkdir()
            (source / "expected.tsv").write_text(
                "id\tRDA1\tRDA2\na\t1\t4\nb\t2\t1\nc\t4\t2\n",
                encoding="utf-8",
            )
            (output / "actual.tsv").write_text(
                "id\tRDA1\tRDA2\na\t-1\t4\nb\t-2\t1\nc\t-4\t2\n",
                encoding="utf-8",
            )
            detail = validation_review.compare_check(
                source,
                output,
                {
                    "type": "tsv_sign_invariant",
                    "path": "actual.tsv",
                    "expected": "expected.tsv",
                    "key_columns": ["id"],
                    "numeric_columns": ["RDA1", "RDA2"],
                    "minimum_absolute_correlation": 0.999999,
                },
            )
            self.assertIn("abs(correlation)", detail)

    def test_private_researcher_id_is_not_fixture_data(self) -> None:
        for path in (ROOT / "validation/packs").rglob("*"):
            if path.is_file():
                self.assertNotIn("0000-0003-4883-2538", path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
