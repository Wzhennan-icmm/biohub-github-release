from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from tools import validate_release


ROOT = Path(__file__).resolve().parents[2]


class ReleaseGateTests(unittest.TestCase):
    def run_gate(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(ROOT / "tools/validate_release.py"), *arguments],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_review_register_is_consistent(self) -> None:
        result = self.run_gate()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("reviews=48", result.stdout)
        self.assertIn("pending=11", result.stdout)

    def test_approved_record_binds_reviewer_date_and_inventory(self) -> None:
        review = next(
            row
            for row in validate_release.read_tsv(ROOT / "validation/reviews.tsv")
            if row["inventory_id"] == "014"
        )
        evidence = ROOT / review["evidence"]
        validate_release.validate_approval_record(evidence, review)
        with tempfile.TemporaryDirectory() as temporary:
            altered = Path(temporary) / "review.md"
            altered.write_text(
                evidence.read_text(encoding="utf-8").replace(
                    "- Reviewer: Zhennan Wang", "- Reviewer: Someone Else"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "reviewer differs"):
                validate_release.validate_approval_record(altered, review)

    def test_release_requires_final_changelog(self) -> None:
        result = self.run_gate("--release", "--tag", "v0.4.0")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Unreleased", result.stderr)

    def test_pending_reviews_block_after_metadata_is_finalized(self) -> None:
        with self.assertRaisesRegex(ValueError, "release blocked by pending validation"):
            validate_release.validate_reviews(ROOT, release=True)

    def test_formal_release_requires_tag(self) -> None:
        result = self.run_gate("--release")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires --tag", result.stderr)

    def test_tag_must_match_package_version(self) -> None:
        result = self.run_gate("--tag", "v9.9.9")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("differs from package version", result.stderr)

    def test_publication_orcid_checksum_is_valid(self) -> None:
        self.assertTrue(validate_release.valid_orcid("0000-0003-4883-2538"))
        self.assertTrue(
            validate_release.valid_orcid("https://orcid.org/0000-0003-4883-2538")
        )
        self.assertFalse(validate_release.valid_orcid("0000-0003-4883-2539"))

    def test_formal_citation_requires_release_date(self) -> None:
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "date-released"):
            validate_release.validate_citation_metadata(citation, True, "2026-08-25")


if __name__ == "__main__":
    unittest.main()
