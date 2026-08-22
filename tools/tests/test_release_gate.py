from __future__ import annotations

from pathlib import Path
import subprocess
import sys
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
        self.assertIn("pending=14", result.stdout)

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


if __name__ == "__main__":
    unittest.main()
