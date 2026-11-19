"""A valid filename/ZIP alone must never authorize taxonomy cache reuse."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import taxonomy_provenance as provenance


class TaxonomyCacheTests(unittest.TestCase):
    model_uuid = "b7986601-1851-4fa2-ace5-98bacf660fa3"
    reads_uuid = "b6fa8776-2166-423d-91f1-1130812b65e9"
    taxonomy_uuid = "a12b1865-41c7-4f82-a90e-7cd5fd891fc9"
    other_uuid = "49ccfb0a-155d-404b-80b4-818d2aeb53b2"

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.model = self.root / "silva-144-full-length-sklearn-1.7.1.qza"
        self.reads = self.root / "orientedRepSeqs.qza"
        self.taxonomy = self.root / "silva144Taxonomy.qza"
        self.artifact(self.model, self.model_uuid, "TaxonomicClassifier")
        self.artifact(self.reads, self.reads_uuid, "FeatureData[Sequence]")
        self.write_taxonomy()

    def artifact(self, path, artifact_uuid, semantic_type, action=None):
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr(f"{artifact_uuid}/metadata.yaml",
                             f"uuid: {artifact_uuid}\ntype: {semantic_type}\n")
            if action is not None:
                archive.writestr(f"{artifact_uuid}/provenance/action/action.yaml", action)

    def write_taxonomy(self, *, model=None, reads=None, confidence="0.7"):
        action = (
            "action:\n    type: method\n"
            "    plugin: !ref 'environment:plugins:feature-classifier'\n"
            "    action: classify_sklearn\n    inputs:\n"
            f"    -   reads: {reads or self.reads_uuid}\n"
            f"    -   classifier: {model or self.model_uuid}\n"
            f"    parameters:\n    -   confidence: {confidence}\n"
        )
        self.artifact(self.taxonomy, self.taxonomy_uuid, "FeatureData[Taxonomy]", action)

    def check(self, confidence=0.7):
        return provenance.check_taxonomy_provenance(
            self.taxonomy, self.model, reads_path=self.reads, confidence=confidence)

    def invoke(self):
        return subprocess.run(
            [sys.executable, ROOT / "scripts/taxonomy_provenance.py",
             "--taxonomy", self.taxonomy, "--classifier", self.model,
             "--reads", self.reads, "--confidence", ".7"],
            capture_output=True, text=True, timeout=10)

    def test_correct_model_reads_and_confidence_allow_reuse(self):
        self.assertEqual(self.check(), {
            "taxonomy_uuid": self.taxonomy_uuid, "classifier_uuid": self.model_uuid,
            "reads_uuid": self.reads_uuid, "confidence": 0.7})
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout + result.stderr, "")

    def test_valid_taxonomy_from_another_model_cannot_hide_under_expected_filename(self):
        self.write_taxonomy(model=self.other_uuid)
        with zipfile.ZipFile(self.taxonomy) as archive:
            self.assertIsNone(archive.testzip())
        with self.assertRaisesRegex(ValueError, "classifier UUID mismatch"):
            self.check()
        result = self.invoke()
        self.assertEqual(result.returncode, 1)
        self.assertIn("Cached taxonomy rejected", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_result_from_other_reads_is_rejected_even_with_correct_model(self):
        self.write_taxonomy(reads=self.other_uuid)
        with self.assertRaisesRegex(ValueError, "reads UUID mismatch"):
            self.check()

    def test_confidence_must_match_numerically_and_be_finite(self):
        self.write_taxonomy(confidence="0.70")
        self.assertEqual(self.check()["confidence"], 0.7)
        with self.assertRaisesRegex(ValueError, "confidence mismatch"):
            self.check(0.9)
        for invalid in ("nan", "inf", "-0.1", "1.1", "null", "[0.7]"):
            with self.subTest(value=invalid):
                self.write_taxonomy(confidence=invalid)
                with self.assertRaisesRegex(ValueError, "confidence"):
                    self.check()
        self.write_taxonomy()
        for invalid in (float("nan"), float("inf")):
            with self.subTest(requested=invalid), self.assertRaisesRegex(ValueError, "confidence"):
                self.check(invalid)

    def test_malformed_or_missing_provenance_rejects_cache_without_traceback(self):
        for action in (None, "action: [", "action: []", "action:\n    type: import\n"):
            with self.subTest(action=action):
                self.artifact(self.taxonomy, self.taxonomy_uuid, "FeatureData[Taxonomy]", action)
                result = self.invoke()
                self.assertEqual(result.returncode, 1)
                self.assertIn("Cached taxonomy rejected", result.stderr)
                self.assertNotIn("Traceback", result.stderr)

    def test_renamed_table_is_not_a_reusable_taxonomy(self):
        self.artifact(self.taxonomy, self.taxonomy_uuid, "FeatureTable[Frequency]")
        with self.assertRaisesRegex(ValueError, r"expected FeatureData\[Taxonomy\]"):
            self.check()


if __name__ == "__main__":
    unittest.main()
