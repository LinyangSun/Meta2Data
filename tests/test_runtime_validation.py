"""Regression coverage for the release gate and changed DADA2 CLI contract."""
import importlib.util
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "runtime_validation", ROOT / "tests/integration/validate_runtime.py")
runtime = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runtime)


class RuntimeGateTests(unittest.TestCase):
    def test_equal_counts_cannot_hide_mismatched_feature_ids(self):
        with self.assertRaisesRegex(RuntimeError, "IDs differ"):
            runtime.check_identifiers(["feature-a"], ["feature-b"])

    def test_tree_and_taxonomy_must_contain_every_retained_feature(self):
        with self.assertRaisesRegex(RuntimeError, "taxonomy"):
            runtime.check_identifiers(["a", "b"], ["a", "b"], ["a"])
        with self.assertRaisesRegex(RuntimeError, "tree"):
            runtime.check_identifiers(["a", "b"], ["a", "b"], ["a", "b"], ["a"])
        # A SEPP tree also contains its reference tips; taxonomy can include
        # oriented features that were subsequently removed by tree placement.
        runtime.check_identifiers(["a"], ["a"], ["a", "removed"], ["a", "reference"])

    def test_zero_feature_results_are_not_a_successful_fixture(self):
        with self.assertRaisesRegex(RuntimeError, "No features"):
            runtime.check_identifiers([], [])

    def test_previous_validation_evidence_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            evidence = path / "validation.json"
            evidence.write_text('{"status":"passed"}\n')
            with self.assertRaisesRegex(RuntimeError, "never overwritten"):
                runtime.Validation(path, "environment", 1)
            self.assertEqual(evidence.read_text(), '{"status":"passed"}\n')


class TaxonomyIdentityTests(unittest.TestCase):
    classifier_uuid = "b7986601-1851-4fa2-ace5-98bacf660fa3"
    other_uuid = "49ccfb0a-155d-404b-80b4-818d2aeb53b2"
    taxonomy_uuid = "a12b1865-41c7-4f82-a90e-7cd5fd891fc9"

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.classifier = self.root / "silva-144-full-length-sklearn-1.7.1.qza"
        self.taxonomy = self.root / "silva144Taxonomy.qza"
        self.write_artifact(self.classifier, self.classifier_uuid, "TaxonomicClassifier")

    def write_artifact(self, path, artifact_uuid, semantic_type, action=None, extra=None):
        # Metadata-only fixtures exercise archive provenance without QIIME or
        # loading a classifier model. Full artifact validation is separate.
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr(f"{artifact_uuid}/metadata.yaml",
                             f"uuid: {artifact_uuid}\ntype: {semantic_type}\n")
            if action is not None:
                archive.writestr(f"{artifact_uuid}/provenance/action/action.yaml", action)
            for name, content in (extra or {}).items():
                archive.writestr(f"{artifact_uuid}/provenance/{name}", content)

    def action(self, classifier_uuid=None):
        return (
            "action:\n"
            "    type: method\n"
            "    plugin: !ref 'environment:plugins:feature-classifier'\n"
            "    action: classify_sklearn\n"
            "    inputs:\n"
            "    -   reads: b6fa8776-2166-423d-91f1-1130812b65e9\n"
            f"    -   classifier: {classifier_uuid or self.classifier_uuid}\n"
            "    parameters:\n"
            "    -   confidence: 0.7\n"
            "    output-name: classification\n"
            "    citations:\n"
            "    - !cite 'action|feature-classifier:2026.7.0|method:classify_sklearn|0'\n"
        )

    def test_classifier_identity_is_read_without_loading_model(self):
        self.write_artifact(self.taxonomy, self.taxonomy_uuid, "FeatureData[Taxonomy]", self.action())
        self.assertEqual(runtime.check_taxonomy_classifier(self.taxonomy, self.classifier),
                         {"taxonomy_uuid": self.taxonomy_uuid, "classifier_uuid": self.classifier_uuid})

    def test_correct_filename_and_ancestor_cannot_hide_wrong_classifier(self):
        self.write_artifact(
            self.taxonomy, self.taxonomy_uuid, "FeatureData[Taxonomy]", self.action(self.other_uuid),
            extra={f"artifacts/{self.classifier_uuid}/metadata.yaml":
                   f"uuid: {self.classifier_uuid}\ntype: TaxonomicClassifier\n"})
        with self.assertRaisesRegex(RuntimeError, "classifier UUID mismatch"):
            runtime.check_taxonomy_classifier(self.taxonomy, self.classifier)

    def test_missing_or_ambiguous_classifier_input_is_rejected(self):
        for action in (
                self.action().replace(f"    -   classifier: {self.classifier_uuid}\n", ""),
                self.action().replace("    parameters:",
                                      f"    -   classifier: {self.other_uuid}\n    parameters:")):
            with self.subTest(action=action):
                self.write_artifact(self.taxonomy, self.taxonomy_uuid, "FeatureData[Taxonomy]", action)
                with self.assertRaisesRegex(RuntimeError, "classifier UUID mismatch"):
                    runtime.check_taxonomy_classifier(self.taxonomy, self.classifier)

    def test_imported_taxonomy_or_missing_provenance_is_not_a_classification_result(self):
        self.write_artifact(self.taxonomy, self.taxonomy_uuid, "FeatureData[Taxonomy]",
                            "action:\n    type: import\n")
        with self.assertRaisesRegex(RuntimeError, "not produced directly"):
            runtime.check_taxonomy_classifier(self.taxonomy, self.classifier)
        self.write_artifact(self.taxonomy, self.taxonomy_uuid, "FeatureData[Taxonomy]")
        with self.assertRaisesRegex(RuntimeError, "Cannot read classifier provenance"):
            runtime.check_taxonomy_classifier(self.taxonomy, self.classifier)

    def test_renamed_nonclassifier_input_is_rejected(self):
        self.write_artifact(self.classifier, self.classifier_uuid, "FeatureData[Sequence]")
        self.write_artifact(self.taxonomy, self.taxonomy_uuid, "FeatureData[Taxonomy]", self.action())
        with self.assertRaisesRegex(RuntimeError, "expected TaxonomicClassifier"):
            runtime.check_taxonomy_classifier(self.taxonomy, self.classifier)


class Silva144RankTests(unittest.TestCase):
    lineage = ("d__Bacteria;k__Pseudomonadati;p__Pseudomonadota;c__Alphaproteobacteria;"
               "o__Rickettsiales;f__Mitochondria;g__Incertae_Sedis")

    def test_full_truncated_and_unassigned_lineages_are_valid(self):
        ranks = self.lineage.split(";")
        lineages = ["; ".join(ranks[:count]) for count in range(1, 8)]
        report = runtime.check_silva144_ranks([*lineages, "Unassigned"])
        self.assertEqual(report["observed_rank_counts"], list(range(8)))
        self.assertEqual(report["unassigned_features"], 1)
        # A confidence threshold may leave no fully resolved lineage.
        self.assertEqual(runtime.check_silva144_ranks(["Unassigned"])["observed_rank_counts"], [0])

    def test_missing_kingdom_reordered_extra_or_unprefixed_ranks_fail(self):
        for lineage in (
                self.lineage.replace("k__Pseudomonadati;", ""),
                "d__Bacteria;p__Pseudomonadota;k__Pseudomonadati",
                self.lineage + ";s__species",
                "Bacteria;Pseudomonadati",
                "d__Bacteria;Unassigned",
                "d__Bacteria;",
                "", None, float("nan")):
            with self.subTest(lineage=lineage), self.assertRaisesRegex(RuntimeError, "SILVA 144"):
                runtime.check_silva144_ranks([lineage])


class Dada2OutputContractTests(unittest.TestCase):
    def capture(self, action, qza_type, platform):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = root / "dataset"
            dataset.mkdir()
            capture = root / "arguments.txt"
            environment = dict(
                os.environ, dataset_path=str(dataset), SCRIPTS=str(ROOT / "scripts"),
                CAPTURE=str(capture), QZA_TYPE=qza_type, platform=platform,
                cpu="1", ION_MAXEE="1.0", DADA2_PACBIO_MIN_LENGTH="1000",
                DADA2_PACBIO_MAX_LENGTH="1600", primer_front="AGAGTTTGATCCTGGCTCAG",
                primer_adapter="", READ_COUNTS_STATE="")
            function = ("Amplicon_Pacbio_DenosingDada2" if action == "denoise-ccs"
                        else "Amplicon_Illumina_DenosingDada2 -s 0 -e 240")
            script = (
                "source " + shlex.quote(str(ROOT / "scripts/AmpliconFunction.sh")) + "\n"
                'Audit_Dada2() { return 0; }\n'
                'qiime() {\n'
                '  if [[ "$1 $2" == "tools peek" ]]; then\n'
                '    printf "Type: %s\\n" "$QZA_TYPE"\n'
                '  else\n'
                '    printf "%s\\n" "$@" > "$CAPTURE"\n'
                '  fi\n'
                '}\n' + function + "\n")
            result = subprocess.run(["bash", "-c", script], env=environment,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            arguments = capture.read_text().splitlines()
            self.assertEqual(arguments[:2], ["dada2", action])
            outputs = {}
            for name in ("table", "representative-sequences", "denoising-stats", "base-transition-stats"):
                flag = "--o-" + name
                self.assertEqual(arguments.count(flag), 1, arguments)
                outputs[name] = arguments[arguments.index(flag) + 1]
            self.assertEqual(len(set(outputs.values())), 4, "DADA2 outputs must not overwrite each other")
            self.assertTrue(outputs["base-transition-stats"].endswith("/dataset-base-transition-stats.qza"))

    def test_each_supported_dada2_entrypoint_supplies_four_distinct_outputs(self):
        for action, qza_type, platform in (
                ("denoise-paired", "SampleData[PairedEndSequencesWithQuality]", "ILLUMINA"),
                ("denoise-single", "SampleData[SequencesWithQuality]", "ILLUMINA"),
                ("denoise-pyro", "SampleData[SequencesWithQuality]", "ION_TORRENT"),
                ("denoise-ccs", "SampleData[SequencesWithQuality]", "PACBIO_SMRT")):
            with self.subTest(action=action):
                self.capture(action, qza_type, platform)


if __name__ == "__main__":
    unittest.main()
