"""Regression coverage for the release gate and changed DADA2 CLI contract."""
import importlib.util
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

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
