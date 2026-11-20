"""Exercise the public TAXA launcher with real input/cache logic and a stub pipeline."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from test_reference_resources import make_artifact

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))
from reference_resources import RESOURCES
import pip_state
import taxa_inputs


class TaxaCliTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.package = self.root / "package"
        (self.package / "bin").mkdir(parents=True)
        (self.package / "scripts").mkdir()
        self.launcher = self.package / "bin/Meta2Data-AmpliconTAXA"
        shutil.copyfile(REPOSITORY / "bin/Meta2Data-AmpliconTAXA", self.launcher)
        for name in ("parameters.py", "resource_profile.sh", "taxa_inputs.py"):
            shutil.copyfile(REPOSITORY / "scripts" / name, self.package / "scripts" / name)
        shutil.copyfile(REPOSITORY / "scripts/reference_resources.py",
                        self.package / "scripts/reference_resources_impl.py")
        # The launcher exercises real download/cache validation with small test
        # artifacts. Production checksums and UUIDs identify the full official
        # models, so only this copied test package uses unpinned fixtures.
        (self.package / "scripts/reference_resources.py").write_text(
            "from dataclasses import replace\n"
            "from reference_resources_impl import RESOURCES, main\n"
            "RESOURCES.update({key: replace(value, expected_uuid=None, md5=None, sha256=None)\n"
            "                  for key, value in RESOURCES.items()})\n"
            "if __name__ == '__main__':\n"
            "    main()\n")
        (self.package / "scripts/taxonomy.sh").write_text(
            "#!/bin/bash\npython3 - <<'PYTHON'\n"
            "import json, os\n"
            "from pathlib import Path\n"
            "names = ('ORIENT_REF', 'NB_CLASSIFIER', 'SEPP_REF', 'DB_LABEL', 'FINAL_DIR', 'INPUT_DIR')\n"
            "Path(os.environ['CAPTURE']).write_text(json.dumps({key: os.environ[key] for key in names}))\n"
            "PYTHON\n")
        self.work = self.root / "project work"
        self.work.mkdir()
        self.output = self.work / "results/taxa"
        self.inputs = self.work / "pip"
        make_artifact(self.inputs / "online/PRJNA1-vsearch-final-table.qza", "FeatureTable[Frequency]")
        make_artifact(self.inputs / "online/PRJNA1-vsearch-final-rep-seqs.qza")
        self.capture = self.root / "capture.json"
        no_network = self.root / "no-network"
        no_network.mkdir()
        curl = no_network / "curl"
        curl.write_text("#!/bin/sh\necho 'UNEXPECTED NETWORK ACCESS' >&2\nexit 99\n")
        curl.chmod(0o755)
        self.env = dict(os.environ, PATH=f"{no_network}:{os.environ['PATH']}", CAPTURE=str(self.capture))
        for variable in ("M2D_PROFILE_DIR", "META2DATA_SEPP_REF", "CONDA_PREFIX", "PREFIX"):
            self.env.pop(variable, None)

    def invoke(self, *args):
        return subprocess.run(["bash", str(self.launcher), *map(str, args)],
                              cwd=self.work, env=self.env, capture_output=True, text=True, timeout=20)

    def cache(self, *kinds):
        for kind in kinds:
            resource = RESOURCES[kind]
            target = self.work / "results/db" / resource.filename
            if not target.exists():
                make_artifact(target, resource.semantic_type, sklearn_version=resource.sklearn_version)

    def test_default_classifier_and_default_tree_use_shared_cache(self):
        self.cache("gg2-sequences", "greengenes-classifier", "sepp")
        result = self.invoke("-i", self.inputs)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        record = json.loads(self.capture.read_text())
        self.assertEqual(record["DB_LABEL"], "gg2")
        state = json.loads((Path(record["FINAL_DIR"]) / "taxa-run-state.json").read_text())
        self.assertEqual(state["classifier_id"], "gg2-2024.09")
        for variable in ("ORIENT_REF", "NB_CLASSIFIER", "SEPP_REF"):
            self.assertEqual(Path(record[variable]).parent, self.work / "results/db")
        self.assertEqual(Path(record["FINAL_DIR"]), self.output / "final-vsearch-multiV")
        self.assertFalse((self.output / "db").exists())

    def test_versioned_classifiers_and_legacy_aliases_select_the_expected_resource(self):
        selections = (
            ("gg2-2024.09", "gg2-2024.09", "gg2", "greengenes-classifier"),
            ("greengenes", "gg2-2024.09", "gg2", "greengenes-classifier"),
            ("silva-138.1", "silva-138.1", "silva", "silva-classifier"),
            ("silva", "silva-138.1", "silva", "silva-classifier"),
            ("silva-144", "silva-144", "silva144", "silva144-classifier"),
        )
        fingerprints = {}
        for selection, canonical, label, kind in selections:
            with self.subTest(selection=selection):
                self.cache("gg2-sequences", kind)
                result = self.invoke("--input", self.inputs, "--classifier", selection, "--notree")
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                record = json.loads(self.capture.read_text())
                self.assertEqual(record["DB_LABEL"], label)
                self.assertEqual(Path(record["NB_CLASSIFIER"]).name, RESOURCES[kind].filename)
                state = json.loads((Path(record["FINAL_DIR"]) / "taxa-run-state.json").read_text())
                self.assertEqual(state["classifier_id"], canonical)
                self.assertEqual(state["taxonomy_dependencies"][label]["classifier_id"], canonical)
                fingerprint = state["taxonomy_dependencies"][label]["fingerprint"]
                if canonical in fingerprints:
                    self.assertEqual(fingerprint, fingerprints[canonical])
                    self.assertNotIn("Recompute affected TAXA steps:", result.stdout)
                fingerprints[canonical] = fingerprint
                self.assertEqual(record["SEPP_REF"], "")

    def test_no_paths_uses_project_input_and_output_defaults(self):
        default_input = self.work / "results/pip"
        default_input.parent.mkdir()
        shutil.move(self.inputs, default_input)
        self.cache("gg2-sequences", "greengenes-classifier", "sepp")
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        record = json.loads(self.capture.read_text())
        self.assertEqual(Path(record["INPUT_DIR"]), default_input)
        self.assertEqual(Path(record["FINAL_DIR"]), self.work / "results/taxa/final-vsearch-multiV")

    def test_explicit_input_does_not_change_fixed_output(self):
        self.cache("gg2-sequences", "greengenes-classifier")
        result = self.invoke("--input", self.inputs, "--notree")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        record = json.loads(self.capture.read_text())
        self.assertEqual(Path(record["INPUT_DIR"]), self.inputs)
        self.assertEqual(Path(record["FINAL_DIR"]), self.work / "results/taxa/final-vsearch-notree")
        self.assertFalse((self.inputs / "final-vsearch-notree").exists())

    def test_removed_output_flags_fail_before_any_output_or_reference_writes(self):
        elsewhere = self.root / "elsewhere/output"
        for args in (("-o", str(elsewhere)), ("--output", str(elsewhere)),
                     ("-o",), ("--output",), (f"--output={elsewhere}",)):
            with self.subTest(args=args):
                result = self.invoke("--input", self.inputs, *args)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(f"{args[0].split('=')[0]} was removed", result.stderr)
                self.assertIn("./results/taxa", result.stderr)
                self.assertIn("./results/db", result.stderr)
                self.assertFalse((self.work / "results").exists())
                self.assertFalse(elsewhere.exists())
                self.assertFalse(self.capture.exists())
                self.assertNotIn("UNEXPECTED NETWORK ACCESS", result.stderr)

    def test_physical_launch_directory_anchors_output_and_references(self):
        self.cache("gg2-sequences", "greengenes-classifier")
        alias = self.root / "project-link"
        alias.symlink_to(self.work, target_is_directory=True)
        env = dict(self.env, PWD=str(alias))
        result = subprocess.run(["bash", str(self.launcher), "--input", str(self.inputs), "--notree"],
                                cwd=alias, env=env, capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        record = json.loads(self.capture.read_text())
        self.assertEqual(Path(record["FINAL_DIR"]), self.work / "results/taxa/final-vsearch-notree")
        self.assertEqual(Path(record["ORIENT_REF"]).parent, self.work / "results/db")

    def test_missing_default_input_fails_before_output_or_download(self):
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Default input directory", result.stderr)
        self.assertIn(str(self.work / "results/pip"), result.stderr)
        self.assertIn("--input DIR", result.stderr)
        self.assertFalse((self.work / "results").exists())
        self.assertNotIn("UNEXPECTED NETWORK ACCESS", result.stderr)

    def test_notree_and_single_region_do_not_prepare_sepp(self):
        self.cache("gg2-sequences", "silva-classifier")
        for flag, suffix in (("--notree", "notree"), ("--singleV", "singleV")):
            with self.subTest(mode=flag):
                result = self.invoke("-i", self.inputs, "--classifier", "silva", flag)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                record = json.loads(self.capture.read_text())
                self.assertEqual(record["DB_LABEL"], "silva")
                self.assertEqual(record["SEPP_REF"], "")
                self.assertEqual(Path(record["FINAL_DIR"]), self.work / "results/taxa" / f"final-vsearch-{suffix}")
        self.assertFalse((self.work / "results/db" / RESOURCES["sepp"].filename).exists())
        self.assertFalse((self.work / "results/db" / RESOURCES["greengenes-classifier"].filename).exists())

    def test_bundled_resources_prepare_automatically_without_download_flags(self):
        for kind in ("gg2-sequences", "greengenes-classifier", "sepp"):
            resource = RESOURCES[kind]
            make_artifact(self.package / "references" / resource.filename, resource.semantic_type,
                          sklearn_version=resource.sklearn_version)
        result = self.invoke("-i", self.inputs)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((self.work / "results/db" / RESOURCES["sepp"].filename).is_file())

    def test_removed_flags_have_migration_message(self):
        for flag in ("--db", "--db-type", "--dl"):
            with self.subTest(flag=flag):
                result = self.invoke(flag)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(f"{flag} was removed", result.stderr)
                self.assertIn("--classifier", result.stderr)
        self.assertFalse((self.work / "results/db").exists())

    def test_invalid_or_missing_classifier_fails_before_preparation(self):
        for args in (("--classifier",), ("--classifier", "other"), ("--classifier", "--notree"),
                     ("--classifier", "/tmp/custom-classifier.qza")):
            with self.subTest(args=args):
                result = self.invoke(*args)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("--classifier", result.stdout + result.stderr)
        self.assertFalse((self.work / "results/db").exists())

    def test_invalid_input_fails_before_preparing_references(self):
        empty = self.work / "empty"
        empty.mkdir()
        result = self.invoke("-i", empty)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No valid complete result pairs", result.stderr)
        self.assertFalse((self.work / "results/db").exists())


class TaxaStateRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.dataset = self.root / "pip/dataset"
        make_artifact(self.dataset / "dataset-vsearch-final-table.qza", "FeatureTable[Frequency]")
        make_artifact(self.dataset / "dataset-vsearch-final-rep-seqs.qza")
        self.collection = self.root / "collection.json"
        taxa_inputs.write_json(self.collection, taxa_inputs.collect(self.root / "pip", "vsearch"))
        self.orient = self.root / "orient.qza"
        self.classifier = self.root / "classifier.qza"
        self.orient.write_bytes(b"orientation-reference")
        self.classifier.write_bytes(b"classification-reference")
        self.final = self.root / "taxa"

    def prepare(self):
        return taxa_inputs.prepare(self.collection, self.final, self.orient,
                                   self.classifier, .7, notree=True)

    def managed_outputs(self):
        paths = [self.final / relative for paths in taxa_inputs.CACHE_OUTPUTS.values()
                 for relative in paths]
        for path in paths:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"cached-result")
        return paths

    def test_corrupt_state_recomputes_owned_outputs_and_preserves_other_files(self):
        self.prepare()
        state_path = self.final / "taxa-run-state.json"
        valid = json.loads(state_path.read_text())
        malformed_classifier = dict(valid, taxonomy_dependencies={"gg2": []})
        invalid_states = ('{"dependencies":', '[]', '{"dependencies": []}',
                          json.dumps(malformed_classifier))
        user_file = self.final / "notes.txt"
        user_file.write_text("keep this user file")
        original = (self.dataset / "dataset-vsearch-final-table.qza").read_bytes()
        for invalid in invalid_states:
            with self.subTest(state=invalid):
                outputs = self.managed_outputs()
                state_path.write_text(invalid)
                self.assertCountEqual(self.prepare(), taxa_inputs.CACHE_OUTPUTS)
                self.assertTrue(all(not path.exists() for path in outputs))
                self.assertEqual(user_file.read_text(), "keep this user file")
                self.assertEqual((self.dataset / "dataset-vsearch-final-table.qza").read_bytes(), original)
                self.assertEqual(json.loads(state_path.read_text())["tree_mode"], "none")
                self.assertEqual(json.loads((self.final / "collection.json").read_text()),
                                 json.loads(self.collection.read_text()))

    def test_valid_unchanged_state_preserves_cached_outputs(self):
        self.prepare()
        outputs = self.managed_outputs()
        self.assertEqual(self.prepare(), [])
        self.assertTrue(all(path.read_bytes() == b"cached-result" for path in outputs))

    def prepare_versions(self):
        """Create real state transitions and distinct outputs for each built-in database."""
        self.prepare()
        shared = [self.final / relative
                  for stage, paths in taxa_inputs.CACHE_OUTPUTS.items() if stage != "taxonomy"
                  for relative in paths]
        for path in shared:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"shared-stage-result")
        classifiers = {"gg2": self.classifier}
        outputs = {}
        for label, canonical in taxa_inputs.CLASSIFIER_IDS.items():
            if label != "gg2":
                classifiers[label] = self.root / f"{canonical}.qza"
                classifiers[label].write_bytes(canonical.encode())
            invalidated = taxa_inputs.prepare(
                self.collection, self.final, self.orient, classifiers[label], .7,
                db_label=label, notree=True, classifier_id=canonical)
            self.assertEqual(invalidated, [] if label == "gg2" else ["taxonomy"])
            outputs[label] = self.final / f"{label}Taxonomy.qza"
            outputs[label].write_bytes(canonical.encode())
            self.assertTrue(all(path.read_bytes() == b"shared-stage-result" for path in shared))
        return classifiers, outputs, shared

    def test_database_versions_cache_independently_and_keep_shared_stages(self):
        classifiers, outputs, shared = self.prepare_versions()
        for label, canonical in taxa_inputs.CLASSIFIER_IDS.items():
            with self.subTest(database=canonical):
                invalidated = taxa_inputs.prepare(
                    self.collection, self.final, self.orient, classifiers[label], .7,
                    db_label=label, notree=True, classifier_id=canonical)
                self.assertEqual(invalidated, [])
                state = json.loads((self.final / "taxa-run-state.json").read_text())
                self.assertEqual(state["classifier_id"], canonical)
                self.assertEqual(state["db_label"], label)
                self.assertEqual(set(state["taxonomy_dependencies"]), set(taxa_inputs.CLASSIFIER_IDS))
        for label, path in outputs.items():
            self.assertEqual(path.read_bytes(), taxa_inputs.CLASSIFIER_IDS[label].encode())
        self.assertTrue(all(path.read_bytes() == b"shared-stage-result" for path in shared))

    def test_classifier_replacement_invalidates_only_its_taxonomy(self):
        classifiers, outputs, shared = self.prepare_versions()
        classifiers["silva144"].write_bytes(b"updated-runtime-compatible-classifier")
        invalidated = taxa_inputs.prepare(
            self.collection, self.final, self.orient, classifiers["silva144"], .7,
            db_label="silva144", notree=True, classifier_id="silva-144")
        self.assertEqual(invalidated, ["taxonomy"])
        self.assertFalse(outputs["silva144"].exists())
        self.assertEqual(outputs["gg2"].read_bytes(), b"gg2-2024.09")
        self.assertEqual(outputs["silva"].read_bytes(), b"silva-138.1")
        self.assertTrue(all(path.read_bytes() == b"shared-stage-result" for path in shared))

    def test_orientation_change_invalidates_every_database_taxonomy(self):
        classifiers, outputs, _ = self.prepare_versions()
        self.orient.write_bytes(b"different-orientation-reference")
        invalidated = taxa_inputs.prepare(
            self.collection, self.final, self.orient, classifiers["silva144"], .7,
            db_label="silva144", notree=True, classifier_id="silva-144")
        self.assertCountEqual(invalidated, ["orient", "taxonomy", "tree"])
        self.assertTrue(all(not path.exists() for path in outputs.values()))
        state = json.loads((self.final / "taxa-run-state.json").read_text())
        self.assertEqual(set(state["taxonomy_dependencies"]), {"silva144"})
        for relative in taxa_inputs.CACHE_OUTPUTS["merge"]:
            self.assertEqual((self.final / relative).read_bytes(), b"shared-stage-result")
        for stage in ("orient", "tree"):
            for relative in taxa_inputs.CACHE_OUTPUTS[stage]:
                self.assertFalse((self.final / relative).exists())

    def test_legacy_state_without_classifier_id_recomputes_only_taxonomy(self):
        self.prepare()
        outputs = self.managed_outputs()
        state_path = self.final / "taxa-run-state.json"
        state = json.loads(state_path.read_text())
        state.pop("classifier_id")
        state.pop("db_label")
        entry = state["taxonomy_dependencies"]["gg2"]
        entry.pop("classifier_id")
        legacy_fingerprint = taxa_inputs.fingerprint({
            "orient": state["dependencies"]["orient"],
            "classifier": entry["classifier"], "confidence": entry["confidence"],
        })
        state["dependencies"]["taxonomy"] = legacy_fingerprint
        entry["fingerprint"] = legacy_fingerprint
        taxa_inputs.write_json(state_path, state)
        self.assertEqual(self.prepare(), ["taxonomy"])
        for path in outputs:
            if path.name == "gg2Taxonomy.qza":
                self.assertFalse(path.exists())
            else:
                self.assertEqual(path.read_bytes(), b"cached-result")
        refreshed = json.loads(state_path.read_text())
        self.assertEqual(refreshed["classifier_id"], "gg2-2024.09")

    def test_mismatched_classifier_identity_and_output_label_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "does not match taxonomy label"):
            taxa_inputs.prepare(self.collection, self.final, self.orient, self.classifier, .7,
                                db_label="silva", notree=True, classifier_id="silva-144")
        self.assertFalse(self.final.exists())

    def test_failed_json_publication_keeps_previous_complete_state(self):
        self.prepare()
        state_path = self.final / "taxa-run-state.json"
        before = state_path.read_bytes()
        with patch("taxa_inputs.os.replace", side_effect=OSError("simulated interruption")):
            with self.assertRaisesRegex(OSError, "simulated interruption"):
                taxa_inputs.write_json(state_path, {"replacement": "new-state"})
        self.assertEqual(state_path.read_bytes(), before)
        self.assertEqual(list(self.final.glob(".taxa-run-state.json.*")), [])

    def test_source_change_rejects_other_method_without_deleting_artifacts(self):
        listing = self.dataset / "dataset_sra.txt"
        listing.write_text("SRR_OLD\n")
        with patch.dict(os.environ, {}, clear=True), \
                patch("pip_state.code_fingerprint", return_value="stable-code"):
            pip_state.prepare(self.dataset, "vsearch")
            make_artifact(self.dataset / "dataset-vsearch-final-table.qza", "FeatureTable[Frequency]")
            make_artifact(self.dataset / "dataset-vsearch-final-rep-seqs.qza")
            self.assertEqual(len(taxa_inputs.collect(self.root / "pip", "vsearch")["datasets"]), 1)
            listing.write_text("SRR_NEW\n")
            pip_state.prepare(self.dataset, "dada2")
        original = {path: path.read_bytes() for path in self.dataset.glob("*-vsearch-final-*.qza")}
        audit = taxa_inputs.collect(self.root / "pip", "vsearch")
        self.assertEqual(audit["datasets"], [])
        self.assertEqual(len(audit["errors"]), 1)
        self.assertIn("stale vsearch results", audit["errors"][0])
        self.assertIn("--vsearch", audit["errors"][0])
        self.assertTrue(all(path.read_bytes() == value for path, value in original.items()))

    def test_external_artifacts_without_pipeline_state_remain_accepted(self):
        audit = taxa_inputs.collect(self.root / "pip", "vsearch")
        self.assertEqual(audit["errors"], [])
        self.assertEqual(len(audit["datasets"]), 1)
        # A per-method record copied along with external results is not enough
        # to establish that a newer input source replaced them.
        (self.dataset / "dataset-vsearch-run.json").write_text(
            json.dumps({"input_fingerprint": "a" * 64}))
        audit = taxa_inputs.collect(self.root / "pip", "vsearch")
        self.assertEqual(audit["errors"], [])
        self.assertEqual(len(audit["datasets"]), 1)


if __name__ == "__main__":
    unittest.main()
