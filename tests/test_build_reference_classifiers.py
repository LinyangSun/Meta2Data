"""Source identity and environment gates for the maintainer-only model rebuild."""
import contextlib
import copy
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import uuid

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import build_reference_classifiers as builder
import reference_resources as refs
from test_reference_resources import make_artifact


class BuildReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_legacy_silva_recipe_retains_the_actual_138_1_source_and_filters(self):
        plan = builder.recipe("silva-138.1")
        self.assertEqual(plan["original_classifier_uuid"],
                         "70b4b5f4-8fce-40bd-b508-afacbc12a5ed")
        self.assertTrue(all("/release_138.1/" in source["url"]
                            for source in plan["sources"].values()))
        steps = {step["action"]: step for step in plan["steps"]}
        taxonomy = steps["parse_silva_taxonomy"]["parameters"]
        self.assertFalse(taxonomy["include_species_labels"])
        self.assertEqual(taxonomy["ranks"],
                         ["domain", "phylum", "class", "order", "family", "genus"])
        self.assertEqual(steps["filter_seqs_length_by_taxon"]["parameters"]["min_lens"],
                         [900, 1200, 1400])
        self.assertEqual(steps["dereplicate"]["parameters"]["mode"], "uniq")
        # These are JSON-encoded strings in both original provenance and the
        # QIIME action signature; Python None fails QIIME's Str validation.
        fit = steps["fit_classifier_naive_bayes"]["parameters"]
        for name in ("classify__class_prior", "feat_ext__preprocessor",
                     "feat_ext__stop_words", "feat_ext__strip_accents", "feat_ext__tokenizer"):
            self.assertEqual(fit[name], "null")
        with self.assertRaises(ValueError):
            builder.recipe("silva-138.2")

    def test_original_training_input_identity_is_required(self):
        payload = b">known\nACGT\n"
        artifact_uuid = str(uuid.uuid4())
        path = make_artifact(self.root / "source.qza", payload=payload,
                             artifact_uuid=artifact_uuid)
        source = {"type": "FeatureData[Sequence]", "uuid": artifact_uuid,
                  "payload_md5": hashlib.md5(payload).hexdigest()}
        builder.validate_source(path, source)
        with self.assertRaisesRegex(ValueError, "UUID"):
            builder.validate_source(path, {**source, "uuid": str(uuid.uuid4())})
        with self.assertRaisesRegex(ValueError, "original classifier input"):
            builder.validate_source(path, {**source, "payload_md5": "0" * 32})

    def test_empty_discarded_sequences_are_allowed_but_empty_references_are_not(self):
        path = make_artifact(self.root / "empty.qza", payload=b"")
        self.assertFalse(refs.valid_artifact(path, "FeatureData[Sequence]"))
        self.assertTrue(refs.valid_artifact(path, "FeatureData[Sequence]",
                                           require_nonempty_data=False))
        outputs = [entry for step in builder.recipe("silva-138.1")["steps"]
                   for entry in step["outputs"] if entry.get("allow_empty_data")]
        self.assertEqual([entry["id"] for entry in outputs], ["discarded"])

    def test_both_compressed_and_decompressed_sources_are_verified(self):
        raw = b"known SILVA source\n"
        compressed = gzip.compress(raw)
        source = {"url": "https://example.invalid/fixed.gz", "filename": "fixed.gz",
                  "md5": hashlib.md5(compressed).hexdigest(),
                  "uncompressed_md5": hashlib.md5(raw).hexdigest()}
        calls = []
        def downloader(url, destination):
            calls.append(url)
            Path(destination).write_bytes(compressed)
        path, record = builder.acquire_source(source, self.root / "ok", downloader)
        self.assertEqual(path.read_bytes(), raw)
        self.assertEqual(record["download_sha256"], hashlib.sha256(compressed).hexdigest())
        self.assertEqual(record["uncompressed_sha256"], hashlib.sha256(raw).hexdigest())
        builder.acquire_source(source, self.root / "ok", downloader)
        self.assertEqual(len(calls), 1)
        with self.assertRaisesRegex(ValueError, "original import provenance"):
            builder.acquire_source({**source, "uncompressed_md5": "0" * 32},
                                   self.root / "wrong", downloader)
        self.assertFalse((self.root / "wrong/fixed").exists())
        self.assertFalse(list((self.root / "wrong").glob(".unpack-*")))
        with self.assertRaisesRegex(ValueError, "published MD5"):
            builder.acquire_source({**source, "md5": "0" * 32}, self.root / "bad", downloader)
        self.assertFalse((self.root / "bad/fixed.gz").exists())
        self.assertFalse(list((self.root / "bad").glob(".download-*")))

    def test_build_environment_cannot_generate_an_old_sklearn_model(self):
        versions = {"qiime2": "2026.7.0", "q2-feature-classifier": "2026.7.0",
                    "rescript": "2026.7.0", "scikit-learn": "1.7.1",
                    "numpy": "2.2.0", "scipy": "1.15.0"}
        with patch.object(builder.importlib.metadata, "version", side_effect=versions.__getitem__):
            self.assertEqual(builder.environment_versions()["scikit-learn"], "1.7.1")
            versions["scikit-learn"] = "1.4.2"
            with self.assertRaisesRegex(RuntimeError, "1.7.1"):
                builder.environment_versions()
            versions["scikit-learn"] = "1.7.1"
            versions["qiime2"] = "2024.10.1"
            with self.assertRaisesRegex(RuntimeError, "supported 2026.7"):
                builder.environment_versions()

    def test_changed_recipe_or_runtime_cannot_reuse_the_same_build_workspace(self):
        plan = builder.recipe("gg2-2024.09")
        versions = {"qiime2": "2026.7.0", "scikit-learn": "1.7.1"}
        first = builder.recipe_digest(plan, versions)
        changed = copy.deepcopy(plan)
        changed["steps"][-1]["parameters"]["classify__alpha"] = 1.0
        self.assertNotEqual(first, builder.recipe_digest(changed, versions))
        self.assertNotEqual(first, builder.recipe_digest(plan, {**versions, "scikit-learn": "1.4.2"}))
        self.assertEqual(first, builder.recipe_digest(copy.deepcopy(plan), dict(versions)))

    def test_plan_is_read_only_and_does_not_require_qiime(self):
        destination = self.root / "must-not-exist"
        output = io.StringIO()
        with patch.object(sys, "argv", ["build_reference_classifiers.py", "--classifier", "all",
                          "--output-dir", str(destination), "--plan"]), \
             patch.object(builder, "environment_versions", side_effect=AssertionError("environment read")), \
             patch.object(refs, "download", side_effect=AssertionError("download")), \
             contextlib.redirect_stdout(output):
            builder.main()
        plans = json.loads(output.getvalue())
        self.assertEqual([item["classifier_id"] for item in plans], ["gg2-2024.09", "silva-138.1"])
        self.assertFalse(destination.exists())


if __name__ == "__main__":
    unittest.main()
