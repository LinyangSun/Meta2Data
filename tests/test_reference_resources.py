"""Offline tests for validated, shared reference preparation."""
import multiprocessing
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import uuid
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import reference_resources as refs


def make_artifact(path, kind="FeatureData[Sequence]", payload=b">reference\nACGT\n"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    root = str(uuid.uuid4())
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(f"{root}/metadata.yaml", f"uuid: {root}\ntype: {kind}\nformat: TestFormat\n")
        archive.writestr(f"{root}/VERSION", "QIIME 2\narchive: 5\nframework: 2024.10\n")
        archive.writestr(f"{root}/data/reference", payload)
    return path


def concurrent_prepare(work, fixture, calls, start, replies):
    def downloader(url, destination):
        with Path(calls).open("a") as handle:
            handle.write("download\n")
        shutil.copyfile(fixture, destination)
    try:
        start.wait(5)
        path = refs.prepare_resource("gg2-sequences", work_dir=work,
                                     builtin_paths=[], downloader=downloader)
        replies.put((str(path), None))
    except Exception as error:
        replies.put((None, str(error)))


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.resource = refs.RESOURCES["gg2-sequences"]
        self.target = self.root / "results/db" / self.resource.filename

    def test_valid_cache_precedes_bundle_and_download(self):
        make_artifact(self.target)
        cached = self.target.read_bytes()
        bundle = make_artifact(self.root / "bundle.qza")
        def forbidden(*args):
            self.fail("valid cache must not download")
        result = refs.prepare_resource("gg2-sequences", self.root,
                                       builtin_paths=[bundle], downloader=forbidden)
        self.assertEqual(result, self.target)
        self.assertEqual(result.read_bytes(), cached)

    def test_valid_bundle_is_copied_into_working_directory(self):
        bundle = make_artifact(self.root / "bundle.qza")
        def forbidden(*args):
            self.fail("valid bundle must not download")
        result = refs.prepare_resource("gg2-sequences", self.root,
                                       builtin_paths=[bundle], downloader=forbidden)
        self.assertEqual(result.read_bytes(), bundle.read_bytes())
        self.assertFalse(result.is_symlink())
        bundle.unlink()
        self.assertTrue(refs.valid_artifact(result, self.resource.semantic_type))

    def test_wrong_type_cache_is_replaced_and_invalid_bundle_skipped(self):
        make_artifact(self.target, "TaxonomicClassifier")
        bad_bundle = make_artifact(self.root / "wrong.qza", "FeatureTable[Frequency]")
        fetched = []
        def downloader(url, destination):
            fetched.append(url)
            make_artifact(destination)
        result = refs.prepare_resource("gg2-sequences", self.root,
                                       builtin_paths=[bad_bundle], downloader=downloader)
        self.assertEqual(fetched, [self.resource.url])
        self.assertTrue(refs.valid_artifact(result, self.resource.semantic_type))
        self.assertEqual(list(self.target.parent.glob("*.partial")), [])

    def test_download_failure_does_not_publish_partial_resource(self):
        def downloader(url, destination):
            Path(destination).write_bytes(b"partial response")
            raise OSError("disconnected")
        with self.assertRaisesRegex(RuntimeError, "Could not prepare gg2-sequences"):
            refs.prepare_resource("gg2-sequences", self.root,
                                  builtin_paths=[], downloader=downloader)
        self.assertFalse(self.target.exists())
        self.assertFalse(any(p.suffix == ".partial" for p in self.target.parent.iterdir()))

    def test_corrupt_or_empty_download_is_rejected(self):
        for payload in (b"html error page", b""):
            with self.subTest(payload=payload):
                def downloader(url, destination):
                    Path(destination).write_bytes(payload)
                with self.assertRaises(RuntimeError):
                    refs.prepare_resource("gg2-sequences", self.root,
                                          builtin_paths=[], downloader=downloader)
                self.assertFalse(self.target.exists())
        empty_artifact = make_artifact(self.root / "empty.qza", payload=b"")
        self.assertFalse(refs.valid_artifact(empty_artifact, self.resource.semantic_type))

    def test_bad_uuid_and_corrupt_payload_are_not_valid_cache(self):
        wrong_uuid = self.root / "wrong-uuid.qza"
        root = str(uuid.uuid4())
        with zipfile.ZipFile(wrong_uuid, "w") as archive:
            archive.writestr(f"{root}/metadata.yaml", f"uuid: {uuid.uuid4()}\ntype: FeatureData[Sequence]\n")
            archive.writestr(f"{root}/data/reference", b"nonempty")
        self.assertFalse(refs.valid_artifact(wrong_uuid, self.resource.semantic_type))
        corrupt = make_artifact(self.root / "corrupt.qza", payload=b"unique-reference-payload")
        raw = corrupt.read_bytes()
        corrupt.write_bytes(raw.replace(b"unique-reference-payload", b"broken-reference-payload"))
        self.assertFalse(refs.valid_artifact(corrupt, self.resource.semantic_type))

    def test_sepp_deployment_reference_is_first_builtin_candidate(self):
        bundle = make_artifact(self.root / "sepp.qza", "SeppReferenceDatabase")
        candidates = list(refs.bundled_candidates("sepp", package_root=self.root,
                          environ={"META2DATA_SEPP_REF": str(bundle)}))
        self.assertEqual(candidates[0], bundle)
        sequence_candidates = list(refs.bundled_candidates("gg2-sequences", package_root=self.root,
                                   environ={"META2DATA_SEPP_REF": str(bundle)}))
        self.assertNotIn(bundle, sequence_candidates)

    def test_concurrent_prepare_downloads_once(self):
        fixture = make_artifact(self.root / "source.qza")
        calls = self.root / "calls.txt"
        context = multiprocessing.get_context("fork")
        start = context.Event()
        replies = context.Queue()
        processes = [context.Process(target=concurrent_prepare,
                     args=(self.root, fixture, calls, start, replies)) for _ in range(4)]
        try:
            for process in processes:
                process.start()
            start.set()
            for process in processes:
                process.join(10)
                self.assertEqual(process.exitcode, 0)
            results = [replies.get(timeout=2) for _ in processes]
            self.assertEqual(results, [(str(self.target), None)] * len(processes))
            self.assertEqual(calls.read_text().splitlines(), ["download"])
        finally:
            for process in processes:
                if process.is_alive():
                    process.terminate()
                    process.join()
            replies.close()


if __name__ == "__main__":
    unittest.main()
