"""Isolated tests for the pinned upstream SRA installer, without downloads."""
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import install_sra_toolkit as installer


class InstallSRATests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.prefix = self.root / "environment"
        self.prefix.mkdir()
        self.archive = self.root / "toolkit.tar.gz"
        self.contents = {
            "bin/sratools.3.4.1": b"test dispatcher\n",
            "bin/prefetch-orig.3.4.1": b"test prefetch\n",
            "bin/fasterq-dump-orig.3.4.1": b"test fasterq\n",
            "schema/insdc/insdc.vschema": b"test schema\n",
        }
        self.write_archive()
        for name, value in {
            "ARCHIVE_SHA256": installer.digest(self.archive),
            "ARCHIVE_MD5": installer.digest(self.archive, "md5"),
            "ARCHIVE_BYTES": self.archive.stat().st_size,
            "PAYLOAD_SHA256": {name: hashlib.sha256(value).hexdigest()
                               for name, value in self.contents.items() if name.startswith("bin/")},
        }.items():
            mocked = patch.object(installer, name, value)
            mocked.start()
            self.addCleanup(mocked.stop)
        platform = patch.object(installer, "check_platform")
        platform.start()
        self.addCleanup(platform.stop)
        self.downloads = []

    def write_archive(self, extra=None):
        with tarfile.open(self.archive, "w:gz") as archive:
            for name, value in self.contents.items():
                entry = tarfile.TarInfo(installer.DIRECTORY + "/" + name)
                entry.size, entry.mode = len(value), 0o755
                archive.addfile(entry, io.BytesIO(value))
            for name in installer.TOOL_FLAGS:
                entry = tarfile.TarInfo(installer.DIRECTORY + "/bin/" + name)
                entry.type, entry.linkname = tarfile.SYMTYPE, "sratools.3.4.1"
                archive.addfile(entry)
            if extra:
                archive.addfile(extra)

    def download(self, destination):
        self.downloads.append(destination)
        shutil.copyfile(self.archive, destination)

    def runner(self, command, **kwargs):
        tool = Path(command[0]).name
        text = (f"{tool} : {installer.VERSION}\n\n" if command[1] == "--version"
                else " ".join(installer.TOOL_FLAGS[tool]))
        return subprocess.CompletedProcess(command, 0, stdout=text, stderr="")

    def install(self, **kwargs):
        return installer.install_toolkit(self.prefix, downloader=self.download,
                                         runner=kwargs.get("runner", self.runner))

    def test_installs_full_layout_exposes_two_commands_and_reuses_verified_bundle(self):
        manifest = self.install()
        evidence = json.loads(manifest.read_text())
        self.assertEqual(evidence["archive"]["sha256"], installer.ARCHIVE_SHA256)
        self.assertEqual(evidence["tools"]["prefetch"]["version"], "3.4.1")
        self.assertEqual(set(path.name for path in (self.prefix / "bin").iterdir()),
                         {"prefetch", "fasterq-dump"})
        for name in installer.TOOL_FLAGS:
            self.assertTrue((self.prefix / "bin" / name).resolve().is_relative_to(self.prefix))
        self.assertTrue((manifest.parent / "schema/insdc/insdc.vschema").is_file())
        self.assertEqual(self.install(), manifest)
        self.assertEqual(len(self.downloads), 1)
        (manifest.parent / "schema/insdc/insdc.vschema").write_text("changed")
        with self.assertRaisesRegex(ValueError, "files changed"):
            self.install()

    def test_existing_commands_are_not_overwritten_or_downloaded_over(self):
        (self.prefix / "bin").mkdir()
        existing = self.prefix / "bin/prefetch"
        existing.write_text("existing command")
        with self.assertRaisesRegex(FileExistsError, "overwrite"):
            self.install()
        self.assertEqual(existing.read_text(), "existing command")
        self.assertFalse(self.downloads)

    def test_bad_download_never_publishes_a_toolkit(self):
        with patch.object(installer, "ARCHIVE_SHA256", "0" * 64):
            with self.assertRaisesRegex(ValueError, "checksum"):
                self.install()
        self.assertFalse((self.prefix / "opt" / installer.DIRECTORY).exists())
        self.assertFalse(list((self.prefix / "opt").glob(".sra-toolkit-*")))
        self.assertFalse(list((self.prefix / "bin").iterdir()))

    def test_escaping_paths_links_and_duplicate_members_are_rejected(self):
        traversal = tarfile.TarInfo(installer.DIRECTORY + "/../../escape")
        symlink = tarfile.TarInfo(installer.DIRECTORY + "/bad-link")
        symlink.type, symlink.linkname = tarfile.SYMTYPE, "../../escape"
        duplicate = tarfile.TarInfo(installer.DIRECTORY + "/bin/sratools.3.4.1")
        for member in (traversal, symlink, duplicate):
            with self.subTest(name=member.name):
                self.write_archive(member)
                with self.assertRaises(ValueError):
                    installer.extract_archive(self.archive, self.root / "unpack")
        self.assertFalse((self.root / "escape").exists())

    def test_wrong_runtime_version_does_not_publish_the_bundle(self):
        def wrong_version(command, **kwargs):
            return subprocess.CompletedProcess(command, 0, stdout="tool : 3.2.1\n", stderr="")
        with self.assertRaisesRegex(ValueError, "version"):
            self.install(runner=wrong_version)
        self.assertFalse((self.prefix / "opt" / installer.DIRECTORY).exists())
        self.assertFalse(list((self.prefix / "bin").iterdir()))

    def test_unregistered_existing_bundle_is_preserved(self):
        target = self.prefix / "opt" / installer.DIRECTORY
        target.mkdir(parents=True)
        (target / "keep").write_text("untouched")
        with self.assertRaisesRegex(FileExistsError, "unregistered"):
            self.install()
        self.assertEqual((target / "keep").read_text(), "untouched")
        self.assertFalse(self.downloads)

    def test_changed_dispatcher_link_is_rejected_before_executing(self):
        manifest = self.install()
        entry = manifest.parent / "bin/prefetch"
        entry.unlink()
        entry.symlink_to("prefetch-orig.3.4.1")
        def must_not_run(*args, **kwargs):
            self.fail("changed installation was executed")
        with self.assertRaisesRegex(ValueError, "links changed"):
            self.install(runner=must_not_run)


if __name__ == "__main__":
    unittest.main()
