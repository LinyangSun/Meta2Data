"""Regression tests for minimal runtime locking and installer isolation."""
from contextlib import redirect_stdout
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import create_conda_environment as installer


class EnvironmentInstallerTests(unittest.TestCase):
    def test_minimal_environment_has_supported_channels_and_no_conflicting_sra_packages(self):
        channels, pins = installer.read_freeze(ROOT / "env.yml")
        self.assertEqual(channels, [
            "conda-forge", "bioconda",
            "https://packages.qiime2.org/qiime2/2026.7/qiime2/released"])
        self.assertIn("q2-feature-classifier=2026.7.0", pins)
        self.assertIn("scikit-learn=1.7.1", pins)
        specifications = installer.read_lock(ROOT / "environment-linux-64.lock",
                                             ROOT / "env.yml", channels, pins)
        for package in ("q2-fondue", "q2-vsearch", "sra-tools", "ossuuid"):
            self.assertFalse(any(pin.startswith(package + "=") for pin in pins))
            self.assertFalse(any("/" + package + "-" in line for line in specifications))

    def test_dry_run_uses_explicit_lock_and_never_installs_sra(self):
        channels, pins = installer.read_freeze(ROOT / "env.yml")
        def fake_run(command, *, env):
            self.assertEqual(command[:2], ["/conda/bin/conda", "create"])
            self.assertIn("--override-channels", command)
            self.assertIn("--no-default-packages", command)
            self.assertIn("--dry-run", command)
            self.assertNotIn("--format", command)
            selected = [command[index + 1] for index, item in enumerate(command)
                        if item == "--channel"]
            self.assertEqual(selected, channels)
            specs = Path(command[command.index("--file") + 1])
            self.assertEqual(specs, ROOT / "environment-linux-64.lock")
            self.assertEqual(env["CONDA_CHANNEL_PRIORITY"], "flexible")
            return subprocess.CompletedProcess(command, 0)
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(installer.shutil, "which", return_value="/conda/bin/conda"), \
                    mock.patch.object(installer.subprocess, "run", side_effect=fake_run) as run, \
                    mock.patch.dict(os.environ, CONDA_CHANNEL_PRIORITY="strict"), \
                    redirect_stdout(io.StringIO()):
                self.assertEqual(installer.create(Path(directory) / "env", "conda", dry_run=True), 0)
                run.assert_called_once()
                self.assertEqual(os.environ["CONDA_CHANNEL_PRIORITY"], "strict")

    def test_existing_environment_is_preserved_before_conda_is_called(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory)
            sentinel = prefix / "existing-data"
            sentinel.write_text("preserve me")
            with mock.patch.object(installer.subprocess, "run") as run:
                with self.assertRaisesRegex(ValueError, "absent or empty"):
                    installer.create(prefix, "conda")
                run.assert_not_called()
            self.assertEqual(sentinel.read_text(), "preserve me")

    def test_nonflat_or_ambiguous_requirements_are_rejected(self):
        cases = (
            ("- pip:\n  - requests\n", "single-version"),
            ("- numpy=2.4.2\n- numpy=1.26.4\n", "duplicate package"),
        )
        with tempfile.TemporaryDirectory() as directory:
            freeze = Path(directory) / "env.yml"
            for dependencies, message in cases:
                with self.subTest(dependencies=dependencies):
                    freeze.write_text("channels:\n- conda-forge\n- nodefaults\n"
                                      "dependencies:\n" + dependencies)
                    with self.assertRaisesRegex(ValueError, message):
                        installer.read_freeze(freeze)

    def test_installer_failure_is_returned_without_retry_or_clobber(self):
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(installer.shutil, "which", return_value="/conda/bin/conda"), \
                    mock.patch.object(installer.subprocess, "run",
                                      return_value=subprocess.CompletedProcess([], 1)) as run, \
                    redirect_stdout(io.StringIO()):
                self.assertEqual(installer.create(Path(directory) / "env", "conda"), 1)
                run.assert_called_once()
                self.assertNotIn("--clobber", run.call_args.args[0])

    def test_successful_conda_creation_requires_successful_sra_setup(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "env"
            with mock.patch.object(installer.shutil, "which", return_value="/conda/bin/conda"), \
                    mock.patch.object(installer.subprocess, "run", side_effect=[
                        subprocess.CompletedProcess([], 0), subprocess.CompletedProcess([], 7)]) as run, \
                    redirect_stdout(io.StringIO()):
                self.assertEqual(installer.create(prefix, "conda"), 7)
                self.assertEqual(run.call_count, 2)
                self.assertEqual(run.call_args.args[0], [
                    str(prefix / "bin/python"), str(ROOT / "scripts/install_sra_toolkit.py"),
                    "--prefix", str(prefix)])

    def test_stale_or_modified_lock_is_rejected(self):
        channels, pins = installer.read_freeze(ROOT / "env.yml")
        original = (ROOT / "environment-linux-64.lock").read_text()
        cases = (
            (original.replace("# env.yml sha256:", "# outdated hash:"), "does not match"),
            (original.replace("# release constraints sha256:", "# outdated constraints:"),
             "does not match release constraints"),
            (original.replace("https://conda.anaconda.org/bioconda/", "https://example.org/"),
             "Unsupported package source"),
            ("\n".join(line for line in original.splitlines() if "/python-3.12.13-" not in line),
             "does not provide python"),
            (original.replace("/bibtexparser-1.4.4-", "/bibtexparser-2.1.0-"),
             "differs from the exact release"),
        )
        with tempfile.TemporaryDirectory() as directory:
            lock = Path(directory) / "explicit.txt"
            for contents, message in cases:
                with self.subTest(message=message):
                    lock.write_text(contents)
                    with self.assertRaisesRegex(ValueError, message):
                        installer.read_lock(lock, ROOT / "env.yml", channels, pins)


if __name__ == "__main__":
    unittest.main()
