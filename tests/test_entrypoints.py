"""Symlinked launchers must locate their own package without host-home fallback."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
COMMANDS = ('Meta2Data', 'Meta2Data-MetaDL', 'Meta2Data-AmpliconPIP', 'Meta2Data-AmpliconTAXA')
MODULES = ('MetaDL', 'AmpliconPIP', 'AmpliconTAXA')


class EntrypointSymlinkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.links = self.root / 'entry links' / 'bin'
        self.links.mkdir(parents=True)
        self.work = self.root / 'unrelated working directory'
        self.work.mkdir()
        home = self.root / 'empty home'
        home.mkdir()
        self.environment = {key: value for key, value in os.environ.items()
                            if not key.startswith('M2D_PROFILE_')}
        self.environment.update(HOME=str(home), CONDA_PREFIX=str(self.root / 'missing conda'),
                                PREFIX=str(self.root / 'missing prefix'))

    def invoke(self, executable, *args, environment=None):
        result = subprocess.run([str(executable), *args], cwd=self.work,
                                env=environment or self.environment,
                                capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('Usage:', result.stdout)
        self.assertFalse((self.work / 'results').exists(), 'Help must not start an analysis')
        return result

    def test_absolute_symlinks_support_main_and_each_standalone_module(self):
        # No scripts/share tree exists beside these links; all must find the target package.
        self.assertFalse((self.links.parent / 'scripts').exists())
        self.assertFalse((self.links.parent / 'share').exists())
        for name in COMMANDS:
            link = self.links / name
            link.symlink_to(ROOT / 'bin' / name)
            with self.subTest(command=name):
                result = self.invoke(link, '--help')
                self.assertIn('Meta2Data', result.stdout)

    def test_main_symlink_dispatches_without_sibling_module_links(self):
        main = self.links / 'Meta2Data'
        main.symlink_to(ROOT / 'bin' / 'Meta2Data')
        for module in MODULES:
            with self.subTest(module=module):
                self.assertFalse((self.links / ('Meta2Data-' + module)).exists())
                result = self.invoke(main, module, '--help')
                self.assertIn(module, result.stdout)

    def test_relative_symlink_chains_work_when_found_through_path(self):
        intermediate = self.root / 'relative chain'
        intermediate.mkdir()
        for name in COMMANDS:
            target = ROOT / 'bin' / name
            middle = intermediate / name
            middle.symlink_to(os.path.relpath(target, intermediate))
            (self.links / name).symlink_to(os.path.relpath(middle, self.links))
        environment = dict(self.environment,
                           PATH=str(self.links) + os.pathsep + self.environment.get('PATH', '/usr/bin:/bin'))
        for name in COMMANDS:
            with self.subTest(command=name):
                self.invoke(name, '--help', environment=environment)
        for module in MODULES:
            with self.subTest(dispatch=module):
                self.invoke('Meta2Data', module, '--help', environment=environment)

    def test_prefix_share_installation_layout_still_works_through_symlinks(self):
        prefix = self.root / 'installed package'
        install_bin = prefix / 'bin'
        install_bin.mkdir(parents=True)
        share = prefix / 'share' / 'Meta2Data'
        share.mkdir(parents=True)
        (share / 'scripts').symlink_to(ROOT / 'scripts', target_is_directory=True)
        for name in COMMANDS:
            shutil.copy2(ROOT / 'bin' / name, install_bin / name)
            (self.links / name).symlink_to(install_bin / name)
        self.assertFalse((prefix / 'scripts').exists())
        for name in COMMANDS:
            with self.subTest(command=name):
                self.invoke(self.links / name, '--help')
        for module in MODULES:
            with self.subTest(dispatch=module):
                self.invoke(self.links / 'Meta2Data', module, '--help')


if __name__ == '__main__':
    unittest.main()
