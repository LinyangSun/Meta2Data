"""Online SRA biological read numbers must not relax local pairing checks."""
import contextlib
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import download_integrity
import read_layout


class OnlineReadLayoutTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dataset = Path(self.temp.name) / 'PRJNA283199'
        self.input = self.dataset / 'ori_fastq'
        self.input.mkdir(parents=True)
        self.output = self.dataset / 'read_layout.json'
        self.sample = 'PRJNA283199_SRR2011463'
        self.record = '@SRR2011463.1\nACGT\n+\nIIII\n'

    def fastq(self, suffix):
        path = self.input / (self.sample + suffix + '.fastq.gz')
        with gzip.open(path, 'wt') as stream:
            stream.write(self.record)
        return path

    def normalize(self, expected=None):
        argv = ['read_layout.py', 'normalize', '--input', str(self.input),
                '--output', str(self.output)]
        if expected:
            argv += ['--expected-layout', expected]
        with patch.object(sys, 'argv', argv):
            read_layout.main()

    def assert_normalized_single(self, source):
        rows = json.loads(self.output.read_text())
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['sample'], self.sample)
        self.assertEqual(rows[0]['layout'], 'SE')
        self.assertEqual(rows[0]['r2'], '')
        archived = self.dataset / 'downloaded_fastq' / source.name
        staged = self.input / (self.sample + '.fastq.gz')
        self.assertEqual(Path(rows[0]['r1']), archived)
        self.assertTrue(staged.is_symlink())
        self.assertEqual(staged.resolve(), archived.resolve())
        with gzip.open(staged, 'rt') as stream:
            self.assertEqual(stream.read(), self.record)

    def assert_rejected_without_moving(self, expected=None, message='FASTQ input error'):
        original = {p.name: p.read_bytes() for p in self.input.iterdir()}
        error = io.StringIO()
        with contextlib.redirect_stderr(error), self.assertRaises(SystemExit) as caught:
            self.normalize(expected)
        self.assertEqual(caught.exception.code, 2)
        self.assertIn(message, error.getvalue())
        self.assertEqual({p.name: p.read_bytes() for p in self.input.iterdir()}, original)
        self.assertFalse((self.dataset / 'downloaded_fastq').exists())
        self.assertFalse(self.output.exists())

    def test_known_single_end_sra_read_two_is_staged_without_mate_suffix(self):
        # Real 454 SRR2011463 has technical read 1 and biological read 2.
        source = self.fastq('_2')
        self.normalize('SE')
        self.assert_normalized_single(source)

    def test_known_single_end_sra_read_one_remains_supported(self):
        source = self.fastq('_1')
        self.normalize('SE')
        self.assert_normalized_single(source)

    def test_existing_online_read_one_behavior_remains_supported(self):
        source = self.fastq('_1')
        self.normalize()
        self.assert_normalized_single(source)

    def test_known_single_end_canonical_filename_remains_supported(self):
        source = self.fastq('')
        self.normalize('SE')
        self.assert_normalized_single(source)

    def test_online_read_two_requires_known_single_end_context(self):
        self.fastq('_2')
        self.assert_rejected_without_moving(message='Unmatched or ambiguous')

    def test_known_paired_input_cannot_accept_read_two_only(self):
        self.fastq('_2')
        self.assert_rejected_without_moving('PE')

    def test_known_paired_input_cannot_accept_read_one_only(self):
        self.fastq('_1')
        self.assert_rejected_without_moving('PE', 'Expected PE FASTQ layout, found SE')

    def test_complete_online_pair_keeps_both_mates_and_sample_id(self):
        first, second = self.fastq('_1'), self.fastq('_2')
        self.normalize()
        rows = json.loads(self.output.read_text())
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['sample'], self.sample)
        self.assertEqual(rows[0]['layout'], 'PE')
        for direction, source in [('r1', first), ('r2', second)]:
            self.assertEqual(Path(rows[0][direction]).name, source.name)
            self.assertTrue((self.input / source.name).is_symlink())
            with gzip.open(self.input / source.name, 'rt') as stream:
                self.assertEqual(stream.read(), self.record)

    def test_known_single_end_platform_rejects_complete_paired_files(self):
        self.fastq('_1')
        self.fastq('_2')
        self.assert_rejected_without_moving('SE', 'Expected SE FASTQ layout, found PE')

    def test_local_numbered_singletons_remain_rejected(self):
        for suffix in ('_1', '_2'):
            with self.subTest(suffix=suffix):
                source = self.fastq(suffix)
                with self.assertRaisesRegex(ValueError, 'Unmatched or ambiguous'):
                    read_layout.discover(self.input)
                source.unlink()

    def test_archive_known_paired_guard_rejects_missing_either_mate(self):
        # This check also runs before NCBI publishes its complete-download
        # manifest, so the platform-specific SE allowance cannot bypass it.
        for suffix in ('_1.fastq', '_2.fastq'):
            with self.subTest(suffix=suffix):
                with self.assertRaisesRegex(ValueError, 'declares PAIRED'):
                    download_integrity.selected([('SRR2011463' + suffix, '')],
                                                'SRR2011463', 'PAIRED')


if __name__ == '__main__':
    unittest.main()
