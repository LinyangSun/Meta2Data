"""Small real FASTQ checks around local staging and guarded fastp publication."""
import gzip
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import fastp_checked
import read_layout


class GuardedLocalFastqTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        source_dir = self.work / 'local_data'
        source_dir.mkdir()
        self.source = source_dir / 'sample.FASTQ.GZ'
        with gzip.open(self.source, 'wt') as stream:
            stream.write('@read\nACGT\n+\nIIII\n')
        self.original = self.source.read_bytes()
        staged = self.work / 'staged'
        read_layout.stage(read_layout.discover(source_dir), staged)
        self.args = SimpleNamespace(
            sample_id='sample', in1=str(staged / 'sample.fastq.gz'), in2=None,
            out1=str(self.work / 'processed/sample.fastq.gz'), out2=None,
            threads=1, audit_dir=str(self.work / 'audit'), work_dir=str(self.work / 'scratch'),
            db_manifest=str(self.work / 'reference.json'), cache_dir=str(self.work / 'cache'),
            report_json=str(self.work / 'fastp.json'), report_html=str(self.work / 'fastp.html'),
            min_length=50, min_identity=98, min_coverage=95)

    def external_tool(self, command, **kwargs):
        def value(option):
            return command[command.index(option) + 1]
        if command[0] == 'fastp':
            shutil.copyfile(value('-i'), value('-o'))
            Path(value('-j')).write_text(json.dumps({
                'summary': {'after_filtering': {'total_reads': 1, 'total_bases': 4}}}))
            Path(value('-h')).write_text('fastp report')
        elif Path(command[1]).name == 'fastp_guard.py':
            Path(value('--output-dir'), 'decision.json').write_text(json.dumps({
                'action': 'accept', 'status': 'no_candidate'}))
        else:
            self.fail(f'Unexpected external command: {command}')

    def test_guard_accepts_uppercase_gzip_source_after_local_staging(self):
        self.assertEqual(list(fastp_checked.records(self.source)), [('@read', 4)])
        with patch.object(fastp_checked.subprocess, 'run', side_effect=self.external_tool) as invoked:
            fastp_checked.process(self.args)
        self.assertEqual(invoked.call_count, 2)
        with gzip.open(self.args.out1, 'rt') as output:
            self.assertEqual(output.read(), '@read\nACGT\n+\nIIII\n')
        decision = json.loads((self.work / 'audit/decision.json').read_text())
        self.assertEqual(decision['publication_status'], 'complete')
        self.assertEqual(decision['accepted_counts']['pipeline_reads'], 1)
        self.assertEqual(self.source.read_bytes(), self.original)

    def test_guard_still_rejects_gzip_to_plain_output_before_running_tools(self):
        self.args.out1 = str(self.work / 'processed/sample.fastq')
        with patch.object(fastp_checked.subprocess, 'run') as invoked:
            with self.assertRaisesRegex(ValueError, 'compression suffixes must agree'):
                fastp_checked.process(self.args)
        invoked.assert_not_called()
        self.assertFalse((self.work / 'processed').exists())
        self.assertFalse((self.work / 'audit').exists())
        self.assertEqual(self.source.read_bytes(), self.original)


if __name__ == '__main__':
    unittest.main()
