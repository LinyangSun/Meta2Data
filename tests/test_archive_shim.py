"""Exercise archive fixtures through the same PATH wrappers as integration runs."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SHIM = Path(__file__).resolve().parent / 'integration' / 'archive_shim.py'


class ArchiveShimTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.bin = self.root / 'fixture-bin'
        self.archive = self.root / 'archive'
        self.archive.mkdir()
        files = []
        for mate in (1, 2):
            name = 'SRR999990001_' + str(mate) + '.fastq.gz'
            payload = gzip.compress(b'@read/1\nACGT\n+\nIIII\n', mtime=0)
            (self.archive / name).write_bytes(payload)
            files.append(dict(path='archive/' + name, md5=hashlib.md5(payload).hexdigest(),
                              url='ftp.sra.ebi.ac.uk/fixtures/' + name))
        self.manifest = self.root / 'public_fixture.json'
        self.data = dict(schema_version=1, synthetic=True,
                         projects={'PRJNA999990001': {'runs': ['SRR999990001']}},
                         runs={'SRR999990001': dict(platform='ILLUMINA', layout='PAIRED', files=files)})
        self.manifest.write_text(json.dumps(self.data))
        self.trace = self.root / 'archive-trace.jsonl'
        subprocess.run([sys.executable, str(SHIM), 'install', '--directory', str(self.bin),
                        '--interpreter', sys.executable], check=True)
        self.env = dict(os.environ, PATH=str(self.bin) + os.pathsep + os.environ.get('PATH', ''),
                        M2D_TEST_ARCHIVE_MANIFEST=str(self.manifest),
                        M2D_TEST_ARCHIVE_TRACE=str(self.trace), M2D_TEST_REAL_PYTHON=sys.executable)

    def run_tool(self, tool, *args):
        return subprocess.run([tool] + list(args), env=self.env, cwd=self.root,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    def events(self):
        return [json.loads(line) for line in self.trace.read_text().splitlines()]

    def test_ena_probe_filereport_and_verified_file_copy(self):
        self.assertEqual(self.run_tool('wget', '-q', '--spider', '--timeout=10',
                                      'https://www.ebi.ac.uk/ena/portal/api/').returncode, 0)
        url = ('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA999990001'
               '&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,library_layout&format=tsv')
        report = self.root / 'report.tsv'
        result = self.run_tool('wget', '-q', '--timeout=30', url, '-O', str(report))
        self.assertEqual(result.returncode, 0, result.stderr)
        with report.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual(rows[0]['run_accession'], 'SRR999990001')
        self.assertEqual(rows[0]['library_layout'], 'PAIRED')
        self.assertEqual(rows[0]['fastq_md5'], ';'.join(item['md5'] for item in self.data['runs']['SRR999990001']['files']))
        url = 'ftp://' + rows[0]['fastq_ftp'].split(';')[0]
        self.assertEqual(self.run_tool('wget', '-q', '--spider', url).returncode, 0)
        destination = self.root / 'download.part'
        result = self.run_tool('wget', '--tries=1', '--timeout=60', url, '-O', str(destination))
        self.assertEqual(result.returncode, 0, result.stderr)
        original = self.root / self.data['runs']['SRR999990001']['files'][0]['path']
        self.assertEqual(destination.read_bytes(), original.read_bytes())
        self.assertEqual([event['action'] for event in self.events()],
                         ['ena_api_probe', 'ena_filereport', 'archive_file_probe', 'archive_file_copy'])

    def test_platform_actions_only_use_declared_runs(self):
        pairs = self.root / 'pairs.tsv'
        pairs.write_text('PRJNA999990001\tSRR999990001\n')
        result = self.run_tool('python', '/unused/py_16s.py', 'batch_get_sequencing_platforms', '--pairs_file', str(pairs))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, 'PRJNA999990001\tILLUMINA\n')
        result = self.run_tool('python3', '/unused/py_16s.py', 'get_sequencing_platform', '--srr_id', 'SRR999990001')
        self.assertEqual(result.stdout, 'ILLUMINA\n')
        result = self.run_tool('python3', '/unused/py_16s.py', 'get_sequencing_platform', '--srr_id', 'SRR404')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, '')
        self.assertIn('Undeclared Run', result.stderr)
        self.assertEqual(self.events()[-1]['status'], 'error')

    def test_unknown_network_and_project_are_rejected_without_output(self):
        destination = self.root / 'must-not-exist'
        for url in ('https://example.org/unlisted.fastq.gz',
                    'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA404'
                    '&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,library_layout&format=tsv'):
            with self.subTest(url=url):
                result = self.run_tool('wget', url, '-O', str(destination))
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(destination.exists())
                self.assertEqual(self.events()[-1]['status'], 'error')

    def test_corrupt_fixture_cannot_fake_successful_download(self):
        entry = self.data['runs']['SRR999990001']['files'][0]
        (self.root / entry['path']).write_bytes(b'corrupt')
        destination = self.root / 'must-not-exist'
        result = self.run_tool('wget', 'https://' + entry['url'], '-O', str(destination))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('MD5 mismatch', result.stderr)
        self.assertFalse(destination.exists())

    def test_other_python_actions_exec_real_interpreter_and_preserve_status(self):
        helper = self.root / 'py_16s.py'
        helper.write_text('import sys\nprint("real helper: " + sys.argv[1])\nsys.exit(7)\n')
        result = self.run_tool('python3', str(helper), 'sanitize_fastq')
        self.assertEqual(result.stdout, 'real helper: sanitize_fastq\n')
        self.assertEqual(result.returncode, 7)
        result = self.run_tool('python', '-c', 'print("real Python")')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, 'real Python\n')
        self.assertTrue(all(event['action'] == 'python_passthrough' for event in self.events()))
        self.assertEqual(sorted(path.name for path in self.bin.iterdir()), ['python', 'python3', 'wget'])


if __name__ == '__main__':
    unittest.main()
