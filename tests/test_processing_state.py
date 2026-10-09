"""Regressions for retaining denoised sequences and invalidating old results."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import pip_state


class ProcessingStateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dataset = Path(self.temp.name) / 'PRJNA123'
        self.dataset.mkdir()
        (self.dataset / 'PRJNA123_sra.txt').write_text('SRR123\n')
        self.environment = patch.dict(os.environ, {}, clear=True)
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.code = patch('pip_state.code_fingerprint', return_value='current-code')
        self.code.start()
        self.addCleanup(self.code.stop)

    def prepare(self, forward='', reverse=''):
        pip_state.prepare(self.dataset, 'dada2', forward=forward, reverse=reverse)
        return json.loads((self.dataset / '.checkpoint.json').read_text())

    def results(self):
        paths = [self.dataset / f'PRJNA123-dada2-final-{suffix}.qza'
                 for suffix in ('table', 'rep-seqs')]
        for path in paths:
            path.write_bytes(b'completed-result')
        return paths

    def test_unchanged_primers_reuse_completed_results(self):
        before = self.prepare('ACGT', 'TGCA')
        paths = self.results()
        after = self.prepare('ACGT', 'TGCA')
        self.assertEqual(before['fingerprint'], after['fingerprint'])
        self.assertTrue(all(path.exists() for path in paths))
        self.assertEqual(after['primer_source'], 'explicit')

    def test_changed_project_primers_invalidate_processing_but_keep_raw(self):
        self.prepare('ACGT', 'TGCA')
        paths = self.results()
        raw = self.dataset / 'downloaded_fastq'
        raw.mkdir()
        (raw / 'SRR123.fastq').write_text('@r\nACGT\n+\nIIII\n')
        (self.dataset / 'tmp').mkdir()
        self.prepare('ACGA', 'TGCA')
        self.assertTrue(all(not path.exists() for path in paths))
        self.assertFalse((self.dataset / 'tmp').exists())
        self.assertTrue((raw / 'SRR123.fastq').exists())

    def test_reverse_primer_alone_changes_fingerprint(self):
        before = self.prepare('ACGT', 'TGCA')
        paths = self.results()
        after = self.prepare('ACGT', 'TGCT')
        self.assertNotEqual(before['fingerprint'], after['fingerprint'])
        self.assertTrue(all(not path.exists() for path in paths))

    def test_automatic_to_explicit_does_not_reuse_results(self):
        self.assertEqual(self.prepare()['primer_source'], 'automatic')
        paths = self.results()
        state = self.prepare('ACGT')
        self.assertEqual(state['primer_source'], 'explicit')
        self.assertTrue(all(not path.exists() for path in paths))

    def test_old_processing_code_cannot_reuse_truncated_pacbio_results(self):
        with patch('pip_state.code_fingerprint', return_value='legacy-v3v4-code'):
            self.prepare('ACGT')
        paths = self.results()
        self.prepare('ACGT')
        self.assertTrue(all(not path.exists() for path in paths))

    def test_changing_one_project_keeps_other_project_results(self):
        self.prepare('ACGT')
        other = Path(self.temp.name) / 'PRJNA456'
        other.mkdir()
        (other / 'PRJNA456_sra.txt').write_text('SRR456\n')
        pip_state.prepare(other, 'dada2', forward='TGCA')
        result = other / 'PRJNA456-dada2-final-table.qza'
        result.write_bytes(b'other-project')
        self.prepare('AAAA')
        pip_state.prepare(other, 'dada2', forward='TGCA')
        self.assertEqual(result.read_bytes(), b'other-project')

    def test_ion_report_is_kept_on_reuse_and_removed_when_settings_change(self):
        self.prepare()
        report = self.dataset / 'ion_quality-dada2.json'
        report.write_text('{"maxee": 5.123826426002262}')
        other = self.dataset / 'ion_quality-vsearch.json'
        other.write_text('{"maxee": 2}')
        self.prepare()
        self.assertTrue(report.exists())
        with patch.dict(os.environ, {'DADA2_ION_MAXEE': '2'}):
            self.prepare()
        self.assertFalse(report.exists())
        self.assertTrue(other.exists())


    def test_454_reports_follow_parameter_invalidation(self):
        pip_state.prepare(self.dataset, 'vsearch')
        reports = [self.dataset / name for name in
                   ('ls454_quality-vsearch.json', 'ls454_members-vsearch.json')]
        for path in reports:
            path.write_text('{}')
        pip_state.prepare(self.dataset, 'vsearch')
        self.assertTrue(all(path.exists() for path in reports))
        with patch.dict(os.environ, {'LS454_LENGTH_FRACTION': '0.7'}):
            pip_state.prepare(self.dataset, 'vsearch')
        self.assertTrue(all(not path.exists() for path in reports))


    def test_454_n_limit_change_invalidates_results_and_reports(self):
        pip_state.prepare(self.dataset, 'vsearch')
        outputs = [self.dataset / f'PRJNA123-vsearch-final-{suffix}.qza'
                   for suffix in ('table', 'rep-seqs')]
        reports = [self.dataset / name for name in
                   ('ls454_quality-vsearch.json', 'ls454_members-vsearch.json')]
        for path in outputs + reports:
            path.write_text('old result')
        with patch.dict(os.environ, {'LS454_MAX_N': '0'}):
            pip_state.prepare(self.dataset, 'vsearch')
        self.assertTrue(all(not path.exists() for path in outputs + reports))


class PacBioOutputTests(unittest.TestCase):
    def test_full_length_representatives_and_table_are_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            work = Path(temp)
            sequence = ''.join((ROOT / 'docs/ecoli_16S_J01859.fasta').read_text().splitlines()[1:])
            self.assertGreater(len(sequence), 1400)
            rep = work / 'full-length.qza'
            table = work / 'table.qza'
            with zipfile.ZipFile(rep, 'w') as out:
                out.writestr('fixture/data/dna-sequences.fasta', '>feature1\n' + sequence + '\n')
            with zipfile.ZipFile(table, 'w') as out:
                out.writestr('fixture/data/table.json', '{"feature1": {"sample1": 10}}')
            dataset = work / 'pacbio'
            dataset.mkdir()
            env = os.environ.copy()
            env.update(FUNCTIONS=str(ROOT / 'scripts/AmpliconFunction.sh'),
                       FULL_REP=str(rep), FULL_TABLE=str(table), dataset_path=str(dataset),
                       MODE='dada2', primer_front='AGAGTTTGATCMTGGCTCAG', primer_adapter='',
                       DADA2_PACBIO_MIN_LENGTH='1000', DADA2_PACBIO_MAX_LENGTH='1600', cpu='1')
            script = '''
source "$FUNCTIONS"
Audit_Dada2() { :; }
qiime() {
    [[ "$1" == dada2 && "$2" == denoise-ccs ]] || return 73
    shift 2
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --o-representative-sequences) cp "$FULL_REP" "$2"; shift 2 ;;
            --o-table) cp "$FULL_TABLE" "$2"; shift 2 ;;
            *) shift ;;
        esac
    done
}
Amplicon_Pacbio_DenosingDada2
Amplicon_Common_FinalFilesCleaning
'''
            result = subprocess.run(['bash', '-c', script], env=env, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertEqual((dataset / 'pacbio-dada2-final-rep-seqs.qza').read_bytes(), rep.read_bytes())
            self.assertEqual((dataset / 'pacbio-dada2-final-table.qza').read_bytes(), table.read_bytes())


if __name__ == '__main__':
    unittest.main()
