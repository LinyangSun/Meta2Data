"""Functions must propagate partial-output failures even from conditional calls."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ChainFailureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.dataset = self.root / 'dataset'
        self.work = self.dataset / 'tmp' / 'step_06_vsearch_cli'
        self.work.mkdir(parents=True)
        self.env = dict(os.environ, dataset_path=str(self.dataset), dataset_name='dataset',
                        SCRIPTS=str(ROOT / 'scripts'), MODE='vsearch',
                        VSEARCH_MIN_FREQUENCY='2', VSEARCH_CLUSTER_IDENTITY='0.97',
                        fastq_path=str(self.root / 'reads'), READ_COUNTS_STATE='',
                        PARTIAL=str(self.root / 'partial'), LATE_AUDIT=str(self.root / 'late_audit'))

    def run_conditional(self, function, stubs):
        script = ('source ' + shlex.quote(str(ROOT / 'scripts' / 'AmpliconFunction.sh')) + '\n'
                  + 'Audit_Fasta() { touch "$LATE_AUDIT"; };\n'
                  + 'Audit_Table() { touch "$LATE_AUDIT"; };\n'
                  + 'Audit_Counts() { touch "$LATE_AUDIT"; };\n'
                  + 'Audit_Fastq() { return 0; };\n' + stubs + '\n'
                  + f'if {function}; then exit 0; else exit "$?"; fi\n')
        return subprocess.run(['bash', '-c', script], env=self.env, capture_output=True, text=True)

    def test_partial_tool_failures_propagate_before_later_audit(self):
        cases = [
            ('Amplicon_Vsearch_ClusterFast97 1', 'vsearch() { echo partial > "$PARTIAL"; return 17; };'),
            ('Amplicon_Vsearch_MapBack vsearch_preprocessed_reads 1',
             'python3() { echo partial > "$PARTIAL"; return 17; };'),
            ('Amplicon_Vsearch_MapBack vsearch_preprocessed_reads 1',
             'python3() { return 0; }; vsearch() { echo partial > "$PARTIAL"; return 17; };'),
            ('Amplicon_Vsearch_ImportResults', 'python3() { echo partial > "$PARTIAL"; return 17; };'),
            ('Amplicon_LS454_FilterLowFreqOTUs 1', 'qiime() { echo partial > "$PARTIAL"; return 17; };'),
            ('Amplicon_LS454_FilterLowFreqOTUs 1',
             'qiime() { [[ "$2" != filter-seqs ]] || { echo partial > "$PARTIAL"; return 17; }; };')]
        for function, stubs in cases:
            with self.subTest(function=function, stubs=stubs):
                result = self.run_conditional(function, stubs)
                self.assertEqual(result.returncode, 17, result.stdout + result.stderr)
                self.assertTrue((self.root / 'partial').exists())
                self.assertFalse((self.root / 'late_audit').exists())

    def test_final_publication_failure_preserves_intermediate_files(self):
        for suffix in ['table-vsearch.qza', 'rep-seqs-vsearch.qza']:
            (self.dataset / ('dataset-' + suffix)).write_text('valid artifact placeholder')
        result = self.run_conditional('Amplicon_Common_FinalFilesCleaning',
            'mv() { echo partial > "$2"; return 17; };')
        self.assertEqual(result.returncode, 17, result.stderr)
        self.assertTrue(self.work.exists())
        self.assertTrue((self.dataset / 'dataset-rep-seqs-vsearch.qza').exists())

    def test_missing_final_representatives_does_not_publish_table(self):
        (self.dataset / 'dataset-table-vsearch.qza').write_text('table')
        result = self.run_conditional('Amplicon_Common_FinalFilesCleaning', '')
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertFalse((self.dataset / 'dataset-vsearch-final-table.qza').exists())
        self.assertTrue(self.work.exists())

    def test_failed_denoising_copy_does_not_remove_intermediates(self):
        denoise = self.dataset / 'tmp' / 'step_05_denoise'
        denoise.mkdir()
        (denoise / 'dataset-table-denoising.qza').write_text('table')
        (denoise / 'dataset-rep-seqs-denoising.qza').write_text('sequences')
        result = self.run_conditional('Amplicon_Common_FinalFilesCleaning',
            'cp() { echo partial > "$2"; return 17; };')
        self.assertEqual(result.returncode, 17, result.stderr)
        self.assertTrue(denoise.exists())


if __name__ == '__main__':
    unittest.main()
