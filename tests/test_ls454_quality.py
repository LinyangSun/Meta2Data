"""454 preprocessing uses per-sample read-weighted length limits, then N removal."""
import gzip
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import ls454_quality


class LS454QualityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'input'
        self.source.mkdir()
        self.output = self.root / 'output'
        self.report = self.root / 'quality.json'

    def write(self, name, sequences, quality='I', newline='\n', final_newline=True):
        path = self.source / name
        lines = []
        for index, sequence in enumerate(sequences):
            lines.extend([f'@read_{index} original description', sequence,
                          '+', quality * len(sequence)])
        data = newline.join(lines)
        if lines and final_newline:
            data += newline
        if name.lower().endswith('.gz'):
            path.write_bytes(gzip.compress(data.encode('ascii')))
        else:
            path.write_bytes(data.encode('ascii'))
        return path

    def run_filter(self, fraction=0.5, max_n=1):
        return ls454_quality.preprocess(self.source, self.output, self.report, fraction, max_n)

    def text(self, name):
        path = self.output / name
        opener = gzip.open if name.endswith('.gz') else open
        with opener(path, 'rt', encoding='ascii', newline='') as stream:
            return stream.read()

    def assert_counts(self, result, **expected):
        expected.setdefault('kept_with_n', 0)
        self.assertEqual(result['totals'], expected)
        self.assertEqual(expected['input'], expected['removed_short']
                         + expected['removed_n'] + expected['kept'])

    def test_boundary_256_is_kept_and_255_is_removed(self):
        self.write('sample.fastq', ['A' * 255, 'C' * 256] + ['G' * 512] * 3)
        result = self.run_filter()
        sample = result['files'][0]
        self.assertEqual(sample['median_length'], 512)
        self.assertEqual(sample['min_length'], 256)
        self.assertNotIn('@read_0 ', self.text('sample.fastq'))
        self.assertIn('@read_1 ', self.text('sample.fastq'))
        self.assert_counts(result, input=5, removed_short=1, removed_n=0,
                           removed_both=0, kept=4)

    def test_even_median_is_fractional_and_threshold_uses_ceiling(self):
        self.write('sample.fastq', ['A' * 511, 'C' * 514])
        result = self.run_filter()
        self.assertEqual(result['files'][0]['median_length'], 512.5)
        self.assertEqual(result['files'][0]['min_length'], 257)

    def test_median_uses_every_read_not_unique_sequence_lengths(self):
        self.write('sample.fastq', ['A' * 100] * 7 + ['C' * 400, 'G' * 900])
        result = self.run_filter()
        self.assertEqual(result['files'][0]['median_length'], 100)
        self.assertEqual(result['files'][0]['min_length'], 50)

    def test_each_sample_gets_its_own_threshold(self):
        self.write('long.fastq', ['A' * 255] + ['A' * 512] * 3)
        self.write('short.fastq', ['C' * 100] + ['C' * 200] * 3)
        result = self.run_filter()
        samples = {entry['sample']: entry for entry in result['files']}
        self.assertEqual(samples['long']['min_length'], 256)
        self.assertEqual(samples['long']['kept'], 3)
        self.assertEqual(samples['short']['min_length'], 100)
        self.assertEqual(samples['short']['kept'], 4)

    def test_zero_n_override_keeps_primary_reasons_mutually_exclusive(self):
        self.write('sample.fastq', ['N' * 10, 'n' + 'A' * 99, 'N' + 'C' * 99,
                                    'G' * 100, 'T' * 100])
        result = self.run_filter(max_n=0)
        self.assertEqual(result['max_n'], 0)
        self.assertEqual(result['files'][0]['min_length'], 50)
        self.assert_counts(result, input=5, removed_short=1, removed_n=2,
                           removed_both=1, kept=2)
        self.assertNotIn('@read_0 ', self.text('sample.fastq'))
        self.assertNotIn('@read_1 ', self.text('sample.fastq'))
        self.assertNotIn('@read_2 ', self.text('sample.fastq'))

    def test_default_one_n_limit_counts_both_cases_and_keeps_the_boundary(self):
        self.write('sample.fastq', ['N' + 'A' * 99, 'n' + 'C' * 99,
                                    'NN' + 'G' * 98, 'Nn' + 'T' * 98,
                                    'nn' + 'A' * 98, 'C' * 100])
        result = self.run_filter()
        self.assertEqual(result['max_n'], 1)
        self.assertEqual(result['filter_order'], ['short_length', 'N_count_exceeds_max_n'])
        self.assert_counts(result, input=6, removed_short=0, removed_n=3,
                           removed_both=0, kept=3, kept_with_n=2)
        self.assertEqual(result['files'][0]['kept_with_n'], 2)
        output = self.text('sample.fastq')
        for index in (0, 1, 5):
            self.assertIn(f'@read_{index} ', output)
        for index in (2, 3, 4):
            self.assertNotIn(f'@read_{index} ', output)
        self.assertIn('N' + 'A' * 99, output)
        self.assertIn('n' + 'C' * 99, output)

    def test_length_has_priority_and_removed_both_requires_excess_n(self):
        self.write('sample.fastq', ['N' + 'A' * 9, 'Nn' + 'C' * 8,
                                    'N' + 'G' * 99, 'nn' + 'T' * 98,
                                    'A' * 100, 'C' * 100, 'G' * 100])
        result = self.run_filter()
        self.assertEqual(result['files'][0]['min_length'], 50)
        self.assert_counts(result, input=7, removed_short=2, removed_n=1,
                           removed_both=1, kept=4, kept_with_n=1)

    def test_invalid_n_limits_cannot_create_or_change_outputs(self):
        source = self.write('sample.fastq', ['A' * 100])
        before = source.read_bytes()
        for value in (-1, 1.0, 0.5, True, False, None, '1', float('nan'), float('inf')):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'nonnegative integer'):
                self.run_filter(max_n=value)
        self.assertFalse(self.output.exists())
        self.assertEqual(source.read_bytes(), before)

    def test_nonnegative_integer_n_limit_can_be_overridden(self):
        self.write('sample.fastq', ['Nn' + 'A' * 98, 'NNn' + 'C' * 97, 'G' * 100])
        result = self.run_filter(max_n=2)
        self.assertEqual(result['max_n'], 2)
        self.assert_counts(result, input=3, removed_short=0, removed_n=1,
                           removed_both=0, kept=2, kept_with_n=1)

    def test_n_reads_are_included_in_the_median(self):
        self.write('sample.fastq', ['A' * 100, 'N' * 600, 'n' * 600])
        with self.assertRaisesRegex(ValueError, 'All 454 reads'):
            self.run_filter()
        report = json.loads(self.report.read_text())
        self.assertEqual(report['files'][0]['median_length'], 600)
        self.assertEqual(report['files'][0]['min_length'], 300)
        self.assert_counts(report, input=3, removed_short=1, removed_n=2,
                           removed_both=0, kept=0)

    def test_phred_scores_do_not_filter_or_modify_the_reads(self):
        original = self.write('sample.fastq', ['ACGT' * 100], quality='!')
        before = original.read_bytes()
        result = self.run_filter()
        self.assertFalse(result['quality_filter'])
        self.assertEqual(result['files'][0]['kept'], 1)
        self.assertEqual((self.output / 'sample.fastq').read_bytes(), before)
        self.assertEqual(original.read_bytes(), before)

    def test_compressed_fq_suffix_is_normalized_without_changing_sample_id(self):
        source = self.write('SRR2011463.FQ.GZ', ['acgt' * 100], newline='\r\n',
                            final_newline=False)
        result = self.run_filter()
        self.assertEqual(result['files'][0]['sample'], 'SRR2011463')
        self.assertEqual(self.text('SRR2011463.fastq.gz'),
                         gzip.decompress(source.read_bytes()).decode('ascii'))

    def test_empty_file_beside_nonempty_sample_is_reported_and_retained_as_empty(self):
        self.write('empty.fastq.gz', [])
        self.write('valid.fastq', ['A' * 100])
        result = self.run_filter()
        sample = next(entry for entry in result['files'] if entry['sample'] == 'empty')
        self.assertEqual(sample['status'], 'empty_input')
        self.assertIsNone(sample['median_length'])
        self.assertIsNone(sample['min_length'])
        self.assertEqual(self.text('empty.fastq.gz'), '')
        self.assert_counts(result, input=1, removed_short=0, removed_n=0,
                           removed_both=0, kept=1)

    def test_one_all_filtered_sample_does_not_discard_other_samples(self):
        self.write('all_n.fastq', ['N' * 100])
        self.write('valid.fastq', ['A' * 100])
        result = self.run_filter()
        self.assertEqual(result['status'], 'completed')
        self.assertEqual(result['files'][0]['status'], 'all_filtered')
        self.assertEqual(self.text('all_n.fastq'), '')
        self.assertEqual(result['totals']['kept'], 1)

    def test_all_empty_inputs_fail_with_complete_saved_report(self):
        self.write('empty.fastq', [])
        self.write('empty2.fastq.gz', [])
        with self.assertRaisesRegex(ValueError, 'All input FASTQ files are empty'):
            self.run_filter()
        report = json.loads(self.report.read_text())
        self.assertEqual(report['status'], 'empty_input')
        self.assertEqual(len(report['files']), 2)
        self.assertEqual(report['totals']['input'], 0)

    def test_no_input_files_fail_with_saved_report(self):
        with self.assertRaisesRegex(ValueError, 'No FASTQ files'):
            self.run_filter()
        self.assertEqual(json.loads(self.report.read_text())['status'], 'failed')
        self.assertFalse(self.output.exists())

    def test_malformed_inputs_fail_before_publishing_any_fastq(self):
        valid = self.write('a_valid.fastq', ['A' * 100])
        before = valid.read_bytes()
        invalid = self.source / 'z_bad.fastq'
        records = {
            'header': b'read\nACGT\n+\nIIII\n',
            'missing_id': b'@\nACGT\n+\nIIII\n',
            'blank_id': b'@ description\nACGT\n+\nIIII\n',
            'separator': b'@read\nACGT\nminus\nIIII\n',
            'unequal_quality': b'@read\nACGT\n+\nIII\n',
            'missing_quality': b'@read\nACGT\n+\n',
            'empty_read': b'@read\n\n+\n\n',
            'sequence_space': b'@read\nAC T\n+\nIIII\n',
            'quality_space': b'@read\nACGT\n+\nIII \n',
            'non_ascii': b'@read\nACGT\n+\nIII\xff\n',
        }
        for label, record in records.items():
            invalid.write_bytes(record)
            with self.subTest(label=label), self.assertRaises(ValueError):
                self.run_filter()
            self.assertFalse(self.output.exists())
            self.assertEqual(valid.read_bytes(), before)
            self.assertEqual(json.loads(self.report.read_text())['status'], 'failed')

    def test_truncated_gzip_fails_without_output_publication(self):
        path = self.write('sample.fastq.gz', ['A' * 100])
        path.write_bytes(path.read_bytes()[:-8])
        with self.assertRaises((ValueError, OSError, EOFError)):
            self.run_filter()
        self.assertFalse(self.output.exists())

    def test_invalid_length_fractions_are_rejected(self):
        self.write('sample.fastq', ['A' * 100])
        for value in (0, -0.1, 1.1, True, None, '0.5', float('nan'), float('inf')):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.run_filter(value)
        self.assertFalse(self.output.exists())

    def test_custom_fraction_one_keeps_the_exact_median_boundary(self):
        self.write('sample.fastq', ['A' * 99, 'C' * 100, 'G' * 101])
        result = self.run_filter(1.0)
        self.assertEqual(result['files'][0]['min_length'], 100)
        self.assertEqual(result['totals']['kept'], 2)

    def test_input_and_output_directory_overlap_is_rejected(self):
        source = self.write('sample.fastq', ['A' * 100])
        before = source.read_bytes()
        for output in (self.source, self.source / 'nested', self.root):
            with self.subTest(output=output), self.assertRaisesRegex(ValueError, 'overlap'):
                ls454_quality.preprocess(self.source, output, self.report)
        self.assertEqual(source.read_bytes(), before)

    def test_symlinked_output_into_input_is_rejected(self):
        self.write('sample.fastq', ['A' * 100])
        self.output.symlink_to(self.source, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'overlap'):
            self.run_filter()

    def test_input_symlink_inside_output_is_rejected(self):
        self.output.mkdir()
        actual = self.output / 'actual.fastq'
        actual.write_text('@read\nACGT\n+\nIIII\n')
        (self.source / 'alias.fastq').symlink_to(actual)
        with self.assertRaisesRegex(ValueError, 'resolves inside'):
            self.run_filter()

    def test_report_cannot_overwrite_an_input_or_output_fastq(self):
        source = self.write('sample.fastq', ['A' * 100])
        before = source.read_bytes()
        for report in (source, self.source / 'report.json', self.output / 'sample.fastq'):
            with self.subTest(report=report), self.assertRaises(ValueError):
                ls454_quality.preprocess(self.source, self.output, report)
        self.assertEqual(source.read_bytes(), before)

    def test_existing_outputs_are_not_deleted_or_mixed_with_new_results(self):
        self.write('sample.fastq', ['A' * 100])
        self.output.mkdir()
        old = self.output / 'old.fastq'
        old.write_text('old result')
        with self.assertRaisesRegex(ValueError, 'empty'):
            self.run_filter()
        self.assertEqual(old.read_text(), 'old result')

    def test_duplicate_sample_names_and_resolved_sources_are_rejected(self):
        source = self.write('sample.fastq', ['A' * 100])
        duplicate = self.write('sample.fq.gz', ['C' * 200])
        with self.assertRaisesRegex(ValueError, 'duplicate sample'):
            self.run_filter()
        duplicate.unlink()
        alias = self.source / 'alias.fastq'
        alias.symlink_to(source)
        with self.assertRaisesRegex(ValueError, 'Duplicate FASTQ'):
            self.run_filter()
        alias.unlink()
        os.link(source, alias)
        with self.assertRaisesRegex(ValueError, 'Duplicate FASTQ'):
            self.run_filter()

    def test_non_fastq_and_nested_files_are_ignored(self):
        self.write('sample.fastq', ['A' * 100])
        (self.source / 'primer_info.json').write_text('{}')
        nested = self.source / 'nested'
        nested.mkdir()
        (nested / 'broken.fastq').write_text('not a FASTQ')
        result = self.run_filter()
        self.assertEqual(len(result['files']), 1)
        self.assertEqual(result['totals']['kept'], 1)


if __name__ == '__main__':
    unittest.main()
