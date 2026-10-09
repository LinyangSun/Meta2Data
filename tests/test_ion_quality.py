"""Ion EE estimation must use all post-primer reads and never mask bad inputs."""
import gzip
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import ion_quality
import parameters


class IonQualityMathTests(unittest.TestCase):
    def test_known_q19_conversions_without_rounding_the_processing_threshold(self):
        # Independent numeric anchors, not a second copy of the helper formula.
        self.assertAlmostEqual(
            ion_quality.maxee_from_median_length(400), 5.03570164717667, places=12)
        self.assertAlmostEqual(
            ion_quality.maxee_from_median_length(407), 5.123826426002262, places=12)
        self.assertNotEqual(ion_quality.maxee_from_median_length(400), 5.04)

    def test_q20_and_fractional_median_lengths(self):
        self.assertAlmostEqual(
            ion_quality.maxee_from_median_length(400, qscore=20), 4.0)
        self.assertAlmostEqual(
            ion_quality.maxee_from_median_length(150.5, qscore=20), 1.505)
        self.assertEqual(ion_quality.maxee_from_median_length(400, qscore=0), 400.0)

    def test_invalid_length_cannot_generate_an_automatic_threshold(self):
        for value in (0, -1, True, False, float('nan'), float('inf'),
                      -float('inf'), None, '400'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                ion_quality.maxee_from_median_length(value)

    def test_invalid_target_quality_is_rejected(self):
        for value in (-1, True, False, float('nan'), float('inf'), None, '19'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                ion_quality.maxee_from_median_length(400, qscore=value)


class IonQualityInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.input = self.root / 'post_primer'
        self.input.mkdir()

    def write_fastq(self, name, lengths, quality='I', newline='\n', final_newline=True):
        path = self.input / name
        parts = []
        for index, length in enumerate(lengths):
            parts.extend(['@read_{}'.format(index), 'A' * length, '+', quality * length])
        content = newline.join(parts)
        if final_newline and parts:
            content += newline
        if path.name.lower().endswith('.gz'):
            path.write_bytes(gzip.compress(content.encode('ascii')))
        else:
            path.write_bytes(content.encode('ascii'))
        return path

    def test_median_is_weighted_by_reads_instead_of_files(self):
        # Per-file medians are 100 and 400, but the pooled median is 400, not 250.
        self.write_fastq('small.fastq', [100])
        self.write_fastq('large.fastq.gz', [400, 400, 400])
        result = ion_quality.estimate_maxee(self.input)
        self.assertEqual(result['read_count'], 4)
        self.assertEqual(result['median_length'], 400)
        self.assertAlmostEqual(result['automatic_maxee'], 5.03570164717667, places=12)
        self.assertEqual(result['maxee'], result['automatic_maxee'])
        self.assertEqual(result['qscore'], 19.0)
        self.assertEqual(result['source'], 'automatic')
        self.assertEqual(len(result['files']), 2)

    def test_odd_count_and_even_count_have_exact_medians(self):
        source = self.write_fastq('sample.fq', [100, 401, 200])
        self.assertEqual(ion_quality.estimate_maxee(self.input)['median_length'], 200)
        source.unlink()
        self.write_fastq('sample.fq', [201, 100])
        result = ion_quality.estimate_maxee(self.input, qscore=20)
        self.assertEqual(result['median_length'], 150.5)
        self.assertAlmostEqual(result['maxee'], 1.505)

    def test_all_supported_suffixes_are_read_without_recursive_discovery(self):
        for name, length in [('a.fastq', 100), ('b.fq', 200),
                             ('c.fastq.gz', 300), ('d.fq.gz', 400)]:
            self.write_fastq(name, [length])
        (self.input / 'notes.txt').write_text('Not FASTQ\n')
        nested = self.input / 'nested'
        nested.mkdir()
        (nested / 'ignored.fastq').write_text('Broken, not a direct input\n')
        result = ion_quality.estimate_maxee(self.input)
        self.assertEqual(result['read_count'], 4)
        self.assertEqual(result['median_length'], 250)
        self.assertEqual(len(result['files']), 4)

    def test_quality_strings_do_not_change_the_length_based_q19_target(self):
        source = self.write_fastq('sample.fastq', [400], quality='I')
        first = ion_quality.estimate_maxee(self.input)
        source.unlink()
        self.write_fastq('sample.fastq', [400], quality='5')
        second = ion_quality.estimate_maxee(self.input)
        self.assertEqual(first['maxee'], second['maxee'])
        self.assertEqual(second['qscore'], 19.0)

    def test_override_zero_is_a_real_override_not_automatic(self):
        self.write_fastq('sample.fastq', [400])
        for override in (0, 0.0, 2, 3.75):
            with self.subTest(override=override):
                result = ion_quality.estimate_maxee(self.input, override=override)
                self.assertEqual(result['source'], 'override')
                self.assertEqual(result['maxee'], override)
                self.assertEqual(result['median_length'], 400)
                self.assertEqual(result['read_count'], 1)
                self.assertAlmostEqual(result['automatic_maxee'], 5.03570164717667)

    def test_invalid_override_is_not_silently_replaced_by_automatic(self):
        self.write_fastq('sample.fastq', [400])
        for override in (-1, True, False, float('nan'), float('inf'), '2'):
            with self.subTest(override=override), self.assertRaises(ValueError):
                ion_quality.estimate_maxee(self.input, override=override)

    def test_crlf_and_missing_final_newline_are_valid(self):
        self.write_fastq('windows.fastq', [100], newline='\r\n')
        self.write_fastq('last.fq.gz', [200], final_newline=False)
        result = ion_quality.estimate_maxee(self.input)
        self.assertEqual(result['read_count'], 2)
        self.assertEqual(result['median_length'], 150)

    def test_empty_dataset_is_an_error_even_with_an_override(self):
        for override in (None, 2):
            with self.subTest(override=override), self.assertRaises(ValueError):
                ion_quality.estimate_maxee(self.input, override=override)

    def test_empty_file_is_not_ignored_beside_a_valid_sample(self):
        self.write_fastq('valid.fastq', [400])
        for name in ('empty.fastq', 'empty.fq.gz'):
            path = self.write_fastq(name, [])
            with self.subTest(name=name), self.assertRaises(ValueError):
                ion_quality.estimate_maxee(self.input)
            path.unlink()

    def test_empty_read_is_not_counted_as_length_zero(self):
        self.write_fastq('sample.fastq', [400, 0])
        with self.assertRaises(ValueError):
            ion_quality.estimate_maxee(self.input)

    def test_broken_fastq_after_valid_record_cannot_produce_partial_estimate(self):
        source = self.write_fastq('sample.fastq', [400])
        prefix = source.read_bytes()
        malformed = {
            'header': b'read2\nACGT\n+\nIIII\n',
            'plus': b'@read2\nACGT\nminus\nIIII\n',
            'unequal_lengths': b'@read2\nACGT\n+\nIII\n',
            'missing_quality': b'@read2\nACGT\n+\n',
            'incomplete_header': b'@read2\n',
        }
        for name, tail in malformed.items():
            source.write_bytes(prefix + tail)
            with self.subTest(name=name), self.assertRaises(ValueError):
                ion_quality.estimate_maxee(self.input, override=2)

    def test_truncated_gzip_cannot_return_a_partial_estimate(self):
        path = self.write_fastq('sample.fq.gz', [400, 401])
        path.write_bytes(path.read_bytes()[:-8])
        with self.assertRaises((ValueError, OSError, EOFError)):
            ion_quality.estimate_maxee(self.input)

    def test_repeated_resolved_file_is_not_double_counted(self):
        path = self.write_fastq('sample.fastq', [400])
        (self.input / 'alias.fastq').symlink_to(path)
        with self.assertRaises(ValueError):
            ion_quality.estimate_maxee(self.input)

    def test_inputs_are_unchanged(self):
        self.write_fastq('a.fastq', [100, 200])
        self.write_fastq('b.fq.gz', [300, 400])
        before = {path.name: path.read_bytes() for path in self.input.iterdir()}
        ion_quality.estimate_maxee(self.input)
        self.assertEqual(before, {path.name: path.read_bytes() for path in self.input.iterdir()})


class IonQualityParameterTests(unittest.TestCase):
    KEYS = ('vsearch.ion_maxee', 'dada2.ion_maxee')

    def test_both_ion_methods_default_to_automatic(self):
        values = parameters.load()
        for key in self.KEYS:
            with self.subTest(key=key):
                self.assertIsNone(values[key])

    def test_null_and_nonnegative_overrides_survive_json_loading(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'parameters.json'
            for value in (None, 0, 0.0, 2, 5.123):
                path.write_text(json.dumps({'vsearch': {'ion_maxee': value},
                                            'dada2': {'ion_maxee': value}}))
                values = parameters.load(path)
                for key in self.KEYS:
                    with self.subTest(value=value, key=key):
                        self.assertEqual(values[key], value)

    def test_negative_bool_and_nonfinite_overrides_are_rejected(self):
        for key in self.KEYS:
            for value in (-1, True, False, float('nan'), float('inf'), '2'):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    parameters.validate(key, value)

    def test_null_environment_and_zero_environment_remain_distinct(self):
        for value, expected in [('null', None), ('0', 0), ('2.5', 2.5)]:
            environment = {parameters.SPEC[key][1]: value for key in self.KEYS}
            with self.subTest(value=value), patch.dict(os.environ, environment, clear=True):
                result = parameters.load(environment=True)
                for key in self.KEYS:
                    self.assertEqual(result[key], expected)

    def test_nullable_ion_setting_does_not_relax_other_numeric_parameters(self):
        with self.assertRaises(ValueError):
            parameters.validate('vsearch.maxee', None)


if __name__ == '__main__':
    unittest.main()
