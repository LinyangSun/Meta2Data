"""Lossless ordering and failure safety before ONT's order-sensitive tools."""
from collections import Counter
import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import ont_deterministic as ordering


def fastq_text(records):
    return ''.join(f'{header}\n{sequence}\n{plus}\n{quality}\n'
                   for sequence, quality, header, plus in records)


def simple_fasta(path):
    lines = path.read_text().splitlines()
    return list(zip(lines[::2], lines[1::2]))


class OntDeterministicTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory(prefix='ONT space ')
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)
        # Repeated original names and exact duplicate records are intentional.
        self.records = [
            ('ACGT', 'IIII', '@same description\tfield', '+same'),
            ('TTAA', 'HHHH', '@same', '+'),
            ('ACGT', 'JJJJ', '@other', '+other details'),
            ('ACGT', 'IIII', '@same description\tfield', '+same'),
            ('aCgN', '!!!!', '@case', '+'),
            ('TTCG', 'ABCD', '@forward', '+'),
            ('CGAA', 'ABCD', '@reverse-complement', '+'),
        ]

    def file(self, name, value):
        path = self.root / name
        path.write_text(value)
        return path

    def assert_no_temporary_files(self):
        self.assertEqual(list(self.root.glob('.ont-sort-*')), [])

    def test_shuffled_fastq_produces_identical_unique_names_without_losing_reads(self):
        first = self.file('first.fastq', fastq_text(self.records))
        second = self.file('second.fastq', fastq_text(self.records[::-1]))
        a, b = self.root / 'a.fa', self.root / 'b.fa'
        ordering.fastq_to_fasta(first, a, 'sample with spaces')
        ordering.fastq_to_fasta(second, b, 'sample with spaces')
        self.assertEqual(a.read_bytes(), b.read_bytes())
        records = simple_fasta(a)
        self.assertEqual(len({header for header, _ in records}), len(self.records))
        self.assertEqual(Counter(sequence for _, sequence in records),
                         Counter(record[0] for record in self.records))
        self.assert_no_temporary_files()

    def test_sample_namespaces_are_unique_without_lossy_sanitizing(self):
        source = self.file('reads.fastq', fastq_text(self.records))
        a, b = self.root / 'a.fa', self.root / 'b.fa'
        ordering.fastq_to_fasta(source, a, 'sample A')
        ordering.fastq_to_fasta(source, b, 'sample_A')
        self.assertTrue({h for h, _ in simple_fasta(a)}.isdisjoint(
            {h for h, _ in simple_fasta(b)}))

    def test_long_sample_names_fit_sam_qname_limit(self):
        source = self.file('reads.fastq', fastq_text(self.records))
        a, b = self.root / 'long-a.fa', self.root / 'long-b.fa'
        ordering.fastq_to_fasta(source, a, 'sample_' + 'x' * 240)
        ordering.fastq_to_fasta(source, b, 'sample_' + 'x' * 239 + 'y')
        first = {header[1:] for header, _ in simple_fasta(a)}
        second = {header[1:] for header, _ in simple_fasta(b)}
        self.assertTrue(first.isdisjoint(second))
        self.assertTrue(all(len(name.encode('ascii')) <= 254 for name in first | second))

    def test_sequence_hash_collision_cannot_merge_records_or_duplicate_names(self):
        source = self.file('reads.fastq', fastq_text(self.records))
        destination = self.root / 'out.fa'
        with patch.object(ordering.hashlib, 'sha256',
                          return_value=Mock(hexdigest=lambda: '0' * 64)):
            ordering.fastq_to_fasta(source, destination, 's')
        records = simple_fasta(destination)
        self.assertEqual(len({header for header, _ in records}), len(self.records))
        self.assertEqual(Counter(sequence for _, sequence in records),
                         Counter(record[0] for record in self.records))

    def test_fastq_sort_preserves_all_fields_and_is_independent_of_input_order(self):
        first = self.file('first.fastq', fastq_text(self.records))
        second = self.file('second.fastq', fastq_text(self.records[::-1]))
        with patch.dict(os.environ, {'LC_ALL': 'POSIX'}):
            ordering.sort_fastq(first, first)
            ordering.sort_fastq(second, second)
            self.assertEqual(os.environ['LC_ALL'], 'POSIX')
        self.assertEqual(first.read_bytes(), second.read_bytes())
        with first.open() as stream:
            self.assertEqual(Counter(ordering._fastq_records(stream)), Counter(self.records))
        self.assert_no_temporary_files()

    def test_wrapped_fastq_is_validated_and_written_as_complete_records(self):
        source = self.file('wrapped.fastq', '@r comment\nAC\nGT\n+r\nII\nJJ\n')
        ordering.sort_fastq(source, source)
        self.assertEqual(source.read_text(), '@r comment\nACGT\n+r\nIIJJ\n')

    def test_fasta_preserves_headers_sizes_case_orientation_and_duplicates(self):
        records = [('>b;size=13;sample=x; extra\tdata', 'TTCG'),
                   ('>a;size=1; reads=8', 'CGAA'),
                   ('>c;size=7; sample=y', 'aCgN'),
                   ('>a;size=1; reads=8', 'CGAA')]
        a = self.file('a.fa', ''.join(f'{h}\n{s}\n' for h, s in records))
        b = self.file('b.fa', ''.join(f'{h}\n{s}\n' for h, s in records[::-1]))
        ordering.sort_fasta(a, a)
        ordering.sort_fasta(b, b)
        self.assertEqual(a.read_bytes(), b.read_bytes())
        self.assertEqual(Counter(simple_fasta(a)), Counter(records))
        saved = a.read_bytes()
        ordering.sort_fasta(a, a)
        self.assertEqual(a.read_bytes(), saved)
        self.assert_no_temporary_files()

    def test_gzip_input_and_paths_with_spaces(self):
        source = self.root / 'reads with spaces.fastq.gz'
        with gzip.open(source, 'wt') as stream:
            stream.write(fastq_text(self.records))
        destination = self.root / 'result with spaces.fa'
        ordering.fastq_to_fasta(source, destination, 'sample with spaces')
        self.assertEqual(len(simple_fasta(destination)), len(self.records))

    def test_corrupt_gzip_preserves_existing_output(self):
        source = self.root / 'corrupt.fastq.gz'
        source.write_bytes(gzip.compress(fastq_text(self.records).encode())[:-6])
        destination = self.file('out.fa', 'previous output\n')
        with self.assertRaises((EOFError, OSError)):
            ordering.fastq_to_fasta(source, destination, 'sample')
        self.assertEqual(destination.read_text(), 'previous output\n')
        self.assert_no_temporary_files()

    def test_invalid_fastq_fails_atomically(self):
        cases = [
            'r\nACGT\n+\nIIII\n',
            '@\nACGT\n+\nIIII\n',
            '@r\nACGT\n',
            '@r\nACGT\n+\nIII\n',
            '@r\nACGT\n+\nIIIII\n',
            '@r\nACGT\n+other\nIIII\n',
            '@r\nAC?T\n+\nIIII\n',
            '@r\n\n+\n\n',
            '@r\nACGT\n+\nI II\n',
            '@r\nACGT\n+\nIIII\ntrailing garbage\n',
        ]
        for malformed in cases:
            with self.subTest(malformed=malformed):
                source = self.file('bad.fastq', malformed)
                destination = self.file('out.fa', 'previous output\n')
                with self.assertRaises(ValueError):
                    ordering.fastq_to_fasta(source, destination, 'sample')
                self.assertEqual(destination.read_text(), 'previous output\n')
                with self.assertRaises(ValueError):
                    ordering.sort_fastq(source, source)
                self.assertEqual(source.read_text(), malformed)
                self.assert_no_temporary_files()

    def test_invalid_fasta_fails_atomically(self):
        for malformed in ('ACGT\n', '>r\n', '>\nACGT\n', '>r\nA?GT\n'):
            with self.subTest(malformed=malformed):
                source = self.file('bad.fa', malformed)
                with self.assertRaises(ValueError):
                    ordering.sort_fasta(source, source)
                self.assertEqual(source.read_text(), malformed)
                self.assert_no_temporary_files()

    def test_external_sort_failure_preserves_output_and_cleans_temporary_files(self):
        source = self.file('reads.fastq', fastq_text(self.records))
        destination = self.file('out.fa', 'previous output\n')
        with patch.object(ordering.subprocess, 'run',
                          side_effect=subprocess.CalledProcessError(2, ['sort'])):
            with self.assertRaises(subprocess.CalledProcessError):
                ordering.fastq_to_fasta(source, destination, 'sample')
        self.assertEqual(destination.read_text(), 'previous output\n')
        self.assert_no_temporary_files()

    def test_atomic_publish_failure_preserves_output_and_cleans_temporary_files(self):
        source = self.file('reads.fastq', fastq_text(self.records))
        destination = self.file('out.fa', 'previous output\n')
        with patch.object(ordering.os, 'replace', side_effect=OSError('disk full')):
            with self.assertRaises(OSError):
                ordering.fastq_to_fasta(source, destination, 'sample')
        self.assertEqual(destination.read_text(), 'previous output\n')
        self.assert_no_temporary_files()

    def test_conversion_cannot_overwrite_fastq_through_same_path_or_link(self):
        original = fastq_text(self.records)
        source = self.file('reads.fastq', original)
        linked = self.root / 'linked.fa'
        linked.hardlink_to(source)
        for destination in (source, linked):
            with self.assertRaises(ValueError):
                ordering.fastq_to_fasta(source, destination, 'sample')
        self.assertEqual(source.read_text(), original)

    def test_compressed_output_is_rejected_instead_of_writing_plaintext_as_gzip(self):
        source = self.file('reads.fastq', fastq_text(self.records))
        with self.assertRaises(ValueError):
            ordering.sort_fastq(source, self.root / 'out.fastq.gz')
        self.assertFalse((self.root / 'out.fastq.gz').exists())

    def test_empty_inputs_are_valid_and_produce_empty_outputs(self):
        source = self.file('empty', '')
        for function, extra in ((ordering.sort_fastq, ()),
                                (ordering.sort_fasta, ()),
                                (ordering.fastq_to_fasta, ('sample',))):
            destination = self.file('out', 'previous output\n')
            function(source, destination, *extra)
            self.assertEqual(destination.read_bytes(), b'')
        self.assert_no_temporary_files()

    def test_cli_accepts_paths_with_spaces_and_returns_failure_for_malformed_input(self):
        source = self.file('reads with spaces.fastq', fastq_text(self.records))
        destination = self.root / 'out with spaces.fa'
        command = [sys.executable, str(ROOT / 'scripts' / 'ont_deterministic.py'),
                   'to-fasta', '--input', str(source), '--output', str(destination),
                   '--sample', 'sample with spaces']
        success = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(success.returncode, 0, success.stderr)
        saved = destination.read_bytes()
        source.write_text('@truncated\nACGT\n+\n')
        failure = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(failure.returncode, 0)
        self.assertIn('truncated', failure.stderr)
        self.assertEqual(destination.read_bytes(), saved)
        self.assert_no_temporary_files()


if __name__ == '__main__':
    unittest.main()
