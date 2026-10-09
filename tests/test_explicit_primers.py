"""Real cutadapt checks for explicit single-end and paired-end primer handling."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / 'scripts' / 'explicit_primers.py'
FORWARD = 'AGAGTTTGATCCTGGCTCAG'
REVERSE = 'TACGGYTACCTTGTTAYGACTT'
REVERSE_BASES = REVERSE.replace('Y', 'C')
INSERT = 'TTGCAACCGGATGATCCGTA' * 5


def rc(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def quality(length, offset=0):
    return ''.join(chr(40 + (i + offset) % 35) for i in range(length))


@unittest.skipUnless(shutil.which('cutadapt'), 'Real cutadapt is required')
class ExplicitPrimerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.input = self.root / 'input'
        self.output = self.root / 'output'
        self.input.mkdir()
        self.originals = {}

    def write_reads(self, filename, reads):
        path = self.input / filename
        path.write_text(''.join(f'@{name}\n{seq}\n+\n{qual}\n'
                                for name, seq, qual in reads))
        self.originals[path] = path.read_bytes()

    def run_helper(self, *options, reverse=True):
        command = [sys.executable, str(HELPER), '--input', str(self.input),
                   '--output', str(self.output), '--forward', FORWARD]
        if reverse:
            command += ['--reverse', REVERSE]
        result = subprocess.run(command + list(options), text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        for source, original in self.originals.items():
            self.assertEqual(source.read_bytes(), original, 'Input FASTQ was modified')
        return json.loads((self.output / 'primer_info.json').read_text())

    def read_reads(self, filename):
        lines = (self.output / filename).read_text().splitlines()
        self.assertEqual(len(lines) % 4, 0)
        reads = []
        for offset in range(0, len(lines), 4):
            name, sequence, separator, qualities = lines[offset:offset + 4]
            self.assertTrue(name.startswith('@'))
            self.assertEqual(separator, '+')
            self.assertEqual(len(sequence), len(qualities))
            reads.append((name[1:], sequence, qualities))
        return reads

    def test_single_end_removes_both_primers(self):
        sequence = FORWARD + INSERT + rc(REVERSE_BASES)
        qualities = quality(len(sequence))
        self.write_reads('sample.fastq', [('forward original description', sequence, qualities)])
        info = self.run_helper()
        self.assertEqual(self.read_reads('sample.fastq'), [
            ('forward original description', INSERT,
             qualities[len(FORWARD):-len(REVERSE_BASES)])])
        self.assertIsNone(info['reverse_primer']['trim_length'])
        self.assertEqual(info['reverse_primer']['trim_method'], 'cutadapt_variable')
        self.assertFalse(info['reverse_complement_search'])

    def test_short_single_end_without_reverse_site_keeps_insert(self):
        sequence = FORWARD + INSERT
        qualities = quality(len(sequence))
        self.write_reads('sample.fastq', [('short', sequence, qualities)])
        self.run_helper()
        self.assertEqual(self.read_reads('sample.fastq'), [
            ('short', INSERT, qualities[len(FORWARD):])])

    def test_missing_forward_site_preserves_read_without_filtering(self):
        sequence = INSERT + rc(REVERSE_BASES)
        reads = [('no_forward', sequence, quality(len(sequence))),
                 ('no_primers', INSERT, quality(len(INSERT), offset=3))]
        self.write_reads('sample.fastq', reads)
        self.run_helper()
        self.assertEqual(self.read_reads('sample.fastq'), reads)

    def test_mixed_orientation_trims_both_ends_and_preserves_ids_and_quality(self):
        sequence = FORWARD + INSERT + rc(REVERSE_BASES)
        forward_quality = quality(len(sequence))
        reverse_quality = quality(len(sequence), offset=7)
        self.write_reads('sample.fastq', [
            ('forward exact header', sequence, forward_quality),
            ('reverse exact header', rc(sequence), reverse_quality)])
        info = self.run_helper('--mixed-orientation')
        self.assertEqual(self.read_reads('sample.fastq'), [
            ('forward exact header', INSERT,
             forward_quality[len(FORWARD):-len(REVERSE_BASES)]),
            ('reverse exact header', INSERT,
             reverse_quality[::-1][len(FORWARD):-len(REVERSE_BASES)])])
        self.assertTrue(info['reverse_complement_search'])

    def test_mixed_forward_only_never_infers_reverse_primer(self):
        sequence = FORWARD + INSERT + rc(REVERSE_BASES)
        forward_quality = quality(len(sequence))
        reverse_quality = quality(len(sequence), offset=11)
        self.write_reads('sample.fastq', [
            ('forward', sequence, forward_quality),
            ('reverse', rc(sequence), reverse_quality)])
        info = self.run_helper('--mixed-orientation', reverse=False)
        self.assertEqual(self.read_reads('sample.fastq'), [
            ('forward', INSERT + rc(REVERSE_BASES), forward_quality[len(FORWARD):]),
            ('reverse', INSERT + rc(REVERSE_BASES), reverse_quality[::-1][len(FORWARD):])])
        self.assertFalse(info['reverse_primer']['detected'])
        self.assertEqual(info['reverse_primer']['consensus'], '')
        self.assertEqual(info['reverse_primer']['trim_length'], 0)

    def test_paired_end_keeps_existing_per_mate_front_trimming(self):
        first = FORWARD + INSERT + rc(REVERSE_BASES)
        second = REVERSE_BASES + rc(INSERT) + rc(FORWARD)
        first_quality = quality(len(first))
        second_quality = quality(len(second), offset=5)
        self.write_reads('sample_R1.fastq', [('pair/1 note', first, first_quality)])
        self.write_reads('sample_R2.fastq', [('pair/2 note', second, second_quality)])
        info = self.run_helper('--mixed-orientation')
        self.assertEqual(self.read_reads('sample_R1.fastq'), [
            ('pair/1 note', INSERT + rc(REVERSE_BASES), first_quality[len(FORWARD):])])
        self.assertEqual(self.read_reads('sample_R2.fastq'), [
            ('pair/2 note', rc(INSERT) + rc(FORWARD), second_quality[len(REVERSE_BASES):])])
        self.assertEqual(info['layout'], 'PE')
        self.assertFalse(info['reverse_complement_search'])

    def test_detect_only_copies_reads_and_records_both_primers_without_applying(self):
        sequence = FORWARD + INSERT + rc(REVERSE_BASES)
        self.write_reads('sample.fastq', [
            ('forward description', sequence, quality(len(sequence))),
            ('reverse description', rc(sequence), quality(len(sequence), offset=9))])
        info = self.run_helper('--detect-only', '--mixed-orientation')
        self.assertEqual((self.output / 'sample.fastq').read_bytes(),
                         self.originals[self.input / 'sample.fastq'])
        for name, sequence in [('forward_primer', FORWARD), ('reverse_primer', REVERSE)]:
            self.assertEqual(info[name]['consensus'], sequence)
            self.assertEqual(info[name]['trim_length'], 0)
            self.assertEqual(info[name]['trim_method'], 'none')
        self.assertTrue(info['detect_only'])
        self.assertFalse(info['reverse_complement_search'])


if __name__ == '__main__':
    unittest.main()
