"""454 membership conservation, chimera exclusion and actual vsearch UC contract."""
from argparse import Namespace
import csv
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import ls454_members as members


def write_fasta(path, records):
    path.write_text(''.join(f'>{name};size={count};\n{sequence}\n'
                            for name, (sequence, count) in records.items()))


def write_fastq(path, sequences):
    path.write_text(''.join(f'@read{i};size=99; original\n{seq}\n+\n{"I" * len(seq)}\n'
                            for i, seq in enumerate(sequences)))


def uc_row(kind, query, target='*', strand='+'):
    return '\t'.join([kind, '0', '120', '99.0', strand, '0', '0', '*', query, target]) + '\n'


class MemberTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'qc'
        self.source.mkdir()
        names = ('dereplicated', 'preclustered', 'precluster_uc', 'denoised',
                 'denoise_uc', 'nonchimeras', 'chimeras', 'borderline', 'report',
                 'members_tsv', 'excluded_fasta')
        self.args = Namespace(input=self.source, output_dir=self.root / 'mapping', minsize=2,
                              **{name: self.root / name for name in names})
        self.seq = {'a': 'ACGT' * 30, 'a1': 'ACGT' * 29 + 'ACGA',
                    'b': 'GATT' * 30, 'b1': 'GATT' * 29 + 'GATC',
                    'c': 'TGCA' * 30, 'd': 'CTAG' * 30, 'e': 'CGAT' * 30}
        self.derep = {name: (seq, count) for (name, seq), count in
                      zip(self.seq.items(), [4, 1, 4, 1, 2, 1, 2])}
        self.pre = {name: (self.seq[name], count) for name, count in
                    [('a', 5), ('b', 5), ('c', 2), ('d', 1), ('e', 2)]}
        self.denoise = {name: (self.seq[name], count) for name, count in
                        [('a', 5), ('b', 7), ('e', 2)]}
        write_fasta(self.args.dereplicated, self.derep)
        write_fasta(self.args.preclustered, self.pre)
        write_fasta(self.args.denoised, self.denoise)
        self.args.precluster_uc.write_text(''.join(uc_row('S', name) for name in self.pre)
                                           + uc_row('H', 'a1', 'a') + uc_row('H', 'b1', 'b'))
        self.args.denoise_uc.write_text(''.join(uc_row('S', name) for name in self.denoise)
                                        + uc_row('H', 'c', 'b'))
        write_fasta(self.args.nonchimeras, {'a': self.denoise['a']})
        write_fasta(self.args.chimeras, {'b': self.denoise['b']})
        write_fasta(self.args.borderline, {'e': self.denoise['e']})
        write_fastq(self.source / 'good.fastq', [self.seq[name] for name in ['a'] * 4 + ['a1', 'd', 'e', 'e']])
        write_fastq(self.source / 'chimeric.fastq', [self.seq[name] for name in ['b'] * 4 + ['b1', 'c', 'c']])
        write_fastq(self.source / 'empty.fastq', [])

    def test_members_trace_both_clusters_and_conserve_each_sample(self):
        result = members.filter_members(self.args)
        self.assertEqual(result['status'], 'completed')
        self.assertEqual(result['totals'], dict(input=15, excluded_chimera=7, kept=8))
        self.assertEqual(result['samples']['chimeric'], dict(input=7, excluded_chimera=7, kept=0))
        self.assertEqual(result['samples']['empty'], dict(input=0, excluded_chimera=0, kept=0))
        self.assertEqual((self.args.output_dir / 'chimeric.fastq').read_text(), '')
        self.assertEqual((self.args.output_dir / 'empty.fastq').read_text(), '')
        rows = list(csv.DictReader(self.args.members_tsv.open(), delimiter='\t'))
        by_id = {row['derep_id']: row for row in rows}
        self.assertEqual(by_id['a1']['classification'], 'nonchimera')
        self.assertEqual(by_id['b1']['denoised_id'], 'b')
        self.assertEqual(by_id['c']['classification'], 'chimera')
        self.assertEqual(by_id['d']['classification'], 'below_unoise_minsize')
        self.assertEqual(by_id['e']['classification'], 'borderline')
        self.assertEqual(result['members']['below_unoise_minsize_reads'], 1)
        self.assertEqual(set(members.fasta(self.args.excluded_fasta)), {'b', 'b1', 'c'})

    def test_original_size_labels_cannot_inflate_dereplication(self):
        output = self.root / 'combined.fasta'
        self.assertEqual(members.to_fasta(self.source, output), 15)
        self.assertEqual(output.read_text().count('>LS454_read_'), 15)
        self.assertNotIn('size=', output.read_text())

    def test_duplicate_adjacent_size_labels_are_rejected(self):
        self.args.dereplicated.write_text('>a;size=2;size=3;\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'Invalid, repeated'):
            members.fasta(self.args.dereplicated)

    def test_invalid_uc_or_abundance_fails_without_mapping_output(self):
        original = self.args.precluster_uc.read_text()
        cases = [original.replace(uc_row('H', 'a1', 'a'), ''),
                 original + uc_row('S', 'a'),
                 original + uc_row('S', 'unknown'),
                 original.replace(uc_row('H', 'a1', 'a'), uc_row('H', 'a1', 'missing')),
                 original.replace(uc_row('H', 'a1', 'a'), uc_row('H', 'a1', 'a', '-'))]
        for text in cases:
            with self.subTest(uc=text):
                self.args.precluster_uc.write_text(text)
                with self.assertRaises(ValueError):
                    members.filter_members(self.args)
                self.assertFalse(self.args.output_dir.exists())
                self.assertEqual(json.loads(self.args.report.read_text())['status'], 'failed')
        self.args.precluster_uc.write_text(original)
        write_fasta(self.args.preclustered, dict(self.pre, a=(self.seq['a'], 99)))
        with self.assertRaisesRegex(ValueError, 'does not conserve'):
            members.filter_members(self.args)

    def test_unaccounted_high_abundance_unoise_input_fails(self):
        self.args.denoise_uc.write_text(self.args.denoise_uc.read_text().replace(uc_row('H', 'c', 'b'), ''))
        with self.assertRaisesRegex(ValueError, 'missing eligible sequence c'):
            members.filter_members(self.args)

    def test_qc_sequences_and_abundances_must_match_dereplication(self):
        good = self.source / 'good.fastq'
        original = good.read_text()
        for content, error in [(original + '@extra\nACGT\n+\nIIII\n', 'absent from dereplication'),
                               (original + original, 'do not conserve')]:
            with self.subTest(error=error):
                good.write_text(content)
                with self.assertRaisesRegex(ValueError, error):
                    members.filter_members(self.args)
                self.assertFalse(self.args.output_dir.exists())

    def test_all_chimera_reads_fail_without_publishing_mapping(self):
        # Remove the below-minsize read so all remaining members can be classified.
        write_fasta(self.args.dereplicated, {k: v for k, v in self.derep.items() if k != 'd'})
        write_fasta(self.args.preclustered, {k: v for k, v in self.pre.items() if k != 'd'})
        self.args.precluster_uc.write_text(self.args.precluster_uc.read_text().replace(uc_row('S', 'd'), ''))
        write_fastq(self.source / 'good.fastq', [self.seq[name] for name in ['a'] * 4 + ['a1', 'e', 'e']])
        write_fasta(self.args.nonchimeras, {})
        write_fasta(self.args.borderline, {})
        write_fasta(self.args.chimeras, self.denoise)
        with self.assertRaisesRegex(ValueError, 'All QC reads'):
            members.filter_members(self.args)
        self.assertFalse(self.args.output_dir.exists())
        report = json.loads(self.args.report.read_text())
        self.assertEqual(report['totals'], dict(input=14, excluded_chimera=14, kept=0))
        self.assertEqual(report['status'], 'failed')

    def test_unsafe_destinations_never_overwrite_inputs_or_evidence(self):
        protected = [self.source / 'good.fastq', self.args.dereplicated, self.args.precluster_uc]
        before = {path: path.read_bytes() for path in protected}
        for flag in ['report', 'members_tsv', 'excluded_fasta']:
            original = getattr(self.args, flag)
            for destination in protected + [self.args.output_dir / 'evidence.txt']:
                with self.subTest(flag=flag, destination=destination):
                    setattr(self.args, flag, destination)
                    with self.assertRaisesRegex(ValueError, 'Output'):
                        members.filter_members(self.args)
                    for path, data in before.items():
                        self.assertEqual(path.read_bytes(), data)
                    self.assertFalse(self.args.output_dir.exists())
            setattr(self.args, flag, original)
        self.assertFalse(self.args.report.exists())

    @unittest.skipUnless(shutil.which('vsearch'), 'vsearch is required for the real UC contract')
    def test_actual_vsearch_uc_sizes_singletons_one_n_and_empty_sample(self):
        rng = random.Random(454)
        a = ''.join(rng.choices('ACGT', k=200))
        near = a[:75] + ('A' if a[75] != 'A' else 'C') + a[76:]
        one_n = a[:120] + 'N' + a[121:]
        b = ''.join(rng.choices('ACGT', k=200))
        low = ''.join(rng.choices('ACGT', k=200))
        for path in self.source.iterdir():
            path.unlink()
        write_fastq(self.source / 'valid.fastq', [a, near, one_n] + [b] * 12 + [low])
        write_fastq(self.source / 'empty.fastq', [])
        combined = self.root / 'combined.fasta'
        members.to_fasta(self.source, combined)
        commands = [
            ['--derep_fulllength', combined, '--output', self.args.dereplicated,
             '--sizeout', '--minuniquesize', '1', '--relabel', 'LS454_'],
            ['--cluster_size', self.args.dereplicated, '--id', '0.99', '--strand', 'plus',
             '--sizein', '--sizeout', '--centroids', self.args.preclustered, '--uc', self.args.precluster_uc],
            ['--cluster_unoise', self.args.preclustered, '--sizein', '--sizeout', '--minsize', '2',
             '--centroids', self.args.denoised, '--uc', self.args.denoise_uc],
            ['--uchime3_denovo', self.args.denoised, '--sizein', '--sizeout',
             '--nonchimeras', self.args.nonchimeras, '--chimeras', self.args.chimeras,
             '--borderline', self.args.borderline]]
        for command in commands:
            length_args = [] if command[0] == '--uchime3_denovo' else ['--minseqlength', '1']
            result = subprocess.run(['vsearch'] + list(map(str, command))
                                    + length_args + ['--threads', '1'],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        # Exact dereplication must keep an N variant separate from its compatible
        # A/C/G/T template; 99% clustering may then join their tracked members.
        derep = members.fasta(self.args.dereplicated)
        self.assertEqual({v['sequence']: v['abundance'] for v in derep.values()},
                         {a: 1, near: 1, one_n: 1, b: 12, low: 1})
        result = members.filter_members(self.args)
        self.assertEqual(result['totals'], dict(input=16, excluded_chimera=0, kept=16))
        self.assertEqual(result['samples']['empty'], dict(input=0, excluded_chimera=0, kept=0))
        self.assertEqual(result['members']['below_unoise_minsize_reads'], 1)
        self.assertEqual(sorted(v['abundance'] for v in members.fasta(self.args.preclustered).values()), [1, 3, 12])
        self.assertEqual(sorted(v['abundance'] for v in members.fasta(self.args.denoised).values()), [3, 12])
        kept = list(members.fastq_records(self.args.output_dir / 'valid.fastq'))
        self.assertEqual(sum(sequence == one_n for _, sequence, _, _ in kept), 1)
        self.assertEqual(len(kept), 16)
        with self.args.members_tsv.open() as stream:
            tracked = {row['derep_id']: row for row in csv.DictReader(stream, delimiter='\t')}
        ids = {record['sequence']: identifier for identifier, record in derep.items()}
        self.assertEqual(tracked[ids[one_n]]['precluster_id'], tracked[ids[a]]['precluster_id'])
        self.assertEqual(tracked[ids[one_n]]['classification'], 'nonchimera')
        self.assertEqual(tracked[ids[one_n]]['excluded'], 'False')


if __name__ == '__main__':
    unittest.main()
