"""Project-level Illumina routing: unit, shell, ledger and real-tool regressions."""
import copy
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import illumina_pe as pe
import entropy_primer_detect as primer
import parameters
import primer_sync
import read_counts


def rec(name, sequence='ACGT' * 20, quality=None):
    return '@' + name, sequence, '+', quality if quality is not None else 'I' * len(sequence)


def put(path, rows):
    with pe.open_fastq(path, 'wt') as stream:
        for row in rows:
            pe.write_record(stream, row)


class FakeVsearch:
    """Only external commands are substituted; FASTQ, routing and audits are real."""
    def __init__(self, merges=(), fail_at=None, ee=True):
        self.merges = list(merges)
        self.fail_at = fail_at
        self.ee = ee
        self.commands = []
        self.merge_index = 0

    def __call__(self, command, log):
        command = list(map(str, command))
        self.commands.append(command)
        if '--version' in command:
            Path(log).write_text('vsearch test double\n')
            return
        if self.fail_at and self.fail_at in command:
            Path(log).write_text('simulated tool execution failure\n')
            raise subprocess.CalledProcessError(2, command)
        def arg(flag):
            return command[command.index(flag) + 1]
        if '--fastq_mergepairs' in command:
            rows = list(pe.records(arg('--fastq_mergepairs')))
            amount = self.merges[self.merge_index]
            self.merge_index += 1
            put(arg('--fastqout'), rows[:amount])
            Path(log).write_text(f'{len(rows) - amount}  overlap too short\n')
        else:
            self.assert_filter_after_merges()
            rows = list(pe.records(arg('--fastq_filter')))
            limit = float(arg('--fastq_maxee'))
            filtered = [r for r in rows if not self.ee or
                        sum(10 ** (-(ord(c) - 33) / 10) for c in r[3]) <= limit]
            put(arg('--fastqout'), filtered)
            Path(log).write_text(f'{len(filtered)} reads kept\n')

    def assert_filter_after_merges(self):
        if self.merge_index != len(self.merges):
            raise AssertionError('EE filtering started before all project merge attempts finished')


class RoutingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'input'
        self.source.mkdir()
        self.output = self.root / 'output'
        self.report = self.root / 'reports/project.json'
        env = patch.dict(os.environ, {}, clear=True)
        env.start()
        self.addCleanup(env.stop)

    def sample(self, name, n, bad_r2=0, quality=None):
        first = [rec(f'{name}-{i}/1', quality=quality) for i in range(n)]
        second = [rec(f'{name}-{i}/2', quality=('!' * 80 if i < bad_r2 else None))
                  for i in range(n)]
        put(self.source / f'{name}_R1.fastq', first)
        put(self.source / f'{name}_R2.fastq', second)

    def run_process(self, merges, **kwargs):
        runner = FakeVsearch(merges)
        rc = pe.process(self.source, self.output, self.report, runner=runner, **kwargs)
        return rc, json.loads(self.report.read_text()), runner

    def test_tail_stops_at_first_q20_and_leaves_internal_low_quality(self):
        quality = ''.join(chr(q+33) for q in [30, 30, 5, 20, 5, 5])
        row = rec('read', 'ACGTAC', quality)
        result = pe.trim_tail(row)
        self.assertEqual(result[1], 'ACGT')
        self.assertEqual(result[3], quality[:4])
        self.assertEqual(pe.trim_tail(rec('all_low', 'ACG', '!!!'))[1], '')
        # The requested rule tests Q only, not N or polyG content.
        self.assertEqual(pe.trim_tail(rec('highq_n', 'ACGN'))[1], 'ACGN')

    def test_exact_half_selects_merged(self):
        self.assertEqual(pe.decide_route(5, 10, 12),
                         ('merged', 0.5, 'merge_fraction_at_or_above_threshold'))
        self.assertEqual(pe.decide_route(4, 10, 12)[0], 'R1')

    def test_weighted_fraction_not_mean_of_sample_fractions(self):
        self.sample('A', 10)
        self.sample('B', 2)
        rc, report, tool = self.run_process([6, 0])
        self.assertEqual(rc, 0)
        self.assertEqual(report['merge_fraction'], 0.5)
        self.assertEqual(report['route'], 'merged')
        self.assertEqual(report['samples']['B']['quality_filtered_reads'], 0)
        self.assertEqual(report['samples']['B']['exclusion_reason'], 'no_reads_in_project_route')
        self.assertFalse((self.output / 'B.fastq').exists())
        self.assertEqual(len(list(pe.records(self.output / 'A.fastq'))), 6)
        self.assertTrue(all('--fastq_maxee' not in c for c in tool.commands if '--fastq_mergepairs' in c))

    def test_fallback_includes_successful_pairs_and_r2_orphans(self):
        self.sample('A', 10, bad_r2=2)
        rc, report, _ = self.run_process([3])
        self.assertEqual(rc, 0)
        self.assertEqual(report['merge_input_pairs'], 8)
        self.assertEqual(report['route'], 'R1')
        self.assertEqual(report['samples']['A']['orphan_r1'], 2)
        self.assertEqual(len(list(pe.records(self.output / 'A.fastq'))), 10)
        self.assertEqual(report['samples']['A']['selected_input_reads'], 10)

    def test_no_pairs_but_good_r1_falls_back_with_na_fraction(self):
        self.sample('A', 4, bad_r2=4)
        rc, report, tool = self.run_process([])
        self.assertEqual(rc, 0)
        self.assertIsNone(report['merge_fraction'])
        self.assertEqual(report['route_reason'], 'no_usable_pairs')
        self.assertEqual(report['route'], 'R1')
        self.assertEqual(report['quality_filtered_reads'], 4)
        self.assertFalse(any('--fastq_mergepairs' in c for c in tool.commands))

    def test_valid_empty_fastq_is_not_a_software_failure(self):
        self.sample('A', 0)
        rc, report, _ = self.run_process([])
        self.assertEqual(rc, 98)
        self.assertEqual(report['status'], 'skipped')
        self.assertEqual(report['reason'], 'no_usable_reads')

    def test_filter_uses_1point5_after_route_and_never_truncates(self):
        # Interior Q0 contributes one expected error; good terminal base stops trimming.
        q = '!' + 'I' * 79
        self.sample('A', 2, quality=q)
        rc, report, tool = self.run_process([0])
        self.assertEqual(rc, 0)
        self.assertEqual(report['parameters']['maxee'], 1.5)
        self.assertTrue(any('--fastq_maxee' in c and c[c.index('--fastq_maxee')+1] == '1.5'
                            for c in tool.commands))
        self.assertEqual([len(r[1]) for r in pe.records(self.output / 'A.fastq')], [80, 80])
        self.assertFalse(any('--fastq_truncee' in c for c in tool.commands))

    def test_merged_route_does_not_switch_back_after_ee_loss(self):
        self.sample('A', 2, quality='!!' + 'I' * 78)
        rc, report, _ = self.run_process([2])
        self.assertEqual(rc, 98)
        self.assertEqual(report['route'], 'merged')
        self.assertEqual(report['quality_filtered_reads'], 0)
        self.assertEqual(report['reason'], 'no_reads_after_ee_filter')

    def test_tool_execution_failure_never_becomes_r1_fallback(self):
        self.sample('A', 2)
        with self.assertRaises(subprocess.CalledProcessError):
            pe.process(self.source, self.output, self.report,
                       runner=FakeVsearch([0], fail_at='--fastq_mergepairs'))
        report = json.loads(self.report.read_text())
        self.assertEqual(report['status'], 'failed')
        self.assertIsNone(report['samples']['A']['merged_pairs'])
        self.assertEqual(report['samples']['A']['merge_status'], 'failed')
        self.assertIsNone(report['route'])
        self.assertEqual(list(self.output.iterdir()), [])

    def test_corrupt_or_misaligned_fastq_fails(self):
        self.sample('A', 2)
        put(self.source / 'A_R2.fastq', [rec('different/2'), rec('A-1/2')])
        with self.assertRaisesRegex(ValueError, 'ID/order mismatch'):
            self.run_process([0])
        self.assertEqual(json.loads(self.report.read_text())['status'], 'failed')

    def test_missing_r2_fails_instead_of_zero_merge(self):
        put(self.source / 'A_R1.fastq', [rec('A/1')])
        with self.assertRaises(ValueError):
            self.run_process([])
        self.assertEqual(json.loads(self.report.read_text())['status'], 'failed')

    def test_truncated_fastq_and_gzip_raise(self):
        path = self.source / 'bad.fastq'
        path.write_text('@read\nACGT\n+\n')
        with self.assertRaisesRegex(ValueError, 'Truncated'):
            list(pe.records(path))
        zipped = self.source / 'bad.fastq.gz'
        zipped.write_bytes(b'not-gzip')
        with self.assertRaises(OSError):
            list(pe.records(zipped))

    def test_existing_outputs_cannot_be_silently_reused(self):
        self.sample('A', 2)
        self.output.mkdir()
        old = self.output / 'old.fastq'
        old.write_text('old data')
        with self.assertRaisesRegex(ValueError, 'refusing to reuse stale'):
            self.run_process([2])
        self.assertEqual(old.read_text(), 'old data')

    def test_merge_reasons_and_tail_metrics_are_saved(self):
        self.sample('A', 3, bad_r2=1)
        rc, report, _ = self.run_process([1])
        self.assertEqual(rc, 0)
        row = report['samples']['A']
        self.assertEqual(row['merge_failure_reasons'], {'overlap too short': 1})
        self.assertEqual(row['tail_trimming']['r2']['removed_bases'], 80)
        self.assertTrue(Path(row['merge_log']).exists())

    def test_read_count_ledger_and_final_table_are_consistent(self):
        self.sample('A', 10, bad_r2=2)
        dataset = self.root / 'PRJ_TEST'
        dataset.mkdir()
        (dataset / 'PRJ_TEST_sra.txt').write_text('RUN_A\tA\n')
        state_path = read_counts.begin(dataset, self.root, 'vsearch')
        state = json.loads(Path(state_path).read_text())
        read_counts.record(state, 'RawReads', {'A': 10}, 'raw')
        read_counts.record(state, 'primer_trimmed_reads', {'A': 10}, 'primer', 'RawReads')
        with patch.dict(os.environ, {'READ_COUNTS_STATE': state_path}):
            rc, report, _ = self.run_process([3])
        self.assertEqual(rc, 0)
        state = json.loads(Path(state_path).read_text())
        self.assertEqual(state['stages']['vsearch_merge_input_pairs']['total'], 8)
        self.assertEqual(state['stages']['vsearch_preprocess_input_reads']['sample_sources'], {'A': 'R1'})
        table = self.root / 'final.qza'
        data = dict(columns=[{'id': 'A'}], matrix_type='sparse', data=[[0, 0, 7]])
        with zipfile.ZipFile(table, 'w') as archive:
            archive.writestr('uuid/data/feature-table.biom', json.dumps(data))
        pe.finalize(self.report, table)
        report = json.loads(self.report.read_text())
        self.assertEqual(report['final_retained_samples'], 1)
        self.assertEqual(report['final_reads'], 7)
        self.assertEqual(report['samples']['A']['raw_pairs'], 10)
        self.assertEqual(report['samples']['A']['final_reads'], 7)


class PrimerSyncTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source, self.output = self.root / 'in', self.root / 'out'
        self.source.mkdir(); self.output.mkdir()
        put(self.source / 'A_R1.fastq', [rec('x/1')])
        put(self.source / 'A_R2.fastq', [rec('x/2')])
        self.settings = parameters.load()

    def decision(self, length=0, decision='no_primer'):
        result = primer.detect_for_reads([], 'test', [], self.settings)
        result.update(decision=decision, primer_length=length, detected=length > 0,
                      primer_name='test', reads_used=1, message=decision)
        return result

    def process(self, first, second, **kwargs):
        with patch.object(primer, '_filtered_reads', return_value=[('ACGT' * 20, 'I' * 80)]), \
             patch.object(primer, 'detect_mixed_orientation', return_value=None), \
             patch.object(primer, 'detect_for_reads', side_effect=[first, second]):
            return primer.process_dataset(str(self.source), str(self.output), self.settings,
                                           sync_r2_to_r1=True, **kwargs)

    def test_authorized_r1_cut_overrides_r2_no_cut_or_other_length_once(self):
        for reverse_length in (0, 5):
            with self.subTest(reverse_length=reverse_length):
                status = self.process(self.decision(20, 'known_database_match'), self.decision(reverse_length))
                self.assertEqual(status, 0)
                self.assertEqual(len(next(pe.records(self.output / 'A_R1.fastq'))[1]), 60)
                self.assertEqual(len(next(pe.records(self.output / 'A_R2.fastq'))[1]), 60)
                info = json.loads((self.output / 'primer_info.json').read_text())
                self.assertEqual(info['reverse_primer']['independent_trim_length'], reverse_length)
                self.assertEqual(info['reverse_primer']['trim_length'], 20)

    def test_r1_no_cut_preserves_independent_r2_cut(self):
        self.assertEqual(self.process(self.decision(), self.decision(5, 'known_database_match')), 0)
        self.assertEqual(len(next(pe.records(self.output / 'A_R1.fastq'))[1]), 80)
        self.assertEqual(len(next(pe.records(self.output / 'A_R2.fastq'))[1]), 75)

    def test_unknown_skip_happens_before_sync_and_writes_no_fastq(self):
        self.settings['primer.skip_unknown'] = True
        self.assertEqual(self.process(self.decision(20, 'unknown_primer'), self.decision()), 98)
        self.assertEqual(list(self.output.glob('*.fastq')), [])
        self.assertEqual(json.loads((self.output / 'primer_info.json').read_text())['forward_primer']['trim_length'], 0)

    def test_unknown_allowed_actual_r1_trim_is_synchronized(self):
        self.assertEqual(self.process(self.decision(20, 'unknown_primer'), self.decision()), 0)
        self.assertEqual(len(next(pe.records(self.output / 'A_R2.fastq'))[1]), 60)

    def test_detect_only_never_trims_or_synchronizes(self):
        self.assertEqual(self.process(self.decision(20, 'known_database_match'), self.decision(), detect_only=True), 0)
        self.assertEqual((self.output / 'A_R2.fastq').read_bytes(), (self.source / 'A_R2.fastq').read_bytes())

    def test_mixed_orientation_uses_each_original_r1_group(self):
        trims = {'r1': {'forward': 20, 'reverse': 0}, 'r2': {'forward': 7, 'reverse': 0}}
        info = {'layout': 'PE', 'original_ends': {end: {group: {} for group in groups}
                                                for end, groups in trims.items()}}
        primer_sync.synchronize_plan(info, 20, 0, trims)
        self.assertEqual(trims, {'r1': {'forward': 20, 'reverse': 0},
                                 'r2': {'forward': 7, 'reverse': 20}})

    def test_other_branches_keep_their_existing_default(self):
        config = parameters.load()
        self.assertEqual(config['vsearch.illumina_pe_maxee'], 1.5)
        self.assertEqual(config['vsearch.maxee'], 1.0)
        self.assertEqual(config['vsearch.ion_maxee'], None)


class ShellScopeTests(unittest.TestCase):
    def test_primer_sync_flag_only_reaches_the_target_branch(self):
        source = (ROOT / 'scripts/run.sh').read_text()
        function = source[source.index('_trim_primers() {'):source.index('_local_register_dataset() {')]
        cases = [('vsearch', 'ILLUMINA', 'paired', 'normal', True),
                 ('dada2', 'ILLUMINA', 'paired', 'normal', False),
                 ('vsearch', 'ILLUMINA', 'single', 'normal', False),
                 ('vsearch', 'ILLUMINA', 'paired', 'degraded_binned', False),
                 ('vsearch', 'PACBIO_SMRT', 'single', 'normal', False)]
        for method, platform, layout, quality, expected in cases:
            for explicit in ('', 'ACGT'):
                with self.subTest(method=method, platform=platform, layout=layout,
                                  quality=quality, explicit=bool(explicit)), tempfile.TemporaryDirectory() as temp:
                    root = Path(temp)
                    script = function + r"""
python3() { printf '%s\n' "$@" >> "$CALLS"; }
Audit_Fastq() { :; }
_log_status() { :; }
_trim_primers "$dataset_path/input" "$dataset_path/output"
"""
                    env = dict(os.environ, MODE=method, platform=platform, sequence_type=layout,
                               quality_status=quality, PRIMER_FWD=explicit, PRIMER_REV='',
                               dataset_path=str(root), dataset_ID='test', SCRIPTS=str(ROOT/'scripts'),
                               READ_COUNTS_REPORTS=str(root), CALLS=str(root/'calls'))
                    result = subprocess.run(['bash', '-c', script], env=env, capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual('--sync-r2-to-r1' in (root/'calls').read_text().splitlines(), expected)

    def test_preprocessor_routes_normal_pe_to_new_helper(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            script = r"""
source "$SCRIPTS/AmpliconFunction.sh"
python3() { printf '%s\n' "$@" >> "$CALLS"; }
Amplicon_Illumina_Vsearch_Preprocess
[[ "$sequence_type" == single ]]
[[ -n "$ILLUMINA_PE_REPORT" ]]
"""
            env = dict(os.environ, SCRIPTS=str(ROOT/'scripts'), CALLS=str(root/'calls'),
                       dataset_path=str(root), fastp_path=str(root/'input'),
                       quality_status='normal', sequence_type='paired', VSEARCH_MAXEE='1.0')
            env.pop('ILLUMINA_PE_MAXEE', None)
            result = subprocess.run(['bash', '-c', script], env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            args = (root/'calls').read_text().splitlines()
            self.assertTrue(args[0].endswith('illumina_pe.py'))
            self.assertEqual(args[args.index('--maxee') + 1], '1.5')

    def test_existing_degraded_branch_uses_existing_helper(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            script = r"""
source "$SCRIPTS/AmpliconFunction.sh"
python3() { printf '%s\n' "$@" >> "$CALLS"; }
Amplicon_Illumina_Vsearch_Preprocess
"""
            env = dict(os.environ, SCRIPTS=str(ROOT/'scripts'), CALLS=str(root/'calls'),
                       dataset_path=str(root), fastp_path=str(root/'input'),
                       quality_status='degraded_binned', sequence_type='paired')
            env.pop('READ_COUNTS_STATE', None)
            result = subprocess.run(['bash', '-c', script], env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            args = (root/'calls').read_text().splitlines()
            self.assertIn('degraded_quality_preprocess', args)
            self.assertFalse(any('illumina_pe.py' in a for a in args))


@unittest.skipUnless(shutil.which('vsearch'), 'Real VSEARCH is required')
class RealVsearchTests(unittest.TestCase):
    def test_real_merge_and_report(self):
        rng = random.Random(804)
        insert = ''.join(rng.choices('ACGT', k=170))
        first, second = insert[:120], insert[50:].translate(str.maketrans('ACGT','TGCA'))[::-1]
        with tempfile.TemporaryDirectory() as temp, patch.dict(os.environ, {}, clear=True):
            root = Path(temp); source = root / 'in'; source.mkdir()
            put(source / 'A_R1.fastq', [rec('r/1', first + 'ACT', 'I'*120 + '!!!')])
            put(source / 'A_R2.fastq', [rec('r/2', second)])
            report = root / 'report/project.json'
            self.assertEqual(pe.process(source, root / 'out', report), 0)
            self.assertEqual(json.loads(report.read_text())['route'], 'merged')
            self.assertEqual(next(pe.records(root / 'out/A.fastq'))[1], insert)

    def test_real_filter_accepts_ee_between_one_and_onepointfive_without_truncation(self):
        rng = random.Random(77)
        insert = ''.join(rng.choices('ACGT', k=80))
        with tempfile.TemporaryDirectory() as temp, patch.dict(os.environ, {}, clear=True):
            root = Path(temp); source = root / 'in'; source.mkdir()
            put(source / 'A_R1.fastq', [rec('r/1', insert, '!' + 'I'*79)])
            put(source / 'A_R2.fastq', [rec('r/2', insert, '!'*80)])
            report = root / 'report/project.json'
            self.assertEqual(pe.process(source, root / 'out', report), 0)
            self.assertEqual(json.loads(report.read_text())['route'], 'R1')
            self.assertEqual(next(pe.records(root / 'out/A.fastq'))[1], insert)


@unittest.skipUnless(shutil.which('cutadapt'), 'Real cutadapt is required')
class RealExplicitSyncTests(unittest.TestCase):
    def test_explicit_matches_decide_once_and_preserve_uncut_r1_behavior(self):
        forward, reverse = 'AGAGTTTGATCCTGGCTCAG', 'TACGGCTACCTTGTTACGACTT'
        body = 'TTGCAACCGGATGATCCGTA' * 5
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); source = root / 'in'; source.mkdir()
            put(source / 'A_R1.fastq', [rec('yes/1', forward+body), rec('no/1', body)])
            put(source / 'A_R2.fastq', [rec('yes/2', reverse+body), rec('no/2', reverse+body)])
            result = subprocess.run([sys.executable, str(ROOT/'scripts/explicit_primers.py'),
                                     '--input', str(source), '--output', str(root/'out'),
                                     '--forward', forward, '--reverse', reverse, '--sync-r2-to-r1'],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            first = list(pe.records(root/'out/A_R1.fastq'))
            second = list(pe.records(root/'out/A_R2.fastq'))
            self.assertEqual([r[1] for r in first], [body, body])
            self.assertEqual([r[1] for r in second], [(reverse+body)[len(forward):], body])
            info = json.loads((root/'out/primer_info.json').read_text())
            self.assertEqual(info['paired_sync']['A_R1.fastq']['synchronized_pairs'], 1)


if __name__ == '__main__':
    unittest.main()
