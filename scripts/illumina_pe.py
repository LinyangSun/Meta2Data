"""Q20-tail/project-route preprocessing for normal Illumina PE VSEARCH input.

Keep all usable R1 reads separately, finish every merge attempt before deciding
one project route, then filter the selected input by expected errors. No EE
truncation, extra N filter, or new minimum-length threshold is introduced.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
from itertools import zip_longest
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

TAIL_QUALITY = 20
MERGE_REASONS = (
    'reads too short', 'reads too long', 'too many N', 'too many differences',
    'too high percentage of differences', 'overlap too short',
    'expected error too high', 'merged fragment too short',
    'merged fragment too long', 'staggered read pairs', 'too few kmers',
    'multiple potential alignments', 'alignment score too low',
)


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(data, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def open_fastq(path, mode='rt'):
    return gzip.open(path, mode) if str(path).lower().endswith('.gz') else open(path, mode)


def records(path):
    """Read complete four-line FASTQ records, including valid empty records."""
    with open_fastq(path) as stream:
        index = 0
        while True:
            header = stream.readline()
            if not header:
                return
            index += 1
            lines = [stream.readline() for _ in range(3)]
            if any(line == '' for line in lines):
                raise ValueError(f'Truncated FASTQ at record {index}: {path}')
            record = tuple(line.rstrip('\r\n') for line in [header, *lines])
            head, seq, plus, qual = record
            if (not head.startswith('@') or not head[1:].split() or
                    not plus.startswith('+') or len(seq) != len(qual)):
                raise ValueError(f'Invalid FASTQ at record {index}: {path}')
            if qual and (min(qual) < '!' or max(qual) > '~'):
                raise ValueError(f'Invalid Phred+33 quality at record {index}: {path}')
            yield record


def read_id(header):
    name = header.split()[0]
    if name.startswith('@'):
        name = name[1:]
    return name[:-2] if name.endswith(('/1', '/2')) else name


def paired_records(first, second):
    for index, pair in enumerate(zip_longest(records(first), records(second)), 1):
        r1, r2 = pair
        if r1 is None or r2 is None:
            raise ValueError(f'Unequal mate counts at pair {index}: {first}, {second}')
        if read_id(r1[0]) != read_id(r2[0]):
            raise ValueError(f'Mate ID/order mismatch at pair {index}: {r1[0]}, {r2[0]}')
        yield r1, r2


def write_record(stream, record):
    stream.write('\n'.join(record) + '\n')


def trim_tail(record):
    head, seq, plus, qual = record
    end = len(qual)
    while end and ord(qual[end - 1]) - 33 < TAIL_QUALITY:
        end -= 1
    return head, seq[:end], plus, qual[:end]


def decide_route(merged, pairs, r1_count, threshold=0.5):
    if not 0 <= merged <= pairs <= r1_count:
        raise ValueError('Invalid project read-count balance')
    if not math.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError('Invalid merge threshold')
    if pairs == 0:
        return ('R1', None, 'no_usable_pairs') if r1_count else (None, None, 'no_usable_reads')
    fraction = merged / pairs
    if fraction >= threshold:
        return 'merged', fraction, 'merge_fraction_at_or_above_threshold'
    return 'R1', fraction, 'merge_fraction_below_threshold'


def run_command(command, log):
    with open(log, 'w') as stream:
        subprocess.run(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT, check=True)


def parse_merge_reasons(path):
    """Keep the full tool log even when a version uses unrecognized reason names."""
    reasons = {}
    for line in Path(path).read_text(errors='replace').splitlines():
        match = re.match(r'^\s*([0-9][0-9,]*)\s+(.+?)\s*$', line)
        if match and any(reason.lower() in match[2].lower() for reason in MERGE_REASONS):
            reasons[match[2]] = int(match[1].replace(',', ''))
    return reasons


def audit(report, fields):
    state_path = os.environ.get('READ_COUNTS_STATE')
    if not state_path:
        return
    import read_counts
    state = json.loads(Path(state_path).read_text())
    for stage, field, parent in fields:
        details = {}
        if stage == 'vsearch_preprocess_input_reads':
            details['sample_sources'] = dict.fromkeys(report['samples'], report['route'])
        read_counts.record(state, stage,
                           {sample: row[field] for sample, row in report['samples'].items()},
                           report['report_path'], parent, **details)


def save_report(report, report_path):
    atomic_json(report_path, report)
    fields = ['sample', 'raw_pairs', 'primer_trimmed_pairs', 'r1_candidates',
              'r2_candidates', 'orphan_r1', 'merge_input_pairs', 'merged_pairs',
              'merge_fraction', 'route', 'selected_input_reads',
              'quality_filtered_reads', 'final_reads', 'exclusion_reason']
    with Path(report_path).with_name('sample_counts.tsv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for sample, row in report['samples'].items():
            values = {key: row.get(key) for key in fields}
            values.update(sample=sample, route=report['route'])
            writer.writerow({key: 'NA' if value is None else value for key, value in values.items()})


def process(input_dir, output_dir, report_path, threads=1, maxee=1.5,
            merge_min=0.5, rows=None, runner=run_command):
    """Return 98 only for valid biological empty-input/output outcomes."""
    if threads < 1 or not math.isfinite(maxee) or maxee < 0:
        raise ValueError('Invalid threads or expected-error threshold')
    decide_route(0, 0, 0, merge_min)
    source, output = Path(input_dir).resolve(), Path(output_dir).resolve()
    if source == output or source in output.parents or output in source.parents:
        raise ValueError('Input and output directories must not overlap')
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError(f'Output must be empty; refusing to reuse stale files: {output}')
    report_path = Path(report_path).resolve()
    report = dict(schema_version=1, status='running', report_path=str(report_path),
                  input=str(source), output=str(output), route=None, reason=None,
                  merge_fraction=None, samples={}, commands=[],
                  parameters=dict(tail_quality=TAIL_QUALITY, tail_window=1,
                                  input_min_length=1, maxee=maxee, merge_min_fraction=merge_min),
                  final_reads=None, final_retained_samples=None)
    logs = report_path.parent / 'logs'
    logs.mkdir(parents=True, exist_ok=True)

    def execute(command, log):
        report['commands'].append(list(map(str, command)))
        save_report(report, report_path)
        runner(command, log)

    try:
        raw_counts, expected_samples = {}, None
        state_path = os.environ.get('READ_COUNTS_STATE')
        if state_path:
            state = json.loads(Path(state_path).read_text())
            raw_counts = state.get('stages', {}).get('RawReads', {}).get('counts', {})
            expected_samples = set(state['samples'])
        if rows is None:
            from read_layout import discover, layout
            rows = discover(str(source))
            if layout(rows) != 'PE':
                raise ValueError('This preprocessor accepts paired-end inputs only')
        grouped = defaultdict(list)
        for row in rows:
            sample = row['sample']
            if (not sample or sample in ('.', '..') or
                    any(char in sample for char in '/\\\t\r\n') or not row.get('r2')):
                raise ValueError(f'Invalid paired sample: {sample!r}')
            grouped[sample].append(row)
        if not grouped:
            raise ValueError('No paired FASTQ files found')
        if expected_samples is not None and set(grouped) != expected_samples:
            raise ValueError('Paired FASTQ sample set differs from the registered dataset')
        execute(['vsearch', '--version'], logs / 'vsearch-version.log')
        report['vsearch_version'] = (logs / 'vsearch-version.log').read_text(errors='replace').strip()
        # Fresh owned staging area; publish no stable outputs until all merges and
        # filters succeed. Temporary files never become implicit checkpoints.
        with tempfile.TemporaryDirectory(prefix='.illumina-pe-', dir=output.parent) as temporary:
            work, paths = Path(temporary), {}
            for index, sample in enumerate(sorted(grouped)):
                report['current_sample'] = sample
                folder = work / str(index)
                folder.mkdir()
                paths[sample] = {name: folder / (name + '.fastq')
                                 for name in ('r1', 'pair_r1', 'pair_r2', 'merged', 'filtered')}
                metrics = dict(raw_pairs=raw_counts.get(sample), primer_trimmed_pairs=0,
                               r1_candidates=0, r2_candidates=0, orphan_r1=0,
                               merge_input_pairs=0, merged_pairs=None, merge_fraction=None,
                               merge_status='not_started',
                               selected_input_reads=None, quality_filtered_reads=None,
                               final_reads=None, exclusion_reason=None)
                report['samples'][sample] = metrics
                histograms = {end: dict(before=Counter(), after=Counter(),
                                       trimmed_reads=0, removed_bases=0, empty_reads=0)
                              for end in ('r1', 'r2')}
                with open(paths[sample]['r1'], 'w') as all_r1, \
                     open(paths[sample]['pair_r1'], 'w') as paired_r1, \
                     open(paths[sample]['pair_r2'], 'w') as paired_r2:
                    for row in grouped[sample]:
                        for first, second in paired_records(row['r1'], row['r2']):
                            trimmed = (trim_tail(first), trim_tail(second))
                            metrics['primer_trimmed_pairs'] += 1
                            for end, original, result in zip(('r1', 'r2'), (first, second), trimmed):
                                before, after = len(original[1]), len(result[1])
                                h = histograms[end]
                                h['before'][before] += 1
                                h['after'][after] += 1
                                h['trimmed_reads'] += before != after
                                h['removed_bases'] += before - after
                                h['empty_reads'] += after == 0
                            good1, good2 = bool(trimmed[0][1]), bool(trimmed[1][1])
                            metrics['r1_candidates'] += good1
                            metrics['r2_candidates'] += good2
                            metrics['orphan_r1'] += good1 and not good2
                            if good1:
                                write_record(all_r1, trimmed[0])
                            if good1 and good2:
                                write_record(paired_r1, trimmed[0])
                                write_record(paired_r2, trimmed[1])
                                metrics['merge_input_pairs'] += 1
                metrics['tail_trimming'] = histograms
                log = logs / f'{index:05d}-merge.log'
                metrics['merge_log'] = str(log)
                if metrics['merge_input_pairs']:
                    metrics['merge_status'] = 'running'
                    execute(['vsearch', '--fastq_mergepairs', paths[sample]['pair_r1'],
                             '--reverse', paths[sample]['pair_r2'], '--fastq_qmax', '93',
                             '--fastq_qmaxout', '93', '--fastqout', paths[sample]['merged'],
                             '--threads', str(threads)], log)
                    metrics['merged_pairs'] = sum(1 for _ in records(paths[sample]['merged']))
                    if metrics['merged_pairs'] > metrics['merge_input_pairs']:
                        raise ValueError('Merge output exceeds its paired input')
                    metrics['merge_fraction'] = metrics['merged_pairs'] / metrics['merge_input_pairs']
                    metrics['merge_failure_reasons'] = parse_merge_reasons(log)
                    metrics['merge_status'] = 'completed'
                else:
                    paths[sample]['merged'].touch()
                    metrics['merged_pairs'] = 0
                    metrics['merge_status'] = 'no_usable_pairs'
                    log.write_text('No nonempty paired reads; merge not attempted.\n')
                    metrics['merge_failure_reasons'] = {}
                metrics['unmerged_pairs'] = metrics['merge_input_pairs'] - metrics['merged_pairs']
                save_report(report, report_path)

            total_pairs = sum(r['merge_input_pairs'] for r in report['samples'].values())
            total_merged = sum(r['merged_pairs'] for r in report['samples'].values())
            total_r1 = sum(r['r1_candidates'] for r in report['samples'].values())
            report['route'], report['merge_fraction'], report['reason'] = decide_route(
                total_merged, total_pairs, total_r1, merge_min)
            report.update(route_reason=report['reason'], merge_input_pairs=total_pairs,
                          merged_pairs=total_merged, r1_candidates=total_r1,
                          raw_pairs=sum(raw_counts.values()) if raw_counts else None)
            report['merged_fraction_of_raw_pairs'] = (
                total_merged / report['raw_pairs'] if report['raw_pairs'] else None)
            audit(report, [('vsearch_tail_r1_reads', 'r1_candidates', 'primer_trimmed_reads'),
                           ('vsearch_tail_r2_reads', 'r2_candidates', 'primer_trimmed_reads'),
                           ('vsearch_merge_input_pairs', 'merge_input_pairs', 'vsearch_tail_r1_reads'),
                           ('vsearch_merge_attempt_reads', 'merged_pairs', 'vsearch_merge_input_pairs')])
            if report['route'] is None:
                report['status'] = 'skipped'
                save_report(report, report_path)
                return 98
            # No EE filtering starts until every sample's merge attempt succeeds.
            for index, (sample, metrics) in enumerate(report['samples'].items()):
                report['current_sample'] = sample
                key = 'merged' if report['route'] == 'merged' else 'r1'
                metrics['selected_input_reads'] = metrics['merged_pairs'] if key == 'merged' else metrics['r1_candidates']
                log = logs / f'{index:05d}-filter.log'
                metrics['filter_log'] = str(log)
                if metrics['selected_input_reads']:
                    execute(['vsearch', '--fastq_filter', paths[sample][key], '--fastq_maxee', str(maxee),
                             '--fastq_qmax', '93', '--fastqout', paths[sample]['filtered'],
                             '--threads', str(threads)], log)
                    metrics['quality_filtered_reads'] = sum(1 for _ in records(paths[sample]['filtered']))
                    if metrics['quality_filtered_reads'] > metrics['selected_input_reads']:
                        raise ValueError('Filter output exceeds its selected input')
                else:
                    paths[sample]['filtered'].touch()
                    log.write_text('Selected project route has no input for this sample.\n')
                    metrics['quality_filtered_reads'] = 0
                if metrics['quality_filtered_reads'] == 0:
                    metrics['exclusion_reason'] = ('no_reads_in_project_route' if not metrics['selected_input_reads']
                                                   else 'no_reads_after_ee_filter')
            parent = 'vsearch_merge_attempt_reads' if report['route'] == 'merged' else 'vsearch_tail_r1_reads'
            audit(report, [('vsearch_preprocess_input_reads', 'selected_input_reads', parent),
                           ('vsearch_quality_filtered_reads', 'quality_filtered_reads', 'vsearch_preprocess_input_reads')])
            report['quality_filtered_reads'] = sum(r['quality_filtered_reads'] for r in report['samples'].values())
            report['samples_after_ee'] = sum(r['quality_filtered_reads'] > 0 for r in report['samples'].values())
            if not report['quality_filtered_reads']:
                report.update(status='skipped', reason='no_reads_after_ee_filter')
                save_report(report, report_path)
                return 98
            for sample, metrics in report['samples'].items():
                if metrics['quality_filtered_reads']:
                    os.replace(paths[sample]['filtered'], output / (sample + '.fastq'))
            report.update(status='preprocessed', current_sample=None)
            save_report(report, report_path)
        print(f"Illumina PE: route={report['route']}; merged={total_merged}/{total_pairs}; "
              f"maxee={maxee}; retained={report['quality_filtered_reads']}")
        return 0
    except Exception as error:
        report.update(status='failed', error=str(error))
        current = report['samples'].get(report.get('current_sample'), {})
        if current.get('merge_status') == 'running':
            current['merge_status'] = 'failed'
        save_report(report, report_path)
        raise


def finalize(report_path, table_path):
    from read_counts import table_counts
    report = json.loads(Path(report_path).read_text())
    if report['status'] not in ('preprocessed', 'completed'):
        raise ValueError('Cannot finalize an incomplete preprocessing attempt')
    counts = table_counts(table_path)
    if set(counts) - set(report['samples']):
        raise ValueError('Final table contains samples outside this preprocessing attempt')
    for sample, row in report['samples'].items():
        row['final_reads'] = counts.get(sample, 0)
        if row['final_reads'] > row['quality_filtered_reads']:
            raise ValueError('Final abundance exceeds quality-filtered input')
        if not row['final_reads'] and not row['exclusion_reason']:
            row['exclusion_reason'] = 'no_reads_after_downstream_processing'
    report.update(status='completed', final_reads=sum(counts.values()),
                  final_retained_samples=sum(value > 0 for value in counts.values()),
                  final_table=str(Path(table_path).resolve()))
    save_report(report, report_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    preprocess = sub.add_parser('preprocess')
    preprocess.add_argument('--input', required=True)
    preprocess.add_argument('--output', required=True)
    preprocess.add_argument('--report', required=True)
    preprocess.add_argument('--threads', type=int, default=1)
    preprocess.add_argument('--maxee', type=float, default=1.5)
    preprocess.add_argument('--merge-min', type=float, default=0.5)
    finish = sub.add_parser('finalize')
    finish.add_argument('--report', required=True)
    finish.add_argument('--table', required=True)
    args = parser.parse_args()
    try:
        if args.action == 'finalize':
            finalize(args.report, args.table)
            return 0
        return process(args.input, args.output, args.report, args.threads, args.maxee, args.merge_min)
    except (OSError, EOFError, ValueError, subprocess.CalledProcessError) as error:
        print(f'Illumina PE preprocessing failed: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
