"""Filter 454 reads by each sample's median length and remaining N bases."""
import argparse
from collections import Counter
import gzip
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

from fastp_checked import save_json
from read_layout import EXT, files


_ALLOWED_BASES = frozenset('ACGTRYSWKMBDHVNacgtryswkmbdhvn')


def fastq_records(path):
    """Yield validated four-line records without altering IDs, bases or scores."""
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii', newline='') as stream:
        number = 0
        while True:
            header = stream.readline()
            if not header:
                return
            number += 1
            sequence, plus, quality = (stream.readline() for _ in range(3))
            seq = sequence.rstrip('\r\n')
            qual = quality.rstrip('\r\n')
            label = f'{path}, record {number}'
            if not sequence or not plus or not quality:
                raise ValueError(f'Truncated FASTQ: {label}')
            if (not header.startswith('@') or not header[1:].strip()
                    or header[1].isspace()):
                raise ValueError(f'Invalid FASTQ header: {label}')
            if not plus.startswith('+'):
                raise ValueError(f'Invalid FASTQ separator: {label}')
            if not seq or any(base not in _ALLOWED_BASES for base in seq):
                raise ValueError(f'Empty or invalid FASTQ sequence: {label}')
            if len(seq) != len(qual) or any(not 33 <= ord(base) <= 126 for base in qual):
                raise ValueError(f'Invalid FASTQ quality: {label}')
            yield (header, sequence, plus, quality), seq


def median_length(histogram):
    """Exact read-weighted median; the two middle lengths are averaged."""
    total = sum(histogram.values())
    if total == 0:
        return None
    lower_rank, upper_rank = (total - 1) // 2, total // 2
    seen = 0
    lower = None
    for length, count in sorted(histogram.items()):
        seen += count
        if lower is None and seen > lower_rank:
            lower = length
        if seen > upper_rank:
            return (lower + length) / 2.0
    raise ValueError('Invalid length histogram')


def _overlap(first, second):
    return first == second or first in second.parents or second in first.parents


def _signature(path):
    info = path.stat()
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns


def _prepare_paths(input_dir, output_dir, report_path):
    source = Path(input_dir).resolve()
    target = Path(output_dir).resolve()
    report = Path(report_path).resolve()
    if not source.is_dir():
        raise ValueError(f'Input FASTQ directory does not exist: {source}')
    if _overlap(source, target):
        raise ValueError('Input and output FASTQ directories must not overlap')
    paths = files(source)
    resolved = [path.resolve() for path in paths]
    if any(path == target or target in path.parents for path in resolved):
        raise ValueError('An input FASTQ resolves inside the output directory')
    identities = [(path.stat().st_dev, path.stat().st_ino) for path in paths]
    if len(set(identities)) != len(paths):
        raise ValueError('Duplicate FASTQ input paths')
    if report == source or source in report.parents or report in resolved:
        raise ValueError('Report must not overwrite input reads or reside in the input directory')
    if report == target or report in target.parents:
        raise ValueError('Report path conflicts with the output directory')
    if report.exists() and report.is_dir():
        raise ValueError('Report path must be a file')
    if report.exists() and (report.stat().st_dev, report.stat().st_ino) in identities:
        raise ValueError('Report must not overwrite input reads')
    names = [EXT.sub('', path.name) + '.fastq'
             + ('.gz' if path.name.lower().endswith('.gz') else '') for path in paths]
    stems = [EXT.sub('', path.name) for path in paths]
    if len(set(stems)) != len(stems) or len(set(names)) != len(names):
        raise ValueError('FASTQ files have duplicate sample names')
    if any(report == target / name for name in names):
        raise ValueError('Report path conflicts with a FASTQ output')
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        raise ValueError('Output directory must be empty to avoid mixing previous results')
    return paths, target, report, names


def preprocess(input_dir, output_dir, report_path, length_fraction=0.5, max_n=1):
    """Two passes per file: estimate its length limit, then filter every read.

    Length is checked first. ``removed_short`` and ``removed_n`` are exclusive;
    ``removed_both`` is an informational subset of ``removed_short`` and must
    not be added to the loss total. N and n are counted together; reads with
    at most max_n are retained. No Phred or expected-error filter is used.
    """
    if (type(length_fraction) not in (int, float) or not math.isfinite(length_fraction)
            or not 0 < length_fraction <= 1):
        raise ValueError('Length fraction must be finite and in (0, 1]')
    if type(max_n) is not int or max_n < 0:
        raise ValueError('Maximum N count must be a nonnegative integer')
    paths, target, report_path, names = _prepare_paths(input_dir, output_dir, report_path)
    report = {
        'status': 'running', 'input_stage': 'post_primer_and_adaptive_tail_trim',
        'length_fraction': length_fraction, 'max_n': max_n,
        'n_rule': 'count(N) + count(n) <= max_n',
        'length_rule': 'length >= ceil(per_sample_median_length * length_fraction)',
        'median_basis': 'all input reads in each sample, including reads containing N',
        'filter_order': ['short_length', 'N_count_exceeds_max_n'],
        'quality_filter': False,
        'counting': 'input = removed_short + removed_n + kept; removed_both is the short-and-N-over-limit subset of removed_short; kept_with_n is a subset of kept',
        'files': [], 'totals': {key: 0 for key in
                              ('input', 'removed_short', 'removed_n', 'removed_both', 'kept', 'kept_with_n')},
    }
    temporary = None
    try:
        if not paths:
            raise ValueError('No FASTQ files found in the input directory')
        # Validate every file before publishing outputs or creating the output directory.
        signatures = []
        for path, name in zip(paths, names):
            before = _signature(path)
            lengths = Counter(len(seq) for _, seq in fastq_records(path))
            if _signature(path) != before:
                raise ValueError(f'FASTQ changed while being read: {path}')
            signatures.append(before)
            median = median_length(lengths)
            report['files'].append({
                'sample': EXT.sub('', path.name), 'path': str(path),
                'output': str(target / name), 'input': sum(lengths.values()),
                'median_length': median,
                'min_length': math.ceil(median * length_fraction) if median is not None else None,
                'removed_short': 0, 'removed_n': 0, 'removed_both': 0, 'kept': 0,
                'kept_with_n': 0,
                'status': 'empty_input' if median is None else 'pending',
            })
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = Path(tempfile.mkdtemp(prefix='.ls454-quality-', dir=target.parent))
        for path, name, before, entry in zip(paths, names, signatures, report['files']):
            if _signature(path) != before:
                raise ValueError(f'FASTQ changed between filtering passes: {path}')
            opener = gzip.open if name.endswith('.gz') else open
            counted = 0
            with opener(temporary / name, 'wt', encoding='ascii', newline='') as output:
                for raw_record, seq in fastq_records(path):
                    counted += 1
                    n_count = seq.count('N') + seq.count('n')
                    exceeds_n = n_count > max_n
                    if len(seq) < entry['min_length']:
                        entry['removed_short'] += 1
                        entry['removed_both'] += int(exceeds_n)
                    elif exceeds_n:
                        entry['removed_n'] += 1
                    else:
                        output.writelines(raw_record)
                        entry['kept'] += 1
                        entry['kept_with_n'] += int(n_count > 0)
            if counted != entry['input'] or _signature(path) != before:
                raise ValueError(f'FASTQ changed between filtering passes: {path}')
            if entry['input'] != entry['removed_short'] + entry['removed_n'] + entry['kept']:
                raise ValueError(f'Filtering count mismatch: {path}')
            if entry['input']:
                entry['status'] = 'completed' if entry['kept'] else 'all_filtered'
            for key in report['totals']:
                report['totals'][key] += entry[key]
            print(f"454 sample {entry['sample']}: input={entry['input']}; "
                  f"median={entry['median_length']}; min_length={entry['min_length']}; "
                  f"removed_short={entry['removed_short']}; removed_n={entry['removed_n']}; "
                  f"max_n={max_n}; kept={entry['kept']}; kept_with_n={entry['kept_with_n']}",
                  file=sys.stderr)
        if target.exists() and any(target.iterdir()):
            raise ValueError('Output directory changed during processing; refusing to publish')
        target.mkdir(parents=True, exist_ok=True)
        for name in names:
            os.replace(temporary / name, target / name)
        if report['totals']['input'] == 0:
            report['status'] = 'empty_input'
            raise ValueError('All input FASTQ files are empty; no 454 reads to process')
        if report['totals']['kept'] == 0:
            report['status'] = 'all_filtered'
            raise ValueError('All 454 reads were removed by the length/N filters')
        report['status'] = 'completed'
        save_json(report_path, report)
        return report
    except (ValueError, OSError, EOFError) as error:
        if report['status'] == 'running':
            report['status'] = 'failed'
        report['error'] = str(error)
        save_json(report_path, report)
        raise
    finally:
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest='command', required=True)
    command = subparsers.add_parser('preprocess', help='Filter each sample by relative length and N bases')
    command.add_argument('--input', required=True)
    command.add_argument('--output-dir', required=True)
    command.add_argument('--report', required=True)
    command.add_argument('--length-fraction', type=float, default=0.5)
    command.add_argument('--max-n', type=int, default=1,
                         help='Maximum combined N/n count per read (default: 1)')
    args = parser.parse_args()
    try:
        preprocess(args.input, args.output_dir, args.report, args.length_fraction, args.max_n)
    except (ValueError, OSError, EOFError) as error:
        parser.exit(2, f'454 quality error: {error}\n')


if __name__ == '__main__':
    main()
