"""Choose one Ion Torrent EE limit from post-primer read lengths and Q19."""
import argparse
from collections import Counter
import math
import sys

from fastp_checked import records, save_json
from read_layout import files


def maxee_from_median_length(median_length, qscore=19.0):
    """Convert a median length and Phred error probability to an EE limit."""
    if (type(median_length) not in (int, float) or not math.isfinite(median_length)
            or median_length <= 0):
        raise ValueError('Median read length must be positive and finite')
    if type(qscore) not in (int, float) or not math.isfinite(qscore) or not 0 <= qscore <= 93:
        raise ValueError('Q score must be between 0 and 93')
    return median_length * 10.0 ** (-qscore / 10.0)


def estimate_maxee(input_dir, qscore=19.0, override=None):
    """Read every FASTQ once; compute the exact, read-weighted length median.

    This is a dataset-wide constant EE threshold, not a per-read Q cutoff.
    Input is the common post-primer stage, before method-specific trimming.
    """
    maxee_from_median_length(1, qscore)
    if override is not None and (type(override) not in (int, float)
            or not math.isfinite(override) or override < 0):
        raise ValueError('Manual maxee must be nonnegative and finite')
    paths = files(input_dir)
    if not paths:
        raise ValueError(f'No FASTQ files found in {input_dir}')
    if len({path.resolve() for path in paths}) != len(paths):
        raise ValueError('Duplicate FASTQ input paths')
    lengths = Counter()
    inputs = []
    for path in paths:
        count = 0
        for _, length in records(path):
            if length == 0:
                raise ValueError(f'Empty read in FASTQ: {path}')
            lengths[length] += 1
            count += 1
        if count == 0:
            raise ValueError(f'Empty FASTQ file: {path}')
        inputs.append({'path': str(path.absolute()), 'read_count': count})
    total = sum(lengths.values())
    lower_rank, upper_rank = (total - 1) // 2, total // 2
    seen = 0
    lower = None
    for length, count in sorted(lengths.items()):
        seen += count
        if lower is None and seen > lower_rank:
            lower = length
        if seen > upper_rank:
            median = (lower + length) / 2.0
            break
    automatic = maxee_from_median_length(median, qscore)
    return {
        'input_stage': 'post_primer_before_method_specific_trimming',
        'read_count': total, 'median_length': median, 'qscore': qscore,
        'error_probability': 10.0 ** (-qscore / 10.0),
        'automatic_maxee': automatic,
        'maxee': automatic if override is None else override,
        'source': 'automatic' if override is None else 'override',
        'files': inputs,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--method', choices=('vsearch', 'dada2'), required=True)
    parser.add_argument('--maxee', type=float, help='Optional fixed EE override')
    args = parser.parse_args()
    try:
        report = estimate_maxee(args.input, override=args.maxee)
        report['method'] = args.method
        save_json(args.output, report)
    except (ValueError, OSError, EOFError) as error:
        parser.exit(2, f'Ion quality error: {error}\n')
    print(f"Ion Torrent: {report['read_count']} reads; median length "
          f"{report['median_length']:g} bp; Q19-derived EE "
          f"{report['automatic_maxee']:.6f}; using {report['maxee']:.6f} "
          f"({report['source']})", file=sys.stderr)
    print(report['maxee'])


if __name__ == '__main__':
    main()
