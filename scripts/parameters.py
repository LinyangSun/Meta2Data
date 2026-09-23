"""Validated defaults shared by the shell entry points and Python helpers."""
import argparse
import json
import math
import os
from pathlib import Path
import shlex

# section.key: default, environment variable, lower bound, upper bound
SPEC = {
    'adapter_guard.enabled': (False, 'ADAPTER_GUARD_ENABLED', None, None),
    'adapter_guard.min_length': (50, 'M2D_ADAPTER_GUARD_MIN_LENGTH', 1, None),
    'adapter_guard.min_identity': (98.0, 'M2D_ADAPTER_GUARD_MIN_IDENTITY', 0, 100),
    'adapter_guard.min_coverage': (95.0, 'M2D_ADAPTER_GUARD_MIN_COVERAGE', 0, 100),
    'primer.window': (20, 'PRIMER_WINDOW', 1, 100),
    'primer.fold_threshold': (16, 'PRIMER_FOLD_THRESHOLD', 1, None),
    'primer.support_frequency': (0.10, 'PRIMER_SUPPORT_FREQUENCY', 0.001, 0.25),
    'primer.database_identity': (0.85, 'PRIMER_DATABASE_IDENTITY', 0, 1),
    'primer.informative_fraction': (0.50, 'PRIMER_INFORMATIVE_FRACTION', 0.01, 1),
    'primer.unknown_trim_length': (20, 'PRIMER_UNKNOWN_TRIM_LENGTH', 1, 100),
    'primer.skip_unknown': (False, 'SKIP_UNKNOWN_PRIMERS', None, None),
    'primer.min_length': (50, 'PRIMER_MIN_LENGTH', 1, None),
    'primer.min_average_quality': (20.0, 'PRIMER_MIN_AVERAGE_QUALITY', 0, 93),
    'primer.min_complexity': (0.3, 'PRIMER_MIN_COMPLEXITY', 0, 1),
    'primer.min_entropy': (1.0, 'PRIMER_MIN_ENTROPY', 0, 2),
    'vsearch.maxee': (1.0, 'VSEARCH_MAXEE', 0, None),
    'vsearch.merge_min_fraction': (0.5, 'VSEARCH_MERGE_MIN', 0, 1),
    'vsearch.precluster_identity': (0.99, 'VSEARCH_PRECLUSTER_IDENTITY', 0.01, 1),
    'vsearch.cluster_identity': (0.97, 'VSEARCH_CLUSTER_IDENTITY', 0.01, 1),
    'vsearch.minsize': (2, 'VSEARCH_MINSIZE', 1, None),
    'vsearch.min_frequency': (2, 'VSEARCH_MIN_FREQUENCY', 1, None),
    'vsearch.degraded_trim_left': (0, 'VSEARCH_DEGRADED_TRIM_LEFT', 0, None),
    'vsearch.min_length': (50, 'VSEARCH_MIN_LENGTH', 1, None),
    'vsearch.max_n': (1, 'VSEARCH_MAX_N', 0, None),
    'vsearch.ion_maxee': (2.0, 'ION_VSEARCH_MAXEE', 0, None),
    'vsearch.ion_trim_left': (0, 'ION_VSEARCH_STRIPLEFT', 0, None),
    'vsearch.pacbio_maxee_rate': (0.01, 'VSEARCH_MAXEE_RATE', 0, 1),
    'vsearch.pacbio_min_length': (1000, 'PACBIO_VSEARCH_MINLEN', 1, None),
    'vsearch.pacbio_max_length': (2000, 'PACBIO_VSEARCH_MAXLEN', 1, None),
    'dada2.ion_trim_left': (0, 'DADA2_ION_TRIM_LEFT', 0, None),
    'dada2.pacbio_min_length': (1000, 'DADA2_PACBIO_MIN_LENGTH', 1, None),
    'dada2.pacbio_max_length': (1600, 'DADA2_PACBIO_MAX_LENGTH', 1, None),
    'dada2.quality_trim_score': (25, 'DADA2_QUALITY_TRIM_SCORE', 0, 93),
    'ont.quality': (10, 'ONT_QUALITY', 0, 93),
    'ont.length_tolerance': (0.15, 'ONT_LENGTH_TOLERANCE', 0, 1),
    'ont.length_floor': (200, 'ONT_LENGTH_FLOOR', 1, None),
    'ont.cluster_identity': (0.97, 'ONT_VSEARCH_IDENTITY', 0.01, 1),
    'ont.map_identity': (0.90, 'ONT_MAP_IDENTITY', 0.01, 1),
    'taxa.confidence': (0.7, 'CONFIDENCE', 0, 1),
    'taxa.singleV': (False, 'SINGLE_V', None, None),
    'taxa.notree': (False, 'NOTREE', None, None),
}


def validate(key, value):
    default, _, lo, hi = SPEC[key]
    if isinstance(default, bool):
        valid = isinstance(value, bool)
    elif isinstance(default, int):
        valid = type(value) is int
    else:
        valid = type(value) in (int, float) and math.isfinite(value)
    if not valid or (lo is not None and value < lo) or (hi is not None and value > hi):
        raise ValueError(f'Invalid parameter {key}: {value!r}')
    return value


def load(path=None, environment=False):
    values = {key: spec[0] for key, spec in SPEC.items()}
    if path:
        data = json.loads(Path(path).read_text())
        if not isinstance(data, dict):
            raise ValueError('Parameter file must contain a JSON object')
        for section, entries in data.items():
            if not isinstance(entries, dict) or section not in {k.split('.')[0] for k in SPEC}:
                raise ValueError(f'Unknown parameter section: {section}')
            for name, value in entries.items():
                key = f'{section}.{name}'
                if key not in SPEC:
                    raise ValueError(f'Unknown parameter: {key}')
                values[key] = validate(key, value)
    if environment:
        for key, (default, env, _, _) in SPEC.items():
            if env in os.environ:
                values[key] = validate(key, json.loads(os.environ[env]))
    for group in ('vsearch.pacbio', 'dada2.pacbio'):
        if values[group + '_min_length'] > values[group + '_max_length']:
            raise ValueError(f'{group}: min_length exceeds max_length')
    if values['taxa.singleV'] and values['taxa.notree']:
        raise ValueError('taxa.singleV and taxa.notree are mutually exclusive')
    return values


def nested(values):
    out = {}
    for key, value in values.items():
        section, name = key.split('.')
        out.setdefault(section, {})[name] = value
    return out


def current():
    return load(os.environ.get('PARAMETER_FILE') or None, environment=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['shell', 'save', 'defaults'])
    parser.add_argument('--parameter')
    parser.add_argument('--output')
    args = parser.parse_args()
    try:
        values = load(args.parameter, environment=args.action == 'save')
        if args.action == 'shell':
            for key, value in values.items():
                print(f'export {SPEC[key][1]}={shlex.quote(json.dumps(value))}')
        else:
            content = json.dumps(nested(values), indent=2) + '\n'
            if args.output:
                Path(args.output).write_text(content)
            else:
                print(content, end='')
    except (ValueError, OSError) as error:
        parser.exit(2, f'Parameter error: {error}\n')
