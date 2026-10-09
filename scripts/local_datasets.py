"""Read local dataset metadata and validate inputs before registering outputs."""
import argparse
import csv
import json
from pathlib import Path
import re
import sys

from project_primers import _atomic_write, PRIMER
from read_layout import discover, layout, validate_local_source

NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')
RESERVED = {'logs', 'tmp', 'datasets_ID.txt', 'datasets.log', 'summary.csv',
            'pip_dataset_read_counts.csv', 'per_dataset_summary.tsv',
            'local_datasets.json', 'selected_metadata.csv', 'project_primers.json'}
PLATFORMS = {'ILLUMINA', 'LS454', 'ION_TORRENT', 'PACBIO_SMRT', 'OXFORD_NANOPORE'}


def validate_name(name):
    if not isinstance(name, str) or not NAME.fullmatch(name) or name in RESERVED:
        raise ValueError(f'Invalid dataset name {name!r}; use letters, digits, underscores, hyphens or dots, and avoid pipeline filenames')


def _read_json(path):
    try:
        value = json.loads(path.read_text())
        if not isinstance(value, dict):
            raise ValueError('expected an object')
        return value
    except (OSError, ValueError) as error:
        raise ValueError(f'Cannot verify existing dataset ownership in {path}: {error}') from error


def ownership_records(target):
    marker = target / '.local-source.json'
    if marker.exists():
        record = _read_json(marker)
        if record.get('source_kind') != 'local' or not isinstance(record.get('source_directory'), str):
            raise ValueError(f'Invalid local dataset ownership record: {marker}')
        yield record
    for path in [target / '.checkpoint.json', *sorted(target.glob('*-run.json'))]:
        if path.exists():
            record = _read_json(path)
            kind = record.get('source_kind')
            if kind is None and 'local_source' in record:
                source = record['local_source']
                if (not isinstance(record.get('inputs'), str) or not isinstance(source, list)
                        or any(not isinstance(row, list) or len(row) != 3
                               or not isinstance(row[0], str) or type(row[1]) is not int
                               or type(row[2]) is not int for row in source)):
                    raise ValueError(f'Invalid legacy dataset source record: {path}')
                kind = 'local' if source else 'archive'
                record = dict(record, source_kind=kind)
            if kind in ('local', 'archive'):
                yield record
            elif kind is not None:
                raise ValueError(f'Invalid dataset source kind in {path}')


def check_ownership(target, source=None):
    """Allow reruns of one source, but never reuse a name for another source."""
    if target.is_symlink():
        raise ValueError(f'Dataset output must not be a symlink: {target}')
    if target.exists() and not target.is_dir():
        raise ValueError(f'Dataset output is already a file: {target}')
    records = list(ownership_records(target)) if target.is_dir() else []
    if source is None:
        if any(record['source_kind'] == 'local' for record in records):
            raise ValueError(f'Dataset name {target.name} already belongs to local data; rename the local dataset or run from another working directory')
        return
    source = Path(source).resolve()
    recorded_directories = [Path(record['source_directory']).resolve() for record in records
                            if record.get('source_directory')]
    for record in records:
        if record['source_kind'] != 'local':
            raise ValueError(f'Local dataset name {target.name} conflicts with an existing online project')
        previous = record.get('source_directory')
        if previous:
            same = Path(previous).resolve() == source
        else:
            # A stored directory is authoritative for symlinked FASTQ inputs.
            # Without one, all legacy paths must demonstrate the same parent;
            # one shared file cannot establish ownership of an entire dataset.
            previous_files = {str(Path(item[0]).resolve()) for item in record.get('local_source', [])
                              if isinstance(item, list) and item and isinstance(item[0], str)}
            same = (all(directory == source for directory in recorded_directories)
                    if recorded_directories else
                    bool(previous_files) and all(Path(path).parent == source for path in previous_files))
        if not same:
            raise ValueError(f'Local dataset name {target.name} already belongs to a different input directory')
    if target.is_dir() and any(target.iterdir()) and not records:
        raise ValueError(f'Existing output {target} has no verifiable source record; choose another dataset name or run from another working directory')


def _source_path(path):
    path = Path(path).resolve()
    if '\t' in str(path) or '\n' in str(path) or '\r' in str(path):
        raise ValueError('Local input paths must not contain tabs or newlines')
    return path


def _primers(name, forward, reverse):
    if reverse and not forward:
        raise ValueError(f'Dataset {name}: a reverse primer requires a forward primer in the same row')
    for value in (forward, reverse):
        if value and not PRIMER.fullmatch(value):
            raise ValueError(f'Dataset {name}: invalid primer sequence {value!r}; use IUPAC DNA bases, not names or paths')
    return dict(forward=forward.upper(), reverse=reverse.upper())


def _online_inputs(metadata, col_bioproject, col_sra=None):
    projects, samples = {}, set()
    with Path(metadata).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError('Online metadata has duplicate column names')
        if col_sra and col_bioproject == col_sra:
            raise ValueError('Online project and Run column names must differ')
        for column in (col_bioproject, col_sra) if col_sra else (col_bioproject,):
            if not column or column not in (reader.fieldnames or []):
                raise ValueError(f'Online metadata is missing column {column!r}')
        for row in reader:
            if None in row:
                raise ValueError(f'Online metadata row {reader.line_num} has more values than column names')
            # Match py_16s._clean_id_series when constructing online sample IDs.
            name = re.sub('[ \n\t]', '', row.get(col_bioproject) or '')
            if not name or name.lower() in {'na', 'nan', 'none', 'null'}:
                continue
            validate_name(name)
            projects[name] = None
            if col_sra:
                run = re.sub('[ \n\t]', '', row.get(col_sra) or '')
                if run and run.lower() not in {'na', 'nan', 'none', 'null'}:
                    samples.add(f'{name}_{run}')
    return list(projects), samples


def plan(metadata, output, col_datasets='datasets', col_path='path', col_platform='platform',
         col_primer_f=None, col_primer_r=None, online_metadata=None, col_bioproject=None,
         col_sra=None):
    metadata = _source_path(metadata)
    output = Path(output).resolve()
    for source_metadata in [metadata, *([_source_path(online_metadata)] if online_metadata else [])]:
        if output == source_metadata or output in source_metadata.parents:
            raise ValueError('Input metadata must be outside the PIP output directory to prevent generated files from overwriting it')
    if col_primer_r and not col_primer_f:
        raise ValueError('--local-primer-r-colNAME requires --local-primer-f-colNAME')
    columns = [col_datasets, col_path, col_platform]
    if not all(columns):
        raise ValueError('Local dataset, path and platform column names must not be empty')
    columns += [column for column in (col_primer_f, col_primer_r) if column]
    if len(columns) != len(set(columns)):
        raise ValueError('Local column names for different roles must differ')
    datasets = {}
    with metadata.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError('Local metadata has duplicate column names')
        for column in columns:
            if column not in fields:
                raise ValueError(f'Local metadata is missing column {column!r}')
        for row in reader:
            if None in row:
                raise ValueError(f'Local metadata row {reader.line_num} has more values than column names')
            if not any((value or '').strip() for value in row.values()):
                continue
            name, source, platform = [(row.get(column) or '').strip() for column in
                                      (col_datasets, col_path, col_platform)]
            if not all((name, source, platform)):
                raise ValueError(f'Local metadata row {reader.line_num}: dataset, path and platform are required')
            validate_name(name)
            if name in datasets:
                raise ValueError(f'Duplicate local dataset name {name!r}')
            platform = platform.upper()
            if platform not in PLATFORMS:
                raise ValueError(f'Dataset {name}: unsupported platform {platform!r}; choose one of {", ".join(sorted(PLATFORMS))}')
            source = Path(source).expanduser()
            source = _source_path(source if source.is_absolute() else metadata.parent / source)
            if not source.is_dir():
                raise ValueError(f'Dataset {name}: FASTQ directory does not exist: {source}')
            forward = (row.get(col_primer_f) or '').strip() if col_primer_f else ''
            reverse = (row.get(col_primer_r) or '').strip() if col_primer_r else ''
            datasets[name] = dict(source=str(source), platform=platform,
                                  **_primers(name, forward, reverse))
    if not datasets:
        raise ValueError('Local metadata contains no datasets')
    online_projects, online_samples = ([], set())
    if online_metadata:
        online_projects, online_samples = _online_inputs(online_metadata, col_bioproject, col_sra)
        for name in online_projects:
            if name in datasets:
                raise ValueError(f'Local dataset name {name!r} conflicts with an online project in the same command')
            check_ownership(output / name)
    destinations = [output / name for name in [*datasets, *online_projects]]
    sample_owners, source_owners = {}, {}
    for name, item in datasets.items():
        source = Path(item['source'])
        if source == output or source in output.parents:
            raise ValueError('Local output must be outside the input directory')
        if source in source_owners:
            raise ValueError(f'Datasets {source_owners[source]} and {name} refer to the same input directory')
        source_owners[source] = name
        for destination in destinations:
            validate_local_source(source, destination)
            if source == destination.resolve() or source in destination.resolve().parents:
                raise ValueError(f'Local input {source} overlaps pipeline output {destination}')
        check_ownership(output / name, source)
        rows = discover(source)
        read_layout = layout(rows)
        if item['platform'] != 'ILLUMINA' and read_layout == 'PE':
            raise ValueError(f'Dataset {name}: paired-end data are unsupported for platform {item["platform"]}')
        samples = sorted({row['sample'] for row in rows})
        for sample in samples:
            if sample in sample_owners:
                raise ValueError(f'Duplicate sample ID {sample!r} in datasets {sample_owners[sample]} and {name}; rename FASTQ files so sample IDs are unique')
            if sample in online_samples:
                raise ValueError(f'Local sample ID {sample!r} conflicts with an online sample in the same command; rename the local FASTQ files')
            sample_owners[sample] = name
        item.update(samples=samples, layout=read_layout)
    return dict(schema_version=2, source_metadata=str(metadata), output_directory=str(output),
                columns=dict(datasets=col_datasets, path=col_path, platform=col_platform,
                             primer_f=col_primer_f, primer_r=col_primer_r),
                dataset_order=list(datasets), datasets=datasets)


def prepare(metadata, output, **columns):
    record = plan(metadata, output, **columns)
    directory = Path(record['output_directory'])
    directory.mkdir(parents=True, exist_ok=True)
    # Only write after every local row and online source has passed preflight.
    for name in record['dataset_order']:
        target = directory / name
        target.mkdir(exist_ok=True)
        marker = dict(schema_version=1, source_kind='local',
                      source_directory=record['datasets'][name]['source'])
        _atomic_write(target / '.local-source.json', lambda handle, value=marker: handle.write(json.dumps(value, indent=2) + '\n'))
    manifest = directory / 'local_datasets.json'
    _atomic_write(manifest, lambda handle: handle.write(json.dumps(record, indent=2) + '\n'))
    return manifest


def read_manifest(path):
    record = _read_json(Path(path))
    if record.get('schema_version') != 2:
        raise ValueError('Unsupported local dataset manifest version; rerun with --local-m to regenerate it')
    order, datasets = record.get('dataset_order'), record.get('datasets')
    if (not isinstance(order, list) or not order or not isinstance(datasets, dict)
            or any(not isinstance(name, str) for name in order)
            or len(set(order)) != len(order) or set(order) != set(datasets)):
        raise ValueError('Invalid local dataset manifest')
    for name in order:
        validate_name(name)
        item = datasets[name]
        item['source'] = str(_source_path(item['source']))
        if item['platform'] not in PLATFORMS:
            raise ValueError(f'Dataset {name}: unsupported platform in manifest')
        item.update(_primers(name, item['forward'], item['reverse']))
    return record


def check_online(metadata, col_bioproject, output):
    projects, _ = _online_inputs(metadata, col_bioproject)
    for name in projects:
        check_ownership(Path(output).resolve() / name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prepare_parser = commands.add_parser('prepare')
    prepare_parser.add_argument('--metadata', required=True)
    prepare_parser.add_argument('--output', required=True)
    prepare_parser.add_argument('--col-datasets', default='datasets')
    prepare_parser.add_argument('--col-path', default='path')
    prepare_parser.add_argument('--col-platform', default='platform')
    prepare_parser.add_argument('--col-primer-f')
    prepare_parser.add_argument('--col-primer-r')
    prepare_parser.add_argument('--online-metadata')
    prepare_parser.add_argument('--col-bioproject')
    prepare_parser.add_argument('--col-sra')
    emit_parser = commands.add_parser('emit')
    emit_parser.add_argument('--manifest', required=True)
    online_parser = commands.add_parser('check-online')
    online_parser.add_argument('--metadata', required=True)
    online_parser.add_argument('--col-bioproject', required=True)
    online_parser.add_argument('--output', required=True)
    args = parser.parse_args()
    try:
        if args.command == 'prepare':
            print(prepare(args.metadata, args.output, col_datasets=args.col_datasets,
                          col_path=args.col_path, col_platform=args.col_platform,
                          col_primer_f=args.col_primer_f, col_primer_r=args.col_primer_r,
                          online_metadata=args.online_metadata, col_bioproject=args.col_bioproject,
                          col_sra=args.col_sra))
        elif args.command == 'emit':
            record = read_manifest(args.manifest)
            for name in record['dataset_order']:
                item = record['datasets'][name]
                print(name, item['source'], item['platform'], item['forward'], item['reverse'], sep='\t')
        else:
            check_online(args.metadata, args.col_bioproject, args.output)
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as error:
        parser.exit(2, f'Local dataset input error: {error}\n')


if __name__ == '__main__':
    main()
