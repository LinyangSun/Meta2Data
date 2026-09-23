"""One FASTQ naming rule for discovery, staging, manifests and read counts."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import re
import shutil

EXT = re.compile(r'\.(?:fastq|fq)(?:\.gz)?$', re.I)
PAIR = re.compile(r'^(?P<sample>.+?)(?:_L(?P<lane>\d{3}))?_(?:R)?(?P<read>[12])(?:_(?P<chunk>\d{3}))?$', re.I)
HINT = 'Rename paired files to sample_R1.fastq.gz / sample_R2.fastq.gz (or sample_1_001.fastq.gz / sample_2_001.fastq.gz); use sample.fastq.gz for single-end data.'


def files(directory):
    return sorted(p for p in Path(directory).iterdir() if p.is_file() and EXT.search(p.name))


def validate_local_source(source, dataset):
    """Keep original reads outside directories owned by the pipeline."""
    source = Path(source).resolve()
    dataset = Path(dataset).resolve()
    candidates = [source] + [path.resolve() for path in files(source)]
    if any(path == dataset or dataset in path.parents for path in candidates):
        raise ValueError('Local input overlaps the pipeline dataset directory; choose a separate output folder outside the original reads.')


def discover(directory, allow_single_r1=False, paths=None):
    groups = {}
    for path in files(directory) if paths is None else sorted(map(Path, paths)):
        stem = EXT.sub('', path.name)
        match = PAIR.fullmatch(stem)
        if match:
            sample, lane, read, chunk = (match.group(k) or '' for k in ('sample', 'lane', 'read', 'chunk'))
        else:
            if re.search(r'(?:^|[_.-])R[12](?:[_.-]|$)', stem, re.I):
                raise ValueError(f'Filename does not match a supported pairing pattern: {path.name}. {HINT}')
            sample, lane, read, chunk = stem, '', '', ''
        if any(c.isspace() for c in sample):
            raise ValueError(f'Sample name contains whitespace: {sample!r}. Rename the file using underscores.')
        key = (sample, lane, chunk)
        group = groups.setdefault(key, {})
        if read in group:
            raise ValueError(f'Ambiguous filenames for {sample}: {group[read]} and {path}. {HINT}')
        group[read] = str(path.absolute())
    if not groups:
        raise ValueError(f'No FASTQ files found in {directory}')
    rows = []
    for (sample, lane, chunk), group in sorted(groups.items()):
        if set(group) == {'1', '2'}:
            layout, r1, r2 = 'PE', group['1'], group['2']
        elif set(group) == {''} or (allow_single_r1 and set(group) == {'1'}):
            layout, r1, r2 = 'SE', next(iter(group.values())), ''
        else:
            raise ValueError(f'Unmatched or ambiguous FASTQ filenames: {list(group.values())}. {HINT}')
        rows.append(dict(sample=sample, lane=lane, chunk=chunk, layout=layout, r1=r1, r2=r2))
    for sample in {r['sample'] for r in rows}:
        if len({r['layout'] for r in rows if r['sample'] == sample}) != 1:
            raise ValueError(f'Sample {sample} has conflicting single/paired-end files. {HINT}')
    return rows


def layout(rows):
    modes = {r['layout'] for r in rows}
    if len(modes) != 1:
        raise ValueError('Mixed single/paired-end samples in one dataset are unsupported. Separate them into folders.')
    return next(iter(modes))


def open_fastq(path, mode='rt'):
    return gzip.open(path, mode) if str(path).lower().endswith('.gz') else open(path, mode)


def stage(rows, destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    for sample in sorted({r['sample'] for r in rows}):
        chunks = [r for r in rows if r['sample'] == sample]
        mode = chunks[0]['layout']
        for direction in ('r1', 'r2') if mode == 'PE' else ('r1',):
            sources = [Path(r[direction]) for r in chunks]
            stem = sample + ('_' + direction[-1] if mode == 'PE' else '')
            compressed = len(sources) > 1 or sources[0].name.lower().endswith('.gz')
            out = destination / (stem + '.fastq' + ('.gz' if compressed else ''))
            if out.exists() or out.is_symlink():
                raise ValueError(f'Staging filename collision: {out}')
            if len(sources) == 1:
                out.symlink_to(sources[0].resolve())
            else:
                with gzip.open(out, 'wb') as dst:
                    for src in sources:
                        with open_fastq(src, 'rb') as fh:
                            shutil.copyfileobj(fh, dst)


def manifest(paths, output, paired, sample_ids=None):
    if not paths:
        raise ValueError('No processed FASTQ files available for the manifest; check pairing and filtering logs.')
    if paired:
        rows = discover(Path(paths[0]).parent, paths=paths)
        if layout(rows) != 'PE' or len({r['sample'] for r in rows}) != len(rows):
            raise ValueError('Paired manifest requires staged, complete sample pairs.')
        values = [[r['sample'], r['r1'], r['r2']] for r in rows]
        header = ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']
    else:
        ids = sample_ids if sample_ids is not None else [EXT.sub('', Path(path).name) for path in paths]
        values = [[sample, str(Path(path).absolute())] for sample, path in zip(ids, paths)]
        header = ['sample-id', 'absolute-filepath']
    if len({v[0] for v in values}) != len(values):
        raise ValueError('Duplicate sample IDs in processed FASTQ files; check filenames.')
    with open(output, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['register', 'stage', 'normalize', 'layout', 'manifest', 'count', 'pairs', 'validate-local'])
    parser.add_argument('--input', required=True)
    parser.add_argument('--output')
    parser.add_argument('--samples')
    parser.add_argument('--list')
    parser.add_argument('--paired', action='store_true')
    args = parser.parse_args()
    try:
        if args.action == 'validate-local':
            validate_local_source(args.input, args.output)
        elif args.action == 'register':
            rows = discover(args.input)
            layout(rows)
            Path(args.output).write_text(json.dumps(rows, indent=2) + '\n')
            Path(args.samples).write_text(''.join(f'{s}\t{s}\n' for s in sorted({r['sample'] for r in rows})))
        elif args.action == 'stage':
            stage(json.loads(Path(args.input).read_text()), args.output)
        elif args.action == 'normalize':
            rows = discover(args.input, allow_single_r1=True)
            layout(rows)
            # Download paths must remain available: move raw downloads aside and stage links.
            raw = Path(args.input).with_name('downloaded_fastq')
            if raw.exists():
                shutil.rmtree(raw)
            Path(args.input).rename(raw)
            for row in rows:
                for key in ('r1', 'r2'):
                    if row[key]:
                        row[key] = str(raw / Path(row[key]).name)
            stage(rows, args.input)
            Path(args.output).write_text(json.dumps(rows, indent=2) + '\n')
        elif args.action == 'pairs':
            rows = discover(args.input)
            if layout(rows) != 'PE':
                raise ValueError('Paired-end input required')
            for row in rows:
                print(row['sample'] + '\t' + row['r1'] + '\t' + row['r2'])
        elif args.action == 'layout':
            print('paired' if layout(discover(args.input)) == 'PE' else 'single')
        elif args.action == 'manifest':
            manifest([str(p) for p in files(args.input)], args.output, args.paired)
        else:
            rows = discover(args.input)
            by_sample = {}
            for row in rows:
                count = 0
                for key in ('r1', 'r2'):
                    if row[key]:
                        with open_fastq(row[key]) as fh:
                            count += sum(1 for _ in fh) // 4
                by_sample[row['sample']] = by_sample.get(row['sample'], 0) + count
            with open(args.output, 'w') as out:
                for line in Path(args.samples).read_text().splitlines():
                    run, sample = line.split('\t')[:2]
                    out.write(f'{run}\t{sample}\t{by_sample.get(sample, 0)}\n')
    except (ValueError, OSError) as error:
        parser.exit(2, f'FASTQ input error: {error}\n')


if __name__ == '__main__':
    main()
