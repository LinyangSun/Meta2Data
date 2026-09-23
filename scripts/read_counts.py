"""Durable read-count ledger. No QIIME process is needed to read QZA tables.

The wide CSVs are materialized views; immutable event JSON and tool reports live
outside tmp, grouped by execution attempt. NA means unobserved/not applicable,
never zero. Pooled centroid abundances are not per-sample read assignments.
"""
import argparse
import csv
from datetime import datetime, timezone
import fcntl
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile
import uuid
import zipfile

from read_layout import discover, files, open_fastq

COMMON = ['RawReads', 'fastp_reads', 'primer_trimmed_reads', 'sanitized_reads']
DADA2 = ['dada2_qc_input_reads', 'dada2_quality_filtered_reads',
         'dada2_input_reads', 'dada2_primer_matched_reads', 'dada2_filtered_reads',
         'dada2_denoised_reads', 'dada2_denoised_forward_reads',
         'dada2_denoised_reverse_reads', 'dada2_merged_reads',
         'dada2_nonchimeric_reads', 'dada2_final_reads']
VSEARCH = ['vsearch_merge_attempt_reads', 'vsearch_preprocess_input_reads',
           'vsearch_quality_filtered_reads', 'vsearch_length_n_filtered_reads',
           'vsearch_adaptive_trimmed_reads', 'vsearch_chopper_reads',
           'vsearch_preprocessed_reads', 'vsearch_dereplicated_reads',
           'vsearch_abundance_filtered_reads', 'vsearch_preclustered_reads',
           'vsearch_denoised_reads', 'vsearch_nonchimeric_reads',
           'vsearch_ont_denoised_reads', 'vsearch_ont_polished_reads',
           'vsearch_ont_relabeled_reads', 'vsearch_clustered_reads',
           'vsearch_mapping_input_reads', 'vsearch_mapped_reads',
           'vsearch_imported_reads', 'vsearch_final_reads']
PIP_COLUMNS = COMMON + DADA2 + VSEARCH + ['FinalReads']
TAXA_COLUMNS = ['taxa_raw_reads', 'taxa_oriented_reads', 'taxa_classified_reads',
                'taxa_tree_placed_reads', 'taxa_final_reads']
META = ['Mode', 'CountUnit', 'ProcessingPath', 'Attempt', 'Status']


def atomic(path, write):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix='.' + path.name)
    try:
        with os.fdopen(fd, 'w', newline='') as stream:
            write(stream)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_json(path, data):
    atomic(path, lambda stream: json.dump(data, stream, indent=2, sort_keys=True))


def write_csv(path, rows, columns):
    def write(stream):
        writer = csv.DictWriter(stream, columns, lineterminator='\n', restval='NA')
        writer.writeheader()
        writer.writerows(rows)
    atomic(path, write)


def replace_group(path, rows, columns, dataset, mode):
    """Replace this dataset/method only; preserve other writers and old columns."""
    path = Path(path)
    with open(str(path) + '.lock', 'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        old = []
        if path.exists():
            with path.open() as stream:
                reader = csv.DictReader(stream)
                columns = list(dict.fromkeys(columns + (reader.fieldnames or [])))
                old = [r for r in reader if not (r['BioProject'] == dataset and
                       r.get('Mode') in (None, '', 'NA', mode))]
        write_csv(path, old + rows, columns)


def integer(value):
    number = float(value)
    if not math.isfinite(number) or number < 0 or number != int(number):
        raise ValueError(f'Invalid read abundance: {value}')
    return int(number)


def table_counts(path):
    with zipfile.ZipFile(path) as archive:
        names = [n for n in archive.namelist() if n.endswith('/data/feature-table.biom')]
        if len(names) != 1:
            raise ValueError(f'Expected exactly one feature table in {path}')
        payload = archive.read(names[0])
        if payload.lstrip().startswith(b'{'):
            data = json.loads(payload)
            samples = [column['id'] for column in data['columns']]
            if len(set(samples)) != len(samples):
                raise ValueError('Duplicate sample IDs in BIOM')
            counts = dict.fromkeys(samples, 0)
            if data['matrix_type'] == 'sparse':
                for _, column, value in data['data']:
                    counts[samples[column]] += integer(value)
            elif data['matrix_type'] == 'dense':
                for row in data['data']:
                    if len(row) != len(samples):
                        raise ValueError('Ragged BIOM matrix')
                    for sample, value in zip(samples, row):
                        counts[sample] += integer(value)
            else:
                raise ValueError('Unknown BIOM matrix type')
            return counts
        from biom import load_table
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / 'table.biom'
            target.write_bytes(payload)
            table = load_table(str(target))
    # Validate individual cells, not just the total (fractional/negative values).
    for value in table.matrix_data.data:
        integer(value)
    # BIOM Table.data rejects a zero-feature table even when sample IDs remain
    # (e.g. SEPP removed zero features). Sparse columns correctly sum to zero.
    return {str(s): sum(integer(v) for v in table.matrix_data.getcol(i).data)
            for i, s in enumerate(table.ids(axis='sample'))}


def fq_metrics(path):
    count, bases = 0, 0
    with open_fastq(path) as stream:
        while stream.readline():
            seq, plus, qual = stream.readline(), stream.readline(), stream.readline()
            if not seq or not plus or not qual:
                raise ValueError(f'Truncated FASTQ: {path}')
            count += 1
            bases += len(seq.rstrip('\r\n'))
    return count, bases


def fq_count(path):
    return fq_metrics(path)[0]


def fastq_counts(directory, metrics=None):
    result, units = {}, set()
    if not Path(directory).is_dir():
        raise FileNotFoundError(directory)
    if not files(directory):
        return {}, None
    for row in discover(directory, allow_single_r1=True):
        count, bases1 = fq_metrics(row['r1'])
        count2, bases2 = fq_metrics(row['r2']) if row['r2'] else (0, 0)
        if row['r2'] and count2 != count:
            raise ValueError(f'Unequal mate counts for {row["sample"]}')
        if metrics is not None:
            metrics.append(dict(SampleName=row['sample'], lane=row['lane'], chunk=row['chunk'],
                layout=row['layout'], r1_path=row['r1'], r2_path=row['r2'],
                r1_bytes=Path(row['r1']).stat().st_size,
                r2_bytes=Path(row['r2']).stat().st_size if row['r2'] else 0,
                r1_compressed=row['r1'].lower().endswith('.gz'),
                r2_compressed=row['r2'].lower().endswith('.gz') if row['r2'] else False,
                raw_reads=count + count2, raw_read_pairs=count if row['r2'] else 'NA',
                raw_fragments=count, total_bases=bases1 + bases2,
                r1_mean_length=bases1 / count if count else 'NA',
                r2_mean_length=bases2 / count2 if count2 else 'NA',
                measurement_source='actual FASTQ traversal'))
        result[row['sample']] = result.get(row['sample'], 0) + count
        units.add('read_pairs' if row['r2'] else 'reads')
    return result, next(iter(units)) if len(units) == 1 else 'fragments'


def fasta_abundance(path):
    total, features, missing = 0, 0, False
    with open(path) as stream:
        for line in stream:
            if line.startswith('>'):
                match = re.search(r';size=(\d+)(?:;|\s|$)', line)
                if not match:
                    missing = True
                else:
                    total += int(match[1])
                features += 1
    return ('NA' if missing else total), features


def otu_counts(path):
    with open(path) as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next(reader)
        if len(set(header[1:])) != len(header[1:]):
            raise ValueError('Duplicate samples in OTU table')
        counts = dict.fromkeys(header[1:], 0)
        for row in reader:
            if len(row) != len(header):
                raise ValueError('Ragged OTU table')
            for sample, value in zip(header[1:], row[1:]):
                counts[sample] += integer(value)
    return counts


def begin(dataset, output, mode, pipeline='pip'):
    root = Path(dataset) / 'read_counts' / mode
    attempt = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '-' + uuid.uuid4().hex[:10]
    directory = root / attempt
    directory.mkdir(parents=True)
    previous = json.loads((root / 'latest.json').read_text()) if (root / 'latest.json').exists() else None
    samples = {}
    if pipeline == 'pip':
        mapping = Path(dataset) / (Path(dataset).name + '_sra.txt')
        for line in mapping.read_text().splitlines():
            run, sample = line.split('\t')[:2]
            if sample in samples:
                raise ValueError(f'Duplicate sample in SRA map: {sample}')
            samples[sample] = run
    state = dict(dataset=Path(dataset).name, dataset_path=str(Path(dataset).resolve()),
                 output=str(Path(output).resolve()), mode=mode, pipeline=pipeline,
                 attempt=attempt, directory=str(directory.resolve()), previous=previous,
                 samples=samples, projects={}, stages={}, status='running', unit='NA', events=0,
                 execution_id=os.environ.get('M2D_EXECUTION_ID', ''))
    checkpoint = Path(dataset) / '.checkpoint.json'
    state['fingerprint'] = json.loads(checkpoint.read_text()).get('fingerprint') if checkpoint.exists() else None
    if checkpoint.exists():
        shutil.copy2(checkpoint, directory / 'input-parameters.json')
    write_json(directory / 'state.json', state)
    write_json(root / 'latest.json', str((directory / 'state.json').resolve()))
    publish(state)
    return str((directory / 'state.json').resolve())


def publish(state):
    stages = state['stages']
    taxa = state['pipeline'] == 'taxa'
    columns = TAXA_COLUMNS if taxa else PIP_COLUMNS
    ids = ['BioProject', 'SampleName'] if taxa else ['BioProject', 'Run', 'SampleName']
    rows = []
    for sample, run in state['samples'].items():
        row = dict(BioProject=state['projects'].get(sample, state['dataset']), SampleName=sample)
        if not taxa:
            row['Run'] = run
        for column in columns:
            event = stages.get(column)
            row[column] = event['counts'].get(sample, 'NA') if event else 'NA'
        if not taxa:
            row['FinalReads'] = row.get(state['mode'] + '_final_reads', 'NA')
        selected = stages.get('dada2_input_reads', {})
        route = selected.get('execution_layout', 'NA')
        if state['mode'] == 'vsearch':
            route = stages.get('vsearch_preprocess_input_reads', {}).get('sample_sources', {}).get(sample, 'NA')
        row.update(Mode=state['mode'], CountUnit=state['unit'], ProcessingPath=route,
                   Attempt=state['attempt'], Status=state['status'])
        rows.append(row)
    out = Path(state['output'])
    if taxa:
        write_csv(out / 'taxa_read_counts.csv', rows, ids + columns + META)
    else:
        replace_group(out / 'summary.csv', rows, ids + columns + META, state['dataset'], state['mode'])
    aggregate = {'BioProject': state['dataset']}
    for column in columns:
        aggregate[column] = stages[column]['total'] if column in stages else 'NA'
    if not taxa:
        aggregate['FinalReads'] = aggregate.get(state['mode'] + '_final_reads', 'NA')
    aggregate.update(Mode=state['mode'], CountUnit=state['unit'], ProcessingPath='see_sample_rows',
                     Attempt=state['attempt'], Status=state['status'])
    if taxa:
        write_csv(out / 'taxa_total_read_counts.csv', [aggregate], ['BioProject'] + columns + META)
    else:
        replace_group(out / 'pip_dataset_read_counts.csv', [aggregate], ['BioProject'] + columns + META,
                      state['dataset'], state['mode'])
    # Full attempt views survive future runs as well as tmp cleanup.
    write_csv(Path(state['directory']) / 'summary.csv', rows, ids + columns + META)
    write_csv(Path(state['directory']) / 'dataset_summary.csv', [aggregate], ['BioProject'] + columns + META)
    losses = []
    for column, event in stages.items():
        parent = stages.get(event['input_stage'])
        if parent is None or event.get('basis') != parent.get('basis'):
            continue
        for sample in list(state['samples']) + ['__dataset_total__']:
            before = parent['total'] if sample == '__dataset_total__' else parent['counts'].get(sample)
            after = event['total'] if sample == '__dataset_total__' else event['counts'].get(sample)
            if not isinstance(before, int) or not isinstance(after, int):
                continue
            if after > before:
                raise ValueError(f'Count increase across filter {event["input_stage"]} -> {column}: {before} -> {after}')
            losses.append(dict(SampleName=sample, InputStage=event['input_stage'], OutputStage=column,
                               InputReads=before, OutputReads=after, LostReads=before-after,
                               LossFraction=(before-after)/before if before else 'NA', Basis=event['basis']))
    loss_columns = ['SampleName', 'InputStage', 'OutputStage', 'InputReads', 'OutputReads', 'LostReads', 'LossFraction', 'Basis']
    write_csv(Path(state['directory']) / 'read_losses.csv', losses, loss_columns)
    if taxa:
        write_csv(out / 'taxa_read_losses.csv', losses, loss_columns)


def save(state):
    write_json(Path(state['directory']) / 'state.json', state)
    publish(state)


def record(state, stage, counts, source, parent='', basis='fragments', total=None, **details):
    if stage not in PIP_COLUMNS + TAXA_COLUMNS:
        raise ValueError(f'Unknown read-count stage: {stage}')
    normalized = {}
    for sample, count in counts.items():
        if sample not in state['samples']:
            clean = sample.removesuffix('.fastq')
            if clean not in state['samples']:
                raise ValueError(f'Unmapped sample {sample!r} at {stage}')
            sample = clean
        if sample in normalized:
            raise ValueError(f'Ambiguous sample {sample!r} at {stage}')
        normalized[sample] = integer(count)
    # Pooled stages deliberately have NO per-sample assignments.
    if total is None:
        normalized = {s: normalized.get(s, 0) for s in state['samples']}
        total = sum(normalized.values())
    if not parent and stage in state['stages']:
        parent = state['stages'][stage]['input_stage']
    if stage == 'dada2_final_reads' and 'dada2_nonchimeric_reads' in state['stages']:
        parent = 'dada2_nonchimeric_reads'
        if normalized != state['stages'][parent]['counts']:
            raise ValueError('Final DADA2 table disagrees with non-chimeric statistics')
    event = dict(stage=stage, counts=normalized, total=integer(total) if total != 'NA' else 'NA', source=str(source),
                 execution_id=os.environ.get('M2D_EXECUTION_ID', state.get('execution_id', '')),
                 resource_invocation_id=os.environ.get('M2D_PROFILE_INVOCATION_ID', ''),
                 input_stage=parent, basis=basis, time=datetime.now(timezone.utc).isoformat(), **details)
    predecessor = state['stages'].get(parent)
    if predecessor and predecessor['basis'] == basis:
        for sample, count in normalized.items():
            before = predecessor['counts'].get(sample)
            if before is not None and count > before:
                raise ValueError(f'Read-count balance failure at {stage}, {sample}: {before} -> {count}')
        if isinstance(event['total'], int) and isinstance(predecessor['total'], int) and event['total'] > predecessor['total']:
            raise ValueError(f'Read-count balance failure at {stage}: output exceeds input')
    state['stages'][stage] = event
    state['events'] += 1
    write_json(Path(state['directory']) / f'event-{state["events"]:05d}.json', event)
    save(state)


def stats(state, path, kind):
    with zipfile.ZipFile(path) as archive:
        names = [n for n in archive.namelist() if re.search(r'/data/stats\.(tsv|csv)$', n)]
        if len(names) != 1:
            raise ValueError(f'Expected one statistics file in {path}')
        content = archive.read(names[0]).decode()
    report = Path(state['directory']) / f'{kind}-{state["events"]:05d}'
    report.mkdir()
    shutil.copy2(path, report / Path(path).name)
    (report / Path(names[0]).name).write_text(content)
    reader = csv.DictReader(io.StringIO(content), delimiter='\t' if names[0].endswith('.tsv') else ',')
    rows = [r for r in reader if not next(iter(r.values())).startswith('#')]
    fields = reader.fieldnames or []
    if kind == 'dada2':
        if not {'input', 'filtered', 'non-chimeric'} <= set(fields):
            raise ValueError(f'Unknown DADA2 statistics fields: {fields}')
        # A forward-only retry replaces selected stats, keeping old immutable events.
        for key in DADA2[2:]:
            state['stages'].pop(key, None)
        mapping = [('input', 'input'), ('primer-matched', 'primer_matched'),
                   ('primer-removed', 'primer_matched'),
                   ('filtered', 'filtered'), ('denoised', 'denoised'),
                   ('denoised-f', 'denoised_forward'), ('denoised-r', 'denoised_reverse'),
                   ('merged', 'merged'), ('non-chimeric', 'nonchimeric')]
        parent = ''
        for field, name in mapping:
            if field not in fields:
                continue
            stage = f'dada2_{name}_reads'
            # Merging requires BOTH directional denoisers; no single directional
            # count is the unique predecessor. Record no fabricated loss edge.
            edge = parent
            if name == 'input':
                edge = next((s for s in ('dada2_quality_filtered_reads', 'sanitized_reads',
                             'primer_trimmed_reads', 'fastp_reads') if s in state['stages']), '')
            if name == 'merged' and 'denoised-f' in fields:
                edge = ''
            if name.startswith('denoised_'):
                edge = 'dada2_filtered_reads'
            record(state, stage, {r[fields[0]]: r[field] for r in rows}, report, edge,
                   execution_layout='paired' if 'merged' in fields else 'single')
            if name != 'denoised_reverse':
                parent = stage
    else:
        aliases = [('total_input_reads', 'dada2_qc_input_reads'),
                   ('total-input-reads', 'dada2_qc_input_reads'),
                   ('total_retained_reads', 'dada2_quality_filtered_reads'),
                   ('total-retained-reads', 'dada2_quality_filtered_reads')]
        found = set()
        for field, stage in aliases:
            if field in fields:
                record(state, stage, {r[fields[0]]: r[field] for r in rows}, report,
                       'dada2_qc_input_reads' if 'quality' in stage else '')
                found.add(stage)
        if len(found) != 2:
            raise ValueError(f'Unknown quality-filter statistics fields: {fields}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['begin', 'recover', 'fastq', 'table', 'fasta', 'fasta-directory', 'otu', 'stats', 'status', 'inherit', 'value', 'taxa-init', 'check-removed', 'export-stats'])
    parser.add_argument('--state')
    parser.add_argument('--dataset')
    parser.add_argument('--output')
    parser.add_argument('--mode')
    parser.add_argument('--pipeline', default='pip')
    parser.add_argument('--stage')
    parser.add_argument('--input')
    parser.add_argument('--parent', default='')
    parser.add_argument('--sample')
    parser.add_argument('--value', type=int)
    parser.add_argument('--basis', default='fragments')
    parser.add_argument('--suffix', default='_polished.fasta')
    args = parser.parse_args()
    if args.action == 'export-stats':
        with zipfile.ZipFile(args.input) as archive:
            names = [n for n in archive.namelist() if re.search(r'/data/stats\.(tsv|csv)$', n)]
            if len(names) != 1:
                raise ValueError('Expected one statistics file')
            content = archive.read(names[0]).decode()
        rows = list(csv.reader(io.StringIO(content), delimiter='\t' if names[0].endswith('.tsv') else ','))
        atomic(args.output, lambda stream: csv.writer(stream, delimiter='\t', lineterminator='\n').writerows(rows))
        return
    if args.action == 'begin':
        print(begin(args.dataset, args.output, args.mode, args.pipeline))
        return
    if args.action == 'recover':
        latest = Path(args.dataset) / 'read_counts' / args.mode / 'latest.json'
        path = json.loads(latest.read_text()) if latest.exists() else begin(args.dataset, args.output, args.mode)
        state = json.loads(Path(path).read_text())
        # A crash after copying the final table must not leave a failed/partial
        # summary forever merely because the runner skips completed datasets.
        record(state, args.mode + '_final_reads', table_counts(args.input), args.input)
        state['status'] = 'success'
        save(state)
        return
    state = json.loads(Path(args.state).read_text())
    if args.action == 'status':
        state['status'] = args.input
        save(state)
    elif args.action == 'inherit':
        previous = state['previous']
        allowed = {'RawReads', 'fastp_reads'}
        if args.input in ('primer', 'chopper'):
            allowed.add('primer_trimmed_reads')
        if args.input == 'chopper':
            allowed.add('vsearch_chopper_reads')
        found = {}
        while previous:
            old = json.loads(Path(previous).read_text())
            if old.get('fingerprint') != state.get('fingerprint') or old['samples'] != state['samples']:
                break
            for stage in allowed - found.keys():
                if stage in old['stages']:
                    found[stage] = (old['stages'][stage], old['attempt'])
            if old['unit'] != 'NA':
                state['unit'] = old['unit']
            previous = old['previous']
        for stage in COMMON + ['vsearch_chopper_reads']:
            if stage in found:
                event, origin = found[stage]
                record(state, stage, event['counts'], event['source'], event['input_stage'],
                       event['basis'], event['total'], inherited_from=origin)
    elif args.action == 'fastq':
        if args.sample:
            counts, unit = {args.sample: fq_count(args.input)}, None
            existing = state['stages'].get(args.stage, {}).get('counts', {})
            # Do not claim unfinished samples were counted: preserve partial NA.
            counts = {**existing, **counts}
            sources = dict(state['stages'].get(args.stage, {}).get('sample_sources', {}))
            sources[args.sample] = str(args.input)
            record(state, args.stage, counts, args.input, args.parent, args.basis,
                   sum(counts.values()), sample_sources=sources)
        else:
            metrics = [] if args.stage == 'RawReads' else None
            counts, unit = fastq_counts(args.input, metrics)
            if args.stage == 'RawReads':
                state['unit'] = unit or 'reads'
                for row in metrics:
                    row.update(BioProject=state['dataset'], Run=state['samples'].get(row['SampleName'], ''),
                               attempt=state['attempt'], execution_id=state.get('execution_id', ''))
                if metrics:
                    write_csv(Path(state['directory']) / 'sample_sizes.csv', metrics, list(metrics[0]))
            record(state, args.stage, counts, args.input, args.parent, args.basis)
    elif args.action == 'table':
        record(state, args.stage, table_counts(args.input), args.input, args.parent, args.basis)
    elif args.action == 'otu':
        record(state, args.stage, otu_counts(args.input), args.input, args.parent, args.basis)
    elif args.action == 'fasta':
        total, features = fasta_abundance(args.input)
        # Racon can discard size annotations. Do not interpret downstream
        # VSEARCH's default size=1 as recovered original-read abundance.
        upstream = {'vsearch_ont_relabeled_reads': 'vsearch_ont_polished_reads',
                    'vsearch_clustered_reads': 'vsearch_ont_relabeled_reads'}.get(args.stage)
        reported_size_sum = total
        if upstream in state['stages'] and state['stages'][upstream]['total'] == 'NA':
            total = 'NA'
        # Archive labels/sizes as evidence without retaining full FASTQ data.
        evidence = Path(state['directory']) / f'{args.stage}-{state["events"]:05d}.headers.txt'
        with open(args.input) as stream, evidence.open('w') as out:
            out.writelines(line for line in stream if line.startswith('>'))
        record(state, args.stage, {}, evidence, args.parent, 'centroid_abundance', total,
               features=features, reported_size_sum=reported_size_sum,
               unavailable_reason='Missing upstream read-abundance annotations' if total == 'NA' else '')
    elif args.action == 'fasta-directory':
        counts, features, missing = {}, 0, []
        evidence = Path(state['directory']) / f'{args.stage}-{state["events"]:05d}.headers.txt'
        with evidence.open('w') as out:
            for path in sorted(Path(args.input).glob('*' + args.suffix)):
                count, number = fasta_abundance(path)
                sample = path.name[:-len(args.suffix)]
                features += number
                if count == 'NA':
                    missing.append(sample)
                else:
                    counts[sample] = count
                with path.open() as stream:
                    out.writelines(line for line in stream if line.startswith('>'))
        record(state, args.stage, counts, evidence, args.parent, 'centroid_abundance',
               'NA' if missing else sum(counts.values()), features=features,
               missing_abundance_samples=missing)
    elif args.action == 'stats':
        stats(state, args.input, args.stage)
    elif args.action == 'value':
        counts = dict(state['stages'].get(args.stage, {}).get('counts', {}))
        counts[args.sample] = args.value
        record(state, args.stage, counts, args.input or '', args.parent, args.basis, sum(counts.values()))
    elif args.action == 'taxa-init':
        counts = table_counts(args.input)
        collection = json.loads((Path(state['output']) / 'collection.json').read_text())
        projects = {}
        expected = {}
        for dataset in collection['datasets']:
            for sample, count in table_counts(dataset['table']['path']).items():
                if sample in projects:
                    raise ValueError(f'Duplicate TAXA sample: {sample}')
                projects[sample], expected[sample] = dataset['dataset'], count
        if counts != expected:
            raise ValueError('Merged TAXA input does not match dataset sample counts')
        state['samples'] = {s: s for s in counts}
        state['projects'] = projects
        state['unit'] = 'feature_table_abundance'
        record(state, 'taxa_raw_reads', counts, args.input)
        for filename in ('collection.json', 'taxa-run-state.json'):
            path = Path(state['output']) / filename
            if path.exists():
                shutil.copy2(path, Path(state['directory']) / filename)
    elif args.action == 'check-removed':
        removed = table_counts(args.input)
        before, after = state['stages'][args.parent]['counts'], state['stages'][args.stage]['counts']
        if set(removed) - set(before) or any(before[s] != after[s] + removed.get(s, 0) for s in before):
            raise ValueError('SEPP retained + removed counts do not equal input')
        write_json(Path(state['directory']) / 'sepp-removed-counts.json', removed)


if __name__ == '__main__':
    main()
