#!/usr/bin/env python3
"""Trace 454 cluster members and prevent identified chimera reads from remapping.

All membership is taken from the actual 99% and UNOISE UC files. Reads omitted
by UNOISE's abundance threshold remain eligible for mapping unless their
precluster belongs to a representative explicitly classified as a chimera.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile

from read_layout import discover, layout, open_fastq

SIZE = re.compile(r';size=(\d+)(?=;|$)')


def label(value):
    parts = value.split()
    return parts[0].split(';', 1)[0] if parts else ''


def fasta(path, allow_empty=False):
    records = {}
    current, chunks = None, []

    def add():
        if current is None:
            return
        identifier = label(current)
        matches = SIZE.findall(current.split()[0])
        sequence = ''.join(chunks).upper().replace('U', 'T')
        if not identifier or identifier in records or len(matches) != 1:
            raise ValueError(f'Invalid, repeated or unannotated FASTA label in {path}: {current}')
        abundance = int(matches[0])
        if abundance < 1 or not sequence or any(base.isspace() for base in sequence):
            raise ValueError(f'Invalid FASTA sequence/abundance in {path}: {current}')
        records[identifier] = dict(sequence=sequence, abundance=abundance)

    with Path(path).open() as stream:
        for line in stream:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                add()
                current, chunks = line[1:], []
            elif line:
                if current is None:
                    raise ValueError(f'Sequence before FASTA header: {path}')
                chunks.append(line)
    add()
    if not records and not allow_empty:
        raise ValueError(f'No representative sequences in {path}')
    return records


def fastq_records(path):
    with open_fastq(path) as stream:
        while True:
            header = stream.readline()
            if not header:
                return
            seq, plus, qual = [stream.readline().rstrip('\r\n') for _ in range(3)]
            header = header.rstrip('\r\n')
            if (not header.startswith('@') or not plus.startswith('+') or not seq
                    or len(seq) != len(qual) or any(base.isspace() for base in seq)):
                raise ValueError(f'Malformed FASTQ record in {path}: {header}')
            yield header, seq, plus, qual


def single_samples(directory):
    rows = discover(directory)
    if layout(rows) != 'SE' or len({row['sample'] for row in rows}) != len(rows):
        raise ValueError('454 member tracking requires one staged single-end file per sample')
    return rows


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')
    temporary.replace(path)


def to_fasta(input_dir, output):
    """Give every QC read a new label, ignoring any original size annotation."""
    rows = single_samples(input_dir)
    output = Path(output)
    sources = {Path(row['r1']).resolve() for row in rows}
    if output.resolve() in sources:
        raise ValueError('Combined FASTA must not overwrite a QC input')
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        total = 0
        with tempfile.NamedTemporaryFile('w', dir=output.parent, delete=False) as stream:
            temporary = Path(stream.name)
            for row in rows:
                for _, sequence, _, _ in fastq_records(row['r1']):
                    total += 1
                    stream.write(f'>LS454_read_{total}\n{sequence}\n')
        if not total:
            raise ValueError('No QC reads remain for 454 dereplication')
        temporary.replace(output)
        return total
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def uc_members(path, sources, centroids, minsize=None):
    """Return source-to-centroid assignments, validating every eligible source."""
    assignments = {}
    with Path(path).open() as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 10 or fields[0] not in ('S', 'H', 'N', 'C'):
                raise ValueError(f'Invalid UC row {number} in {path}')
            kind, query = fields[0], label(fields[8])
            if kind == 'C':
                if query not in centroids:
                    raise ValueError(f'Unknown UC cluster summary {query} in {path}')
                continue
            if query not in sources or query in assignments:
                raise ValueError(f'Unknown or repeated UC query {query} in {path}')
            if kind == 'N':
                assignments[query] = None
                continue
            target = query if kind == 'S' else label(fields[9])
            if target not in centroids:
                raise ValueError(f'Unknown UC centroid {target} in {path}')
            if kind == 'H' and fields[4] not in ('+', '*'):
                raise ValueError(f'Unexpected reverse-strand member in 454 plus-strand UC: {query}')
            assignments[query] = target
    for query, record in sources.items():
        if query not in assignments:
            assignments[query] = None
        if assignments[query] is None:
            if minsize is None or record['abundance'] >= minsize:
                raise ValueError(f'UC is missing eligible sequence {query} in {path}')
    totals = Counter()
    for query, target in assignments.items():
        if target is not None:
            totals[target] += sources[query]['abundance']
    for target, centroid in centroids.items():
        if assignments.get(target) != target or sources[target]['sequence'] != centroid['sequence']:
            raise ValueError(f'UC centroid identity is inconsistent for {target} in {path}')
        if totals[target] != centroid['abundance']:
            raise ValueError(f'UC member abundance does not conserve centroid {target} in {path}: '
                             f"members={totals[target]}, centroid={centroid['abundance']}")
    return assignments


def trace_members(dereplicated, preclustered, precluster_uc, denoised, denoise_uc,
                  nonchimeras, chimeras, borderline, minsize):
    if minsize < 1:
        raise ValueError('UNOISE minsize must be positive')
    derep = fasta(dereplicated)
    if len({record['sequence'] for record in derep.values()}) != len(derep):
        raise ValueError('Dereplicated FASTA contains a repeated sequence')
    pre = fasta(preclustered)
    denoise = fasta(denoised)
    first = uc_members(precluster_uc, derep, pre)
    second = uc_members(denoise_uc, pre, denoise, minsize)
    classified = {}
    for category, filename in [('nonchimera', nonchimeras), ('chimera', chimeras), ('borderline', borderline)]:
        for identifier, record in fasta(filename, allow_empty=True).items():
            if identifier in classified or identifier not in denoise or record != denoise[identifier]:
                raise ValueError(f'Inconsistent chimera classification for {identifier}: {filename}')
            classified[identifier] = category
    if set(classified) != set(denoise):
        raise ValueError('Chimera output files do not account for every denoised representative')
    members = []
    for identifier, record in derep.items():
        precluster = first[identifier]
        representative = second[precluster]
        category = classified[representative] if representative is not None else 'below_unoise_minsize'
        members.append(dict(derep_id=identifier, precluster_id=precluster,
                            denoised_id=representative or '', classification=category,
                            excluded=category == 'chimera', abundance=record['abundance'],
                            sequence_sha256=hashlib.sha256(record['sequence'].encode()).hexdigest()))
    return derep, members


def validate_output_paths(args):
    """Validate write destinations before even writing a failure report."""
    source_dir = Path(args.input).resolve()
    destination = Path(args.output_dir).resolve()
    evidence = {Path(getattr(args, name)).resolve() for name in
                ('dereplicated', 'preclustered', 'precluster_uc', 'denoised',
                 'denoise_uc', 'nonchimeras', 'chimeras', 'borderline')}
    outputs = [Path(getattr(args, name)).resolve() for name in
               ('report', 'members_tsv', 'excluded_fasta')]
    if len(set(outputs)) != len(outputs):
        raise ValueError('Report and member evidence outputs must be distinct')
    if (destination == source_dir or destination in source_dir.parents
            or source_dir in destination.parents):
        raise ValueError('Mapping output must be separate from the QC input directory')
    for path in outputs:
        if (path == source_dir or source_dir in path.parents or path in source_dir.parents
                or path in evidence or path == destination or destination in path.parents
                or path in destination.parents):
            raise ValueError(f'Output must not overwrite QC, input evidence or mapping reads: {path}')
        if any(path in item.parents or item in path.parents for item in evidence):
            raise ValueError(f'Output conflicts with an input evidence path: {path}')
    if any(destination == item or destination in item.parents or item in destination.parents
           for item in evidence):
        raise ValueError('Mapping output must be separate from input evidence')
    for path in outputs:
        if any(path != other and (path in other.parents or other in path.parents) for other in outputs):
            raise ValueError('Report and member evidence outputs must not contain one another')
        if path.exists() and not path.is_file():
            raise ValueError(f'Output must be a file: {path}')


def filter_members(args):
    validate_output_paths(args)
    report = dict(schema_version=1, status='failed', mode='LS454/vsearch',
                  exclusion_policy='Exclude only members of representatives explicitly classified as chimeric; '
                                   'retain borderline and below-minsize reads for mapping.',
                  input=str(Path(args.input).resolve()), output=str(Path(args.output_dir).resolve()),
                  minsize=args.minsize, samples={})
    staging = None
    try:
        rows = single_samples(args.input)
        destination = Path(args.output_dir).resolve()
        input_dir = Path(args.input).resolve()
        if (destination == input_dir or destination in input_dir.parents or input_dir in destination.parents):
            raise ValueError('Mapping output must be separate from the QC input directory')
        if destination.exists() and any(destination.iterdir()):
            raise ValueError('Mapping output directory must be absent or empty')
        derep, members = trace_members(args.dereplicated, args.preclustered, args.precluster_uc,
            args.denoised, args.denoise_uc, args.nonchimeras, args.chimeras, args.borderline, args.minsize)
        by_sequence = {record['sequence']: identifier for identifier, record in derep.items()}
        excluded = {row['derep_id'] for row in members if row['excluded']}
        observed = Counter()
        destination.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix='.ls454-mapping-', dir=destination.parent))
        for row in rows:
            target = staging / Path(row['r1']).name
            count, removed = 0, 0
            with open_fastq(target, 'wt') as out:
                for header, sequence, plus, quality in fastq_records(row['r1']):
                    identifier = by_sequence.get(sequence.upper().replace('U', 'T'))
                    if identifier is None:
                        raise ValueError(f'QC FASTQ contains sequence absent from dereplication: {row["sample"]}')
                    observed[identifier] += 1
                    count += 1
                    if identifier in excluded:
                        removed += 1
                    else:
                        out.write(f'{header}\n{sequence}\n{plus}\n{quality}\n')
            report['samples'][row['sample']] = dict(input=count, excluded_chimera=removed, kept=count - removed)
        expected = {identifier: record['abundance'] for identifier, record in derep.items()}
        if dict(observed) != expected:
            raise ValueError('QC read counts do not conserve dereplicated abundances')
        report['totals'] = {key: sum(row[key] for row in report['samples'].values())
                            for key in ('input', 'excluded_chimera', 'kept')}
        assert report['totals']['input'] == report['totals']['excluded_chimera'] + report['totals']['kept']
        report['members'] = dict(dereplicated=len(members), excluded_chimera=len(excluded),
            below_unoise_minsize=sum(row['classification'] == 'below_unoise_minsize' for row in members),
            below_unoise_minsize_reads=sum(row['abundance'] for row in members
                                          if row['classification'] == 'below_unoise_minsize'),
            borderline=sum(row['classification'] == 'borderline' for row in members))
        for filename in (args.members_tsv, args.excluded_fasta):
            Path(filename).parent.mkdir(parents=True, exist_ok=True)
        with Path(args.members_tsv).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(members[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(members)
        with Path(args.excluded_fasta).open('w') as stream:
            for identifier in sorted(excluded):
                record = derep[identifier]
                stream.write(f'>{identifier};size={record["abundance"]};\n{record["sequence"]}\n')
        report['evidence'] = {name: str(Path(getattr(args, name)).resolve()) for name in
            ('precluster_uc', 'denoise_uc', 'chimeras', 'borderline', 'nonchimeras', 'members_tsv', 'excluded_fasta')}
        if report['totals']['kept'] == 0:
            raise ValueError('All QC reads belong to identified chimeras; no reads remain for mapping')
        if destination.exists():
            destination.rmdir()
        staging.replace(destination)
        staging = None
        report['status'] = 'completed'
        save_json(args.report, report)
        return report
    except (OSError, ValueError) as error:
        report.update(status='failed', error=str(error))
        save_json(args.report, report)
        raise
    finally:
        if staging is not None:
            shutil.rmtree(staging)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    combined = commands.add_parser('to-fasta', help='Convert QC reads with unique, unit-abundance labels')
    combined.add_argument('--input', required=True)
    combined.add_argument('--output', required=True)
    filtering = commands.add_parser('filter', help='Remove traced chimera members before mapping')
    for name in ('input', 'output-dir', 'dereplicated', 'preclustered', 'precluster-uc',
                 'denoised', 'denoise-uc', 'nonchimeras', 'chimeras', 'borderline',
                 'report', 'members-tsv', 'excluded-fasta'):
        filtering.add_argument('--' + name, required=True)
    filtering.add_argument('--minsize', type=int, required=True)
    args = parser.parse_args()
    try:
        if args.command == 'to-fasta':
            print(f'LS454_QC_READS={to_fasta(args.input, args.output)}')
        else:
            report = filter_members(args)
            print('LS454_MAPPING_READS=' + str(report['totals']['kept']))
            print('LS454_EXCLUDED_CHIMERA_READS=' + str(report['totals']['excluded_chimera']))
    except (OSError, ValueError) as error:
        parser.exit(2, f'454 member tracking failed: {error}\n')


if __name__ == '__main__':
    main()
