"""Validate complete archive FASTQ sets before automatic layout detection."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile


def metadata(map_path, run):
    if map_path and Path(map_path).is_file():
        with open(map_path) as source:
            for row in csv.reader(source, delimiter='\t'):
                if row and row[0] == run:
                    return (row + [''] * 4)[:4]
    return [run, '', '', '']


def destination(name, run, prefix):
    if '/' in prefix or not name.startswith(run):
        raise ValueError(f'Unsupported archive filename or sample prefix: {name}, {prefix}')
    suffix = name[len(run):].replace('_subreads.fastq', '.fastq')
    if suffix not in ('.fastq', '_1.fastq', '_2.fastq', '.fastq.gz', '_1.fastq.gz', '_2.fastq.gz'):
        raise ValueError(f'Unsupported FASTQ filename: {name}')
    return prefix + suffix


def selected(rows, run, layout):
    # ENA may list an unpaired/orphan file alongside the two paired files.
    names = {Path(row[0]).name for row in rows}
    paired = any(n == run + '_1.fastq' or n == run + '_1.fastq.gz' for n in names) and any(
        n == run + '_2.fastq' or n == run + '_2.fastq.gz' for n in names)
    # LibraryLayout describes the submitted experiment, not necessarily the
    # released FASTQ representation (some paired experiments release only
    # merged/forward reads in a complete canonical singleton). Do not invent
    # a missing mate when the authoritative archive lists one canonical file.
    canonical = len(rows) == 1 and names <= {run + suffix for suffix in
        ('.fastq', '.fastq.gz', '_subreads.fastq', '_subreads.fastq.gz')}
    if layout.upper() == 'PAIRED' and not paired and not canonical:
        raise ValueError(f'{run}: archive declares PAIRED but the required R1/R2 files are incomplete')
    if paired:
        rows = [row for row in rows if Path(row[0]).name not in (run + '.fastq', run + '.fastq.gz')]
    if not rows:
        raise ValueError(f'{run}: no required FASTQ files found')
    if len({Path(row[0]).name for row in rows}) != len(rows):
        raise ValueError(f'{run}: duplicate FASTQ names in archive manifest')
    return rows


def ena_manifest(map_path, run, prefix):
    _, urls, checksums, layout = metadata(map_path, run)
    if not urls:
        raise ValueError(f'{run}: no ENA FASTQ URLs available')
    urls = urls.split(';')
    checksums = checksums.split(';') if checksums else [''] * len(urls)
    if len(checksums) != len(urls):
        raise ValueError(f'{run}: ENA URL/MD5 list lengths differ')
    rows = selected(list(zip(urls, checksums)), run, layout)
    return [(destination(Path(url).name, run, prefix), url, md5) for url, md5 in rows]


def signature(path):
    info = path.stat()
    return [info.st_size, info.st_mtime_ns, info.st_ctime_ns, info.st_ino, info.st_dev]


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as out:
            json.dump(value, out, sort_keys=True)
            out.write('\n')
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def receipt_path(path, receipts):
    return Path(receipts) / (hashlib.sha256(str(path.absolute()).encode()).hexdigest() + '.json')


def verify(path, md5='', receipts=None, require_gzip=False):
    path = Path(path)
    require_gzip = require_gzip or str(path).endswith('.gz')
    before = signature(path)
    if before[0] == 0:
        raise ValueError(f'{path}: empty file')
    expected = md5.lower()
    if expected and (len(expected) != 32 or any(c not in '0123456789abcdef' for c in expected)):
        raise ValueError(f'{path}: invalid expected MD5 value')
    record_path = receipt_path(path, receipts) if receipts else None
    if record_path and record_path.exists():
        record = json.loads(record_path.read_text())
        if record.get('schema') == 2 and record.get('signature') == before and (not expected or record.get('md5') == expected) and (not require_gzip or record.get('compressed')):
            return record
    actual = ''
    if expected:
        digest = hashlib.md5()
        with path.open('rb') as source:
            for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
                digest.update(block)
        actual = digest.hexdigest()
        if actual != expected:
            raise ValueError(f'{path}: MD5 mismatch: expected {expected}, found {actual}')
    # Read through EOF: gzip CRC/trailer errors must not be hidden by a pipe.
    with path.open('rb') as probe:
        compressed = probe.read(2) == b'\x1f\x8b'
    if require_gzip and not compressed:
        raise ValueError(f'{path}: expected gzip content but gzip header is missing')
    opener = gzip.open if compressed else open
    reads = 0
    previous_id, previous_mate = None, None
    interleaved_pairs = adjacent_duplicate_ids = 0
    with opener(path, 'rb') as source:
        while True:
            header = source.readline()
            if not header:
                break
            sequence = source.readline().rstrip(b'\r\n')
            plus = source.readline()
            quality_line = source.readline()
            quality = quality_line.rstrip(b'\r\n')
            if not header.startswith(b'@') or not plus.startswith(b'+') or not sequence or not quality_line or len(sequence) != len(quality):
                raise ValueError(f'{path}: invalid FASTQ record {reads + 1}')
            fields = header.split(None, 2)
            read_id, mate = fields[0], None
            if read_id.endswith((b'/1', b'/2')):
                read_id, mate = read_id[:-2], read_id[-1:]
            elif len(fields) > 1 and fields[1].startswith((b'1:', b'2:')):
                mate = fields[1][:1]
            if read_id == previous_id:
                adjacent_duplicate_ids += 1
                if previous_mate == b'1' and mate == b'2':
                    interleaved_pairs += 1
            previous_id, previous_mate = read_id, mate
            reads += 1
    if not reads:
        raise ValueError(f'{path}: no FASTQ records')
    if signature(path) != before:
        raise ValueError(f'{path}: file changed during integrity verification')
    record = dict(schema=2, path=str(path.absolute()), signature=before, reads=reads, md5=actual,
                  compressed=compressed, explicit_interleaved_pairs=interleaved_pairs,
                  adjacent_duplicate_read_ids=adjacent_duplicate_ids)
    if record_path:
        atomic_json(record_path, record)
    return record


def verify_set(rows, target, receipts):
    if not rows:
        raise ValueError('Empty required FASTQ manifest')
    records = {}
    for name, _, md5 in rows:
        if Path(name).name != name:
            raise ValueError(f'Invalid destination in FASTQ manifest: {name}')
        records[name] = verify(Path(target) / name, md5, receipts)
    if len(records) == 1:
        name, record = next(iter(records.items()))
        if record.get('explicit_interleaved_pairs'):
            raise ValueError(f'{name}: interleaved R1/R2 records detected in one FASTQ; '
                             'split into paired files before layout detection')
    for name, record in records.items():
        if '_1.fastq' in name:
            mate = name.replace('_1.fastq', '_2.fastq')
            if mate in records and records[mate]['reads'] != record['reads']:
                raise ValueError(f'Paired FASTQ read counts differ: {name}, {mate}')
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['ena-manifest', 'verify-ena', 'verify-list', 'verify-file', 'publish', 'ncbi-manifest'])
    parser.add_argument('--map')
    parser.add_argument('--run')
    parser.add_argument('--prefix')
    parser.add_argument('--target')
    parser.add_argument('--receipts')
    parser.add_argument('--file')
    parser.add_argument('--destination')
    parser.add_argument('--md5', default='')
    args = parser.parse_args()
    try:
        if args.action == 'ena-manifest':
            csv.writer(sys.stdout, delimiter='\t', lineterminator='\n').writerows(ena_manifest(args.map, args.run, args.prefix))
        elif args.action == 'verify-ena':
            verify_set(ena_manifest(args.map, args.run, args.prefix), args.target, args.receipts)
        elif args.action == 'verify-list':
            with open(args.file) as source:
                rows = [(row + ['', ''])[:3] for row in csv.reader(source, delimiter='\t') if row]
            verify_set(rows, args.target, args.receipts)
        elif args.action in ('verify-file', 'publish'):
            record = verify(args.file, args.md5, args.receipts, bool(args.destination and args.destination.endswith('.gz')))
            if args.action == 'publish':
                destination_path = Path(args.destination)
                os.replace(args.file, destination_path)
                record.update(path=str(destination_path.absolute()), signature=signature(destination_path))
                if args.receipts:
                    atomic_json(receipt_path(destination_path, args.receipts), record)
        else:
            archive = metadata(args.map, args.run)
            rows = selected([(str(p), '') for p in sorted(Path(args.target).glob(args.run + '*.fastq'))],
                            args.run, archive[3])
            if len(rows) == 1 and archive[1]:
                archive_names = {Path(url).name for url in archive[1].split(';')}
                if {args.run + '_1.fastq.gz', args.run + '_2.fastq.gz'} <= archive_names:
                    raise ValueError(f'{args.run}: ENA explicitly lists R1/R2 but NCBI extraction produced one file')
            counts = {}
            for raw, _ in rows:
                record = verify(raw)
                if len(rows) == 1 and record.get('explicit_interleaved_pairs'):
                    raise ValueError(f'{raw}: interleaved R1/R2 records detected in singleton extraction')
                counts[Path(raw).name] = record['reads']
            left, right = args.run + '_1.fastq', args.run + '_2.fastq'
            if left in counts and right in counts and counts[left] != counts[right]:
                raise ValueError(f'{args.run}: fasterq-dump produced unequal paired read counts')
            csv.writer(sys.stdout, delimiter='\t', lineterminator='\n').writerows(
                (raw, destination(Path(raw).name, args.run, args.prefix) + '.gz') for raw, _ in rows)
    except (OSError, EOFError, ValueError) as error:
        parser.exit(1, f'FASTQ integrity error: {error}\n')


if __name__ == '__main__':
    main()
