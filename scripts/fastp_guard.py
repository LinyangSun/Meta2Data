#!/usr/bin/env python3
"""Check fastp's de novo adapter candidates against an immutable 16S reference.

A strong, single-HSP 16S match requests a sample-local fastp retry. The absence
of a match is not evidence that an inferred sequence is a genuine adapter.
This module records decisions and never modifies FASTQ files.
"""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import tempfile
import zipfile

SCHEMA_VERSION = 1
OUTFMT = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
MAX_EVALUE = 1e-5


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def object_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


@contextmanager
def locked(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def blast_version():
    output = subprocess.run(['blastn', '-version'], check=True, text=True,
                            capture_output=True).stdout.strip()
    if not output.startswith('blastn:'):
        raise ValueError('Could not identify the blastn version')
    return output.splitlines()[0]


def copy_reference(source, target):
    """Copy one validated FASTA from a plain/gzip FASTA or Sequence QZA.

    No archive member is extracted by pathname. All names are checked anyway
    so an unsafe or ambiguous QZA is rejected before its sequence data is used.
    """
    source, target = Path(source), Path(target)
    if zipfile.is_zipfile(source):
        with zipfile.ZipFile(source) as archive:
            names = []
            for member in archive.infolist():
                name = PurePosixPath(member.filename)
                if (name.is_absolute() or '..' in name.parts or '\\' in member.filename
                        or stat.S_ISLNK(member.external_attr >> 16)):
                    raise ValueError('Unsafe QZA member path or symbolic link')
                names.append(member.filename)
            if len(names) != len(set(names)):
                raise ValueError('Duplicate QZA archive members')
            metadata = [name for name in names if len(PurePosixPath(name).parts) == 2
                        and name.endswith('/metadata.yaml')]
            if len(metadata) != 1:
                raise ValueError('Expected exactly one QZA root metadata file')
            metadata_text = archive.read(metadata[0]).decode('utf-8')
            artifact_type = re.search(r'^type:\s*[\'"]?([^\r\n\'"]+)', metadata_text, re.MULTILINE)
            if not artifact_type or artifact_type.group(1).strip() != 'FeatureData[Sequence]':
                raise ValueError('Reference QZA must have type FeatureData[Sequence]')
            root = PurePosixPath(metadata[0]).parts[0]
            candidates = [name for name in names if name.startswith(root + '/data/')
                          and PurePosixPath(name).suffix.lower() in ('.fasta', '.fa', '.fna')]
            if len(candidates) != 1:
                raise ValueError('Expected exactly one nucleotide FASTA in the reference QZA')
            with archive.open(candidates[0]) as reader, target.open('wb') as writer:
                shutil.copyfileobj(reader, writer, 1024 * 1024)
    else:
        if source.suffix.lower() == '.qza':
            raise ValueError('Reference QZA is not a valid ZIP archive')
        opener = gzip.open if source.suffix.lower() == '.gz' else open
        with opener(source, 'rb') as reader, target.open('wb') as writer:
            shutil.copyfileobj(reader, writer, 1024 * 1024)
    sequences = bases = current_bases = 0
    with target.open('rt', encoding='ascii') as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if sequences and not current_bases:
                    raise ValueError('Reference FASTA has an empty sequence')
                if not line[1:].strip():
                    raise ValueError('Reference FASTA has an empty identifier')
                sequences += 1
                current_bases = 0
            else:
                if not sequences or re.fullmatch('[ACGTRYSWKMBDHVNUacgtryswkmbdhvnu]+', line) is None:
                    raise ValueError('Reference FASTA contains invalid nucleotide sequence data')
                current_bases += len(line)
                bases += len(line)
    if not sequences or not current_bases:
        raise ValueError('Reference FASTA is empty or ends with an empty sequence')
    return {'sequences': sequences, 'bases': bases}


def file_record(path):
    path = Path(path)
    return {'name': path.name, 'bytes': path.stat().st_size, 'sha256': sha256(path)}


def load_manifest(path, verify_hashes=False):
    path = Path(path).resolve()
    manifest = json.loads(path.read_text())
    if not isinstance(manifest, dict) or manifest.get('schema_version') != SCHEMA_VERSION:
        raise ValueError('Unsupported adapter reference manifest')
    if not re.fullmatch('[0-9a-f]{64}', str(manifest.get('reference_sha256', ''))):
        raise ValueError('Missing reference SHA256 in adapter manifest')
    if not isinstance(manifest.get('blast_version'), str) or not manifest['blast_version'].startswith('blastn:'):
        raise ValueError('Missing BLAST version in adapter manifest')
    prefix = Path(manifest.get('database_prefix', ''))
    if not prefix.is_absolute() or prefix.parent != path.parent:
        raise ValueError('BLAST database prefix must be beside its manifest')
    files = manifest.get('database_files')
    if not isinstance(files, list) or not files:
        raise ValueError('Adapter reference manifest has no database files')
    names = set()
    for item in files:
        if not isinstance(item, dict):
            raise ValueError('Invalid database file record')
        name = item.get('name', '')
        if not isinstance(name, str) or Path(name).name != name or not name.startswith(prefix.name + '.') or name in names:
            raise ValueError('Invalid or duplicate database file name')
        names.add(name)
        index = prefix.parent / name
        if not index.is_file() or index.stat().st_size != item.get('bytes'):
            raise ValueError(f'Missing or incomplete BLAST database file: {index}')
        if verify_hashes and sha256(index) != item.get('sha256'):
            raise ValueError(f'BLAST database checksum mismatch: {index}')
    if not all(any(name.endswith(suffix) for name in names) for suffix in ('.nhr', '.nin', '.nsq')):
        raise ValueError('Incomplete nucleotide BLAST database index')
    return manifest


def prepare(reference, cache_dir):
    reference = Path(reference).resolve(strict=True)
    cache = Path(cache_dir).resolve()
    staging_root = cache / 'staging'
    staging_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='reference-', dir=staging_root))
    try:
        fasta = staging / 'reference.fasta'
        reference_stats = copy_reference(reference, fasta)
        reference_sha = sha256(fasta)
        version = blast_version()
        key = object_hash({'schema_version': SCHEMA_VERSION, 'reference_sha256': reference_sha,
                           'blast_version': version})
        database_dir = cache / 'databases' / key
        manifest_path = database_dir / 'manifest.json'
        with locked(cache / 'locks' / ('database-' + key + '.lock')):
            if database_dir.exists():
                cached = load_manifest(manifest_path, verify_hashes=True)
                if cached['reference_sha256'] != reference_sha or cached['blast_version'] != version:
                    raise ValueError('Cached reference manifest identity mismatch')
                if sha256(database_dir / 'reference.fasta') != reference_sha:
                    raise ValueError('Cached reference FASTA checksum mismatch')
                return manifest_path
            prefix = staging / 'reference'
            command = ['makeblastdb', '-in', str(fasta), '-dbtype', 'nucl',
                       '-parse_seqids', '-out', str(prefix)]
            with (staging / 'makeblastdb.log').open('w') as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            index_files = [p for p in sorted(staging.glob('reference.*')) if p != fasta]
            # BLAST v5 may produce extra index files; list all of them permanently.
            if not all(any(p.name.endswith(suffix) for p in index_files)
                       for suffix in ('.nhr', '.nin', '.nsq')):
                raise ValueError('makeblastdb did not produce a complete nucleotide index')
            manifest = dict(schema_version=SCHEMA_VERSION, reference_sha256=reference_sha,
                            source_path=str(reference), source_sha256=sha256(reference),
                            reference_fasta=str(database_dir / 'reference.fasta'),
                            blast_version=version, database_prefix=str(database_dir / 'reference'),
                            database_files=[file_record(p) for p in index_files],
                            reference_stats=reference_stats,
                            created_at=datetime.now(timezone.utc).isoformat(),
                            build_command=command)
            save_json(staging / 'manifest.json', manifest)
            database_dir.parent.mkdir(parents=True, exist_ok=True)
            staging.replace(database_dir)
            load_manifest(manifest_path, verify_hashes=True)
        return manifest_path
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def candidate(report, end):
    section = report.get('adapter_cutting', {})
    if not isinstance(section, dict):
        raise ValueError('fastp adapter_cutting must be an object')
    name = f'{end}_adapter_sequence'
    if name not in section:
        return None
    value = section[name]
    if not isinstance(value, str):
        raise ValueError(f'fastp {name} must be a string')
    value = value.strip()
    if not value or value.lower() == 'unspecified':
        return None
    if not re.fullmatch('[ACGTacgt]+', value):
        raise ValueError(f'fastp {name} contains invalid nucleotide symbols')
    return value.upper()


def parameters(min_length=50, min_identity=98, min_coverage=95):
    if (not isinstance(min_length, int) or isinstance(min_length, bool) or min_length < 1
            or not math.isfinite(min_identity) or not math.isfinite(min_coverage)
            or not 0 <= min_identity <= 100 or not 0 <= min_coverage <= 100):
        raise ValueError('Invalid adapter guard thresholds')
    return dict(min_length=min_length, min_identity=float(min_identity), min_coverage=float(min_coverage),
                max_evalue=MAX_EVALUE, task='blastn-short', strand='both', dust='no',
                soft_masking=False, word_size=7, max_target_seqs=100, max_hsps=5,
                coverage_rule='single_hsp_query_span', outfmt=OUTFMT)


def parse_hits(path, query_length, params):
    hits = []
    with Path(path).open() as stream:
        for line in stream:
            if not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 14 or fields[0] != 'candidate':
                raise ValueError('Invalid adapter BLAST result row')
            qstart, qend, qlen = int(fields[6]), int(fields[7]), int(fields[12])
            pident, evalue, bitscore = float(fields[2]), float(fields[10]), float(fields[11])
            if (qlen != query_length or not 1 <= qstart <= qlen or not 1 <= qend <= qlen
                    or not all(math.isfinite(v) for v in (pident, evalue, bitscore))
                    or not 0 <= pident <= 100 or evalue < 0):
                raise ValueError('Inconsistent adapter BLAST query length or coordinates')
            coverage = 100 * (abs(qend - qstart) + 1) / qlen
            strong = (query_length >= params['min_length'] and pident >= params['min_identity']
                      and coverage >= params['min_coverage'] and evalue <= params['max_evalue'])
            hits.append(dict(subject_id=fields[1], identity_percent=pident, alignment_length=int(fields[3]),
                             mismatches=int(fields[4]), gap_opens=int(fields[5]), query_start=qstart,
                             query_end=qend, subject_start=int(fields[8]), subject_end=int(fields[9]),
                             evalue=evalue, bitscore=bitscore, query_length=qlen, subject_length=int(fields[13]),
                             query_coverage_percent=coverage, qualifies=strong))
    hits.sort(key=lambda hit: (hit['qualifies'], hit['bitscore'], hit['query_coverage_percent']), reverse=True)
    return hits


def query_reference(sequence, manifest, cache_dir, params, threads):
    """Publish a complete result once, then share it across identical candidates."""
    cache = Path(cache_dir).resolve()
    identity = dict(schema_version=SCHEMA_VERSION, sequence=sequence,
                    reference_sha256=manifest['reference_sha256'], blast_version=manifest['blast_version'],
                    parameters=params)
    key = object_hash(identity)
    folder = cache / 'queries' / key
    with locked(cache / 'locks' / ('query-' + key + '.lock')):
        if folder.exists():
            result = json.loads((folder / 'result.json').read_text())
            if not isinstance(result, dict) or result.get('identity') != identity or result.get('cache_key') != key:
                raise ValueError('Adapter query cache identity mismatch')
            for artifact in result.get('artifacts', []):
                filename = artifact.get('name', '')
                if Path(filename).name != filename or sha256(folder / filename) != artifact.get('sha256'):
                    raise ValueError('Adapter query cache checksum mismatch')
            if {artifact['name'] for artifact in result.get('artifacts', [])} != {'query.fasta', 'blast.tsv', 'command.json', 'blast.stderr.log'}:
                raise ValueError('Adapter query cache has incomplete evidence')
            hits = parse_hits(folder / 'blast.tsv', len(sequence), params)
            expected = dict(status='protected_16s' if any(hit['qualifies'] for hit in hits) else 'no_strong_hit',
                            hsp_count=len(hits), qualifying_hsp_count=sum(hit['qualifies'] for hit in hits),
                            best_hit=hits[0] if hits else None)
            if any(result.get(name) != value for name, value in expected.items()):
                raise ValueError('Adapter query cache decision disagrees with its BLAST evidence')
            return result, folder, True
        (cache / 'staging').mkdir(parents=True, exist_ok=True)
        stage = Path(tempfile.mkdtemp(prefix='query-', dir=cache / 'staging'))
        try:
            query = stage / 'query.fasta'
            query.write_text('>candidate\n' + sequence + '\n')
            output = stage / 'blast.tsv'
            command = ['blastn', '-task', params['task'], '-query', str(query),
                       '-db', manifest['database_prefix'], '-out', str(output), '-outfmt', OUTFMT,
                       '-evalue', str(params['max_evalue']), '-word_size', str(params['word_size']),
                       '-dust', 'no', '-soft_masking', 'false', '-strand', 'both',
                       '-perc_identity', str(params['min_identity']),
                       '-max_target_seqs', str(params['max_target_seqs']), '-max_hsps', str(params['max_hsps']),
                       '-num_threads', str(threads)]
            save_json(stage / 'command.json', command)
            with (stage / 'blast.stderr.log').open('w') as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            hits = parse_hits(output, len(sequence), params)
            result = dict(identity=identity, cache_key=key, status='protected_16s' if any(h['qualifies'] for h in hits) else 'no_strong_hit',
                          hsp_count=len(hits), qualifying_hsp_count=sum(h['qualifies'] for h in hits),
                          best_hit=hits[0] if hits else None,
                          created_at=datetime.now(timezone.utc).isoformat(),
                          artifacts=[file_record(stage / name) for name in
                                     ('query.fasta', 'blast.tsv', 'command.json', 'blast.stderr.log')])
            save_json(stage / 'result.json', result)
            folder.parent.mkdir(parents=True, exist_ok=True)
            stage.replace(folder)
            return result, folder, False
        finally:
            if stage.exists():
                shutil.rmtree(stage)


def check(report_path, db_manifest, cache_dir, output_dir, sample_id, threads=1,
          min_length=50, min_identity=98, min_coverage=95, auto_r1=False, auto_r2=False):
    if not isinstance(threads, int) or threads < 1:
        raise ValueError('BLAST thread count must be positive')
    if not sample_id or '/' in sample_id or sample_id in ('.', '..'):
        raise ValueError('Invalid sample ID')
    params = parameters(min_length, min_identity, min_coverage)
    manifest = load_manifest(db_manifest)
    if blast_version() != manifest['blast_version']:
        raise ValueError('BLAST version differs from the prepared reference manifest; prepare again')
    report_path = Path(report_path).resolve(strict=True)
    report = json.loads(report_path.read_text())
    if not isinstance(report, dict):
        raise ValueError('fastp report must be a JSON object')
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    per_end = {}
    for end, enabled in (('read1', auto_r1), ('read2', auto_r2)):
        evidence = {'auto_inference': bool(enabled), 'status': 'not_checked', 'sequence': None}
        per_end[end] = evidence
        if not enabled:
            evidence['reason'] = 'de_novo_inference_not_enabled_for_this_end'
            continue
        sequence = candidate(report, end)
        evidence.update(sequence=sequence, length=len(sequence) if sequence else 0)
        if sequence is None:
            evidence.update(status='no_candidate', reason='no_inferred_adapter_sequence_in_report')
        elif len(sequence) < min_length:
            evidence.update(status='uncertain', reason='candidate_shorter_than_min_length')
        else:
            result, folder, cache_hit = query_reference(sequence, manifest, cache_dir, params, threads)
            sample_evidence = out / end
            sample_evidence.mkdir(exist_ok=True)
            for name in ('query.fasta', 'blast.tsv', 'command.json', 'blast.stderr.log', 'result.json'):
                shutil.copyfile(folder / name, sample_evidence / name)
            evidence.update(status=result['status'], cache_hit=cache_hit, cache_key=result['cache_key'],
                            cache_directory=str(folder), evidence_directory=str(sample_evidence),
                            hsp_count=result['hsp_count'], qualifying_hsp_count=result['qualifying_hsp_count'],
                            best_hit=result['best_hit'])
    states = {end['status'] for end in per_end.values()}
    status = next((value for value in ('protected_16s', 'uncertain', 'no_strong_hit') if value in states), 'no_candidate')
    decision = dict(schema_version=SCHEMA_VERSION, sample_id=sample_id, status=status,
                    action='fallback' if status == 'protected_16s' else 'accept', per_end=per_end,
                    report_path=str(report_path), report_sha256=sha256(report_path),
                    reference_manifest=str(Path(db_manifest).resolve()),
                    reference_sha256=manifest['reference_sha256'], blast_version=manifest['blast_version'],
                    parameters=params, created_at=datetime.now(timezone.utc).isoformat(),
                    interpretation='A strong 16S match triggers protection; a missing match does not establish adapter identity.')
    decision_path = out / 'decision.json'
    save_json(decision_path, decision)
    return decision_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prepare_parser = commands.add_parser('prepare', help='Validate the reference and publish a cached BLAST database')
    prepare_parser.add_argument('--reference', required=True)
    prepare_parser.add_argument('--cache-dir', required=True)
    check_parser = commands.add_parser('check', help='Evaluate one sample report without changing FASTQs')
    for name in ('report', 'db-manifest', 'cache-dir', 'output-dir', 'sample-id'):
        check_parser.add_argument('--' + name, required=True)
    check_parser.add_argument('--threads', type=int, default=1)
    check_parser.add_argument('--min-length', type=int, default=50)
    check_parser.add_argument('--min-identity', type=float, default=98)
    check_parser.add_argument('--min-coverage', type=float, default=95)
    check_parser.add_argument('--auto-r1', action='store_true')
    check_parser.add_argument('--auto-r2', action='store_true')
    args = parser.parse_args()
    try:
        if args.command == 'prepare':
            result = prepare(args.reference, args.cache_dir)
        else:
            result = check(args.report, args.db_manifest, args.cache_dir, args.output_dir, args.sample_id,
                           args.threads, args.min_length, args.min_identity, args.min_coverage,
                           args.auto_r1, args.auto_r2)
        print(result)
    except (OSError, ValueError, KeyError, UnicodeError, zipfile.BadZipFile, subprocess.CalledProcessError) as error:
        parser.exit(2, f'fastp adapter guard failed: {error}\n')


if __name__ == '__main__':
    main()
