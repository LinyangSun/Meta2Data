"""Invalidate processing checkpoints without discarding reusable source reads."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import uuid
from parameters import current, nested
from read_layout import files, validate_local_source

SCHEMA_VERSION = 2
RAW_PATHS = ('downloaded_fastq', 'ori_fastq', 'download_integrity', 'download_source.tsv')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def atomic_json(path, value):
    """Publish complete JSON so interrupted writes cannot create a checkpoint."""
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as out:
            json.dump(value, out, indent=2)
            out.write('\n')
            out.flush()
            os.fsync(out.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def read_state(path):
    if not path.exists():
        return {}, 'missing'
    try:
        value = json.loads(path.read_text())
        if not isinstance(value, dict):
            raise ValueError('Checkpoint is not an object')
        return value, 'valid'
    except (ValueError, OSError):
        return {}, 'invalid'


def input_identity(state):
    """Derive the same input identity from current and pre-schema checkpoints."""
    if not isinstance(state.get('inputs'), str) or not isinstance(state.get('local_source'), list):
        return None
    source = state['local_source']
    if any(not isinstance(row, list) or len(row) != 3 or
           not isinstance(row[0], str) or type(row[1]) is not int or
           type(row[2]) is not int for row in source):
        return None
    kind = state.get('source_kind', 'local' if source else 'archive')
    if kind not in ('local', 'archive'):
        return None
    identity = dict(inputs=state['inputs'], local_source=source, source_kind=kind)
    if state.get('schema_version') == SCHEMA_VERSION and state.get('input_fingerprint') != digest(identity):
        return None
    return identity


def code_fingerprint():
    scripts = Path(__file__).parent
    code = hashlib.sha256()
    paths = list(scripts.glob('*.py')) + list(scripts.glob('*.sh'))
    paths += list((scripts.parent / 'docs').glob('*.fas'))
    for path in sorted(paths):
        code.update(path.name.encode())
        code.update(path.read_bytes())
    return code.hexdigest()


def guard_identity(settings):
    """Fingerprint biological decisions, excluding cache locations/timestamps."""
    enabled = os.environ.get('M2D_ADAPTER_GUARD_ENABLED', '0')
    if enabled not in ('0', '1'):
        raise ValueError('M2D_ADAPTER_GUARD_ENABLED must be 0 or 1')
    if enabled == '0':
        return {'enabled': False}
    location = os.environ.get('M2D_ADAPTER_GUARD_DB_MANIFEST', '')
    if not location:
        raise ValueError('Enabled adapter guard requires M2D_ADAPTER_GUARD_DB_MANIFEST')
    manifest, status = read_state(Path(location))
    if status != 'valid':
        raise ValueError(f'Adapter guard manifest is missing or invalid: {location}')
    reference = manifest.get('reference_sha256')
    version = manifest.get('blast_version')
    schema = manifest.get('schema_version')
    prefix = manifest.get('database_prefix')
    if (not isinstance(reference, str) or re.fullmatch(r'[0-9a-fA-F]{64}', reference) is None or
            not isinstance(version, str) or not version.strip() or
            type(schema) is not int or schema < 1 or
            not isinstance(prefix, str) or not prefix.strip()):
        raise ValueError(f'Adapter guard manifest has invalid identity fields: {location}')
    return dict(enabled=True, reference_sha256=reference.lower(), blast_version=version,
                schema_version=schema,
                min_length=settings['adapter_guard.min_length'],
                min_identity=settings['adapter_guard.min_identity'],
                min_coverage=settings['adapter_guard.min_coverage'])


def remove_owned(path):
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.exists():
        shutil.rmtree(path)


def archive_raw(dataset, previous, old_identity, new_identity, reason):
    present = [dataset / name for name in RAW_PATHS
               if (dataset / name).exists() or (dataset / name).is_symlink()]
    if not present:
        return
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    previous_hash = digest(old_identity)[:16] if old_identity is not None else 'unknown-input'
    archive = dataset / 'raw_cache_archive' / f'{previous_hash}-{stamp}-{uuid.uuid4().hex[:8]}'
    archive.mkdir(parents=True)
    atomic_json(archive / 'archive.json', dict(reason=reason, previous_state=previous,
                previous_input_identity=old_identity, next_input_identity=new_identity,
                archived_paths=[path.name for path in present]))
    for path in present:
        path.rename(archive / path.name)


def prepare(dataset, method, local_source='', platform='', forward='', reverse=''):
    dataset = Path(dataset)
    name = dataset.name
    if local_source:
        validate_local_source(local_source, dataset)
    settings = {k: v for k, v in current().items() if not k.startswith('taxa.')}
    inputs = (dataset / f'{name}_sra.txt').read_text()
    source = []
    if local_source:
        for path in files(local_source):
            stat = path.stat()
            source.append([str(path.resolve()), stat.st_size, stat.st_mtime_ns])
    identity = dict(inputs=inputs, local_source=source,
                    source_kind='local' if local_source else 'archive')
    state = dict(schema_version=SCHEMA_VERSION, method=method, parameters=nested(settings),
                 **identity, input_fingerprint=digest(identity), platform=platform,
                 primer_fwd=forward, primer_rev=reverse, code=code_fingerprint(),
                 adapter_guard=guard_identity(settings))
    token = digest(state)
    completed = dataset / f'{name}-{method}-run.json'
    checkpoint = dataset / '.checkpoint.json'
    old, old_status = read_state(completed)
    previous, previous_status = read_state(checkpoint)
    # A missing shared checkpoint can be migrated from its old per-method record.
    # A corrupt shared checkpoint is never rescued this way: input ownership is uncertain.
    if previous_status == 'missing' and old_status == 'valid':
        previous, previous_status = old, old_status
    previous_identity = input_identity(previous) if previous_status == 'valid' else None
    same_inputs = previous_identity == identity
    if not same_inputs:
        reason = 'input_identity_changed' if previous_identity is not None else 'input_identity_unknown'
        archive_raw(dataset, previous, previous_identity, identity, reason)
        (dataset / 'platform.txt').unlink(missing_ok=True)
        if not local_source:
            (dataset / 'read_layout.json').unlink(missing_ok=True)
    if old.get('schema_version') != SCHEMA_VERSION or old.get('fingerprint') != token or not same_inputs:
        for suffix in ('table.qza', 'rep-seqs.qza'):
            (dataset / f'{name}-{method}-final-{suffix}').unlink(missing_ok=True)
        (dataset / f'{name}-{method}-primer_info.json').unlink(missing_ok=True)
    if (previous.get('schema_version') != SCHEMA_VERSION or
            previous.get('fingerprint') != token or not same_inputs):
        # Before normalization, ori_fastq may still contain the only raw download.
        # Move that cache out of the staging directory before removing checkpoints.
        original = dataset / 'ori_fastq'
        downloaded = dataset / 'downloaded_fastq'
        if same_inputs and not local_source and original.exists() and not downloaded.exists():
            original.rename(downloaded)
        for folder in ('tmp', 'ori_fastq', 'working_fastq'):
            remove_owned(dataset / folder)
        for suffix in ('quality_status.txt', 'raw_read_counts.tsv'):
            (dataset / f'{name}_{suffix}').unlink(missing_ok=True)
    state['fingerprint'] = token
    atomic_json(completed, state)
    atomic_json(checkpoint, state)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--method', choices=['dada2', 'vsearch'], required=True)
    parser.add_argument('--local-source', default='')
    parser.add_argument('--platform', default='')
    parser.add_argument('--forward', default='')
    parser.add_argument('--reverse', default='')
    args = parser.parse_args()
    try:
        prepare(args.dataset, args.method, args.local_source, args.platform, args.forward, args.reverse)
    except (ValueError, OSError) as error:
        parser.exit(2, f'Pipeline input error: {error}\n')
