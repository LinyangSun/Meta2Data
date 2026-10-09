#!/usr/bin/python3
"""Test-only archive boundary; never install these wrappers into a product PATH."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import shlex
import sys
from urllib.parse import parse_qs, urlsplit

PLATFORM_ACTIONS = {'batch_get_sequencing_platforms', 'get_sequencing_platform'}
DEFAULT_REAL_PYTHON = '/opt/conda/envs/qiime2-amplicon-2024.10/bin/python3'


def trace(tool, argv, action, status, **details):
    target = os.environ.get('M2D_TEST_ARCHIVE_TRACE')
    if not target:
        return
    record = dict(time_utc=datetime.now(timezone.utc).isoformat(), pid=os.getpid(),
                  tool=tool, argv=argv, action=action, status=status, **details)
    payload = (json.dumps(record, sort_keys=True) + '\n').encode()
    fd = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    try:
        os.write(fd, payload)
    finally:
        os.close(fd)


def url_key(url):
    parsed = urlsplit(url if '://' in url else '//' + url)
    if parsed.scheme not in ('', 'ftp', 'http', 'https') or not parsed.netloc:
        raise ValueError('Unsupported archive URL: ' + url)
    if parsed.query or parsed.fragment or parsed.username or parsed.password:
        raise ValueError('Unexpected archive URL components: ' + url)
    return parsed.netloc + parsed.path


class Archive:
    def __init__(self):
        manifest = Path(os.environ['M2D_TEST_ARCHIVE_MANIFEST']).resolve()
        self.data = json.loads(manifest.read_text())
        if self.data.get('synthetic') is not True or self.data.get('schema_version') != 1:
            raise ValueError('Expected a schema_version=1 synthetic fixture manifest')
        self.projects = self.data['projects']
        self.runs = self.data['runs']
        self.files = {}
        for run in self.runs.values():
            for entry in run['files']:
                path = Path(entry['path'])
                entry['_path'] = path if path.is_absolute() else manifest.parent / path
                entry['_url'] = url_key(entry.get('url', 'ftp.sra.ebi.ac.uk/fixtures/' + path.name))
                if entry['_url'] in self.files:
                    raise ValueError('Duplicate fixture file URL: ' + entry['_url'])
                self.files[entry['_url']] = entry

    def platform(self, run, project=None):
        if run not in self.runs:
            raise ValueError('Undeclared Run: ' + run)
        if project and (project not in self.projects or run not in self.projects[project]['runs']):
            raise ValueError('Undeclared project/Run mapping: ' + project + '/' + run)
        platform = self.runs[run]['platform']
        if not platform or any(char.isspace() for char in platform):
            raise ValueError('Invalid fixture platform for ' + run)
        return platform

    def verified_bytes(self, entry):
        payload = entry['_path'].read_bytes()
        actual = hashlib.md5(payload).hexdigest()
        if actual != entry['md5'].lower():
            raise ValueError('Fixture MD5 mismatch: ' + str(entry['_path']))
        return payload

    def filereport(self, project):
        if project not in self.projects:
            raise ValueError('Undeclared BioProject: ' + project)
        output = io.StringIO()
        writer = csv.writer(output, delimiter='\t', lineterminator='\n')
        writer.writerow(['run_accession', 'fastq_ftp', 'fastq_md5', 'library_layout'])
        for accession in self.projects[project]['runs']:
            self.platform(accession, project)
            run = self.runs[accession]
            for entry in run['files']:
                self.verified_bytes(entry)
            writer.writerow([accession, ';'.join(entry['_url'] for entry in run['files']),
                             ';'.join(entry['md5'].lower() for entry in run['files']), run['layout']])
        return output.getvalue().encode()


def parse_wget(argv):
    output, spider, urls = None, False, []
    index = 0
    while index < len(argv):
        arg = argv[index]
        if arg in ('-O', '--output-document'):
            index += 1
            if index == len(argv):
                raise ValueError('Missing wget output path')
            output = argv[index]
        elif arg.startswith('--output-document='):
            output = arg.split('=', 1)[1]
        elif arg.startswith('-O') and len(arg) > 2:
            output = arg[2:]
        elif arg.startswith('-qO') and len(arg) > 3:
            output = arg[3:]
        elif arg == '--spider':
            spider = True
        elif arg in ('-q', '--quiet', '-nv', '--no-verbose'):
            pass
        elif arg.startswith(('--timeout=', '--tries=', '--user-agent=')):
            pass
        elif arg in ('--timeout', '--tries', '--user-agent'):
            index += 1
            if index == len(argv):
                raise ValueError('Missing wget option value')
        elif arg.startswith('-'):
            raise ValueError('Unsupported test wget option: ' + arg)
        else:
            urls.append(arg)
        index += 1
    if len(urls) != 1:
        raise ValueError('Test wget requires exactly one declared archive URL')
    if not spider and output is None:
        raise ValueError('Test wget downloads require -O/--output-document')
    return urls[0], output, spider


def wget(argv):
    url, output, spider = parse_wget(argv)
    archive = Archive()
    parsed = urlsplit(url)
    if parsed.scheme == 'https' and parsed.netloc == 'www.ebi.ac.uk' and parsed.path.rstrip('/') == '/ena/portal/api':
        if not spider or parsed.query:
            raise ValueError('Only the ENA API reachability probe is declared')
        return 'ena_api_probe', dict(url=url)
    if parsed.scheme == 'https' and parsed.netloc == 'www.ebi.ac.uk' and parsed.path == '/ena/portal/api/filereport':
        query = parse_qs(parsed.query)
        required = {'result': ['read_run'], 'format': ['tsv'],
                    'fields': ['run_accession,fastq_ftp,fastq_md5,library_layout']}
        if any(query.get(key) != value for key, value in required.items()) or len(query.get('accession', [])) != 1:
            raise ValueError('Unsupported ENA filereport query')
        payload = archive.filereport(query['accession'][0])
        action = 'ena_filereport'
    else:
        key = url_key(url)
        if key not in archive.files:
            raise ValueError('Undeclared archive URL: ' + url)
        payload = archive.verified_bytes(archive.files[key])
        action = 'archive_file_probe' if spider else 'archive_file_copy'
    if not spider:
        if output == '-':
            sys.stdout.buffer.write(payload)
        else:
            Path(output).write_bytes(payload)
    return action, dict(url=url, output=output, bytes=len(payload))


def platform_action(argv):
    # Product invocations pass the script directly; also tolerate ordinary
    # interpreter flags without redirecting unrelated -c/-m Python programs.
    index = 0
    while index < len(argv) and argv[index] in ('-u', '-B', '-E', '-s', '-S', '-I', '-O', '-OO'):
        index += 1
    if len(argv) > index + 1 and Path(argv[index]).name == 'py_16s.py' and argv[index + 1] in PLATFORM_ACTIONS:
        return argv[index + 1], argv[index + 2:]
    return None, []


def platform(argv):
    action, args = platform_action(argv)
    parser = argparse.ArgumentParser(prog='archive-shim ' + action)
    if action == 'get_sequencing_platform':
        parser.add_argument('--srr_id', required=True)
        parser.add_argument('--bioproject_id')
    else:
        parser.add_argument('--pairs_file', required=True)
    parsed = parser.parse_args(args)
    archive = Archive()
    if action == 'get_sequencing_platform':
        print(archive.platform(parsed.srr_id, parsed.bioproject_id))
        return action, dict(run=parsed.srr_id)
    lines = []
    for line in Path(parsed.pairs_file).read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) not in (2, 3) or not parts[0] or not parts[1]:
            raise ValueError('Malformed platform query pair: ' + line)
        lines.append(parts[0] + '\t' + archive.platform(parts[1], parts[2] if len(parts) == 3 else None))
    # Validate the entire request before publishing any synthetic results.
    if lines:
        print('\n'.join(lines))
    return action, dict(pairs_file=parsed.pairs_file, datasets=len(lines))


def install(directory, interpreter):
    interpreter = Path(interpreter)
    if not interpreter.is_absolute() or not interpreter.is_file():
        raise ValueError('Wrapper interpreter must be an existing absolute Python path')
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    script = str(Path(__file__).resolve())
    for tool in ('wget', 'python', 'python3'):
        wrapper = directory / tool
        wrapper.write_text('#!/bin/sh\nexec ' + ' '.join(map(shlex.quote, (str(interpreter), script, tool))) + ' "$@"\n')
        wrapper.chmod(0o755)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'install':
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('install')
        parser.add_argument('--directory', required=True)
        parser.add_argument('--interpreter', default='/usr/bin/python3')
        args = parser.parse_args()
        install(args.directory, args.interpreter)
        return 0
    if len(sys.argv) < 2 or sys.argv[1] not in ('wget', 'python', 'python3'):
        raise SystemExit('Usage: archive_shim.py {wget|python|python3} [arguments]')
    tool, argv = sys.argv[1], sys.argv[2:]
    try:
        if tool == 'wget':
            action, details = wget(argv)
        elif platform_action(argv)[0]:
            action, details = platform(argv)
        else:
            real = os.environ.get('M2D_TEST_REAL_PYTHON', DEFAULT_REAL_PYTHON)
            if not Path(real).is_absolute() or not Path(real).is_file():
                raise ValueError('M2D_TEST_REAL_PYTHON must name an existing absolute interpreter')
            trace(tool, argv, 'python_passthrough', 'exec', executable=real)
            os.execv(real, [real] + argv)
        trace(tool, argv, action, 'success', **details)
        return 0
    except SystemExit as error:
        code = int(error.code or 0)
        trace(tool, argv, 'platform_arguments', 'error' if code else 'success', exit_code=code)
        return code
    except (OSError, KeyError, ValueError) as error:
        trace(tool, argv, 'rejected', 'error', reason=str(error))
        print('archive-shim: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
