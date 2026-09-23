"""Run and adjudicate one sample's fastp output before publishing stable filenames."""
import argparse
import gzip
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def save_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False) as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')
        temporary = Path(stream.name)
    temporary.replace(path)


def records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as stream:
        while True:
            header = stream.readline()
            if not header:
                return
            sequence, plus, quality = (stream.readline().rstrip('\r\n') for _ in range(3))
            if not header.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                raise ValueError(f'Invalid FASTQ output: {path}')
            name = header.split()[0]
            if name.endswith(('/1', '/2')):
                name = name[:-2]
            yield name, len(sequence)


def validate_outputs(paths, report):
    count = bases = 0
    if len(paths) == 2:
        for first, second in itertools.zip_longest(records(paths[0]), records(paths[1])):
            if first is None or second is None or first[0] != second[0]:
                raise ValueError('Accepted fastp paired outputs have unmatched records')
            count += 2
            bases += first[1] + second[1]
    else:
        for _, length in records(paths[0]):
            count += 1
            bases += length
    expected = report['summary']['after_filtering']
    if count != expected['total_reads'] or bases != expected['total_bases']:
        raise ValueError('Accepted FASTQ counts do not match the fastp JSON report')
    return dict(individual_reads=count, bases=bases,
                count_unit='read_pairs' if len(paths) == 2 else 'reads',
                pipeline_reads=count // len(paths))


def process(args):
    if bool(args.in2) != bool(args.out2):
        raise ValueError('--in2 and --out2 must be provided together')
    if args.threads < 1 or not args.sample_id or '/' in args.sample_id or args.sample_id in ('.', '..'):
        raise ValueError('Invalid thread count or sample ID')
    inputs = [Path(args.in1).resolve()] + ([Path(args.in2).resolve()] if args.in2 else [])
    outputs = [Path(args.out1).absolute()] + ([Path(args.out2).absolute()] if args.out2 else [])
    if any(not source.is_file() for source in inputs):
        raise ValueError('An input FASTQ is missing')
    if len(set(outputs)) != len(outputs) or any(output.resolve() in inputs for output in outputs):
        raise ValueError('Outputs must be distinct and must not overwrite input reads')
    for source, target in zip(inputs, outputs):
        # The staging input may be a link whose resolved name is not canonical.
        if str(source).endswith('.gz') != str(target).endswith('.gz'):
            raise ValueError('Input/output compression suffixes must agree')
    report_dir = Path(args.audit_dir).absolute()
    report_dir.mkdir(parents=True, exist_ok=True)
    work_root = Path(args.work_dir).absolute()
    work_root.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='sample-', dir=work_root))
    decision_path = report_dir / 'decision.json'
    environment = os.environ.copy()
    environment['M2D_PROFILE_SAMPLE'] = args.sample_id
    environment['M2D_PROFILE_DECISION_REPORT'] = str(decision_path)
    commands = []

    def run_fastp(attempt):
        folder = work / attempt
        folder.mkdir()
        temporary_outputs = [folder / output.name for output in outputs]
        json_path = report_dir / f'{attempt}.json'
        html_path = report_dir / f'{attempt}.html'
        command = ['fastp', '-i', str(inputs[0]), '-o', str(temporary_outputs[0])]
        if len(inputs) == 2:
            command += ['-I', str(inputs[1]), '-O', str(temporary_outputs[1])]
            if attempt == 'initial':
                command += ['--detect_adapter_for_pe']
            # The PE fallback retains overlap-based trimming, without de novo inference.
        elif attempt == 'fallback':
            command += ['--disable_adapter_trimming']
        command += ['--disable_quality_filtering', '--disable_length_filtering',
                    '-w', str(args.threads), '-j', str(json_path), '-h', str(html_path)]
        commands.append(dict(attempt=attempt, command=command))
        save_json(report_dir / 'commands.json', commands)
        run_environment = dict(environment, M2D_PROFILE_STAGE='fastp_' + attempt)
        subprocess.run(command, env=run_environment, check=True)
        report = json.loads(json_path.read_text())
        if not isinstance(report.get('summary', {}).get('after_filtering'), dict):
            raise ValueError(f'Missing fastp summary in {json_path}')
        return temporary_outputs, report, json_path, html_path

    try:
        initial = run_fastp('initial')
        command = ['python3', str(Path(__file__).with_name('fastp_guard.py')), 'check',
                   '--report', str(initial[2]), '--db-manifest', args.db_manifest,
                   '--cache-dir', args.cache_dir, '--output-dir', str(report_dir),
                   '--sample-id', args.sample_id, '--threads', str(args.threads),
                   '--min-length', str(args.min_length), '--min-identity', str(args.min_identity),
                   '--min-coverage', str(args.min_coverage), '--auto-r1']
        if len(inputs) == 2:
            command += ['--auto-r2']
        commands.append(dict(attempt='guard', command=command))
        save_json(report_dir / 'commands.json', commands)
        check_env = dict(environment, M2D_PROFILE_STAGE='fastp_adapter_guard')
        subprocess.run(command, env=check_env, check=True)
        decision = json.loads(decision_path.read_text())
        if decision.get('action') not in ('accept', 'fallback'):
            raise ValueError('Guard did not produce a valid processing action')
        accepted_attempt = 'fallback' if decision['action'] == 'fallback' else 'initial'
        accepted = run_fastp('fallback') if accepted_attempt == 'fallback' else initial
        counts = validate_outputs(accepted[0], accepted[1])
        decision.update(sample_id=args.sample_id, accepted_attempt=accepted_attempt,
                        input_files=list(map(str, inputs)), final_outputs=list(map(str, outputs)),
                        initial_report=str(initial[2]), accepted_report=str(accepted[2]),
                        initial_summary=initial[1]['summary'], accepted_summary=accepted[1]['summary'],
                        accepted_counts=counts, publication_status='validated',
                        fallback_policy=('paired_overlap_only' if len(inputs) == 2 else 'single_no_adapter_trimming')
                        if accepted_attempt == 'fallback' else None)
        save_json(decision_path, decision)
        # Downstream processing starts only after this process and the dataset loop finish.
        for source, target in zip(accepted[0], outputs):
            target.parent.mkdir(parents=True, exist_ok=True)
            os.replace(source, target)
        for source, target in ((accepted[2], Path(args.report_json)), (accepted[3], Path(args.report_html))):
            target.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as stream:
                staged = Path(stream.name)
            shutil.copyfile(source, staged)
            staged.replace(target)
        decision['publication_status'] = 'complete'
        save_json(decision_path, decision)
        shutil.rmtree(work)
        print(f"[fastp guard] {args.sample_id}: {decision['status']}; accepted={accepted_attempt}; "
              f"{counts['pipeline_reads']} {counts['count_unit']}")
    except Exception as error:
        failure = dict(sample_id=args.sample_id, status='error', error=str(error),
                       input_files=list(map(str, inputs)), work_directory=str(work))
        save_json(report_dir / 'failure.json', failure)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('sample-id', 'in1', 'out1', 'report-json', 'report-html', 'audit-dir',
                 'work-dir', 'db-manifest', 'cache-dir'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--in2')
    parser.add_argument('--out2')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--min-length', type=int, default=50)
    parser.add_argument('--min-identity', type=float, default=98)
    parser.add_argument('--min-coverage', type=float, default=95)
    args = parser.parse_args()
    try:
        process(args)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        parser.exit(2, f'fastp guard failed for {args.sample_id}: {error}\n')


if __name__ == '__main__':
    main()
