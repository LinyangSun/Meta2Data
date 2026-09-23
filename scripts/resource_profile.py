"""Optional command-level resource measurements; raw events are the authority.

CPU is Linux wait4 accounting (inclusive of waited-for descendants). Command
hierarchies must not be summed. Process RSS is sampled and may count shared pages
twice; cgroup observations describe the enclosing job/step, not a private command
cgroup. No counters are reset. Unsupported values are null/NA, never invented 0s.
"""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid

TOOLS = ('python', 'python3', 'qiime', 'vsearch', 'fastp', 'cutadapt', 'seqkit',
         'chopper', 'minimap2', 'racon', 'mafft', 'FastTree', 'fasttree',
         'prefetch', 'fasterq-dump', 'wget', 'curl', 'gzip', 'pigz', 'blastn', 'makeblastdb')


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name('.' + path.name + '.' + uuid.uuid4().hex)
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def stage_for(argv, module, method):
    tool = Path(argv[0]).name
    args = argv[1:]
    prefix = method + '_' if method in ('dada2', 'vsearch') else ''
    if tool in ('wget', 'curl', 'prefetch', 'fasterq-dump'):
        if '--spider' in args or any('/ena/portal/api/' in a for a in args):
            return 'download_source_probe', 'metadata'
        return ('database_download' if module == 'taxa' else 'read_download'), 'download'
    if tool == 'qiime' and len(args) > 1:
        action = '_'.join(args[:2]).replace('-', '_')
        mapping = {'feature_table_merge': 'taxa_merge_tables',
                   'feature_table_merge_seqs': 'taxa_merge_sequences',
                   'rescript_orient_seqs': 'taxa_orientation',
                   'feature_classifier_classify_sklearn': 'taxa_classification',
                   'fragment_insertion_sepp': 'taxa_sepp',
                   'fragment_insertion_filter_features': 'taxa_tree_filter'}
        if module == 'taxa':
            if action == 'feature_table_filter_features':
                return 'taxa_orientation_table_filter', 'compute'
            return mapping.get(action, 'taxa_' + action), ('reporting' if 'summar' in action else 'compute')
        return (action if action.startswith(prefix) else prefix + action), ('reporting' if 'summar' in action else 'compute')
    if tool == 'vsearch':
        actions = {'--fastq_mergepairs': 'merge_pairs', '--fastq_filter': 'quality_filter',
                   '--derep_fulllength': 'dereplicate', '--fastx_uniques': 'dereplicate',
                   '--sortbysize': 'abundance_filter', '--cluster_size': 'precluster',
                   '--cluster_unoise': 'denoise', '--uchime3_denovo': 'chimera_filter',
                   '--uchime_denovo': 'chimera_filter', '--cluster_fast': 'cluster',
                   '--usearch_global': 'map', '--fastx_filter': 'length_filter'}
        return 'vsearch_' + next((v for k, v in actions.items() if k in args), 'other'), 'compute'
    if tool == 'blastn':
        return 'reference_blast', 'auditing'
    if tool == 'makeblastdb':
        return 'blast_database_build', 'preparation'
    if tool in ('gzip', 'pigz'):
        return ('fastq_integrity_check', 'auditing') if '-t' in args or '--test' in args else ('data_compression', 'preparation')
    if tool in ('python', 'python3') and args:
        script = Path(args[0]).name
        action = args[1] if len(args) > 1 else ''
        if script == 'fastp_checked.py':
            return 'fastp_checked', 'compute'
        if script == 'fastp_guard.py':
            return ('fastp_adapter_guard_prepare', 'preparation') if action == 'prepare' else ('fastp_adapter_guard', 'auditing')
        if script == 'download_integrity.py':
            if action in ('ena-manifest', 'ncbi-manifest'):
                return 'download_file_manifest', 'metadata'
            return 'fastq_integrity_check', 'auditing'
        if script == 'read_counts.py':
            return 'read_count_audit_' + action, 'auditing'
        if script == 'py_16s.py':
            category = 'compute'
            if action in ('GenerateDatasetsIDsFile', 'GenerateSRAsFile', 'get_sequencing_platform', 'batch_get_platforms', 'batch_get_sequencing_platforms'):
                category = 'metadata'
            elif action == 'build_per_dataset_summary':
                category = 'reporting'
            return prefix + action, category
        if script == 'read_layout.py':
            return 'read_layout_' + action, 'auditing' if action == 'count' else 'preparation'
        if script.endswith('.py'):
            return prefix + script[:-3], 'compute'
        return 'python_helper', 'orchestration'
    return {'fastp': 'fastp', 'cutadapt': 'primer_trim'}.get(tool, prefix + tool), 'compute'


def bootstrap(directory):
    """PATH shims also see tools invoked from Python, unlike shell functions."""
    import psutil  # Fail before starting work if the optional dependency is absent.
    target = Path(directory).resolve()
    target.mkdir(parents=True, exist_ok=True)
    python = str(Path(sys.executable).absolute())
    script = str(Path(__file__).resolve())
    resolved = {}
    for tool in TOOLS:
        real = shutil.which(tool)
        if not real:
            continue
        if Path(real).parent == target:
            raise ValueError('Profiler bootstrap must use the original PATH')
        resolved[tool] = real
        shim = target / tool
        # Absolute interpreter prevents the shim from invoking itself.
        shim.write_text('#!/bin/bash\nexec ' + ' '.join(map(shlex.quote, (
            python, script, 'run', '--', real))) + ' "$@"\n')
        shim.chmod(0o755)
    write_json(target / 'executables.json', dict(python=python, scripts=script,
                                                psutil=psutil.__version__, tools=resolved))


def state_context(argv):
    path = os.environ.get('READ_COUNTS_STATE', '')
    if '--state' in argv:
        path = argv[argv.index('--state') + 1]
    state = {}
    if path and Path(path).is_file():
        state = json.loads(Path(path).read_text())
    return path, state


def cgroup_path():
    try:
        for line in Path('/proc/self/cgroup').read_text().splitlines():
            if line.startswith('0::'):
                path = Path('/sys/fs/cgroup') / line[3:].lstrip('/')
                if (path / 'memory.current').exists():
                    return path
    except OSError:
        pass
    return None


def cgroup_snapshot(path):
    if path is None:
        return {}
    result = {}
    try:
        result['job_cgroup_bytes'] = int((path / 'memory.current').read_text())
        mem = dict(line.split() for line in (path / 'memory.stat').read_text().splitlines())
        result['job_cgroup_anon_bytes'] = int(mem['anon'])
        cpu = dict(line.split() for line in (path / 'cpu.stat').read_text().splitlines())
        result['job_cgroup_cpu_usec'] = int(cpu['usage_usec'])
    except (OSError, ValueError, KeyError):
        pass
    return result


def option_value(argv, *names):
    for index, token in enumerate(argv[1:], 1):
        if token in names and index + 1 < len(argv):
            return argv[index + 1]
        for name in names:
            if token.startswith(name + '='):
                return token[len(name) + 1:]
    return None


# BLAST creates a family of index files rather than a file at the -db/-out
# prefix. Only known suffixes count; an unrelated prefix.fasta is not an index.
BLAST_INDEX_SUFFIX = re.compile(r'(?:\.\d+)?\.(?:nhr|nin|nsq|ndb|not|ntf|nto|nog|nos|nsd|nsi|nhd|nhi|nal|phr|pin|psq|pdb|pot|ptf|pto|pog|pos|psd|psi|phd|phi|pal)$')


def file_arguments(argv):
    """Explicit files and BLAST index prefixes; no recursive filesystem scans."""
    found = {}
    tool = Path(argv[0]).name
    output_flags = {'-o', '-O', '--output', '--out', '--centroids',
                    '--fastaout', '--fastqout', '--notmerged_fwd', '--notmerged_rev',
                    '--nonchimeras', '--chimeras', '--otutabout', '--uc', '--output_csv',
                    '--out1', '--out2', '--report-json', '--report-html'}
    input_flags = {'-i', '-I', '--input', '--input-path', '--FilePath', '--reference', '--db',
                   '--in1', '--in2', '--report', '--db-manifest'}
    if tool == 'fastp':
        output_flags.update(('-h', '-j', '--html', '--json', '--out1', '--out2',
                             '--unpaired1', '--unpaired2', '--merged_out', '--failed_out'))
    if tool == 'blastn':
        output_flags.add('-out')
        input_flags.update(('-query', '-subject', '-gilist', '-seqidlist'))
    if tool == 'makeblastdb':
        input_flags.update(('-in', '-taxid_map'))
    for i, token in enumerate(argv[1:], 1):
        if token.startswith('-') or token in ('/dev/null', '/dev/stdout', '/dev/stderr'):
            continue
        prev = argv[i - 1]
        if '\n' in token or len(token) > 4096:
            continue
        path = Path(token)
        if tool == 'blastn' and prev == '-db':
            found[str(path.absolute())] = 'input_index_prefix'
            continue
        if tool == 'makeblastdb' and prev == '-out':
            found[str(path.absolute())] = 'output_index_prefix'
            continue
        try:
            exists = path.is_file()
        except OSError:
            continue
        is_output = prev.startswith('--o-') or prev in output_flags
        if exists or is_output or prev in input_flags or prev.startswith('--i-'):
            found[str(path.absolute())] = 'output' if is_output else 'input'
    if tool == 'makeblastdb' and not option_value(argv, '-out'):
        prefix = option_value(argv, '-in')
        if prefix:
            # The input FASTA can also be the implicit output index prefix.
            found[str(Path(prefix).absolute()) + '/.blast-index-prefix'] = 'output_index_prefix'
    return found


def sizes(paths, kind):
    observations, seen = {}, set()
    for name, role in paths.items():
        if role == kind:
            candidates = [Path(name)]
        elif role == kind + '_index_prefix':
            path = Path(name)
            if path.name == '.blast-index-prefix':
                path = path.parent
            try:
                candidates = [item for item in path.parent.iterdir()
                              if item.name.startswith(path.name) and
                              BLAST_INDEX_SUFFIX.fullmatch(item.name[len(path.name):])]
            except OSError:
                candidates = []
        else:
            continue
        for path in candidates:
            try:
                resolved = path.resolve()
                if path.is_file() and resolved not in seen:
                    observations[str(path.absolute())] = path.stat().st_size
                    seen.add(resolved)
            except OSError:
                pass
    return observations


def report_identity(path):
    try:
        stat = Path(path).stat()
        return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    except OSError:
        return None


def fastp_report_snapshot(command, events, invocation, prior_identity=None):
    """Capture the invocation's own report before the caller can replace it.

    fastp JSON counts individual reads, including for PE input. Convert to pairs
    only for a synchronized paired-output command with valid even totals. These
    are attempt counts, and can differ from the accepted read-count ledger.
    """
    if Path(command[0]).name != 'fastp':
        return {}
    report_name = option_value(command, '-j', '--json') or 'fastp.json'
    path = Path(report_name).absolute()
    result = {'fastp_report_path': str(path), 'fastp_report_status': 'missing',
              'fastp_counts_scope': 'this invocation only; not necessarily accepted downstream'}
    try:
        data = path.read_bytes()
        result['fastp_report_sha256'] = hashlib.sha256(data).hexdigest()
        snapshot = events / (invocation + '.fastp-report.json')
        snapshot.write_bytes(data)
        result['fastp_report_snapshot'] = str(snapshot)
        if prior_identity is not None and report_identity(path) == prior_identity:
            result['fastp_report_status'] = 'unchanged_preexisting_report'
            return result
        report = json.loads(data)
        summary = report['summary']
        values = {}
        for label, section in (('input', 'before_filtering'), ('output', 'after_filtering')):
            for name, field in (('individual_reads', 'total_reads'), ('bases', 'total_bases')):
                value = summary[section][field]
                if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                    raise ValueError(f'Invalid {section}.{field}')
                values[label + '_' + name] = value
        result.update({label + '_mean_read_length': (values[label + '_bases'] / values[label + '_individual_reads'] if values[label + '_individual_reads'] else None) for label in ('input', 'output')})
        result.update(values, fastp_report_status='valid',
                      fastp_report_summary=summary,
                      read_link_backend='immutable fastp invocation JSON snapshot')
        paired = option_value(command, '-I', '--in2') is not None
        paired_output = option_value(command, '-O', '--out2') is not None
        unsupported = any(flag in command for flag in ('--interleaved_in', '--merge', '-m',
                          '--unpaired1', '--unpaired2', '--merged_out', '--include_unmerged'))
        if not paired and not paired_output and not unsupported:
            divisor, unit = 1, 'reads'
        elif paired and paired_output and not unsupported and all(values[key] % 2 == 0 for key in ('input_individual_reads', 'output_individual_reads')):
            # If mate summaries are supplied, they must agree with a paired interpretation.
            for section, key in (('before_filtering', 'input_individual_reads'), ('after_filtering', 'output_individual_reads')):
                for mate in ('read1_', 'read2_'):
                    mate_reads = report.get(mate + section, {}).get('total_reads')
                    if mate_reads is not None and (isinstance(mate_reads, bool) or not isinstance(mate_reads, int) or mate_reads != values[key] // 2):
                        raise ValueError('fastp mate counts do not describe synchronized pairs')
            divisor, unit = 2, 'read_pairs'
        else:
            result['fastp_pair_conversion_status'] = 'unsupported_or_inconsistent_layout; individual reads retained'
            return result
        result.update(input_reads=values['input_individual_reads'] // divisor,
                      output_reads=values['output_individual_reads'] // divisor,
                      read_unit=unit, read_output_stage='fastp_attempt',
                      fastp_pair_conversion_status='paired' if divisor == 2 else 'single')
    except FileNotFoundError:
        pass
    except (OSError, ValueError, KeyError, TypeError) as error:
        result.update(fastp_report_status='invalid', fastp_report_error=str(error))
    return result


def attach_fastp_decision(row):
    """Attach acceptance only when the current report matches the captured bytes."""
    if 'fastp_report_path' not in row:
        return
    row['fastp_attempt_decision'] = 'unresolved'
    decision_path = row.get('decision_report')
    if not decision_path or row.get('fastp_report_status') != 'valid' or row.get('status') != 'success':
        return
    try:
        decision = json.loads(Path(decision_path).read_text())
        current = Path(row['fastp_report_path']).read_bytes()
        if hashlib.sha256(current).hexdigest() != row.get('fastp_report_sha256'):
            return
        if decision.get('publication_status') != 'complete':
            return
        if decision.get('sample_id') != row.get('SampleName'):
            return
        accepted = str(Path(decision['accepted_report']).absolute())
        if accepted == row['fastp_report_path']:
            row['fastp_attempt_decision'] = 'accepted'
        elif str(Path(decision['initial_report']).absolute()) == row['fastp_report_path']:
            row['fastp_attempt_decision'] = 'rejected'
    except (OSError, ValueError, KeyError, TypeError):
        return


def run(command, stage=None, scope='command'):
    import psutil
    root = Path(os.environ['M2D_PROFILE_DIR']).resolve()
    events = root / 'events'
    events.mkdir(parents=True, exist_ok=True)
    invocation = uuid.uuid4().hex
    parent = os.environ.get('M2D_PROFILE_INVOCATION_ID', '')
    module = os.environ.get('M2D_PROFILE_MODULE', '')
    method = os.environ.get('M2D_PROFILE_METHOD', os.environ.get('MODE', ''))
    mapped_stage, category = stage_for(command, module, method)
    path, state = state_context(command)
    args_paths = file_arguments(command)
    inputs = sizes(args_paths, 'input')
    interval = float(os.environ.get('M2D_PROFILE_INTERVAL', '1'))
    if not 0.05 <= interval <= 60:
        raise ValueError('M2D_PROFILE_INTERVAL must be between 0.05 and 60 seconds')
    sample_ids = [s for s in state.get('samples', {}) if any(s in a for a in command)]
    explicit_sample = os.environ.get('M2D_PROFILE_SAMPLE') or option_value(command, '--sample-id')
    stage_override = os.environ.get('M2D_PROFILE_STAGE')
    record = dict(schema_version=1, execution_id=os.environ.get('M2D_EXECUTION_ID', root.name),
                  invocation_id=invocation, parent_invocation_id=parent,
                  SlurmJobID=os.environ.get('SLURM_JOB_ID', ''), hostname=os.uname().nodename,
                  BioProject=os.environ.get('M2D_PROFILE_PROJECT') or state.get('dataset', ''),
                  SampleName=explicit_sample or (sample_ids[0] if len(sample_ids) == 1 else ''),
                  module=module, method=method, stage=stage or stage_override or mapped_stage,
                  decision_report=os.environ.get('M2D_PROFILE_DECISION_REPORT', ''),
                  category=scope if scope in ('module', 'dataset') else module if module in ('auditing', 'reporting') else category,
                  scope=scope, attempt=state.get('attempt', ''), read_counts_state=path,
                  start_utc=utc(), command=command, cwd=os.getcwd(),
                  allocated_cpus=os.environ.get('M2D_PROFILE_CPUS') or os.environ.get('SLURM_CPUS_PER_TASK', ''),
                  slurm_allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK', ''),
                  sampling_interval_s=interval, input_files=inputs,
                  input_bytes=sum(inputs.values()) if inputs else None,
                  input_size_scope='explicit file arguments and BLAST index-prefix files only', status='running',
                  cpu_backend='Linux wait4; inclusive waited-for descendants',
                  memory_backend='wait4 maxrss + sampled process-tree RSS + enclosing cgroup context',
                  io_backend='wait4 block I/O + sampled root /proc I/O including waited-for children (lower bound, no PID sums)',
                  cache_state='cached_result_reused' if mapped_stage == 'read_count_audit_recover' else
                  'checkpoint_counts_inherited' if mapped_stage == 'read_count_audit_inherit' else 'executed')
    write_json(events / (invocation + '.start.json'), record)
    env = dict(os.environ, M2D_PROFILE_INVOCATION_ID=invocation)
    # An override describes only this process; nested tools need their own stage.
    env.pop('M2D_PROFILE_STAGE', None)
    group = cgroup_path()
    before = cgroup_snapshot(group)
    record['cgroup_path'] = str(group) if group else None
    record['cgroup_scope'] = 'enclosing job/step context; includes concurrent siblings, cache and profiler'
    record['memory_caveat'] = 'sampled RSS sums may double-count shared pages and miss short peaks; maxrss is not simultaneous tree RSS'
    prior_report = report_identity(option_value(command, '-j', '--json') or 'fastp.json') if Path(command[0]).name == 'fastp' else None
    started = time.monotonic()
    try:
        # Native dataset workers use fd 3 for console milestones while their
        # stdout/stderr go to the project log. Preserve only that explicit fd.
        pass_fds = ()
        if scope == 'dataset':
            try:
                os.fstat(3)
                pass_fds = (3,)
            except OSError:
                pass
        child = subprocess.Popen(command, env=env, start_new_session=True,
                                 pass_fds=pass_fds)
    except OSError as exc:
        record.update(status='launch_failed', exit_code=127, error=str(exc),
                      end_utc=utc(), wall_s=time.monotonic() - started)
        write_json(events / (invocation + '.json'), record)
        return 127
    stop = threading.Event()
    observations = []
    # Linux /proc/PID/io rolls waited-for child counters into the parent.
    # Summing historical per-PID maxima would count those bytes twice. Observe
    # only the command root here, and retain wait4 block-I/O as the exit measure.
    root_io = [None, None]
    monitor_errors = []
    def monitor():
        try:
            proc = psutil.Process(child.pid)
            with gzip.open(events / (invocation + '.samples.jsonl.gz'), 'wt') as stream:
                while True:
                    obs = dict(t_monotonic_s=time.monotonic() - started, **cgroup_snapshot(group))
                    rss, processes = 0, 0
                    try:
                        tree = [proc] + proc.children(recursive=True)
                    except psutil.Error:
                        tree = []
                    for item in tree:
                        try:
                            rss += item.memory_info().rss
                            processes += 1
                            if item.pid == child.pid:
                                counters = item.io_counters()
                                root_io[0] = max(root_io[0] or 0, counters.read_bytes)
                                root_io[1] = max(root_io[1] or 0, counters.write_bytes)
                                obs.update(root_process_read_bytes=counters.read_bytes,
                                           root_process_write_bytes=counters.write_bytes)
                        except (psutil.Error, OSError, AttributeError):
                            pass
                    obs.update(process_tree_rss_bytes=rss, process_count=processes)
                    observations.append(obs)
                    stream.write(json.dumps(obs) + '\n')
                    stream.flush()
                    if stop.wait(interval):
                        break
        except Exception as exc:
            monitor_errors.append(str(exc))
    worker = threading.Thread(target=monitor, daemon=True)
    old_handlers = {}
    def forward(signum, frame):
        try:
            os.killpg(child.pid, signum)
        except ProcessLookupError:
            pass
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        old_handlers[sig] = signal.signal(sig, forward)
    worker.start()
    try:
        _, status, usage = os.wait4(child.pid, 0)
        child.returncode = os.waitstatus_to_exitcode(status)
    finally:
        stop.set()
        worker.join(timeout=5)
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
    elapsed = time.monotonic() - started
    code = child.returncode if child.returncode >= 0 else 128 - child.returncode
    after = cgroup_snapshot(group)
    record.update(fastp_report_snapshot(command, events, invocation, prior_report))
    if isinstance(record.get('input_reads'), int) and elapsed > 0:
        record['input_reads_per_s'] = record['input_reads'] / elapsed
    outputs = sizes(args_paths, 'output')
    def peak(key):
        values = [v[key] for v in observations if key in v]
        return max(values) if values else None
    cpu_delta = (after['job_cgroup_cpu_usec'] - before['job_cgroup_cpu_usec']) / 1e6 if 'job_cgroup_cpu_usec' in before and 'job_cgroup_cpu_usec' in after else None
    record.update(end_utc=utc(), wall_s=elapsed, exit_code=code,
                  status='success' if code == 0 else 'failed',
                  cpu_user_s=usage.ru_utime, cpu_system_s=usage.ru_stime,
                  process_maxrss_kib=usage.ru_maxrss,
                  sampled_peak_process_tree_rss_bytes=peak('process_tree_rss_bytes'),
                  sampled_peak_job_cgroup_bytes=peak('job_cgroup_bytes'),
                  sampled_peak_job_cgroup_anon_bytes=peak('job_cgroup_anon_bytes'),
                  job_cgroup_cpu_delta_s=cpu_delta,
                  rusage_block_read_bytes=usage.ru_inblock * 512,
                  rusage_block_write_bytes=usage.ru_oublock * 512,
                  sampled_process_read_bytes=root_io[0],
                  sampled_process_write_bytes=root_io[1],
                  output_files=outputs, output_bytes=sum(outputs.values()) if outputs else None,
                  monitoring_samples=len(observations), monitoring_errors=monitor_errors)
    write_json(events / (invocation + '.json'), record)
    return code


def export(root):
    root = Path(root)
    records = []
    for start in sorted((root / 'events').glob('*.start.json')):
        final = start.with_name(start.name.replace('.start.json', '.json'))
        row = json.loads((final if final.exists() else start).read_text())
        if not final.exists():
            row['status'] = 'incomplete'
        records.append(row)
    records.sort(key=lambda r: (r['start_utc'], r['invocation_id']))
    children = {r['parent_invocation_id'] for r in records}
    for row in records:
        row['has_recorded_children'] = row['invocation_id'] in children
        row['additivity'] = 'inclusive; do not sum ancestors and descendants'
        attach_fastp_decision(row)
    # Read counts remain authoritative in the existing ledger. Resource rows
    # link to the next matching observation in the SAME attempt. Missing,
    # failed, inherited and unobserved stages never acquire fabricated counts.
    mapping = {
        'fastp': 'fastp_reads', 'fastp_checked': 'fastp_reads', 'primer_trim': 'primer_trimmed_reads',
        'vsearch_dereplicate': 'vsearch_dereplicated_reads',
        'vsearch_abundance_filter': 'vsearch_abundance_filtered_reads',
        'vsearch_precluster': 'vsearch_preclustered_reads',
        'vsearch_denoise': 'vsearch_denoised_reads',
        'vsearch_chimera_filter': 'vsearch_nonchimeric_reads',
        'vsearch_cluster': 'vsearch_clustered_reads',
        'vsearch_map': 'vsearch_mapped_reads',
        'vsearch_merge_pairs': 'vsearch_merge_attempt_reads',
        'vsearch_quality_filter': 'vsearch_quality_filtered_reads',
        'vsearch_sanitize_fastq': 'sanitized_reads',
        'dada2_sanitize_fastq': 'sanitized_reads',
        'taxa_orientation_table_filter': 'taxa_oriented_reads',
        'taxa_classification': 'taxa_classified_reads',
        'taxa_tree_filter': 'taxa_tree_placed_reads',
        'dada2_denoise_paired': 'dada2_nonchimeric_reads',
        'dada2_denoise_single': 'dada2_nonchimeric_reads',
        'dada2_denoise_pyro': 'dada2_nonchimeric_reads',
        'dada2_denoise_ccs': 'dada2_nonchimeric_reads'}
    ledgers = {}
    for row in records:
        state_path = row.get('read_counts_state')
        if not state_path or state_path in ledgers or not Path(state_path).is_file():
            continue
        ledgers[state_path] = [(p, json.loads(p.read_text())) for p in sorted(Path(state_path).parent.glob('event-*.json'))]
    for row in records:
        if row.get('status') != 'success' or row.get('read_counts_state') not in ledgers:
            continue
        # A rejected fastp attempt must never inherit the accepted final audit.
        # Every fastp invocation retains its own immutable report counts instead.
        if 'fastp_report_path' in row:
            continue
        target = mapping.get(row['stage'])
        if not target:
            continue
        candidates = ledgers[row['read_counts_state']]
        sample = row.get('SampleName')
        for path, event in candidates:
            if event['stage'] != target or event.get('inherited_from') or event['time'] < row['end_utc']:
                continue
            if sample and sample not in event['counts']:
                continue
            # Per-sample FASTQ events can summarize many commands. Require an
            # actual output file under the audited directory when possible.
            if target in ('fastp_reads', 'primer_trimmed_reads'):
                source = Path(event['source'])
                if not any(Path(p) == source or source in Path(p).parents for p in row.get('output_files', {})):
                    continue
            value = event['counts'].get(sample, 'NA') if sample else event['total']
            row.update(output_reads=value, read_unit=event['basis'], read_output_stage=target,
                       read_event_path=str(path), read_link_backend='same attempt, stage mapping, next audit; command is not itself a read counter')
            parent_stage = 'dada2_input_reads' if row['stage'].startswith('dada2_denoise_') else event['input_stage']
            prior = [e for _, e in candidates if e['stage'] == parent_stage and e['time'] <= event['time']]
            if prior and prior[-1]['basis'] == event['basis']:
                row['input_reads'] = prior[-1]['counts'].get(sample, 'NA') if sample else prior[-1]['total']
                if isinstance(row['input_reads'], int) and row.get('wall_s', 0) > 0:
                    row['input_reads_per_s'] = row['input_reads'] / row['wall_s']
            break
    preferred = ['execution_id', 'invocation_id', 'parent_invocation_id', 'SlurmJobID',
                 'BioProject', 'SampleName', 'module', 'method', 'stage', 'category',
                 'scope', 'attempt', 'start_utc', 'end_utc', 'wall_s', 'cpu_user_s',
                 'cpu_system_s', 'allocated_cpus', 'process_maxrss_kib',
                 'sampled_peak_process_tree_rss_bytes', 'input_bytes', 'output_bytes',
                 'status', 'exit_code']
    fields = preferred + sorted({k for r in records for k in r} - set(preferred))
    def table(name, rows):
        with (root / name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fields, restval='NA')
            writer.writeheader()
            for row in rows:
                writer.writerow({k: 'NA' if v is None else json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in row.items()})
    table('step_resources.csv', records)
    table('dataset_resources.csv', [r for r in records if r['module'] == 'pip' and
          (r['scope'] == 'dataset' or (r['scope'] == 'module' and r.get('BioProject') and
                                     r['stage'] == 'pip_total'))])
    table('pipeline_resources.csv', [r for r in records if r['scope'] == 'module' and r['module'] == 'pip'])
    table('taxa_resources.csv', [r for r in records if r['scope'] == 'module' and r['module'] == 'taxa'])
    write_json(root / 'resource_export.json', dict(records=len(records),
               incomplete=sum(r['status'] == 'incomplete' for r in records), exported_at_utc=utc()))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest='action', required=True)
    init = subs.add_parser('bootstrap')
    init.add_argument('--directory', required=True)
    execute = subs.add_parser('run')
    execute.add_argument('--stage')
    execute.add_argument('--scope', default='command', choices=['command', 'module', 'dataset'])
    execute.add_argument('command', nargs=argparse.REMAINDER)
    finish = subs.add_parser('export')
    finish.add_argument('--directory', required=True)
    args = parser.parse_args()
    if args.action == 'bootstrap':
        bootstrap(args.directory)
    elif args.action == 'export':
        export(args.directory)
    else:
        command = args.command[1:] if args.command[:1] == ['--'] else args.command
        if not command:
            parser.error('run requires a command after --')
        sys.exit(run(command, args.stage, args.scope))


if __name__ == '__main__':
    main()
