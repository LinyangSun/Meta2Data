"""Offline integration checks: real CLI/runner with metadata and tool boundary stubs."""
import csv
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]

METADATA_STUB = '''import csv, json, sys
from pathlib import Path
args = sys.argv[1:]
def value(flag): return args[args.index(flag)+1]
command = args[0]
if command in ('GenerateDatasetsIDsFile', 'GenerateSRAsFile', 'subset_meta_for_test'):
    with open(value('--FilePath')) as handle:
        reader = csv.DictReader(handle); fields = reader.fieldnames; rows = list(reader)
    column = value('--Bioproject'); out = Path(value('--OutputDir')); out.mkdir(parents=True, exist_ok=True)
    projects = list(dict.fromkeys(row[column] for row in rows))
    if command == 'GenerateDatasetsIDsFile':
        (out/'datasets_ID.txt').write_text('\\n'.join(projects)+'\\n')
    elif command == 'GenerateSRAsFile':
        for project in projects:
            folder=out/project; folder.mkdir(exist_ok=True)
            (folder/(project+'_sra.txt')).write_text(''.join(row[value('--SRA_Number')]+'\\t'+project+'_'+row[value('--SRA_Number')]+'\\n' for row in rows if row[column]==project))
    else:
        selected=[]; counts={}
        for row in rows:
            project=row[column]; counts[project]=counts.get(project,0)+1
            if counts[project]<=2: selected.append(row)
        dest=out/'test_subset_meta.csv'
        with dest.open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields); writer.writeheader(); writer.writerows(selected)
        print(dest)
elif command == 'batch_get_sequencing_platforms':
    Path(value('--pairs_file')).with_name('platform_queries.txt').write_text(Path(value('--pairs_file')).read_text())
    for line in Path(value('--pairs_file')).read_text().splitlines():
        print(line.split('\\t')[0]+'\\tUNSUPPORTED_TEST_PLATFORM')
'''
READ_COUNTS_STUB = '''import json,os,sys
from pathlib import Path
args=sys.argv[1:]
if args[0]=='begin':
    dataset=Path(args[args.index('--dataset')+1])
    (dataset/'worker_primers.json').write_text(json.dumps({'forward':os.environ.get('PRIMER_FWD',''),'reverse':os.environ.get('PRIMER_REV','')}))
    (dataset/'worker_threads.txt').write_text(os.environ.get('THREADS_PER_DATASET',''))
    (dataset/'worker_context.json').write_text(json.dumps({'local':os.environ.get('LOCAL_MODE',''),'platform':os.environ.get('LOCAL_PLATFORM','')}))
    print(dataset/'read_counts/state.json')
'''


class PipCliTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.repo = self.work / 'repo'
        self.scripts = self.repo / 'scripts'
        self.scripts.mkdir(parents=True)
        (self.repo / 'bin').mkdir()
        shutil.copy(ROOT / 'bin/Meta2Data-AmpliconPIP', self.repo / 'bin')
        for name in ('parameters.py', 'project_primers.py', 'read_layout.py', 'pip_state.py', 'run.sh', 'local_datasets.py'):
            shutil.copy(ROOT / 'scripts' / name, self.scripts / name)
        (self.scripts / 'resource_profile.sh').write_text('')
        (self.scripts / 'AmpliconFunction.sh').write_text('Audit_Exit() { return 0; }\n')
        (self.scripts / 'py_16s.py').write_text(METADATA_STUB)
        (self.scripts / 'read_counts.py').write_text(READ_COUNTS_STUB)
        self.metadata = self.work / 'metadata.csv'
        self.metadata.write_text('Bioproject,Run\nPRJNA1,SRR1\nPRJNA2,ERR2\nPRJNA1,SRR3\nPRJNA1,SRR4\nPRJNA9,DRR9\n')
        self.output = self.work / 'results/pip'
        self.local_columns = ['--local-datasets-colNAME', 'datasets', '--local-path-colNAME', 'path', '--local-platform-colNAME', 'platform']
        self.env = {key: value for key, value in os.environ.items()
                    if not key.startswith(('M2D_', 'PRIMER_', 'ADAPTER_', 'VSEARCH_', 'TAXA_'))}
        self.base = ['--public-m', str(self.metadata), '--public-bioproject-colNAME', 'Bioproject', '--public-sra-colNAME', 'Run',
                     '--vsearch', '-t', '2', '--no-adapter-guard']

    def run_cli(self, options, success=True):
        result = subprocess.run(['bash', str(self.repo / 'bin/Meta2Data-AmpliconPIP'), *options],
                                cwd=self.work, env=self.env, text=True, capture_output=True, timeout=30)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_parallel_workers_and_state_receive_own_primers(self):
        self.run_cli([*self.base, '--public-bioprojectIDs', 'PRJNA2', 'PRJNA1',
                      '--public-primer-fwd', 'AAAA', 'CCCC', '--public-primer-rev', 'TTTT', 'GGGG'])
        for project, forward, reverse in [('PRJNA1', 'CCCC', 'GGGG'), ('PRJNA2', 'AAAA', 'TTTT')]:
            worker = json.loads((self.output / project / 'worker_primers.json').read_text())
            state = json.loads((self.output / project / f'{project}-vsearch-run.json').read_text())
            self.assertEqual(worker, {'forward': forward, 'reverse': reverse})
            self.assertEqual((state['primer_fwd'], state['primer_rev']), (forward, reverse))
        self.assertFalse((self.output / 'PRJNA9').exists())
        self.assertEqual(len((self.output / 'PRJNA1/PRJNA1_sra.txt').read_text().splitlines()), 3)

    def test_selection_before_test_subset_and_no_reverse_autofill(self):
        self.run_cli([*self.base, '--test', '--public-bioprojectIDs', 'PRJNA1', '--public-primer-fwd', 'ACGT'])
        with (self.output / 'selected_metadata.csv').open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 3)
        with (self.output / 'test_subset_meta.csv').open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 2)
        worker = json.loads((self.output / 'PRJNA1/worker_primers.json').read_text())
        self.assertEqual(worker, {'forward': 'ACGT', 'reverse': ''})

    def test_batch_keeps_automatic_primers(self):
        self.run_cli(self.base)
        self.assertFalse((self.output / 'project_primers.json').exists())
        for project in ('PRJNA1', 'PRJNA2', 'PRJNA9'):
            worker = json.loads((self.output / project / 'worker_primers.json').read_text())
            self.assertEqual(worker, {'forward': '', 'reverse': ''})

    def test_invalid_scopes_fail_before_launch(self):
        for extra, message in [(['--public-primer-fwd', 'ACGT'], 'require --public-bioprojectIDs'),
                               (['--public-bioprojectIDs', 'PRJNA1', 'PRJNA2', '--public-primer-fwd', 'ACGT'], 'one primer per'),
                               (['--local', '--public-bioprojectIDs', 'PRJNA1', '--public-primer-fwd', 'ACGT'], 'Use --local-m local_metadata.csv'),
                               (['--public-bioprojectIDs', 'PRJNA1', 'PRJNA1', '--public-primer-fwd', 'AC', 'TG'], 'duplicate'),
                               (['--public-bioprojectIDs', 'PRJNA404', '--public-primer-fwd', 'ACGT'], 'not present')]:
            with self.subTest(extra=extra):
                result = self.run_cli([*self.base, *extra], success=False)
                self.assertIn(message, result.stdout + result.stderr)
                self.assertFalse((self.output / 'datasets_ID.txt').exists())

    def test_selection_and_test_cannot_overwrite_original_metadata(self):
        self.output.mkdir(parents=True)
        source = self.output / 'test_subset_meta.csv'
        original = self.metadata.read_bytes()
        source.write_bytes(original)
        options = [str(source) if argument == str(self.metadata) else argument for argument in self.base]
        result = self.run_cli([*options, '--test', '--public-bioprojectIDs', 'PRJNA1', '--public-primer-fwd', 'ACGT'], success=False)
        self.assertIn('would be overwritten', result.stderr)
        self.assertEqual(source.read_bytes(), original)
        self.assertFalse((self.output / 'datasets_ID.txt').exists())

    def test_removed_output_options_fail_without_writes(self):
        requested = self.work / 'requested_output'
        for option in ('-o', '--output'):
            for options in ([option], [*self.base, option, str(requested)]):
                with self.subTest(options=options):
                    result = self.run_cli(options, success=False)
                    self.assertIn(f'{option} has been removed', result.stderr)
                    self.assertIn('results/pip/', result.stderr)
                    self.assertFalse(requested.exists())
                    self.assertFalse((self.work / 'results').exists())
        result = self.run_cli([*self.base, f'--output={requested}'], success=False)
        self.assertIn('has been removed', result.stderr)
        self.assertFalse(requested.exists())
        self.assertFalse((self.work / 'results').exists())

    def test_renamed_options_explain_migration_without_writes(self):
        replacements = {
            '-m': '--public-m', '--metadata': '--public-m',
            '--col-bioproject': '--public-bioproject-colNAME',
            '--col-sra': '--public-sra-colNAME',
            '--bioprojectIDs': '--public-bioprojectIDs',
            '--primer-fwd': '--public-primer-fwd',
            '--primer-rev': '--public-primer-rev',
            '--local-col-datasets': '--local-datasets-colNAME',
            '--local-col-path': '--local-path-colNAME',
            '--local-col-platform': '--local-platform-colNAME',
            '--local-col-primer-f': '--local-primer-f-colNAME',
            '--local-col-primer-r': '--local-primer-r-colNAME',
        }
        for old, new in replacements.items():
            for options in ([old], [*self.base, old, 'value'], [f'{old}=value']):
                with self.subTest(options=options):
                    result = self.run_cli(options, success=False)
                    self.assertIn(f'{old} has been removed. Use {new} instead', result.stderr)
                    self.assertFalse((self.work / 'results').exists())

    def test_value_options_fail_before_configuration_or_output_when_value_missing(self):
        self.configure_reference_stubs()
        scalar_options = [
            '--public-m', '--public-bioproject-colNAME', '--public-sra-colNAME',
            '--local-m', '--local-datasets-colNAME', '--local-path-colNAME',
            '--local-platform-colNAME', '--local-primer-f-colNAME', '--local-primer-r-colNAME',
            '-t', '--threads', '--max-parallel', '--parameter',
        ]
        list_options = ['--public-bioprojectIDs', '--public-primer-fwd', '--public-primer-rev']
        for option in scalar_options + list_options:
            expected = 'requires a value' if option in scalar_options else 'requires one or more values'
            for suffix in ([], [''], ['--vsearch']):
                with self.subTest(option=option, suffix=suffix):
                    result = self.run_cli([option, *suffix], success=False)
                    self.assertIn(f'{option} {expected}', result.stderr)
                    self.assertFalse((self.work / 'results').exists())
                    self.assertFalse((self.work / 'reference_called').exists())

    def test_removed_db_options_explain_migration(self):
        for option in ('--db', '--dl'):
            result = self.run_cli([option], success=False)
            self.assertIn('removed', result.stderr)

    def test_explicit_trimming_never_calls_automatic_detection(self):
        runner = (ROOT / 'scripts/run.sh').read_text()
        function = '_trim_primers() {' + runner.split('_trim_primers() {', 1)[1].split('_local_register_dataset()', 1)[0]
        (self.scripts / 'explicit_primers.py').write_text("""import json,os,sys
from pathlib import Path
args=sys.argv[1:]; output=Path(args[args.index('--output')+1])
(output/'explicit.json').write_text(json.dumps({'forward':args[args.index('--forward')+1], 'reverse':args[args.index('--reverse')+1]}))
sys.exit(int(os.environ.get('TEST_PRIMER_STATUS','0')))
""")
        (self.scripts / 'entropy_primer_detect.py').write_text("raise RuntimeError('automatic detection must not run')\n")
        for status in (0, 7):
            folder = self.work / f'trim-{status}'; folder.mkdir()
            env = dict(self.env, SCRIPTS=str(self.scripts), PRIMER_FWD='ACGT', PRIMER_REV='',
                       dataset_path=str(folder), dataset_ID='PRJNA1', MODE='vsearch',
                       READ_COUNTS_REPORTS=str(folder), TEST_PRIMER_STATUS=str(status))
            script = function + '\n_trim_primers "$dataset_path/input" "$dataset_path/trimmed" --detect-only\n'
            result = subprocess.run(['bash', '-c', script], env=env, text=True, capture_output=True, timeout=10)
            self.assertEqual(result.returncode, status, result.stdout + result.stderr)
            record = json.loads((folder / 'trimmed/explicit.json').read_text())
            self.assertEqual(record, {'forward': 'ACGT', 'reverse': ''})
            self.assertNotIn('automatic detection must not run', result.stderr)

    def configure_reference_stubs(self):
        (self.scripts / 'reference_resources.py').write_text("""import sys
from pathlib import Path
args=sys.argv[1:]
assert args[0]=='prepare'
work=Path(args[args.index('--work-dir')+1])
assert args[args.index('--kind')+1]=='gg2-sequences'
(work/'reference_called').write_text(str(work))
ref=work/'results/db/2024.09.backbone.full-length.fna.qza'
ref.parent.mkdir(parents=True,exist_ok=True); ref.write_text('stub-reference')
print(ref)
""")
        (self.scripts / 'fastp_guard.py').write_text("""import json,sys
from pathlib import Path
args=sys.argv[1:]
cache=Path(args[args.index('--cache-dir')+1]); cache.mkdir(parents=True,exist_ok=True)
manifest=cache/'manifest.json'
manifest.write_text(json.dumps({'reference':args[args.index('--reference')+1], 'reference_sha256':'a'*64, 'blast_version':'blastn: test', 'schema_version':1, 'database_prefix':str(cache/'reference')}))
print(manifest)
""")
        (self.scripts / 'run.sh').write_text('exit 0\n')

    def test_default_guard_uses_fixed_gg2_and_cache_for_all_input_modes(self):
        self.configure_reference_stubs()
        self.local_folder('dataset_A', 'sample_A')
        local_metadata = self.local_csv()
        # Old environment variable names cannot redirect the fixed resources.
        inherited_cache = self.work / 'custom_cache'
        self.env.update(ADAPTER_REF=str(self.work / 'custom.fasta'), ADAPTER_CACHE=str(inherited_cache))
        cases = [
            ('public', self.base),
            ('local', self.local_options(local_metadata)),
            ('mixed', [*self.base, '--local-m', str(local_metadata), *self.local_columns]),
        ]
        for mode, options in cases:
            with self.subTest(mode=mode):
                called = self.work / 'reference_called'
                if called.exists():
                    called.unlink()
                self.run_cli([argument for argument in options if argument != '--no-adapter-guard'])
                self.assertEqual(called.read_text(), str(self.work))
                manifest = json.loads((self.output / 'adapter_guard_reference.json').read_text())
                self.assertEqual(manifest['reference'], str(self.work / 'results/db/2024.09.backbone.full-length.fna.qza'))
                self.assertEqual(manifest['database_prefix'], str(self.work / 'results/db/adapter_guard/reference'))
                self.assertEqual(json.loads((self.work / 'results/db/adapter_guard/manifest.json').read_text()), manifest)
                self.assertFalse(inherited_cache.exists())
                self.assertFalse((self.output / 'db').exists())
                if mode != 'public':
                    local_manifest = json.loads((self.output / 'local_datasets.json').read_text())
                    self.assertEqual(local_manifest['dataset_order'], ['dataset_A'])
                    self.assertTrue((self.output / 'dataset_A/.local-source.json').is_file())

    def test_guard_config_and_project_primers_reach_actual_worker_state(self):
        self.configure_reference_stubs()
        shutil.copy(ROOT / 'scripts/run.sh', self.scripts / 'run.sh')
        configuration = self.work / 'parameters.json'
        configuration.write_text(json.dumps({'adapter_guard': {'enabled': False, 'min_identity': 99.2}}))
        self.run_cli([*self.base, '--parameter', str(configuration), '--adapter-guard',
                      '--public-bioprojectIDs', 'PRJNA2', 'PRJNA1', '--public-primer-fwd', 'AAAA', 'CCCC'])
        for project, forward in [('PRJNA1', 'CCCC'), ('PRJNA2', 'AAAA')]:
            state = json.loads((self.output / project / f'{project}-vsearch-run.json').read_text())
            self.assertEqual(state['primer_fwd'], forward)
            self.assertTrue(state['adapter_guard']['enabled'])
            self.assertEqual(state['adapter_guard']['reference_sha256'], 'a' * 64)
            self.assertEqual(state['adapter_guard']['min_identity'], 99.2)
        effective = json.loads((self.output / 'effective-parameters-vsearch.json').read_text())
        self.assertTrue(effective['adapter_guard']['enabled'])

    def test_removed_adapter_reference_and_cache_options_fail_without_writes(self):
        self.configure_reference_stubs()
        custom = self.work / 'custom.fasta'
        original = '>reference\nACGT\n'
        custom.write_text(original)
        requested_cache = self.work / 'requested_cache'
        for option, value in [('--adapter-ref', custom), ('--adapter-cache', requested_cache)]:
            for options in ([option], [option, str(value)], [f'{option}={value}'],
                            [*self.base, '--adapter-guard', option, str(value)],
                            [*self.base, f'{option}={value}']):
                with self.subTest(options=options):
                    result = self.run_cli(options, success=False)
                    self.assertIn(f'{option} has been removed', result.stderr)
                    self.assertIn('automatically uses GG2', result.stderr)
                    self.assertIn('fixed at results/db/adapter_guard/', result.stderr)
                    self.assertFalse((self.work / 'results').exists())
                    self.assertFalse((self.work / 'reference_called').exists())
                    self.assertFalse(requested_cache.exists())
                    self.assertEqual(custom.read_text(), original)
        help_output = self.run_cli(['--help']).stdout
        self.assertNotIn('--adapter-ref', help_output)
        self.assertNotIn('--adapter-cache', help_output)
        self.assertIn('--adapter-guard', help_output)
        self.assertIn('--no-adapter-guard', help_output)

    def local_folder(self, name, sample):
        folder = self.work / 'local_data' / name
        folder.mkdir(parents=True)
        (folder / f'{sample}.fastq').write_text('@r\nACGT\n+\nIIII\n')
        return folder

    def local_csv(self, rows=None, filename='local_metadata.csv', fields=None):
        if rows is None:
            rows = [{'datasets': folder.name, 'path': str(folder), 'platform': 'LS454'}
                    for folder in sorted((self.work / 'local_data').iterdir()) if folder.is_dir()]
        destination = self.work / filename
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
        return destination

    def local_options(self, metadata=None):
        # LS454/ONT with dada2 deliberately skips unsupported analysis after real
        # worker initialization. No sequencing tools or network are used here.
        return ['--local-m', str(metadata or self.local_csv()), *self.local_columns, '--dada2', '-t', '8', '--no-adapter-guard']

    def test_local_csv_only_automatic_rows_never_start_online_lookup(self):
        self.local_folder('dataset_A', 'sample_A')
        self.local_folder('dataset_B', 'sample_B')
        result = self.run_cli(self.local_options())
        output = self.work / 'results/pip'
        self.assertIn('2 parallel datasets, 4 threads per dataset', result.stdout)
        for name in ('dataset_A', 'dataset_B'):
            self.assertEqual(json.loads((output / name / 'worker_primers.json').read_text()),
                             {'forward': '', 'reverse': ''})
            self.assertEqual(json.loads((output / name / 'worker_context.json').read_text()),
                             {'local': '1', 'platform': 'LS454'})
            self.assertEqual((output / name / 'worker_threads.txt').read_text(), '4')
        self.assertFalse((output / 'platform_queries.txt').exists())
        self.assertEqual((output / 'datasets_ID.txt').read_text().splitlines(), ['dataset_A', 'dataset_B'])

    def test_local_csv_rows_control_names_platforms_and_optional_primers(self):
        a = self.local_folder('folder_A', 'sample_A')
        b = self.local_folder('folder_B', 'sample_B')
        c = self.local_folder('folder_C', 'sample_C')
        metadata = self.local_csv([
            {'datasets': 'named_B', 'path': str(b), 'platform': 'OXFORD_NANOPORE', 'fwd': 'AAAA', 'rev': 'TTTT'},
            {'datasets': 'named_A', 'path': str(a), 'platform': 'LS454', 'fwd': '', 'rev': ''},
            {'datasets': 'named_C', 'path': str(c), 'platform': 'LS454', 'fwd': 'CCCC', 'rev': ''},
        ])
        original = metadata.read_bytes()
        result = self.run_cli([*self.local_options(metadata), '--max-parallel', '6',
                              '--local-primer-f-colNAME', 'fwd', '--local-primer-r-colNAME', 'rev'])
        output = self.work / 'results/pip'
        self.assertIn('3 parallel datasets, 2 threads per dataset', result.stdout)
        self.assertEqual((output / 'datasets_ID.txt').read_text().splitlines(), ['named_B', 'named_A', 'named_C'])
        for name, platform, forward, reverse, source in [
                ('named_B', 'OXFORD_NANOPORE', 'AAAA', 'TTTT', b),
                ('named_A', 'LS454', '', '', a), ('named_C', 'LS454', 'CCCC', '', c)]:
            self.assertEqual(json.loads((output / name / 'worker_primers.json').read_text()),
                             {'forward': forward, 'reverse': reverse})
            self.assertEqual(json.loads((output / name / 'worker_context.json').read_text()),
                             {'local': '1', 'platform': platform})
            state = json.loads((output / name / f'{name}-dada2-run.json').read_text())
            self.assertEqual((state['primer_fwd'], state['primer_rev'], state['platform']), (forward, reverse, platform))
            self.assertEqual(Path(state['local_source'][0][0]).parent, source)
            self.assertEqual(state['source_kind'], 'local')
        self.assertEqual(metadata.read_bytes(), original)
        self.assertFalse((output / 'folder_A').exists())

    def test_mixed_queue_keeps_online_and_local_source_context_separate(self):
        a = self.local_folder('folder_A', 'sample_A')
        b = self.local_folder('folder_B', 'sample_B')
        metadata = self.local_csv([
            {'datasets': 'local_B', 'path': str(b), 'platform': 'OXFORD_NANOPORE', 'fwd': 'AAAA', 'rev': 'TTTT'},
            {'datasets': 'local_A', 'path': str(a), 'platform': 'LS454', 'fwd': '', 'rev': ''},
        ])
        online = [argument for argument in self.base if argument != '--vsearch']
        result = self.run_cli([*online, '--dada2', '--local-m', str(metadata), *self.local_columns,
                              '--local-primer-f-colNAME', 'fwd', '--local-primer-r-colNAME', 'rev',
                              '-t', '10', '--max-parallel', '6'])
        self.assertIn('5 parallel datasets, 2 threads per dataset', result.stdout)
        self.assertEqual((self.output / 'datasets_ID.txt').read_text().splitlines(),
                         ['PRJNA1', 'PRJNA2', 'PRJNA9', 'local_B', 'local_A'])
        queries = (self.output / 'platform_queries.txt').read_text().splitlines()
        self.assertEqual([line.split('\t')[0] for line in queries], ['PRJNA1', 'PRJNA2', 'PRJNA9'])
        for name in ('PRJNA1', 'PRJNA2', 'PRJNA9'):
            self.assertEqual(json.loads((self.output / name / 'worker_context.json').read_text()),
                             {'local': '0', 'platform': ''})
            self.assertEqual(json.loads((self.output / name / 'worker_primers.json').read_text()),
                             {'forward': '', 'reverse': ''})
            state = json.loads((self.output / name / f'{name}-dada2-run.json').read_text())
            self.assertEqual(state['source_kind'], 'archive')
            self.assertEqual(state['local_source'], [])
        self.assertEqual(json.loads((self.output / 'local_B/worker_primers.json').read_text()),
                         {'forward': 'AAAA', 'reverse': 'TTTT'})
        self.assertEqual(json.loads((self.output / 'local_A/worker_context.json').read_text()),
                         {'local': '1', 'platform': 'LS454'})
        self.assertEqual((self.output / 'datasets.log').read_text().count('# === RUN'), 1)

    def test_local_custom_columns_and_csv_relative_paths(self):
        folder = self.local_folder('folder A', 'sample_A')
        metadata = self.local_csv([{'ID': 'named_A', 'Directory': '../local_data/folder A',
                                    'Technology': 'LS454', 'F': 'ACGT'}], filename='tables/local.csv')
        result = self.run_cli([*self.local_options(metadata), '--local-datasets-colNAME', 'ID',
                              '--local-path-colNAME', 'Directory', '--local-platform-colNAME', 'Technology',
                              '--local-primer-f-colNAME', 'F'])
        output = self.work / 'results/pip/named_A'
        self.assertIn('1 parallel datasets, 8 threads per dataset', result.stdout)
        self.assertEqual(json.loads((output / 'worker_primers.json').read_text()), {'forward': 'ACGT', 'reverse': ''})
        state = json.loads((output / 'named_A-dada2-run.json').read_text())
        self.assertEqual(Path(state['local_source'][0][0]).parent, folder)

    def test_unmapped_primer_columns_are_not_silently_enabled(self):
        folder = self.local_folder('folder_A', 'sample_A')
        metadata = self.local_csv([{'datasets': 'named_A', 'path': str(folder), 'platform': 'LS454', 'fwd': 'ACGT'}])
        self.run_cli(self.local_options(metadata))
        record = json.loads((self.work / 'results/pip/named_A/worker_primers.json').read_text())
        self.assertEqual(record, {'forward': '', 'reverse': ''})

    def test_legacy_local_options_give_csv_migration_message(self):
        for option in ('--local', '--localinput', '--input', '--localIDs', '--platform'):
            with self.subTest(option=option):
                result = self.run_cli([option], success=False)
                self.assertIn('has been removed', result.stderr)
                self.assertIn('Use --local-m local_metadata.csv', result.stderr)

    def test_local_csv_argument_scope_validation(self):
        for options, message in [(['--local-m'], '--local-m requires a value'),
                                 (['--local-path-colNAME', 'directory'], 'require --local-m'),
                                 (['--local-m', 'local.csv', '--local-primer-r-colNAME', 'R'], 'requires --local-primer-f-colNAME'),
                                 (['--local-m', 'local.csv', '--test'], 'cannot be combined'),
                                 (['--local-m', 'local.csv', '--public-primer-fwd', 'ACGT'], 'local primers must be mapped'),
                                 (['--local-m', 'local.csv', '--public-bioprojectIDs', 'PRJNA1'], 'cannot be combined'),
                                 (['--local-m', 'local.csv', '--public-bioproject-colNAME', 'Project'], 'require online')]:
            with self.subTest(options=options):
                result = self.run_cli(options, success=False)
                self.assertIn(message, result.stderr)

    def test_all_local_column_mappings_are_required_before_any_writes(self):
        self.local_folder('dataset_A', 'sample_A')
        metadata = self.local_csv()
        before = metadata.read_bytes()
        self.configure_reference_stubs()
        # Standard column names in the CSV do not make the CLI flags optional.
        for online in (False, True):
            cases = [([], '--local-datasets-colNAME')]
            for index in range(0, len(self.local_columns), 2):
                cases.append((self.local_columns[:index] + self.local_columns[index + 2:],
                              self.local_columns[index]))
            for columns, missing in cases:
                with self.subTest(online=online, missing=missing):
                    options = [argument for argument in self.base if argument != '--no-adapter-guard'] if online else ['--dada2']
                    result = self.run_cli([*options, '--local-m', str(metadata), *columns], success=False)
                    self.assertIn(f'--local-m requires {missing} NAME', result.stderr)
                    self.assertFalse(self.output.exists())
                    self.assertFalse((self.work / 'results').exists())
                    self.assertFalse((self.work / 'reference_called').exists())
                    self.assertEqual(metadata.read_bytes(), before)

    def test_invalid_local_csv_fails_before_guard_or_online_registration(self):
        folder = self.local_folder('folder_A', 'sample_A')
        self.configure_reference_stubs()
        rows = [
            {'datasets': 'named_A', 'path': str(folder), 'platform': 'INVALID', 'fwd': '', 'rev': ''},
            {'datasets': 'named_A', 'path': 'missing_folder', 'platform': 'LS454', 'fwd': '', 'rev': ''},
            {'datasets': 'named_A', 'path': str(folder), 'platform': 'LS454', 'fwd': '', 'rev': 'TTTT'},
        ]
        for row in rows:
            metadata = self.local_csv([row])
            options = [argument for argument in self.base if argument != '--no-adapter-guard']
            self.run_cli([*options, '--local-m', str(metadata), *self.local_columns, '--local-primer-f-colNAME', 'fwd',
                          '--local-primer-r-colNAME', 'rev'], success=False)
            self.assertFalse((self.work / 'reference_called').exists())
            self.assertFalse((self.output / 'datasets_ID.txt').exists())
            self.assertFalse((self.output / 'PRJNA1').exists())
            self.assertFalse((self.output / 'local_datasets.json').exists())

    def test_mixed_same_batch_name_and_sample_collisions_fail_before_writes(self):
        folder = self.local_folder('folder_A', 'PRJNA1_SRR1')
        for name in ('PRJNA1', 'different_dataset'):
            metadata = self.local_csv([{'datasets': name, 'path': str(folder), 'platform': 'LS454'}])
            self.run_cli([*self.base, '--local-m', str(metadata), *self.local_columns], success=False)
            self.assertFalse((self.output / 'PRJNA1').exists())
            self.assertFalse((self.output / 'local_datasets.json').exists())
            self.assertFalse((self.output / 'datasets_ID.txt').exists())

    def test_mixed_custom_online_run_column_is_validated_before_local_markers(self):
        self.local_folder('folder_A', 'sample_A')
        local = self.local_csv()
        self.metadata.write_text('Project,Accession\nPRJNA1,SRR1\n')
        options = ['--public-m', str(self.metadata), '--public-bioproject-colNAME', 'Project', '--public-sra-colNAME', 'Missing',
                   '--local-m', str(local), *self.local_columns, '--vsearch', '--no-adapter-guard']
        self.run_cli(options, success=False)
        self.assertFalse((self.work / 'results/pip/folder_A').exists())
        options[options.index('Missing')] = 'Accession'
        (self.scripts / 'run.sh').write_text('exit 0\n')
        self.run_cli(options)
        self.assertTrue((self.work / 'results/pip/folder_A/.local-source.json').exists())

    def test_online_default_output_is_startup_workdir(self):
        metadata_dir = self.work / 'other'; metadata_dir.mkdir()
        source = metadata_dir / 'metadata.csv'; source.write_bytes(self.metadata.read_bytes())
        self.run_cli(['--public-m', str(source), '--public-bioproject-colNAME', 'Bioproject', '--public-sra-colNAME', 'Run',
                      '--vsearch', '--no-adapter-guard'])
        self.assertTrue((self.work / 'results/pip/PRJNA1/worker_primers.json').exists())
        self.assertFalse((metadata_dir / 'PRJNA1').exists())

    def test_builtin_test_retains_unprefixed_flag_and_online_defaults(self):
        (self.repo / 'test').mkdir()
        shutil.copy(self.metadata, self.repo / 'test/ampliconpiptest.csv')
        self.run_cli(['--test', '--vsearch', '--no-adapter-guard'])
        self.assertTrue((self.output / 'test_subset_meta.csv').is_file())
        queries = (self.output / 'platform_queries.txt').read_text().splitlines()
        self.assertEqual([line.split('\t')[0] for line in queries], ['PRJNA1', 'PRJNA2', 'PRJNA9'])
        self.assertEqual(len((self.output / 'PRJNA1/PRJNA1_sra.txt').read_text().splitlines()), 2)
        self.assertFalse((self.output / 'local_datasets.json').exists())

    def test_test_mode_default_output_matches_other_modes(self):
        self.run_cli(['--public-m', str(self.metadata), '--public-bioproject-colNAME', 'Bioproject', '--public-sra-colNAME', 'Run',
                      '--test', '--vsearch', '--no-adapter-guard'])
        self.assertTrue((self.work / 'results/pip/test_subset_meta.csv').is_file())
        self.assertFalse((self.work / 'test_subset_meta.csv').exists())

    def test_online_dataset_cannot_be_overwritten_by_same_named_local_dataset(self):
        self.run_cli(self.base)
        state = self.output / 'PRJNA1/PRJNA1-vsearch-run.json'
        original = state.read_bytes()
        folder = self.local_folder('folder_A', 'local_sample')
        local = self.local_csv([{'datasets': 'PRJNA1', 'path': str(folder), 'platform': 'LS454'}])
        self.run_cli(self.local_options(local), success=False)
        self.assertEqual(state.read_bytes(), original)

    def test_local_dataset_cannot_be_overwritten_by_same_named_online_project(self):
        self.local_folder('PRJNA1', 'local_sample')
        self.run_cli(self.local_options())
        state = self.output / 'PRJNA1/PRJNA1-dada2-run.json'
        original = state.read_bytes()
        self.run_cli(self.base, success=False)
        self.assertEqual(state.read_bytes(), original)

    def test_default_output_online_local_online_keeps_all_dataset_state(self):
        online = ['--public-m', str(self.metadata), '--public-bioproject-colNAME', 'Bioproject', '--public-sra-colNAME', 'Run',
                  '--vsearch', '--no-adapter-guard']
        self.run_cli(online)
        output = self.work / 'results/pip'
        original_online = (output / 'PRJNA1/PRJNA1-vsearch-run.json').read_bytes()
        self.local_folder('dataset_A', 'sample_A')
        self.run_cli(self.local_options())
        local_state = output / 'dataset_A/dataset_A-dada2-run.json'
        original_local = local_state.read_bytes()
        self.assertEqual((output / 'PRJNA1/PRJNA1-vsearch-run.json').read_bytes(), original_online)
        self.run_cli(online)
        self.assertEqual(local_state.read_bytes(), original_local)
        self.assertTrue((output / 'dataset_A/.local-source.json').exists())
        self.assertTrue((output / 'PRJNA2/PRJNA2-vsearch-run.json').exists())
        self.assertEqual((output / 'datasets.log').read_text().count('# === RUN'), 3)

    def test_output_lock_blocks_writes_and_allows_run_after_release(self):
        self.output.mkdir(parents=True)
        artifacts = {'selected_metadata.csv': 'unchanged input subset',
                     'project_primers.json': 'unchanged project mapping',
                     'effective-parameters-vsearch.json': 'unchanged parameters',
                     'local_datasets.json': 'unchanged local mapping'}
        for filename, content in artifacts.items():
            (self.output / filename).write_text(content)
        with (self.output / '.pip.lock').open('a') as handle:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            result = self.run_cli([*self.base, '--public-bioprojectIDs', 'PRJNA1', '--public-primer-fwd', 'ACGT'], success=False)
            self.assertIn('another AmpliconPIP command is using output directory', result.stderr)
            self.assertFalse((self.output / 'datasets_ID.txt').exists())
            for filename, content in artifacts.items():
                self.assertEqual((self.output / filename).read_text(), content)
        self.run_cli([*self.base, '--public-bioprojectIDs', 'PRJNA1', '--public-primer-fwd', 'ACGT'])
        self.assertTrue((self.output / 'PRJNA1/worker_primers.json').exists())

    def test_wrapper_holds_output_lock_after_python_helper_exits(self):
        ready = self.work / 'worker_ready'
        release = self.work / 'worker_release'
        self.env['TEST_WORKER_READY'] = str(ready)
        self.env['TEST_WORKER_RELEASE'] = str(release)
        (self.scripts / 'run.sh').write_text("""python3 - <<'WAIT_PY'
import os, time
from pathlib import Path
Path(os.environ['TEST_WORKER_READY']).write_text('ready')
for attempt in range(200):
    if Path(os.environ['TEST_WORKER_RELEASE']).exists():
        break
    time.sleep(0.05)
else:
    raise SystemExit('test worker timed out')
WAIT_PY
""")
        command = ['bash', str(self.repo / 'bin/Meta2Data-AmpliconPIP'), *self.base]
        worker = subprocess.Popen(command, cwd=self.work, env=self.env, text=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            for attempt in range(500):
                if ready.exists() or worker.poll() is not None:
                    break
                time.sleep(0.01)
            self.assertTrue(ready.exists(), 'first wrapper did not reach the running worker')
            result = self.run_cli(self.base, success=False)
            self.assertIn('another AmpliconPIP command is using output directory', result.stderr)
            release.write_text('finish')
            stdout, stderr = worker.communicate(timeout=10)
            self.assertEqual(worker.returncode, 0, stdout + stderr)
            self.run_cli(self.base)
        finally:
            release.write_text('finish')
            if worker.poll() is None:
                worker.terminate()
            worker.communicate(timeout=10)



if __name__ == '__main__':
    unittest.main()
