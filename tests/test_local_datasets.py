"""CSV dataset planning validates all rows before touching prior results."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import local_datasets as local


class LocalDatasetTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.input = self.work / 'local_data'
        self.output = self.work / 'results/pip'
        self.fastq(self.input / 'dataset_A', 'A01')
        self.fastq(self.input / 'dataset_B', 'B01')
        self.metadata = self.work / 'local_metadata.csv'
        self.write_csv([['dataset_A', 'local_data/dataset_A', 'ILLUMINA', '', ''],
                        ['dataset_B', 'local_data/dataset_B', 'ILLUMINA', '', '']])

    def fastq(self, folder, sample):
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / (sample + '.fastq')
        path.write_text('@read\nACGT\n+\nIIII\n')
        return path

    def write_csv(self, rows, headers=None, path=None):
        path = path or self.metadata
        with path.open('w', encoding='utf-8-sig', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(headers or ['datasets', 'path', 'platform', 'primer f', 'primer r'])
            writer.writerows(rows)
        return path

    def single(self, source, platform='ILLUMINA', name='dataset_A', forward='', reverse=''):
        return self.write_csv([[name, str(source), platform, forward, reverse]])

    def prepare(self, **kwargs):
        return local.prepare(self.metadata, self.output, **kwargs)

    def test_csv_order_explicit_names_and_sources_preserve_inputs(self):
        self.fastq(self.input / 'nested/deeper', 'not_selected')
        self.write_csv([['named_B', 'local_data/dataset_B', 'LS454'],
                        ['named_A', 'local_data/dataset_A', 'ILLUMINA']])
        before = {p: p.read_bytes() for p in self.input.rglob('*.fastq')}
        record = local.read_manifest(self.prepare())
        self.assertEqual(record['dataset_order'], ['named_B', 'named_A'])
        self.assertEqual(record['datasets']['named_B']['platform'], 'LS454')
        self.assertEqual(record['datasets']['named_A']['forward'], '')
        self.assertEqual(Path(record['datasets']['named_A']['source']), self.input / 'dataset_A')
        self.assertTrue((self.output / 'named_A/.local-source.json').is_file())
        self.assertEqual(before, {p: p.read_bytes() for p in self.input.rglob('*.fastq')})

    def test_explicit_primers_and_empty_rows_use_their_own_platform(self):
        self.write_csv([['B', 'local_data/dataset_B', 'pacbio_smrt', ' cctn ', 'GGRY'],
                        ['A', 'local_data/dataset_A', 'ILLUMINA', '', '']])
        record = local.read_manifest(self.prepare(col_primer_f='primer f', col_primer_r='primer r'))
        self.assertEqual(record['datasets']['B']['forward'], 'CCTN')
        self.assertEqual(record['datasets']['B']['reverse'], 'GGRY')
        self.assertEqual(record['datasets']['B']['platform'], 'PACBIO_SMRT')
        self.assertEqual(record['datasets']['A']['forward'], '')
        self.assertEqual(record['datasets']['A']['reverse'], '')

    def test_unmapped_primer_columns_are_ignored(self):
        self.single(self.input / 'dataset_A', forward='primer-name', reverse='bad path')
        record = local.read_manifest(self.prepare())
        self.assertEqual(record['datasets']['dataset_A']['forward'], '')
        self.assertEqual(record['datasets']['dataset_A']['reverse'], '')

    def test_forward_only_mapping_does_not_read_unmapped_reverse(self):
        self.single(self.input / 'dataset_A', forward='acgtn', reverse='ignore me')
        record = local.read_manifest(self.prepare(col_primer_f='primer f'))
        self.assertEqual(record['datasets']['dataset_A']['forward'], 'ACGTN')
        self.assertEqual(record['datasets']['dataset_A']['reverse'], '')

    def test_invalid_row_values_do_not_register_any_dataset(self):
        bad_rows = [['', 'local_data/dataset_B', 'ILLUMINA', '', ''],
                    ['B', '', 'ILLUMINA', '', ''], ['B', 'local_data/dataset_B', '', '', ''],
                    ['../B', 'local_data/dataset_B', 'ILLUMINA', '', ''],
                    ['dataset_A', 'local_data/dataset_B', 'ILLUMINA', '', ''],
                    ['B', 'missing', 'ILLUMINA', '', ''],
                    ['B', 'local_data/dataset_B', 'ILLUMINA', 'primer.fasta', ''],
                    ['B', 'local_data/dataset_B', 'ILLUMINA', '', 'ACGT']]
        for row in bad_rows:
            with self.subTest(row=row):
                self.write_csv([['dataset_A', 'local_data/dataset_A', 'ILLUMINA', 'ACGT', 'TGCA'], row])
                with self.assertRaises(ValueError):
                    self.prepare(col_primer_f='primer f', col_primer_r='primer r')
                self.assertFalse(self.output.exists())

    def test_invalid_column_mappings_and_headers_are_rejected(self):
        for options in [dict(col_primer_r='primer r'), dict(col_primer_f='missing'),
                        dict(col_path='datasets'), dict(col_platform='missing')]:
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.prepare(**options)
            self.assertFalse(self.output.exists())
        self.write_csv([['A', 'local_data/dataset_A', 'ILLUMINA']], ['datasets', 'path', 'path'])
        with self.assertRaisesRegex(ValueError, 'duplicate column names'):
            self.prepare()

    def test_custom_headers_paths_relative_to_csv_and_quoted_fields(self):
        folder = self.work / 'metadata folder'
        folder.mkdir()
        metadata = self.write_csv([['A', '../local_data/dataset_A', 'ILLUMINA', 'ACGN']],
                                  ['name', 'path to datasets', 'technology', 'fwd'],
                                  folder / 'local.csv')
        record = local.read_manifest(local.prepare(metadata, self.output, col_datasets='name',
            col_path='path to datasets', col_platform='technology', col_primer_f='fwd'))
        self.assertEqual(record['datasets']['A']['source'], str(self.input / 'dataset_A'))
        self.assertEqual(record['datasets']['A']['forward'], 'ACGN')
        comma_folder = self.work / 'reads, batch'
        self.fastq(comma_folder, 'C01')
        self.single(comma_folder, name='C')
        self.assertEqual(local.read_manifest(self.prepare())['datasets']['C']['source'], str(comma_folder))

    def test_row_lists_one_fastq_directory_without_recursion(self):
        self.fastq(self.input, 'at_root')
        self.single(self.input)
        record = local.read_manifest(self.prepare())
        self.assertEqual(record['datasets']['dataset_A']['samples'], ['at_root'])

    def test_empty_metadata_and_ragged_rows_are_rejected(self):
        self.write_csv([])
        with self.assertRaisesRegex(ValueError, 'contains no datasets'):
            self.prepare()
        self.write_csv([['A', 'local_data/dataset_A', 'ILLUMINA', '', '', 'extra']])
        with self.assertRaisesRegex(ValueError, 'more values'):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_mixed_online_names_samples_and_source_overlap_fail_preflight(self):
        online = self.work / 'online.csv'
        options = dict(online_metadata=online, col_bioproject='project', col_sra='run')
        online.write_text('project,run\ndataset_A,SRR1\n')
        with self.assertRaisesRegex(ValueError, 'same command'):
            self.prepare(**options)
        self.assertFalse(self.output.exists())
        online.write_text('project,run\nPRJNA1,SRR1\n')
        self.fastq(self.input / 'dataset_A', 'PRJNA1_SRR1')
        with self.assertRaisesRegex(ValueError, 'sample ID.*conflicts'):
            self.prepare(**options)
        self.assertFalse(self.output.exists())
        (self.input / 'dataset_A/PRJNA1_SRR1.fastq').unlink()
        online_source = self.output / 'PRJNA1'
        self.fastq(online_source, 'kept')
        self.single(online_source)
        with self.assertRaisesRegex(ValueError, 'overlaps'):
            self.prepare(**options)
        self.assertTrue((online_source / 'kept.fastq').exists())
        self.assertFalse((self.output / 'dataset_A').exists())

    def test_mixed_valid_inputs_and_bad_online_columns(self):
        online = self.work / 'online.csv'
        online.write_text('project,run\nPRJNA1,SRR1\n')
        with self.assertRaisesRegex(ValueError, 'missing column'):
            self.prepare(online_metadata=online, col_bioproject='project', col_sra='wrong')
        self.assertFalse(self.output.exists())
        self.prepare(online_metadata=online, col_bioproject='project', col_sra='run')
        self.assertFalse((self.output / 'PRJNA1').exists())

    def test_duplicate_sample_ids_fail_before_any_dataset_is_registered(self):
        self.fastq(self.input / 'dataset_B', 'A01')
        with self.assertRaisesRegex(ValueError, 'Duplicate sample ID.*A01'):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_invalid_layout_in_one_dataset_fails_entire_preflight(self):
        self.fastq(self.input / 'dataset_B', 'unpaired_R2')
        with self.assertRaises(ValueError):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_input_output_overlap_and_duplicate_sources_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'outside the input'):
            local.prepare(self.single(self.input / 'dataset_A'), self.input / 'dataset_A/results/pip')
        (self.input / 'alias').symlink_to(self.input / 'dataset_A', target_is_directory=True)
        self.write_csv([['A', 'local_data/dataset_A', 'ILLUMINA'],
                        ['alias', 'local_data/alias', 'ILLUMINA']])
        with self.assertRaisesRegex(ValueError, 'same input directory'):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_symlinked_fastq_inside_output_is_rejected(self):
        target = self.fastq(self.output / 'dataset_B', 'hidden_source')
        (self.input / 'dataset_A/linked.fastq').symlink_to(target)
        with self.assertRaisesRegex(ValueError, 'overlaps'):
            self.prepare()
        self.assertTrue(target.exists())
        self.assertFalse((self.output / 'dataset_A').exists())

    def test_same_source_reruns_but_another_source_cannot_reuse_name(self):
        self.prepare()
        result = self.output / 'dataset_A/dataset_A-dada2-final-table.qza'
        result.write_bytes(b'existing')
        self.fastq(self.input / 'dataset_A', 'A02')
        self.prepare()
        self.assertEqual(result.read_bytes(), b'existing')
        other = self.work / 'other/dataset_A'
        self.fastq(other, 'A99')
        with self.assertRaisesRegex(ValueError, 'different input directory'):
            local.prepare(self.single(other), self.output)
        self.assertEqual(result.read_bytes(), b'existing')

    def test_local_cannot_overwrite_online_or_unknown_existing_results(self):
        dataset = self.output / 'dataset_A'
        dataset.mkdir(parents=True)
        result = dataset / 'dataset_A-dada2-final-table.qza'
        result.write_bytes(b'existing')
        with self.assertRaisesRegex(ValueError, 'no verifiable source'):
            self.prepare()
        (dataset / '.checkpoint.json').write_text(json.dumps({'source_kind': 'archive'}))
        with self.assertRaisesRegex(ValueError, 'existing online project'):
            self.prepare()
        self.assertEqual(result.read_bytes(), b'existing')
        self.assertFalse((self.output / 'dataset_B').exists())

    def test_online_cannot_overwrite_local_but_same_online_project_can_repeat(self):
        self.prepare()
        metadata = self.work / 'metadata.csv'
        metadata.write_text('Bioproject,Run\ndataset_A,SRR1\n')
        with self.assertRaisesRegex(ValueError, 'already belongs to local'):
            local.check_online(metadata, 'Bioproject', self.output)
        metadata.write_text('Bioproject,Run\nPRJNA1,SRR1\n')
        online = self.output / 'PRJNA1'
        online.mkdir()
        (online / '.checkpoint.json').write_text(json.dumps({'source_kind': 'archive'}))
        local.check_online(metadata, 'Bioproject', self.output)

    def test_old_local_checkpoint_migrates_for_the_same_source(self):
        target = self.output / 'dataset_A'
        target.mkdir(parents=True)
        source = self.input / 'dataset_A/A01.fastq'
        (target / '.checkpoint.json').write_text(json.dumps({'source_kind': 'local', 'local_source': [[str(source), 16, 1]]}))
        self.prepare()
        self.assertTrue((target / '.local-source.json').is_file())

    def test_pre_schema_local_record_blocks_online_and_migrates_same_directory(self):
        target = self.output / 'dataset_A'
        target.mkdir(parents=True)
        source = self.input / 'dataset_A/A01.fastq'
        (target / '.checkpoint.json').write_text(json.dumps({'inputs': 'A01\n', 'local_source': [[str(source), 16, 1]]}))
        metadata = self.work / 'metadata.csv'
        metadata.write_text('Bioproject,Run\ndataset_A,SRR1\n')
        with self.assertRaisesRegex(ValueError, 'already belongs to local'):
            local.check_online(metadata, 'Bioproject', self.output)
        self.prepare()
        self.assertTrue((target / '.local-source.json').is_file())

    def test_shared_file_does_not_let_another_directory_claim_legacy_results(self):
        target = self.output / 'dataset_A'
        target.mkdir(parents=True)
        source = self.input / 'dataset_A/A01.fastq'
        (target / '.checkpoint.json').write_text(json.dumps({'source_kind': 'local', 'local_source': [[str(source), 16, 1]]}))
        other = self.work / 'different/dataset_A'
        self.fastq(other, 'A02')
        (other / 'A01.fastq').symlink_to(source)
        with self.assertRaisesRegex(ValueError, 'different input directory'):
            local.prepare(self.single(other), self.output)
        self.assertFalse((target / '.local-source.json').exists())

    def test_directory_marker_allows_same_source_with_symlinked_fastq(self):
        external = self.fastq(self.work / 'external', 'A99')
        source = self.input / 'dataset_A'
        (source / 'A99.fastq').symlink_to(external)
        self.prepare()
        target = self.output / 'dataset_A'
        old_files = [[str(p.resolve()), 16, 1] for p in source.glob('*.fastq')]
        (target / '.checkpoint.json').write_text(json.dumps({'source_kind': 'local', 'local_source': old_files}))
        self.prepare()
        self.assertEqual(json.loads((target / '.local-source.json').read_text())['source_directory'], str(source))

    def test_platform_and_tsv_path_validation(self):
        self.single(self.input / 'dataset_A', 'UNKNOWN')
        with self.assertRaisesRegex(ValueError, 'unsupported platform'):
            self.prepare()
        path = self.work / 'bad\tpath'
        self.fastq(path, 'sample')
        self.single(path)
        with self.assertRaisesRegex(ValueError, 'tabs or newlines'):
            self.prepare()
        paired = self.work / 'paired'
        self.fastq(paired, 'sample_R1')
        self.fastq(paired, 'sample_R2')
        self.single(paired, 'PACBIO_SMRT')
        with self.assertRaisesRegex(ValueError, 'paired-end data are unsupported'):
            self.prepare()

    def test_empty_input_has_clear_error(self):
        empty = self.work / 'empty'
        empty.mkdir()
        self.single(empty)
        with self.assertRaisesRegex(ValueError, 'No FASTQ files'):
            self.prepare()

    def test_metadata_in_generated_dataset_directory_is_not_overwritten(self):
        target = self.output / 'dataset_A'
        target.mkdir(parents=True)
        metadata = target / 'metadata.csv'
        self.write_csv([['dataset_A', str(self.input / 'dataset_A'), 'ILLUMINA']], path=metadata)
        before = metadata.read_bytes()
        with self.assertRaisesRegex(ValueError, 'metadata must be outside'):
            local.prepare(metadata, self.output)
        self.assertEqual(metadata.read_bytes(), before)
        self.assertFalse((target / '.local-source.json').exists())

    def test_both_metadata_files_are_kept_outside_generated_output(self):
        self.output.mkdir(parents=True)
        for filename in ['local_datasets.json', 'summary.csv', 'effective-parameters-dada2.json',
                         'datasets_ID.txt']:
            with self.subTest(filename=filename):
                online = self.output / filename
                online.write_text('project,run\nPRJNA1,SRR1\n')
                before = online.read_bytes()
                with self.assertRaisesRegex(ValueError, 'metadata must be outside'):
                    self.prepare(online_metadata=online, col_bioproject='project', col_sra='run')
                self.assertEqual(online.read_bytes(), before)
                self.assertFalse((self.output / 'dataset_A').exists())

    def test_online_preflight_uses_same_id_normalization_and_unambiguous_headers(self):
        online = self.work / 'online.csv'
        options = dict(online_metadata=online, col_bioproject='project', col_sra='run')
        self.fastq(self.input / 'dataset_A', 'PRJNA1_SRR1')
        online.write_text('project,run\nPRJNA1,SRR 1\n')
        with self.assertRaisesRegex(ValueError, 'sample ID.*conflicts'):
            self.prepare(**options)
        online.write_text('project,run,run\nPRJNA1,SRR1,SRR2\n')
        with self.assertRaisesRegex(ValueError, 'duplicate column names'):
            self.prepare(**options)
        online.write_text('project,run\nPRJNA1,SRR1,extra\n')
        with self.assertRaisesRegex(ValueError, 'more values'):
            self.prepare(**options)
        self.assertFalse(self.output.exists())


if __name__ == '__main__':
    unittest.main()
