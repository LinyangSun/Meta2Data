"""Offline regressions for current-run metadata outputs and archive platforms."""
import contextlib
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
DEPENDENCIES = ('pandas', 'numpy', 'Bio', 'requests')
MISSING = [name for name in DEPENDENCIES if importlib.util.find_spec(name) is None]
if not MISSING:
    import pandas as pd
    import metadata_downloader as metadata
    import py_16s


@unittest.skipIf(MISSING, 'MetaDL tests require: ' + ', '.join(MISSING))
class MetadataResultsTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.stack = contextlib.ExitStack()
        self.addCleanup(self.stack.close)
        self.stack.enter_context(contextlib.redirect_stdout(io.StringIO()))

    def test_keyword_mode_reads_only_generated_accessions_and_applies_api_key(self):
        searched = self.root / 'searched_keywords'
        searched.mkdir()
        (searched / 'bioproject_ids.txt').write_text('PRJNA1\nPRJNA2\n')
        (searched / 'search_summary.txt').write_text('Total queries: 1\nTime: today\n')
        instance = mock.Mock()
        instance.search_and_download_batch.return_value = {'dirs': {'results': searched}}
        instance.combine_batch_results.return_value = pd.DataFrame({'accession': ['PRJNA1', 'PRJNA2']})
        received = []

        def run(input_path, *args):
            received.extend(metadata.read_input_ids(input_path))
            return {'input_ids': received}

        def create_downloader():
            self.assertEqual(metadata.Entrez.api_key, 'test-key')
            return instance

        argv = ['metadata_downloader.py', '--keywords', '--field', '16S',
                '--organism', 'Apis', '-o', str(self.root), '-k', 'test-key']
        with mock.patch.object(sys, 'argv', argv), \
                mock.patch.object(metadata.Entrez, 'api_key', None, create=True), \
                mock.patch.object(metadata, 'BioProjectDownloader', side_effect=create_downloader), \
                mock.patch.object(metadata, 'run_unified_pipeline', side_effect=run), \
                self.assertRaises(SystemExit) as exit_code:
            metadata.main()
        self.assertEqual(exit_code.exception.code, 0)
        self.assertEqual(received, ['PRJNA1', 'PRJNA2'])

    def test_input_directory_keeps_all_id_files_and_deduplicates(self):
        (self.root / 'a.txt').write_text('# comment\nPRJNA1\nPRJNA2\n')
        (self.root / 'b.txt').write_text('PRJNA2\nSAMN3\n')
        self.assertEqual(metadata.read_input_ids(self.root), ['PRJNA1', 'PRJNA2', 'SAMN3'])

    def test_keyword_result_without_optional_annotations_is_usable(self):
        xml = self.root / 'result.xml'
        xml.write_text('<Package><Project><ProjectID><ArchiveID accession="PRJNA1" />'
                       '</ProjectID></Project></Package>')
        result = {'queries': [{'files': {'ncbi': xml}}], 'dirs': {'results': self.root}}
        combined = metadata.BioProjectDownloader().combine_batch_results(result)
        self.assertEqual(combined['accession'].tolist(), ['PRJNA1'])
        self.assertIn('description', combined.columns)
        self.assertEqual((self.root / 'bioproject_ids.txt').read_text(), 'PRJNA1\n')

    def test_empty_keyword_results_replace_previous_search_outputs(self):
        (self.root / 'combined_results.csv').write_text('accession\nPRJNA999\n')
        (self.root / 'bioproject_ids.txt').write_text('PRJNA999\n')
        result = {'queries': [], 'dirs': {'results': self.root}}
        combined = metadata.BioProjectDownloader().combine_batch_results(result)
        self.assertTrue(combined.empty)
        self.assertTrue(pd.read_csv(self.root / 'combined_results.csv').empty)
        self.assertEqual((self.root / 'bioproject_ids.txt').read_text(), '')

    def test_pipeline_merges_only_requested_checkpoint_not_old_project_or_group(self):
        inputs = self.root / 'ids'
        inputs.mkdir()
        (inputs / 'projects.txt').write_text('PRJNA1\n')
        output = self.root / 'metadata'
        tmp = output / 'tmp'
        state = metadata.StateManager(tmp)
        current = pd.DataFrame({'Run': ['SRR1'], 'Bioproject': ['PRJNA1']})
        old = pd.DataFrame({'Run': ['SRR999'], 'Bioproject': ['PRJNA999']})
        current_path = tmp / 'PRJNA1.processed.csv'
        current.to_csv(current_path, index=False)
        state.mark_complete('PRJNA1', current_path, 1, 'NCBI')
        for name in ('PRJNA999', 'BIOSAMPLE_INPUT', 'SRA_INPUT'):
            old.to_csv(tmp / (name + '.processed.csv'), index=False)
        with mock.patch.object(metadata, 'fetch_bioproject_description', return_value=''), \
                mock.patch.object(metadata, 'build_bioproject_absdesc'), \
                mock.patch.object(metadata, 'download_ncbi_metadata',
                                  side_effect=AssertionError('checkpoint should be reused')):
            result = metadata.run_unified_pipeline(inputs, output, max_workers=1)
        self.assertEqual(result['final_df']['Run'].tolist(), ['SRR1'])
        self.assertEqual(pd.read_csv(output / 'all_metadata_merged.csv')['Run'].tolist(), ['SRR1'])
        self.assertTrue((tmp / 'PRJNA999.processed.csv').exists())

    def test_pipeline_summary_counts_biosample_ids_not_grouped_result_files(self):
        for include_project in (False, True):
            with self.subTest(include_project=include_project):
                case = self.root / ('mixed' if include_project else 'biosamples')
                inputs = case / 'ids'
                inputs.mkdir(parents=True)
                ids = ['SAMN1', 'SAMN2'] + (['PRJNA3'] if include_project else [])
                (inputs / 'input.txt').write_text('\n'.join(ids) + '\n')
                output = case / 'metadata'
                tmp = output / 'tmp'
                state = metadata.StateManager(tmp)
                # Unrelated historical states must not enter this run's totals.
                state.mark_status('SAMN999', metadata.STATUS_HAS_DATA)
                state.mark_status('SAMN998', metadata.STATUS_NO_DATA)
                if include_project:
                    project = pd.DataFrame({'Run': ['SRR3'], 'Bioproject': ['PRJNA3']})
                    checkpoint = tmp / 'PRJNA3.processed.csv'
                    project.to_csv(checkpoint, index=False)
                    state.mark_complete('PRJNA3', checkpoint, 1, 'NCBI')
                biosamples = pd.DataFrame({'Run': ['SRR1', 'SRR2'],
                                          'Biosample': ['SAMN1', 'SAMN2'],
                                          'Bioproject': ['PRJNA1', 'PRJNA2']})
                captured = io.StringIO()
                with contextlib.redirect_stdout(captured), \
                        mock.patch.object(metadata, 'download_ncbi_metadata_from_biosamples',
                                          return_value=biosamples), \
                        mock.patch.object(metadata, '_assign_bioproject_descriptions'), \
                        mock.patch.object(metadata, 'fetch_bioproject_description', return_value=''), \
                        mock.patch.object(metadata, 'build_bioproject_absdesc'), \
                        mock.patch.object(metadata, 'download_ncbi_metadata',
                                          side_effect=AssertionError('checkpoint should be reused')):
                    result = metadata.run_unified_pipeline(inputs, output, max_workers=1)
                total = len(ids)
                self.assertEqual(len(result['results']), 1 + int(include_project))
                self.assertEqual(len(pd.read_csv(tmp / 'BIOSAMPLE_INPUT.processed.csv')), 2)
                self.assertEqual(dict(zip(result['status_df']['ID'], result['status_df']['Status'])),
                                 {accession: metadata.STATUS_HAS_DATA for accession in ids})
                self.assertEqual(set(result['final_df']['Run']),
                                 {'SRR1', 'SRR2', 'SRR3'} if include_project else {'SRR1', 'SRR2'})
                summary = captured.getvalue().split('PIPELINE COMPLETE', 1)[1]
                self.assertIn(f'Total input IDs: {total}\n', summary)
                self.assertIn(f'  With valid Run data: {total}\n', summary)
                self.assertIn('  No data/No Run info: 0\n', summary)

    def test_pipeline_summary_separates_missing_data_from_invalid_ids(self):
        inputs = self.root / 'ids'
        inputs.mkdir()
        (inputs / 'input.txt').write_text('SAMN1\nSAMN2\nnot-an-accession\n')
        output = self.root / 'metadata'
        biosamples = pd.DataFrame({'Run': ['SRR1'], 'Biosample': ['SAMN1'],
                                  'Bioproject': ['PRJNA1']})
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured), \
                mock.patch.object(metadata, 'download_ncbi_metadata_from_biosamples',
                                  return_value=biosamples), \
                mock.patch.object(metadata, '_assign_bioproject_descriptions'), \
                mock.patch.object(metadata, 'build_bioproject_absdesc'):
            result = metadata.run_unified_pipeline(inputs, output, max_workers=1)
        self.assertEqual(dict(zip(result['status_df']['ID'], result['status_df']['Status'])),
                         {'SAMN1': metadata.STATUS_HAS_DATA, 'SAMN2': metadata.STATUS_NO_DATA,
                          'not-an-accession': metadata.STATUS_INVALID_FORMAT})
        summary = captured.getvalue().split('PIPELINE COMPLETE', 1)[1]
        self.assertIn('Total input IDs: 3\n', summary)
        self.assertIn('  With valid Run data: 1\n', summary)
        self.assertIn('  No data/No Run info: 1\n', summary)
        self.assertIn('  invalid_format: 1\n', summary)

    def test_empty_current_results_replace_previous_final_reports(self):
        (self.root / 'all_metadata_merged.csv').write_text('Run,Bioproject\nSRR999,PRJNA999\n')
        (self.root / 'column_description.tsv').write_text('stale\n')
        (self.root / 'RecordWithoutRUNinfo.csv').write_text('stale\n')
        result = metadata.merge_all_results([], self.root)
        self.assertTrue(result.empty)
        self.assertTrue(pd.read_csv(self.root / 'all_metadata_merged.csv').empty)
        columns = pd.read_csv(self.root / 'column_description.tsv', sep='\t')
        self.assertIn('Run', columns['ColumnName'].tolist())
        self.assertFalse((self.root / 'RecordWithoutRUNinfo.csv').exists())

    def test_without_run_report_reflects_current_results_only(self):
        frame = pd.DataFrame({'Run': ['SRR1', None], 'Bioproject': ['PRJNA1', 'PRJNA1']})
        metadata.merge_all_results([{'df': frame}], self.root)
        self.assertEqual(len(pd.read_csv(self.root / 'RecordWithoutRUNinfo.csv')), 1)
        metadata.merge_all_results([{'df': frame.iloc[:1]}], self.root)
        self.assertFalse((self.root / 'RecordWithoutRUNinfo.csv').exists())


@unittest.skipIf(MISSING, 'Platform tests require: ' + ', '.join(MISSING))
class PlatformLookupTests(unittest.TestCase):
    XML = '''<EXPERIMENT_PACKAGE_SET>
      <EXPERIMENT_PACKAGE><EXPERIMENT><PLATFORM><ILLUMINA /></PLATFORM></EXPERIMENT>
        <RUN_SET><RUN accession="SRR1" /><RUN accession="SRR2" /></RUN_SET>
      </EXPERIMENT_PACKAGE>
      <EXPERIMENT_PACKAGE><EXPERIMENT><PLATFORM><PACBIO_SMRT /></PLATFORM></EXPERIMENT>
        <RUN_SET><RUN accession="SRR3" /></RUN_SET>
      </EXPERIMENT_PACKAGE>
    </EXPERIMENT_PACKAGE_SET>'''

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.stack = contextlib.ExitStack()
        self.addCleanup(self.stack.close)
        self.stack.enter_context(contextlib.redirect_stderr(io.StringIO()))
        self.stack.enter_context(mock.patch.object(py_16s, '_configure_entrez'))
        self.stack.enter_context(mock.patch.object(metadata.Entrez, 'esearch', side_effect=lambda **kw: io.StringIO('')))
        self.stack.enter_context(mock.patch.object(metadata.Entrez, 'read', return_value={'IdList': ['1', '2']}))
        self.fetch = self.stack.enter_context(mock.patch.object(metadata.Entrez, 'efetch',
                                                       side_effect=lambda **kw: io.StringIO(self.XML)))

    def test_single_lookup_selects_matching_run_not_first_experiment(self):
        self.assertEqual(py_16s._get_platform_from_ncbi('SRR3'), 'PACBIO_SMRT')
        self.assertEqual(self.fetch.call_args.kwargs['id'], '1,2')

    def test_missing_run_does_not_borrow_another_platform(self):
        self.assertIsNone(py_16s._get_platform_from_ncbi('SRR404'))
        self.assertIsNone(py_16s._parse_cncb_platform_response('Run,Platform\nCRR1,Illumina\n', 'CRR404'))
        self.assertIsNone(py_16s._parse_cncb_platform_response('Other,Platform\nCRR1,Illumina\n', 'CRR1'))

    def test_batch_lookup_includes_second_run_and_shared_query(self):
        pairs = self.root / 'pairs.tsv'
        pairs.write_text('project_A\tSRR2\nproject_B\tSRR3\nproject_C\tSRR2\n')
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            py_16s.batch_get_sequencing_platforms(pairs)
        self.assertEqual(set(out.getvalue().splitlines()),
                         {'project_A\tILLUMINA', 'project_B\tPACBIO_SMRT', 'project_C\tILLUMINA'})

    def test_cncb_matching_run_uses_its_platform(self):
        data = 'Run,Platform\nCRR1,Illumina\nCRR2,PacBio Sequel\n'
        self.assertEqual(py_16s._parse_cncb_platform_response(data, 'CRR2'), 'PACBIO_SMRT')


if __name__ == '__main__':
    unittest.main()
