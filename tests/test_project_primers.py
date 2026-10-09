import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from project_primers import prepare, read_mapping, validate_primers


class ProjectPrimersTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.metadata = self.root / 'input.csv'
        self.metadata.write_text('Project,Accession,Other\nPRJNA1,SRR1,a\nPRJNA2,ERR2,b\nPRJNA1,SRR3,c\nPRJNA3,,d\n')
        self.original = self.metadata.read_bytes()

    def prepare(self, projects=None, forward=None, reverse=None):
        return prepare(self.metadata, self.root / 'output', 'Project', 'Accession',
                       projects or ['PRJNA2', 'PRJNA1'], forward or ['acgt', 'NNRY'], reverse)

    def test_subset_preserves_all_selected_runs_and_original(self):
        subset, mapping = self.prepare(reverse=['TTAA', 'CCGG'])
        with subset.open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([row['Accession'] for row in rows], ['SRR1', 'ERR2', 'SRR3'])
        self.assertEqual(self.metadata.read_bytes(), self.original)
        self.assertEqual(read_mapping(mapping), {'PRJNA2': {'forward': 'ACGT', 'reverse': 'TTAA'},
                                                'PRJNA1': {'forward': 'NNRY', 'reverse': 'CCGG'}})
        self.assertEqual(json.loads(mapping.read_text())['project_order'], ['PRJNA2', 'PRJNA1'])

    def test_forward_only_never_fills_reverse(self):
        _, mapping = self.prepare()
        self.assertTrue(all(item['reverse'] == '' for item in read_mapping(mapping).values()))

    def test_validation_rejects_ambiguous_or_invalid_lists(self):
        cases = [(['P1', 'P1'], ['AC', 'TG'], None),
                 (['P1', 'P2'], ['AC'], None),
                 (['P1', 'P2'], ['AC', 'TG'], ['AC']),
                 (['P1'], ['primer-file.fa'], None),
                 (['../P1'], ['AC'], None),
                 (['P1'], [''], None)]
        for args in cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                validate_primers(*args)

    def test_missing_or_empty_project_fails_without_writing(self):
        for project, message in [('PRJNA404', 'not present'), ('PRJNA3', 'no usable Run')]:
            with self.subTest(project=project), self.assertRaisesRegex(ValueError, message):
                self.prepare([project], ['ACGT'])
        self.assertFalse((self.root / 'output').exists())

    def test_invalid_run_fails_before_output(self):
        self.metadata.write_text('Project,Accession\nPRJNA1,not-a-run\n')
        with self.assertRaisesRegex(ValueError, 'Invalid Run'):
            self.prepare(['PRJNA1'], ['ACGT'])
        self.assertFalse((self.root / 'output').exists())

    def test_bad_columns_and_overwriting_input_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'missing column'):
            prepare(self.metadata, self.root / 'output', 'wrong', 'Accession', ['PRJNA1'], ['AC'])
        source = self.root / 'selected_metadata.csv'
        source.write_bytes(self.original)
        with self.assertRaisesRegex(ValueError, 'Input metadata must differ'):
            prepare(source, self.root, 'Project', 'Accession', ['PRJNA1'], ['AC'])
        self.assertEqual(source.read_bytes(), self.original)

    def test_corrupt_mapping_cannot_trigger_automatic_detection(self):
        _, mapping = self.prepare()
        data = json.loads(mapping.read_text())
        data['projects']['PRJNA1']['forward'] = ''
        mapping.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'Invalid primer'):
            read_mapping(mapping)


if __name__ == '__main__':
    unittest.main()
