"""PacBio fixtures derive their known insert from native primer binding sites."""
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('fixture_generator',
                                            ROOT / 'tests/integration/generate_fixture.py')
fixture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fixture)

# The insert has no primer sequence; expected products are defined independently
# of either cutadapt or the pipeline output.
FORWARD = 'AGAGTTTGATCCTGGCTCAG'
REVERSE_TARGET = 'AAGTCGTAACAAGGTAGCCGTA'
INSERT = 'GCTAACCGAT' * 142
DOWNSTREAM = 'CTGGAAGGTGCGGCTG'


def rc(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


class PacbioFixtureConstructionTests(unittest.TestCase):
    def test_complete_native_reverse_site_replaces_site_and_discards_downstream(self):
        reference = FORWARD + INSERT + REVERSE_TARGET + DOWNSTREAM
        result = fixture.build_pacbio_template(reference)
        self.assertEqual(result['pacbio'], FORWARD + INSERT + REVERSE_TARGET)
        self.assertEqual(result['pacbio_insert'], INSERT)
        binding = result['pacbio_binding']
        self.assertEqual(binding['reverse_start'], len(FORWARD) + len(INSERT))
        self.assertEqual(binding['reverse_covered_end'], len(reference) - len(DOWNSTREAM))
        self.assertEqual(binding['reverse_reference_bases'], 22)
        self.assertEqual(binding['downstream_reference_bases_removed'], len(DOWNSTREAM))
        self.assertEqual(reference[binding['insert_start']:binding['insert_end']], INSERT)

    def test_terminal_partial_reverse_site_is_completed_once(self):
        for covered in (19, 20, 21):
            with self.subTest(covered=covered):
                reference = FORWARD + INSERT + REVERSE_TARGET[:covered]
                result = fixture.build_pacbio_template(reference)
                self.assertEqual(result['pacbio'], FORWARD + INSERT + REVERSE_TARGET)
                self.assertEqual(result['pacbio_insert'], INSERT)
                self.assertEqual(result['pacbio_binding']['reverse_reference_bases'], covered)
                self.assertEqual(result['pacbio_binding']['downstream_reference_bases_removed'], 0)

    def test_reference_forward_variants_are_replaced_without_changing_insert(self):
        reference = 'AGAGTTCGATCCTGGCTCAG' + INSERT + REVERSE_TARGET
        result = fixture.build_pacbio_template(reference)
        self.assertEqual(result['pacbio'], FORWARD + INSERT + REVERSE_TARGET)
        self.assertEqual(result['pacbio_insert'], INSERT)

    def test_missing_reverse_site_is_rejected(self):
        with self.assertRaisesRegex(ValueError, '1492R.*unique'):
            fixture.build_pacbio_template(FORWARD + INSERT)

    def test_repeated_native_reverse_sites_are_rejected(self):
        with self.assertRaisesRegex(ValueError, '1492R.*unique'):
            fixture.build_pacbio_template(FORWARD + INSERT + REVERSE_TARGET * 2)

    def test_reverse_site_far_from_reference_end_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'terminal 100 bp'):
            fixture.build_pacbio_template(FORWARD + INSERT + REVERSE_TARGET + 'G' * 101)

    def test_select_templates_reports_skipped_missing_reverse_site(self):
        # This record passes the short-read binding and reference-length gates.
        forward = 'CCTACGGGAGGCAGCAG'
        reverse = rc('GACTACAGGGGTATCTAATCC')
        prefix = FORWARD + 'A' * (300 - len(FORWARD))
        reference = prefix + forward + 'C' * 400 + reverse
        reference += 'G' * (1470 - len(reference))
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / 'reference.fasta'
            fasta.write_text('>no_1492R\n' + reference + '\n')
            with self.assertRaisesRegex(ValueError, 'PacBio binding-site rejections:.*1492R'):
                fixture.select_templates(fasta, 1)

    @unittest.skipUnless(shutil.which('cutadapt'), 'Real cutadapt is required')
    def test_complete_and_partial_fixtures_roundtrip_both_orientations_exactly(self):
        references = [FORWARD + INSERT + REVERSE_TARGET + DOWNSTREAM,
                      FORWARD + INSERT + REVERSE_TARGET[:19]]
        for reference in references:
            with self.subTest(reference_length=len(reference)), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source = root / 'input'
                target = root / 'output'
                source.mkdir()
                amplicon = fixture.build_pacbio_template(reference)['pacbio']
                forward_quality = ''.join(chr(40 + index % 35) for index in range(len(amplicon)))
                reverse_quality = ''.join(chr(40 + (index + 9) % 35) for index in range(len(amplicon)))
                contents = (f'@forward unchanged description\n{amplicon}\n+\n{forward_quality}\n'
                            f'@reverse unchanged description\n{rc(amplicon)}\n+\n{reverse_quality}\n')
                path = source / 'sample.fastq'
                path.write_text(contents)
                command = [sys.executable, str(ROOT / 'scripts/explicit_primers.py'),
                           '--input', str(source), '--output', str(target),
                           '--forward', FORWARD, '--reverse', fixture.PACBIO_REVERSE,
                           '--mixed-orientation']
                result = subprocess.run(command, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual((target / 'sample.fastq').read_text(),
                    f'@forward unchanged description\n{INSERT}\n+\n'
                    f'{forward_quality[len(FORWARD):-len(REVERSE_TARGET)]}\n'
                    f'@reverse unchanged description\n{INSERT}\n+\n'
                    f'{reverse_quality[::-1][len(FORWARD):-len(REVERSE_TARGET)]}\n')
                self.assertEqual(path.read_text(), contents)


if __name__ == '__main__':
    unittest.main()
