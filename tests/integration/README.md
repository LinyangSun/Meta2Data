# Synthetic archive boundary

On a SLURM cluster, place fixture generation, image builds and test commands in
a batch script and submit it with `sbatch`; do not run them on the login node.

`archive_shim.py` is an integration-test helper. Its explicit test PATH contains
only `wget`, `python` and `python3` wrappers:

- `wget` serves declared ENA probes, four-column filereports and FASTQ files.
  File bytes must match the manifest's MD5. Unlisted URLs and projects fail;
  the wrapper never calls the network or a real `wget`.
- Python intercepts only `py_16s.py batch_get_sequencing_platforms` and
  `py_16s.py get_sequencing_platform`. Unknown Runs fail. Every other Python
  invocation executes the real interpreter and preserves its exit status.
- `curl`, FASTQ validation, primer handling, adapter protection and all biological
  tools remain real. Reference downloads can therefore still use the network.

The fixture generator supplies `public_fixture.json` with
`schema_version: 1`, `synthetic: true`, `projects[project].runs` and
`runs[run].{platform,layout,files}`. Each file has `path`, `md5` and an optional
`url`; relative paths resolve against the manifest's directory. The default URL
is `ftp.sra.ebi.ac.uk/fixtures/<filename>`. Extra generation metadata is retained
but ignored by the shim.

`M2D_TEST_ARCHIVE_MANIFEST` selects the fixture. `M2D_TEST_ARCHIVE_TRACE` receives
append-only JSONL records for accepted requests, rejected requests and Python
passthroughs. `M2D_TEST_REAL_PYTHON` must be an absolute real interpreter path,
not one of the wrappers. The wrappers themselves use `/usr/bin/python3` so
adding their directory to PATH cannot recurse.

PacBio fixtures are PCR amplicons bounded by known primer sites. The generator
locates the native reverse-primer binding region near the reference end and
replaces that region with the experimental primer; it must not append another
primer after the complete reference. The original binding coordinates and the
expected primer-free inserts are recorded before reads are simulated. Both
orientations are generated, and validation compares representative sequences
against the expected inserts exactly.

Example with a complete SIF and fixtures already generated:

```bash
M2D_SOURCE=/absolute/path/to/Meta2Data-main
M2D_FIXTURE=/absolute/path/to/fixtures
M2D_RUN=/absolute/path/to/online-test
M2D_SIF=/absolute/path/to/Meta2Data.sif
mkdir -p "$M2D_RUN"
/usr/bin/python3 "$M2D_SOURCE/tests/integration/archive_shim.py" install \
  --directory "$M2D_RUN/fixture-bin"
M2D_CONTAINER_PATH="$(apptainer exec --cleanenv "$M2D_SIF" /usr/bin/printenv PATH)"

apptainer exec --cleanenv \
  --bind "$M2D_SOURCE:$M2D_SOURCE:ro" \
  --bind "$M2D_FIXTURE:$M2D_FIXTURE:ro" --bind "$M2D_RUN:$M2D_RUN" \
  --pwd "$M2D_RUN" \
  --env "PATH=$M2D_RUN/fixture-bin:$M2D_CONTAINER_PATH" \
  --env "M2D_TEST_ARCHIVE_MANIFEST=$M2D_FIXTURE/public_fixture.json" \
  --env "M2D_TEST_ARCHIVE_TRACE=$M2D_RUN/archive-trace.jsonl" \
  --env M2D_TEST_REAL_PYTHON=/opt/conda/envs/qiime2-amplicon-2024.10/bin/python3 \
  "$M2D_SIF" Meta2Data AmpliconPIP \
  --public-m "$M2D_FIXTURE/online_metadata.csv" \
  --public-bioproject-colNAME Bioproject --public-sra-colNAME Run \
  --vsearch -t 4
```

The source-directory bind above exposes this **test helper**, while the product
command executes from the SIF. Do not add `fixture-bin` to a normal analysis
PATH. This simulated archive run exercises the real downstream pipeline but
does not validate live ENA/NCBI connectivity or SRA Toolkit fallback.

Run the portable shim checks without a SIF or network:

```bash
python3 -m unittest discover -s tests -p 'test_archive_shim.py' -q
```
