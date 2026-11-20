# Final QIIME 2 2026.7 release gates

Current tested source: `b9c9b70c7da9c5c8f0b8d2cad189ce7e4edbc208`.
The [Build SIF workflow, run 44640635399](https://github.com/LinyangSun/Meta2Data/actions/runs/44640635399)
uses that exact revision for both jobs. Both the fresh native gate and the
final SIF gate passed. The image was published after the packaged tests
completed successfully.

## Fresh native installation

**Passed.** The native job installed the runtime in a new environment and ran
the complete compatibility suite. Its validation started at
`2026-11-19T23:50:58.790570+00:00` and finished at
`2026-11-20T00:18:32.260252+00:00`.

| Check | Observed result |
|---|---|
| Runtime | QIIME 2 2026.7.0 and scikit-learn 1.7.1 |
| Installed Conda packages | All 485 package URLs match the repository lock; no extra or missing packages |
| Unit tests | 317 passed, zero failures, errors or skips |
| PIP | 14 dataset/method cases passed: 6 DADA2 and 8 VSEARCH |
| TAXA | Five cases passed: three classifiers with DADA2/no-tree, plus GG2 with VSEARCH/SEPP and VSEARCH/de novo |
| Execution steps | All 40 recorded subprocess steps returned zero |

The PIP cases cover explicit and detected paired/single Illumina input, Ion
Torrent and PacBio for both methods, plus 454 and ONT for VSEARCH. Each case
retained five features, two samples and nonzero sample counts. These are small
generated fixtures processed by the real tools. Fixture sizes and validation
thread counts are not production analysis limits.

The taxonomy artifacts use the registered model UUIDs and the expected
`gg2Taxonomy.qza`, `silvaTaxonomy.qza` and `silva144Taxonomy.qza` names. The
SILVA 144 results contain four, five or seven ranks; the validator checked
the `d__`, `k__`, `p__`, `c__`, `o__`, `f__`, `g__` prefix order while allowing
assignments to stop at a higher rank. Table, representative-sequence and
taxonomy IDs, tree tips, and the new feature-table summary outputs also
passed the release checks.

Evidence:

- [Current native validation report](release-native-rebuild-validation.json)
- [Unedited current native unit-test output](release-native-rebuild-unit-tests.txt)
- [Actual installed Conda package URLs](release-native-environment-linux-64.lock)
- [Complete current native logs artifact](https://github.com/LinyangSun/Meta2Data/actions/runs/44640635399/artifacts/14587893546)

The installed lock is byte-for-byte identical to the retained lock above;
no duplicate copy was added. The downloaded native artifact ZIP has SHA-256
`fa18c2101fc6c859cc2006a4c19ea6cb5976a3234047205d7417d2ed72c3acb2`.
The complete Actions artifact has a 14-day retention period; the report,
unit-test output and installed lock above are retained in this repository.

### Earlier native pass

[Run 44633113721](https://github.com/LinyangSun/Meta2Data/actions/runs/44633113721)
passed the native gate on source
`3915439a5322a4e197e6373374776d9803d4c038`, from
`2026-11-19T22:58:29.540912+00:00` to
`2026-11-19T23:35:27.852543+00:00`. Its original evidence remains available:

- [Earlier native validation report](release-native-validation.json)
- [Unedited earlier native unit-test output](release-native-unit-tests.txt)
- [Earlier complete native logs artifact](https://github.com/LinyangSun/Meta2Data/actions/runs/44633113721/artifacts/14584783589)

The earlier artifact ZIP SHA-256 was
`30f2936cad625a2630de3f7a7215919157b34e08a256bbcabc780a01f8b9d150`, matching
GitHub's uploaded artifact digest. Its Actions artifact also has 14-day retention.

## Final SIF

**Passed and published.** The corrected fresh build, its `%test` environment
check, and the tests inside the sealed image all passed in
[job 135984600907](https://github.com/LinyangSun/Meta2Data/actions/runs/44640635399/job/135984600907).
The packaged validation started at `2026-11-20T00:27:42.306234+00:00` and
finished at `2026-11-20T00:53:29.611783+00:00`. All 30 recorded subprocess
steps returned zero.

| Packaged workflow | Observed result |
|---|---|
| PE DADA2 and PE VSEARCH | Both passed; three features and two samples each |
| Three classifiers with DADA2/no-tree | All passed with the registered classifier UUIDs and expected taxonomy filenames |
| GG2 with VSEARCH/SEPP | Passed; filtered table, sequences and tree tips agree |
| GG2 with VSEARCH/de novo | Passed; oriented table, sequences and rooted tree tips agree |
| SILVA 144 rank handling | Four- and seven-rank results passed the `d__`, `k__`, `p__`, `c__`, `o__`, `f__`, `g__` prefix-order checks |
| Summary outputs | Existing visualizations and the required feature/sample frequency artifacts were produced |

This is the final image's smoke suite: two PIP and five TAXA cases. The full
317-unit-test and 14-PIP suite was run in the native job, not repeated inside
the SIF. The sealed-image command uses `apptainer exec --cleanenv`; the
validator asserts that the entrypoint resolves to
`/opt/meta2data-current/bin/Meta2Data` and that the analysis tools come from
`/opt/conda/envs/meta2data/bin`. Only the working directory is bound from the
host for fixtures, results and logs.

Evidence:

- [Unedited final SIF validation report](release-sif-validation.json)
- [Complete SIF logs artifact](https://github.com/LinyangSun/Meta2Data/actions/runs/44640635399/artifacts/14590945171)

The downloaded artifact ZIP has SHA-256
`f63091def05d86a10e11bfcaf306aa44c52fd795ec8a2cf8cdfe11bd612e5273`, matching
GitHub's artifact digest. The complete Actions artifact has 14-day retention;
the report is retained here, and the complete logs were also downloaded to
`/workspaces/meta2data-validation/release-gates-44640635399/sif` in the
existing Codespace. Both native-run artifacts were likewise saved there in
their respective `release-gates-<run-id>/native` directories.

### First build failure and correction

The first build stopped during `%test`: `qiime info` attempted to create
`/opt/conda/envs/meta2data/etc/qiime2_config.toml` after the image filesystem
had become read-only. Environment installation and plugin-cache generation
had completed. The [unedited failure excerpt](release-sif-initial-failure.txt)
records the source revision, workflow and build job.

The failure was specific to image initialization. The corrected image passed
both `%test` and the final packaged workflow suite before publication.

The correction runs `qiime info > /dev/null` in `%post`, while the environment
is writable. The public command writes and reads the complete upstream
default configuration. Its nonzero exit still stops the build. No runtime
write permission, overlay, analysis parameter or validation gate is changed.
The behavior is defined by the pinned upstream
[configuration loader](https://github.com/rachis-org/rachis/blob/2026.7.0/src/rachis/sdk/parallel_config.py)
and [info command](https://github.com/rachis-org/rachis-cli/blob/2026.7.0/q2cli/builtin/info.py).

That default configuration contains Parsl worker counts calculated on the
build host. Meta2Data does not enable QIIME's `--parallel` or
`--parallel-config` execution mode. It passes the user thread budget directly
to the plugin parameters for DADA2, classification, SEPP and de novo trees;
those normal calls are not capped by the Parsl values. Direct, explicitly
parallel QIIME calls outside Meta2Data would use those defaults unless given
a separate parallel configuration. See the upstream
[CLI execution branches](https://github.com/rachis-org/rachis-cli/blob/2026.7.0/q2cli/commands.py).

## Published SIF v1.1

The release is publicly available as
`oras://ghcr.io/linyangsun/meta2data:v1.1`.

```bash
apptainer pull Meta2Data-v1.1.sif oras://ghcr.io/linyangsun/meta2data:v1.1
```

| Identity | Value |
|---|---|
| SIF source revision | `b9c9b70c7da9c5c8f0b8d2cad189ce7e4edbc208` |
| SIF bytes | 2,147,151,872 |
| SIF SHA-256 | `84f0093b37cb976cab8974195a20149d291994b677c8d37da16431c86d62dce1` |
| OCI manifest SHA-256 | `e07f288ce27d52fb2aa39271856cd8b43f3df0905039d3aa6aa868c65b6b976b` |

The successful build first published `qiime2026.7-silva144-v1`. The requested
short release name was then published by
[run 44701077801](https://github.com/LinyangSun/Meta2Data/actions/runs/44701077801)
at `2026-11-20T08:15:29Z`. Both tags refer to the exact same manifest and SIF
bytes; the image was not rebuilt or modified for the name change. The job
verified the original build's successful conclusion and source revision,
the fixed manifest and SIF hashes, and byte-identical anonymous retrieval of
the new tag. It refused to replace any different existing target manifest.

Evidence is retained in the [publication receipt](release-sif-v1.1-publication.json)
and [original OCI manifest](release-sif-v1.1-manifest.json). The complete
[publication artifact](https://github.com/LinyangSun/Meta2Data/actions/runs/44701077801/artifacts/14611929183)
has ZIP SHA-256
`a30fd77fe0adf4523f1eb94a56019c52f3a9edc063c8dd3b09e7cfcf67c1d37e`.
This verifies the public manifest without downloading the 2.15 GB SIF again;
the SIF file itself was hashed and exercised in the completed build job.

The earlier tag-publication attempt,
[run 44700765831](https://github.com/LinyangSun/Meta2Data/actions/runs/44700765831),
stopped on an HTTP 403 while reading the source run through the unauthenticated
GitHub API, before writing the target tag. The successful attempt used the
job's authenticated, read-only Actions access. The publication-only workflow
is preserved on `release/sif-v1.1`; it is not part of the upgrade's production
branch or the normal fresh-build workflow.

## Final source documentation

Changes after the tested source are documentation, retained validation
evidence and one lock-file status comment. The production entrypoints,
analysis scripts, installer, container definition, workflows, tests,
`env.yml` and release constraints are unchanged. Every non-comment line in
the lock remains identical. Its raw file SHA changed only because the
candidate-status comment now records completed validation; historical reports
retain the SHA of the exact file they tested.

## Scope

These gates establish software compatibility and output handling. They do
not evaluate database accuracy or claim an exhaustive combination of every
classifier, platform and tree mode. Full-reference training counts and class
counts are recorded separately in the [classifier build audit](classifier-build-audit.md).
