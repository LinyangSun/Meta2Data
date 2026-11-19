# Final QIIME 2 2026.7 release gates

Tested source: `3915439a5322a4e197e6373374776d9803d4c038`.
The [Build SIF workflow, run 44633113721](https://github.com/LinyangSun/Meta2Data/actions/runs/44633113721)
checks that exact revision for both jobs.

## Fresh native installation

**Passed.** The native job installed the runtime in a new environment and ran
the complete compatibility suite. Its validation started at
`2026-11-19T22:58:29.540912+00:00` and finished at
`2026-11-19T23:35:27.852543+00:00`.

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

- [Native validation report](release-native-validation.json)
- [Unedited native unit-test output](release-native-unit-tests.txt)
- [Actual installed Conda package URLs](release-native-environment-linux-64.lock)
- [Complete native logs artifact](https://github.com/LinyangSun/Meta2Data/actions/runs/44633113721/artifacts/14584783589)

The downloaded native artifact ZIP has SHA-256
`30f2936cad625a2630de3f7a7215919157b34e08a256bbcabc780a01f8b9d150`, matching
GitHub's uploaded artifact digest. The complete Actions artifact has a 14-day
retention period; the report, unit-test output and installed lock above are
retained in this repository.

## Final SIF

The first build stopped during `%test`: `qiime info` attempted to create
`/opt/conda/envs/meta2data/etc/qiime2_config.toml` after the image filesystem
had become read-only. Environment installation and plugin-cache generation
had completed. The [unedited failure excerpt](release-sif-initial-failure.txt)
records the source revision, workflow and build job.

The failure is specific to image initialization. The native results above
remain valid. The corrected image must pass `%test` and the final packaged
workflow suite before publication or merging.

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

## Scope

These gates establish software compatibility and output handling. They do
not evaluate database accuracy or claim an exhaustive combination of every
classifier, platform and tree mode. Full-reference training counts and class
counts are recorded separately in the [classifier build audit](classifier-build-audit.md).
