# Parameter, output-name, and workflow audit

Date: 19 November 2026. Environment: QIIME 2 2026.7.0, scikit-learn 1.7.1,
Python 3.12.13, and standalone VSEARCH 2.30.4, using the pinned 485-package
minimal runtime in the 4-core / 16 GB Codespace.

The parameter/file audit, 14 PIP cases, SILVA 144 no-tree and both tree modes,
three real taxonomy-cache scenarios, and the isolated 317-test suite passed.
Compatible replacement GG2 and SILVA 138.1 models and final SIF validation are
still pending. This is not the complete all-classifier release gate.

## Revisions actually exercised

| Revision | Evidence |
|---|---|
| `82d41e0` | All 14 PIP cases and the first complete SILVA 144 no-tree run |
| `86dea0b` | Both SILVA 144 tree workflows and initial artifact-identity validation |
| `550d440` | Production cache-provenance fix, three real resume cases, tree artifact validation, current help, and all 317 unit tests |

The later changes do not alter PIP processing. The production change at
550d440 affects taxonomy reuse; the real resume checks below exercise it.
Driver snapshots in `audit-drivers/` record the paths and commands used in
this Codespace. They are audit evidence, not public product entrypoints.

## Public parameters and filenames

The selector is still `--classifier NAME`. New canonical values identify
releases explicitly; old aliases retain their database identities.

| Selection | Retained alias | Input model in results/db | Final taxonomy output |
|---|---|---|---|
| `gg2-2024.09`, default | `greengenes` | `gg2-2024.09-full-length-sklearn-1.7.1.qza` | `gg2Taxonomy.qza`, unchanged |
| `silva-138.1` | `silva` | `silva-138.1-full-length-sklearn-1.7.1.qza` | `silvaTaxonomy.qza`, unchanged |
| `silva-144` | None | `silva-144-full-length-sklearn-1.7.1.qza` | `silva144Taxonomy.qza`, new |

There is no `--silva144`/`--silva144.x` flag or user-supplied classifier-path
option. The `silva` alias does not select the newest release. The default still
means GG2 2024.09. Confidence remains 0.7, threads remain 4, and the default tree
mode remains SEPP. Existing PIP table/representative-sequence and TAXA table/tree
basenames remain unchanged. Classifier selection remains a CLI option.

Current availability is stated in README and help: SILVA 144 downloads from the
pinned official source; the compatible GG2 and SILVA 138.1 models are not yet
published. Their valid parameters, including the default, stop with NOT READY
unless a compatible maintainer-built resource is installed. Old 1.4.2 models
are not silently substituted.

QIIME's new DADA2 base-transition artifact is an intermediate under
`tmp/step_05_denoise/`; normal successful cleanup removes it unless
`KEEP_INTERMEDIATE=1`. Each TAXA summary retains its old .qzv basename and adds
feature-frequencies and sample-frequencies .qza files alongside it under
`tmp/`. See the exact table in [the upgrade document](../../qiime2-2026.7-upgrade.md#public-parameters-and-output-files).

## Concrete issues corrected

1. README/help documentation now consistently identifies the canonical default,
   exact model-cache filenames, final result filenames, and pending legacy-model
   availability. The parameter reference and pipeline guide match these choices.
2. A valid ZIP placed under an expected taxonomy filename could previously be
   reused if current state was retained, even when the artifact came from another
   model, reads artifact, or confidence. Reuse now checks the direct
   feature-classifier.classify_sklearn action and its classifier UUID, reads UUID,
   and numeric confidence. Missing or mismatched provenance triggers recomputation.
3. The environment installer now prints the exact component-only retry command
   if Conda succeeds but standalone SRA setup fails. Existing-file/provenance
   conflicts still require inspection; overwrite protections remain in place.
4. Container instructions now describe running the candidate branch workflow
   before merging, using the existing workflow_dispatch entrypoint.

## Real PIP workflows

Fixture generation used five real reference templates with --reads 8000 and
--pacbio-reads 2000. Extra platform cases reuse those generated files; per-case
sample counts are recorded in the JSON evidence. The checks require valid artifact types,
matching feature IDs, expected sample identities, positive bounded read counts,
and the exact expected final basenames. The PacBio DADA2 representatives also
match the known primer-free inserts.

| Method/case | Features | Result |
|---|---:|---|
| dada2/ion | 5 | Passed |
| dada2/pacbio | 5 | Passed |
| dada2/pe_auto | 5 | Passed |
| dada2/pe_explicit | 5 | Passed |
| dada2/se_auto | 5 | Passed |
| dada2/se_explicit | 5 | Passed |
| vsearch/ion | 5 | Passed |
| vsearch/ls454 | 5 | Passed |
| vsearch/ont | 5 | Passed |
| vsearch/pacbio | 5 | Passed |
| vsearch/pe_auto | 5 | Passed |
| vsearch/pe_explicit | 5 | Passed |
| vsearch/se_auto | 5 | Passed |
| vsearch/se_explicit | 5 | Passed |

Detailed products, sample counts, commands, and times: [pip-workflow-audit.json](pip-workflow-audit.json).

## Real SILVA 144 outputs

The no-tree run produced 5 features and the expected
`silva144Taxonomy.qza`, with no incorrectly named GG2 or old-SILVA taxonomy.
Its classifier UUID is `b7986601-1851-4fa2-ace5-98bacf660fa3`; default confidence is 0.7.
Observed lineage lengths were [4, 7], with ordered prefixes
`d__/k__/p__/c__/o__/f__/g__`. Higher-rank truncated calls are retained as produced.
All 12 QIIME output archives were validated, including the two
three-file summary bundles. See [silva144-workflow.json](silva144-workflow.json)
and the [exported taxonomy](silva144-workflow-taxonomy.tsv).

| Tree mode | Retained features | Taxonomy features | Tree tips | Result |
|---|---:|---:|---:|---|
| sepp | 5 | 5 | 203457 | Passed |
| singleV | 5 | 5 | 5 | Passed |

Tree checks verify rooted-tree type, retained feature IDs in the tree, table and
sequence consistency, selected classifier/reads/confidence provenance, and actual
summary/frequency artifact types. Exact filenames and UUIDs are in
[tree-workflow-audit.json](tree-workflow-audit.json).

## Real cache regression

All three cases used the actual official SILVA 144 model and the same complete
PIP output; no substitute classifier was trained.

| Case | Observed behavior |
|---|---|
| Repeat confidence 0.7 | Every artifact UUID and SHA256 remained unchanged; taxonomy was reused |
| Change confidence to 0.8 | Only the taxonomy artifact changed; its provenance recorded 0.8 |
| Restore saved 0.7 taxonomy under the 0.8 filename while retaining current state | Cache rejected for confidence mismatch; only taxonomy regenerated, with recorded confidence 0.8 |

See [taxonomy-resume-audit.json](taxonomy-resume-audit.json) and the
[actual rejection/recomputation log](taxonomy-cache-rejection.txt).
Unit regressions also cover a different classifier, different reads, malformed
provenance, and a non-taxonomy artifact under the expected filename.

## Unit tests and the initial parallel timeout

The first run executed 317 tests concurrently with real workflows and reported
three subprocess.TimeoutExpired errors at the existing 30-second CLI-test bound.
There were no assertion failures or skips. These logs are preserved in
[unit-tests-parallel.txt](unit-tests-parallel.txt).

After all real workflows finished, the same source and environment passed all
317 tests in isolation: zero failures, zero errors, and zero skips. No timeout
limit or assertion was relaxed. Shell syntax checks and actual AmpliconTAXA help
also passed. See [unit-audit-isolated.json](unit-audit-isolated.json),
[unit-tests-isolated.txt](unit-tests-isolated.txt), and
[AmpliconTAXA-help.txt](AmpliconTAXA-help.txt).

## Remaining release work

- Build and publish the compatible fixed GG2 2024.09 and SILVA 138.1 classifiers,
  register their verified identities/checksums, and run the all-classifier gate.
- Build the fresh SIF and validate the packaged runtime and resources.
- Review those results before merging the upgrade branch.

Full reference-model training remains paused. These fixtures check execution,
artifact handling, and caching; they do not assess taxonomic accuracy.
