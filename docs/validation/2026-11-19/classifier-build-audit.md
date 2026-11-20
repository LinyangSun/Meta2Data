# Compatible classifier build checkpoint

**Status: BOTH MODELS BUILT, AUDITED AND PUBLISHED.** GG2 completed in run04;
SILVA completed in recovery02. Public download and cache-reuse checks passed
for both models. Fresh native integration and final-SIF validation passed.
The initial image's read-only configuration error was corrected before the
successful build and publication. See the [release gate record](release-gates.md).
Training cleanup is complete.

## Environment and inputs

- Active Codespace: 4 CPUs, 16 GiB RAM, standard Linux machine with 32 GB disk.
  The machine choices exposed during this work offered 8 or 16 GiB RAM.
- Frozen environment: QIIME 2, q2-feature-classifier and RESCRIPt **2026.7.0**;
  scikit-learn **1.7.1**.
- GG2 run04 build code: commit `ec56140`, using training batches of **1,000** references.
- GG2 2024.09 full-length training input: **337,506 sequences** and **28,054
  taxonomy classes**. The build retains all reference sequences.
- SILVA 138.1 full-length training input after the recorded preparation:
  **435,518 sequences** and **10,572 taxonomy classes**, with batches of **1,000**.
- All four raw SILVA 138.1 source files were prefetched through the existing
  `build_reference_classifiers.acquire_source` function. Their compressed and
  uncompressed MD5 checks passed. Those checks verified source acquisition;
  completed SILVA training and model validation are recorded below.

Run04 used `/tmp/meta2data-classifier-build-20261119`; the SILVA workspace was
`silva-138.1-45ed4331f34efbfe`. That build-work directory was absent after the
Codespace restart. The GG2 output and manifest remain under
`/workspaces/meta2data-validation/trained-classifiers`.

## GG2 training attempts

| Attempt | Training batch size | Observed result | Elapsed time | Last sampled process RSS | Last sampled available RAM |
|---|---:|---|---|---:|---:|
| run01 | 20,000 | Exit 143 | 3 min 41.512 s | Not recorded | Not recorded |
| run02 | 20,000 | Terminated | 3 min 18.587 s | 9.646 GiB | 4.009 GiB |
| run03 | 1,000 | Terminated | 3 min 50.920 s | 9.561 GiB | 4.029 GiB |
| run04 | 1,000 | GG2 model and build manifest produced | Not recorded | See raw monitor log | See raw monitor log |

Memory observations for run02 and run03 were sampled every five seconds. Their
last recorded RSS values are not measured peak RSS at termination. No kernel
OOM message was found, and the inspected cgroup OOM counters were zero. The
observed cgroup memory peak, **14,281,216,000 bytes**, covers the group rather
than only the training process. The source of SIGTERM has not been established;
these observations do **not** establish a kernel OOM kill.

### Batch size is a computational setting

Reducing `classify__chunk_size` from 20,000 to 1,000 reduces one dense label
matrix for 28,054 classes from approximately **4.18037 GiB to 0.209 GiB**.
This does not remove the full feature matrix or fitted model arrays, and is not
an estimate of total peak memory. It does not limit the number of training
sequences, change the reference database, or subsample the training set.

### Temporary swap for run04

A **16 GiB** temporary swap file was enabled at
`/tmp/meta2data-training-swap.7gnm2h`, owned by root with permissions `0600`.
Observed swap use reached **1.648 GiB** in the final recorded run04 sample.
This observation does not identify the earlier termination cause.

No permanent `fstab` entry or container change was made. The temporary swap
file was absent after the Codespace restart; removal of that file is no longer
pending.

## Run04 output and interruption

The GG2 artifact `gg2-2024.09-full-length-sklearn-1.7.1.qza` and its
`.build.json` manifest were produced. The artifact was shown as approximately
**172 MiB** by `ls`. The manifest records:

| Field | Recorded value |
|---|---|
| `built_utc` | `2026-11-19T18:29:00.948950+00:00` |
| UUID | `b1cd95ec-3416-4183-9286-1f426af96c81` |
| SHA-256 | `b68e1d5632fc3b6ee3fd815d176e5e2273b3f41343858b007e374e83f7f5f5f6` |
| MD5 | `e3a82923c8c412b5c010247974a139c1` |
| Recipe SHA-256 | `1629f9b0fae80f4e10ce66e7ccf8ccf8fde523ec22d486d2c1935b8fc3824893` |

The run04 log stopped during SILVA's `rescript.reverse_transcribe`. Its last
memory sample, at **18:34:06 UTC**, recorded **100% CPU**, **0.614 GiB process
RSS**, **12.93 GiB available RAM** and **1.648 GiB swap use**. No `run04.exit`
record was present. At the **20:53 UTC** check, the Codespace had restarted,
with no live training process. The old `/tmp` build-work and swap file were
absent; model output and logs under `/workspaces` were preserved. The cause
of the interruption has not been established.

## Scoped GG2 audit

The standalone GG2 audit **passed**. Its report is
[`gg2-model-validation.json`](gg2-model-validation.json), copied from
`/workspaces/meta2data-validation/gg2-model-audit-01/model-validation.json`,
with audit source revision `0a103ba240e3fd885d5dd1a5829c98b5f15cec61`.
The artifact SHA-256 matches the build record above.

| Check or field | Recorded result |
|---|---|
| `requested_classifiers` | `["gg2-2024.09"]` |
| Training sequences | 337,506 |
| Training classes | 28,054 |
| `all_class_counts_match` | `true` |
| Classified features | 5 |
| Taxonomy UUID | `aed9345f-1dda-4503-af3b-ad2a26e73d95` |
| Classifier UUID | `b1cd95ec-3416-4183-9286-1f426af96c81` |
| Reads UUID | `10ce232d-c773-4835-a12f-ee4769842489` |
| Confidence threshold | 0.7 |
| Duration | 84.97 seconds |

The audit checked complete training counts, every class count and real
classification with matching output IDs and provenance. GG2 reference inputs
were restored through the same pinned source-acquisition functions under
`/workspaces/meta2data-validation/classifier-build/gg2-2024.09-1629f9b0fae80f4e/sources`.

This report applies to GG2 only. SILVA passed a separate scoped audit during
recovery02, recorded below. Both reports explicitly record
`requested_classifiers`; the audit driver's default `all` checks both models.

## Published model resources

The GitHub Release
[`classifiers-2026.7-sklearn-1.7.1-v1`](https://github.com/LinyangSun/Meta2Data/releases/tag/classifiers-2026.7-sklearn-1.7.1-v1)
was published at **2026-11-19T22:41:14Z**, with `draft=false` and seven uploaded
assets: two models, both build manifests, both scoped validation reports and
`reference-model-licenses.txt`. GitHub's recorded SHA-256 values for both
models match their build manifests. Public release metadata is retained in
[`classifier-release.json`](classifier-release.json).

| Model download | Bytes |
|---|---:|
| [GG2 2024.09](https://github.com/LinyangSun/Meta2Data/releases/download/classifiers-2026.7-sklearn-1.7.1-v1/gg2-2024.09-full-length-sklearn-1.7.1.qza) | 180,058,060 |
| [SILVA 138.1](https://github.com/LinyangSun/Meta2Data/releases/download/classifiers-2026.7-sklearn-1.7.1-v1/silva-138.1-full-length-sklearn-1.7.1.qza) | 218,310,928 |

Registry URLs, UUIDs and SHA-256 pins were compared with the committed build
manifests. Public download from an empty cache and subsequent reuse **passed**
for both models at revision `ff66625`. Each model downloaded once; repeat
preparation preserved its UUID, SHA-256, size, modification time and inode
without another download. Evidence is in
[`public-model-download-audit.json`](public-model-download-audit.json).
Draft-stage release metadata was saved in
`/workspaces/meta2data-validation/logs/classifier-draft-assets.json`.

### Published-resource regression checks

All **13 resource tests passed**. The first TAXA CLI run had **12 failures
among 22 tests** because its test shim read the removed `CLASSIFIER_ID`
environment variable and raised `KeyError`. The probe was corrected to read
the classifier identity from generated state, preserving selector, alias,
filename and cache assertions. All **22 tests then passed in 5.525 seconds**
in the frozen environment at revision `f2f7762`. No production code or timeout
was changed for this fix.

The summary is in [`published-resource-checks.json`](published-resource-checks.json).
Raw records are retained for the
[resource tests](published-resource-unit-tests.txt),
[initial CLI run](published-taxa-cli-initial.txt) and
[successful CLI rerun](published-taxa-cli-rerun.txt).

## SILVA recovery02 completion

A detached SILVA-only recovery started at **21:09:11 UTC** in the existing
Codespace using `run-silva-recovery-02.sh`, with build work under
`/workspaces/meta2data-validation/classifier-build`. The process was confirmed
alive during `rescript.reverse_transcribe`, after all four imports and
taxonomy parsing. Its launcher records
`logs/silva-recovery-02.log` and `logs/silva-recovery-02.exit` beneath
`/workspaces/meta2data-validation`, and runs the SILVA model audit after a
successful build. The build completed in **87 min 46.360 s**. The full wrapper,
including the model audit, exited **0** at **2026-11-19T22:37:32Z**.

Recovery01 exited with status 127 before building because `/usr/bin/time` was
absent. Its log is retained. Recovery02 uses Bash's built-in `time` instead.
Neither attempt changes the training recipe or limits the input sequences.

Completed intermediates, inputs, models and logs now reside on the persistent
`/workspaces` mount; QIIME temporary extraction remains on the separate `/tmp`
mount. A temporary 16 GiB swap file at
`/tmp/meta2data-training-swap.hXFB9Y` was enabled with root ownership and mode
0600. No `fstab` entry was added. This swap was disabled and removed during
the cleanup recorded below.
`nohup` detaches the process from the terminal; it cannot preserve a running
process across Codespace shutdown. Valid completed intermediates can be reused
after a restart.

No new Codespace was created, and no account timeout setting was changed.

## SILVA model and scoped audit

The completed build is recorded in
[`silva-classifier-build.json`](silva-classifier-build.json). Its full-length
model passed the audit in
[`silva-model-validation.json`](silva-model-validation.json), copied from
`/workspaces/meta2data-validation/silva-model-audit-02/model-validation.json`.
The audit used source revision `d0b7a38f28afb1c5e4ca77c3a547dfd33f8eb7b9`
and requested only `silva-138.1`.

| Field or check | Recorded result |
|---|---|
| Model filename | `silva-138.1-full-length-sklearn-1.7.1.qza` |
| `built_utc` | `2026-11-19T22:36:47.835193+00:00` |
| Classifier UUID | `7bb8da5e-bc83-4d33-8271-97f3c4028973` |
| SHA-256 | `81aa6ceef0fda3cef4e988423832eaeb8c134b67f300f5ba07ea3ac01e1da1d1` |
| MD5 | `2cfc9c00b994e5a6dd58c251c91db7c0` |
| Recipe SHA-256 | `45ed4331f34efbfee4e9540c31c99d8c76ba7ce726d2fd8e0a313ea878947bc6` |
| Training sequences | 435,518 |
| Training classes | 10,572 |
| `all_class_counts_match` | `true` |
| Classified features | 5 |
| Taxonomy UUID | `b6eaea05-c748-4fe0-a8c4-b3827defe2d3` |
| Reads UUID | `10ce232d-c773-4835-a12f-ee4769842489` |
| Confidence threshold | 0.7 |
| Audit duration | 29.61 seconds |

The audit matched the total and every class count against all prepared
reference sequences, then checked real classification, output IDs and
provenance. Its reference artifacts are retained under
`/workspaces/meta2data-validation/classifier-build/silva-138.1-45ed4331f34efbfe/`.
These two scoped model audits are separate from the subsequently completed
[native and final-SIF release gates](release-gates.md).

## Training cleanup

Cleanup completed at **2026-11-19T22:43:00.079695+00:00**. The temporary swap
was disabled and removed; `swapon --show` was empty. The temporary `py-spy`
installation and three raw stack snapshots were removed. Models, prepared
reference intermediates and raw logs were retained; monitors and jobs had
exited. Evidence is retained in [`training-cleanup.json`](training-cleanup.json).

## Raw evidence locations

The original build and memory logs are under
`/workspaces/meta2data-validation/logs/`:

```text
classifier-build-20261119-run01.log
classifier-build-20261119-run02.log
classifier-build-20261119-run03.log
classifier-build-20261119-run04.log
classifier-memory-run02.jsonl
classifier-memory-run03.jsonl
classifier-memory-run04.jsonl
```

The memory logs begin with a `MONITOR_PID` line. Preserve them as raw text when
archiving; they are not strict JSON Lines files from the first line onward.

## Completed release checks

The final production-change audit passed with no findings.

README, help, resource documentation and release status describe the verified
implementation. The fresh native and final-SIF gates passed on the same source
revision. The validated image is published as `v1.1`; image identity and the
complete validation scope are recorded in the [release gate record](release-gates.md).
