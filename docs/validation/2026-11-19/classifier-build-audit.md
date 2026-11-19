# Compatible classifier build checkpoint

**Status: PARTIAL.** Run04 produced the GG2 model and its build manifest, and
the separate GG2 model audit has passed. SILVA preparation stopped at
`rescript.reverse_transcribe`; its replacement model has not been produced.
Neither model has been published, and the upgrade branch has not been merged.

## Environment and inputs

- Active Codespace: 4 CPUs, 16 GiB RAM, standard Linux machine with 32 GB disk.
  The machine choices exposed during this work offered 8 or 16 GiB RAM.
- Frozen environment: QIIME 2, q2-feature-classifier and RESCRIPt **2026.7.0**;
  scikit-learn **1.7.1**.
- Build code: commit `ec56140`, using training batches of **1,000** references.
- GG2 2024.09 full-length training input: **337,506 sequences** and **28,054
  taxonomy classes**. The build retains all reference sequences.
- All four raw SILVA 138.1 source files were prefetched through the existing
  `build_reference_classifiers.acquire_source` function. Their compressed and
  uncompressed MD5 checks passed. This verifies acquisition, not completion of
  SILVA preparation or model training.

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

This pass applies to GG2 only. Reports explicitly record
`requested_classifiers`; the audit driver's default `all` checks both models.
The automatic audit scheduled after both builds has not run.

## Draft asset backup

The draft GitHub Release `classifiers-2026.7-sklearn-1.7.1-v1` contains the
GG2 model, its build manifest, the scoped validation report and
`reference-model-licenses.txt`. All four assets report `uploaded`.
The model is **180,058,060 bytes**; GitHub's recorded SHA-256 matches
`b68e1d5632fc3b6ee3fd815d176e5e2273b3f41343858b007e374e83f7f5f5f6`.
Release metadata is saved in
`/workspaces/meta2data-validation/logs/classifier-draft-assets.json`.
This is a draft backup, not a published model download or completed source
installation test. The registry URLs remain pending until the bundle is ready.

## SILVA recovery launch

A detached SILVA-only recovery started at **21:09:11 UTC** in the existing
Codespace using `run-silva-recovery-02.sh`, with build work under
`/workspaces/meta2data-validation/classifier-build`. The process was confirmed
alive during `rescript.reverse_transcribe`, after all four imports and
taxonomy parsing. Its launcher records
`logs/silva-recovery-02.log` and `logs/silva-recovery-02.exit` beneath
`/workspaces/meta2data-validation`, and runs the SILVA model audit after a
successful build. The model has not yet been produced or audited.

Recovery01 exited with status 127 before building because `/usr/bin/time` was
absent. Its log is retained. Recovery02 uses Bash's built-in `time` instead.
Neither attempt changes the training recipe or limits the input sequences.

Completed intermediates, inputs, models and logs now reside on the persistent
`/workspaces` mount; QIIME temporary extraction remains on the separate `/tmp`
mount. A temporary 16 GiB swap file at
`/tmp/meta2data-training-swap.hXFB9Y` was enabled with root ownership and mode
0600. No `fstab` entry was added. Remove this swap after training and auditing.
`nohup` detaches the process from the terminal; it cannot preserve a running
process across Codespace shutdown. Valid completed intermediates can be reused
after a restart.

No new Codespace was created, and no account timeout setting was changed.

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

## Required remaining checks

1. Complete the SILVA model build and
   retain both models' `.build.json` records and required reference inputs.
2. Run [verify-rebuilt-models.py](audit-drivers/verify-rebuilt-models.py) with
   `--classifier silva-138.1`. The SILVA audit must verify model identity and
   checksums, exact total and per-class training counts against all actual
   training references, and real classification of the existing representative
   sequences with matching output IDs and provenance. The scoped GG2 audit
   has already passed.
3. Publish the verified versioned assets, register their actual download URLs,
   UUIDs and checksums, and verify download from an empty classifier cache and
   subsequent valid cache reuse.
4. Complete the final code simplicity review and synchronize README, help,
   resource documentation and release status with the verified implementation.
5. Pass the existing native and final-SIF release gates, publish the validated
   image tag, and only then complete the merge. Earlier workflow results do not
   substitute for successful execution of these remaining gates.
6. Record the final build outcome, retaining the unsuccessful attempts and
   interruption evidence.
