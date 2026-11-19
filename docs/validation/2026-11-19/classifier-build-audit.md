# Compatible classifier build checkpoint

**Status: PARTIAL.** Run04 produced the GG2 model and its build manifest, but
the model audit has not run. SILVA preparation stopped at
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

These are build records, not an independent model audit. No model audit has
run, including the automatic audit scheduled after both builds.

The run04 log stopped during SILVA's `rescript.reverse_transcribe`. Its last
memory sample, at **18:34:06 UTC**, recorded **100% CPU**, **0.614 GiB process
RSS**, **12.93 GiB available RAM** and **1.648 GiB swap use**. No `run04.exit`
record was present. At the **20:53 UTC** check, the Codespace had restarted,
with no live training process. The old `/tmp` build-work and swap file were
absent; model output and logs under `/workspaces` were preserved. The cause
of the interruption has not been established.

## Raw evidence locations

All logs are under `/workspaces/meta2data-validation/logs/`:

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

1. Complete the SILVA model build, retain both `.build.json` records, and
   restore the reference inputs needed for the model audits.
2. Run [verify-rebuilt-models.py](audit-drivers/verify-rebuilt-models.py). This
   audit driver is ready and reviewed but has **not run**. It must verify model
   identity and checksums, exact total and per-class training counts against
   all actual training references, and real classification of the existing
   representative sequences with matching output IDs and provenance.
   `--classifier gg2-2024.09` permits a separate GG2 audit; the default `all`
   audits both models. Reports explicitly record `requested_classifiers`, so
   a passed scoped audit applies only to the requested model.
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
