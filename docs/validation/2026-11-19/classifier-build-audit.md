# Compatible classifier build checkpoint

**Status: RUNNING.** At this checkpoint, neither replacement model has been
produced or published. The upgrade branch has not been merged. The observations
below describe an ongoing build, not a successful release validation.

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

Build work is under `/tmp/meta2data-classifier-build-20261119`; the SILVA
workspace is `silva-138.1-45ed4331f34efbfe`. Finished models are intended for
`/workspaces/meta2data-validation/trained-classifiers`.

## GG2 training attempts

| Attempt | Training batch size | Observed result | Elapsed time | Last sampled process RSS | Last sampled available RAM |
|---|---:|---|---|---:|---:|
| run01 | 20,000 | Exit 143 | 3 min 41.512 s | Not recorded | Not recorded |
| run02 | 20,000 | Terminated | 3 min 18.587 s | 9.646 GiB | 4.009 GiB |
| run03 | 1,000 | Terminated | 3 min 50.920 s | 9.561 GiB | 4.029 GiB |
| run04 | 1,000 | Still running beyond the earlier failure point | Ongoing | See raw monitor log | See raw monitor log |

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
Observed swap use reached **1.535 GiB** while run04 remained active. Swap use
and survival beyond previous failure times do not yet demonstrate successful
training or identify the earlier termination cause.

No permanent `fstab` entry or container change was made. **Disable this swap
and remove this specific file after training ends and memory permits.** Cleanup
remains outstanding at this checkpoint.

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

1. Complete both full model builds and retain their `.build.json` records.
2. Run [verify-rebuilt-models.py](audit-drivers/verify-rebuilt-models.py). This
   audit driver is ready and reviewed but has **not run**. It must verify model
   identity and checksums, exact total and per-class training counts against
   all actual training references, and real classification of the existing
   representative sequences with matching output IDs and provenance.
3. Publish the verified versioned assets, register their actual download URLs,
   UUIDs and checksums, and verify download from an empty classifier cache and
   subsequent valid cache reuse.
4. Complete the final code simplicity review and synchronize README, help,
   resource documentation and release status with the verified implementation.
5. Pass the existing native and final-SIF release gates, publish the validated
   image tag, and only then complete the merge. Earlier workflow results do not
   substitute for successful execution of these remaining gates.
6. Finish temporary swap cleanup and record the final build outcome, including
   unsuccessful attempts rather than replacing their evidence.
