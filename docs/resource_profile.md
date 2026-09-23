# Optional resource profiling

Enable profiling by setting an absolute `M2D_PROFILE_DIR` before invoking
AmpliconPIP or AmpliconTAXA. `M2D_EXECUTION_ID` identifies the shared run;
`M2D_PROFILE_PROJECT` and `M2D_PROFILE_METHOD` can label a project invocation.
`M2D_PROFILE_INTERVAL` defaults to one second. With Apptainer `--cleanenv`, pass
these variables explicitly using `--env`.

```bash
M2D_PROFILE_DIR=/absolute/run/b4 M2D_EXECUTION_ID=run001 \
  Meta2Data AmpliconPIP --vsearch -m metadata.csv \
  --col-bioproject Bioproject --col-sra Run -o results --max-parallel 2 -t 8
python3 scripts/resource_profile.py export --directory /absolute/run/b4
```

Native AmpliconPIP workers automatically receive a separate `dataset` scope,
with stage `pip_dataset_total`, for each BioProject. With `-t 8 --max-parallel 2`,
at most two datasets run concurrently and each receives four CPU threads. The
worker context includes its BioProject, thread budget and parent invocation ID;
its child commands inherit those values. Actual worker failures are reported
with a nonzero PIP exit status after other workers finish and the summary is saved.

The B3/B4 submission script additionally wraps complete module invocations with
`resource_profile.py run --scope module --stage pip_batch_total|taxa_total -- CMD`.
These scopes observe shell orchestration, unwrapped child tools and waiting as
well as the scientific commands. A unique invocation identifier is propagated to nested
commands. Supported executables are wrapped through a private PATH directory,
so Python subprocess calls using PATH are observed as well. Absolute executable
paths are covered by their enclosing command's accounting, without a separate
row. Arguments, standard streams, environment and exit codes are preserved.
SIGINT/SIGTERM/SIGHUP are forwarded to the child process group.

## Durable outputs

- `events/ID.start.json` is written before execution, `events/ID.json` on exit,
  and `events/ID.samples.jsonl.gz` holds flushed process/cgroup observations.
- `step_resources.csv` contains one row per invocation, including failed work.
  Unfinished start records are exported as `incomplete`, with unavailable end
  metrics left NA. Running jobs should not be interpreted from a final export.
- `dataset_resources.csv` selects inclusive native dataset worker rows, plus
  historical single-project `pip_total` module rows. `pipeline_resources.csv`
  selects complete PIP invocations, including multi-project batches;
  `taxa_resources.csv` selects complete TAXA module invocations. The batch wall
  time is not the sum of overlapping dataset wall times.
- Original read-count `summary.csv` and TAXA tables keep their existing names.
  Read-count events gain `execution_id` and `resource_invocation_id` pointers to
  the exact auditing command. Permanent `sample_sizes.csv` files are written in
  the raw-count attempt directory while FASTQ is already being traversed.

The benchmark assembler copies actual input measurements into a run-level
`sample_sizes.csv`, records raw checksums, exports `artifact_workloads.csv`
(feature/sample counts and representative-sequence lengths), and saves
`read_event_links.csv`. Existing read-count histories remain the authority.
Where the same attempt, explicit stage mapping and next matching audit allow
it, resource rows contain input/output read counts and the evidence path.
Unmapped stages, failed commands and inherited observations do not acquire
fabricated throughput denominators. For per-sample fastp/cutadapt commands,
the actual output path must also match the audited output directory.

## Measurement interpretation

- `wall_s` uses a monotonic clock around the command, including child startup
  and waiting. The outer benchmark also records Apptainer invocation timestamps.
- CPU user/system seconds come from Linux `wait4`, including descendants reaped
  by the command. This catches CPU used by short-lived children between samples.
  Detached/orphan processes are outside that guarantee; the scientific commands
  must wait for their workers. **Do not sum parents and descendants.** Both
  levels are kept so the selected accounting scope remains inspectable.
- `allocated_cpus` is the configured scope budget (four for each dataset in an
  eight-CPU/two-dataset run); `slurm_allocated_cpus` retains the complete task
  allocation. These are capacity settings, not measured CPU seconds.
- `process_maxrss_kib` is the Linux exit-accounting high-water statistic, not
  the sum of simultaneous child peaks. Sampled process-tree RSS can double-count
  shared pages and miss brief peaks. Memory observations are estimates.
- All `job_cgroup_*` values describe the enclosing job/step, including the
  profiler, sibling processes and charged page cache. They are contextual
  observations, not isolated per-command measurements. Lifetime `memory.peak`
  is never reused as a stage peak. No cgroup or host cache counters are reset.
- `rusage_block_*_bytes` are Linux block-I/O accounting (blocks x 512);
  `sampled_process_*_bytes` are sampled command-root I/O lower bounds, including
  waited-for children. Historical child counters are not added again. They are not
  network traffic or FASTQ sizes. Unsupported observations are NA.
- Input/output byte counts cover explicitly named files only. No repeated
  recursive directory scan is performed. Artifact/input inventories are separate
  boundary measurements. The profiler does not claim exact peak scratch usage.
- Download, metadata, auditing, reporting and compute have separate category
  labels. Cached stages are represented by their audit/checkpoint evidence and
  absent computation, never an invented zero-second denoising run.

Profiling is opt-in and requires `psutil` (already available in the production
SIF). Small-fixture off/on overhead is measured during validation; its order/cache
effects are reported and no guessed correction is subtracted from real timings.

References: [Python wait4](https://docs.python.org/3/library/os.html#os.wait4),
[Linux process I/O](https://man7.org/linux/man-pages/man5/proc_pid_io.5.html),
[Slurm sacct](https://slurm.schedmd.com/sacct.html).

## Adapter protection attempts

With `--adapter-guard`, `fastp_initial` and `fastp_fallback` each record their own
wall/CPU/memory cost and a permanent `events/<invocation>.fastp-report.json`
snapshot. `fastp_report_sha256`, `fastp_report_status`, `input/output_individual_reads`,
`input/output_reads`, bases and average lengths identify the exact attempt. PE
normalization uses `read_pairs` only for a consistent paired report. At export,
`fastp_attempt_decision` links accepted/rejected attempts to the permanent
`decision_report`; unverified or incomplete decisions remain unresolved.

`fastp_adapter_guard`, `reference_blast` and `blast_database_build` are also
measured. BLAST input sizes include its index-file prefix. `M2D_PROFILE_STAGE`
is consumed by the next recorded command; sample and decision context persist
into children. Parent costs include children, so use either parent totals or
non-overlapping tool events when comparing efficiency. Refer to
[adapter_guard.md](adapter_guard.md) for scientific scope and threshold limits.
