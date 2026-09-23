# Permanent read-count outputs

AmpliconPIP extends the existing **`summary.csv`**. It does not create a replacement
`pip_read_counts.csv`. Each row is a BioProject / Run / sample / method combination.
The identifier columns come first, then `RawReads`, `fastp_reads`,
`primer_trimmed_reads`, `sanitized_reads`, the ordered `dada2_*` and `vsearch_*`
columns, and the backwards-compatible `FinalReads`. The other branch's columns
are `NA`. `FinalReads` equals the selected branch's final feature-table abundance.
`pip_dataset_read_counts.csv` has the corresponding BioProject totals.

- `0`: the step completed and retained zero reads, including samples absent from
  the resulting feature table. Such samples are never silently removed from the CSV.
- `NA`: the step was not executed, has not completed, or cannot provide this
  sample's read abundance. It is not a zero and must not enter a loss calculation.
- `CountUnit`: original SE inputs use `reads`; original PE inputs use `read_pairs`.
  Each pair contributes one starting fragment, not two. Later merged reads or
  forward-only reads are comparable to this fragment baseline. The legacy
  `<dataset>_raw_read_counts.tsv` continues to contain R1+R2 individual-read totals.
- `ProcessingPath`: DADA2's selected `paired`/`single` execution, or VSEARCH's
  per-sample preprocessing input path (merged FASTQ versus original forward FASTQ).
  Detailed events record each observation's source, basis and predecessor.
- `Status`: `running`, `success`, `failed`, or `skipped`. A partially completed
  dataset has only the steps actually observed; downstream values remain `NA`.

## PIP steps

DADA2 records optional external Q-score filtering, denoiser input, CCS primer
matching/removal when reported, quality filtering, denoising, pair merging, and
non-chimeric output. The actual plugin statistics determine which columns apply.
Some versions expose one `denoised` count; versions exposing `denoised-f` and
`denoised-r` populate separate directional columns. Percentages are not counts.
All original columns remain available in the archived native statistics file.

VSEARCH records pair-merge attempts, the chosen preprocessing input, quality or
combined length/N filtering, adaptive trimming or chopper filtering where used,
preprocessed reads, dereplication, abundance filtering, preclustering, UNOISE,
chimera removal, OTU clustering, mapping input, mapped abundance, imported table,
and final low-frequency-filtered table. ONT also has denoising, polishing and
relabeling columns. A combined tool operation has one output count; this does not
invent separate intermediate outputs for simultaneous filter criteria.

**Pooled VSEARCH stages:** the BioProject table records the sum of `;size=`
annotations, not the number of representative sequences. Per-sample columns at
these stages are `NA` because pooled centroids do not carry unambiguous sample
assignments. The permanent header evidence retains sequence IDs and sizes. Feature
numbers are separately recorded in the events and are never inserted as read
counts. If ONT polishing loses abundance annotations, affected counts remain `NA`
with an explicit reason, including downstream stages whose default size=1 would
otherwise be mistaken for original-read abundance.

**Mapping restarts from all preprocessed reads**, not just surviving centroids.
Therefore `vsearch_mapped_reads` can exceed `vsearch_clustered_reads`. Mapping loss
is calculated against `vsearch_mapping_input_reads`; differences between arbitrary
adjacent columns are not necessarily read losses. SE fallback similarly has its
own input. A successful merge attempt yielding zero reads is recorded as zero;
a failed tool invocation has no invented count.

## TAXA steps

Within each `final-<method>-<region>/` output directory:

- `taxa_read_counts.csv`: per-sample retained abundance, beginning with
  `taxa_raw_reads` from the merged input feature table, followed by
  `taxa_oriented_reads`, `taxa_classified_reads`, `taxa_tree_placed_reads`,
  `taxa_final_reads`.
- `taxa_total_read_counts.csv`: matching combined totals.
- `taxa_read_losses.csv`: input/output counts, lost counts and loss fractions for
  each comparable step and sample, plus a `__dataset_total__` row per step.

The initial sample universe and BioProject attribution are checked against the
collected dataset tables. Taxonomic annotation does not remove unclassified reads,
so classification retains the oriented table's counts. Single-region and `--notree`
runs do not execute SEPP (`taxa_tree_placed_reads=NA`). In `--notree` mode,
`taxa_final_reads` is the oriented/classified table abundance; no placement filter
is applied. The tree mode is also preserved in the attempt's `taxa-run-state.json`.
Alignment/masking/tree construction
operate on representatives and do not constitute additional read-table filters.
Multi-region runs check that retained plus SEPP-unplaced abundance equals the
oriented input, including an empty removed table. Samples with total dropout remain
as zero rows. TAXA counts are **feature-table abundance**, not ASV/OTU richness.

## Persistence, retries and cleanup

Counts are committed after each completed operation. CSV replacement uses an
atomic rename; shared PIP CSV updates are protected by a file lock and keyed by
dataset/method. Re-running one branch preserves the other branch's rows. The
per-dataset region summary also selects the corresponding method.

Each dataset (or TAXA final directory) keeps:

```text
read_counts/<method>/<attempt>/
  state.json                 # latest state of this attempt
  event-00001.json ...        # immutable observations, sources and predecessor edges
  summary.csv                # this attempt's per-sample view
  dataset_summary.csv        # this attempt's aggregate view
  read_losses.csv            # comparable input/output edges only
  reports/fastp/<sample>.*    # distinct per-sample JSON and HTML reports (PIP)
  dada2-*/stats.tsv           # native statistics and source QZA, one set per invocation
  qc-*/stats.csv              # native external QC statistics and source QZA, if used
  *.headers.txt              # pooled abundance evidence, independent of tmp
```

`latest.json` points to the latest attempt's state. A PE-to-SE retry archives both
DADA2 invocations, then replaces the selected DADA2 columns and clears the obsolete
merge count. Valid checkpoints inherit only the appropriate upstream observations,
with the originating attempt recorded and the input/configuration fingerprint
checked. Fresh executions clear current rows, preserving previous attempt files.
A restart that skips an already completed dataset repairs its summary from the
retained ledger and final table. Historical runs without this ledger cannot recover
intermediate counts that were already deleted; no counts are fabricated.

The statistics and tool reports are outside `tmp` and survive normal cleanup.
DADA2 statistics export now recognizes `stats.tsv`; exported QC and DADA2 tables
are normalized to TSV. The existing `summary.csv` filename and `RawReads` /
`FinalReads` names are preserved. Shell helper calls outside the PIP runner do not
create a ledger unless initialized by the runner.

Validation and example outputs: `validation/read-counts-20260921/`.
The pre-existing SIF is not rebuilt by this source-code change; the validation uses
its installed tools to execute the revised scripts from the host source directory.
