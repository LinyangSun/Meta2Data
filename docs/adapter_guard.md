# Per-sample fastp inferred-adapter guard

The guard is enabled by default in AmpliconPIP. It checks whether a sequence
inferred as a sequencing adapter may actually be biological 16S sequence, and
adjusts processing for protected samples. `--no-adapter-guard` disables this
additional check; `--adapter-guard` explicitly enables it.

The GG2 backbone is automatically prepared in
`<launch working directory>/results/db/`, shared with TAXA. Valid cached or
bundled references are reused; missing references are downloaded. BLAST indexes
and query results use the fixed cache directory
`<launch working directory>/results/db/adapter_guard/`; indexes are prepared
before parallel dataset workers start.
No classifier or SEPP reference is needed for this check.

Thresholds are in the `adapter_guard` section of `--parameter` JSON. The bundled
E. coli 16S sequence serves amplified-region reporting, independently of this
GG2-based guard.

Each sample first runs the existing adapter-only fastp command. The guard checks
only inferred adapter fields from read ends for which de novo detection was
requested. An inferred sequence of at least 50 nt with a single BLAST HSP at
least 98% identical and covering at least 95% of the query (E-value <= 1e-5)
triggers protection. These conservative engineering thresholds are configurable;
a hit is evidence of potential biological mis-trimming, not a taxonomic call.
An absent inferred adapter, a short query or lack of a strong hit does not prove
that a sample is adapter-free. In particular this guard does not remove adapters
that fastp failed to infer.

A protected single-end sample reruns from the original FASTQ with adapter
trimming disabled. A protected paired-end sample reruns both original mates
without `--detect_adapter_for_pe`, retaining fastp's overlap-based adapter
trimming. Quality and length filtering remain disabled at this adapter stage,
as in the original pipeline; downstream primer and quality steps remain active.
Only the accepted attempt is published into the usual adapter-removed directory.
Names and sample IDs do not change, and PE order and pairing are validated.
The accepted FASTQ reads and bases must agree with its fastp report.

Permanent evidence lives under each read-count execution's
`reports/fastp_guard/<sample>/`: initial JSON/HTML, optional fallback JSON/HTML,
commands, BLAST evidence and `decision.json`. The original `reports/fastp/`
location contains the accepted report. `summary.csv` records only accepted
`fastp_reads`; the discarded initial count remains in the decision sidecar.
Reference and query caches use sequence/reference/tool/parameter fingerprints,
locks and atomic publication so concurrent dataset workers can share them safely.

B4 keeps separate initial-fastp, reference-check and fallback-fastp events with
sample, project, method, wall time, CPU and peak sampled RSS. Each fastp event
captures its own report counts and SHA. Module/dataset/guard parents contain
nested events: do not sum parent and child CPU or wall times. An accepted read
count is not a measure of total compute cost; retries remain real measured work.
Cache hits must be reported when comparing cold and warm execution costs.

Changing guard settings, reference identity or processing code invalidates
processing checkpoints. Unchanged input accession mappings retain the archived
raw downloads. Changed inputs archive old raw files separately so old samples
cannot silently enter a new run.
