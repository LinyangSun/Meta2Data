# AmpliconPIP Pipeline Flow

Detailed workflow of `Meta2Data AmpliconPIP`. The orchestration lives in
[scripts/run.sh](../scripts/run.sh); platform-specific steps are
implemented as functions in
[scripts/AmpliconFunction.sh](../scripts/AmpliconFunction.sh).

```mermaid
flowchart TD
    Start([Meta2Data AmpliconPIP<br/>-m metadata.csv -t N]):::terminal --> PP1

    PP1[<b>PHASE 1</b><br/>Parse metadata CSV<br/>→ datasets_ID.txt + per-dataset SRA lists]:::phase
    PP1 --> PP15

    PP15[<b>PHASE 1.5</b><br/>Batch Entrez platform query<br/>→ .platform_cache.txt]:::phase
    PP15 --> PP2

    PP2[<b>PHASE 2</b><br/>Per-dataset loop<br/>max-parallel concurrency]:::phase --> DISPATCH

    DISPATCH{Platform?}:::decision
    DISPATCH -->|ILLUMINA| I1
    DISPATCH -->|LS454| L1
    DISPATCH -->|ION_TORRENT| T1
    DISPATCH -->|PACBIO_SMRT| B1
    DISPATCH -->|OXFORD_NANOPORE| SKIP1[Skip → failed_log]:::skip

    %% ====== ILLUMINA branch ======
    I1[Download SRA]:::illumina --> I2[Detect PE/SE layout]:::illumina
    I2 --> I3[fastp adapter removal]:::illumina
    I3 --> I4[entropy_primer_detect.py]:::illumina
    I4 --> IQ{q-score<br/>diverse?}:::decision

    IQ -->|degraded<br/>1 unique value| ID1[degraded_quality_preprocess<br/>trim + length + N filter]:::illumina
    ID1 --> ID2[vsearch dereplicate]:::illumina
    ID2 --> ID3[vsearch UNOISE3 denoise]:::illumina
    ID3 --> ID4[vsearch usearch_global<br/>map reads → ZOTUs]:::illumina
    ID4 --> ID5[Import to QIIME2]:::illumina
    ID5 --> ID6[Filter low-freq OTUs]:::illumina
    ID6 --> ICLEAN[FinalFilesCleaning]:::illumina

    IQ -->|normal| IN1[Import to QIIME2]:::illumina
    IN1 --> IN2[DADA2 denoise-paired/single]:::illumina
    IN2 --> IN3{PE retention<br/>≥ 50%?}:::decision
    IN3 -->|yes| ICLEAN
    IN3 -->|no, PE mode| IN4[Re-import as SE<br/>Re-run DADA2 denoise-single]:::illumina
    IN4 --> IN5{SE retention<br/>≥ 50%?}:::decision
    IN5 -->|yes| ICLEAN
    IN5 -->|no| IN6[Mark low_quality<br/>preserve tmp/]:::skip
    IN6 --> ICLEAN

    ICLEAN --> SUMMARY

    %% ====== LS454 branch ======
    L1[Download SRA]:::ls454 --> L2[fastp adapter removal]:::ls454
    L2 --> L3[entropy_primer_detect.py]:::ls454
    L3 --> L4[adaptive_tail_trim<br/>data-driven 3' N removal]:::ls454
    L4 --> L5[Import to QIIME2]:::ls454
    L5 --> L6[quality-filter q-score<br/>--min-quality 0 + length]:::ls454
    L6 --> L7[vsearch dereplicate]:::ls454
    L7 --> L8[vsearch uchime-denovo]:::ls454
    L8 --> L9[cluster-features-de-novo 97%]:::ls454
    L9 --> L10[filter-features min-freq 2]:::ls454
    L10 --> LCLEAN[FinalFilesCleaning]:::ls454
    LCLEAN --> SUMMARY

    %% ====== ION_TORRENT branch ======
    T1[Download SRA]:::iontorrent --> T2[fastp adapter removal]:::iontorrent
    T2 --> T3[entropy_primer_detect.py]:::iontorrent
    T3 --> T4[Import to QIIME2]:::iontorrent
    T4 --> T5[quality-filter q-score]:::iontorrent
    T5 --> T6[DADA2 denoise-pyro<br/>--p-trim-left 10]:::iontorrent
    T6 --> TCLEAN[FinalFilesCleaning]:::iontorrent
    TCLEAN --> SUMMARY

    %% ====== PACBIO branch ======
    B1[Download SRA]:::pacbio --> B2[fastp adapter removal<br/>no primer detection]:::pacbio
    B2 --> B3[Sample 1000 reads<br/>check > 1400bp ratio]:::pacbio
    B3 --> BQ{> 50% long<br/>reads?}:::decision
    BQ -->|no| BSKIP[Skip → failed_log]:::skip
    BQ -->|yes| B4[Read 27F/1492R primers<br/>from docs/]:::pacbio
    B4 --> B5[Import to QIIME2]:::pacbio
    B5 --> B6[quality-filter q-score]:::pacbio
    B6 --> B7[DADA2 denoise-ccs<br/>with 27F/1492R]:::pacbio
    B7 --> B8[Extract reads]:::pacbio
    B8 --> BCLEAN[FinalFilesCleaning]:::pacbio
    BCLEAN --> SUMMARY

    %% ====== Aggregation ======
    SUMMARY[Append summary.csv row<br/>+ success/failed/skipped/low_quality logs]:::phase --> FIN([Final Summary<br/>N datasets complete]):::terminal

    %% ====== Styling ======
    classDef phase fill:#e1f5ff,stroke:#0277bd,stroke-width:2px,color:#000
    classDef decision fill:#fff3e0,stroke:#e65100,stroke-width:2px,color:#000
    classDef terminal fill:#f3e5f5,stroke:#6a1b9a,stroke-width:2px,color:#000
    classDef skip fill:#ffebee,stroke:#c62828,stroke-width:1px,color:#000
    classDef illumina fill:#e8f5e9,stroke:#2e7d32,color:#000
    classDef ls454 fill:#fff8e1,stroke:#f57f17,color:#000
    classDef iontorrent fill:#fce4ec,stroke:#ad1457,color:#000
    classDef pacbio fill:#e0f2f1,stroke:#00695c,color:#000
```

## Legend

| Color | Meaning |
|---|---|
| 🔵 Phase (blue) | Top-level orchestration stages run by [run.sh](../scripts/run.sh) |
| 🟠 Decision (orange) | Runtime branching based on data characteristics |
| 🟢 ILLUMINA (green) | Illumina pipeline, supports both DADA2 and VSEARCH sub-paths |
| 🟡 LS454 (yellow) | 454 pipeline, always VSEARCH OTU clustering |
| 🩷 ION_TORRENT (pink) | Ion Torrent pipeline, DADA2 `denoise-pyro` with 10bp trim |
| 🟦 PACBIO_SMRT (teal) | PacBio pipeline, full-length 16S CCS only (`denoise-ccs`) |
| 🔴 Skip (red) | Dataset skipped or marked low-quality; preserved for debugging |

## Cross-cutting mechanisms not shown in the diagram

1. **Checkpoint / Resume** — Each branch checks whether
   `tmp/step_02_fastp/` (or PacBio's `tmp/step_01_adapter_removed/`)
   contains the expected number of FASTQ files for the dataset's SRA
   list. If it matches, the branch resumes from that checkpoint and
   skips download + pre-processing. If it does not, `tmp/` and
   `ori_fastq/` are wiped and the branch restarts from scratch.

2. **Parallel execution** — The diagram shows the flow for a single
   dataset; in practice `PHASE 2` runs up to `--max-parallel` datasets
   concurrently, each in its own log file. A dedicated file descriptor
   (`fd 3`) carries milestone messages ("Dataset X: success / failed")
   to the main console while per-dataset stdout goes to
   `logs/<dataset_id>.log`.

3. **Quality-degraded dispatch is cross-platform** — The diagram draws
   the `q-score diverse?` decision inside the ILLUMINA branch, but the
   underlying `check_quality_diversity` helper
   ([run.sh:517](../scripts/run.sh#L517)) is generic. ILLUMINA is the
   only branch where it actually flips the pipeline to VSEARCH at
   runtime, because 454 is hard-coded to VSEARCH and Ion
   Torrent / PacBio have usable q-scores in practice.

4. **Platform detection has two layers** — `PHASE 1.5` does a single
   batched Entrez query for all datasets up-front and caches the
   result in `.platform_cache.txt`, avoiding NCBI's 3 req/s rate limit
   during parallel execution. `PHASE 2` reads from the cache first and
   only falls back to a live API query on cache miss.
