# Total RNA-seq Snakemake Workflow

End-to-end, reproducible RNA-seq pipeline driven by a metadata sheet.

```
FastQC ─► TrimGalore ─► HISAT2 ─► samtools ─► featureCounts
                                                  │
                              process_fc (Python) │
                                                  ▼
                      DESeq2 (R): PCA + per-contrast DEGs + volcanoes
                                                  │
                                              MultiQC
```

Supports **paired-end** and **single-end** libraries via a single config switch.
Conda environments are split into `tools` / `python` / `r` so each step gets a
clean, reproducible interpreter.

---

## Repository layout

```
.
├── Snakefile                  # entry point, rule all, wildcard constraints
├── config/
│   ├── config.yaml            # paths, PE/SE, tool params, thresholds
│   ├── samples.tsv            # sample, rep, condition, time
│   └── contrasts.tsv          # named DESeq2 comparisons
├── profiles/
│   └── slurm/
│       └── config.yaml        # cluster execution profile (you create)
└── workflow/
    ├── rules/
    │   ├── qc.smk             # FastQC + MultiQC
    │   ├── trim.smk           # TrimGalore (PE or SE branch)
    │   ├── align.smk          # HISAT2, samtools view/sort/index/flagstat
    │   ├── counts.smk         # featureCounts + tidy-up
    │   └── deseq2.smk         # PCA + per-contrast DE/volcano
    ├── envs/
    │   ├── tools.yaml         # fastqc, trim-galore, hisat2, samtools, subread, multiqc
    │   ├── python.yaml        # python 3.11 + pandas
    │   └── r.yaml             # R 4.3 + DESeq2/edgeR/EnhancedVolcano/ggplot2
    └── scripts/
        ├── process_fc.py      # cleans featureCounts column names
        └── deseq2.R           # subcommands: pca | contrast
```

---

## 1. Inputs

### Raw FASTQs

Put them under `fastq_dir` (default `data/fastq/`) using the naming convention:

| Layout | Files |
|--------|-------|
| PE     | `{sample}_R1_001.fastq.gz`, `{sample}_R2_001.fastq.gz` |
| SE     | `{sample}_R1_001.fastq.gz` |

If your suffixes differ, override `fastq_suffix_r1` / `fastq_suffix_r2` in
`config/config.yaml`.

### `config/samples.tsv`

Tab-separated, header required:

| sample                | rep | condition | time |
|-----------------------|-----|-----------|------|
| 21161R-74-01_S33_L005 | 1   | DMSO      | 6    |
| 21161R-74-02_S34_L005 | 2   | DMSO      | 6    |
| 21161R-74-03_S35_L005 | 1   | ISD       | 6    |
| ...                   | ... | ...       | NA   |

- `sample` must match the FASTQ basename (without the suffix).
- `time` may be `NA` if not applicable.

### `config/contrasts.tsv`

Defines pairwise DESeq2 comparisons:

| name             | group_col | numerator | denominator | subset_col | subset_val |
|------------------|-----------|-----------|-------------|------------|------------|
| Drug_vs_DMSO_6h   | condition | Drug       | DMSO        | time       | 6          |
| Drug_vs_DMSO_28h  | condition | Drug       | DMSO        | time       | 28         |

Leave `subset_col` / `subset_val` blank to use all samples.

---

## 2. Configure

Edit `config/config.yaml`:

| Key | Description |
|-----|-------------|
| `library_type` | `"PE"` or `"SE"` |
| `fastq_dir` | directory containing raw FASTQs |
| `fastq_suffix_r1` / `fastq_suffix_r2` | suffix patterns |
| `hisat2_index` | HISAT2 index prefix (path passed to `hisat2 -x`) |
| `gtf` | GTF annotation for featureCounts |
| `trim_galore.*`, `hisat2.*`, `samtools.*`, `featurecounts.*` | tool params |
| `deseq2.padj_cutoff`, `deseq2.log2fc_cutoff` | significance thresholds |
| `threads.*` | per-rule thread budgets |

---

## 3. Run

### Local

```bash
mamba create -n snakemake -c bioconda -c conda-forge snakemake=8
conda activate snakemake

snakemake -n                        # dry run: show what will execute
snakemake --use-conda --cores 32    # full run
snakemake --use-conda --cores 32 --rerun-incomplete
```

### SLURM (via profile)

`profiles/slurm/config.yaml` example:

```yaml
executor: slurm
jobs: 50
use-conda: true
software-deployment-method:
  - conda
conda-frontend: mamba

default-resources:
  - slurm_account=your_allocation_for_SLURM
  - slurm_partition=short
  - runtime=240
  - mem_mb=50000
  - cpus_per_task=16

set-threads:
  hisat2_pe: 16
  hisat2_se: 16
  featurecounts: 8
  trim_galore_pe: 8
  trim_galore_se: 8
```

Then:

```bash
snakemake --profile profiles/slurm
```

> With a SLURM profile, set `--jobs N` (max concurrent jobs) — `--cores all`
> only matters for local execution.

---

## 4. Outputs

```
results/
├── fastqc/raw/                              # per-sample FastQC HTML/zip
├── trimmed/                                 # TrimGalore output + .trim.log
├── bam/
│   ├── {sample}.bam, .bam.bai
│   └── {sample}.flagstat
├── counts/
│   ├── raw_counts.txt                       # featureCounts output
│   ├── raw_counts.txt.summary               # featureCounts QC
│   └── counts_processed.tsv                 # cleaned column names
├── deseq2/
│   ├── PCA.pdf
│   ├── dds_all.rds
│   └── {contrast}/
│       ├── results.tsv                      # full DE table
│       ├── significant.tsv                  # padj/log2FC filtered
│       └── volcano.pdf                      # labelled volcano
└── multiqc/multiqc_report.html              # aggregated QC
logs/                                        # per-rule stderr/stdout
```

---

## 5. How rule order is decided

Snakemake does **not** run rules in the order they appear in the file. It
builds a DAG by working backward from `rule all`'s targets:

1. Match each target to the rule whose `output:` produces it.
2. Treat that rule's `input:` as new sub-targets and recurse.
3. Stop when sub-targets exist on disk (raw FASTQs).
4. Execute leaves first; run independent branches in parallel.

Useful introspection commands:

```bash
snakemake --rulegraph | dot -Tpng > rules.png
snakemake --dag       | dot -Tpng > dag.png
snakemake -n -p                 # dry-run with shell commands
snakemake --reason -n           # explain why each job will run
snakemake --summary             # status of every output
```

---

## 6. Troubleshooting

### `Nothing to be done`

Snakemake thinks all `rule all` targets exist. Check:

```bash
ls -R results/ 2>/dev/null | head
snakemake --summary | head -30
```

Force a rerun:

```bash
snakemake --use-conda -F                       # rebuild everything
snakemake --use-conda -R deseq2_contrast       # rerun a specific rule
rm -rf results/ logs/ .snakemake/              # nuclear option
```

If `SAMPLES` is empty (parsing of `samples.tsv` failed), `expand()` produces
nothing and only the static targets remain. Verify:

```bash
python -c "import pandas as pd; \
  print(pd.read_csv('config/samples.tsv', sep='\t', dtype=str))"
```

The header must be exactly `sample<TAB>rep<TAB>condition<TAB>time` — comma
separators or BOMs silently produce 0 rows.

### `ImportError: Unable to import required dependencies: numpy ...`

You're running the system Python (e.g., 3.13 in `~/.local`) instead of the
conda env's Python 3.11. Two fixes:

1. **Make sure `--use-conda` is active.** In your profile:

   ```yaml
   use-conda: true
   software-deployment-method:
     - conda
   ```

   Then verify the env was built: `ls .snakemake/conda/`.

2. **Block `~/.local` from shadowing the env.** The Python rule now sets
   `PYTHONNOUSERSITE=1`, but if you hit a similar issue in R, add to the
   shell command: `R_LIBS_USER=/dev/null Rscript ...`.

### `Permission denied` / `command not found` on cluster

Add `--use-conda` to the profile (above). Without it, the SLURM worker uses
whatever Python/R lives on the node — often broken or ancient.

---

## 7. Adding new comparisons

Append a row to `config/contrasts.tsv` and re-run. The wildcard
`{contrast}` in `rule deseq2_contrast` will pick it up automatically and
produce `results/deseq2/<name>/{results.tsv,significant.tsv,volcano.pdf}`.

## 8. Switching between gene- and transcript-level counts

In `config/config.yaml`:

```yaml
featurecounts:
  feature_type: "exon"
  attribute: "gene_id"        # or "transcript_id"
```

Then `snakemake -F featurecounts process_counts` to refresh.

---

## Citing the tools

| Step | Tool | Reference |
|------|------|-----------|
| QC   | FastQC | Andrews S. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| Trim | TrimGalore / cutadapt | Krueger F.; Martin M., EMBnet J. 2011 |
| Align | HISAT2 | Kim et al., Nat. Biotechnol. 2019 |
| BAM ops | samtools | Danecek et al., GigaScience 2021 |
| Counts | featureCounts (subread) | Liao et al., Bioinformatics 2014 |
| DE    | DESeq2 | Love et al., Genome Biol. 2014 |
| Volcano | EnhancedVolcano | Blighe et al. |
| QC report | MultiQC | Ewels et al., Bioinformatics 2016 |