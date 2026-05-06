# =====================================================================
# RNA-seq Snakemake workflow
# =====================================================================
# Steps: FastQC -> TrimGalore -> HISAT2 -> samtools sort/index ->
#        featureCounts -> process_fc -> DESeq2 (PCA, DEGs, volcanoes) ->
#        MultiQC
#
# Run:
#   snakemake --use-conda --cores 32
# Dry run:
#   snakemake -n
# =====================================================================

import os
import re
import pandas as pd

configfile: "config/config.yaml"


# load dependencies


# shell.prefix(
#     "set -euo pipefail; "
#     "module purge; "
#     "module load fastqc; "
#     "module load cutadapt/4.2; "
#     "module load TrimGalore; "
#     "module load samtools; "
#     "module load subread; "
#     "module load deeptools; "
#     "module load multiqc; "
#     "module load bedtools; "
#     "module load hisat2; "
# )



# 
# Load metadata
# 
samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("sample", drop=False)
SAMPLES = samples.index.tolist()

contrasts = pd.read_csv(config["contrasts"], sep="\t", dtype=str).fillna("")
CONTRASTS = contrasts["name"].tolist()

LIB = config["library_type"].upper()
assert LIB in ("PE", "SE"), "library_type must be PE or SE"
PAIRED = LIB == "PE"

wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES),
    read="R[12]",
    contrast="|".join(re.escape(c) for c in CONTRASTS),

FASTQ_DIR = config["fastq_dir"]
R1_SFX = config["fastq_suffix_r1"]
R2_SFX = config["fastq_suffix_r2"]


# Helpers
# 
def raw_fastq(sample, read):
    sfx = R1_SFX if read == "R1" else R2_SFX
    return os.path.join(FASTQ_DIR, f"{sample}{sfx}")

def trim_inputs(wc):
    if PAIRED:
        return {"r1": raw_fastq(wc.sample, "R1"),
                "r2": raw_fastq(wc.sample, "R2")}
    return {"r1": raw_fastq(wc.sample, "R1")}

def trimmed_fastq(sample, read):
    # Trim Galore default naming
    if PAIRED:
        if read == "R1":
            return f"results/trimmed/{sample}_R1_001_val_1.fq.gz"
        return f"results/trimmed/{sample}_R2_001_val_2.fq.gz"
    return f"results/trimmed/{sample}_R1_001_trimmed.fq.gz"

def hisat2_input(wc):
    if PAIRED:
        return {"r1": trimmed_fastq(wc.sample, "R1"),
                "r2": trimmed_fastq(wc.sample, "R2")}
    return {"r1": trimmed_fastq(wc.sample, "R1")}

#    
# Targets
#    
rule all:
    input:
        # QC
        expand("results/fastqc/raw/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        # Alignments
        expand("results/bam/{sample}.bam.bai", sample=SAMPLES),
        # Counts
        "results/counts/raw_counts.txt",
        "results/counts/counts_processed.tsv",
        # DESeq2
        expand("results/deseq2/{contrast}/results.tsv", contrast=CONTRASTS),
        expand("results/deseq2/{contrast}/volcano.pdf", contrast=CONTRASTS),
        "results/deseq2/PCA.pdf",
        # MultiQC
        "results/multiqc/multiqc_report.html",

#    
# Include rule modules
#    
include: "workflow/rules/qc.smk"
include: "workflow/rules/trim.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/counts.smk"
include: "workflow/rules/deseq2.smk"