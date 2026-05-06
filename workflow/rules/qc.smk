# FastQC on raw reads + MultiQC aggregating everything downstream.

rule fastqc_raw:
    input:
        lambda wc: raw_fastq(wc.sample, wc.read),
    output:
        html="results/fastqc/raw/{sample}_{read}_001_fastqc.html",
        zip="results/fastqc/raw/{sample}_{read}_001_fastqc.zip",
    log:
        "logs/fastqc/{sample}_{read}.log",
    threads: config["threads"]["fastqc"]
    conda:
         "../envs/tools.yaml"
    shell:
        "fastqc -t {threads} -o results/fastqc/raw {input} > {log} 2>&1"


def multiqc_inputs(wc):
    inputs = []
    reads = ["R1", "R2"] if PAIRED else ["R1"]
    for s in SAMPLES:
        for r in reads:
            inputs.append(f"results/fastqc/raw/{s}_{r}_001_fastqc.zip")
        inputs.append(f"results/trimmed/{s}.trim.log")
        inputs.append(f"results/bam/{s}.flagstat")
    inputs.append("results/counts/raw_counts.txt.summary")
    return inputs


rule multiqc:
    input:
        multiqc_inputs,
    output:
        "results/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc -f -o results/multiqc results/ > {log} 2>&1"