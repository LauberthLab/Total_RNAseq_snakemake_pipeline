# featureCounts on all BAMs at once, then tidy column names.

rule featurecounts:
    input:
        bams=expand("results/bam/{sample}.bam", sample=SAMPLES),
        bais=expand("results/bam/{sample}.bam.bai", sample=SAMPLES),
        gtf=config["gtf"],
    output:
        counts="results/counts/raw_counts.txt",
        summary="results/counts/raw_counts.txt.summary",
    log:
        "logs/featurecounts.log",
    threads: config["threads"]["featurecounts"]
    params:
        feature=config["featurecounts"]["feature_type"],
        attr=config["featurecounts"]["attribute"],
        strand=config["featurecounts"]["strandness"],
        paired_flag="-p --countReadPairs" if PAIRED else "",
        extra=config["featurecounts"]["extra"],
    conda:
         "../envs/tools.yaml"
    shell:
        r"""
        featureCounts -T {threads} -t {params.feature} -F GTF \
            -g {params.attr} -s {params.strand} {params.paired_flag} \
            {params.extra} -a {input.gtf} \
            -o {output.counts} {input.bams} > {log} 2>&1
        """


rule process_counts:
    input:
        "results/counts/raw_counts.txt",
    output:
        "results/counts/counts_processed.tsv",
    log:
        "logs/process_counts.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/process_fc.py {input} {output} > {log} 2>&1"