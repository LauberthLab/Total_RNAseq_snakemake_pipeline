# HISAT2 alignment + samtools sort/index.

def hisat2_strand_flag():
    s = config["hisat2"]["strandness"]
    return f"--rna-strandness {s}" if s else ""


if PAIRED:

    rule hisat2_pe:
        input:
            unpack(hisat2_input),
        output:
            sam=temp("results/bam/{sample}.sam"),
        log:
            "logs/hisat2/{sample}.log",
        threads: config["threads"]["align"]
        params:
            idx=config["hisat2_index"],
            strand=hisat2_strand_flag(),
            extra=config["hisat2"]["extra"],
        conda:
             "../envs/tools.yaml"
        shell:
            r"""
            hisat2 -p {threads} -x {params.idx} {params.strand} {params.extra} \
                -1 {input.r1} -2 {input.r2} -S {output.sam} 2> {log}
            """

else:

    rule hisat2_se:
        input:
            unpack(hisat2_input),
        output:
            sam=temp("results/bam/{sample}.sam"),
        log:
            "logs/hisat2/{sample}.log",
        threads: config["threads"]["align"]
        params:
            idx=config["hisat2_index"],
            strand=hisat2_strand_flag(),
            extra=config["hisat2"]["extra"],
        conda:
            "../envs/tools.yaml"
        shell:
            r"""
            hisat2 -p {threads} -x {params.idx} {params.strand} {params.extra} \
                -U {input.r1} -S {output.sam} 2> {log}
            """


rule sam_to_sorted_bam:
    input:
        sam="results/bam/{sample}.sam",
    output:
        bam="results/bam/{sample}.bam",
    log:
        "logs/samtools/{sample}.sort.log",
    threads: config["threads"]["sort"]
    params:
        mapq=config["samtools"]["min_mapq"],
    conda:
        "../envs/tools.yaml"
    shell:
        r"""
        samtools view -@ {threads} -q {params.mapq} -b {input.sam} \
            | samtools sort -@ {threads} -o {output.bam} - 2> {log}
        """


rule index_bam:
    input:
        bam="results/bam/{sample}.bam",
    output:
        bai="results/bam/{sample}.bam.bai",
    threads: config["threads"]["sort"]
    conda:
        "../envs/tools.yaml"
    shell:
        "samtools index -@ {threads} {input.bam}"


rule flagstat:
    input:
        bam="results/bam/{sample}.bam",
    output:
        "results/bam/{sample}.flagstat",
    conda:
        "../envs/tools.yaml"
    shell:
        "samtools flagstat {input.bam} > {output}"