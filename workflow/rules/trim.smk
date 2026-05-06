# TrimGalore — paired-end or single-end depending on config["library_type"].

if PAIRED:

    rule trim_galore_pe:
        input:
            unpack(trim_inputs),
        output:
            r1="results/trimmed/{sample}_R1_001_val_1.fq.gz",
            r2="results/trimmed/{sample}_R2_001_val_2.fq.gz",
            log="results/trimmed/{sample}.trim.log",
        log:
            "logs/trim_galore/{sample}.log",
        threads: config["threads"]["trim"]
        params:
            quality=config["trim_galore"]["quality"],
            adapter=config["trim_galore"]["adapter"],
            extra=config["trim_galore"]["extra"],
            outdir="results/trimmed",
        conda:
             "../envs/tools.yaml"
        shell:
            r"""
            trim_galore -j {threads} {params.adapter} --paired \
                -q {params.quality} {params.extra} \
                -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1
            cp {log} {output.log}
            """

else:

    rule trim_galore_se:
        input:
            unpack(trim_inputs),
        output:
            r1="results/trimmed/{sample}_R1_001_trimmed.fq.gz",
            log="results/trimmed/{sample}.trim.log",
        log:
            "logs/trim_galore/{sample}.log",
        threads: config["threads"]["trim"]
        params:
            quality=config["trim_galore"]["quality"],
            adapter=config["trim_galore"]["adapter"],
            extra=config["trim_galore"]["extra"],
            outdir="results/trimmed",
        conda:
            "../envs/tools.yaml"
        shell:
            r"""
            trim_galore -j {threads} {params.adapter} \
                -q {params.quality} {params.extra} \
                -o {params.outdir} {input.r1} > {log} 2>&1
            cp {log} {output.log}
            """