# DESeq2 — global PCA + per-contrast DE results and volcano plots.

rule deseq2_pca:
    input:
        counts="results/counts/counts_processed.tsv",
        meta=config["samples"],
    output:
        pdf="results/deseq2/PCA.pdf",
        rds="results/deseq2/dds_all.rds",
    log:
        "logs/deseq2/pca.log",
    params:
        min_cpm=config["deseq2"]["min_cpm"],
        min_samples=config["deseq2"]["min_samples_above_cpm"],
    conda:
         "../envs/r.yaml"
    shell:
        r"""
        Rscript workflow/scripts/deseq2.R pca \
            --counts {input.counts} \
            --meta {input.meta} \
            --min-cpm {params.min_cpm} \
            --min-samples {params.min_samples} \
            --out-pdf {output.pdf} \
            --out-rds {output.rds} > {log} 2>&1
        """


rule deseq2_contrast:
    input:
        counts="results/counts/counts_processed.tsv",
        meta=config["samples"],
        contrasts=config["contrasts"],
    output:
        results="results/deseq2/{contrast}/results.tsv",
        sig="results/deseq2/{contrast}/significant.tsv",
        volcano="results/deseq2/{contrast}/volcano.pdf",
    log:
        "logs/deseq2/{contrast}.log",
    params:
        contrast="{contrast}",
        min_cpm=config["deseq2"]["min_cpm"],
        min_samples=config["deseq2"]["min_samples_above_cpm"],
        padj=config["deseq2"]["padj_cutoff"],
        l2fc=config["deseq2"]["log2fc_cutoff"],
        top=config["deseq2"]["top_labels"],
    conda:
        "../envs/r.yaml"
    shell:
        r"""
        Rscript workflow/scripts/deseq2.R contrast \
            --counts {input.counts} \
            --meta {input.meta} \
            --contrasts {input.contrasts} \
            --name {params.contrast} \
            --min-cpm {params.min_cpm} \
            --min-samples {params.min_samples} \
            --padj {params.padj} \
            --log2fc {params.l2fc} \
            --top-labels {params.top} \
            --out-results {output.results} \
            --out-sig {output.sig} \
            --out-volcano {output.volcano} > {log} 2>&1
        """