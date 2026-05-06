#!/usr/bin/env Rscript
# =====================================================================
# DESeq2 analysis driver
# =====================================================================
# Two subcommands:
#   pca       — global PCA across all samples
#   contrast  — pairwise DE + volcano for one contrast row of contrasts.tsv
# =====================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
  library(EnhancedVolcano)
})

# ---------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------
read_counts <- function(path) {
  counts <- read.delim(path, check.names = FALSE)
  colnames(counts) <- gsub("^X", "", colnames(counts))
  counts <- counts[!duplicated(counts$Geneid), ]
  rownames(counts) <- counts$Geneid
  drop <- intersect(c("Geneid", "Length"), colnames(counts))
  counts[, !(colnames(counts) %in% drop), drop = FALSE]
}

read_meta <- function(path) {
  meta <- read.delim(path, stringsAsFactors = FALSE)
  rownames(meta) <- meta$sample
  meta
}

filter_low <- function(counts, min_cpm, min_samples) {
  cpm_data <- edgeR::cpm(counts, normalized.lib.sizes = TRUE)
  keep <- rowSums(cpm_data > min_cpm) >= min_samples
  message(sprintf("Kept %d / %d features (CPM > %.2f in >= %d samples)",
                  sum(keep), nrow(counts), min_cpm, min_samples))
  counts[keep, , drop = FALSE]
}

align_counts_meta <- function(counts, meta) {
  shared <- intersect(colnames(counts), rownames(meta))
  if (length(shared) == 0) {
    stop("No overlap between count column names and metadata sample ids. ",
         "Counts cols: ", paste(head(colnames(counts)), collapse = ","),
         " | Meta rows: ", paste(head(rownames(meta)), collapse = ","))
  }
  list(counts = counts[, shared, drop = FALSE],
       meta   = meta[shared, , drop = FALSE])
}

# ---------------------------------------------------------------------
# PCA subcommand
# ---------------------------------------------------------------------
run_pca <- function(args) {
  counts <- read_counts(args$counts)
  meta   <- read_meta(args$meta)

  aligned <- align_counts_meta(counts, meta)
  counts <- aligned$counts
  meta   <- aligned$meta

  counts <- filter_low(counts, args$`min-cpm`, args$`min-samples`)

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~1)
  rl  <- rlog(dds, blind = TRUE)

  intgroup <- intersect(c("condition", "time"), colnames(meta))
  if (length(intgroup) == 0) intgroup <- colnames(meta)[1]

  pca_df <- plotPCA(rl, intgroup = intgroup, returnData = TRUE)
  pct    <- round(100 * attr(pca_df, "percentVar"))

  aes_args <- list(x = quote(PC1), y = quote(PC2))
  if ("condition" %in% intgroup) aes_args$colour <- quote(condition)
  if ("time" %in% intgroup) {
    pca_df$time <- as.factor(pca_df$time)
    aes_args$shape <- quote(time)
  }

  p <- ggplot(pca_df, do.call(aes, aes_args)) +
    geom_point(size = 4) +
    xlab(paste0("PC1 (", pct[1], "%)")) +
    ylab(paste0("PC2 (", pct[2], "%)")) +
    theme_minimal()

  dir.create(dirname(args$`out-pdf`), showWarnings = FALSE, recursive = TRUE)
  ggsave(args$`out-pdf`, p, width = 7, height = 5)
  saveRDS(dds, args$`out-rds`)
  message("Wrote ", args$`out-pdf`)
}

# ---------------------------------------------------------------------
# Contrast subcommand
# ---------------------------------------------------------------------
run_contrast <- function(args) {
  counts    <- read_counts(args$counts)
  meta      <- read_meta(args$meta)
  contrasts <- read.delim(args$contrasts, stringsAsFactors = FALSE)

  row <- contrasts[contrasts$name == args$name, , drop = FALSE]
  if (nrow(row) != 1) stop("Contrast not found in contrasts file: ", args$name)

  group_col <- row$group_col
  num       <- row$numerator
  den       <- row$denominator
  subset_col <- if ("subset_col" %in% colnames(row)) row$subset_col else ""
  subset_val <- if ("subset_val" %in% colnames(row)) row$subset_val else ""

  if (!group_col %in% colnames(meta)) {
    stop("group_col '", group_col, "' not in metadata columns: ",
         paste(colnames(meta), collapse = ","))
  }

  if (nzchar(subset_col) && nzchar(subset_val)) {
    if (!subset_col %in% colnames(meta)) {
      stop("subset_col '", subset_col, "' not in metadata columns")
    }
    keep_samples <- rownames(meta)[as.character(meta[[subset_col]]) == subset_val]
    meta <- meta[keep_samples, , drop = FALSE]
    message(sprintf("Subset to %s == %s : %d samples",
                    subset_col, subset_val, nrow(meta)))
  }

  meta <- meta[meta[[group_col]] %in% c(num, den), , drop = FALSE]
  if (nrow(meta) < 2) stop("Fewer than 2 samples remain after subsetting.")

  aligned <- align_counts_meta(counts, meta)
  counts  <- aligned$counts
  meta    <- aligned$meta

  counts <- filter_low(counts, args$`min-cpm`, args$`min-samples`)

  meta[[group_col]] <- factor(meta[[group_col]], levels = c(den, num))
  design_formula <- as.formula(paste("~", group_col))

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData   = meta,
                                design    = design_formula)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c(group_col, num, den))
  res <- res[order(res$padj), ]

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

  dir.create(dirname(args$`out-results`), showWarnings = FALSE, recursive = TRUE)
  write.table(res_df, args$`out-results`, sep = "\t",
              quote = FALSE, row.names = FALSE)

  sig <- res_df[!is.na(res_df$padj) &
                  res_df$padj <= args$padj &
                  abs(res_df$log2FoldChange) >= args$log2fc, ]
  write.table(sig, args$`out-sig`, sep = "\t",
              quote = FALSE, row.names = FALSE)
  message(sprintf("Significant: %d (padj<=%.3g, |log2FC|>=%.3f)",
                  nrow(sig), args$padj, args$log2fc))

  # ----- Volcano -----
  res_plot <- res_df[!is.na(res_df$padj), ]
  up   <- sig[sig$log2FoldChange > 0, ]
  down <- sig[sig$log2FoldChange < 0, ]
  top_up   <- head(up$gene[order(up$padj)],   args$`top-labels`)
  top_down <- head(down$gene[order(down$padj)], args$`top-labels`)

  keyvals <- ifelse(
    res_plot$log2FoldChange >  args$log2fc & res_plot$padj < args$padj, "red",
    ifelse(
      res_plot$log2FoldChange < -args$log2fc & res_plot$padj < args$padj, "blue",
      "gray"))
  names(keyvals)[keyvals == "red"]  <- "Upregulated"
  names(keyvals)[keyvals == "blue"] <- "Downregulated"
  names(keyvals)[keyvals == "gray"] <- "Not significant"

  xmax <- max(abs(res_plot$log2FoldChange), na.rm = TRUE)

  v <- EnhancedVolcano(
    res_plot,
    title           = paste0(num, " vs ", den),
    subtitle        = args$name,
    lab             = res_plot$gene,
    selectLab       = c(top_up, top_down),
    x               = "log2FoldChange",
    y               = "padj",
    pCutoffCol      = "padj",
    pCutoff         = args$padj,
    FCcutoff        = args$log2fc,
    colCustom       = keyvals,
    drawConnectors  = TRUE,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendLabSize   = 12,
    labSize         = 4,
    xlim            = c(-xmax, xmax)
  )

  ggsave(args$`out-volcano`, v, width = 8, height = 7)
  message("Wrote ", args$`out-volcano`)
}

# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
parse_args_for <- function(sub) {
  opts <- list(
    make_option("--counts", type = "character"),
    make_option("--meta", type = "character"),
    make_option("--min-cpm", type = "double", default = 0.3),
    make_option("--min-samples", type = "integer", default = 2)
  )
  if (sub == "pca") {
    opts <- c(opts, list(
      make_option("--out-pdf", type = "character"),
      make_option("--out-rds", type = "character")
    ))
  } else {
    opts <- c(opts, list(
      make_option("--contrasts", type = "character"),
      make_option("--name", type = "character"),
      make_option("--padj", type = "double", default = 0.05),
      make_option("--log2fc", type = "double", default = 0.5849),
      make_option("--top-labels", type = "integer", default = 10),
      make_option("--out-results", type = "character"),
      make_option("--out-sig", type = "character"),
      make_option("--out-volcano", type = "character")
    ))
  }
  parse_args(OptionParser(option_list = opts),
             args = commandArgs(trailingOnly = TRUE)[-1])
}

main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  if (length(argv) < 1) {
    stop("Usage: deseq2.R {pca|contrast} [options]")
  }
  sub <- argv[1]
  args <- parse_args_for(sub)
  if (sub == "pca")          run_pca(args)
  else if (sub == "contrast") run_contrast(args)
  else stop("Unknown subcommand: ", sub)
}

main()