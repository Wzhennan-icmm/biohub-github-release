#!/usr/bin/env Rscript

parse_args <- function(values) {
  if (length(values) %% 2 != 0) stop("arguments must be --name value pairs", call. = FALSE)
  result <- list()
  for (index in seq(1, length(values), by = 2)) {
    key <- sub("^--", "", values[[index]])
    if (key == values[[index]] || key %in% names(result)) {
      stop(sprintf("invalid or repeated argument: %s", values[[index]]), call. = FALSE)
    }
    result[[key]] <- values[[index + 1]]
  }
  result
}

require_args <- function(args, names) {
  missing <- names[!names %in% names(args)]
  if (length(missing) > 0) stop(sprintf("missing arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

write_tsv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c(
  "counts", "samples", "contrasts", "design", "minimum-count", "minimum-samples",
  "alpha", "minimum-absolute-log2-fold-change", "independent-filtering", "differential",
  "contrast-summary", "normalized-counts", "sample-qc", "pca", "log"
)
require_args(args, required)

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("R package DESeq2 is required", call. = FALSE)
}

minimum_count <- as.integer(args[["minimum-count"]])
minimum_samples <- as.integer(args[["minimum-samples"]])
alpha <- as.numeric(args[["alpha"]])
lfc_threshold <- as.numeric(args[["minimum-absolute-log2-fold-change"]])
independent_filtering <- tolower(args[["independent-filtering"]]) == "true"
if (any(is.na(c(minimum_count, minimum_samples, alpha, lfc_threshold))) || minimum_count < 0 || minimum_samples < 1 || alpha <= 0 || alpha > 1 || lfc_threshold < 0) stop("invalid numeric parameter", call. = FALSE)

counts_frame <- read.delim(args$counts, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts_frame) < 3 || names(counts_frame)[[1]] != "gene_id") {
  stop("counts matrix needs gene_id and at least two sample columns", call. = FALSE)
}
if (anyDuplicated(counts_frame$gene_id) || any(counts_frame$gene_id == "")) stop("gene_id values must be non-empty and unique", call. = FALSE)
count_values <- as.matrix(counts_frame[, -1, drop = FALSE])
suppressWarnings(storage.mode(count_values) <- "numeric")
if (anyNA(count_values) || any(!is.finite(count_values)) || any(count_values < 0) || any(abs(count_values - round(count_values)) > 1e-8)) {
  stop("count values must be finite non-negative integers", call. = FALSE)
}
storage.mode(count_values) <- "integer"
rownames(count_values) <- counts_frame$gene_id

samples <- read.delim(args$samples, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample_id" %in% names(samples) || anyNA(samples$sample_id) || anyDuplicated(samples$sample_id) || any(samples$sample_id == "")) {
  stop("samples_tsv needs unique non-empty sample_id values", call. = FALSE)
}
if (!setequal(colnames(count_values), samples$sample_id)) {
  missing_metadata <- setdiff(colnames(count_values), samples$sample_id)
  missing_counts <- setdiff(samples$sample_id, colnames(count_values))
  stop(sprintf("sample mismatch; missing metadata=[%s], missing counts=[%s]", paste(missing_metadata, collapse = ","), paste(missing_counts, collapse = ",")), call. = FALSE)
}
samples <- samples[match(colnames(count_values), samples$sample_id), , drop = FALSE]
rownames(samples) <- samples$sample_id

design_formula <- as.formula(args$design)
design_variables <- all.vars(design_formula)
if (length(design_variables) == 0 || any(!design_variables %in% names(samples))) {
  stop("design must name columns present in samples_tsv", call. = FALSE)
}
if (any(!complete.cases(samples[, design_variables, drop = FALSE]))) stop("design variables contain missing values", call. = FALSE)
for (name in design_variables) {
  if (is.character(samples[[name]]) || is.logical(samples[[name]])) samples[[name]] <- factor(samples[[name]])
}

contrasts <- read.delim(args$contrasts, check.names = FALSE, stringsAsFactors = FALSE)
contrast_columns <- c("contrast_id", "factor", "numerator", "denominator")
if (nrow(contrasts) == 0 || any(!contrast_columns %in% names(contrasts))) stop("contrasts manifest needs contrast_id, factor, numerator and denominator", call. = FALSE)
if (anyNA(contrasts[, contrast_columns, drop = FALSE]) || anyDuplicated(contrasts$contrast_id) || any(!grepl("^[A-Za-z0-9._-]+$", contrasts$contrast_id))) stop("contrast rows must be complete and contrast_id values must be unique safe identifiers", call. = FALSE)
for (index in seq_len(nrow(contrasts))) {
  factor_name <- contrasts$factor[[index]]
  if (!factor_name %in% design_variables) stop(sprintf("contrast factor is not in design: %s", factor_name), call. = FALSE)
  if (!is.factor(samples[[factor_name]])) stop(sprintf("contrast factor must be categorical: %s", factor_name), call. = FALSE)
  levels_present <- levels(samples[[factor_name]])
  requested <- c(contrasts$numerator[[index]], contrasts$denominator[[index]])
  if (any(!requested %in% levels_present) || requested[[1]] == requested[[2]]) {
    stop(sprintf("invalid levels for contrast %s; available=[%s]", contrasts$contrast_id[[index]], paste(levels_present, collapse = ",")), call. = FALSE)
  }
}

keep <- rowSums(count_values >= minimum_count) >= minimum_samples
if (!any(keep)) stop("no genes remain after count filter", call. = FALSE)
filtered_counts <- count_values[keep, , drop = FALSE]
dds <- DESeq2::DESeqDataSetFromMatrix(countData = filtered_counts, colData = samples, design = design_formula)
design_matrix <- model.matrix(design_formula, data = as.data.frame(SummarizedExperiment::colData(dds)))
if (qr(design_matrix)$rank < ncol(design_matrix)) stop("design matrix is not full rank", call. = FALSE)
dds <- DESeq2::DESeq(dds, quiet = TRUE)

normalized <- as.data.frame(DESeq2::counts(dds, normalized = TRUE), check.names = FALSE)
normalized <- cbind(gene_id = rownames(normalized), normalized)
write_tsv(normalized, args[["normalized-counts"]])

all_results <- list()
summary_rows <- list()
for (index in seq_len(nrow(contrasts))) {
  row <- contrasts[index, ]
  result <- DESeq2::results(
    dds,
    contrast = c(row$factor, row$numerator, row$denominator),
    alpha = alpha,
    independentFiltering = independent_filtering
  )
  result_frame <- as.data.frame(result)
  result_frame$gene_id <- rownames(result_frame)
  result_frame$contrast_id <- row$contrast_id
  result_frame$factor <- row$factor
  result_frame$numerator <- row$numerator
  result_frame$denominator <- row$denominator
  result_frame$significant <- !is.na(result_frame$padj) & result_frame$padj <= alpha & abs(result_frame$log2FoldChange) >= lfc_threshold
  result_frame <- result_frame[, c("contrast_id", "factor", "numerator", "denominator", "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significant")]
  all_results[[index]] <- result_frame
  summary_rows[[index]] <- data.frame(
    contrast_id = row$contrast_id,
    factor = row$factor,
    numerator = row$numerator,
    denominator = row$denominator,
    tested_genes = sum(!is.na(result_frame$pvalue)),
    adjusted_genes = sum(!is.na(result_frame$padj)),
    significant_genes = sum(result_frame$significant),
    up_genes = sum(result_frame$significant & result_frame$log2FoldChange > 0),
    down_genes = sum(result_frame$significant & result_frame$log2FoldChange < 0),
    stringsAsFactors = FALSE
  )
}
write_tsv(do.call(rbind, all_results), args$differential)
write_tsv(do.call(rbind, summary_rows), args[["contrast-summary"]])

raw_library <- colSums(count_values)
detected <- colSums(count_values > 0)
log_norm <- log2(t(as.matrix(normalized[, -1, drop = FALSE])) + 1)
pca_scores <- matrix(NA_real_, nrow = nrow(samples), ncol = 2, dimnames = list(rownames(samples), c("PC1", "PC2")))
if (nrow(log_norm) >= 2 && ncol(log_norm) >= 2) {
  pca <- prcomp(log_norm, center = TRUE, scale. = FALSE)
  dimensions <- min(2, ncol(pca$x))
  pca_scores[, seq_len(dimensions)] <- pca$x[, seq_len(dimensions), drop = FALSE]
}
sample_qc <- data.frame(sample_id = rownames(samples), raw_library_size = raw_library[rownames(samples)], detected_genes = detected[rownames(samples)], PC1 = pca_scores[, 1], PC2 = pca_scores[, 2], stringsAsFactors = FALSE)
write_tsv(sample_qc, args[["sample-qc"]])

dir.create(dirname(args$pca), recursive = TRUE, showWarnings = FALSE)
pdf(args$pca, width = 7, height = 6)
if (all(is.finite(pca_scores[, 1:2]))) {
  plot(pca_scores[, 1], pca_scores[, 2], xlab = "PC1", ylab = "PC2", main = "RNA-seq sample PCA", pch = 19)
  text(pca_scores[, 1], pca_scores[, 2], labels = rownames(samples), pos = 3, cex = 0.7)
} else {
  plot.new(); text(0.5, 0.5, "PCA unavailable for current matrix")
}
dev.off()

log_lines <- c(
  sprintf("DESeq2_version\t%s", as.character(utils::packageVersion("DESeq2"))),
  sprintf("input_genes\t%d", nrow(count_values)),
  sprintf("retained_genes\t%d", nrow(filtered_counts)),
  sprintf("samples\t%d", ncol(filtered_counts)),
  sprintf("design\t%s", args$design),
  sprintf("contrasts\t%d", nrow(contrasts)),
  sprintf("alpha\t%.17g", alpha),
  sprintf("minimum_absolute_log2_fold_change\t%.17g", lfc_threshold),
  sprintf("independent_filtering\t%s", independent_filtering)
)
dir.create(dirname(args$log), recursive = TRUE, showWarnings = FALSE)
writeLines(log_lines, args$log)
