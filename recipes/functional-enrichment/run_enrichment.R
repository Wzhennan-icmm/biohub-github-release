#!/usr/bin/env Rscript

parse_args <- function(values) {
  result <- list(); index <- 1
  while (index <= length(values)) {
    key <- values[[index]]
    if (!startsWith(key, "--") || index == length(values)) stop(sprintf("invalid argument: %s", key), call. = FALSE)
    key <- sub("^--", "", key)
    if (key %in% names(result)) stop(sprintf("repeated argument: --%s", key), call. = FALSE)
    result[[key]] <- values[[index + 1]]; index <- index + 2
  }
  result
}

write_tsv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("foreground", "background", "associations", "sources", "minimum-term-size", "maximum-term-size", "minimum-overlap", "fdr", "adjustment-scope", "require-foreground-in-background", "plot-top-terms", "all", "significant", "summary", "plot", "log")
missing <- required[!required %in% names(args)]
if (length(missing)) stop(sprintf("missing arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)

minimum_term_size <- as.integer(args[["minimum-term-size"]])
maximum_term_size <- as.integer(args[["maximum-term-size"]])
minimum_overlap <- as.integer(args[["minimum-overlap"]])
fdr <- as.numeric(args$fdr)
plot_top <- as.integer(args[["plot-top-terms"]])
scope <- args[["adjustment-scope"]]
sources <- unique(strsplit(args$sources, ",", fixed = TRUE)[[1]])
require_in_background <- tolower(args[["require-foreground-in-background"]]) == "true"
if (anyNA(c(minimum_term_size, maximum_term_size, minimum_overlap, fdr, plot_top)) || minimum_term_size < 1 || maximum_term_size < 1 || minimum_overlap < 1 || fdr <= 0 || fdr > 1 || plot_top < 1) stop("invalid numeric parameter", call. = FALSE)
sources <- trimws(sources)
if (minimum_term_size > maximum_term_size) stop("minimum_term_size exceeds maximum_term_size", call. = FALSE)
if (!length(sources) || any(sources == "")) stop("sources must be non-empty", call. = FALSE)

foreground <- read.delim(args$foreground, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("set_id", "gene_id") %in% names(foreground)) || nrow(foreground) == 0) stop("foreground_tsv needs set_id and gene_id", call. = FALSE)
if (anyNA(foreground[, c("set_id", "gene_id"), drop = FALSE]) || any(!grepl("^[A-Za-z0-9._-]+$", foreground$set_id)) || any(foreground$gene_id == "")) stop("invalid set_id or empty foreground gene_id", call. = FALSE)
foreground <- unique(foreground[, c("set_id", "gene_id")])

background_table <- read.delim(args$background, header = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
background <- unique(trimws(background_table[[1]]))
background <- background[background != ""]
if (!length(background)) stop("background gene list is empty", call. = FALSE)

outside <- setdiff(foreground$gene_id, background)
if (require_in_background && length(outside)) stop(sprintf("foreground genes outside background: %s", paste(head(sort(outside), 20), collapse = ",")), call. = FALSE)

associations <- read.delim(args$associations, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c("gene_id", "term_id", "source")
if (!all(needed %in% names(associations)) || nrow(associations) == 0) stop("associations_tsv needs gene_id, term_id and source", call. = FALSE)
if (anyNA(associations[, needed, drop = FALSE]) || any(associations$gene_id == "") || any(associations$term_id == "") || any(associations$source == "")) stop("association keys must be complete and non-empty", call. = FALSE)
if (!"term_name" %in% names(associations)) associations$term_name <- associations$term_id
associations$term_name[is.na(associations$term_name) | associations$term_name == ""] <- associations$term_id[is.na(associations$term_name) | associations$term_name == ""]
associations <- unique(associations[, c("gene_id", "term_id", "source", "term_name")])
associations <- associations[associations$gene_id %in% background & associations$source %in% sources, , drop = FALSE]
if (!nrow(associations)) stop("no selected associations overlap background", call. = FALSE)
missing_sources <- setdiff(sources, unique(associations$source))
if (length(missing_sources)) stop(sprintf("selected sources have no background associations: %s", paste(missing_sources, collapse = ",")), call. = FALSE)

name_counts <- aggregate(term_name ~ source + term_id, associations, function(x) length(unique(x)))
if (any(name_counts$term_name > 1)) stop("a term_id/source pair has conflicting term_name values", call. = FALSE)

sets <- sort(unique(foreground$set_id))
rows <- list(); row_index <- 1
coverage_rows <- list()
N <- length(background)
for (set_id in sets) {
  original <- unique(foreground$gene_id[foreground$set_id == set_id])
  selected <- intersect(original, background)
  annotated_selected <- intersect(selected, unique(associations$gene_id))
  coverage_rows[[length(coverage_rows) + 1]] <- data.frame(
    set_id = set_id,
    input_genes = length(original),
    genes_in_background = length(selected),
    genes_outside_background = length(setdiff(original, background)),
    annotated_genes = length(annotated_selected),
    stringsAsFactors = FALSE
  )
  if (!length(selected)) next
  for (source in sources) {
    source_data <- associations[associations$source == source, , drop = FALSE]
    for (term_id in sort(unique(source_data$term_id))) {
      term_data <- source_data[source_data$term_id == term_id, , drop = FALSE]
      term_genes <- unique(term_data$gene_id)
      K <- length(term_genes); k_genes <- sort(intersect(selected, term_genes)); k <- length(k_genes)
      if (K < minimum_term_size || K > maximum_term_size || k < minimum_overlap) next
      p_value <- phyper(k - 1, K, N - K, length(selected), lower.tail = FALSE)
      rows[[row_index]] <- data.frame(
        set_id = set_id,
        source = source,
        term_id = term_id,
        term_name = term_data$term_name[[1]],
        background_size = N,
        term_size = K,
        foreground_size = length(selected),
        overlap_size = k,
        fold_enrichment = (k / length(selected)) / (K / N),
        pvalue = p_value,
        overlap_genes = paste(k_genes, collapse = ","),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1
    }
  }
}

columns <- c("set_id", "source", "term_id", "term_name", "background_size", "term_size", "foreground_size", "overlap_size", "fold_enrichment", "pvalue", "padj", "significant", "overlap_genes")
if (!length(rows)) {
  all_results <- as.data.frame(setNames(replicate(length(columns), character(0), simplify = FALSE), columns), stringsAsFactors = FALSE)
} else {
  all_results <- do.call(rbind, rows)
  all_results$padj <- NA_real_
  groups <- switch(
    scope,
    set_source = interaction(all_results$set_id, all_results$source, drop = TRUE),
    set = factor(all_results$set_id),
    global = factor(rep("global", nrow(all_results))),
    stop("invalid adjustment_scope", call. = FALSE)
  )
  for (group in levels(groups)) {
    selected <- which(groups == group)
    all_results$padj[selected] <- p.adjust(all_results$pvalue[selected], method = "BH")
  }
  all_results$significant <- all_results$padj <= fdr
  all_results <- all_results[order(all_results$padj, all_results$pvalue, all_results$set_id, all_results$source, all_results$term_id), columns]
}
write_tsv(all_results, args$all)
significant <- all_results[!is.na(all_results$padj) & all_results$padj <= fdr, , drop = FALSE]
write_tsv(significant, args$significant)

coverage <- do.call(rbind, coverage_rows)
coverage$tested_terms <- vapply(coverage$set_id, function(x) sum(all_results$set_id == x), integer(1))
coverage$significant_terms <- vapply(coverage$set_id, function(x) sum(significant$set_id == x), integer(1))
write_tsv(coverage, args$summary)

dir.create(dirname(args$plot), recursive = TRUE, showWarnings = FALSE)
pdf(args$plot, width = 10, height = 7)
for (set_id in sets) {
  selected <- all_results[all_results$set_id == set_id, , drop = FALSE]
  selected <- head(selected[order(selected$padj, selected$pvalue), , drop = FALSE], plot_top)
  if (!nrow(selected)) {
    plot.new(); title(main = sprintf("Enrichment: %s", set_id)); text(0.5, 0.5, "No tested terms")
  } else {
    score <- -log10(pmax(selected$padj, .Machine$double.xmin))
    labels <- paste(selected$source, selected$term_name, sep = ": ")
    par(mar = c(5, max(8, min(24, max(nchar(labels)) * 0.55)), 4, 2))
    barplot(rev(score), names.arg = rev(labels), horiz = TRUE, las = 1, xlab = "-log10(BH-adjusted p)", main = sprintf("Enrichment: %s", set_id))
  }
}
dev.off()

log_lines <- c(
  sprintf("background_genes\t%d", length(background)),
  sprintf("foreground_sets\t%d", length(sets)),
  sprintf("foreground_genes_outside_background\t%d", length(outside)),
  sprintf("selected_sources\t%s", paste(sources, collapse = ",")),
  sprintf("tested_terms\t%d", nrow(all_results)),
  sprintf("significant_terms\t%d", nrow(significant)),
  sprintf("adjustment_scope\t%s", scope),
  sprintf("fdr\t%.17g", fdr)
)
dir.create(dirname(args$log), recursive = TRUE, showWarnings = FALSE)
writeLines(log_lines, args$log)
