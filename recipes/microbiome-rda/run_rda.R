#!/usr/bin/env Rscript

parse_args <- function(values) {
  if (length(values) %% 2 != 0) stop("arguments must be --name value pairs", call. = FALSE)
  result <- list()
  for (index in seq(1, length(values), by = 2)) {
    key <- sub("^--", "", values[[index]])
    if (key == values[[index]] || key %in% names(result)) stop(sprintf("invalid or repeated argument: %s", values[[index]]), call. = FALSE)
    result[[key]] <- values[[index + 1]]
  }
  result
}

write_tsv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

anova_frame <- function(value, test_name) {
  frame <- as.data.frame(value, check.names = FALSE)
  frame$row <- rownames(frame)
  rownames(frame) <- NULL
  frame$test <- test_name
  frame[, c("test", "row", setdiff(names(frame), c("test", "row"))), drop = FALSE]
}

score_frame <- function(model, display, scaling) {
  value <- vegan::scores(model, display = display, choices = seq_len(min(2, model$CCA$rank)), scaling = scaling)
  if (is.null(value) || !length(value)) {
    return(data.frame(id = character(0), RDA1 = numeric(0), RDA2 = numeric(0), stringsAsFactors = FALSE))
  }
  value <- as.matrix(value)
  if (ncol(value) < 2) value <- cbind(value, RDA2 = NA_real_)
  frame <- data.frame(id = rownames(value), RDA1 = value[, 1], RDA2 = value[, 2], stringsAsFactors = FALSE)
  rownames(frame) <- NULL
  frame
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("features", "metadata", "constraints", "conditions", "transform", "minimum-prevalence", "minimum-total-abundance", "drop-incomplete-samples", "permutations", "random-seed", "scaling", "model-summary", "features-kept", "overall-test", "term-tests", "axis-tests", "site-scores", "feature-scores", "biplot-scores", "plot", "log")
missing <- required[!required %in% names(args)]
if (length(missing)) stop(sprintf("missing arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
if (!requireNamespace("vegan", quietly = TRUE)) stop("R package vegan is required", call. = FALSE)

split_names <- function(value) if (value == "") character(0) else strsplit(value, ",", fixed = TRUE)[[1]]
constraints <- split_names(args$constraints)
conditions <- split_names(args$conditions)
if (!length(constraints)) stop("at least one constraint is required", call. = FALSE)
if (length(intersect(constraints, conditions))) stop("constraint and condition_covariate sets must not overlap", call. = FALSE)
minimum_prevalence <- as.numeric(args[["minimum-prevalence"]])
minimum_total <- as.numeric(args[["minimum-total-abundance"]])
permutations <- as.integer(args$permutations)
seed <- as.integer(args[["random-seed"]])
scaling <- as.integer(args$scaling)
drop_incomplete <- tolower(args[["drop-incomplete-samples"]]) == "true"
if (anyNA(c(minimum_prevalence, minimum_total, permutations, seed, scaling)) || minimum_prevalence < 0 || minimum_prevalence > 1 || minimum_total < 0 || permutations < 1 || seed < 0 || !scaling %in% c(1, 2)) stop("invalid numeric parameter", call. = FALSE)

feature_frame <- read.delim(args$features, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(feature_frame) < 3 || names(feature_frame)[[1]] != "feature_id") stop("feature table needs feature_id and at least two samples", call. = FALSE)
if (anyDuplicated(feature_frame$feature_id) || any(feature_frame$feature_id == "")) stop("feature_id values must be non-empty and unique", call. = FALSE)
feature_values <- as.matrix(feature_frame[, -1, drop = FALSE])
suppressWarnings(storage.mode(feature_values) <- "numeric")
if (anyNA(feature_values) || any(!is.finite(feature_values)) || any(feature_values < 0)) stop("feature abundances must be finite non-negative numbers", call. = FALSE)
rownames(feature_values) <- feature_frame$feature_id

metadata <- read.delim(args$metadata, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample_id" %in% names(metadata) || anyNA(metadata$sample_id) || anyDuplicated(metadata$sample_id) || any(metadata$sample_id == "")) stop("metadata needs unique non-empty sample_id values", call. = FALSE)
if (!setequal(colnames(feature_values), metadata$sample_id)) stop("feature and metadata sample sets must match exactly", call. = FALSE)
metadata <- metadata[match(colnames(feature_values), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id
model_variables <- c(constraints, conditions)
if (any(!model_variables %in% names(metadata))) stop("constraint or condition column missing from metadata", call. = FALSE)
complete <- complete.cases(metadata[, model_variables, drop = FALSE])
if (any(!complete) && !drop_incomplete) stop(sprintf("metadata has incomplete model rows: %s", paste(metadata$sample_id[!complete], collapse = ",")), call. = FALSE)
metadata <- metadata[complete, , drop = FALSE]
feature_values <- feature_values[, metadata$sample_id, drop = FALSE]
if (nrow(metadata) < 3) stop("fewer than three complete samples remain", call. = FALSE)
for (name in model_variables) if (is.character(metadata[[name]]) || is.logical(metadata[[name]])) metadata[[name]] <- factor(metadata[[name]])

prevalence <- rowMeans(feature_values > 0)
total <- rowSums(feature_values)
keep <- prevalence >= minimum_prevalence & total >= minimum_total
if (sum(keep) < 2) stop("fewer than two features remain after filtering", call. = FALSE)
kept_table <- data.frame(feature_id = rownames(feature_values), prevalence = prevalence, total_abundance = total, retained = keep, stringsAsFactors = FALSE)
write_tsv(kept_table, args[["features-kept"]])
feature_values <- feature_values[keep, , drop = FALSE]

community <- t(feature_values)
sample_totals <- rowSums(community)
if (args$transform %in% c("relative", "hellinger") && any(sample_totals <= 0)) stop("relative and Hellinger transforms require positive sample totals", call. = FALSE)
community <- switch(
  args$transform,
  relative = community / sample_totals,
  hellinger = sqrt(community / sample_totals),
  log1p = log1p(community),
  none = community,
  stop("invalid transform", call. = FALSE)
)

rhs <- paste(constraints, collapse = " + ")
if (length(conditions)) rhs <- paste(rhs, sprintf("Condition(%s)", paste(conditions, collapse = " + ")), sep = " + ")
model_formula <- as.formula(paste("community ~", rhs))
environment(model_formula) <- environment()
model <- vegan::rda(model_formula, data = metadata)
if (is.null(model$CCA) || model$CCA$rank < 1) stop("RDA has no constrained axis; inspect constraints and feature matrix", call. = FALSE)

set.seed(seed)
overall_test <- anova_frame(vegan::anova.cca(model, permutations = permutations), "overall")
set.seed(seed)
term_tests <- anova_frame(vegan::anova.cca(model, by = "term", permutations = permutations), "term")
set.seed(seed)
axis_tests <- anova_frame(vegan::anova.cca(model, by = "axis", permutations = permutations), "axis")
write_tsv(overall_test, args[["overall-test"]])
write_tsv(term_tests, args[["term-tests"]])
write_tsv(axis_tests, args[["axis-tests"]])

sites <- score_frame(model, "sites", scaling)
species <- score_frame(model, "species", scaling)
biplot <- score_frame(model, "bp", scaling)
write_tsv(sites, args[["site-scores"]])
write_tsv(species, args[["feature-scores"]])
write_tsv(biplot, args[["biplot-scores"]])

constrained_eigen <- model$CCA$eig
total_inertia <- model$tot.chi
summary_table <- data.frame(
  metric = c("input_samples", "analyzed_samples", "input_features", "retained_features", "total_inertia", "constrained_inertia", "constrained_fraction", "constrained_axes"),
  value = c(ncol(feature_frame) - 1, nrow(metadata), nrow(feature_frame), sum(keep), total_inertia, model$CCA$tot.chi, model$CCA$tot.chi / total_inertia, model$CCA$rank),
  stringsAsFactors = FALSE
)
write_tsv(summary_table, args[["model-summary"]])

dir.create(dirname(args$plot), recursive = TRUE, showWarnings = FALSE)
pdf(args$plot, width = 8, height = 7)
plot(model, scaling = scaling, type = "n", main = "Constrained ordination (RDA)")
points(model, display = "sites", scaling = scaling, pch = 19, col = "#2c7fb8")
text(model, display = "bp", scaling = scaling, col = "#d7301f", cex = 0.8)
dev.off()

log_lines <- c(
  sprintf("vegan_version\t%s", as.character(utils::packageVersion("vegan"))),
  sprintf("formula\t%s", paste(deparse(model_formula), collapse = "")),
  sprintf("transform\t%s", args$transform),
  sprintf("samples_analyzed\t%d", nrow(metadata)),
  sprintf("features_retained\t%d", sum(keep)),
  sprintf("permutations\t%d", permutations),
  sprintf("random_seed\t%d", seed),
  sprintf("scaling\t%d", scaling),
  sprintf("constrained_eigenvalues\t%s", paste(format(constrained_eigen, digits = 17), collapse = ","))
)
dir.create(dirname(args$log), recursive = TRUE, showWarnings = FALSE)
writeLines(log_lines, args$log)
