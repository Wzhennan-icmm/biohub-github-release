#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: biohub run psmc-plot --input <merged.tsv> --output <plot.pdf|png>",
  "[--x-scale linear|log10] [--y-scale linear|log10] [--stages <stages.tsv>]",
  "[--width <inches>] [--height <inches>] [--force]"
)

if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
  cat(usage, "\n")
  quit(status = if (length(args) == 0) 2 else 0)
}

known_values <- c("--input", "--output", "--x-scale", "--y-scale", "--stages", "--width", "--height")
known_flags <- "--force"
values <- list(); force <- FALSE; index <- 1
while (index <= length(args)) {
  flag <- args[[index]]
  if (flag %in% known_flags) {
    if (force) stop("--force specified more than once", call. = FALSE)
    force <- TRUE; index <- index + 1; next
  }
  if (!flag %in% known_values || index == length(args)) stop(sprintf("unknown option or missing value: %s", flag), call. = FALSE)
  if (flag %in% names(values)) stop(sprintf("option specified more than once: %s", flag), call. = FALSE)
  values[[flag]] <- args[[index + 1]]; index <- index + 2
}

if (is.null(values[["--input"]]) || is.null(values[["--output"]])) stop("--input and --output are required", call. = FALSE)
input <- values[["--input"]]; output <- values[["--output"]]
x_scale <- if (is.null(values[["--x-scale"]])) "log10" else values[["--x-scale"]]
y_scale <- if (is.null(values[["--y-scale"]])) "log10" else values[["--y-scale"]]
width <- if (is.null(values[["--width"]])) 8 else as.numeric(values[["--width"]])
height <- if (is.null(values[["--height"]])) 5 else as.numeric(values[["--height"]])
if (!x_scale %in% c("linear", "log10") || !y_scale %in% c("linear", "log10")) stop("scales must be linear or log10", call. = FALSE)
if (!is.finite(width) || !is.finite(height) || width <= 0 || height <= 0) stop("width and height must be positive", call. = FALSE)
if (!file.exists(input)) stop(sprintf("input does not exist: %s", input), call. = FALSE)
if (file.exists(output) && !force) stop(sprintf("output exists: %s; pass --force to replace", output), call. = FALSE)

data <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("Sample", "Time", "Ne")
if (nrow(data) == 0 || any(!required %in% names(data))) stop("input needs Sample, Time and Ne columns", call. = FALSE)
data$Time <- suppressWarnings(as.numeric(data$Time)); data$Ne <- suppressWarnings(as.numeric(data$Ne))
if (anyNA(data$Sample) || any(!is.finite(data$Time)) || any(!is.finite(data$Ne)) || any(data$Sample == "")) stop("Sample, Time and Ne must be non-empty finite values", call. = FALSE)
if (x_scale == "log10" && any(data$Time <= 0)) stop("log10 x scale requires positive Time values", call. = FALSE)
if (y_scale == "log10" && any(data$Ne <= 0)) stop("log10 y scale requires positive Ne values", call. = FALSE)

transform_axis <- function(value, scale) if (scale == "log10") log10(value) else value
data$x <- transform_axis(data$Time, x_scale); data$y <- transform_axis(data$Ne, y_scale)
samples <- unique(data$Sample)
colors <- setNames(grDevices::hcl.colors(length(samples), "Dark 3"), samples)

stages <- NULL
if (!is.null(values[["--stages"]])) {
  if (!file.exists(values[["--stages"]])) stop("stages file does not exist", call. = FALSE)
  stages <- read.delim(values[["--stages"]], check.names = FALSE, stringsAsFactors = FALSE)
  stage_required <- c("label", "start", "end", "color")
  if (nrow(stages) == 0 || any(!stage_required %in% names(stages))) stop("stages needs label, start, end and color columns", call. = FALSE)
  stages$start <- suppressWarnings(as.numeric(stages$start)); stages$end <- suppressWarnings(as.numeric(stages$end))
  if (anyNA(stages[, c("label", "color"), drop = FALSE]) || any(stages$label == "") || any(stages$color == "") || any(!is.finite(stages$start)) || any(!is.finite(stages$end)) || any(stages$end <= stages$start)) stop("stage rows need labels/colors and finite intervals with end > start", call. = FALSE)
  if (x_scale == "log10" && any(stages$start <= 0)) stop("log10 x scale requires positive stage starts", call. = FALSE)
  invisible(lapply(stages$color, grDevices::col2rgb))
}

extension <- tolower(tools::file_ext(output))
if (extension == "pdf") {
  pdf(output, width = width, height = height)
} else if (extension == "png") {
  png(output, width = round(width * 220), height = round(height * 220), res = 220)
} else {
  stop("output extension must be .pdf or .png", call. = FALSE)
}

x_label <- if (x_scale == "log10") "Time (log10 scale)" else "Time"
y_label <- if (y_scale == "log10") "Effective population size (log10 scale)" else "Effective population size"
plot(range(data$x), range(data$y), type = "n", xlab = x_label, ylab = y_label, main = "PSMC demographic trajectories")
if (!is.null(stages)) {
  bounds <- par("usr")
  for (row in seq_len(nrow(stages))) {
    left <- transform_axis(stages$start[[row]], x_scale); right <- transform_axis(stages$end[[row]], x_scale)
    rect(left, bounds[[3]], right, bounds[[4]], col = grDevices::adjustcolor(stages$color[[row]], alpha.f = 0.12), border = NA)
    text((left + right) / 2, bounds[[4]], labels = stages$label[[row]], adj = c(0.5, 1.2), cex = 0.65)
  }
}
for (sample in samples) {
  selected <- data[data$Sample == sample, , drop = FALSE]
  selected <- selected[order(selected$Time), , drop = FALSE]
  lines(selected$x, selected$y, col = colors[[sample]], lwd = 1.2)
}
legend("topright", legend = samples, col = colors[samples], lty = 1, lwd = 1.2, bty = "n", cex = 0.75)
dev.off()
