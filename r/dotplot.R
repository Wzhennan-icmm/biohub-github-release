#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: biohub r run dotplot --input <paf-or-coords> --output <pdf-or-png> [--format paf|coords]"

value_for <- function(flag, required = TRUE) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) {
    if (required) stop(sprintf("missing value for %s", flag), call. = FALSE)
    return(NULL)
  }
  args[[index + 1]]
}

if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
  cat(usage, "\n")
  quit(status = if (length(args) == 0) 2 else 0)
}

input <- value_for("--input")
output <- value_for("--output")
format <- value_for("--format", required = FALSE)
if (!file.exists(input)) stop(sprintf("input does not exist: %s", input), call. = FALSE)
if (file.exists(output) && !any(args == "--force")) {
  stop(sprintf("output exists: %s; pass --force to replace", output), call. = FALSE)
}
if (is.null(format)) {
  format <- if (grepl("\\.paf$", input, ignore.case = TRUE)) "paf" else "coords"
}
if (!format %in% c("paf", "coords")) stop("--format must be paf or coords", call. = FALSE)

rows <- readLines(input, warn = FALSE)
rows <- rows[nzchar(rows) & !startsWith(rows, "#")]
if (length(rows) == 0) stop("input has no alignment rows", call. = FALSE)

points <- lapply(rows, function(row) {
  fields <- if (format == "paf") {
    strsplit(row, "\\t", fixed = FALSE)[[1]]
  } else {
    strsplit(trimws(row), "[[:space:]|]+")[[1]]
  }
  if (format == "paf") {
    if (length(fields) < 9) return(NULL)
    list(
      x = (as.numeric(fields[[3]]) + as.numeric(fields[[4]])) / 2,
      y = (as.numeric(fields[[8]]) + as.numeric(fields[[9]])) / 2,
      strand = fields[[5]]
    )
  } else {
    numeric_fields <- suppressWarnings(as.numeric(fields))
    numeric_fields <- numeric_fields[!is.na(numeric_fields)]
    if (length(numeric_fields) < 4) return(NULL)
    list(
      x = (numeric_fields[[1]] + numeric_fields[[2]]) / 2,
      y = (numeric_fields[[3]] + numeric_fields[[4]]) / 2,
      strand = "+"
    )
  }
})
points <- Filter(Negate(is.null), points)
if (length(points) == 0) stop("no parseable alignment rows", call. = FALSE)

x <- vapply(points, `[[`, numeric(1), "x")
y <- vapply(points, `[[`, numeric(1), "y")
strand <- vapply(points, `[[`, character(1), "strand")
colors <- ifelse(strand == "-", "#d95f02", "#1b9e77")

extension <- tolower(tools::file_ext(output))
if (extension == "pdf") {
  pdf(output, width = 8, height = 8)
} else if (extension == "png") {
  png(output, width = 1800, height = 1800, res = 220)
} else {
  stop("output extension must be .pdf or .png", call. = FALSE)
}
plot(x, y, pch = 16, cex = 0.35, col = colors,
     xlab = "Query coordinate", ylab = "Target coordinate",
     main = sprintf("BioHub dot plot (%s)", toupper(format)))
legend("topright", legend = c("forward", "reverse"), pch = 16,
       col = c("#1b9e77", "#d95f02"), bty = "n")
dev.off()
