# R backend contract

R-backed commands live here after migration. Each entrypoint must:

- accept explicit input and output arguments through `optparse`;
- never call `install.packages`, `install_github`, `file.choose`, or `setwd`;
- write only declared output paths;
- return non-zero on validation or processing failure;
- declare package requirements in its command-registry metadata.

Native users provide R and required packages. Docker supplies `Rscript`; image
extensions add command-specific CRAN and Bioconductor packages. Use
`biohub doctor` before running an R-backed command.
