# BioHub Slurm profile

Requires Snakemake 8.6 or newer plus `snakemake-executor-plugin-slurm`.
Packaged profile submits one job at a time. Copy directory and add cluster-approved
`jobs`, `default-resources`, account, partition, memory, and runtime settings before
production use. Pass copied directory with `--profile PATH`.

Reference: <https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html>
