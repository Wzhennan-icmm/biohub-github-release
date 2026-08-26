# BioHub domain validation record: annotation-coordinates

- Inventory IDs: `014`, `015`, `016`
- Commands: `convert-gemoma-gff3`, `convert-gene-annotation-contigs2chr-PASA`, `convert-gene-annotation-scaffold2chr-nextgenomics`
- Reviewer: Zhennan Wang
- Reviewer account: `wangzhennan`
- Reviewer affiliation: 中国中医科学院中药研究所 / Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences
- Reviewer ORCID: https://orcid.org/0000-0003-4883-2538
- Review date: 2026-08-26
- Decision: approved
- BioHub head commit: `005958bc61285d7bf5b3d27ca1ce909c50d2e909`
- CI merge commit recorded by evidence: `3e6bb854d2f07488a88e504b18ea01fb58d331cd`
- Evidence pack and schema version: `annotation-coordinates`, schema `1`
- GitHub Actions run: https://github.com/Wzhennan-icmm/biohub-github-release/actions/runs/32813865044
- Artifact: `domain-review-annotation-coordinates` (ID `9550726788`)
- Input-manifest SHA256: `35bd30a667f689f934c723d82893f2aca8e32e8a0de76dead831760abb58f31f`
- Output-manifest SHA256: `baf1c1a542c7b71c6e2d0b42b01b95ac738a9d999c95f7474ed14ddb418687ed`
- Evidence-manifest SHA256: `2a8384230efa50158865e8d5cfa58cc6821b38b95bc3e62abccba051806028e0`
- Reference software and version: BioHub `0.4.0`; Python `3.12.3`; exact byte comparisons executed by `validation_review.py` schema `1`
- Reference command and parameters: `diff -u` for inventory IDs `014` and `015`; `diff -ru` for inventory ID `016`; commands and parameters retained in artifact `reference-commands.txt`
- Numeric or visual acceptance tolerance: exact byte equality for all seven declared output and audit files; coordinate semantics reviewed as 1-based closed intervals
- Differences found and resolution: none reported; all seven automated comparisons passed
- Data/license confirmation: public CC0-1.0 synthetic fixtures only; no private genomes, sample identifiers, or unpublished results
- Approval source: explicit maintainer instruction `批准 annotation-coordinates` on 2026-08-26

## Accepted scope

- Inventory `014`: GFF3 record order, feature types, `ID`/`Parent` hierarchy,
  phase, strand, and coordinates match reviewed expected output.
- Inventory `015`: forward and reverse mappings preserve interval length under
  documented 1-based closed-coordinate arithmetic.
- Inventory `016`: split-scaffold, out-of-range, and unmapped records remain
  visible in declared outputs or audit logs; no silent record loss was observed.

## Limits

Approval covers committed CC0 fixtures, documented input contracts, and tested
boundary cases. It does not claim correctness for untested GFF dialects,
coordinate conventions, malformed records, species, or mapping-table layouts.
