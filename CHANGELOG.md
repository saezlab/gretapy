# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

## 0.0.2

### Added

- `tl.cre_to_tss_distance`: computes the distance from each CRE in a GRN to the promoter window of its target gene, using the Promoters database. Accepts a single GRN DataFrame or a dictionary of GRNs.
- `pl.cre_to_tss_distance`: plots CRE-to-TSS distance distributions as horizontal boxplots, with an optional vertical threshold line (default 250,000 bp).

## 0.0.1

### Added

- First release
