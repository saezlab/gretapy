# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

## 0.0.2

### Added

- `ds.read_imaginary_metrics`: samples one row per unique benchmark configuration from the metrics table and relabels it as `'ImaginaryMethod'`, useful for baseline comparisons or testing visualizations.
- `get_terms`: retrieves the filtering terms for a given dataset as a dictionary mapping evaluation database names to lists of terms.
- `tl.cre_to_tss_distance`: computes the distance from each CRE in a GRN to the promoter window of its target gene, using the Promoters database. Accepts a single GRN DataFrame or a dictionary of GRNs.
- `pl.cre_to_tss_distance`: plots CRE-to-TSS distance distributions as horizontal boxplots, with an optional vertical threshold line (default 250,000 bp).
- Tutorial notebooks.
- Added `scanpy` and `xgboost` as explicit dependencies; removed the `[full]` extra from `decoupler`.

### Changed

- `pl.links`: legend is now displayed once on the middle TF panel (with multi-column layout for many GRNs) instead of being repeated on every subplot.
- `mt.correlation` now filters for promoters first and then computes correlations, reducing the number of compute.

## 0.0.1

### Added

- First release
