# gretapy - Evaluation and analysis of Gene Regulatory Networks (GRNs)
<img src="https://raw.githubusercontent.com/saezlab/gretapy/refs/heads/main/docs/_static/images/logo.png" align="right" width="120" class="no-scaled-link" alt='GRETA logo' />

[![Tests][badge-tests]][tests]
[![Documentation][badge-docs]][documentation]

[![Issues][badge-issues]][issue tracker]
[![Coverage][badge-coverage]][codecoverage]
[![Stars][badge-stars]](https://github.com/scverse/anndata/stargazers)

[![PyPI][badge-pypi]][pypi]
[![Downloads month][badge-mdown]][down]
[![Downloads all][badge-adown]][down]

[![Conda version][badge-condav]][conda]
[![Conda downloads][badge-condad]][conda]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/saezlab/gretapy/test.yaml?branch=main
[badge-docs]: https://img.shields.io/readthedocs/gretapy
[badge-condav]: https://img.shields.io/conda/vn/conda-forge/gretapy.svg
[badge-condad]: https://img.shields.io/conda/dn/conda-forge/gretapy.svg
[badge-issues]: https://img.shields.io/github/issues/saezlab/gretapy
[badge-coverage]: https://codecov.io/gh/saezlab/gretapy/branch/main/graph/badge.svg
[badge-pypi]: https://img.shields.io/pypi/v/gretapy.svg
[badge-mdown]: https://static.pepy.tech/badge/gretapy/month
[badge-adown]: https://static.pepy.tech/badge/gretapy
[badge-stars]: https://img.shields.io/github/stars/saezlab/gretapy?style=flat&logo=github&color=yellow

`gretapy` is a comprehensive framework for benchmarking and evaluating gene regulatory networks (GRNs) inferred from single-cell multiome (RNA+ATAC) data. It provides a systematic evaluation across four complementary dimensions: prior knowledge validation (TF markers, known TF-TF interactions, reference networks), genomic annotations (TF binding sites, cis-regulatory elements, chromatin-gene links), predictive performance (pathway enrichment, expression correlation), and mechanistic validation (perturbation forecasting, Boolean network simulations). The package includes built-in GRN inference methods, curated benchmark datasets, and visualization tools to facilitate rigorous comparison of network inference approaches.

<img src="https://raw.githubusercontent.com/saezlab/gretapy/refs/heads/main/docs/_static/images/gabstract.png" align="center" width="100%" alt='GRETA graphical abstract' />

## Getting started

Please refer to the [documentation][],
in particular, the [API documentation][].

## Installation

You need to have Python 3.11 or newer installed on your system.
If you don't have Python installed, we recommend installing [uv][].

There are several alternative options to install gretapy:

1. Install the latest stable release from [PyPI][pypi] with minimal dependancies:

```bash
pip install gretapy
```

2. Install the latest stable full release from [PyPI][pypi] with extra dependancies:

```bash
pip install gretapy[full]
```

3. Install the latest stable version from [conda-forge][conda] using mamba or conda:

```bash
mamba create -n=greta conda-forge::gretapy
```

4. Install the latest development version:

```bash
pip install git+https://github.com/saezlab/gretapy.git@main
```

## Release notes

See the [changelog][].

## Contact

For questions and help requests, you can reach out in the [scverse discourse][].
If you found a bug, please use the [issue tracker][].

## Citation

> Badia-i-Mompel P., Casals-Franch R., Wessels L., Müller-Dott S., Trimbour R., Yang Y., Ramirez Flores R.O., Saez-Rodriguez J. 2024. Comparison and evaluation of methods to infer gene regulatory networks from multimodal single-cell data. bioRxiv. https://doi.org/10.1101/2024.12.20.629764

[uv]: https://github.com/astral-sh/uv
[scverse discourse]: https://discourse.scverse.org/
[issue tracker]: https://github.com/saezlab/gretapy/issues
[tests]: https://github.com/saezlab/gretapy/actions/workflows/test.yaml
[documentation]: https://gretapy.readthedocs.io
[changelog]: https://gretapy.readthedocs.io/en/latest/changelog.html
[api documentation]: https://gretapy.readthedocs.io/en/latest/api.html
[pypi]: https://pypi.org/project/gretapy
[down]: https://pepy.tech/project/gretapy
[conda]: https://anaconda.org/conda-forge/gretapy
[codecoverage]: https://codecov.io/gh/saezlab/gretapy
