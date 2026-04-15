# Frequently Asked Questions

## How to evaluate a non-multimodal GRN (only TF-Gene interactions)?
While GRETA is designed to evaluate multimodal GRNs, unimodal GRNs (i.e., containing only TF-Gene information and no cis-regulatory element (CRE) information) can still be evaluated.

To do so, apply the same filtering strategy used in the baselines of the GRETA manuscript: TF-Gene interactions are filtered out if the promoter region of the target gene is not accessible (i.e., there is no accessibility peak in that genomic region). Run the following line once per dataset:

```python
mdata  # Multiome dataset as a MuData object
org  # Organism name, either 'hg38' or 'mm10'
grn  # GRN as a pandas DataFrame with columns: `source`, `target`, and optionally `weight`

mgrn = gt.mt.lit_grn(mdata=mdata, organism=org, grn=grn)
```

This function will: (1) download the promoter regions for the given organism, (2) check whether those regions are accessible in the multiome data, (3) filter TF-Gene interactions whose target gene promoter is inaccessible, and (4) annotate the remaining interactions with their promoter region as the CRE, returning a TF-CRE-Gene interaction network.

## How can I change where gretapy downloads datasets or databases?
By default, `gretapy` stores all downloaded data at `./gretapy_data`. This path is defined in the config and can be changed before running any gretapy functions:

```python
import gretapy as gt

# Update config before running gretapy functions
gt.config.PATH_DATA = 'path/to/new/location'  # Instead of './gretapy_data'

# Run gretapy code afterwards
```

## How much disk space does gretapy's data occupy?
When running the full benchmark, gretapy downloads a total of ~40 GB of data (datasets and evaluation databases combined). Make sure you have sufficient disk space before starting.

## For which organisms can I run the benchmark?
Currently gretapy supports human (`hg38`) and mouse (`mm10`). Other organisms are not yet supported, but this may change in the future. We are open to collaborations to add more model organisms.

## How can I run the benchmark with my custom multiome dataset?
Some evaluation databases are filtered by cell type and tissue terms. For a custom dataset, you will need to define those terms manually and provide them as a dictionary when calling `gt.tl.eval_grn_dataset` or `gt.tl.benchmark`. To see all available terms per database, run:

```python
all_terms = gt.show_terms(organism='hg38')
```

## Why do you use the overlap coefficient as a similarity metric when comparing GRNs instead of the Jaccard coefficient?

GRN inference methods produce networks of vastly different sizes, and the Jaccard coefficient is sensitive to this size difference, making it a poor choice for comparing GRNs. The overlap coefficient handles this more gracefully.

The two metrics are defined as follows. Given two sets of edges A and B:

- **Jaccard coefficient:** |A ∩ B| / |A ∪ B|
- **Overlap coefficient:** |A ∩ B| / min(|A|, |B|)

The overlap coefficient asks: *"What fraction of the smaller network's edges are recovered in the larger network?"* This is a more meaningful question when comparing methods that differ greatly in the number of edges they predict.

**Example:** Consider two methods applied to the same dataset:
- Method A (conservative): 1,000 edges
- Method B (permissive): 10,000 edges
- Shared edges: 900

| Metric | Calculation | Result |
|--------|-------------|--------|
| Jaccard | 900 / (1,000 + 10,000 − 900) | ≈ 0.09 |
| Overlap | 900 / min(1,000, 10,000) | = 0.90 |

The Jaccard score of 0.09 suggests almost no agreement, but this low value is driven entirely by the 10-fold size difference between the networks — not by any real disagreement in the regulatory logic they capture. The overlap score of 0.90 is more informative: 90% of the smaller network's edges are recovered in the larger one, indicating strong agreement in the interactions both methods consider most confident.
