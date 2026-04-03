import anndata as ad
import numpy as np
import pandas as pd
import pyranges as pr
import pytest

import gretapy as gp


def get_tf_names():
    """Return common TF names used across fixtures (matches toy() defaults)."""
    return ["PAX5", "GATA3", "SPI1"]


def get_target_genes():
    """Return common target gene names used across fixtures (matches toy() defaults)."""
    # PAX5 targets: CD19, MS4A1, CD79A, BCL2, IRF4
    # GATA3 targets: CD3E, IL7R, TCF7, FOXP3, RUNX1
    # SPI1 targets: CD14, CD19, MS4A1, IRF4, RUNX1 (overlaps with PAX5/GATA3)
    return ["CD19", "MS4A1", "CD79A", "BCL2", "IRF4", "CD3E", "IL7R", "TCF7", "FOXP3", "RUNX1", "CD14"]


def get_all_genes():
    """Return all unique gene names (TFs + targets) used across fixtures."""
    return list(dict.fromkeys(get_tf_names() + get_target_genes()))


@pytest.fixture
def toy_data():
    """Base toy data fixture using gretapy.ds.toy()."""
    return gp.ds.toy(seed=42)


@pytest.fixture
def mudata_with_celltype(toy_data):
    """MuData object with RNA and ATAC modalities from toy data."""
    mdata, _ = toy_data
    return mdata


@pytest.fixture
def simple_grn(toy_data):
    """GRN DataFrame from toy data."""
    _, grn = toy_data
    return grn


@pytest.fixture
def adata(mudata_with_celltype):
    """Basic AnnData object extracted from MuData RNA modality."""
    rna = mudata_with_celltype.mod["rna"].copy()
    rna.obs["celltype"] = mudata_with_celltype.obs["celltype"].values
    return rna


@pytest.fixture
def reference_grn_db():
    """Reference GRN database DataFrame with shared gene names."""
    return pd.DataFrame(
        {
            "source": ["PAX5", "PAX5", "GATA3", "SPI1", "SPI1"],
            "target": ["CD19", "MS4A1", "IL7R", "CD14", "RUNX1"],
        }
    )


@pytest.fixture
def tfm_db():
    """TF markers database DataFrame (two columns: TF name and cell type)."""
    return pd.DataFrame(
        {
            0: ["PAX5", "GATA3", "SPI1", "RUNX1"],
            1: ["B cell", "T cell", "Monocyte", "Stem cell"],
        }
    )


@pytest.fixture
def tfp_db():
    """TF pairs database DataFrame."""
    return pd.DataFrame(
        {
            0: ["PAX5", "PAX5", "GATA3"],
            1: ["GATA3", "SPI1", "SPI1"],
        }
    )


@pytest.fixture
def pyranges_db():
    """PyRanges database for genomic tests with coordinates overlapping peak_names."""
    return pr.PyRanges(
        pd.DataFrame(
            {
                "Chromosome": [
                    "chr16",
                    "chr16",
                    "chr11",
                    "chr11",
                    "chr11",
                    "chr5",
                    "chr5",
                    "chr19",
                    "chr18",
                    "chr6",
                    "chr1",
                    "chr5",
                    "chrX",
                    "chr21",
                ],
                "Start": [
                    28931100,
                    28932100,
                    60223100,
                    60224100,
                    118209100,
                    35871100,
                    140013100,
                    41879100,
                    63318100,
                    411100,
                    207940100,
                    134110100,
                    49250100,
                    34799100,
                ],
                "End": [
                    28931400,
                    28932400,
                    60223400,
                    60224400,
                    118209400,
                    35871400,
                    140013400,
                    41879400,
                    63318400,
                    411400,
                    207940400,
                    134110400,
                    49250400,
                    34799400,
                ],
                "Name": [
                    "CD19",
                    "CD19",
                    "MS4A1",
                    "MS4A1",
                    "CD3E",
                    "IL7R",
                    "CD14",
                    "CD79A",
                    "BCL2",
                    "IRF4",
                    "CD34",
                    "TCF7",
                    "FOXP3",
                    "RUNX1",
                ],
                "Score": [
                    "B cell",
                    "B cell",
                    "B cell",
                    "B cell",
                    "T cell",
                    "T cell",
                    "Monocyte",
                    "B cell",
                    "B cell",
                    "B cell",
                    "Stem cell",
                    "T cell",
                    "T cell",
                    "Stem cell",
                ],
            }
        )
    )


@pytest.fixture
def gene_set_db():
    """Gene set database for pathway enrichment.

    Contains genes that are both present and absent in adata to test filtering.
    Each pathway has enough genes to pass decoupler's tmin threshold.
    """
    return pd.DataFrame(
        {
            "source": [
                # B cell activation pathway (7 genes: 6 present, 1 absent)
                "B_CELL_ACTIVATION",
                "B_CELL_ACTIVATION",
                "B_CELL_ACTIVATION",
                "B_CELL_ACTIVATION",
                "B_CELL_ACTIVATION",
                "B_CELL_ACTIVATION",
                "B_CELL_ACTIVATION",
                # T cell activation pathway (7 genes: 6 present, 1 absent)
                "T_CELL_ACTIVATION",
                "T_CELL_ACTIVATION",
                "T_CELL_ACTIVATION",
                "T_CELL_ACTIVATION",
                "T_CELL_ACTIVATION",
                "T_CELL_ACTIVATION",
                "T_CELL_ACTIVATION",
                # Myeloid pathway (6 genes: 5 present, 1 absent)
                "MYELOID_DIFFERENTIATION",
                "MYELOID_DIFFERENTIATION",
                "MYELOID_DIFFERENTIATION",
                "MYELOID_DIFFERENTIATION",
                "MYELOID_DIFFERENTIATION",
                "MYELOID_DIFFERENTIATION",
            ],
            "target": [
                # B cell genes (CD19, MS4A1, CD79A, PAX5, BCL2, IRF4 present; EBF1 absent)
                "CD19",
                "MS4A1",
                "CD79A",
                "PAX5",
                "BCL2",
                "IRF4",
                "EBF1",
                # T cell genes (CD3E, IL7R, GATA3, TCF7, FOXP3, RUNX1 present; LEF1 absent)
                "CD3E",
                "IL7R",
                "GATA3",
                "TCF7",
                "FOXP3",
                "RUNX1",
                "LEF1",
                # Myeloid genes (CD14, SPI1, IRF4, MS4A1, RUNX1 present; CSF1R absent)
                "CD14",
                "SPI1",
                "IRF4",
                "MS4A1",
                "RUNX1",
                "CSF1R",
            ],
            "weight": [1.0] * 20,
        }
    )


@pytest.fixture
def knocktf_db():
    """Mock KnockTF database for TF scoring and forecasting tests."""
    np.random.seed(42)
    all_genes = get_all_genes()
    n_experiments = 5
    exp_names = [f"Exp{i}" for i in range(n_experiments)]

    X = np.random.randn(n_experiments, len(all_genes)).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.var_names = all_genes
    adata.obs_names = exp_names
    adata.obs["TF"] = ["PAX5", "GATA3", "PAX5", "SPI1", "GATA3"]
    adata.obs["logFC"] = [-1.0, -0.8, -1.2, -0.6, -0.9]
    adata.obs["Tissue.Type"] = ["Blood", "Blood", "Brain", "Blood", "Brain"]

    return adata
