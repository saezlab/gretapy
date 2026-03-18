from __future__ import annotations

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd

mu.set_options(pull_on_update=False)

# Hardcoded GRN with hg38-style coordinates
# Some TF-target pairs have multiple CREs (e.g., PAX5->CD19, GATA3->CD3E, SPI1->CD14)
_GRN = pd.DataFrame(
    {
        "source": [
            # PAX5 targets (B cell) - CD19 has 2 CREs
            "PAX5",
            "PAX5",
            "PAX5",
            "PAX5",
            "PAX5",
            "PAX5",
            # EBF1 targets (B cell) - CD19 has 2 CREs
            "EBF1",
            "EBF1",
            "EBF1",
            "EBF1",
            "EBF1",
            "EBF1",
            # GATA3 targets (T cell) - CD3E has 2 CREs
            "GATA3",
            "GATA3",
            "GATA3",
            "GATA3",
            "GATA3",
            "GATA3",
            # TCF7 targets (T cell) - CD3E has 2 CREs
            "TCF7",
            "TCF7",
            "TCF7",
            "TCF7",
            "TCF7",
            "TCF7",
            # SPI1 targets (Monocyte) - CD14 has 2 CREs
            "SPI1",
            "SPI1",
            "SPI1",
            "SPI1",
            "SPI1",
            "SPI1",
            # CEBPA targets (Monocyte) - CD14 has 2 CREs
            "CEBPA",
            "CEBPA",
            "CEBPA",
            "CEBPA",
            "CEBPA",
            "CEBPA",
            # New hallmark-overlapping targets
            # B cell
            "PAX5",  # CD79B
            "EBF1",  # CD79B
            "PAX5",  # CCND2
            "EBF1",  # CCND2
            "PAX5",  # BANK1
            # T cell
            "GATA3",  # GZMA
            "TCF7",  # GZMA
            "GATA3",  # IL2RA
            "TCF7",  # IL2RA
            "GATA3",  # CD38
            # Monocyte
            "SPI1",  # IL1B
            "CEBPA",  # IL1B
            "SPI1",  # CXCL10
            "CEBPA",  # CXCL10
            "SPI1",  # PIM1
            "CEBPA",  # PIM1
        ],
        "target": [
            # PAX5 targets (CD19 with 2 CREs)
            "CD19",
            "CD19",
            "MS4A1",
            "CD79A",
            "BCL2",
            "IRF4",
            # EBF1 targets (CD19 with 2 CREs)
            "CD19",
            "CD19",
            "MS4A1",
            "CD79A",
            "IRF4",
            "VPREB1",
            # GATA3 targets (CD3E with 2 CREs)
            "CD3E",
            "CD3E",
            "IL7R",
            "LEF1",
            "RUNX1",
            "IRF4",
            # TCF7 targets (CD3E with 2 CREs)
            "CD3E",
            "CD3E",
            "IL7R",
            "LEF1",
            "RUNX1",
            "BCL11B",
            # SPI1 targets (CD14 with 2 CREs)
            "CD14",
            "CD14",
            "CD68",
            "CSF1R",
            "IRF4",
            "RUNX1",
            # CEBPA targets (CD14 with 2 CREs)
            "CD14",
            "CD14",
            "CD68",
            "CSF1R",
            "RUNX1",
            "LYZ",
            # New hallmark-overlapping targets
            # B cell
            "CD79B",
            "CD79B",
            "CCND2",
            "CCND2",
            "BANK1",
            # T cell
            "GZMA",
            "GZMA",
            "IL2RA",
            "IL2RA",
            "CD38",
            # Monocyte
            "IL1B",
            "IL1B",
            "CXCL10",
            "CXCL10",
            "PIM1",
            "PIM1",
        ],
        "cre": [
            # PAX5 CREs (2 for CD19)
            "chr16-28926950-28927450",
            "chr16-28936000-28936500",
            "chr11-60443809-60444309",
            "chr19-41827233-41827733",
            "chr18-63073346-63073846",
            "chr6-366739-367239",
            # EBF1 CREs (2 for CD19, some shared with PAX5)
            "chr16-28926950-28927450",
            "chr16-28936000-28936500",
            "chr11-60443809-60444309",
            "chr19-41832233-41832733",
            "chr6-366739-367239",
            "chr22-22497505-22498005",
            # GATA3 CREs (2 for CD3E)
            "chr11-118336316-118336816",
            "chr11-118296000-118296500",
            "chr5-35861274-35861774",
            "chr4-108037038-108037538",
            "chr21-34772801-34773301",
            "chr6-421739-422239",
            # TCF7 CREs (2 for CD3E, some shared with GATA3)
            "chr11-118336316-118336816",
            "chr11-118296000-118296500",
            "chr5-35861274-35861774",
            "chr4-108042038-108042538",
            "chr21-34772801-34773301",
            "chr14-99159874-99160374",
            # SPI1 CREs (2 for CD14)
            "chr5-140628728-140629228",
            "chr5-140638000-140638500",
            "chr17-7573628-7574128",
            "chr5-150043291-150043791",
            "chr6-431739-432239",
            "chr21-34777801-34778301",
            # CEBPA CREs (2 for CD14, some shared with SPI1)
            "chr5-140628728-140629228",
            "chr5-140638000-140638500",
            "chr17-7573628-7574128",
            "chr5-150048291-150048791",
            "chr21-34777801-34778301",
            "chr12-69338350-69338850",
            # New hallmark-overlapping CREs
            # B cell (shared CREs for PAX5+EBF1 targets)
            "chr17-63940000-63940500",  # CD79B (PAX5)
            "chr17-63940000-63940500",  # CD79B (EBF1)
            "chr12-4280000-4280500",  # CCND2 (PAX5)
            "chr12-4280000-4280500",  # CCND2 (EBF1)
            "chr4-101420000-101420500",  # BANK1 (PAX5)
            # T cell
            "chr5-55112000-55112500",  # GZMA (GATA3)
            "chr5-55112000-55112500",  # GZMA (TCF7)
            "chr10-6052000-6052500",  # IL2RA (GATA3)
            "chr10-6052000-6052500",  # IL2RA (TCF7)
            "chr4-15788000-15788500",  # CD38 (GATA3)
            # Monocyte (shared CREs for SPI1+CEBPA targets)
            "chr2-112846000-112846500",  # IL1B (SPI1)
            "chr2-112846000-112846500",  # IL1B (CEBPA)
            "chr4-76033000-76033500",  # CXCL10 (SPI1)
            "chr4-76033000-76033500",  # CXCL10 (CEBPA)
            "chr6-37180000-37180500",  # PIM1 (SPI1)
            "chr6-37180000-37180500",  # PIM1 (CEBPA)
        ],
        "score": [
            0.85,
            0.72,
            0.78,
            0.72,
            0.68,
            0.81,
            0.82,
            0.70,
            0.75,
            0.69,
            0.79,
            0.71,
            0.88,
            0.74,
            0.76,
            0.73,
            0.67,
            0.74,
            0.86,
            0.71,
            0.74,
            0.70,
            0.65,
            0.77,
            0.84,
            0.73,
            0.79,
            0.71,
            0.76,
            0.69,
            0.81,
            0.69,
            0.77,
            0.68,
            0.66,
            0.73,
            # New hallmark-overlapping scores
            0.74,  # PAX5->CD79B
            0.71,  # EBF1->CD79B
            0.69,  # PAX5->CCND2
            0.67,  # EBF1->CCND2
            0.72,  # PAX5->BANK1
            0.75,  # GATA3->GZMA
            0.70,  # TCF7->GZMA
            0.77,  # GATA3->IL2RA
            0.73,  # TCF7->IL2RA
            0.68,  # GATA3->CD38
            0.80,  # SPI1->IL1B
            0.76,  # CEBPA->IL1B
            0.78,  # SPI1->CXCL10
            0.74,  # CEBPA->CXCL10
            0.71,  # SPI1->PIM1
            0.69,  # CEBPA->PIM1
        ],
    }
)

# TF-celltype mapping for expression patterns
_TF_CELLTYPE = {
    "PAX5": "B cell",
    "EBF1": "B cell",
    "GATA3": "T cell",
    "TCF7": "T cell",
    "SPI1": "Monocyte",
    "CEBPA": "Monocyte",
}

# Target genes per TF for expression patterns
_TF_TARGETS = {
    "PAX5": ["CD19", "MS4A1", "CD79A", "BCL2", "IRF4", "CD79B", "CCND2", "BANK1"],
    "EBF1": ["CD19", "MS4A1", "CD79A", "IRF4", "VPREB1", "CD79B", "CCND2"],
    "GATA3": ["CD3E", "IL7R", "LEF1", "RUNX1", "IRF4", "GZMA", "IL2RA", "CD38"],
    "TCF7": ["CD3E", "IL7R", "LEF1", "RUNX1", "BCL11B", "GZMA", "IL2RA"],
    "SPI1": ["CD14", "CD68", "CSF1R", "IRF4", "RUNX1", "IL1B", "CXCL10", "PIM1"],
    "CEBPA": ["CD14", "CD68", "CSF1R", "RUNX1", "LYZ", "IL1B", "CXCL10", "PIM1"],
}

# All genes (TFs + unique targets)
_TFS = ["PAX5", "EBF1", "GATA3", "TCF7", "SPI1", "CEBPA"]
_TARGETS = [
    "CD19",
    "MS4A1",
    "CD79A",
    "BCL2",
    "IRF4",
    "VPREB1",
    "CD3E",
    "IL7R",
    "LEF1",
    "RUNX1",
    "BCL11B",
    "CD14",
    "CD68",
    "CSF1R",
    "LYZ",
    # Hallmark-overlapping targets
    "CD79B",
    "CCND2",
    "BANK1",
    "GZMA",
    "IL2RA",
    "CD38",
    "IL1B",
    "CXCL10",
    "PIM1",
]
_GENES = _TFS + _TARGETS

# Hardcoded peaks: GRN CREs + promoters for all genes
_PEAKS = [
    # GRN CREs (unique)
    "chr16-28926950-28927450",
    "chr16-28936000-28936500",
    "chr11-60443809-60444309",
    "chr19-41827233-41827733",
    "chr18-63073346-63073846",
    "chr6-366739-367239",
    "chr19-41832233-41832733",
    "chr22-22497505-22498005",
    "chr11-118336316-118336816",
    "chr11-118296000-118296500",
    "chr5-35861274-35861774",
    "chr4-108037038-108037538",
    "chr21-34772801-34773301",
    "chr6-421739-422239",
    "chr4-108042038-108042538",
    "chr14-99159874-99160374",
    "chr5-140628728-140629228",
    "chr5-140638000-140638500",
    "chr17-7573628-7574128",
    "chr5-150043291-150043791",
    "chr6-431739-432239",
    "chr21-34777801-34778301",
    "chr5-150048291-150048791",
    "chr12-69338350-69338850",
    # New hallmark-overlapping CREs (unique)
    "chr17-63940000-63940500",  # CD79B
    "chr12-4280000-4280500",  # CCND2
    "chr4-101420000-101420500",  # BANK1
    "chr5-55112000-55112500",  # GZMA
    "chr10-6052000-6052500",  # IL2RA
    "chr4-15788000-15788500",  # CD38
    "chr2-112846000-112846500",  # IL1B
    "chr4-76033000-76033500",  # CXCL10
    "chr6-37180000-37180500",  # PIM1
    # Promoters for TFs (within 500bp upstream of TSS)
    "chr9-37034268-37034768",  # PAX5 (TSS=37034268, -)
    "chr5-159099916-159100416",  # EBF1 (TSS=159099916, -)
    "chr10-8044878-8045378",  # GATA3 (TSS=8045378, +)
    "chr5-134114181-134114681",  # TCF7 (TSS=134114681, +)
    "chr11-47378547-47379047",  # SPI1 (TSS=47378547, -)
    "chr19-33302534-33303034",  # CEBPA (TSS=33302534, -)
    # Promoters for targets
    "chr16-28931465-28931965",  # CD19 (TSS=28931965, +)
    "chr11-60455346-60455846",  # MS4A1 (TSS=60455846, +)
    "chr19-41876779-41877279",  # CD79A (TSS=41877279, +)
    "chr18-63320128-63320628",  # BCL2 (TSS=63320128, -)
    "chr6-391239-391739",  # IRF4 (TSS=391739, +)
    "chr22-22244280-22244780",  # VPREB1 (TSS=22244780, +)
    "chr11-118304230-118304730",  # CD3E (TSS=118304730, +)
    "chr5-35852195-35852695",  # IL7R (TSS=35852695, +)
    "chr4-108168956-108169456",  # LEF1 (TSS=108168956, -)
    "chr21-36004667-36005167",  # RUNX1 (TSS=36004667, -)
    "chr14-99272197-99272697",  # BCL11B (TSS=99272197, -)
    "chr5-140633700-140634200",  # CD14 (TSS=140633700, -)
    "chr17-7578991-7579491",  # CD68 (TSS=7579491, +)
    "chr5-150113372-150113872",  # CSF1R (TSS=150113372, -)
    "chr12-69347881-69348381",  # LYZ (TSS=69348381, +)
    # Promoters for new hallmark-overlapping targets
    "chr17-63932629-63933129",  # CD79B (TSS=63932629, -)
    "chr12-4269271-4269771",  # CCND2 (TSS=4269771, +)
    "chr4-101410786-101411286",  # BANK1 (TSS=101411286, +)
    "chr5-55102146-55102646",  # GZMA (TSS=55102646, +)
    "chr10-6062370-6062870",  # IL2RA (TSS=6062370, -)
    "chr4-15777775-15778275",  # CD38 (TSS=15778275, +)
    "chr2-112836816-112837316",  # IL1B (TSS=112836816, -)
    "chr4-76023497-76023997",  # CXCL10 (TSS=76023497, -)
    "chr6-37169652-37170152",  # PIM1 (TSS=37170152, +)
]

# Map CREs to celltypes for accessibility patterns
_CRE_CELLTYPE = {
    # B cell CREs
    "chr16-28926950-28927450": "B cell",
    "chr16-28936000-28936500": "B cell",
    "chr11-60443809-60444309": "B cell",
    "chr19-41827233-41827733": "B cell",
    "chr18-63073346-63073846": "B cell",
    "chr6-366739-367239": "B cell",
    "chr19-41832233-41832733": "B cell",
    "chr22-22497505-22498005": "B cell",
    # T cell CREs
    "chr11-118336316-118336816": "T cell",
    "chr11-118296000-118296500": "T cell",
    "chr5-35861274-35861774": "T cell",
    "chr4-108037038-108037538": "T cell",
    "chr21-34772801-34773301": "T cell",
    "chr6-421739-422239": "T cell",
    "chr4-108042038-108042538": "T cell",
    "chr14-99159874-99160374": "T cell",
    # Monocyte CREs
    "chr5-140628728-140629228": "Monocyte",
    "chr5-140638000-140638500": "Monocyte",
    "chr17-7573628-7574128": "Monocyte",
    "chr5-150043291-150043791": "Monocyte",
    "chr6-431739-432239": "Monocyte",
    "chr21-34777801-34778301": "Monocyte",
    "chr5-150048291-150048791": "Monocyte",
    "chr12-69338350-69338850": "Monocyte",
    # New hallmark-overlapping CREs
    "chr17-63940000-63940500": "B cell",  # CD79B
    "chr12-4280000-4280500": "B cell",  # CCND2
    "chr4-101420000-101420500": "B cell",  # BANK1
    "chr5-55112000-55112500": "T cell",  # GZMA
    "chr10-6052000-6052500": "T cell",  # IL2RA
    "chr4-15788000-15788500": "T cell",  # CD38
    "chr2-112846000-112846500": "Monocyte",  # IL1B
    "chr4-76033000-76033500": "Monocyte",  # CXCL10
    "chr6-37180000-37180500": "Monocyte",  # PIM1
}

# Map gene promoters to celltypes
_PROMOTER_CELLTYPE = {
    # B cell TFs
    "chr9-37034268-37034768": "B cell",  # PAX5
    "chr5-159099916-159100416": "B cell",  # EBF1
    # B cell targets
    "chr16-28931465-28931965": "B cell",  # CD19
    "chr11-60455346-60455846": "B cell",  # MS4A1
    "chr19-41876779-41877279": "B cell",  # CD79A
    "chr18-63320128-63320628": "B cell",  # BCL2
    "chr22-22244280-22244780": "B cell",  # VPREB1
    # T cell TFs
    "chr10-8044878-8045378": "T cell",  # GATA3
    "chr5-134114181-134114681": "T cell",  # TCF7
    # T cell targets
    "chr11-118304230-118304730": "T cell",  # CD3E
    "chr5-35852195-35852695": "T cell",  # IL7R
    "chr4-108168956-108169456": "T cell",  # LEF1
    "chr14-99272197-99272697": "T cell",  # BCL11B
    # Monocyte TFs
    "chr11-47378547-47379047": "Monocyte",  # SPI1
    "chr19-33302534-33303034": "Monocyte",  # CEBPA
    # Monocyte targets
    "chr5-140633700-140634200": "Monocyte",  # CD14
    "chr17-7578991-7579491": "Monocyte",  # CD68
    "chr5-150113372-150113872": "Monocyte",  # CSF1R
    "chr12-69347881-69348381": "Monocyte",  # LYZ
    # New hallmark-overlapping target promoters
    "chr17-63932629-63933129": "B cell",  # CD79B
    "chr12-4269271-4269771": "B cell",  # CCND2
    "chr4-101410786-101411286": "B cell",  # BANK1
    "chr5-55102146-55102646": "T cell",  # GZMA
    "chr10-6062370-6062870": "T cell",  # IL2RA
    "chr4-15777775-15778275": "T cell",  # CD38
    "chr2-112836816-112837316": "Monocyte",  # IL1B
    "chr4-76023497-76023997": "Monocyte",  # CXCL10
    "chr6-37169652-37170152": "Monocyte",  # PIM1
    # Hub genes (accessible in multiple celltypes)
    "chr6-391239-391739": None,  # IRF4
    "chr21-36004667-36005167": None,  # RUNX1
}


def toy(n_cells: int = 60, seed: int = 42) -> tuple[mu.MuData, pd.DataFrame]:
    """
    Generate synthetic MuData and GRN for testing and demonstration.

    Creates biologically structured test data with RNA and ATAC modalities,
    along with a gene regulatory network (GRN) DataFrame.

    Parameters
    ----------
    n_cells : int
        Number of cells to generate. Default is 60.
    seed : int
        Random seed for reproducibility. Default is 42.

    Returns
    -------
    tuple[mu.MuData, pd.DataFrame]
        A tuple containing:
        - MuData with 'rna' and 'atac' modalities
        - GRN DataFrame with columns: source, target, cre, score
    """
    np.random.seed(seed)

    celltypes = ["B cell", "T cell", "Monocyte"]
    n_celltypes = len(celltypes)

    # Create cell type assignments
    cells_per_type = n_cells // n_celltypes
    remainder = n_cells % n_celltypes
    celltype_list = []
    for i, ct in enumerate(celltypes):
        count = cells_per_type + (1 if i < remainder else 0)
        celltype_list.extend([ct] * count)

    # Create RNA expression matrix
    n_genes = len(_GENES)
    X_rna = np.random.exponential(scale=1.0, size=(n_cells, n_genes)).astype(np.float32)
    X_rna = np.clip(X_rna, 0, 4.0)

    # Add TF and target expression patterns per celltype
    for tf, ct in _TF_CELLTYPE.items():
        tf_idx = _GENES.index(tf)
        cell_mask = np.array([c == ct for c in celltype_list])
        X_rna[cell_mask, tf_idx] += np.random.uniform(3.0, 6.0, size=cell_mask.sum())

        for target in _TF_TARGETS[tf]:
            target_idx = _GENES.index(target)
            X_rna[cell_mask, target_idx] += np.random.uniform(2.0, 5.0, size=cell_mask.sum())

    X_rna = np.clip(X_rna, 0, 12.0)

    # Create RNA AnnData
    rna = ad.AnnData(X=X_rna)
    rna.var_names = _GENES
    rna.obs_names = [f"Cell{i}" for i in range(n_cells)]

    # Create ATAC accessibility matrix
    n_peaks = len(_PEAKS)
    X_atac = np.random.exponential(scale=0.8, size=(n_cells, n_peaks)).astype(np.float32)
    X_atac = np.clip(X_atac, 0, 3.0)

    peak_idx_map = {p: i for i, p in enumerate(_PEAKS)}

    # Add accessibility patterns for CREs
    for cre, ct in _CRE_CELLTYPE.items():
        if cre in peak_idx_map:
            peak_idx = peak_idx_map[cre]
            cell_mask = np.array([c == ct for c in celltype_list])
            X_atac[cell_mask, peak_idx] += np.random.uniform(2.0, 5.0, size=cell_mask.sum())

    # Add accessibility patterns for promoters
    for promoter, ct in _PROMOTER_CELLTYPE.items():
        if promoter in peak_idx_map:
            peak_idx = peak_idx_map[promoter]
            if ct is None:
                # Hub gene promoters accessible in all celltypes
                X_atac[:, peak_idx] += np.random.uniform(1.5, 3.0, size=n_cells)
            else:
                cell_mask = np.array([c == ct for c in celltype_list])
                X_atac[cell_mask, peak_idx] += np.random.uniform(2.0, 5.0, size=cell_mask.sum())

    X_atac = np.clip(X_atac, 0, 12.0)

    # Create ATAC AnnData
    atac = ad.AnnData(X=X_atac)
    atac.var_names = _PEAKS
    atac.obs_names = rna.obs_names.copy()

    # Create MuData
    mdata = mu.MuData({"rna": rna, "atac": atac})
    mdata.obs["celltype"] = celltype_list

    return mdata, _GRN.copy()
