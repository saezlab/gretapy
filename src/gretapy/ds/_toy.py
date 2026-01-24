from __future__ import annotations

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from decoupler._download import _log

# Default TF-celltype associations for biologically structured expression
_TF_CELLTYPE = {
    "PAX5": "B cell",
    "EBF1": "B cell",  # Co-regulates B cell genes with PAX5
    "GATA3": "T cell",
    "TCF7": "T cell",  # Co-regulates T cell genes with GATA3
    "SPI1": "Monocyte",
    "CEBPA": "Monocyte",  # Co-regulates myeloid genes with SPI1
}

# Default target genes per TF (immune-relevant, with overlap between TFs)
# TFs that regulate the same celltype share target genes
# Some "hub" genes (IRF4, RUNX1) are regulated by TFs across celltypes
_TF_TARGETS = {
    "PAX5": ["CD19", "MS4A1", "CD79A", "BCL2", "IRF4"],
    "EBF1": ["CD19", "MS4A1", "CD79A", "IRF4", "VPREB1"],  # Overlaps with PAX5, shares IRF4
    "GATA3": ["CD3E", "IL7R", "TCF7", "RUNX1", "IRF4"],  # IRF4 is a hub gene
    "TCF7": ["CD3E", "IL7R", "LEF1", "RUNX1", "BCL11B"],  # Overlaps with GATA3
    "SPI1": ["CD14", "CD68", "CSF1R", "IRF4", "RUNX1"],  # IRF4, RUNX1 are hubs
    "CEBPA": ["CD14", "CD68", "CSF1R", "RUNX1", "IRF4"],  # Overlaps with SPI1, shares hubs
}

# CREs that are shared between TFs of the same celltype (co-binding)
# Format: {target_gene: {celltype: CRE_offset_from_TSS}}
# TFs of the same celltype will share this CRE in addition to their unique ones
_SHARED_CRES = {
    "CD19": {"B cell": -5000},  # PAX5 and EBF1 co-bind near CD19 promoter
    "MS4A1": {"B cell": -12000},  # Shared B cell enhancer
    "CD3E": {"T cell": -8000},  # GATA3 and TCF7 co-bind
    "IL7R": {"T cell": 15000},  # Shared T cell enhancer
    "CD14": {"Monocyte": -3000},  # SPI1 and CEBPA co-bind at promoter
    "CD68": {"Monocyte": -20000},  # Shared myeloid enhancer
    "IRF4": {"B cell": -25000, "T cell": 30000, "Monocyte": -40000},  # Hub gene, multiple enhancers
    "RUNX1": {"T cell": -15000, "Monocyte": 25000},  # Hub gene
}

# TSS positions for genes (hg38 approximate coordinates)
# Format: (chromosome, TSS position)
_GENE_TSS = {
    # B cell genes
    "CD19": ("chr16", 28931950),
    "MS4A1": ("chr11", 60455809),
    "CD79A": ("chr19", 41877233),
    "BCL2": ("chr18", 63123346),
    "IRF4": ("chr6", 391739),
    "VPREB1": ("chr22", 22547505),
    "IGLL1": ("chr22", 22355167),
    # T cell genes
    "CD3E": ("chr11", 118344316),
    "IL7R": ("chr5", 35876274),
    "TCF7": ("chr5", 134117803),
    "FOXP3": ("chrX", 49250436),
    "RUNX1": ("chr21", 34787801),
    "LEF1": ("chr4", 108047038),
    "BCL11B": ("chr14", 99169874),
    # Myeloid genes
    "CD14": ("chr5", 140631728),
    "CD68": ("chr17", 7593628),
    "CSF1R": ("chr5", 150053291),
    "MPO": ("chr17", 58269861),
    "LYZ": ("chr12", 69348350),
    # TFs (for reference)
    "PAX5": ("chr9", 36833275),
    "EBF1": ("chr5", 158523083),
    "GATA3": ("chr10", 8045378),
    "SPI1": ("chr11", 47380115),
    "CEBPA": ("chr19", 33299934),
    "CD34": ("chr1", 207928564),
}


def toy(
    n_cells: int = 60,
    n_tfs: int = 3,
    n_targets_per_tf: int = 5,
    n_peaks_per_target: int = 1,
    celltypes: list[str] | None = None,
    seed: int = 42,
    verbose: bool = False,
    max_expr: float = 12.0,
) -> tuple[mu.MuData, pd.DataFrame]:
    """
    Generate synthetic MuData and GRN for testing and demonstration.

    Creates biologically structured test data with RNA and ATAC modalities,
    along with a gene regulatory network (GRN) DataFrame. CREs are placed
    within +/- 100kb of target gene TSS for biological realism. Some genes
    (e.g., IRF4, RUNX1) act as "hub" genes regulated by 4+ TFs across celltypes.

    Parameters
    ----------
    n_cells : int
        Number of cells to generate. Should be divisible by number of celltypes
        for equal distribution. Default is 60.
    n_tfs : int
        Number of transcription factors. Must be between 1 and 6. TFs are paired
        by celltype (PAX5/EBF1 for B cells, GATA3/TCF7 for T cells, SPI1/CEBPA
        for monocytes). TFs of the same celltype can co-bind to shared CREs,
        while also having their own unique CREs. Default is 3.
    n_targets_per_tf : int
        Number of target genes per TF. Must be between 1 and 5. Default is 5.
    n_peaks_per_target : int
        Number of CREs (peaks) per target gene. Each CRE is placed within
        +/- 100kb of the gene's TSS. Default is 1.
    celltypes : list[str] | None
        Custom celltype names. If None, uses default celltypes based on n_tfs:
        ["B cell", "T cell", "Monocyte"]. Default is None.
    seed : int
        Random seed for reproducibility. Default is 42.
    verbose : bool
        Whether to log progress messages. Default is False.
    max_expr : float
        Maximum value for expression/accessibility values. Values are scaled
        to resemble log-normalized counts (0 to max_expr). Default is 12.0.

    Returns
    -------
    tuple[mu.MuData, pd.DataFrame]
        A tuple containing:
        - MuData with 'rna' and 'atac' modalities
        - GRN DataFrame with columns: source, target, cre, score
    """
    np.random.seed(seed)

    # Validate parameters
    if n_tfs < 1 or n_tfs > 6:
        raise ValueError("n_tfs must be between 1 and 6")
    if n_targets_per_tf < 1 or n_targets_per_tf > 5:
        raise ValueError("n_targets_per_tf must be between 1 and 5")
    if n_peaks_per_target < 1:
        raise ValueError("n_peaks_per_target must be at least 1")

    # Select TFs and their properties
    tf_names = list(_TF_CELLTYPE.keys())[:n_tfs]
    tf_celltypes = [_TF_CELLTYPE[tf] for tf in tf_names]

    # Set up celltypes (use unique celltypes from selected TFs)
    if celltypes is None:
        # Preserve order while getting unique celltypes
        seen_ct = set()
        celltypes = []
        for ct in tf_celltypes:
            if ct not in seen_ct:
                celltypes.append(ct)
                seen_ct.add(ct)
    n_celltypes = len(celltypes)

    if n_cells % n_celltypes != 0:
        _log(
            f"n_cells={n_cells} not divisible by {n_celltypes} celltypes, distribution will be uneven",
            level="warning",
            verbose=verbose,
        )

    _log(f"Generating toy data with {n_cells} cells, {n_tfs} TFs", level="info", verbose=verbose)

    # Build gene list: TFs + unique targets (preserve order)
    target_genes = []
    seen = set()
    for tf in tf_names:
        for target in _TF_TARGETS[tf][:n_targets_per_tf]:
            if target not in seen:
                target_genes.append(target)
                seen.add(target)
    all_genes = tf_names + target_genes

    # Build GRN DataFrame with TSS-proximal CREs
    # TFs get unique CREs, but TFs of the same celltype can share CREs (co-binding)
    grn_records = []
    peak_set = set()  # Track unique peaks
    peak_names = []
    tf_target_peaks = {}  # Map (tf, target) -> list of peak names for ATAC patterns
    cre_to_celltype = {}  # Track which celltype each CRE is accessible in

    # First pass: generate shared CREs for co-binding TFs of the same celltype
    shared_cre_map = {}  # Map (target, celltype) -> CRE name
    for target, celltype_offsets in _SHARED_CRES.items():
        if target not in _GENE_TSS:
            continue
        chrom, tss = _GENE_TSS[target]
        for celltype, offset in celltype_offsets.items():
            start = max(1, tss + offset)
            end = start + 500
            cre = f"{chrom}-{start}-{end}"
            shared_cre_map[(target, celltype)] = cre
            if cre not in peak_set:
                peak_set.add(cre)
                peak_names.append(cre)
                cre_to_celltype[cre] = celltype

    # Second pass: generate unique CREs per TF and add shared CREs where applicable
    for tf_idx, tf in enumerate(tf_names):
        tf_celltype = _TF_CELLTYPE[tf]
        targets = _TF_TARGETS[tf][:n_targets_per_tf]

        for target in targets:
            # Get TSS info for this target
            if target in _GENE_TSS:
                chrom, tss = _GENE_TSS[target]
            else:
                # Fallback for unknown genes
                chrom, tss = "chr1", 10000000

            tf_target_peaks[(tf, target)] = []

            # Add shared CRE if this TF's celltype has one for this target
            shared_key = (target, tf_celltype)
            if shared_key in shared_cre_map:
                shared_cre = shared_cre_map[shared_key]
                tf_target_peaks[(tf, target)].append(shared_cre)
                score = np.random.uniform(0.5, 0.95)  # Shared CREs tend to have good scores
                grn_records.append(
                    {
                        "source": tf,
                        "target": target,
                        "cre": shared_cre,
                        "score": round(score, 2),
                    }
                )

            # Add unique CREs for this TF
            for peak_idx in range(n_peaks_per_target):
                # Place CRE within +/- 100kb of TSS
                # Use tf_idx and peak_idx to ensure different TFs get different CREs
                offset = np.random.randint(-100000, 100000)
                # Add TF-specific shift to ensure different TFs get different CREs
                offset += tf_idx * 5000 + peak_idx * 1000
                start = max(1, tss + offset)
                end = start + 500

                cre = f"{chrom}-{start}-{end}"

                # Track peaks for this TF-target pair
                tf_target_peaks[(tf, target)].append(cre)

                # Only add to peak list if not already present
                if cre not in peak_set:
                    peak_set.add(cre)
                    peak_names.append(cre)
                    cre_to_celltype[cre] = tf_celltype

                score = np.random.uniform(0.4, 0.95)
                grn_records.append(
                    {
                        "source": tf,
                        "target": target,
                        "cre": cre,
                        "score": round(score, 2),
                    }
                )

    grn = pd.DataFrame(grn_records)

    # Create cell type assignments
    cells_per_type = n_cells // n_celltypes
    remainder = n_cells % n_celltypes
    celltype_list = []
    for i, ct in enumerate(celltypes):
        count = cells_per_type + (1 if i < remainder else 0)
        celltype_list.extend([ct] * count)

    # Create RNA expression matrix with biological structure
    # Values resemble log-normalized counts (sparse, with most values low)
    n_genes = len(all_genes)
    # Use exponential distribution for sparse-like pattern, then scale
    X_rna = np.random.exponential(scale=1.0, size=(n_cells, n_genes)).astype(np.float32)
    X_rna = np.clip(X_rna, 0, max_expr * 0.4)  # Base expression is low

    # Add TF-specific expression patterns
    for i, tf in enumerate(tf_names):
        tf_idx = all_genes.index(tf)
        ct = tf_celltypes[i]
        # Find cells of this celltype
        cell_mask = np.array([c == ct for c in celltype_list])
        # Add higher expression for TF in its celltype
        X_rna[cell_mask, tf_idx] += np.random.uniform(3.0, 6.0, size=cell_mask.sum())

        # Also upregulate target genes in the corresponding celltype
        targets = _TF_TARGETS[tf][:n_targets_per_tf]
        for target in targets:
            if target in all_genes:
                target_idx = all_genes.index(target)
                X_rna[cell_mask, target_idx] += np.random.uniform(2.0, 5.0, size=cell_mask.sum())

    # Clip to max expression value
    X_rna = np.clip(X_rna, 0, max_expr)

    # Create RNA AnnData
    rna = ad.AnnData(X=X_rna)
    rna.var_names = all_genes
    rna.obs_names = [f"Cell{i}" for i in range(n_cells)]

    # Create ATAC accessibility matrix with biological structure
    n_peaks = len(peak_names)
    # Use exponential distribution for sparse-like pattern
    X_atac = np.random.exponential(scale=0.8, size=(n_cells, n_peaks)).astype(np.float32)
    X_atac = np.clip(X_atac, 0, max_expr * 0.3)  # Base accessibility is low

    # Add TF-specific accessibility patterns
    # CREs associated with a TF should be more accessible in cells where that TF is active
    peak_name_to_idx = {p: i for i, p in enumerate(peak_names)}
    for i, tf in enumerate(tf_names):
        ct = tf_celltypes[i]
        cell_mask = np.array([c == ct for c in celltype_list])
        targets = _TF_TARGETS[tf][:n_targets_per_tf]

        for target in targets:
            # Get CREs associated with this TF-target pair
            cre_list = tf_target_peaks.get((tf, target), [])
            for cre in cre_list:
                if cre in peak_name_to_idx:
                    peak_idx = peak_name_to_idx[cre]
                    # Increase accessibility in cells of the corresponding celltype
                    X_atac[cell_mask, peak_idx] += np.random.uniform(2.0, 5.0, size=cell_mask.sum())

    # Clip to max value
    X_atac = np.clip(X_atac, 0, max_expr)

    # Create ATAC AnnData
    atac = ad.AnnData(X=X_atac)
    atac.var_names = peak_names
    atac.obs_names = rna.obs_names.copy()

    # Create MuData
    mdata = mu.MuData({"rna": rna, "atac": atac})
    mdata.obs["celltype"] = celltype_list

    _log(
        f"Created MuData with {n_cells} cells, {n_genes} genes, {n_peaks} peaks",
        level="info",
        verbose=verbose,
    )
    _log(f"Created GRN with {len(grn)} edges", level="info", verbose=verbose)

    return mdata, grn
