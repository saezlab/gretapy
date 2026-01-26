import decoupler as dc
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyranges as pr

import gretapy as gt


def _norm_score(x: np.ndarray, axis: int | None = None) -> np.ndarray:
    """Normalize values to 0-1 range."""
    if np.all(x == x[0]):
        return x
    return (x - x.min(axis=axis)) / (x.max(axis=axis) - x.min(axis=axis))


def _get_tss_window(gannot: pr.PyRanges, gene: str, w_size: int) -> pr.PyRanges:
    """Get genomic window around TSS of a gene."""
    gene_gr = gannot[gannot.df["Name"] == gene]
    strand = gene_gr.df["Strand"].values[0]
    if strand == "+":
        tss = gene_gr.df["Start"].values[0]
    else:
        tss = gene_gr.df["End"].values[0]
    tss_window = pr.from_dict(
        {
            "Chromosome": gene_gr.Chromosome,
            "Start": tss - w_size,
            "End": tss + w_size,
        }
    )
    return tss_window


def _get_gannot_data(gannot: pr.PyRanges, target: str, w_size: int, atac_var_names: pd.Index):
    """Extract genomic annotation data for plotting."""
    target_gr = gannot[gannot.df["Name"] == target]
    wind = _get_tss_window(gannot, target, w_size)
    x_min = wind.df["Start"].values[0]
    x_max = wind.df["End"].values[0]
    gs_gr = gannot.join(wind).sort()

    chromosome = gs_gr.df["Chromosome"].values[0]
    strand = target_gr.df["Strand"].values[0]
    if strand == "+":
        tss = target_gr.df["Start"].values[0]
    else:
        tss = target_gr.df["End"].values[0]

    atac_vars = atac_var_names[atac_var_names.str.startswith(chromosome)]
    cres_gr = []
    for cre in atac_vars:
        cre_chr, cre_start, cre_end = cre.split("-")
        if cre_chr == chromosome:
            cres_gr.append([cre_chr, cre_start, cre_end, cre])
    cres_gr = pr.PyRanges(pd.DataFrame(cres_gr, columns=["Chromosome", "Start", "End", "Name"]))

    return x_min, x_max, gs_gr, tss, chromosome, strand, cres_gr


def _plot_links(
    links: pd.DataFrame,
    tf: str,
    tss: int,
    strand: str,
    palette: dict,
    ax: plt.Axes,
    show_legend: bool = True,
) -> None:
    """Plot arcs connecting CREs to TSS for a given TF."""
    grns = links["grn"].sort_values().unique()
    links = links.copy()
    links = links[links["tf"] == tf]
    is_empty = []
    for grn in grns:
        link = links[links["grn"] == grn]
        if link.shape[0] == 0:
            is_empty.append(True)
        else:
            is_empty.append(False)
        for _, row in link.iterrows():
            cre, score = row["cre"], row["score"]
            if score > 0.00:
                _, cre_start, cre_end = cre.split("-")
                cre_start, cre_end = float(cre_start), float(cre_end)
                if strand == "+":
                    if cre_end < tss:
                        arc_start = cre_end
                        arc_width = tss - cre_end
                    elif cre_start > tss:
                        arc_start = tss
                        arc_width = cre_start - tss
                    elif cre_start < tss < cre_end:
                        arc_start = cre_start
                        arc_width = tss - cre_start
                else:
                    if cre_end < tss:
                        arc_start = cre_end
                        arc_width = tss - cre_end
                    elif cre_start > tss:
                        arc_start = tss
                        arc_width = cre_start - tss
                    elif cre_start < tss < cre_end:
                        arc_start = tss
                        arc_width = cre_end - tss
                center_x = (arc_start + arc_start + arc_width) // 2
                arc = patches.Arc(
                    (center_x, 0), arc_width, score * 2, angle=0, theta1=0, theta2=180, color=palette[grn], linewidth=1
                )
                ax.add_patch(arc)
    ax.set_ylim(0, 1.05)
    if show_legend:
        handles = [
            plt.Line2D([0], [0], marker="o", color="w", label=grn, markerfacecolor=palette[grn], markersize=10)
            for grn in grns
        ]
        handles = [h for i, h in enumerate(handles) if not is_empty[i]]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1, 0.5), title="", frameon=False)
    ax.set_ylabel(tf)
    yticks = np.arange(0.25, 1, 0.25)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)


def _plot_gannot(gs_gr: pr.PyRanges, x_min: int, x_max: int, target: str, ax: plt.Axes) -> None:
    """Plot genomic annotations (genes as arrows)."""

    def plot_gene(row, i, ax):
        row = row[["Chromosome", "Start", "End", "Strand", "Name"]]
        chrm, strt, end, strnd, name = row
        ax.plot([strt, end], [i, i], lw=1, zorder=0, color="#440a82")
        if strnd == "+":
            s_marker = ">"
        else:
            s_marker = "<"
        ax.scatter([strt, end], [i, i], marker=s_marker, color="#440a82", s=10)
        ax.text(x_max, i, name, ha="left", va="center")

    # Sort genes by unstranded start site (TSS) descending, with target gene at top
    df = gs_gr.df.copy()
    df["_tss"] = np.where(df["Strand"] == "+", df["Start"], df["End"])
    df = df.sort_values("_tss", ascending=False).reset_index(drop=True)
    # Move target gene to the end so it appears at the top of the plot
    target_row = df[df["Name"] == target]
    other_rows = df[df["Name"] != target]
    df = pd.concat([other_rows, target_row]).reset_index(drop=True)

    for i, row in df.iterrows():
        plot_gene(row, i, ax)

    ax.set_xlim(x_min, x_max)
    ax.tick_params(axis="y", labelleft=False, length=0)
    w_size = x_max - x_min
    xticks = np.arange(x_min, x_max + 1, w_size // 2)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks)
    pad = gs_gr.df.shape[0] * 0.1
    ax.set_ylim(0 - pad, gs_gr.df.shape[0] + pad - 1)


def _plot_omic(
    adata,
    feat_gr: pr.PyRanges,
    x_min: int,
    x_max: int,
    cmap,
    mode: str,
    ax: plt.Axes,
) -> None:
    """Plot omics data (expression heatmap or accessibility peaks)."""
    for y, ctype in enumerate(adata.obs_names):
        if mode == "heatmap":
            y += 0.5
            ax.text(x_max, y, ctype, ha="left", va="center")
            for s, e, feat_name in zip(feat_gr.df["Start"], feat_gr.df["End"], feat_gr.df["Name"], strict=True):
                if feat_name in adata.var_names:
                    val = adata[ctype, feat_name].X.ravel()[0]
                    rect = patches.Rectangle((s, y - 0.5), width=e - s, height=0.95, color=cmap(val))
                    ax.add_patch(rect)
        elif mode == "peaks":
            x_coord = [x_min]
            y_coord = [y]
            for s, e, feat_name in zip(feat_gr.df["Start"], feat_gr.df["End"], feat_gr.df["Name"], strict=True):
                if feat_name in adata.var_names:
                    if s >= x_min and e <= x_max:
                        val = adata[ctype, feat_name].X.ravel()[0]
                        x_coord.extend([s, e])
                        y_coord.extend([val + y, y])
            x_coord.append(x_max)
            y_coord.append(y)
            y_coord = (_norm_score(np.array(y_coord)) * 0.75 + y) + 0.05
            ax.step(x_coord, y_coord, where="post", color="gold")

    ax.grid(axis="y", lw=0.25)
    yticks = np.arange(0, adata.shape[0] + 1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ax.tick_params(axis="y", labelleft=False, length=0)
    ax.set_ylim(0, y + 1)


def links(
    mdata,
    grn: pd.DataFrame | dict[str, pd.DataFrame],
    target: str,
    tfs: list[str],
    gannot: pr.PyRanges | str = "hg38",
    w_size: int = 100_000,
    sample_col: str = "celltype",
    agg_mode: str = "mean",
    palette: dict[str, str] | None = None,
    expr_cmap: None | str | list[str] = None,
    figsize: tuple[float, float] | None = None,
    dpi: int = 125,
    return_fig: bool = False,
) -> plt.Figure | None:
    """
    Plot CRE-to-gene links for TFs in a genomic region.

    Creates a multi-panel figure showing:
    - Top panels: Arcs connecting CREs to the target gene TSS for each TF
    - Middle panel: Gene expression heatmap and chromatin accessibility peaks
    - Bottom panel: Genomic annotations (genes as arrows)

    Parameters
    ----------
    mdata
        MuData object with 'rna' and 'atac' modalities. Can be either single-cell
        data (which will be automatically aggregated by cell type) or pre-computed
        pseudobulk data. The 'rna' modality should have genes as variables, and
        'atac' should have CREs in format "chr-start-end".
    grn
        GRN predictions. Either a DataFrame with columns 'source' (TF), 'target' (gene),
        'cre', and 'score', or a dict mapping GRN names to DataFrames.
    target
        Target gene name to plot.
    tfs
        List of TF names to plot (one panel per TF).
    gannot
        Gene annotations. If a string (e.g., 'hg38'), reads the Promoters database
        for that organism using `gretapy.ds.read_db`. If a PyRanges object, uses it
        directly. Default is 'hg38'.
    w_size
        Window size (bp) around target gene TSS. Default is 250000.
    sample_col
        Column in `mdata.obs` containing cell type labels for aggregation.
        If the column exists, data will be aggregated by cell type using
        `decoupler.get_pseudobulk`. Set to None to skip aggregation.
        Default is 'celltype'.
    agg_mode
        Aggregation mode for pseudobulk computation. One of 'mean', 'sum',
        or 'median'. Default is 'mean'.
    palette
        Color palette mapping GRN names to colors. If None, uses default colors.
    expr_cmap
        Colormap for expression heatmap. Default is white to purple.
    figsize
        Figure size (width, height). If None, auto-calculated.
    dpi
        Figure DPI. Default is 150.
    return_fig
        Whether to return the figure. Default is False.

    Returns
    -------
    plt.Figure or None
        Figure if return_fig is True.

    Examples
    --------
    >>> import gretapy as gt
    >>> # Load data and generate GRN
    >>> mdata, grn = gt.ds.toy()
    >>> # Plot links - data is automatically aggregated by cell type
    >>> gt.pl.links(mdata, grn, target="CD19", tfs=["PAX5", "EBF1"])
    >>> # Plot links for multiple GRNs
    >>> grns = {"granie": grn1, "celloracle": grn2}
    >>> gt.pl.links(mdata, grns, target="AREG", tfs=["JUND", "FOSL2"])
    """
    # Extract RNA and ATAC modalities from MuData
    rna = mdata.mod["rna"].copy()
    atac = mdata.mod["atac"].copy()

    # Aggregate by cell type if sample_col is provided and exists
    if sample_col is not None and sample_col in mdata.obs.columns:
        rna.obs = mdata.obs[[sample_col]].copy()
        atac.obs = mdata.obs[[sample_col]].copy()
        rna = dc.pp.pseudobulk(adata=rna, sample_col=sample_col, groups_col=None, mode=agg_mode)
        atac = dc.pp.pseudobulk(adata=atac, sample_col=sample_col, groups_col=None, mode=agg_mode)

    # Load gene annotations if string provided
    if isinstance(gannot, str):
        gannot = gt.show_genome_annotation(organism=gannot)

    # Normalize GRN input to dict
    if isinstance(grn, pd.DataFrame):
        grn_dict = {"grn": grn}
    else:
        grn_dict = grn

    # Build links DataFrame
    links_list = []
    for grn_name, grn_df in grn_dict.items():
        df = grn_df.copy()
        df = df.rename(columns={"source": "tf", "target": "gene", "score": "raw_score"})
        df["grn"] = grn_name
        df = df[df["gene"] == target]
        df["score"] = df["raw_score"].abs().rank(method="average", pct=True)
        links_list.append(df)
    links_df = pd.concat(links_list)

    # Set default palette
    if palette is None:
        default_colors = plt.cm.tab10.colors
        grn_names = sorted(grn_dict.keys())
        palette = {name: default_colors[i % len(default_colors)] for i, name in enumerate(grn_names)}

    # Determine if legend should be shown (only for multiple GRNs)
    show_legend = len(grn_dict) > 1

    # Extract genomic annotation data (before figure creation to know counts)
    x_min, x_max, gs_gr, tss, chromosome, strand, cres_gr = _get_gannot_data(gannot, target, w_size, atac.var_names)

    # Calculate dynamic heights based on content
    n_tfs = len(tfs)
    n_celltypes = rna.shape[0]
    n_genes = gs_gr.df.shape[0]

    # Height per cell type row (minimum 0.4 inches per row)
    omic_height = max(n_celltypes * 0.4, 1.5)
    # Height per gene row (minimum 0.3 inches per row)
    gannot_height = max(n_genes * 0.3, 1.0)

    # Set up figure with dynamic sizing
    if figsize is None:
        fig_height = n_tfs * 1.0 + omic_height + gannot_height
        figsize = (3, fig_height)
    height_ratios = [1] * n_tfs + [omic_height, gannot_height]
    fig, axes = plt.subplots(2 + n_tfs, 1, figsize=figsize, dpi=dpi, sharex=True, height_ratios=height_ratios)
    axes = axes.ravel()

    # Expression colormap
    if expr_cmap is None:
        expr_cmap = ["white", "#3f007d"]
    if isinstance(expr_cmap, list):
        cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", expr_cmap)
    else:
        cmap = plt.get_cmap(expr_cmap)

    # Plot TF links
    for i, tf in enumerate(tfs):
        ax = axes[i]
        _plot_links(links_df, tf, tss, strand, palette, ax, show_legend=show_legend)

    # Plot omics data
    ax = axes[-2]
    _plot_omic(rna, gs_gr, x_min, x_max, cmap, "heatmap", ax)
    _plot_omic(atac, cres_gr, x_min, x_max, cmap, "peaks", ax)

    # Plot gene annotations
    ax = axes[-1]
    _plot_gannot(gs_gr, x_min, x_max, target, ax)
    ax.set_xlabel(chromosome)

    fig.subplots_adjust(wspace=0, hspace=0.0)

    if return_fig:
        plt.close()
        return fig
