import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from decoupler._Plotter import Plotter
from matplotlib.gridspec import GridSpec
from scipy.stats import rankdata

from gretapy.tl._ranking import CLASS_ORDER, _class_mean
from gretapy.tl._ranking import ranking as _tl_ranking

# Default color palette
PALETTE = {
    "CellOracle": "#1e76b4",
    "CREMA": "#aec6e7",
    "Dictys": "#ff7f0e",
    "DIRECT-NET": "#ffba78",
    "FigR": "#2ca02c",
    "GRaNIE": "#97de8a",
    "HuMMuS": "#d52727",
    "Inferelator3.0": "#ff9796",
    "Pando": "#9467bc",
    "scDoRI": "#c4afd4",
    "SCENIC+": "#8c554a",
    "scMTNI": "#c49c94",
    "CollecTRI": "#e276c1",
    "DoRothEA": "#f6b6d2",
    "GRNBoost2": "#bcbd22",
    "SCENIC": "#dada8d",
    "Pearson": "#16becf",
    "Spearman": "#9edae4",
    "scGPT": "#7e7e7e",
    "Random": "black",
}


def _compute_db_aggregation(df):
    """Compute mean F0.1 per method per (db, task) and the class mapping."""
    mean = df.groupby(["name", "db", "task"])["f01"].mean().reset_index()
    pivot = mean.pivot(index="name", columns=["db", "task"], values="f01")
    hierarchy = df[["db", "task", "class"]].drop_duplicates()
    return pivot, hierarchy


def _compute_ranks(values):
    """Rank each column: rank 1 = best (highest value). NaN stays NaN."""
    ranks = pd.DataFrame(index=values.index, columns=values.columns)
    for col in values.columns:
        v = values[col].values
        mask = ~np.isnan(v)
        if mask.sum() > 0:
            r = np.full(len(v), np.nan)
            r[mask] = mask.sum() + 1 - rankdata(v[mask], method="average")
            ranks[col] = r
    return ranks


def _build_task_list(hierarchy, col_order):
    """Build list of (class, task, [(db, task), ...]) grouped by class and task."""
    task_list = []
    for cls in CLASS_ORDER:
        cls_rows = hierarchy[hierarchy["class"] == cls]
        if cls_rows.empty:
            continue
        cls_cols = [(d, t) for d, t in col_order if any((cls_rows["db"] == d) & (cls_rows["task"] == t))]
        current_task = None
        task_cols = []
        for d, t in cls_cols:
            if t != current_task:
                if current_task is not None:
                    task_list.append((cls, current_task, task_cols))
                current_task = t
                task_cols = []
            task_cols.append((d, t))
        if current_task is not None:
            task_list.append((cls, current_task, task_cols))
    return task_list


def _rank_text(rank_val):
    """Format rank value as string."""
    if rank_val == int(rank_val):
        return f"{int(rank_val)}"
    return f"{rank_val:.1f}"


def _draw_barplot(ax, method_order, overall_mean, palette, n_methods):
    """Draw the horizontal barplot on the given axes."""
    y = np.arange(n_methods)
    colors = [palette.get(m, "gray") for m in method_order]
    vals = overall_mean.values
    ax.barh(y, vals, color=colors, height=0.95)
    ax.set_ylim(-0.5, n_methods - 0.5)
    ax.invert_yaxis()
    ax.set_xlabel(r"Mean F$\mathrm{_{0.1}}$", fontsize=10)
    ax.set_yticks([])
    ax.set_xlim(0, vals.max() * 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.axvline(x=np.mean(vals), ls="--", lw=1, c="black")


def _draw_method_names(ax, method_order, n_methods):
    """Draw method name labels on the given axes."""
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, n_methods - 0.5)
    ax.invert_yaxis()
    ax.axis("off")
    for i, m in enumerate(method_order):
        ax.text(1, i, m, ha="right", va="center", fontsize=10, fontweight="bold")


def _ranking_class(df, overall_mean, class_mean, method_order, n_methods, palette, figsize):
    """Create the class-level ranking figure."""
    class_cols = [c for c in CLASS_ORDER if c in class_mean.columns]
    class_mean = class_mean.loc[method_order, class_cols]
    class_ranks = _compute_ranks(class_mean)

    if figsize is None:
        figsize = (6, max(4, n_methods * 0.35))

    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = GridSpec(
        2,
        4,
        figure=fig,
        width_ratios=[2, 2, 0.05, len(class_cols)],
        height_ratios=[20, 0.4],
        wspace=0.02,
        hspace=0.02,
    )

    _draw_method_names(fig.add_subplot(gs[0, 0]), method_order, n_methods)
    _draw_barplot(fig.add_subplot(gs[0, 1]), method_order, overall_mean, palette, n_methods)

    # Heatmap
    ax_hmap = fig.add_subplot(gs[0, 3])
    data = class_mean.values
    ranks = class_ranks.values
    vmin, vmax = 0, np.nanmax(data)
    threshold = vmax / 2
    ax_hmap.imshow(data, aspect="auto", cmap="viridis", vmin=vmin, vmax=vmax)
    for i in range(n_methods):
        for j in range(len(class_cols)):
            rv, cv = ranks[i, j], data[i, j]
            if pd.isna(rv) or pd.isna(cv):
                txt, color = "-", "gray"
            else:
                color = "white" if cv < threshold else "black"
                txt = _rank_text(rv)
            ax_hmap.text(j, i, txt, ha="center", va="center", fontsize=10, color=color)
    ax_hmap.set_xticks(range(len(class_cols)))
    ax_hmap.set_xticklabels(class_cols, fontsize=10, rotation=90)
    ax_hmap.xaxis.set_ticks_position("top")
    ax_hmap.set_yticks([])

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=vmin, vmax=vmax))
    cbar = plt.colorbar(sm, cax=fig.add_subplot(gs[1, 3]), orientation="horizontal")
    cbar.set_label(r"Mean F$\mathrm{_{0.1}}$", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    return fig


def _ranking_task(df, method_order, n_methods, figsize):
    """Create the task-level ranking figure with hierarchical headers."""
    db_mean, hierarchy = _compute_db_aggregation(df)

    # Order columns: CLASS_ORDER, tasks alphabetically, dbs alphabetically
    available = set(db_mean.columns.tolist())
    col_order = []
    for cls in CLASS_ORDER:
        cls_rows = hierarchy[hierarchy["class"] == cls]
        for task in sorted(cls_rows["task"].unique()):
            for db_name in sorted(cls_rows[cls_rows["task"] == task]["db"].unique()):
                if (db_name, task) in available:
                    col_order.append((db_name, task))

    db_mean = db_mean.loc[method_order, col_order]
    db_ranks = _compute_ranks(db_mean)

    # Build task groupings
    task_list = _build_task_list(hierarchy, col_order)

    # Compute wrapped task titles and effective widths
    chars_per_col = 2.5
    task_widths = []
    task_wrapped = {}
    max_lines = 1
    for _cls, task, cols in task_list:
        wrapper = textwrap.TextWrapper(width=12, break_long_words=False, break_on_hyphens=False)
        lines = wrapper.wrap(task)
        task_wrapped[task] = "\n".join(lines)
        max_lines = max(max_lines, len(lines))
        title_width = max(len(l) for l in lines) / chars_per_col if lines else 0
        task_widths.append(max(len(cols), title_width))

    # Build GridSpec width ratios: [names, task1, task2, gap_between_classes, ...]
    gap_width = 0.3
    width_ratios = [1.5]  # names
    gs_col_mapping = []
    current_class = None
    for i, (cls, _task, _cols) in enumerate(task_list):
        if current_class is not None and cls != current_class:
            width_ratios.append(gap_width)  # gap between classes
        width_ratios.append(task_widths[i])
        gs_col_mapping.append(len(width_ratios) - 1)
        current_class = cls

    n_gs_cols = len(width_ratios)
    total_data_cols = sum(len(cols) for _, _, cols in task_list)
    task_header_height = 0.4 * max_lines

    if figsize is None:
        figsize = (max(14, total_data_cols * 0.5 + 4), max(6, n_methods * 0.35))

    fig = plt.figure(figsize=figsize)

    # Rows: class_header, task_header, heatmap, db_labels, colorbar
    gs = GridSpec(
        5,
        n_gs_cols,
        figure=fig,
        height_ratios=[0.4, task_header_height, 10, 1.5, 0.3],
        width_ratios=width_ratios,
        wspace=0.02,
        hspace=0.02,
    )

    # Method names (span the heatmap row)
    _draw_method_names(fig.add_subplot(gs[2, 0]), method_order, n_methods)

    # Shared colormap
    vmin, vmax = 0, np.nanmax(db_mean.values)
    threshold = vmax / 2
    cmap = plt.cm.viridis

    # Track class spans for class headers
    class_spans = {}
    for i, (cls, _task, _cols) in enumerate(task_list):
        gs_col = gs_col_mapping[i]
        if cls not in class_spans:
            class_spans[cls] = [gs_col, gs_col]
        else:
            class_spans[cls][1] = gs_col

    # Class headers (spanning)
    for cls, (start_col, end_col) in class_spans.items():
        ax = fig.add_subplot(gs[0, start_col : end_col + 1])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.add_patch(
            plt.Rectangle(
                (0, 0),
                1,
                1,
                facecolor="white",
                edgecolor="black",
                linewidth=1,
            )
        )
        ax.text(0.5, 0.5, cls, ha="center", va="center", fontsize=10, fontweight="bold")

    im = None
    # Render each task block
    for i, (_cls, task, cols) in enumerate(task_list):
        gs_col = gs_col_mapping[i]
        n_cols = len(cols)

        # Task header
        ax_task = fig.add_subplot(gs[1, gs_col])
        ax_task.set_xlim(0, 1)
        ax_task.set_ylim(0, 1)
        ax_task.axis("off")
        ax_task.add_patch(
            plt.Rectangle(
                (0, 0),
                1,
                1,
                facecolor="white",
                edgecolor="gray",
                linewidth=0.5,
            )
        )
        ax_task.text(0.5, 0.5, task_wrapped[task], ha="center", va="center", fontsize=9)

        # Heatmap for this task
        ax_hm = fig.add_subplot(gs[2, gs_col])
        data = db_mean[cols].values
        im = ax_hm.imshow(data, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        for ri in range(n_methods):
            for cj in range(n_cols):
                rv = db_ranks.loc[method_order[ri], cols[cj]]
                cv = data[ri, cj]
                if pd.isna(rv) or pd.isna(cv):
                    txt, color = "-", "gray"
                else:
                    color = "white" if cv < threshold else "black"
                    txt = _rank_text(rv)
                ax_hm.text(cj, ri, txt, ha="center", va="center", fontsize=8, color=color)
        ax_hm.set_xticks([])
        ax_hm.set_yticks([])

        # Database labels at bottom
        ax_db = fig.add_subplot(gs[3, gs_col])
        ax_db.set_xlim(-0.5, n_cols - 0.5)
        ax_db.set_ylim(0, 1)
        ax_db.axis("off")
        for j, (db_name, _) in enumerate(cols):
            ax_db.text(j, 1, db_name, ha="center", va="top", fontsize=8, rotation=90)

    # Colorbar (centered at bottom)
    cbar_width = 0.15
    cbar_ax = fig.add_axes([0.5 - cbar_width / 2, 0.02, cbar_width, 0.015])
    cbar = plt.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label(r"Mean F$\mathrm{_{0.1}}$", fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    return fig


def ranking(
    df,
    level="class",
    palette=None,
    metric_weights=None,
    **kwargs,
):
    """
    Plot a ranking figure with a barplot of mean F0.1 and a heatmap of rankings.

    Parameters
    ----------
    df : pandas.DataFrame
        Metrics dataframe with columns: class, task, db, dataset, name, f01.
    level : str
        ``'class'`` for summary heatmap at class level (Predictive, Genomic, etc.),
        ``'task'`` for the fine-grained (db, task) heatmap with hierarchical headers.
    palette : dict or None
        Method name -> color mapping. Uses default palette if None.
    metric_weights : dict or None
        Weight for each metric class when computing the final mean F0.1 that
        determines the row ordering. Keys are class names (``predictive``,
        ``genomic``, ``literature``, ``mechanistic``); all must be non-negative.
        Any metric class missing from the dictionary is treated as having a
        weight of 0 (i.e. excluded from the ranking).
        Defaults to ``dict(predictive=1, genomic=1, literature=1, mechanistic=1)``.
        The overall score is the weighted mean across classes, i.e. the sum of
        ``weight * value`` divided by the total of the weights of the classes
        present for that method.
    **kwargs
        Additional arguments passed to ``decoupler.Plotter`` (e.g. ``figsize``,
        ``dpi``, ``return_fig``, ``save``). ``figsize`` defaults to auto-computed
        based on data dimensions. ``dpi`` defaults to ``100``.

    Returns
    -------
    fig : matplotlib.figure.Figure or None
    """
    if palette is None:
        palette = PALETTE

    # Weighted overall mean (validates metric_weights, sorted descending) drives ordering
    overall_mean = _tl_ranking(df, metric_weights)["mean_f01"]
    method_order = overall_mean.index.tolist()
    class_mean = _class_mean(df).loc[method_order]
    n_methods = len(method_order)

    # Extract user-provided figsize before Plotter (sub-functions auto-calc when None)
    user_figsize = kwargs.get("figsize", None)
    kwargs.setdefault("figsize", (4, 3))
    kwargs.setdefault("dpi", 100)
    kwargs["ax"] = None
    bp = Plotter(**kwargs)
    bp.fig.delaxes(bp.ax)
    plt.close(bp.fig)

    if level == "class":
        fig = _ranking_class(df, overall_mean, class_mean, method_order, n_methods, palette, user_figsize)
    elif level == "task":
        fig = _ranking_task(df, method_order, n_methods, user_figsize)
    else:
        raise ValueError(f"level must be 'class' or 'task', got '{level}'")

    bp.fig = fig
    if user_figsize is not None:
        bp.fig.set_figwidth(bp.figsize[0])
        bp.fig.set_figheight(bp.figsize[1])
    bp.fig.set_dpi(bp.dpi)
    return bp._return()
