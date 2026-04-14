import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from decoupler._Plotter import Plotter


_CNAME_DICT = {
    "n_sources": "TFs",
    "n_cres": "CREs",
    "n_targets": "Genes",
    "n_edges": "Edges",
    "mean_regulon_size": "Regulon",
}


def stats(
    df: pd.DataFrame,
    order: list | None = None,
    palette: str | None = None,
    **kwargs,
) -> plt.Figure | None:
    """
    Plot GRN statistics as horizontal bar plots.

    Parameters
    ----------
    df
        Output from tl.stats() with columns: name, n_sources, n_cres,
        n_targets, n_edges, mean_regulon_size.
    order
        Order of GRN names on the y-axis. If None, uses the order in df.
    palette
        Seaborn color palette. If None, uses the default palette.
    **kwargs
        Additional arguments passed to ``decoupler.Plotter`` (e.g. ``figsize``,
        ``dpi``, ``return_fig``, ``save``).

    Returns
    -------
    plt.Figure or None
        Figure if ``return_fig=True``.
    """
    cols = list(_CNAME_DICT.keys())

    kwargs.setdefault("figsize", (9, 2))
    kwargs["ax"] = None
    bp = Plotter(**kwargs)
    bp.fig.delaxes(bp.ax)
    plt.close(bp.fig)

    fig, axes = plt.subplots(
        1, len(cols), sharey=True, figsize=bp.figsize, tight_layout=True
    )

    for ax, col in zip(axes, cols):
        sns.barplot(
            data=df,
            x=col,
            y="name",
            hue="name",
            order=order,
            orient="h",
            palette=palette,
            ax=ax,
            legend=False,
        )
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title(_CNAME_DICT[col])

    bp.fig = fig
    bp.fig.set_figwidth(bp.figsize[0])
    bp.fig.set_figheight(bp.figsize[1])
    bp.fig.set_dpi(bp.dpi)
    return bp._return()
