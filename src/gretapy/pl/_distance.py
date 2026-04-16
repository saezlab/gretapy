import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from decoupler._Plotter import Plotter


def cre_to_tss_distance(
    df: pandas.DataFrame,
    order: list | None = None,
    palette: str | None = None,
    thr_distance: int | None = 250_000,
    **kwargs,
) -> plt.Figure | None:
    """
    Plot CRE-to-TSS distance distributions as horizontal boxplots.

    Parameters
    ----------
    df
        DataFrame with columns: grn, cre, target, distance.
    order
        Order of GRN names on the y-axis. If None, uses the order in df.
    palette
        Seaborn color palette. If None, uses the default palette.
    thr_distance
        Distance threshold to draw as a vertical dashed line. If None, no
        line is drawn. Default is 250_000.
    **kwargs
        Additional arguments passed to ``decoupler.Plotter`` (e.g. ``figsize``,
        ``dpi``, ``return_fig``, ``save``).

    Returns
    -------
    plt.Figure or None
        Figure if ``return_fig=True``.
    """
    kwargs.setdefault("figsize", (4, 3))
    kwargs["ax"] = None
    bp = Plotter(**kwargs)
    bp.fig.delaxes(bp.ax)
    plt.close(bp.fig)

    fig, ax = plt.subplots(1, 1, figsize=bp.figsize, tight_layout=True)

    sns.boxplot(
        data=df,
        x="distance",
        y="grn",
        hue="grn",
        order=order,
        palette=palette,
        ax=ax,
    )
    ax.set_ylabel("")
    ax.set_xlabel("CRE Distance To TSS (bp)")

    if thr_distance is not None:
        ax.axvline(thr_distance, ls="--", c="gray")

    bp.fig = fig
    bp.fig.set_figwidth(bp.figsize[0])
    bp.fig.set_figheight(bp.figsize[1])
    bp.fig.set_dpi(bp.dpi)
    return bp._return()
