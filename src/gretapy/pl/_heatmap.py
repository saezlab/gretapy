import marsilea as ma
import marsilea.plotter as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from decoupler._Plotter import Plotter


def _make_sim_mat(df: pd.DataFrame, col: str) -> pd.DataFrame:
    names = list(set(df["grn_a"].tolist() + df["grn_b"].tolist()))
    mat = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    for _, row in df.iterrows():
        mat.loc[row["grn_a"], row["grn_b"]] = row[col]
        mat.loc[row["grn_b"], row["grn_a"]] = row[col]
    return mat


def heatmap(
    df: pd.DataFrame,
    level: str = "edge",
    order: list | None = None,
    title: str | None = None,
    cmap: str = "Purples",
    vmin: float = 0,
    vmax: float = 1,
    width: float = 2,
    height: float = 2,
    **kwargs,
) -> plt.Figure | None:
    """
    Plot overlap coefficient heatmap.

    Parameters
    ----------
    df
        Output from tl.ocoeff with columns: grn_a, grn_b, source, cre, target, edge.
    level
        Which level to plot: "source", "cre", "target", or "edge". Default is "edge".
    order
        Order of GRN names for rows/columns. If None, uses alphabetical order.
    title
        Title for the heatmap. If None, uses the level name.
    cmap
        Colormap name. Default is "Purples".
    vmin
        Minimum value for colormap. Default is 0.
    vmax
        Maximum value for colormap. Default is 1.
    width
        Width of the heatmap in marsilea units. Default is 2.
    height
        Height of the heatmap in marsilea units. Default is 2.
    **kwargs
        Additional arguments passed to ``decoupler.Plotter`` (e.g. ``figsize``,
        ``dpi``, ``return_fig``, ``save``).

    Returns
    -------
    plt.Figure or None
        Figure if ``return_fig=True``.
    """
    if level not in {"source", "cre", "target", "edge"}:
        raise ValueError(f'level must be "source", "cre", "target", or "edge", got {level}')

    mat = _make_sim_mat(df, col=level)

    if order is not None:
        mat = mat.loc[order, order]
    else:
        mat = mat.sort_index(axis=0).sort_index(axis=1)

    if title is None:
        title = level.capitalize()

    kwargs["ax"] = None
    bp = Plotter(**kwargs)
    bp.fig.delaxes(bp.ax)
    plt.close(bp.fig)

    h = ma.Heatmap(mat, cmap=cmap, width=width, height=height, label="Overlap\nCoefficient", vmin=vmin, vmax=vmax)
    h.add_bottom(mp.Labels(mat.columns))
    h.add_left(mp.Labels(mat.index))
    h.add_top(mp.Title(title))
    h.add_legends()
    h.render()

    bp.fig = h.figure
    bp.fig.set_figwidth(bp.figsize[0])
    bp.fig.set_figheight(bp.figsize[1])
    bp.fig.set_dpi(bp.dpi)
    return bp._return()
