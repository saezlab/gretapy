import marsilea as ma
import marsilea.plotter as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
    return_fig: bool = False,
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
        Width of the heatmap. Default is 2.
    height
        Height of the heatmap. Default is 2.
    return_fig
        Whether to return the figure. Default is False.

    Returns
    -------
    plt.Figure or None
        Figure if return_fig is True.
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

    h = ma.Heatmap(mat, cmap=cmap, width=width, height=height, label="Overlap\nCoefficient", vmin=vmin, vmax=vmax)
    h.add_bottom(mp.Labels(mat.columns))
    h.add_left(mp.Labels(mat.index))
    h.add_top(mp.Title(title))
    h.add_legends()
    h.render()

    if return_fig:
        fig = h.figure
        plt.close()
        return fig
