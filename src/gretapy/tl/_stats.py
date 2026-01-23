from itertools import combinations

import numpy as np
import pandas as pd
import pyranges as pr


def _map_regions(regions_a: np.ndarray, regions_b: np.ndarray) -> pd.DataFrame:
    """Map overlapping genomic regions between two sets of peaks."""

    def to_df(regions):
        data = [r.split("-") for r in regions]
        return pd.DataFrame(data, columns=["Chromosome", "Start", "End"]).assign(
            Start=lambda x: x.Start.astype(int),
            End=lambda x: x.End.astype(int),
            region=regions,
        )

    pr_a = pr.PyRanges(to_df(regions_a))
    pr_b = pr.PyRanges(to_df(regions_b))
    joined = pr_b.join(pr_a, suffix="_a")
    if joined.empty:
        return pd.DataFrame(columns=["region_a", "region_b"])
    return joined.df[["region", "region_a"]].rename(columns={"region": "region_b"})


def _ocoeff(df_a: pd.DataFrame, df_b: pd.DataFrame, on: list, use_overlap: bool = False) -> float:
    tmp_a, tmp_b = df_a.drop_duplicates(on), df_b.drop_duplicates(on)
    a_size, b_size = tmp_a.shape[0], tmp_b.shape[0]
    if (a_size > 0) and (b_size > 0):
        if use_overlap and len(on) == 1 and on[0] == "cre":
            regions_a = tmp_a["cre"].unique()
            regions_b = tmp_b["cre"].unique()
            mapping = _map_regions(regions_a, regions_b)
            if mapping.empty:
                i_size = 0
            else:
                # Count overlapping regions from the smaller set
                if a_size <= b_size:
                    i_size = mapping["region_a"].nunique()
                else:
                    i_size = mapping["region_b"].nunique()
        else:
            inter = pd.merge(tmp_a, tmp_b, on=on, how="inner")
            i_size = inter.shape[0]
        coeff = i_size / np.min([a_size, b_size])
    else:
        coeff = 0.0
    return coeff


def ocoeff(
    grns: dict,
) -> pd.DataFrame:
    """
    Compute pairwise overlap coefficient between GRNs at all levels.

    Parameters
    ----------
    grns
        Dictionary of GRNs with names as keys and DataFrames as values.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: grn_a, grn_b, source, cre, target, edge.
    """
    if len(grns) < 2:
        raise ValueError("At least 2 GRNs are required")

    names = list(grns.keys())
    rows = []
    for name_a, name_b in combinations(names, 2):
        grn_a, grn_b = grns[name_a], grns[name_b]
        source_oc = _ocoeff(grn_a, grn_b, on=["source"])
        cre_oc = (
            _ocoeff(grn_a, grn_b, on=["cre"], use_overlap=True)
            if "cre" in grn_a.columns and "cre" in grn_b.columns
            else np.nan
        )
        target_oc = _ocoeff(grn_a, grn_b, on=["target"])
        edge_oc = _ocoeff(grn_a, grn_b, on=["source", "target"])
        rows.append([name_a, name_b, source_oc, cre_oc, target_oc, edge_oc])

    return pd.DataFrame(rows, columns=["grn_a", "grn_b", "source", "cre", "target", "edge"])


def _get_grn_stats(grn: pd.DataFrame) -> tuple:
    n_s = grn["source"].nunique()
    n_c = grn["cre"].nunique() if "cre" in grn.columns else 0
    n_t = grn["target"].nunique()
    tdf = grn.drop_duplicates(["source", "target"])
    n_e = tdf.shape[0]
    n_r = tdf.groupby("source")["target"].count().mean()
    if np.isnan(n_r):
        n_r = 0.0
    return n_s, n_c, n_t, n_e, n_r


def stats(
    grns: pd.DataFrame | dict,
) -> pd.DataFrame:
    """
    Compute descriptive statistics for GRN(s).

    Parameters
    ----------
    grns
        Single GRN DataFrame or dictionary of GRNs with names as keys.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: name, n_sources, n_cres, n_targets, n_edges, mean_regulon_size.
    """
    if isinstance(grns, pd.DataFrame):
        grns = {"grn": grns}

    rows = []
    for name, grn in grns.items():
        n_s, n_c, n_t, n_e, n_r = _get_grn_stats(grn)
        rows.append([name, n_s, n_c, n_t, n_e, n_r])

    return pd.DataFrame(rows, columns=["name", "n_sources", "n_cres", "n_targets", "n_edges", "mean_regulon_size"])
