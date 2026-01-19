import mudata as mu
import numpy as np
import pandas as pd
import pyranges as pr
import scipy.stats as sts
from decoupler._download import _log

from pygreta.ds._db import read_db
from pygreta.pp._check import _check_organism


def correlation(
    mdata: mu.MuData,
    tfs: np.ndarray | list,
    organism: str = "hg38",
    method: str = "pearson",
    thr_r: float = 0.3,
    min_targets: int = 5,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Infer a GRN based on TF-gene expression correlation.

    Parameters
    ----------
    mdata
        MuData object with "rna" and "atac" modalities.
    tfs
        Array or list of transcription factor names.
    organism
        Which organism to use. Default is "hg38".
    method
        Correlation method: "pearson" or "spearman". Default is "pearson".
    thr_r
        Minimum absolute correlation threshold. Default is 0.3.
    min_targets
        Minimum number of targets required for a TF to be included. Default is 5.
    verbose
        Whether to print progress messages. Default is False.

    Returns
    -------
    pd.DataFrame
        GRN DataFrame with columns: source, cre, target, score.
    """
    _check_organism(organism)
    if not isinstance(mdata, mu.MuData):
        raise ValueError(f"mdata must be a MuData object, got {type(mdata)}")
    if not {"rna", "atac"}.issubset(mdata.mod):
        raise ValueError('MuData must contain "rna" and "atac" modalities')
    if method not in {"pearson", "spearman"}:
        raise ValueError(f'method must be "pearson" or "spearman", got {method}')

    genes = mdata.mod["rna"].var_names.astype("U")
    peaks = mdata.mod["atac"].var_names.astype("U")

    # Filter TFs to those present in the dataset
    tfs = np.array(list(set(genes) & set(tfs)))
    _log(f"Found {len(tfs)} TFs in dataset", level="info", verbose=verbose)

    # Compute correlation
    _log(f"Computing {method} correlation...", level="info", verbose=verbose)
    x = mdata.mod["rna"][:, tfs].X
    y = mdata.mod["rna"].X

    if method == "spearman":
        x = sts.rankdata(x, axis=0)
        y = sts.rankdata(y, axis=0)

    corr = np.corrcoef(x=x, y=y, rowvar=False)
    grn = pd.DataFrame(corr[: tfs.size, tfs.size :], index=tfs, columns=genes)
    grn = grn.reset_index(names="source").melt(id_vars="source", var_name="target", value_name="score")

    # Filter by threshold and remove self-regulation
    grn = grn[grn["score"].abs() > thr_r]
    grn = grn[grn["source"] != grn["target"]]
    _log(f"GRN edges after correlation filtering: {len(grn)}", level="info", verbose=verbose)

    # Transform peaks to PyRanges
    peaks_df = pd.DataFrame(peaks, columns=["cre"])
    peaks_df[["Chromosome", "Start", "End"]] = peaks_df["cre"].str.split("-", n=2, expand=True)
    peaks_df["Start"] = peaks_df["Start"].astype(int)
    peaks_df["End"] = peaks_df["End"].astype(int)
    peaks_pr = pr.PyRanges(peaks_df[["Chromosome", "Start", "End"]])

    # Load promoters and filter by genes
    _log("Downloading promoter annotations...", level="info", verbose=verbose)
    proms = read_db(organism=organism, db_name="Promoters", verbose=verbose)
    proms = proms[proms.Name.astype("U").isin(genes)]

    # Find promoters that overlap with peaks
    proms_nearest = proms.nearest(peaks_pr)
    proms_df = proms_nearest.df
    proms_df = proms_df[proms_df["Distance"] == 0]

    proms_df["cre"] = (
        proms_df["Chromosome"].astype(str) + "-" + proms_df["Start_b"].astype(str) + "-" + proms_df["End_b"].astype(str)
    )
    proms_df = proms_df[["cre", "Name"]].rename(columns={"Name": "target"}).drop_duplicates()

    # Merge GRN with promoter-peak mappings
    grn = pd.merge(grn, proms_df, how="inner")[["source", "cre", "target", "score"]]
    grn = grn.sort_values(["source", "target", "cre"])
    _log(f"GRN edges after peak filtering: {len(grn)}", level="info", verbose=verbose)

    # Filter TFs with less than min_targets targets
    n_targets = grn.groupby(["source"]).size().reset_index(name="counts")
    n_targets = n_targets[n_targets["counts"] >= min_targets]
    grn = grn[grn["source"].isin(n_targets["source"])]
    _log(f"GRN edges after min_targets filtering: {len(grn)}", level="info", verbose=verbose)

    grn = grn.reset_index(drop=True)

    return grn
