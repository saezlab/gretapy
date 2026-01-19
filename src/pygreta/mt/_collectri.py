import mudata as mu
import pandas as pd
import pyranges as pr
from decoupler._download import _log

from pygreta.ds._db import read_db
from pygreta.pp._check import _check_organism


def collectri(
    mdata: mu.MuData,
    organism: str = "hg38",
    min_targets: int = 5,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Generate a GRN from CollecTRI prior knowledge filtered by available peaks and genes.

    This method downloads the CollecTRI reference GRN and promoter annotations,
    then prunes the network based on the peaks and genes available in the input
    MuData object.

    Parameters
    ----------
    mdata
        MuData object with "rna" and "atac" modalities.
    organism
        Which organism to use. Default is "hg38".
    min_targets
        Minimum number of targets required for a TF to be included. Default is 5.
    verbose
        Whether to print progress messages. Default is False.

    Returns
    -------
    pd.DataFrame
        GRN DataFrame with columns: source, cre, target, score.
        - source: Transcription factor name.
        - cre: Cis-regulatory element (peak) in format chrX-start-end.
        - target: Target gene name.
        - score: Edge weight from CollecTRI.

    Examples
    --------
    >>> import mudata as mu
    >>> import pygreta as pg
    >>> mdata = mu.read("my_multiome.h5mu")
    >>> grn = pg.mt.collectri(mdata, organism="hg38")
    >>> grn.head()
    """
    # Validate inputs
    _check_organism(organism)
    if not isinstance(mdata, mu.MuData):
        raise ValueError(f"mdata must be a MuData object, got {type(mdata)}")
    if not {"rna", "atac"}.issubset(mdata.mod):
        raise ValueError('MuData must contain "rna" and "atac" modalities')

    # Extract genes and peaks from mdata
    genes = mdata.mod["rna"].var_names.astype("U")
    peaks = mdata.mod["atac"].var_names.astype("U")

    _log("Downloading CollecTRI GRN...", level="info", verbose=verbose)
    grn = read_db(organism=organism, db_name="CollecTRI", verbose=verbose)

    _log("Downloading promoter annotations...", level="info", verbose=verbose)
    proms = read_db(organism=organism, db_name="Promoters", verbose=verbose)

    # Transform peaks to PyRanges
    peaks_df = pd.DataFrame(peaks, columns=["cre"])
    peaks_df[["Chromosome", "Start", "End"]] = peaks_df["cre"].str.split("-", n=2, expand=True)
    peaks_df["Start"] = peaks_df["Start"].astype(int)
    peaks_df["End"] = peaks_df["End"].astype(int)
    peaks_pr = pr.PyRanges(peaks_df[["Chromosome", "Start", "End"]])

    # Filter GRN by available genes
    grn = grn[grn["source"].astype("U").isin(genes) & grn["target"].astype("U").isin(genes)]
    _log(f"GRN edges after gene filtering: {len(grn)}", level="info", verbose=verbose)

    # Filter promoters by available genes
    proms = proms[proms.Name.astype("U").isin(genes)]

    # Find promoters that overlap with peaks
    proms_nearest = proms.nearest(peaks_pr)
    proms_df = proms_nearest.df
    proms_df = proms_df[proms_df["Distance"] == 0]

    # Create CRE column from overlapping peak coordinates
    proms_df["cre"] = (
        proms_df["Chromosome"].astype(str) + "-" + proms_df["Start_b"].astype(str) + "-" + proms_df["End_b"].astype(str)
    )
    proms_df = proms_df[["cre", "Name"]].rename(columns={"Name": "target"}).drop_duplicates()

    # Merge GRN with promoter-peak mappings
    grn = pd.merge(grn, proms_df, how="inner")[["source", "cre", "target", "weight"]]
    grn = grn.sort_values(["source", "target", "cre"]).rename(columns={"weight": "score"})
    _log(f"GRN edges after peak filtering: {len(grn)}", level="info", verbose=verbose)

    # Filter TFs with less than min_targets targets
    n_targets = grn.groupby(["source"]).size().reset_index(name="counts")
    n_targets = n_targets[n_targets["counts"] >= min_targets]
    grn = grn[grn["source"].isin(n_targets["source"])]
    _log(f"GRN edges after min_targets filtering: {len(grn)}", level="info", verbose=verbose)

    grn = grn.reset_index(drop=True)

    return grn
