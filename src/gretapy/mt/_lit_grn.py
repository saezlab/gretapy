import mudata as mu
import pandas as pd
import pyranges as pr
from decoupler._download import _log

from gretapy.ds._db import read_db
from gretapy.pp._check import _check_organism

VALID_GRN_DATABASES = {"CollecTRI", "DoRothEA"}
REQUIRED_GRN_COLUMNS = {"source", "target", "weight"}


def lit_grn(
    mdata: mu.MuData,
    grn: str | pd.DataFrame = "CollecTRI",
    organism: str = "hg38",
    min_targets: int = 5,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Generate a GRN from literature-based prior knowledge filtered by available peaks and genes.

    This method accepts a reference GRN (CollecTRI, DoRothEA, or a custom DataFrame)
    and promoter annotations, then prunes the network based on the peaks and genes
    available in the input MuData object.

    Parameters
    ----------
    mdata
        MuData object with "rna" and "atac" modalities.
    grn
        Which GRN to use. Can be "CollecTRI", "DoRothEA", or a pandas DataFrame
        with columns: source, target, weight. Default is "CollecTRI".
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
        - score: Edge weight from the reference GRN.

    Examples
    --------
    >>> import mudata as mu
    >>> import gretapy as gt
    >>> mdata = mu.read("my_multiome.h5mu")
    >>> grn = gt.mt.lit_grn(mdata, grn="CollecTRI", organism="hg38")
    >>> grn.head()
    """
    # Validate inputs
    _check_organism(organism)
    if not isinstance(mdata, mu.MuData):
        raise ValueError(f"mdata must be a MuData object, got {type(mdata)}")
    if not {"rna", "atac"}.issubset(mdata.mod):
        raise ValueError('MuData must contain "rna" and "atac" modalities')

    # Resolve GRN
    if isinstance(grn, str):
        if grn not in VALID_GRN_DATABASES:
            raise ValueError(f"Invalid grn string '{grn}'. Must be one of {VALID_GRN_DATABASES} or a pandas DataFrame.")
        _log(f"Downloading {grn} GRN...", level="info", verbose=verbose)
        grn_df = read_db(organism=organism, db_name=grn, verbose=verbose)
    elif isinstance(grn, pd.DataFrame):
        missing = REQUIRED_GRN_COLUMNS - set(grn.columns)
        if missing:
            raise ValueError(f"Custom GRN DataFrame is missing required columns: {missing}")
        grn_df = grn.copy()
    else:
        raise ValueError(f"grn must be a string or pandas DataFrame, got {type(grn)}")

    # Extract genes and peaks from mdata
    genes = mdata.mod["rna"].var_names.astype("U")
    peaks = mdata.mod["atac"].var_names.astype("U")

    _log("Downloading promoter annotations...", level="info", verbose=verbose)
    proms = read_db(organism=organism, db_name="Promoters", verbose=verbose)

    # Transform peaks to PyRanges
    peaks_df = pd.DataFrame(peaks, columns=["cre"])
    peaks_df[["Chromosome", "Start", "End"]] = peaks_df["cre"].str.split("-", n=2, expand=True)
    peaks_df["Start"] = peaks_df["Start"].astype(int)
    peaks_df["End"] = peaks_df["End"].astype(int)
    peaks_pr = pr.PyRanges(peaks_df[["Chromosome", "Start", "End"]])

    # Filter GRN by available genes
    grn_df = grn_df[grn_df["source"].astype("U").isin(genes) & grn_df["target"].astype("U").isin(genes)]
    _log(f"GRN edges after gene filtering: {len(grn_df)}", level="info", verbose=verbose)

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
    grn_df = pd.merge(grn_df, proms_df, how="inner")[["source", "cre", "target", "weight"]]
    grn_df = grn_df.sort_values(["source", "target", "cre"]).rename(columns={"weight": "score"})
    _log(f"GRN edges after peak filtering: {len(grn_df)}", level="info", verbose=verbose)

    # Filter TFs with less than min_targets targets
    n_targets = grn_df.groupby(["source"]).size().reset_index(name="counts")
    n_targets = n_targets[n_targets["counts"] > min_targets]
    grn_df = grn_df[grn_df["source"].isin(n_targets["source"])]
    _log(f"GRN edges after min_targets filtering: {len(grn_df)}", level="info", verbose=verbose)

    grn_df = grn_df.reset_index(drop=True)

    return grn_df
