import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm

from gretapy.tl._utils import _prc_rcl_f01

# TODO: assert that genes and cres in grn are in mdata input
# TODO: add grn.empty check in main function, not in individual metrics
# TODO: actually that at least 5 edges, else many errors in metrics (tfp, omics, etc)
# TODO: add support for when users only have RNA


def _grn_to_pr(
    grn: pd.DataFrame,
    column: str | None = None,
) -> pr.PyRanges:
    if grn.empty:
        return None
    pr_grn = grn.copy()
    pr_grn[["Chromosome", "Start", "End"]] = pr_grn["cre"].str.split("-", n=2, expand=True)
    if column is not None:
        pr_grn = pr_grn[["Chromosome", "Start", "End", column]].rename(columns={column: "Name"})
    else:
        pr_grn = pr_grn[["Chromosome", "Start", "End"]]
    pr_grn = pr_grn.drop_duplicates()
    pr_grn = pr.PyRanges(pr_grn)
    return pr_grn


def _peaks_to_pr(
    peaks: np.ndarray | list,
) -> pr.PyRanges:
    df_peaks = pd.DataFrame(peaks, columns=["cre"])
    df_peaks[["Chromosome", "Start", "End"]] = df_peaks["cre"].str.split("-", n=2, expand=True)
    pr_peaks = pr.PyRanges(df_peaks[["Chromosome", "Start", "End"]])
    return pr_peaks


def _cre_column(
    grn: pd.DataFrame,
    db: pd.DataFrame,
    genes: np.ndarray | list,
    peaks: np.ndarray | list,
    cats: np.ndarray | list | None,
    column: str,
    verbose: bool = True,
) -> tuple:
    # Filter db
    if cats is not None:
        cats_set = set(cats)
        mask = db.df["Score"].apply(lambda x: any(c.strip() in cats_set for c in x.split(",")))
        db = db[mask]
    db = db[db.df["Name"].astype("U").isin(genes)]
    pr_peaks = _peaks_to_pr(peaks)
    db = db.overlap(pr_peaks)
    # Remove features that are in GRN but not in db
    pr_grn = _grn_to_pr(grn, column=column)
    features_grn = pr_grn.df["Name"].unique()
    features_db = db.df["Name"].unique()
    features = np.setdiff1d(features_grn, features_db)
    pr_grn = pr_grn[~pr_grn.df["Name"].isin(features)]
    # Compute overlaps
    tps = 0
    fps = 0
    fns = 0
    for feature in tqdm(features_db, disable=not verbose, bar_format="{l_bar}{bar:20}{r_bar}"):
        f_grn = pr_grn[pr_grn.df["Name"] == feature]
        f_db = db[db.df["Name"] == feature]
        tps += f_grn.overlap(f_db).df.shape[0]
        fps += f_grn.overlap(f_db, invert=True).df.shape[0]
        fns += f_db.overlap(f_grn, invert=True).df.shape[0]
    # Compute metric
    prc, rcl, f01 = _prc_rcl_f01(tps=tps, fps=fps, fns=fns)
    return prc, rcl, f01


def _cre(
    grn: pd.DataFrame,
    db: pd.DataFrame,
    peaks: np.ndarray | list,
    cats: np.ndarray | list | None,
    reverse: bool = False,
    verbose: bool = True,
) -> tuple:
    # Filter db
    if cats is not None:
        cats_set = set(cats)
        mask = db.df["Score"].apply(lambda x: any(c.strip() in cats_set for c in x.split(",")))
        db = db[mask]
    pr_peaks = _peaks_to_pr(peaks)
    db = db.overlap(pr_peaks)
    # Compute
    pr_grn = _grn_to_pr(grn, column=None)
    if not reverse:
        tps = pr_grn.overlap(db).df.shape[0]
        fps = pr_grn.overlap(db, invert=True).df.shape[0]
        fns = db.overlap(pr_grn, invert=True).df.shape[0]
    else:
        tps = pr_grn.overlap(db, invert=True).df.shape[0]
        fps = pr_grn.overlap(db).df.shape[0]
        fns = pr_peaks.overlap(db, invert=True).df.shape[0]
    # Compute metric
    prc, rcl, f01 = _prc_rcl_f01(tps=tps, fps=fps, fns=fns)
    return prc, rcl, f01
