from itertools import combinations

import numpy as np
import pandas as pd
import scipy.stats as sts
from tqdm import tqdm

from gretapy.tl._utils import _prc_rcl_f01


def _grn(
    grn: pd.DataFrame,
    db: pd.DataFrame,
    genes: np.ndarray | list,
    verbose: bool = True,
) -> tuple:
    # Ensure uniqueness
    grn = grn[["source", "target"]].drop_duplicates(["source", "target"])
    # Filter resource by measured genes
    db = db[db["source"].astype("U").isin(genes) & db["target"].astype("U").isin(genes)]
    # Compute
    set_grn = set(grn["source"] + "|" + grn["target"])
    set_db = set(db["source"] + "|" + db["target"])
    tps = len(set_grn & set_db)
    fps = len(set_grn - set_db)
    fns = len(set_db - set_grn)
    # Compute metric
    prc, rcl, f01 = _prc_rcl_f01(tps=tps, fps=fps, fns=fns)
    return prc, rcl, f01


def _tfm(
    grn: pd.DataFrame,
    db: pd.DataFrame,
    genes: np.ndarray | list,
    cats: np.ndarray | list | None,
    verbose: bool = True,
) -> tuple:
    # Ensure uniqueness
    grn = grn[["source", "target"]].drop_duplicates(["source", "target"])
    # Filter db
    db = db[db[0].astype("U").isin(genes)]
    if cats is not None:
        db = db[db[1].str.contains("|".join(cats))]
    # Compute
    y_pred = set(grn["source"])
    y = set(db[0])
    tps = len(y_pred & y)
    fps = len(y_pred - y)
    fns = len(y - y_pred)
    # Compute metric
    prc, rcl, f01 = _prc_rcl_f01(tps=tps, fps=fps, fns=fns)
    return prc, rcl, f01


def _find_pairs(
    grn: pd.DataFrame,
    thr_pval: float,
    verbose: bool = True,
) -> set:
    # Pre-compute target sets per TF and total targets once
    targets_per_tf = grn.groupby("source")["target"].apply(set).to_dict()
    total = set(grn["target"])
    n_total = len(total)
    tfs = list(targets_per_tf)
    df = []
    for tf_a, tf_b in tqdm(list(combinations(tfs, r=2)), disable=not verbose, bar_format="{l_bar}{bar:20}{r_bar}"):
        trg_a = targets_per_tf[tf_a]
        trg_b = targets_per_tf[tf_b]
        a = len(trg_a & trg_b)
        if a > 0:
            b = len(trg_a) - a
            c = len(trg_b) - a
            d = n_total - len(trg_a | trg_b)
            s, p = sts.fisher_exact([[a, b], [c, d]], alternative="greater")
        else:
            s, p = 0.0, 1.0
        df.append([tf_a, tf_b, s, p])
    df = pd.DataFrame(df, columns=["tf_a", "tf_b", "stat", "pval"]).dropna()
    if df.shape[0] > 0:
        df["padj"] = sts.false_discovery_control(df["pval"], method="bh")
        df = df[df["padj"] < thr_pval]
        pairs = {"|".join(sorted([a, b])) for a, b in zip(df["tf_a"], df["tf_b"], strict=True)}
    else:
        pairs = set()
    return pairs


def _tfp(
    grn: pd.DataFrame,
    db: pd.DataFrame,
    genes: np.ndarray | list,
    thr_pval: float = 0.01,
    verbose: bool = True,
) -> tuple:
    # Ensure uniqueness
    grn = grn[["source", "target"]].drop_duplicates(["source", "target"])
    # Filter resource by measured genes
    db = db[db[0].astype("U").isin(genes) & db[1].astype("U").isin(genes)]
    # Filter grn
    tfs = set(db[0]) | set(db[1])
    grn = grn[grn["source"].isin(tfs)]
    db = {"|".join(sorted([a, b])) for a, b in zip(db[0], db[1], strict=True)}
    # Find pairs
    p_grn = _find_pairs(grn=grn, thr_pval=thr_pval, verbose=verbose)
    # Compute
    tps = len(p_grn & db)
    fps = len(p_grn - db)
    fns = len(db - p_grn)
    # Compute metric
    prc, rcl, f01 = _prc_rcl_f01(tps=tps, fps=fps, fns=fns)
    return prc, rcl, f01
