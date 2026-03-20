import anndata as ad
import decoupler as dc
import mudata as mu
import numpy as np
import pandas as pd
import scipy.stats as sts
from sklearn.model_selection import train_test_split
from tqdm import tqdm
from xgboost import XGBRegressor

from gretapy.tl._utils import _f_beta_score, _prc_rcl_f01


def _remove_zeros(
    X: np.ndarray,
    y: np.ndarray,
) -> tuple:
    msk = y != 0.0
    y = y[msk]
    X = X[msk, :]
    return X, y


def _extract_data(
    data: mu.MuData | ad.AnnData,
    train_obs_names: np.ndarray,
    test_obs_names: np.ndarray,
    target: str,
    sources: np.ndarray | list,
    mod_source: str | None,
    mod_target: str | None,
) -> tuple:
    if isinstance(data, mu.MuData):
        train_X = data.mod[mod_source][train_obs_names, sources].X
        train_y = data.mod[mod_target][train_obs_names, target].X.ravel()
        test_X = data.mod[mod_source][test_obs_names, sources].X
        test_y = data.mod[mod_target][test_obs_names, target].X.ravel()
    elif isinstance(data, ad.AnnData):
        train_X = data[train_obs_names, sources].X
        train_y = data[train_obs_names, target].X.ravel()
        test_X = data[test_obs_names, sources].X
        test_y = data[test_obs_names, target].X.ravel()
    train_X, train_y = _remove_zeros(train_X, train_y)
    test_X, test_y = _remove_zeros(test_X, test_y)
    return train_X, train_y, test_X, test_y


def _test_predictability(
    data: mu.MuData | ad.AnnData,
    train_obs_names: np.ndarray,
    test_obs_names: np.ndarray,
    grn: pd.DataFrame,
    col_source: str,
    col_target: str,
    mod_source: str | None,
    mod_target: str | None,
    ntop: int = 5,
    verbose: bool = True,
) -> pd.DataFrame:
    # Make unique and select ntop sources per target
    net = grn.drop_duplicates([col_source, col_target]).copy()
    net["abs_score"] = abs(net["score"])
    net = net.sort_values("abs_score", ascending=False)
    net = net.groupby(col_target)[col_source].apply(lambda x: list(x) if ntop is None else list(x)[:ntop])
    cor = []
    for target in tqdm(net.index, disable=not verbose, bar_format="{l_bar}{bar:20}{r_bar}"):
        sources = net[target]
        sources = [s for s in sources if s != target]  # remove self loop
        if len(sources) >= 1:  # Needed if self loop is removed
            train_X, train_y, test_X, test_y = _extract_data(
                data=data,
                train_obs_names=train_obs_names,
                test_obs_names=test_obs_names,
                target=target,
                sources=sources,
                mod_source=mod_source,
                mod_target=mod_target,
            )
            if test_y.size >= 10:  # Needed if zeros are removed
                reg = XGBRegressor(random_state=0, n_jobs=1).fit(train_X, train_y)
                pred_y = reg.predict(test_X)
                if np.any(pred_y != pred_y[0]):
                    s, p = sts.spearmanr(pred_y, test_y)  # Spearman to control for outliers
                    cor.append([target, pred_y.size, len(sources), s, p])
    # Format to df
    cor = pd.DataFrame(cor, columns=["target", "n_obs", "n_vars", "coef", "pval"])
    if cor.shape[0] > 0:
        cor["padj"] = sts.false_discovery_control(cor["pval"], method="bh")
    else:
        cor["padj"] = pd.Series(dtype=float)
    return cor


def _omics(
    data: mu.MuData | ad.AnnData,
    grn: pd.DataFrame,
    col_source: str,
    col_target: str,
    mod_source: str | None,
    mod_target: str | None,
    test_size: float = 0.33,
    seed: int = 42,
    ntop: int = 5,
    verbose: bool = True,
):
    # Split by train test
    train_obs_names, test_obs_names = train_test_split(
        data.obs_names,
        test_size=test_size,
        random_state=seed,
        stratify=data.obs["celltype"],
    )
    # Compute
    cor = _test_predictability(
        data=data,
        train_obs_names=train_obs_names,
        test_obs_names=test_obs_names,
        grn=grn,
        col_source=col_source,
        col_target=col_target,
        mod_source=mod_source,
        mod_target=mod_target,
        ntop=ntop,
        verbose=verbose,
    )
    sig_cor = cor[(cor["padj"] < 0.05) & (cor["coef"] > 0.05)]
    n_hits = sig_cor.shape[0]
    # Compute metric
    if n_hits > 0:
        if isinstance(data, mu.MuData):
            universe_size = data.mod[mod_target].var_names.size
        else:
            universe_size = data.var_names.size
        prc = n_hits / cor.shape[0]
        rcl = n_hits / universe_size
        f01 = _f_beta_score(prc=prc, rcl=rcl)
    else:
        prc, rcl, f01 = 0.0, 0.0, 0.0
    return prc, rcl, f01


def _ora_overlap_fdr(features, pw_targets, n_bg=20000):
    """Run ORA with FDR applied only to pathways with overlap > 0."""
    features_set = set(features)
    results = []
    for source, targets in pw_targets.items():
        a = len(features_set & targets)
        if a > 0:
            b = len(targets) - a
            c = len(features_set) - a
            d = int(n_bg - a - b - c)
            _, pv = sts.fisher_exact([[a, b], [c, d]], alternative="greater")
            results.append([source, pv])
    df = pd.DataFrame(results, columns=["source", "pval"])
    df["padj"] = sts.false_discovery_control(df["pval"], method="bh")
    return df


def _gset(
    adata: ad.AnnData,
    grn: pd.DataFrame,
    db: pd.DataFrame,
    thr_pval: float = 0.01,
    thr_prop: float = 0.20,
    verbose: bool = True,
) -> tuple:
    # Ensure uniqueness
    grn = grn[["source", "target"]].drop_duplicates(["source", "target"])
    # Infer pathway enrichment scores
    dc.mt.ulm(data=adata, net=db)
    scores = adata.obsm["score_ulm"]
    # Recompute raw p-values from t-values (matching ULM internals)
    df_t = adata.shape[1] - 2
    raw_pvals = 2 * sts.t.sf(np.abs(scores.values), df_t)
    # Apply global FDR correction (old behavior: flatten, correct, reshape)
    flat_padj = sts.false_discovery_control(raw_pvals.ravel(), method="bh")
    padj = pd.DataFrame(flat_padj.reshape(raw_pvals.shape), index=scores.index, columns=scores.columns)
    # Find pathway hits
    hits = ((padj < thr_pval) & (scores > 0)).sum(0)
    hits = hits.sort_values(ascending=False) / padj.shape[0]
    hits = hits[hits > thr_prop].index.values.astype("U")
    # Find pathway hits in grn
    pw_targets = {src: set(tgts) for src, tgts in db.groupby("source")["target"]}
    sig_pws = set()
    for source in tqdm(grn["source"].unique(), disable=not verbose, bar_format="{l_bar}{bar:20}{r_bar}"):
        features = grn[grn["source"] == source]["target"]
        pws = _ora_overlap_fdr(features=features, pw_targets=pw_targets)
        sig_pws.update(pws[pws["padj"] < thr_pval]["source"])
    sig_pws = np.array(list(sig_pws))
    # Compute
    tps = np.intersect1d(sig_pws, hits).size
    fps = np.setdiff1d(sig_pws, hits).size
    fns = np.setdiff1d(hits, sig_pws).size
    # Compute metric
    prc, rcl, f01 = _prc_rcl_f01(tps=tps, fps=fps, fns=fns)
    return prc, rcl, f01
