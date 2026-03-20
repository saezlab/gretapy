import logging
import multiprocessing

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.stats as sts
from sklearn.linear_model import Ridge
from tqdm import tqdm

from gretapy.tl._utils import _f_beta_score


def _get_pyboolnet():
    try:
        from pyboolnet import file_exchange, trap_spaces

        return file_exchange, trap_spaces
    except ImportError as e:
        raise ImportError(
            "pyboolnet is required for Boolean simulation but is not installed. "
            "Install it with: pip install git+https://github.com/hklarner/pyboolnet.git"
        ) from e


def _get_source_markers(
    adata: ad.AnnData,
    sources: np.ndarray | list,
    thr_deg_lfc: float,
    thr_deg_padj: float,
) -> pd.DataFrame:
    # Filter and update
    inter = adata.var_names.intersection(sources)
    sdata = adata[:, inter].copy()
    # Extract DEG tfs
    sc.tl.rank_genes_groups(sdata, groupby="celltype", method="wilcoxon")
    sm_df = sc.get.rank_genes_groups_df(sdata, group=None)
    # Filter results
    sm_df = sm_df[(sm_df["pvals_adj"] < thr_deg_padj) & (sm_df["logfoldchanges"] > thr_deg_lfc)]
    n_group = sm_df.groupby("group", as_index=False, observed=True).size()
    n_group = n_group[n_group["size"] >= 1]  # At least one TF per celltype
    groups = n_group["group"].values
    sm_df["group"] = sm_df["group"].astype(str)
    sm_df = sm_df[sm_df["group"].isin(groups)]
    sm_df = sm_df[["group", "names"]]
    sm_df.columns = ["celltype", "source"]
    return sm_df


def _define_bool_rules(
    grn: pd.DataFrame,
    indegree: int,
) -> str:
    targets = grn["target"].unique()
    rules = ""
    for target in targets:
        t_form = f"{target}, "
        msk = grn["target"] == target
        pos = ""
        neg = ""
        tgrn = grn.loc[msk]
        tgrn = tgrn.loc[tgrn["score"].abs().sort_values(ascending=False).index].head(indegree)
        for source, sign in zip(tgrn["source"], tgrn["score"], strict=True):
            if sign >= 0:
                pos += f"{source} | "
            else:
                neg += f"!{source} & "
        pos, neg = pos[:-3], neg[:-3]
        if (pos != "") and (neg != ""):
            t_form += f"({pos}) & ({neg})"
        elif (pos != "") and (neg == ""):
            t_form += f"{pos}"
        elif (pos == "") and (neg != ""):
            t_form += f"{neg}"
        rules += t_form + "\n"
    return rules


def _fisher_test(
    hits: set,
    sm_set: set,
    sources: set,
) -> float:
    a = hits & sm_set
    b = hits - sm_set
    c = sm_set - hits
    d = sources - a - b - c
    a, b, c, d = len(a), len(b), len(c), len(d)
    _, p = sts.fisher_exact([[a, b], [c, d]], alternative="greater")
    return p


def _get_sim_hits(
    sss: list,
    sm_sets: pd.Series,
    sources: set,
    thr_fisher_padj: float,
    verbose: bool = True,
) -> pd.DataFrame:
    pvals = np.zeros((len(sss), sm_sets.shape[0]))
    for i, pss in enumerate(tqdm(sss, disable=not verbose, bar_format="{l_bar}{bar:20}{r_bar}")):
        hits = set()
        for k in pss:
            if pss[k]:
                hits.add(k)
        for j, sm_set in enumerate(sm_sets):
            pvals[i, j] = _fisher_test(hits=hits, sm_set=sm_set, sources=sources)
    pvals = sts.false_discovery_control(pvals, axis=1)
    hits_df = pd.DataFrame(pvals < thr_fisher_padj, columns=sm_sets.index)
    return hits_df


def _sss_prc_rcl(hits: pd.DataFrame) -> tuple:
    # prc at sss level
    n_ct_per_ss = hits.sum(1)
    tp = (n_ct_per_ss == 1).sum()
    if tp > 0:
        fp = (n_ct_per_ss != 1).sum()
        prc = tp / (tp + fp)
    else:
        prc = 0
    # rcl at celltype level
    n_cts = hits.sum(0)
    tp = (n_cts > 0).sum()
    if tp > 0:
        tn = (n_cts == 0).sum()
        rcl = tp / (tp + tn)
    else:
        rcl = 0
    return prc, rcl


def _run_sim_block(queue, sgrn, indegree, sm_sets, f_sources, thr_fisher_padj, verbose):
    """Run the Boolean simulation block in a subprocess."""
    try:
        file_exchange, trap_spaces = _get_pyboolnet()
        bool_rules = _define_bool_rules(grn=sgrn, indegree=indegree)
        primes = file_exchange.bnet2primes(bool_rules)
        logging.getLogger("pyboolnet.external.potassco").disabled = True
        sss = trap_spaces.compute_steady_states(primes, max_output=100_000)
        hits = _get_sim_hits(
            sss=sss, sm_sets=sm_sets, sources=f_sources, thr_fisher_padj=thr_fisher_padj, verbose=verbose
        )
        prc, rcl = _sss_prc_rcl(hits=hits)
        f01 = _f_beta_score(prc=prc, rcl=rcl)
        queue.put((prc, rcl, f01))
    except (OSError, ValueError, RuntimeError):
        queue.put((np.nan, np.nan, np.nan))


def _sim(
    adata: ad.AnnData,
    grn: pd.DataFrame,
    indegree: int = 10,
    thr_deg_lfc: float = 2,
    thr_deg_padj: float = 2.22e-16,
    thr_fisher_padj: float = 0.01,
    verbose: bool = True,
    timeout: int = 3600,
) -> tuple:
    # Ensure uniqueness but keep score
    grn = grn.groupby(["source", "target"], as_index=False)["score"].mean()
    # Find TF markers
    sources = grn["source"].unique()
    sm_df = _get_source_markers(
        adata=adata,
        sources=sources,
        thr_deg_lfc=thr_deg_lfc,
        thr_deg_padj=thr_deg_padj,
    )
    # Find sets and filter grn
    sm_sets = sm_df.groupby("celltype")["source"].apply(lambda x: set(x))
    f_sources = set(sm_df["source"])
    sgrn = grn[(grn["source"].isin(f_sources)) & (grn["target"].isin(f_sources))]
    # Simulate
    if sgrn.shape[0] >= 5:  # Remove outlier networks
        queue = multiprocessing.Queue()
        proc = multiprocessing.Process(
            target=_run_sim_block,
            args=(queue, sgrn, indegree, sm_sets, f_sources, thr_fisher_padj, verbose),
        )
        proc.start()
        proc.join(timeout=timeout)
        if proc.is_alive():
            proc.terminate()
            proc.join()
            prc, rcl, f01 = np.nan, np.nan, np.nan
        else:
            prc, rcl, f01 = queue.get()
    else:
        prc, rcl, f01 = np.nan, np.nan, np.nan
    return prc, rcl, f01


def _tfa(
    adata: ad.AnnData,
    grn: pd.DataFrame,
    db: ad.AnnData,
    cats: list | None,
    thr_pert_lfc: float = -0.5,
    thr_score_padj: float = 0.05,
    verbose: bool = True,
) -> tuple:
    # Ensure uniqueness but keep score
    grn = grn.groupby(["source", "target"], as_index=False)["score"].mean()
    # Filter db
    msk = (db.obs["source"].isin(adata.var_names)) & (db.obs["logFC"] < thr_pert_lfc)
    if cats is not None:
        msk = msk & db.obs["Tissue.Type"].isin(cats)
    db = db[msk, :].copy()
    # Infer enrichment activity scores
    scores = []
    pvals = []
    for dataset in tqdm(db.obs.index, disable=not verbose, bar_format="{l_bar}{bar:20}{r_bar}"):
        source = db.obs.loc[dataset, "source"]
        source_mat = db[[dataset], :].to_df()
        source_grn = grn[grn["source"] == source].rename(columns={"score": "weight"})
        if source_grn.shape[0] >= 3:
            score, pval = dc.mt.ulm(data=source_mat, net=source_grn, tmin=3)
            score, pval = float(score.values[0, 0]), float(pval.values[0, 0])
            scores.append(score)
            pvals.append(pval)
    scores = np.array(scores)
    pvals = np.array(pvals)
    padj = sts.false_discovery_control(pvals, method="bh")
    # Compute metric
    tp = np.sum((scores < 0) & (padj < thr_score_padj))
    if tp > 0:
        prc = tp / scores.size
        rcl = tp / db.shape[0]
        f01 = _f_beta_score(prc=prc, rcl=rcl)
    else:
        prc, rcl, f01 = 0.0, 0.0, 0.0
    return prc, rcl, f01


def _coefmat(
    adata: ad.AnnData,
    grn: pd.DataFrame,
    alpha: float = 1,
    seed: int = 42,
    smin: int = 3,
    verbose: bool = False,
):
    coefmat = np.zeros((adata.var_names.size, adata.var_names.size))
    coefmat = pd.DataFrame(coefmat, index=adata.var_names, columns=adata.var_names)
    for target in tqdm(adata.var_names, disable=not verbose, bar_format="{l_bar}{bar:20}{r_bar}"):
        t_grn = grn[grn["target"] == target]
        sources = [s for s in t_grn["source"].unique() if s != target]
        if len(sources) >= smin:
            model = Ridge(alpha=alpha, random_state=seed)
            model.fit(adata[:, sources].X, adata[:, target].X)
            coefmat.loc[target, sources] = model.coef_
    return coefmat


def _simulate(
    df: pd.DataFrame,
    coefmat: pd.DataFrame,
    gene: str,
    n_steps: int = 3,
):
    # Initial diff after knockdown
    delta_init = df.copy()
    delta_init[gene] = 0.0
    delta_init = delta_init - df
    # Propagate
    delta_sim = delta_init.copy()
    for _ in range(n_steps):
        delta_sim = delta_sim.dot(coefmat.T)
        delta_sim[delta_init != 0] = delta_init
        gex_tmp = df + delta_sim
        gex_tmp[gex_tmp < 0] = 0.0
        delta_sim = gex_tmp - df
    # Remove invariant genes
    delta_sim = delta_sim.loc[:, delta_sim.abs().sum(0) != 0]
    return delta_sim


def _frc(
    adata: ad.AnnData,
    grn: pd.DataFrame,
    db: ad.AnnData,
    cats: list | None,
    thr_pert_lfc: float = -0.5,
    n_steps: int = 3,
    min_size: int = 10,
    thr_cor_stat: float = 0.05,
    thr_cor_padj: float = 0.05,
    verbose: bool = True,
) -> float:
    # Ensure uniqueness, keep first occurrence
    grn = grn.drop_duplicates(["source", "target"], keep="first")
    grn_sources = set(grn["source"])
    # Filter for grn genes
    g_universe = set(grn["source"]) | set(grn["target"])
    fdata = adata[:, adata.var_names.isin(g_universe)].copy()
    # Filter db
    msk = (db.obs["source"].isin(adata.var_names)) & (db.obs["logFC"] < thr_pert_lfc)
    if cats is not None:
        msk = msk & db.obs["Tissue.Type"].isin(cats)
    g_universe = list(g_universe & set(db.var_names))
    fdb = db[msk][:, g_universe].copy()
    # Fit grn with OLS (alpha=0), allow single-regulator targets
    coefmat = _coefmat(adata=fdata, grn=grn, alpha=0, smin=1)
    # Determine valid TFs (targets with at least 3 non-zero regulators)
    valid_tfs = set(coefmat.index[(coefmat != 0).sum(axis=1) >= 3])
    # Generate seed gex profile
    profile = fdata.to_df().mean(0).to_frame().T
    # Count GRN source experiments for precision denominator
    n_grn_tfs_in_prt = fdb.obs["source"].isin(grn_sources).sum()
    # Simulate perturbations per tf
    coefs = []
    pvals = []
    for dataset in tqdm(fdb.obs_names, disable=not verbose, bar_format="{l_bar}{bar:20}{r_bar}"):
        # Extract real lfc
        tf = fdb.obs.loc[dataset, "source"]
        if tf in valid_tfs:
            y = fdb[[dataset], :].to_df()
            y = y[y != 0].dropna(axis=1)
            # Run perturbation simulation for the current TF
            x = _simulate(profile, coefmat, gene=tf, n_steps=n_steps)
            # Intersect
            inter = np.intersect1d(x.columns, y.columns)
            x, y = x.loc[:, inter].values[0], y.loc[:, inter].values[0]
            # Compute correlation
            if x.size >= min_size:
                r, p = sts.spearmanr(x, y)
            else:
                r, p = 0.0, 1.0
            coefs.append(r)
            pvals.append(p)
    coefs = np.array(coefs)
    pvals = np.array(pvals)
    padj = sts.false_discovery_control(pvals, method="bh")
    # Compute metrics
    tp = np.sum((coefs > thr_cor_stat) & (padj < thr_cor_padj))
    if tp > 0:
        prc = tp / n_grn_tfs_in_prt
        rcl = tp / fdb.obs.shape[0]
        f01 = _f_beta_score(prc=prc, rcl=rcl)
    else:
        prc, rcl, f01 = 0.0, 0.0, 0.0
    return prc, rcl, f01
