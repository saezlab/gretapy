import mudata as mu
import numpy as np
import pandas as pd
import pyranges as pr
from decoupler._download import _log

from gretapy.ds._db import read_db
from gretapy.mt._utils import _trim_grn
from gretapy.pp._check import _check_organism


def _get_cres_pr(peaks: np.ndarray) -> pr.PyRanges:
    cres = [c.split("-") for c in peaks]
    return pr.PyRanges(pd.DataFrame(cres, columns=["Chromosome", "Start", "End"]))


def _get_window(gannot: pr.PyRanges, target: str, w_size: int) -> pr.PyRanges:
    target_gr = gannot[gannot.df["Name"] == target]
    tss = target_gr.df["Start"].values[0] + 1000
    return pr.from_dict(
        {
            "Chromosome": target_gr.Chromosome,
            "Start": tss - w_size,
            "End": tss + w_size,
        }
    )


def _get_overlap_cres(gene: str, gannot: pr.PyRanges, cres_pr: pr.PyRanges, w_size: int) -> np.ndarray | None:
    wnd = _get_window(gannot, target=gene, w_size=w_size)
    o_cres = cres_pr.overlap(wnd)
    if len(o_cres) > 0:
        return o_cres.df.assign(
            name=lambda x: x["Chromosome"].astype(str) + "-" + x["Start"].astype(str) + "-" + x["End"].astype(str)
        )["name"].values
    return None


def random(
    mdata: mu.MuData,
    organism: str = "hg38",
    tfs: np.ndarray | list | None = None,
    g_perc: float = 0.25,
    scale: float = 1.0,
    tf_g_ratio: float = 0.10,
    w_size: int = 250000,
    min_targets: int = 5,
    seed: int = 0,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Generate a random GRN for benchmarking purposes.

    Parameters
    ----------
    mdata
        MuData object with "rna" and "atac" modalities.
    organism
        Which organism to use. Default is "hg38".
    tfs
        Array or list of transcription factor names. If None, uses LambertTFs.
    g_perc
        Percentage of genes to include. Default is 0.25.
    scale
        Scale parameter for exponential distribution sampling. Default is 1.0.
    tf_g_ratio
        Ratio of TFs to genes. Default is 0.10.
    w_size
        Window size around TSS for peak overlap. Default is 250000.
    min_targets
        Minimum number of targets required for a TF to be included. Default is 5.
    seed
        Random seed for reproducibility. Default is 42.
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

    # Extract names from mdata
    genes = mdata.mod["rna"].var_names.values.astype("U")
    peaks = mdata.mod["atac"].var_names.values.astype("U")

    # Shuffle (pre phase) — own rng
    rng_pre = np.random.default_rng(seed=seed)
    genes = rng_pre.choice(genes, genes.size, replace=False)
    peaks = rng_pre.choice(peaks, peaks.size, replace=False)

    # Get TFs
    if tfs is None:
        _log("Loading TFs from LambertTFs...", level="info", verbose=verbose)
        tfs = read_db(organism=organism, db_name="Lambert TFs", verbose=verbose)
    tfs = np.intersect1d(genes, tfs)
    _log(f"Found {len(tfs)} TFs in dataset", level="info", verbose=verbose)

    _log("Downloading promoter annotations...", level="info", verbose=verbose)
    gannot = read_db(organism=organism, db_name="Promoters", verbose=verbose)

    # Filter genes to those in annotation
    g_in_ann = set(gannot.df["Name"].values)
    genes = genes[np.isin(genes, list(g_in_ann))]
    _log(f"Genes in annotation: {len(genes)}", level="info", verbose=verbose)

    # Get peaks as PyRanges
    cres_pr = _get_cres_pr(peaks)

    # p2g phase — own rng
    rng = np.random.default_rng(seed=seed)

    # Randomly sample genes
    n_genes = int(np.round(len(genes) * g_perc))
    sampled_genes = rng.choice(genes, n_genes, replace=False)
    n_cres_per_gene = np.ceil(rng.exponential(scale=scale, size=n_genes)).astype(int)

    # Generate peak-gene connections
    _log("Generating random peak-gene connections...", level="info", verbose=verbose)
    p2g_rows = []
    for i, gene in enumerate(sampled_genes):
        o_cres = _get_overlap_cres(gene, gannot, cres_pr, w_size)
        if o_cres is not None:
            n_cre = min(n_cres_per_gene[i], len(o_cres))
            r_cres = rng.choice(o_cres, n_cre, replace=False)
            for cre in r_cres:
                p2g_rows.append([cre, gene])

    if not p2g_rows:
        _log("No peak-gene connections found", level="warning", verbose=verbose)
        return pd.DataFrame(columns=["source", "cre", "target", "score"])

    p2g = (
        pd.DataFrame(p2g_rows, columns=["cre", "target"])
        .sort_values(["cre", "target"])
        .drop_duplicates(["cre", "target"])
    )

    # tfb phase — own rng
    rng = np.random.default_rng(seed=seed)

    # Generate TF-peak connections
    _log("Generating random TF-peak connections...", level="info", verbose=verbose)
    cres = p2g["cre"].unique()
    n_tfs_per_cre = np.ceil(rng.exponential(scale=scale, size=len(cres))).astype(int)

    tfb_rows = []
    for i, cre in enumerate(cres):
        n_tf = n_tfs_per_cre[i]
        r_tfs = rng.choice(tfs, n_tf)
        for tf in r_tfs:
            tfb_rows.append([cre, tf])

    tfb = (
        pd.DataFrame(tfb_rows, columns=["cre", "source"])
        .sort_values(["cre", "source"])
        .drop_duplicates(["cre", "source"])
    )

    # Merge to create GRN
    grn = pd.merge(tfb[["source", "cre"]], p2g[["cre", "target"]], how="inner", on="cre")[["source", "cre", "target"]]
    grn = grn.sort_values(["source", "cre", "target"])
    grn["score"] = 1.0

    # mdl phase — own rng
    rng = np.random.default_rng(seed=seed)

    # Subsample TFs based on ratio
    unique_tfs = grn["source"].unique().astype("U")
    n_tfs = int(np.round(len(grn["target"].unique()) * tf_g_ratio))
    n_tfs = min(n_tfs, len(unique_tfs))
    selected_tfs = rng.choice(unique_tfs, n_tfs, replace=False)
    grn = grn[grn["source"].astype("U").isin(selected_tfs)]
    _log(f"GRN edges before filtering: {len(grn)}", level="info", verbose=verbose)

    # Filter TFs with less than min_targets targets
    n_targets = grn.groupby(["source"]).size().reset_index(name="counts")
    n_targets = n_targets[n_targets["counts"] > min_targets]
    grn = grn[grn["source"].isin(n_targets["source"])]
    _log(f"GRN edges after min_targets filtering: {len(grn)}", level="info", verbose=verbose)

    # Trim
    grn = _trim_grn(grn)

    return grn.reset_index(drop=True)
