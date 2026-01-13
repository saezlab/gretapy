import anndata as ad
import mudata as mu
import pandas as pd
import pyranges as pr
from decoupler._log import _log

from pygreta.config import DATA, METRIC_CATS
from pygreta.ds._db import _read_db
from pygreta.pp._check import _check_dataset, _check_dts_grn, _check_grn, _check_metrics, _check_terms
from pygreta.tl._genomic import _cre, _cre_column
from pygreta.tl._mechanistic import _frc, _sim, _tfa
from pygreta.tl._predictive import _gset, _omics
from pygreta.tl._prior import _grn, _tfm, _tfp


def eval_grn_dataset(
    organism: str,
    grn: pd.DataFrame,
    dataset: str | mu.MuData | ad.AnnData,
    terms: dict | None,
    metrics: str | list | None = None,
    min_edges: int = 5,
) -> pd.DataFrame:
    """
    Evaluate a GRN against a dataset using multiple metrics.

    Parameters
    ----------
    organism
        Which organism to use (e.g., "hg38", "mm10").
    grn
        GRN DataFrame with columns "source", "target", and optionally "cre" and "score".
    dataset
        Dataset name (str) to load from config, or loaded MuData/AnnData object.
    terms
        Dictionary mapping database names to lists of terms for filtering.
        If None and dataset is str, terms are auto-loaded from config.
        Cannot be None if dataset is MuData/AnnData.
    metrics
        Metric(s) to evaluate. Can be category name, metric type, or database name.
        If None, all available metrics are evaluated.
    min_edges
        Minimum number of edges required in the GRN to run evaluation.
        GRNs with fewer edges will return an empty DataFrame.

    Returns
    -------
    DataFrame with columns: category, metric, db, precision, recall, f01.

    Example
    -------
    .. code-block:: python

        import pygreta as pg
        import pandas as pd

        grn = pd.read_csv("grn.csv")
        results = pg.tl.eval_grn_dataset(
            organism="hg38",
            grn=grn,
            dataset="pbmc10k",
            terms=None,
            metrics=None,
        )
    """
    result_cols = ["category", "metric", "db", "precision", "recall", "f01"]
    # Validate inputs
    metrics_list = _check_metrics(organism=organism, metrics=metrics)
    grn = _check_grn(grn=grn)
    if grn.shape[0] < min_edges:
        _log(
            f"GRN has {grn.shape[0]} edges, minimum required is {min_edges}. Returning empty results.",
            level="warning",
            verbose=True,
        )
        return pd.DataFrame(columns=result_cols)
    dataset = _check_dataset(organism=organism, dataset=dataset)
    terms = _check_terms(organism=organism, dataset=dataset, terms=terms)
    _check_dts_grn(dataset=dataset, grn=grn)
    # Check capabilities
    has_cre = "cre" in grn.columns
    is_mudata = isinstance(dataset, mu.MuData)
    can_run_genomic = has_cre and is_mudata
    if not has_cre:
        _log("GRN does not have 'cre' column. Genomic metrics will be skipped.", level="warning", verbose=True)
    if not is_mudata:
        _log("Dataset is AnnData (no ATAC modality). Genomic metrics will be skipped.", level="warning", verbose=True)
    # Extract data from dataset
    if is_mudata:
        genes, peaks, adata = (
            dataset.mod["rna"].var_names.tolist(),
            dataset.mod["atac"].var_names.tolist(),
            dataset.mod["rna"],
        )
    else:
        genes, peaks, adata = dataset.var_names.tolist(), [], dataset
    # Evaluate metrics
    results = []
    for db_name in metrics_list:
        db_info = DATA[organism]["dbs"].get(db_name)
        if db_info is None:
            continue
        metric_type, category = db_info["metric"], METRIC_CATS.get(db_info["metric"], "Unknown")
        # Handle metrics without file
        if db_info["fname"] is None:
            result = _run_fileless_metric(metric_type, db_name, dataset, grn, adata, is_mudata, has_cre)
            if result is not None:
                results.append([category, metric_type, db_name, *result])
            continue
        # Skip genomic metrics if not possible
        if metric_type in {"TF binding", "CREs", "CRE to gene links"} and not can_run_genomic:
            continue
        # Load database and run metric
        db = _read_db(organism=organism, db_name=db_name)
        cats = terms.get(db_name, None)
        result = _run_metric(metric_type, db_name, grn, db, genes, peaks, cats, adata)
        if result is not None:
            results.append([category, metric_type, db_name, *result])
    return pd.DataFrame(results, columns=result_cols)


def _run_metric(
    metric_type: str,
    db_name: str,
    grn: pd.DataFrame,
    db: pd.DataFrame | pr.PyRanges | ad.AnnData,
    genes: list,
    peaks: list,
    cats: list | None,
    adata: ad.AnnData,
) -> tuple | None:
    """Run a metric that requires a database file."""
    if metric_type == "Reference GRN":
        return _grn(grn=grn, db=db, genes=genes)
    elif metric_type == "TF markers":
        return _tfm(grn=grn, db=db, genes=genes, cats=cats)
    elif metric_type == "TF pairs":
        return _tfp(grn=grn, db=db)
    elif metric_type == "TF binding":
        return _cre_column(grn=grn, db=db, genes=genes, peaks=peaks, cats=cats, column="source")
    elif metric_type == "CREs":
        return _cre(grn=grn, db=db, peaks=peaks, cats=cats, reverse=(db_name == "ENCODE Blacklist"))
    elif metric_type == "CRE to gene links":
        return _cre_column(grn=grn, db=db, genes=genes, peaks=peaks, cats=cats, column="target")
    elif metric_type == "Gene sets":
        return _gset(adata=adata, grn=grn, db=db)
    elif metric_type == "TF scoring":
        return _tfa(adata=db, grn=grn, db=db, cats=cats)
    elif metric_type == "Perturbation forecasting":
        return _frc(adata=db, grn=grn, db=db, cats=cats)
    return None


def _run_fileless_metric(
    metric_type: str,
    db_name: str,
    dataset: mu.MuData | ad.AnnData,
    grn: pd.DataFrame,
    adata: ad.AnnData,
    is_mudata: bool,
    has_cre: bool,
) -> tuple | None:
    """Run metrics that don't require a database file (Omics, Boolean rules)."""
    if metric_type == "Omics":
        return _run_omics_metric(db_name, dataset, grn, is_mudata, has_cre)
    elif metric_type == "Steady state simulation":
        return _sim(adata=adata, grn=grn)
    return None


def _run_omics_metric(
    db_name: str,
    dataset: mu.MuData | ad.AnnData,
    grn: pd.DataFrame,
    is_mudata: bool,
    has_cre: bool,
) -> tuple | None:
    """Run omics metric based on the specific type."""
    if db_name == "gene ~ TFs":
        return _omics(
            data=dataset, grn=grn, col_source="source", col_target="target", mod_source="rna", mod_target="rna"
        )
    elif db_name == "gene ~ CREs" and is_mudata and has_cre:
        return _omics(data=dataset, grn=grn, col_source="cre", col_target="target", mod_source="atac", mod_target="rna")
    elif db_name == "CRE ~ TFs" and is_mudata and has_cre:
        return _omics(data=dataset, grn=grn, col_source="source", col_target="cre", mod_source="rna", mod_target="atac")
    elif db_name in {"gene ~ CREs", "CRE ~ TFs"}:
        _log(
            f"Skipping '{db_name}': requires MuData with ATAC and GRN with 'cre' column.", level="warning", verbose=True
        )
    return None
