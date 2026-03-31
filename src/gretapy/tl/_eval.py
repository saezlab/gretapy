import time

import anndata as ad
import mudata as mu
import pandas as pd
import pyranges as pr
from decoupler._log import _log

from gretapy.config import DATA, METRIC_CATS
from gretapy.ds._db import read_db
from gretapy.pp._check import (
    _check_dataset,
    _check_dts_grn,
    _check_grn,
    _check_metrics,
    _check_organism,
    _check_terms,
)
from gretapy.tl._genomic import _cre, _cre_column
from gretapy.tl._mechanistic import _frc, _sim, _tfa
from gretapy.tl._predictive import _gset, _omics
from gretapy.tl._prior import _grn, _tfm, _tfp

_SEP = "\u2550" * 50


def _format_log_prefix(grn_name: str | None = None, dataset_name: str | None = None) -> str:
    """Build the optional bracket prefix for log messages."""
    parts = []
    if grn_name is not None:
        parts.append(grn_name)
    if dataset_name is not None:
        parts.append(dataset_name)
    if parts:
        return f"[{' | '.join(parts)}] "
    return ""


def _format_label(grn_name: str | None = None, dataset_name: str | None = None) -> str:
    """Build a label string from available names."""
    parts = []
    if grn_name is not None:
        parts.append(grn_name)
    if dataset_name is not None:
        parts.append(dataset_name)
    return " | ".join(parts) if parts else ""


def benchmark(
    grns: dict,
    organism: str | None = None,
    datasets: list | dict | None = None,
    terms: dict | None = None,
    metrics: str | list | None = None,
    min_edges: int = 5,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Run the benchmark for one or multiple GRNs across one or multiple datasets.

    Parameters
    ----------
    grns
        Dictionary mapping GRN names to per-organism per-dataset GRN DataFrames.
        Structure: ``{grn_name: {organism: {dataset_name: DataFrame}}}``.
    organism
        Ignored when organism keys are present in ``grns``.  Kept for clarity
        but organisms are inferred from the second level of ``grns``.
    datasets
        Dataset(s) to evaluate against. Can be:
        - None: Use all datasets present in the grns dict for each organism.
        - list: A whitelist of dataset names (applied across all organisms).
        - dict: A flat dictionary mapping dataset names to pre-loaded MuData/AnnData objects.
    terms
        Optional dictionary specifying filtering terms per organism, dataset, and metric.
        Structure: ``{organism: {dataset_name: {db_name: [terms]}}}``.
        If None, terms are auto-loaded from config for each dataset.
    metrics
        Metric(s) to evaluate. Can be category name, metric type, or database name.
        If None, all available metrics are evaluated.
    min_edges
        Minimum number of edges required in a GRN to run evaluation.
    verbose
        Whether to log progress messages and show progress bars.

    Returns
    -------
    DataFrame with columns: grn, organism, dataset, class, task, db, precision, recall, f01.

    Example
    -------
    .. code-block:: python

        import gretapy as gt
        import pandas as pd

        # Multi-organism GRNs
        grns = {
            "method_a": {
                "hg38": {
                    "PBMC": pd.read_csv("grn_a_pbmc.csv"),
                    "Lung": pd.read_csv("grn_a_lung.csv"),
                },
                "mm10": {
                    "Palate": pd.read_csv("grn_a_palate.csv"),
                },
            },
            "method_b": {
                "hg38": {
                    "PBMC": pd.read_csv("grn_b_pbmc.csv"),
                },
            },
        }
        results = gt.tl.benchmark(grns=grns)

        # With pre-loaded datasets
        results = gt.tl.benchmark(
            grns=grns,
            datasets={"PBMC": mudata_obj, "Lung": mudata_obj2},
        )
    """
    # Validate grns: must be dict[str, dict[str, dict[str, pd.DataFrame]]]
    if not isinstance(grns, dict):
        raise ValueError(f"grns must be dict[str, dict[str, dict[str, DataFrame]]], got {type(grns)}")
    for grn_name, grn_inner in grns.items():
        if not isinstance(grn_inner, dict):
            raise ValueError(
                f"grns['{grn_name}'] must be a dict mapping organism names to dicts, got {type(grn_inner)}"
            )
        for org_key, org_inner in grn_inner.items():
            if not isinstance(org_inner, dict):
                raise ValueError(
                    f"grns['{grn_name}']['{org_key}'] must be a dict mapping dataset names to DataFrames, "
                    f"got {type(org_inner)}"
                )
    grns_dict = grns
    # Extract and validate organisms from grns
    organisms_in_grns = {org for inner in grns_dict.values() for org in inner}
    if not organisms_in_grns:
        raise ValueError("grns is empty or contains no organism keys. Provide at least one organism.")
    if organism is not None:
        _log(
            f"'organism' parameter ('{organism}') is ignored when organisms are encoded in the grns dict. "
            "Organisms are inferred from grns keys.",
            level="warning",
            verbose=verbose,
        )
    for org in organisms_in_grns:
        _check_organism(organism=org)
    # Validate datasets input type
    if not (datasets is None or isinstance(datasets, (list, dict))):
        raise ValueError(f"datasets must be None, list, or dict, got {type(datasets)}")
    datasets_objects = datasets if isinstance(datasets, dict) else None
    # Count pairs for logging
    n_pairs = sum(
        1
        for inner in grns_dict.values()
        for org_inner in inner.values()
        for ds in org_inner
        if datasets is None or (isinstance(datasets, list) and ds in datasets) or (isinstance(datasets, dict) and ds in datasets)
    )
    _log(_SEP, level="info", verbose=verbose)
    _log(
        f"Starting benchmark: {len(grns_dict)} GRN(s), {len(organisms_in_grns)} organism(s), {n_pairs} pair(s)",
        level="info",
        verbose=verbose,
    )
    _log(_SEP, level="info", verbose=verbose)
    t_start_bench = time.time()
    all_results = []
    for grn_name, grn_inner in grns_dict.items():
        for org, org_inner in grn_inner.items():
            # Determine dataset list for this organism
            if datasets is None:
                ds_list = list(org_inner.keys())
            else:
                ds_list = [d for d in (datasets if isinstance(datasets, list) else datasets.keys()) if d in org_inner]
            # Validate metrics per organism
            _check_metrics(organism=org, metrics=metrics)
            for dataset_name in ds_list:
                grn_df = org_inner[dataset_name]
                # Resolve dataset: string name or pre-loaded object
                dataset_arg = datasets_objects[dataset_name] if datasets_objects else dataset_name
                # Resolve terms before eval (new 3-level structure)
                if terms is None:
                    dataset_terms = _check_terms(organism=org, dataset=dataset_name, terms=None)
                else:
                    dataset_terms = terms.get(org, {}).get(dataset_name, {})
                # Warn if no auto-loaded terms for pre-loaded datasets not in config
                if terms is None and datasets_objects is not None and not dataset_terms:
                    _log(
                        f"No terms auto-loaded for dataset '{dataset_name}' (not in config). "
                        "Metrics requiring terms will run unfiltered.",
                        level="warning",
                        verbose=verbose,
                    )
                # Run evaluation
                result = eval_grn_dataset(
                    organism=org,
                    grn=grn_df,
                    dataset=dataset_arg,
                    terms=dataset_terms,
                    metrics=metrics,
                    min_edges=min_edges,
                    grn_name=grn_name,
                    dataset_name=dataset_name,
                    verbose=verbose,
                )
                # Add identifiers
                if not result.empty:
                    result.insert(0, "grn", grn_name)
                    result.insert(1, "organism", org)
                    result.insert(2, "dataset", dataset_name)
                    all_results.append(result)
    elapsed = time.time() - t_start_bench
    _log(_SEP, level="info", verbose=verbose)
    _log(f"Benchmark complete ({len(all_results)} result(s), {elapsed:.1f}s)", level="info", verbose=verbose)
    _log(_SEP, level="info", verbose=verbose)
    if not all_results:
        return pd.DataFrame(columns=["grn", "organism", "dataset", "class", "task", "db", "precision", "recall", "f01"])
    return pd.concat(all_results, ignore_index=True)


def eval_grn_dataset(
    organism: str,
    grn: pd.DataFrame,
    dataset: str | mu.MuData | ad.AnnData,
    terms: dict | None,
    metrics: str | list | None = None,
    min_edges: int = 5,
    grn_name: str | None = None,
    dataset_name: str | None = None,
    verbose: bool = True,
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
    grn_name
        Optional name for the GRN (used in log messages).
    dataset_name
        Optional name for the dataset (used in log messages).
    verbose
        Whether to log progress messages and show progress bars.

    Returns
    -------
    DataFrame with columns: category, metric, db, precision, recall, f01.

    Example
    -------
    .. code-block:: python

        import gretapy as gt
        import pandas as pd

        grn = pd.read_csv("grn.csv")
        results = gt.tl.eval_grn_dataset(
            organism="hg38",
            grn=grn,
            dataset="pbmc10k",
            terms=None,
            metrics=None,
        )
    """
    result_cols = ["class", "task", "db", "precision", "recall", "f01"]
    # Validate inputs
    metrics_list = _check_metrics(organism=organism, metrics=metrics)
    grn = _check_grn(grn=grn)
    if grn.shape[0] < min_edges:
        _log(
            f"GRN has {grn.shape[0]} edges, minimum required is {min_edges}. Returning empty results.",
            level="warning",
            verbose=verbose,
        )
        return pd.DataFrame(columns=result_cols)
    # Resolve dataset_name for logging
    if dataset_name is None and isinstance(dataset, str):
        dataset_name = dataset
    dataset = _check_dataset(organism=organism, dataset=dataset, verbose=verbose)
    terms = _check_terms(organism=organism, dataset=dataset, terms=terms)
    _check_dts_grn(dataset=dataset, grn=grn)
    # Check capabilities
    has_cre = "cre" in grn.columns
    is_mudata = isinstance(dataset, mu.MuData)
    can_run_genomic = has_cre and is_mudata
    if not has_cre:
        _log("GRN does not have 'cre' column. Genomic metrics will be skipped.", level="warning", verbose=verbose)
    if not is_mudata:
        _log(
            "Dataset is AnnData (no ATAC modality). Genomic metrics will be skipped.", level="warning", verbose=verbose
        )
    # Extract data from dataset
    if is_mudata:
        genes, peaks, adata = (
            dataset.mod["rna"].var_names.tolist(),
            dataset.mod["atac"].var_names.tolist(),
            dataset.mod["rna"],
        )
        adata.obs = dataset.obs.copy()
    else:
        genes, peaks, adata = dataset.var_names.tolist(), [], dataset
    # Build log prefix
    prefix = _format_log_prefix(grn_name=grn_name, dataset_name=dataset_name)
    label = _format_label(grn_name=grn_name, dataset_name=dataset_name)
    label_suffix = f": {label}" if label else ""
    _log(_SEP, level="info", verbose=verbose)
    _log(f"Starting evaluation{label_suffix}", level="info", verbose=verbose)
    _log(_SEP, level="info", verbose=verbose)
    t_start_eval = time.time()
    # Evaluate metrics
    results = []
    n_metrics = len(metrics_list)
    for i, db_name in enumerate(metrics_list, 1):
        db_info = DATA[organism]["dbs"].get(db_name)
        if db_info is None:
            continue
        metric_type, category = db_info["metric"], METRIC_CATS.get(db_info["metric"], "Unknown")
        _log(f"{prefix}{category} > {metric_type} > {db_name} [{i} / {n_metrics}]", level="info", verbose=verbose)
        # Handle metrics without file
        if db_info["fname"] is None:
            result = _run_fileless_metric(
                metric_type, db_name, dataset, grn, adata, is_mudata, has_cre, verbose=verbose
            )
            if result is not None:
                results.append([category, metric_type, db_name, *result])
            continue
        # Skip genomic metrics if not possible
        if metric_type in {"TF binding", "CREs", "CRE to gene links"} and not can_run_genomic:
            continue
        # Load database and run metric
        db = read_db(organism=organism, db_name=db_name, verbose=verbose)
        cats = terms.get(db_name, None)
        result = _run_metric(metric_type, db_name, grn, db, genes, peaks, cats, adata, verbose=verbose)
        if result is not None:
            results.append([category, metric_type, db_name, *result])
    elapsed = time.time() - t_start_eval
    _log(_SEP, level="info", verbose=verbose)
    _log(f"Evaluation complete{label_suffix} ({len(results)} metrics, {elapsed:.1f}s)", level="info", verbose=verbose)
    _log(_SEP, level="info", verbose=verbose)
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
    verbose: bool = True,
) -> tuple | None:
    """Run a metric that requires a database file."""
    if metric_type == "Reference GRN":
        return _grn(grn=grn, db=db, genes=genes, verbose=verbose)
    elif metric_type == "TF Markers":
        return _tfm(grn=grn, db=db, genes=genes, cats=cats, verbose=verbose)
    elif metric_type == "TF Pairs":
        return _tfp(grn=grn, db=db, genes=genes, verbose=verbose)
    elif metric_type == "TF Binding":
        return _cre_column(grn=grn, db=db, genes=genes, peaks=peaks, cats=cats, column="source", verbose=verbose)
    elif metric_type == "CREs":
        return _cre(grn=grn, db=db, peaks=peaks, cats=cats, reverse=(db_name == "ENCODE Blacklist"), verbose=verbose)
    elif metric_type == "CRE to gene links":
        return _cre_column(grn=grn, db=db, genes=genes, peaks=peaks, cats=cats, column="target", verbose=verbose)
    elif metric_type == "Gene sets":
        return _gset(adata=adata, grn=grn, db=db, verbose=verbose)
    elif metric_type == "TF Scoring":
        return _tfa(adata=adata, grn=grn, db=db, cats=cats, verbose=verbose)
    elif metric_type == "Perturbation Forecasting":
        return _frc(adata=adata, grn=grn, db=db, cats=cats, verbose=verbose)
    return None


def _run_fileless_metric(
    metric_type: str,
    db_name: str,
    dataset: mu.MuData | ad.AnnData,
    grn: pd.DataFrame,
    adata: ad.AnnData,
    is_mudata: bool,
    has_cre: bool,
    verbose: bool = True,
) -> tuple | None:
    """Run metrics that don't require a database file (Omics, Boolean rules)."""
    if metric_type == "Omics":
        return _run_omics_metric(db_name, dataset, grn, is_mudata, has_cre, verbose=verbose)
    elif metric_type == "Steady state simulation":
        return _sim(adata=adata, grn=grn, verbose=verbose)
    return None


def _run_omics_metric(
    db_name: str,
    dataset: mu.MuData | ad.AnnData,
    grn: pd.DataFrame,
    is_mudata: bool,
    has_cre: bool,
    verbose: bool = True,
) -> tuple | None:
    """Run omics metric based on the specific type."""
    if db_name == "Gene ~ TFs":
        return _omics(
            data=dataset,
            grn=grn,
            col_source="source",
            col_target="target",
            mod_source="rna",
            mod_target="rna",
            verbose=verbose,
        )
    elif db_name == "Gene ~ CREs" and is_mudata and has_cre:
        return _omics(
            data=dataset,
            grn=grn,
            col_source="cre",
            col_target="target",
            mod_source="atac",
            mod_target="rna",
            verbose=verbose,
        )
    elif db_name == "CRE ~ TFs" and is_mudata and has_cre:
        return _omics(
            data=dataset,
            grn=grn,
            col_source="source",
            col_target="cre",
            mod_source="rna",
            mod_target="atac",
            verbose=verbose,
        )
    elif db_name in {"Gene ~ CREs", "CRE ~ TFs"}:
        _log(
            f"Skipping '{db_name}': requires MuData with ATAC and GRN with 'cre' column.",
            level="warning",
            verbose=verbose,
        )
    return None
