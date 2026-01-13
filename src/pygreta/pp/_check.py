import anndata as ad
import mudata as mu
import pandas as pd

from pygreta.config import DATA
from pygreta.ds._utils import show_terms
from pygreta.tl._utils import show_metrics


def _check_metrics(
    organism: str,
    metrics: str | list | None = None,
) -> list:
    """
    Validate and expand metric selection to database names.

    Parameters
    ----------
    organism
        Which organism to use.
    metrics
        Metric(s) to evaluate. Can be category, metric type, or database name.
        If None, all metrics are selected.

    Returns
    -------
    List of database names to evaluate.
    """
    assert isinstance(metrics, (str, list)) or metrics is None, (
        f"metrics must be str, list or None, got {type(metrics)}"
    )
    df_metrics = show_metrics(organism=organism)
    if isinstance(metrics, str):
        metrics = [metrics]
    if metrics is not None:
        # Expand if necessary or raise error
        tmp_metrics = []
        cats = df_metrics["category"].tolist()
        sub_cats = df_metrics["metric"].tolist()
        names = df_metrics["db"].tolist()
        for metric in metrics:
            if metric in cats:
                tmp_metrics.extend(df_metrics[df_metrics["category"] == metric]["db"].to_list())
            elif metric in sub_cats:
                tmp_metrics.extend(df_metrics[df_metrics["metric"] == metric]["db"].to_list())
            elif metric in names:
                tmp_metrics.extend(df_metrics[df_metrics["db"] == metric]["db"].to_list())
            else:
                raise ValueError(
                    f"Invalid metric or database: '{metric}'. View available options: pygreta.tl.show_metrics()"
                )
        metrics = list(set(tmp_metrics))
    else:
        # Do all metrics and dbs
        metrics = df_metrics["db"].tolist()
    return metrics


def _check_grn(
    grn: pd.DataFrame,
) -> pd.DataFrame:
    """
    Validate GRN DataFrame structure.

    Parameters
    ----------
    grn
        GRN DataFrame with at least "source" and "target" columns.

    Returns
    -------
    Validated and deduplicated GRN DataFrame.
    """
    assert isinstance(grn, pd.DataFrame), f"grn must be pd.DataFrame, got {type(grn)}"
    assert {"source", "target"}.issubset(grn.columns), (
        f'grn must contain "source" and "target" column names, got {grn.columns}'
    )
    if "cre" in grn.columns:
        pattern = r"^chr[0-9a-zA-Z]+-\d+-\d+$"
        assert grn["cre"].str.match(pattern).all(), (
            'Region names in column "cre" are not well formatted. Expected format: chrX-start-end (e.g., chr7-87600040-87600540)'
        )
        unique_cols = ["source", "cre", "target"]
    else:
        unique_cols = ["source", "target"]
    # Ensure uniqueness
    if "score" in grn.columns:
        grn = grn.groupby(unique_cols, as_index=False)["score"].mean()
    else:
        grn = grn.drop_duplicates(unique_cols)
    return grn


def _check_dataset(
    organism: str,
    dataset: str | mu.MuData | ad.AnnData,
) -> mu.MuData | ad.AnnData:
    """
    Validate dataset input.

    Parameters
    ----------
    organism
        Which organism to use.
    dataset
        Dataset name (str) or loaded data object (MuData/AnnData).

    Returns
    -------
    Loaded dataset as MuData or AnnData.
    """
    if isinstance(dataset, str):
        assert dataset in DATA[organism]["dts"], (
            f'Dataset "{dataset}" not found in config. Run pygreta.ds.show_datasets() to see available datasets'
        )
        # TODO: download and return dataset
        # dataset = _download_dataset(name=dataset)
        raise NotImplementedError(
            "Downloading datasets by name is not yet implemented. Please load the dataset manually."
        )
    elif isinstance(dataset, mu.MuData):
        assert {"rna", "atac"}.issubset(dataset.mod), 'Modalities "rna" and "atac" missing in dataset.mod'
        assert "celltype" in dataset.obs.columns, 'Column "celltype" not found in dataset.obs'
    elif isinstance(dataset, ad.AnnData):
        assert "celltype" in dataset.obs.columns, 'Column "celltype" not found in dataset.obs'
    else:
        raise ValueError(f"Invalid type for dataset={type(dataset)}. Must be str, mudata.MuData or anndata.AnnData")
    return dataset


def _check_dts_grn(
    dataset: mu.MuData | ad.AnnData,
    grn: pd.DataFrame,
) -> None:
    """
    Validate gene overlap between dataset and GRN.

    Parameters
    ----------
    dataset
        Loaded dataset (MuData or AnnData).
    grn
        GRN DataFrame.

    Raises
    ------
    AssertionError
        If genes in GRN do not exist in dataset.
    """
    genes_grn = set(grn["source"]) | set(grn["target"])
    if isinstance(dataset, ad.AnnData):
        genes_dts = set(dataset.var_names)
    elif isinstance(dataset, mu.MuData):
        genes_dts = set(dataset.mod["rna"].var_names)
    genes_dif = list(genes_grn - genes_dts)
    n_diff = len(genes_dif)
    assert n_diff == 0, f"{n_diff} genes from grn do not exist in dataset.var_names: {genes_dif[:5]}"


def _check_terms(
    organism: str,
    dataset: str | mu.MuData | ad.AnnData,
    terms: dict | None,
) -> dict:
    """
    Validate terms against available options.

    Parameters
    ----------
    organism
        Which organism to use.
    dataset
        Dataset name (str) or loaded data object.
    terms
        Dictionary mapping database names to lists of terms for filtering.
        If None and dataset is str, terms are loaded from config.

    Returns
    -------
    Validated terms dictionary.
    """
    if isinstance(dataset, str):
        # Auto-load terms from config
        if terms is None:
            dts_config = DATA[organism]["dts"].get(dataset, {})
            terms = dts_config.get("terms", {})
    else:
        # terms cannot be None for non-string datasets
        if terms is None:
            raise ValueError("terms cannot be None when dataset is not a string. Please provide a terms dictionary.")

    # Validate terms against available options
    terms_df = show_terms(organism=organism)
    for db in terms:
        assert db in DATA[organism]["dbs"], (
            f'db="{db}" not found in databases. View available options: pygreta.tl.show_metrics()'
        )
        og_db_terms = set(terms_df[terms_df["db_name"] == db]["term"])
        db_terms = set(terms[db])
        diff_terms = list(db_terms - og_db_terms)
        n_diff = len(diff_terms)
        assert n_diff == 0, (
            f"{n_diff} terms do not exist in db={db}: {diff_terms[:5]} View available options: pygreta.ds.show_terms()"
        )
    return terms
