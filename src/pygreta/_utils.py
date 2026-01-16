"""Shared utility functions for pygreta."""

import os
import shutil

import pandas as pd
from decoupler._download import _download

from pygreta.config import DATA, METRIC_CATS, PATH_DATA, URL_END, URL_STR


def show_organisms() -> list:
    """
    Show all available organisms.

    Returns
    -------
    List of available organisms.

    Example
    -------
    .. code-block:: python

        import pygreta as pg

        pg.show_organisms()
    """
    return list(DATA.keys())


def show_datasets(organism: str) -> pd.DataFrame:
    """
    Show all available datasets for an organism.

    Parameters
    ----------
    organism
        Which organism to use (e.g., "hg38", "mm10").

    Returns
    -------
    DataFrame with columns: name, pubmed, geo.

    Example
    -------
    .. code-block:: python

        import pygreta as pg

        pg.show_datasets(organism="hg38")
    """
    organisms = show_organisms()
    assert organism in organisms, f"organism={organism} not available ({organisms})"
    datasets = DATA[organism]["dts"]
    rows = []
    for name, info in datasets.items():
        rows.append(
            {
                "name": name,
                "pubmed": info.get("pubmed", ""),
                "geo": info.get("geo", ""),
            }
        )
    return pd.DataFrame(rows)


def show_terms(organism: str) -> pd.DataFrame:
    """
    Show all available terms for filtering databases.

    Parameters
    ----------
    organism
        Which organism to use (e.g., "hg38", "mm10").

    Returns
    -------
    DataFrame listing all available terms per database.

    Example
    -------
    .. code-block:: python

        import pygreta as pg

        pg.show_terms(organism="hg38")
    """
    organisms = show_organisms()
    assert organism in organisms, f"organism={organism} not available ({organisms})"
    fname_terms = DATA[organism]["terms"]
    path_terms = os.path.join(PATH_DATA, fname_terms)
    if not os.path.isfile(path_terms):
        url = URL_STR + fname_terms + URL_END
        data = _download(url, verbose=False)
        data.seek(0)
        with open(path_terms, "wb") as f:
            shutil.copyfileobj(data, f)
    terms_df = pd.read_csv(path_terms)
    return terms_df


def show_metrics(organism: str | None = None) -> pd.DataFrame:
    """
    Show all available metrics.

    Parameters
    ----------
    organism
        Which organism to use. If None, shows metrics for all organisms.

    Returns
    -------
    DataFrame listing all available metrics per organism.

    Example
    -------
    .. code-block:: python

        import pygreta as pg

        pg.show_metrics(organism="hg38")
    """
    assert isinstance(organism, str) or organism is None
    df = []
    for org in DATA:
        org_dbs = DATA[org]["dbs"]
        for db in org_dbs:
            metric = org_dbs[db]["metric"]
            metric_cat = METRIC_CATS[metric]
            df.append([org, metric_cat, metric, db])
    cols = ["organism", "category", "metric", "db"]
    df = pd.DataFrame(df, columns=cols)
    df = df.sort_values(cols)
    if organism:
        organisms = show_organisms()
        assert organism in organisms, f"organism={organism} not available ({organisms})"
        df = df[df["organism"] == organism].drop(columns="organism")
    else:
        df = df.sort_values(cols)
    df = df.reset_index(drop=True)
    return df
