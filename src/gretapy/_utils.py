"""Shared utility functions for gretapy."""

import os
import shutil

import pandas as pd
import pyranges as pr
from decoupler._download import _download

from gretapy.config import DATA, METRIC_CATS, PATH_DATA, URL_END, URL_STR


def show_organisms() -> list:
    """
    Show all available organisms.

    Returns
    -------
    List of available organisms.

    Example
    -------
    .. code-block:: python

        import gretapy as gt

        gt.show_organisms()
    """
    return list(DATA.keys())


def show_datasets(organism: str | None = None) -> pd.DataFrame:
    """
    Show all available datasets for an organism.

    Parameters
    ----------
    organism
        Which organism to use (e.g., "hg38", "mm10"). If None, shows datasets for all organisms.

    Returns
    -------
    DataFrame with columns: name, pubmed, geo.

    Example
    -------
    .. code-block:: python

        import gretapy as gt

        gt.show_datasets(organism="hg38")
    """
    assert isinstance(organism, str) or organism is None
    rows = []
    for org in DATA:
        datasets = DATA[org]["dts"]
        for name, info in datasets.items():
            rows.append(
                {
                    "organism": org,
                    "name": name,
                    "pubmed": info.get("pubmed", ""),
                    "geo": info.get("geo", ""),
                }
            )
    df = pd.DataFrame(rows)
    df = df.sort_values(["organism", "name"])
    if organism:
        organisms = show_organisms()
        assert organism in organisms, f"organism={organism} not available ({organisms})"
        df = df[df["organism"] == organism].drop(columns="organism")
    df = df.reset_index(drop=True)
    return df


def show_terms(organism: str | None = None) -> pd.DataFrame:
    """
    Show all available terms for filtering databases.

    Parameters
    ----------
    organism
        Which organism to use (e.g., "hg38", "mm10"). If None, shows terms for all organisms.

    Returns
    -------
    DataFrame listing all available terms per database.

    Example
    -------
    .. code-block:: python

        import gretapy as gt

        gt.show_terms(organism="hg38")
    """
    assert isinstance(organism, str) or organism is None
    dfs = []
    for org in DATA:
        fname_terms = DATA[org]["terms"]
        path_terms = os.path.join(PATH_DATA, fname_terms)
        if not os.path.isfile(path_terms):
            url = URL_STR + fname_terms + URL_END
            data = _download(url, verbose=False)
            data.seek(0)
            with open(path_terms, "wb") as f:
                shutil.copyfileobj(data, f)
        terms_df = pd.read_csv(path_terms)
        terms_df.insert(0, "organism", org)
        dfs.append(terms_df)
    df = pd.concat(dfs)
    if organism:
        organisms = show_organisms()
        assert organism in organisms, f"organism={organism} not available ({organisms})"
        df = df[df["organism"] == organism].drop(columns="organism")
    df = df.reset_index(drop=True)
    df["db_name"] = df["db_name"].str.replace("Human Protein Atlas (HPA)", "HPA", regex=False)
    return df


def show_genome_annotation(organism: str) -> pr.PyRanges:
    """
    Show genome annotation for an organism.

    Parameters
    ----------
    organism
        Which organism to use (e.g., "hg38", "mm10").

    Returns
    -------
    PyRanges object with genome annotation.

    Example
    -------
    .. code-block:: python

        import gretapy as gt

        gt.show_genome_annotation(organism="hg38")
    """
    organisms = show_organisms()
    assert organism in organisms, f"organism={organism} not available ({organisms})"
    fname_gann = DATA[organism]["gann"]
    path_gann = os.path.join(PATH_DATA, fname_gann)
    if not os.path.isfile(path_gann):
        url = URL_STR + fname_gann + URL_END
        data = _download(url, verbose=False)
        data.seek(0)
        with open(path_gann, "wb") as f:
            shutil.copyfileobj(data, f)
    gann_pr = pr.read_bed(path_gann)
    return gann_pr


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

        import gretapy as gt

        gt.show_metrics(organism="hg38")
    """
    assert isinstance(organism, str) or organism is None
    df = []
    for org in DATA:
        org_dbs = DATA[org]["dbs"]
        for db in org_dbs:
            metric = org_dbs[db]["metric"]
            if metric is not None:
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
