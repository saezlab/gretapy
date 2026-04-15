import gzip
import os
import shutil
import tempfile

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import pyranges as pr
from decoupler._download import _download, _log

from gretapy.config import DATA, PATH_DATA, URL_END, URL_STR

def _download_metrics(
    verbose: bool = False,
):
    fname = 'metrics.csv.gz'
    path_fname = os.path.join(PATH_DATA, fname)
    if not os.path.isfile(path_fname): 
        url = URL_STR + fname + URL_END
        data = _download(url, verbose=verbose)
        data.seek(0)
        with open(path_fname, "wb") as f:
            shutil.copyfileobj(data, f)
        m = f"Metrics saved in {path_fname}"
        _log(m, level="info", verbose=verbose)
    else:
        m = f"Metrics found in {path_fname}"
        _log(m, level="info", verbose=verbose)
    return path_fname


def read_metrics(
    remove_paired: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Read the GRN benchmark metrics table.

    Downloads the metrics file if not already cached locally, then loads and
    preprocesses it. Column names are standardized,
    and optionally paired datasets are removed.

    Parameters
    ----------
    remove_paired : bool
        Whether to remove paired datasets ('Synthetic Pituitary' and
        'Unpaired Pituitary') from the results. Default is True.
    verbose : bool
        Whether to print progress messages. Default is False.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns: name, organism, dataset, task, db, precision,
        recall, and any other columns present in the source file.
    """
    path_fname = _download_metrics(verbose=verbose)
    df = pd.read_csv(path_fname, compression="gzip").dropna()
    df = df.rename(columns={'org': 'organism', 'dts': 'dataset', 'prc': 'precision', 'rcl': 'recall'})
    if remove_paired:
        df = df[~df['dataset'].isin(['Synthetic Pituitary', 'Unpaired Pituitary'])]
    col = []
    for t, d in zip(df['task'], df['db']):
        if d == 'KnockTF':
            if t == 'Perturbation Forecasting':
                col.append('KnockTF (forecasting)')
            elif t == 'TF Scoring':
                col.append('KnockTF (scoring)')
        else:
            col.append(d)
    df['db'] = col
    return df.reset_index(drop=True)
    

def read_imaginary_metrics(
    seed=None,
    remove_paired: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Read benchmark metrics for an imaginary method.

    Calls :func:`read_metrics` and samples one row per unique benchmark
    configuration (class, task, db, organism, dataset), then sets the method
    name to ``'ImaginaryMethod'``. Useful for baseline comparisons or testing
    visualizations.

    Parameters
    ----------
    seed : int or None
        Seed for :func:`numpy.random.default_rng`. Default is None
        (non-deterministic).
    remove_paired : bool
        Passed to :func:`read_metrics`. Default is True.
    verbose : bool
        Passed to :func:`read_metrics`. Default is False.

    Returns
    -------
    pd.DataFrame
        A DataFrame with the same columns as :func:`read_metrics`, with one
        row per unique benchmark configuration and ``name`` set to
        ``'ImaginaryMethod'``.
    """
    rng = np.random.default_rng(seed)
    df = read_metrics(remove_paired=remove_paired, verbose=verbose)
    rows = []
    for _, group in df.groupby(['class', 'task', 'db', 'organism', 'dataset'], sort=False):
        idx = rng.integers(0, len(group))
        rows.append(group.iloc[idx])
    result = pd.DataFrame(rows)
    result['name'] = 'ImaginaryMethod'
    return result.reset_index(drop=True)


def _download_db(
    organism: str,
    db_name: str,
    verbose: bool = False,
) -> str:
    os.makedirs(PATH_DATA, exist_ok=True)
    assert organism in DATA, f"organism={organism} not available:\n{DATA.keys()}"
    assert db_name in DATA[organism]["dbs"], (
        f"db_name={db_name} not available as a database:\n{DATA[organism]['dbs'].keys()}"
    )
    fname = DATA[organism]["dbs"][db_name]["fname"]
    path_fname = os.path.join(PATH_DATA, fname)
    if not os.path.isfile(path_fname):
        url = URL_STR + fname + URL_END
        data = _download(url, verbose=verbose)
        data.seek(0)
        if not '.h5ad' in fname:
            with open(path_fname, "wb") as f:
                shutil.copyfileobj(data, f)
        else:
            with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
                tmp_path = tmp.name
                with gzip.GzipFile(fileobj=data) as gz:
                    shutil.copyfileobj(gz, tmp)
            adata = ad.read_h5ad(tmp_path)
            adata.write(path_fname)
            os.remove(tmp_path)
        m = f"Database {db_name} saved in {path_fname}"
        _log(m, level="info", verbose=verbose)
    else:
        m = f"Database {db_name} found in {path_fname}"
        _log(m, level="info", verbose=verbose)
    return path_fname


def read_db(organism: str, db_name: str, verbose: bool = False) -> pd.DataFrame:
    """
    Read a database file for a given organism.

    Downloads the database if not already cached locally, then reads it
    based on its file format (bed, tsv, csv, or h5ad).

    Parameters
    ----------
    organism : str
        The organism identifier (e.g., 'hg38' for human).
    db_name : str
        The name of the database to read (e.g., 'Promoters', 'CollecTRI').
    verbose : bool
        Whether to print progress messages. Default is False.

    Returns
    -------
    pd.DataFrame | pr.PyRanges | ad.AnnData
        The loaded database. The return type depends on the file format:
        - bed files return PyRanges objects
        - tsv/csv files return pandas DataFrames
        - h5ad files return AnnData objects
    """
    path_fname = _download_db(organism=organism, db_name=db_name, verbose=verbose)
    f_format = os.path.basename(path_fname).replace(".gz", "").split(".")[-1]
    if f_format == "bed":
        db = pr.read_bed(path_fname)
    elif f_format == "tsv":
        db = pd.read_csv(path_fname, sep="\t", compression="gzip", header=None)
    elif f_format == "csv":
        db = pd.read_csv(path_fname, compression="gzip")
    elif f_format == "h5ad":
        db = ad.read_h5ad(path_fname)
    elif f_format == "txt":
        db = pd.read_csv(path_fname, header=None)[0].tolist()
    return db
