import os
import shutil

import decoupler as dc
import pandas as pd
from decoupler._download import _download, _log

from pygreta.config import DATA, PATH_DATA, URL_END, URL_STR


def show_organisms() -> None:
    """
    Show all available organisms.

    Returns
    -------
    List of available organisms.

    Example
    -------
    .. code-block:: python

        import pygreta as pg

        pg.show_organisms(organism="hg38")
    """
    return list(DATA.keys())


def _download_data(
    organism: str,
    db_name: str,
    verbose: bool = False,
    retries: int = 5,
    wait_time: int = 20,
) -> str:
    os.makedirs(PATH_DATA, exist_ok=True)
    assert organism in DATA, f"organism={organism} not available:\n{DATA.keys()}"
    assert db_name in DATA[organism]["dbs"], f"db_name={db_name} not available as a database:\n{DATA[organism]['dbs']}"
    fname = DATA[organism]["dbs"][db_name]["fname"]
    path_fname = os.path.join(PATH_DATA, fname)
    if not os.path.isfile(path_fname):
        if fname != "hg38_prt_knocktf.h5ad":
            url = URL_STR + fname + URL_END
            data = _download(url, verbose=verbose)
            data.seek(0)  # Move pointer to beginning
            with open(path_fname, "wb") as f:
                shutil.copyfileobj(data, f)
            m = f"Database {db_name} saved in {path_fname}"
            _log(m, level="info", verbose=verbose)
        else:
            adata = dc.ds.knocktf(thr_fc=100_000, verbose=verbose)  # Do not filter here
            adata.write(path_fname)
    else:
        m = f"Database {db_name} found in {path_fname}"
        _log(m, level="info", verbose=verbose)
    return path_fname


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

        pg.ds.show_terms(organism="hg38")
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
