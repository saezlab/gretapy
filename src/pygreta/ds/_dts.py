import os
import shutil

import mudata as mu
import pandas as pd
from decoupler._download import _download, _log

from pygreta.config import DATA, PATH_DATA, URL_END, URL_STR


def _download_dts(
    organism: str,
    dts_name: str,
    verbose: bool = False,
) -> str:
    os.makedirs(PATH_DATA, exist_ok=True)
    assert organism in DATA, f"organism={organism} not available:\n{DATA.keys()}"
    assert dts_name in DATA[organism]["dts"], (
        f"dts_name={dts_name} not available as a dataset:\n{DATA[organism]['dts']}"
    )
    fname = DATA[organism]["dts"][dts_name]["fname"]
    path_fname = os.path.join(PATH_DATA, fname)
    if not os.path.isfile(path_fname):
        url = URL_STR + fname + URL_END
        data = _download(url, verbose=verbose)
        data.seek(0)  # Move pointer to beginning
        with open(path_fname, "wb") as f:
            shutil.copyfileobj(data, f)
        m = f"Dataset {dts_name} saved in {path_fname}"
        _log(m, level="info", verbose=verbose)
    else:
        m = f"Database {dts_name} found in {path_fname}"
        _log(m, level="info", verbose=verbose)
    return path_fname


def _read_dts(organism: str, dts_name: str, verbose: bool = False) -> pd.DataFrame:
    path_fname = _download_dts(organism=organism, dts_name=dts_name, verbose=verbose)
    dts = mu.read(path_fname)
    return dts
