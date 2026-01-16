import os
import shutil

import anndata as ad
import decoupler as dc
import pandas as pd
import pyranges as pr
from decoupler._download import _download, _log

from pygreta.config import DATA, PATH_DATA, URL_END, URL_STR


def _download_db(
    organism: str,
    db_name: str,
    verbose: bool = False,
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


def _read_db(organism: str, db_name: str, verbose: bool = False) -> pd.DataFrame:
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
    return db
