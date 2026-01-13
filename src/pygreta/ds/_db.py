import os

import anndata as ad
import pandas as pd
import pyranges as pr

from pygreta.ds._utils import _download_data


def _read_db(organism: str, db_name: str, verbose: bool = False) -> pd.DataFrame:
    path_fname = _download_data(organism=organism, db_name=db_name, verbose=verbose)
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
