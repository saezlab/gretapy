def _trim_grn(grn, ntop=100_000):
    grn["abs_score"] = grn["score"].abs()
    grn = grn.nlargest(ntop, "abs_score", keep="all").reset_index(drop=True)
    grn = grn.drop(columns=["abs_score"])
    return grn
