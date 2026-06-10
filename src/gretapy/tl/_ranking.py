import pandas as pd

# Metric classes, in display order.
CLASS_ORDER = ["Predictive", "Genomic", "Literature", "Mechanistic"]


def _class_mean(df):
    """Compute hierarchical mean F0.1 collapsed to a (name x class) matrix."""
    s = df.groupby(["name", "class", "task", "db", "dataset"])["f01"].mean()
    s = s.groupby(["name", "class", "task", "db"]).mean()
    s = s.groupby(["name", "class", "task"]).mean()
    return s.groupby(["name", "class"]).mean().unstack()


def _check_weights(metric_weights):
    """Validate ``metric_weights`` and fill the default if None."""
    if metric_weights is None:
        metric_weights = {"predictive": 1, "genomic": 1, "literature": 1, "mechanistic": 1}
    valid = {c.lower() for c in CLASS_ORDER}
    assert set(metric_weights).issubset(valid), (
        f"metric_weights keys must be a subset of {sorted(valid)}, got {sorted(set(metric_weights) - valid)}"
    )
    assert all(w >= 0 for w in metric_weights.values()), "metric_weights must be non-negative"
    return metric_weights


def _weighted_overall(class_mean, metric_weights):
    """Weighted mean across classes: sum(w * v) / sum(w over present classes)."""
    weights = pd.Series({c: metric_weights.get(c.lower(), 0.0) for c in class_mean.columns})
    present_weights = class_mean.notna().mul(weights, axis=1)
    return class_mean.mul(weights, axis=1).sum(axis=1) / present_weights.sum(axis=1)


def ranking(df, metric_weights=None):
    """
    Compute the weighted mean F0.1 per method that determines the ranking.

    Metrics are aggregated hierarchically (dataset -> db -> task -> class) and
    then combined across metric classes using a weighted mean.

    Parameters
    ----------
    df : pandas.DataFrame
        Metrics dataframe with columns: class, task, db, dataset, name, f01.
    metric_weights : dict or None
        Weight for each metric class. Keys are class names (``predictive``,
        ``genomic``, ``literature``, ``mechanistic``); all must be non-negative.
        Any metric class missing from the dictionary is treated as having a
        weight of 0 (i.e. excluded from the ranking).
        Defaults to ``dict(predictive=1, genomic=1, literature=1, mechanistic=1)``.
        The overall score is the weighted mean across classes, i.e. the sum of
        ``weight * value`` divided by the total of the weights of the classes
        present for that method.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by method name (``name``) with a single column
        ``mean_f01`` holding the weighted mean F0.1, sorted in descending order.
    """
    metric_weights = _check_weights(metric_weights)
    overall = _weighted_overall(_class_mean(df), metric_weights)
    out = overall.sort_values(ascending=False).to_frame("mean_f01")
    out.index.name = "name"
    return out
