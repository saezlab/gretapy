import pandas as pd

# Metric classes, in display order.
CLASS_ORDER = ["Predictive", "Genomic", "Literature", "Mechanistic"]


def _class_mean(df):
    """Compute hierarchical mean F0.1 collapsed to a (name x class) matrix."""
    s = df.groupby(["name", "class", "task", "db", "dataset"])["f01"].mean()
    s = s.groupby(["name", "class", "task", "db"]).mean()
    s = s.groupby(["name", "class", "task"]).mean()
    return s.groupby(["name", "class"]).mean().unstack()


def _dataset_class_mean(df):
    """Mean F0.1 collapsed to a (name, dataset) x class matrix, averaging over db then task."""
    s = df.groupby(["name", "dataset", "class", "task", "db"])["f01"].mean()
    s = s.groupby(["name", "dataset", "class", "task"]).mean()  # over db
    s = s.groupby(["name", "dataset", "class"]).mean()  # over task
    return s.unstack()


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


def _weighted_per_dataset(dataset_class_mean, metric_weights):
    """Weighted mean across classes within each dataset.

    For each (name, dataset) the score is ``sum(w * v) / sum(w over present
    classes)``. Returns a Series indexed by (name, dataset).
    """
    weights = pd.Series({c: metric_weights.get(c.lower(), 0.0) for c in dataset_class_mean.columns})
    present_weights = dataset_class_mean.notna().mul(weights, axis=1)
    return dataset_class_mean.mul(weights, axis=1).sum(axis=1) / present_weights.sum(axis=1)


def ranking(df, metric_weights=None):
    """
    Compute the weighted mean F0.1 per method that determines the ranking.

    For each dataset, metric scores are aggregated to the class level (averaging
    over db then task) and combined across metric classes using a weighted mean.
    The final score is the mean of these per-dataset weighted means across
    datasets.

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
        Within each dataset the weighted mean is the sum of ``weight * value``
        divided by the total of the weights of the classes present in that
        dataset; these per-dataset values are then averaged across datasets.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by method name (``name``). The first column
        ``mean_f01`` holds the overall weighted mean F0.1 (averaged across
        datasets), followed by one column per dataset holding that dataset's
        weighted mean F0.1. Rows are sorted by ``mean_f01`` in descending order.
    """
    metric_weights = _check_weights(metric_weights)
    per_dataset = _weighted_per_dataset(_dataset_class_mean(df), metric_weights)
    out = per_dataset.unstack("dataset")
    out.insert(0, "mean_f01", out.mean(axis=1))
    out = out.sort_values("mean_f01", ascending=False)
    out.index.name = "name"
    out.columns.name = None
    return out
