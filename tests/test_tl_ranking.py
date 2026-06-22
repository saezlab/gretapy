"""Tests for gretapy.tl._ranking module."""

import numpy as np
import pandas as pd
import pytest

from gretapy.tl._ranking import _check_weights, _class_mean, _dataset_class_mean, ranking


@pytest.fixture
def sample_ranking_df():
    """Sample metrics DataFrame with all four classes for three methods."""
    np.random.seed(42)
    methods = ["MethodA", "MethodB", "MethodC"]
    rows = (
        [
            {
                "name": m,
                "class": "Predictive",
                "task": "Gene Sets",
                "db": "Hallmarks",
                "dataset": "DatasetX",
                "f01": np.random.uniform(0, 1),
            }
            for m in methods
        ]
        + [
            {
                "name": m,
                "class": "Genomic",
                "task": "TF Binding",
                "db": "ChIP-Atlas",
                "dataset": "DatasetX",
                "f01": np.random.uniform(0, 1),
            }
            for m in methods
        ]
        + [
            {
                "name": m,
                "class": "Literature",
                "task": "TF Markers",
                "db": "HPA",
                "dataset": "DatasetX",
                "f01": np.random.uniform(0, 1),
            }
            for m in methods
        ]
        + [
            {
                "name": m,
                "class": "Mechanistic",
                "task": "TF Scoring",
                "db": "KnockTF",
                "dataset": "DatasetX",
                "f01": np.random.uniform(0, 1),
            }
            for m in methods
        ]
    )
    return pd.DataFrame(rows)


@pytest.fixture
def sample_ranking_df_with_nan():
    """DataFrame where one method has no data for a class → NaN in class_mean."""
    rows = [
        {
            "name": "MethodA",
            "class": "Predictive",
            "task": "Gene Sets",
            "db": "Hallmarks",
            "dataset": "DatasetX",
            "f01": 0.5,
        },
        {
            "name": "MethodA",
            "class": "Literature",
            "task": "TF Markers",
            "db": "HPA",
            "dataset": "DatasetX",
            "f01": 0.6,
        },
        # MethodB only has Predictive data → NaN for Literature class
        {
            "name": "MethodB",
            "class": "Predictive",
            "task": "Gene Sets",
            "db": "Hallmarks",
            "dataset": "DatasetX",
            "f01": 0.4,
        },
    ]
    return pd.DataFrame(rows)


class TestClassMean:
    """Tests for the _class_mean helper."""

    def test_returns_dataframe(self, sample_ranking_df):
        cm = _class_mean(sample_ranking_df)
        assert isinstance(cm, pd.DataFrame)

    def test_methods_in_index(self, sample_ranking_df):
        cm = _class_mean(sample_ranking_df)
        for method in ["MethodA", "MethodB", "MethodC"]:
            assert method in cm.index

    def test_classes_in_columns(self, sample_ranking_df):
        cm = _class_mean(sample_ranking_df)
        for cls in ["Predictive", "Genomic", "Literature", "Mechanistic"]:
            assert cls in cm.columns


class TestDatasetClassMean:
    """Tests for the _dataset_class_mean helper."""

    def test_index_has_name_and_dataset(self, sample_ranking_df):
        dcm = _dataset_class_mean(sample_ranking_df)
        assert list(dcm.index.names) == ["name", "dataset"]

    def test_classes_in_columns(self, sample_ranking_df):
        dcm = _dataset_class_mean(sample_ranking_df)
        for cls in ["Predictive", "Genomic", "Literature", "Mechanistic"]:
            assert cls in dcm.columns


class TestCheckWeights:
    """Tests for the _check_weights validation helper."""

    def test_default_when_none(self):
        assert _check_weights(None) == {"predictive": 1, "genomic": 1, "literature": 1, "mechanistic": 1}

    def test_subset_allowed(self):
        assert _check_weights({"predictive": 2}) == {"predictive": 2}

    def test_invalid_key_raises(self):
        with pytest.raises(AssertionError, match="subset"):
            _check_weights({"foo": 1})

    def test_negative_weight_raises(self):
        with pytest.raises(AssertionError, match="non-negative"):
            _check_weights({"predictive": -1})

    def test_zero_weight_allowed(self):
        assert _check_weights({"predictive": 0}) == {"predictive": 0}


class TestRanking:
    """Tests for the public ranking function."""

    def test_returns_dataframe(self, sample_ranking_df):
        out = ranking(sample_ranking_df)
        assert isinstance(out, pd.DataFrame)

    def test_mean_f01_is_first_column(self, sample_ranking_df):
        out = ranking(sample_ranking_df)
        assert out.columns[0] == "mean_f01"

    def test_per_dataset_columns_present(self, sample_ranking_df):
        out = ranking(sample_ranking_df)
        assert list(out.columns) == ["mean_f01", "DatasetX"]

    def test_mean_f01_is_mean_of_dataset_columns(self, sample_ranking_df):
        out = ranking(sample_ranking_df)
        dataset_cols = out.columns.drop("mean_f01")
        assert np.allclose(out["mean_f01"].values, out[dataset_cols].mean(axis=1).values)

    def test_index_named_name(self, sample_ranking_df):
        out = ranking(sample_ranking_df)
        assert out.index.name == "name"
        for method in ["MethodA", "MethodB", "MethodC"]:
            assert method in out.index

    def test_sorted_descending(self, sample_ranking_df):
        out = ranking(sample_ranking_df)
        assert out["mean_f01"].is_monotonic_decreasing

    def test_default_equals_plain_class_mean(self, sample_ranking_df):
        """Equal default weights reproduce the plain mean across classes."""
        out = ranking(sample_ranking_df)
        plain = _class_mean(sample_ranking_df).mean(axis=1)
        assert np.allclose(out["mean_f01"].values, plain.loc[out.index].values)

    def test_single_class_weight(self, sample_ranking_df):
        """A single weight key reduces the score to that class only."""
        out = ranking(sample_ranking_df, {"predictive": 2})
        cm = _class_mean(sample_ranking_df)
        assert np.allclose(out["mean_f01"].loc[cm.index].values, cm["Predictive"].loc[cm.index].values)

    def test_zero_weight_excludes_class(self, sample_ranking_df):
        """Weight 0 for a class is equivalent to dropping it from the dict."""
        with_zero = ranking(sample_ranking_df, {"predictive": 1, "genomic": 1, "literature": 0, "mechanistic": 0})
        without = ranking(sample_ranking_df, {"predictive": 1, "genomic": 1})
        pd.testing.assert_frame_equal(with_zero, without)

    def test_missing_class_treated_as_zero(self, sample_ranking_df):
        """Omitting a class from the dict gives the same result as weight 0."""
        omitted = ranking(sample_ranking_df, {"predictive": 1})
        explicit = ranking(sample_ranking_df, {"predictive": 1, "genomic": 0, "literature": 0, "mechanistic": 0})
        pd.testing.assert_frame_equal(omitted, explicit)

    def test_nan_class_divides_by_present_weights(self, sample_ranking_df_with_nan):
        """A method missing a class is averaged only over its present classes."""
        out = ranking(sample_ranking_df_with_nan)
        # MethodB only has Predictive (0.4) → weighted mean over present classes is 0.4
        assert np.isclose(out["mean_f01"].loc["MethodB"], 0.4)
        # MethodA has Predictive (0.5) and Literature (0.6) → mean 0.55
        assert np.isclose(out["mean_f01"].loc["MethodA"], 0.55)

    def test_invalid_key_raises(self, sample_ranking_df):
        with pytest.raises(AssertionError, match="subset"):
            ranking(sample_ranking_df, {"foo": 1})

    def test_negative_weight_raises(self, sample_ranking_df):
        with pytest.raises(AssertionError, match="non-negative"):
            ranking(sample_ranking_df, {"predictive": -1})

    def test_weighting_is_per_dataset_then_mean_across_datasets(self):
        """Weighting happens per dataset, then averaged across datasets.

        MethodM has Predictive in two datasets but Mechanistic in only one.
        Per-dataset weighted means (equal weights): D1 = (0.8 + 0.4) / 2 = 0.6,
        D2 = 0.2 / 1 = 0.2. Final = mean(0.6, 0.2) = 0.4. This differs from
        weighting the class means first ((mean(0.8, 0.2) + 0.4) / 2 = 0.45).
        """
        df = pd.DataFrame(
            [
                {
                    "name": "MethodM",
                    "class": "Predictive",
                    "task": "Gene Sets",
                    "db": "Hallmarks",
                    "dataset": "D1",
                    "f01": 0.8,
                },
                {
                    "name": "MethodM",
                    "class": "Mechanistic",
                    "task": "TF Scoring",
                    "db": "KnockTF",
                    "dataset": "D1",
                    "f01": 0.4,
                },
                {
                    "name": "MethodM",
                    "class": "Predictive",
                    "task": "Gene Sets",
                    "db": "Hallmarks",
                    "dataset": "D2",
                    "f01": 0.2,
                },
            ]
        )
        out = ranking(df)
        assert np.isclose(out["mean_f01"].loc["MethodM"], 0.4)
        # Per-dataset columns hold each dataset's weighted mean.
        assert np.isclose(out["D1"].loc["MethodM"], 0.6)
        assert np.isclose(out["D2"].loc["MethodM"], 0.2)
        # Confirm it diverges from the old class-mean-then-weight ordering (0.45).
        assert not np.isclose(out["mean_f01"].loc["MethodM"], 0.45)
