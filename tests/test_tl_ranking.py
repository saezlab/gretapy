"""Tests for gretapy.tl._ranking module."""

import numpy as np
import pandas as pd
import pytest

from gretapy.tl._ranking import _check_weights, _class_mean, ranking


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

    def test_single_mean_f01_column(self, sample_ranking_df):
        out = ranking(sample_ranking_df)
        assert list(out.columns) == ["mean_f01"]

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
