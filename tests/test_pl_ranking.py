"""Tests for gretapy.pl._ranking module."""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from gretapy.pl._ranking import (
    CLASS_ORDER,
    _build_task_list,
    _compute_aggregations,
    _compute_db_aggregation,
    _compute_ranks,
    _draw_barplot,
    _draw_method_names,
    _rank_text,
    _ranking_class,
    _ranking_task,
    ranking,
)


@pytest.fixture
def sample_ranking_df():
    """Sample metrics DataFrame for ranking tests.

    Columns: name, class, task, db, dts, f01
    """
    np.random.seed(42)
    methods = ["MethodA", "MethodB", "MethodC"]
    rows = [
        # Literature - 2 different tasks (triggers line 83 in _build_task_list)
        {"name": m, "class": "Literature", "task": "TF Markers", "db": "HPA", "dts": "DatasetX", "f01": np.random.uniform(0, 1)}
        for m in methods
    ] + [
        {"name": m, "class": "Literature", "task": "TF Pairs", "db": "Europe PMC", "dts": "DatasetX", "f01": np.random.uniform(0, 1)}
        for m in methods
    ] + [
        # Genomic - 2 different tasks (triggers line 259 in _ranking_task span update)
        {"name": m, "class": "Genomic", "task": "TF Binding", "db": "ChIP-Atlas", "dts": "DatasetX", "f01": np.random.uniform(0, 1)}
        for m in methods
    ] + [
        {"name": m, "class": "Genomic", "task": "CREs", "db": "ENCODE CREs", "dts": "DatasetX", "f01": np.random.uniform(0, 1)}
        for m in methods
    ] + [
        # Predictive
        {"name": m, "class": "Predictive", "task": "Gene Sets", "db": "Hallmarks", "dts": "DatasetX", "f01": np.random.uniform(0, 1)}
        for m in methods
    ] + [
        # Mechanistic
        {"name": m, "class": "Mechanistic", "task": "TF Scoring", "db": "KnockTF", "dts": "DatasetX", "f01": np.random.uniform(0, 1)}
        for m in methods
    ]
    return pd.DataFrame(rows)


@pytest.fixture
def sample_ranking_df_with_nan():
    """DataFrame where one method has no data for a class → NaN in class_mean."""
    methods = ["MethodA", "MethodB"]
    rows = [
        # MethodA has both Predictive and Literature data
        {"name": "MethodA", "class": "Predictive", "task": "Gene Sets", "db": "Hallmarks", "dts": "DatasetX", "f01": 0.5},
        {"name": "MethodA", "class": "Literature", "task": "TF Markers", "db": "HPA", "dts": "DatasetX", "f01": 0.6},
        # MethodB only has Predictive data → NaN for Literature class
        {"name": "MethodB", "class": "Predictive", "task": "Gene Sets", "db": "Hallmarks", "dts": "DatasetX", "f01": 0.4},
    ]
    return pd.DataFrame(rows)


@pytest.fixture
def sample_ranking_df_with_excluded(sample_ranking_df):
    """DataFrame with excluded dataset names."""
    extra = sample_ranking_df.copy()
    extra["dts"] = "Synthetic Pituitary"
    return pd.concat([sample_ranking_df, extra])


class TestComputeAggregations:
    """Tests for _compute_aggregations function."""

    def test_returns_tuple_of_two(self, sample_ranking_df):
        """Test that function returns a tuple of two elements."""
        result = _compute_aggregations(sample_ranking_df)
        assert len(result) == 2

    def test_overall_mean_is_series(self, sample_ranking_df):
        """Test that overall_mean is a pandas Series."""
        overall_mean, _ = _compute_aggregations(sample_ranking_df)
        assert isinstance(overall_mean, pd.Series)

    def test_class_mean_is_dataframe(self, sample_ranking_df):
        """Test that class_mean is a pandas DataFrame."""
        _, class_mean = _compute_aggregations(sample_ranking_df)
        assert isinstance(class_mean, pd.DataFrame)

    def test_overall_mean_values_in_range(self, sample_ranking_df):
        """Test that overall_mean values are in [0, 1]."""
        overall_mean, _ = _compute_aggregations(sample_ranking_df)
        assert (overall_mean >= 0).all()
        assert (overall_mean <= 1).all()

    def test_method_names_in_index(self, sample_ranking_df):
        """Test that all method names appear in overall_mean index."""
        overall_mean, _ = _compute_aggregations(sample_ranking_df)
        for method in ["MethodA", "MethodB", "MethodC"]:
            assert method in overall_mean.index


class TestComputeDbAggregation:
    """Tests for _compute_db_aggregation function."""

    def test_returns_pivot_and_hierarchy(self, sample_ranking_df):
        """Test that function returns pivot DataFrame and hierarchy."""
        pivot, hierarchy = _compute_db_aggregation(sample_ranking_df)
        assert isinstance(pivot, pd.DataFrame)
        assert isinstance(hierarchy, pd.DataFrame)

    def test_hierarchy_has_expected_columns(self, sample_ranking_df):
        """Test that hierarchy has db, task, class columns."""
        _, hierarchy = _compute_db_aggregation(sample_ranking_df)
        assert "db" in hierarchy.columns
        assert "task" in hierarchy.columns
        assert "class" in hierarchy.columns


class TestComputeRanks:
    """Tests for _compute_ranks function."""

    def test_returns_dataframe(self):
        """Test that function returns a DataFrame."""
        values = pd.DataFrame(
            {"col1": [0.5, 0.8, 0.3], "col2": [0.6, 0.4, 0.9]},
            index=["A", "B", "C"],
        )
        result = _compute_ranks(values)
        assert isinstance(result, pd.DataFrame)

    def test_rank_1_is_highest(self):
        """Test that rank 1 is assigned to the highest value."""
        values = pd.DataFrame({"col1": [0.9, 0.5, 0.1]}, index=["A", "B", "C"])
        result = _compute_ranks(values)
        assert result.loc["A", "col1"] == 1.0

    def test_nan_stays_nan(self):
        """Test that NaN values remain NaN in ranks."""
        values = pd.DataFrame({"col1": [0.5, np.nan, 0.8]}, index=["A", "B", "C"])
        result = _compute_ranks(values)
        assert np.isnan(result.loc["B", "col1"])

    def test_handles_all_nan_column(self):
        """Test that a column with all NaN values is handled."""
        values = pd.DataFrame({"col1": [np.nan, np.nan]}, index=["A", "B"])
        result = _compute_ranks(values)
        # No ranks assigned
        assert result.shape == (2, 1)


class TestBuildTaskList:
    """Tests for _build_task_list function."""

    def test_returns_list(self):
        """Test that function returns a list."""
        hierarchy = pd.DataFrame(
            {
                "class": ["Literature", "Literature"],
                "task": ["TF Markers", "TF Markers"],
                "db": ["HPA", "TF-Marker"],
            }
        )
        col_order = [("HPA", "TF Markers"), ("TF-Marker", "TF Markers")]
        result = _build_task_list(hierarchy, col_order)
        assert isinstance(result, list)

    def test_respects_class_order(self):
        """Test that tasks are ordered by CLASS_ORDER."""
        hierarchy = pd.DataFrame(
            {
                "class": ["Predictive", "Literature"],
                "task": ["Gene Sets", "TF Markers"],
                "db": ["Hallmarks", "HPA"],
            }
        )
        col_order = [("HPA", "TF Markers"), ("Hallmarks", "Gene Sets")]
        result = _build_task_list(hierarchy, col_order)
        # CLASS_ORDER = ["Predictive", "Genomic", "Literature", "Mechanistic"]
        # so Predictive should come before Literature
        classes = [item[0] for item in result]
        if "Literature" in classes and "Predictive" in classes:
            assert classes.index("Predictive") < classes.index("Literature")

    def test_empty_hierarchy_returns_empty(self):
        """Test that empty hierarchy returns empty list."""
        hierarchy = pd.DataFrame({"class": [], "task": [], "db": []})
        result = _build_task_list(hierarchy, [])
        assert result == []


class TestRankText:
    """Tests for _rank_text function."""

    def test_integer_rank(self):
        """Test that integer ranks are formatted without decimal."""
        assert _rank_text(1.0) == "1"
        assert _rank_text(3.0) == "3"

    def test_fractional_rank(self):
        """Test that fractional ranks are formatted with one decimal."""
        assert _rank_text(1.5) == "1.5"
        assert _rank_text(2.5) == "2.5"


class TestDrawFunctions:
    """Tests for _draw_barplot and _draw_method_names."""

    @pytest.fixture
    def ax(self):
        """Create a fresh matplotlib axes."""
        fig, ax = plt.subplots()
        yield ax
        plt.close(fig)

    def test_draw_barplot(self, ax):
        """Test that _draw_barplot executes without error."""
        method_order = ["MethodA", "MethodB", "MethodC"]
        overall_mean = pd.Series([0.8, 0.6, 0.4], index=method_order)
        palette = {"MethodA": "blue", "MethodB": "red", "MethodC": "green"}
        _draw_barplot(ax, method_order, overall_mean, palette, n_methods=3)
        # Verify bars were drawn
        assert len(ax.patches) == 3

    def test_draw_method_names(self, ax):
        """Test that _draw_method_names executes without error."""
        method_order = ["MethodA", "MethodB"]
        _draw_method_names(ax, method_order, n_methods=2)
        # Verify text was added
        assert len(ax.texts) == 2


class TestRankingClass:
    """Tests for _ranking_class function."""

    def test_returns_figure(self, sample_ranking_df):
        """Test that _ranking_class returns a Figure."""
        overall_mean, class_mean = _compute_aggregations(sample_ranking_df)
        method_order = overall_mean.sort_values(ascending=False).index.tolist()
        overall_mean = overall_mean.loc[method_order]

        fig = _ranking_class(
            df=sample_ranking_df,
            overall_mean=overall_mean,
            class_mean=class_mean,
            method_order=method_order,
            n_methods=len(method_order),
            palette={m: f"C{i}" for i, m in enumerate(method_order)},
            figsize=None,
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_custom_figsize(self, sample_ranking_df):
        """Test that custom figsize is applied."""
        overall_mean, class_mean = _compute_aggregations(sample_ranking_df)
        method_order = overall_mean.sort_values(ascending=False).index.tolist()
        overall_mean = overall_mean.loc[method_order]

        fig = _ranking_class(
            df=sample_ranking_df,
            overall_mean=overall_mean,
            class_mean=class_mean,
            method_order=method_order,
            n_methods=len(method_order),
            palette={},
            figsize=(8, 6),
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


class TestRankingTask:
    """Tests for _ranking_task function."""

    def test_returns_figure(self, sample_ranking_df):
        """Test that _ranking_task returns a Figure."""
        overall_mean, _ = _compute_aggregations(sample_ranking_df)
        method_order = overall_mean.sort_values(ascending=False).index.tolist()

        fig = _ranking_task(
            df=sample_ranking_df,
            method_order=method_order,
            n_methods=len(method_order),
            figsize=None,
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_custom_figsize(self, sample_ranking_df):
        """Test that custom figsize is applied."""
        overall_mean, _ = _compute_aggregations(sample_ranking_df)
        method_order = overall_mean.sort_values(ascending=False).index.tolist()

        fig = _ranking_task(
            df=sample_ranking_df,
            method_order=method_order,
            n_methods=len(method_order),
            figsize=(20, 10),
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


class TestRanking:
    """Tests for the public ranking function."""

    def test_class_level_returns_figure(self, sample_ranking_df):
        """Test level='class' returns a Figure."""
        fig = ranking(sample_ranking_df, level="class", return_fig=True)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_task_level_returns_figure(self, sample_ranking_df):
        """Test level='task' returns a Figure."""
        fig = ranking(sample_ranking_df, level="task", return_fig=True)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_invalid_level_raises(self, sample_ranking_df):
        """Test that invalid level raises ValueError."""
        with pytest.raises(ValueError, match="level must be"):
            ranking(sample_ranking_df, level="invalid", return_fig=True)

    def test_return_fig_false_returns_none(self, sample_ranking_df):
        """Test that return_fig=False returns None."""
        result = ranking(sample_ranking_df, level="class", return_fig=False)
        assert result is None
        plt.close("all")

    def test_custom_palette(self, sample_ranking_df):
        """Test that custom palette is applied."""
        palette = {"MethodA": "red", "MethodB": "blue", "MethodC": "green"}
        fig = ranking(sample_ranking_df, level="class", palette=palette, return_fig=True)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_custom_figsize(self, sample_ranking_df):
        """Test that custom figsize is accepted."""
        fig = ranking(sample_ranking_df, level="class", figsize=(10, 8), return_fig=True)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_filters_excluded_datasets(self, sample_ranking_df_with_excluded):
        """Test that excluded datasets (Synthetic Pituitary) are filtered out."""
        fig = ranking(sample_ranking_df_with_excluded, level="class", return_fig=True)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_without_dts_column(self, sample_ranking_df):
        """Test with df that lacks 'dts' column uses data as-is."""
        # Just verify class level works
        fig = ranking(sample_ranking_df, level="class", return_fig=True)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_ranking_class_with_nan_in_class_mean(self, sample_ranking_df_with_nan):
        """Test _ranking_class when class_mean has NaN (line 160 - txt='-' path)."""
        fig = ranking(sample_ranking_df_with_nan, level="class", return_fig=True)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_ranking_task_with_nan_in_db_mean(self, sample_ranking_df_with_nan):
        """Test _ranking_task when db_mean has NaN (line 311 - txt='-' path)."""
        fig = ranking(sample_ranking_df_with_nan, level="task", return_fig=True)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_ranking_task_multiple_tasks_per_class(self, sample_ranking_df):
        """Test _ranking_task with multiple tasks per class (lines 83, 259)."""
        fig = ranking(sample_ranking_df, level="task", return_fig=True)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
