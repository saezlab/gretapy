"""Tests for pygreta.pl module (plotting functions)."""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

# Use non-interactive backend for testing
matplotlib.use("Agg")

from pygreta.pl._heatmap import _make_sim_mat, heatmap

# ============================================================================
# Fixtures specific to pl module tests
# ============================================================================


@pytest.fixture
def ocoeff_df():
    """DataFrame mimicking output from tl.ocoeff with overlap coefficients."""
    return pd.DataFrame(
        {
            "grn_a": ["GRN1", "GRN1", "GRN2"],
            "grn_b": ["GRN2", "GRN3", "GRN3"],
            "source": [0.8, 0.6, 0.7],
            "cre": [0.5, 0.4, 0.6],
            "target": [0.75, 0.55, 0.65],
            "edge": [0.4, 0.3, 0.5],
        }
    )


@pytest.fixture
def ocoeff_df_large():
    """Larger DataFrame with more GRNs for testing ordering."""
    return pd.DataFrame(
        {
            "grn_a": ["A", "A", "A", "B", "B", "C"],
            "grn_b": ["B", "C", "D", "C", "D", "D"],
            "source": [0.9, 0.8, 0.7, 0.6, 0.5, 0.4],
            "cre": [0.85, 0.75, 0.65, 0.55, 0.45, 0.35],
            "target": [0.8, 0.7, 0.6, 0.5, 0.4, 0.3],
            "edge": [0.7, 0.6, 0.5, 0.4, 0.3, 0.2],
        }
    )


# ============================================================================
# _make_sim_mat helper function tests
# ============================================================================


class TestMakeSimMat:
    """Tests for _make_sim_mat helper function."""

    def test_returns_dataframe(self, ocoeff_df):
        """Test that function returns a DataFrame."""
        result = _make_sim_mat(ocoeff_df, col="edge")

        assert isinstance(result, pd.DataFrame)

    def test_square_matrix(self, ocoeff_df):
        """Test that result is a square matrix."""
        result = _make_sim_mat(ocoeff_df, col="edge")

        assert result.shape[0] == result.shape[1]

    def test_correct_size(self, ocoeff_df):
        """Test that matrix size matches number of unique GRNs."""
        result = _make_sim_mat(ocoeff_df, col="edge")
        unique_grns = set(ocoeff_df["grn_a"]) | set(ocoeff_df["grn_b"])

        assert result.shape[0] == len(unique_grns)

    def test_diagonal_is_one(self, ocoeff_df):
        """Test that diagonal values are 1.0 (self-similarity)."""
        result = _make_sim_mat(ocoeff_df, col="edge")

        for idx in result.index:
            assert result.loc[idx, idx] == 1.0

    def test_symmetric_matrix(self, ocoeff_df):
        """Test that matrix is symmetric."""
        result = _make_sim_mat(ocoeff_df, col="edge")

        np.testing.assert_array_almost_equal(result.values, result.values.T)

    def test_correct_values(self, ocoeff_df):
        """Test that matrix contains correct values from input."""
        result = _make_sim_mat(ocoeff_df, col="edge")

        # Check specific values from ocoeff_df
        assert result.loc["GRN1", "GRN2"] == 0.4
        assert result.loc["GRN2", "GRN1"] == 0.4  # Symmetric
        assert result.loc["GRN1", "GRN3"] == 0.3
        assert result.loc["GRN2", "GRN3"] == 0.5

    def test_different_columns(self, ocoeff_df):
        """Test that different columns produce different matrices."""
        result_edge = _make_sim_mat(ocoeff_df, col="edge")
        result_source = _make_sim_mat(ocoeff_df, col="source")

        # Values should differ (except diagonal)
        assert result_edge.loc["GRN1", "GRN2"] != result_source.loc["GRN1", "GRN2"]

    def test_single_pair(self):
        """Test with single GRN pair."""
        df = pd.DataFrame(
            {
                "grn_a": ["A"],
                "grn_b": ["B"],
                "edge": [0.5],
            }
        )

        result = _make_sim_mat(df, col="edge")

        assert result.shape == (2, 2)
        assert result.loc["A", "B"] == 0.5
        assert result.loc["B", "A"] == 0.5
        assert result.loc["A", "A"] == 1.0
        assert result.loc["B", "B"] == 1.0


# ============================================================================
# heatmap function tests
# ============================================================================


class TestHeatmap:
    """Tests for heatmap function."""

    def test_invalid_level_raises(self, ocoeff_df):
        """Test that invalid level raises ValueError."""
        with pytest.raises(ValueError, match='level must be "source", "cre", "target", or "edge"'):
            heatmap(ocoeff_df, level="invalid")

    def test_valid_level_source(self, ocoeff_df):
        """Test that level='source' works."""
        result = heatmap(ocoeff_df, level="source", return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_valid_level_cre(self, ocoeff_df):
        """Test that level='cre' works."""
        result = heatmap(ocoeff_df, level="cre", return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_valid_level_target(self, ocoeff_df):
        """Test that level='target' works."""
        result = heatmap(ocoeff_df, level="target", return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_valid_level_edge(self, ocoeff_df):
        """Test that level='edge' works."""
        result = heatmap(ocoeff_df, level="edge", return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_return_fig_true(self, ocoeff_df):
        """Test that return_fig=True returns a Figure."""
        result = heatmap(ocoeff_df, return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_return_fig_false(self, ocoeff_df):
        """Test that return_fig=False returns None."""
        result = heatmap(ocoeff_df, return_fig=False)

        assert result is None
        plt.close("all")

    def test_custom_order(self, ocoeff_df_large):
        """Test that custom order is respected."""
        custom_order = ["D", "C", "B", "A"]
        result = heatmap(ocoeff_df_large, order=custom_order, return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_custom_title(self, ocoeff_df):
        """Test that custom title is accepted."""
        result = heatmap(ocoeff_df, title="Custom Title", return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_custom_cmap(self, ocoeff_df):
        """Test that custom colormap is accepted."""
        result = heatmap(ocoeff_df, cmap="viridis", return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_custom_vmin_vmax(self, ocoeff_df):
        """Test that custom vmin/vmax are accepted."""
        result = heatmap(ocoeff_df, vmin=0.2, vmax=0.8, return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_custom_dimensions(self, ocoeff_df):
        """Test that custom width/height are accepted."""
        result = heatmap(ocoeff_df, width=4, height=3, return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_default_title_from_level(self, ocoeff_df):
        """Test that default title is capitalized level name."""
        # We can't easily check the title text, but we verify no error
        result = heatmap(ocoeff_df, level="source", title=None, return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_alphabetical_order_by_default(self, ocoeff_df_large):
        """Test that GRNs are ordered alphabetically by default."""
        result = heatmap(ocoeff_df_large, return_fig=True)

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_all_parameters_combined(self, ocoeff_df):
        """Test with all parameters specified."""
        result = heatmap(
            ocoeff_df,
            level="target",
            order=["GRN1", "GRN2", "GRN3"],
            title="Test Heatmap",
            cmap="Blues",
            vmin=0.1,
            vmax=0.9,
            width=3,
            height=3,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)
