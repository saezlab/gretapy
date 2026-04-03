"""Tests for gretapy.pl module (plotting functions)."""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

# Use non-interactive backend for testing
matplotlib.use("Agg")

import anndata as ad
import mudata as mu
import pyranges as pr

from gretapy.pl._heatmap import _make_sim_mat, heatmap
from gretapy.pl._links import (
    _get_gannot_data,
    _get_tss_window,
    _norm_score,
    links,
)

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

    def test_save_writes_file(self, ocoeff_df, tmp_path):
        """Test that save parameter writes a file."""
        out = tmp_path / "heatmap.png"
        result = heatmap(ocoeff_df, save=str(out))

        assert result is None
        assert out.exists()
        plt.close("all")


# ============================================================================
# Fixtures for links tests
# ============================================================================


@pytest.fixture
def gannot_pyranges():
    """Gene annotations as PyRanges for links tests."""
    return pr.PyRanges(
        pd.DataFrame(
            {
                "Chromosome": ["chr1", "chr1", "chr1"],
                "Start": [1000000, 1200000, 1400000],
                "End": [1050000, 1250000, 1450000],
                "Name": ["GENE_A", "GENE_B", "GENE_C"],
                "Score": [0, 0, 0],
                "Strand": ["+", "-", "+"],
            }
        )
    )


@pytest.fixture
def rna_pseudobulk():
    """RNA pseudobulk AnnData with cell types as obs and genes as var."""
    X = np.array([[1.0, 0.5, 0.2], [0.3, 1.2, 0.8], [0.1, 0.4, 1.5]])
    adata = ad.AnnData(X=X.astype(np.float32))
    adata.obs_names = ["CellType_A", "CellType_B", "CellType_C"]
    adata.var_names = ["GENE_A", "GENE_B", "GENE_C"]
    return adata


@pytest.fixture
def atac_pseudobulk():
    """ATAC pseudobulk AnnData with cell types as obs and CREs as var."""
    # CREs in chr1 region that overlap with gene annotations
    cre_names = [
        "chr1-950000-960000",
        "chr1-1020000-1030000",
        "chr1-1100000-1110000",
        "chr1-1180000-1190000",
        "chr1-1350000-1360000",
    ]
    X = np.array(
        [
            [0.8, 0.3, 0.5, 0.2, 0.7],
            [0.2, 0.9, 0.4, 0.6, 0.3],
            [0.4, 0.1, 0.8, 0.9, 0.5],
        ]
    )
    adata = ad.AnnData(X=X.astype(np.float32))
    adata.obs_names = ["CellType_A", "CellType_B", "CellType_C"]
    adata.var_names = cre_names
    return adata


@pytest.fixture
def grn_single():
    """Single GRN DataFrame for links tests."""
    return pd.DataFrame(
        {
            "source": ["TF1", "TF1", "TF2", "TF2", "TF1"],
            "target": ["GENE_A", "GENE_A", "GENE_A", "GENE_A", "GENE_B"],
            "cre": [
                "chr1-950000-960000",
                "chr1-1020000-1030000",
                "chr1-1100000-1110000",
                "chr1-950000-960000",
                "chr1-1180000-1190000",
            ],
            "score": [0.9, 0.7, 0.5, 0.8, 0.6],
        }
    )


@pytest.fixture
def grn_dict(grn_single):
    """Dictionary of GRNs for links tests."""
    grn2 = grn_single.copy()
    grn2["score"] = grn2["score"] * 0.8
    return {"GRN_Method1": grn_single, "GRN_Method2": grn2}


@pytest.fixture
def mdata_pseudobulk(rna_pseudobulk, atac_pseudobulk):
    """MuData object with rna and atac modalities for links tests."""
    return mu.MuData({"rna": rna_pseudobulk, "atac": atac_pseudobulk})


# ============================================================================
# _norm_score helper function tests
# ============================================================================


class TestNormScore:
    """Tests for _norm_score helper function."""

    def test_normalizes_to_0_1(self):
        """Test that values are normalized to 0-1 range."""
        x = np.array([10.0, 20.0, 30.0, 40.0])
        result = _norm_score(x)

        assert result.min() == 0.0
        assert result.max() == 1.0

    def test_constant_array_unchanged(self):
        """Test that constant arrays are returned unchanged."""
        x = np.array([5.0, 5.0, 5.0])
        result = _norm_score(x)

        np.testing.assert_array_equal(result, x)

    def test_negative_values(self):
        """Test normalization with negative values."""
        x = np.array([-10.0, 0.0, 10.0])
        result = _norm_score(x)

        assert result.min() == 0.0
        assert result.max() == 1.0
        assert result[1] == 0.5


# ============================================================================
# _get_tss_window helper function tests
# ============================================================================


class TestGetTssWindow:
    """Tests for _get_tss_window helper function."""

    def test_plus_strand_tss(self, gannot_pyranges):
        """Test TSS window for plus strand gene."""
        result = _get_tss_window(gannot_pyranges, "GENE_A", w_size=100000)

        # GENE_A is + strand, TSS should be at Start (1000000)
        assert result.df["Start"].values[0] == 1000000 - 100000
        assert result.df["End"].values[0] == 1000000 + 100000

    def test_minus_strand_tss(self, gannot_pyranges):
        """Test TSS window for minus strand gene."""
        result = _get_tss_window(gannot_pyranges, "GENE_B", w_size=100000)

        # GENE_B is - strand, TSS should be at End (1250000)
        assert result.df["Start"].values[0] == 1250000 - 100000
        assert result.df["End"].values[0] == 1250000 + 100000

    def test_window_size(self, gannot_pyranges):
        """Test that window size is correctly applied."""
        result = _get_tss_window(gannot_pyranges, "GENE_A", w_size=50000)
        window_size = result.df["End"].values[0] - result.df["Start"].values[0]

        assert window_size == 100000  # 2 * w_size


# ============================================================================
# _get_gannot_data helper function tests
# ============================================================================


class TestGetGannotData:
    """Tests for _get_gannot_data helper function."""

    def test_returns_correct_tuple_length(self, gannot_pyranges, atac_pseudobulk):
        """Test that function returns tuple with 7 elements."""
        result = _get_gannot_data(gannot_pyranges, "GENE_A", 500000, atac_pseudobulk.var_names)

        assert len(result) == 7

    def test_extracts_chromosome(self, gannot_pyranges, atac_pseudobulk):
        """Test that chromosome is correctly extracted."""
        x_min, x_max, gs_gr, tss, chromosome, strand, cres_gr = _get_gannot_data(
            gannot_pyranges, "GENE_A", 500000, atac_pseudobulk.var_names
        )

        assert chromosome == "chr1"

    def test_extracts_strand(self, gannot_pyranges, atac_pseudobulk):
        """Test that strand is correctly extracted."""
        _, _, _, _, _, strand, _ = _get_gannot_data(gannot_pyranges, "GENE_A", 500000, atac_pseudobulk.var_names)

        assert strand == "+"

    def test_extracts_tss_plus_strand(self, gannot_pyranges, atac_pseudobulk):
        """Test TSS extraction for plus strand."""
        _, _, _, tss, _, strand, _ = _get_gannot_data(gannot_pyranges, "GENE_A", 500000, atac_pseudobulk.var_names)

        assert strand == "+"
        assert tss == 1000000

    def test_extracts_tss_minus_strand(self, gannot_pyranges, atac_pseudobulk):
        """Test TSS extraction for minus strand."""
        _, _, _, tss, _, strand, _ = _get_gannot_data(gannot_pyranges, "GENE_B", 500000, atac_pseudobulk.var_names)

        assert strand == "-"
        assert tss == 1250000

    def test_cres_are_pyranges(self, gannot_pyranges, atac_pseudobulk):
        """Test that CREs are returned as PyRanges."""
        _, _, _, _, _, _, cres_gr = _get_gannot_data(gannot_pyranges, "GENE_A", 500000, atac_pseudobulk.var_names)

        assert isinstance(cres_gr, pr.PyRanges)


# ============================================================================
# links function tests
# ============================================================================


class TestLinks:
    """Tests for links function."""

    def test_returns_figure_when_return_fig_true(self, mdata_pseudobulk, grn_single, gannot_pyranges):
        """Test that return_fig=True returns a Figure."""
        result = links(
            mdata_pseudobulk,
            grn_single,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_returns_none_when_return_fig_false(self, mdata_pseudobulk, grn_single, gannot_pyranges):
        """Test that return_fig=False returns None."""
        result = links(
            mdata_pseudobulk,
            grn_single,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            return_fig=False,
        )

        assert result is None
        plt.close("all")

    def test_single_grn_dataframe(self, mdata_pseudobulk, grn_single, gannot_pyranges):
        """Test with single GRN as DataFrame."""
        result = links(
            mdata_pseudobulk,
            grn_single,
            target="GENE_A",
            tfs=["TF1", "TF2"],
            gannot=gannot_pyranges,
            w_size=500000,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_multiple_grns_dict(self, mdata_pseudobulk, grn_dict, gannot_pyranges):
        """Test with multiple GRNs as dict."""
        result = links(
            mdata_pseudobulk,
            grn_dict,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_multiple_tfs(self, mdata_pseudobulk, grn_single, gannot_pyranges):
        """Test plotting multiple TFs."""
        result = links(
            mdata_pseudobulk,
            grn_single,
            target="GENE_A",
            tfs=["TF1", "TF2"],
            gannot=gannot_pyranges,
            w_size=500000,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        # Should have 4 subplots: 2 TFs + 1 omics + 1 gannot
        assert len(result.axes) == 4
        plt.close(result)

    def test_custom_palette(self, mdata_pseudobulk, grn_dict, gannot_pyranges):
        """Test with custom color palette."""
        palette = {"GRN_Method1": "red", "GRN_Method2": "blue"}
        result = links(
            mdata_pseudobulk,
            grn_dict,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            palette=palette,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_custom_expr_cmap_list(self, mdata_pseudobulk, grn_single, gannot_pyranges):
        """Test with custom expression colormap as list."""
        result = links(
            mdata_pseudobulk,
            grn_single,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            expr_cmap=["white", "red"],
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_custom_expr_cmap_string(self, mdata_pseudobulk, grn_single, gannot_pyranges):
        """Test with custom expression colormap as string."""
        result = links(
            mdata_pseudobulk,
            grn_single,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            expr_cmap="viridis",
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_custom_figsize(self, mdata_pseudobulk, grn_single, gannot_pyranges):
        """Test with custom figure size."""
        result = links(
            mdata_pseudobulk,
            grn_single,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            figsize=(5, 6),
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_custom_dpi(self, mdata_pseudobulk, grn_single, gannot_pyranges):
        """Test with custom DPI."""
        result = links(
            mdata_pseudobulk,
            grn_single,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            dpi=100,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        assert result.dpi == 100
        plt.close(result)

    def test_all_parameters_combined(self, mdata_pseudobulk, grn_dict, gannot_pyranges):
        """Test with all parameters specified."""
        palette = {"GRN_Method1": "green", "GRN_Method2": "orange"}
        result = links(
            mdata_pseudobulk,
            grn_dict,
            target="GENE_A",
            tfs=["TF1", "TF2"],
            gannot=gannot_pyranges,
            w_size=300000,
            palette=palette,
            expr_cmap="Blues",
            figsize=(4, 5),
            dpi=120,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_sample_col_triggers_pseudobulk_aggregation(self, grn_single, gannot_pyranges):
        """Test that sample_col in mdata.obs triggers pseudobulk aggregation (lines 288-291)."""
        import anndata as ad
        import mudata as mu

        np.random.seed(42)
        n_cells = 12

        rna_X = np.abs(np.random.rand(n_cells, 3)).astype(np.float32)
        rna = ad.AnnData(X=rna_X)
        rna.var_names = ["GENE_A", "GENE_B", "GENE_C"]
        rna.obs_names = [f"Cell{i}" for i in range(n_cells)]

        atac_X = np.abs(np.random.rand(n_cells, 5)).astype(np.float32)
        atac = ad.AnnData(X=atac_X)
        atac.var_names = [
            "chr1-950000-960000",
            "chr1-1020000-1030000",
            "chr1-1100000-1110000",
            "chr1-1180000-1190000",
            "chr1-1350000-1360000",
        ]
        atac.obs_names = [f"Cell{i}" for i in range(n_cells)]

        mdata = mu.MuData({"rna": rna, "atac": atac})
        # Add celltype column to mdata.obs so sample_col is in mdata.obs.columns
        celltypes = ["CellType_A"] * 4 + ["CellType_B"] * 4 + ["CellType_C"] * 4
        mdata.obs["celltype"] = celltypes

        result = links(
            mdata,
            grn_single,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            sample_col="celltype",
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_gannot_as_string_loads_annotation(self, mdata_pseudobulk, grn_single):
        """Test that gannot as string triggers show_genome_annotation call (line 295)."""
        from unittest.mock import patch

        mock_gannot = pr.PyRanges(
            pd.DataFrame(
                {
                    "Chromosome": ["chr1", "chr1", "chr1"],
                    "Start": [1000000, 1200000, 1400000],
                    "End": [1050000, 1250000, 1450000],
                    "Name": ["GENE_A", "GENE_B", "GENE_C"],
                    "Score": [0, 0, 0],
                    "Strand": ["+", "-", "+"],
                }
            )
        )

        with patch("gretapy.pl._links.gt.show_genome_annotation", return_value=mock_gannot) as mock_show:
            result = links(
                mdata_pseudobulk,
                grn_single,
                target="GENE_A",
                tfs=["TF1"],
                gannot="hg38",  # string triggers show_genome_annotation
                w_size=500000,
                return_fig=True,
            )

        mock_show.assert_called_once_with(organism="hg38")
        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_plot_links_with_empty_grn_for_one_grn(self, mdata_pseudobulk, gannot_pyranges):
        """Test _plot_links when one GRN has no TF1 links (line 80 - is_empty append True).

        Both GRNs must have links to the target (so they appear in links_df),
        but one GRN has no TF1 links → empty when filtered by TF in _plot_links.
        """
        # GRN1: TF1 links to GENE_A (appears in links_df for target=GENE_A)
        grn_with_tf1 = pd.DataFrame(
            {
                "source": ["TF1"],
                "target": ["GENE_A"],
                "cre": ["chr1-1020000-1030000"],
                "score": [0.9],
            }
        )
        # GRN2: Only TF2 links to GENE_A (no TF1 links, but appears in links_df)
        grn_tf2_only = pd.DataFrame(
            {
                "source": ["TF2"],
                "target": ["GENE_A"],
                "cre": ["chr1-1100000-1110000"],
                "score": [0.5],
            }
        )
        # Plot with TF1 only: GRN2 (grn_tf2_only) has no TF1 links → is_empty.append(True)
        grn_dict_mixed = {"GRN_With_TF1": grn_with_tf1, "GRN_TF2_Only": grn_tf2_only}

        result = links(
            mdata_pseudobulk,
            grn_dict_mixed,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_plot_links_cre_overlapping_tss(self, mdata_pseudobulk, gannot_pyranges):
        """Test arc computation when CRE overlaps TSS (lines 95-97: cre_start < tss < cre_end)."""
        # GENE_A is + strand, TSS at 1000000
        # CRE that overlaps TSS: start < 1000000 < end
        grn_overlapping = pd.DataFrame(
            {
                "source": ["TF1"],
                "target": ["GENE_A"],
                "cre": ["chr1-990000-1010000"],  # overlaps TSS at 1000000
                "score": [0.8],
            }
        )

        result = links(
            mdata_pseudobulk,
            grn_overlapping,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)

    def test_save_writes_file(self, mdata_pseudobulk, grn_single, gannot_pyranges, tmp_path):
        """Test that save parameter writes a file."""
        out = tmp_path / "links.png"
        result = links(
            mdata_pseudobulk,
            grn_single,
            target="GENE_A",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            save=str(out),
        )

        assert result is None
        assert out.exists()
        plt.close("all")

    def test_plot_links_minus_strand_cre_positions(self, mdata_pseudobulk, gannot_pyranges):
        """Test arc computation for minus-strand gene (lines 99-107)."""
        # GENE_B is - strand, TSS at End=1250000
        # Test different CRE positions relative to TSS
        grn_minus = pd.DataFrame(
            {
                "source": ["TF1", "TF1", "TF1"],
                "target": ["GENE_B", "GENE_B", "GENE_B"],
                "cre": [
                    "chr1-1100000-1200000",   # cre_end < tss (1250000)
                    "chr1-1300000-1400000",   # cre_start > tss (1250000)
                    "chr1-1200000-1300000",   # cre_start < tss < cre_end
                ],
                "score": [0.9, 0.8, 0.7],
            }
        )

        result = links(
            mdata_pseudobulk,
            grn_minus,
            target="GENE_B",
            tfs=["TF1"],
            gannot=gannot_pyranges,
            w_size=500000,
            return_fig=True,
        )

        assert isinstance(result, plt.Figure)
        plt.close(result)
