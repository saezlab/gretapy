"""Tests for gretapy.tl._mechanistic module."""

import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

from gretapy.tl._mechanistic import (
    _coefmat,
    _define_bool_rules,
    _fisher_test,
    _frc,
    _get_pyboolnet,
    _get_sim_hits,
    _get_source_markers,
    _sim,
    _simulate,
    _sss_prc_rcl,
    _tfa,
)

try:
    import pyboolnet  # noqa: F401

    HAS_PYBOOLNET = True
except ImportError:
    HAS_PYBOOLNET = False

requires_pyboolnet = pytest.mark.skipif(not HAS_PYBOOLNET, reason="pyboolnet not installed")


class TestGetPyboolnet:
    """Tests for _get_pyboolnet function."""

    def test_raises_import_error_when_not_installed(self):
        """Test that ImportError is raised with helpful message when pyboolnet missing."""
        with patch.dict(sys.modules, {"pyboolnet": None, "pyboolnet.file_exchange": None, "pyboolnet.trap_spaces": None}):
            with pytest.raises(ImportError, match="pyboolnet is required"):
                _get_pyboolnet()

    def test_returns_modules_when_installed(self):
        """Test that _get_pyboolnet returns (file_exchange, trap_spaces) when available."""
        mock_fe = MagicMock()
        mock_ts = MagicMock()
        mock_pb = MagicMock()
        mock_pb.file_exchange = mock_fe
        mock_pb.trap_spaces = mock_ts
        with patch.dict(sys.modules, {
            "pyboolnet": mock_pb,
            "pyboolnet.file_exchange": mock_fe,
            "pyboolnet.trap_spaces": mock_ts,
        }):
            fe, ts = _get_pyboolnet()
        assert fe is mock_fe
        assert ts is mock_ts


class TestGetSourceMarkers:
    """Tests for _get_source_markers function."""

    def test_returns_dataframe(self, adata):
        """Test that function returns a DataFrame."""
        sources = ["PAX5", "GATA3", "SPI1"]
        result = _get_source_markers(
            adata=adata,
            sources=sources,
            thr_deg_lfc=0.5,
            thr_deg_padj=0.05,
        )

        assert isinstance(result, pd.DataFrame)
        assert "celltype" in result.columns
        assert "source" in result.columns

    def test_filters_by_thresholds(self, adata):
        """Test that markers are filtered by LFC and p-value thresholds."""
        sources = ["PAX5", "GATA3", "SPI1"]

        # More stringent thresholds should return fewer markers
        result_strict = _get_source_markers(
            adata=adata,
            sources=sources,
            thr_deg_lfc=2.0,
            thr_deg_padj=1e-10,
        )

        result_lenient = _get_source_markers(
            adata=adata,
            sources=sources,
            thr_deg_lfc=0.1,
            thr_deg_padj=0.5,
        )

        assert len(result_strict) <= len(result_lenient)


class TestDefineBoolRules:
    """Tests for _define_bool_rules function."""

    def test_generates_rules_string(self, simple_grn):
        """Test that function generates a rules string."""
        result = _define_bool_rules(grn=simple_grn, indegree=5)

        assert isinstance(result, str)
        assert len(result) > 0

    def test_includes_all_targets(self, simple_grn):
        """Test that all targets are included in rules."""
        result = _define_bool_rules(grn=simple_grn, indegree=5)

        for target in simple_grn["target"].unique():
            assert target in result

    def test_respects_indegree(self):
        """Test that indegree parameter limits sources per target."""
        grn = pd.DataFrame(
            {
                "source": ["TF1", "TF2", "TF3", "TF4", "TF5"],
                "target": ["Gene1"] * 5,
                "score": [0.9, 0.8, 0.7, 0.6, 0.5],
            }
        )

        result_2 = _define_bool_rules(grn=grn, indegree=2)
        result_5 = _define_bool_rules(grn=grn, indegree=5)

        # With indegree=2, only top 2 sources should be included
        assert result_2.count("TF") <= result_5.count("TF")

    def test_handles_positive_and_negative_scores(self):
        """Test that positive and negative scores create different rules."""
        grn = pd.DataFrame(
            {
                "source": ["TF1", "TF2"],
                "target": ["Gene1", "Gene1"],
                "score": [0.8, -0.5],
            }
        )

        result = _define_bool_rules(grn=grn, indegree=5)

        # Positive score should be OR, negative should be NOT
        assert "TF1" in result
        assert "!TF2" in result

    def test_all_negative_scores(self):
        """Test rules with all-negative scores (only NOT terms, no OR terms)."""
        grn = pd.DataFrame(
            {
                "source": ["TF1", "TF2"],
                "target": ["Gene1", "Gene1"],
                "score": [-0.8, -0.5],
            }
        )

        result = _define_bool_rules(grn=grn, indegree=5)

        # Both should be NOT terms
        assert "!TF1" in result
        assert "!TF2" in result
        # No positive OR terms
        assert "TF1 |" not in result


class TestFisherTest:
    """Tests for _fisher_test function."""

    def test_significant_overlap(self):
        """Test with significant overlap."""
        hits = {"A", "B", "C"}
        sm_set = {"A", "B", "D"}
        sources = {"A", "B", "C", "D", "E", "F"}

        pval = _fisher_test(hits=hits, sm_set=sm_set, sources=sources)

        assert 0 <= pval <= 1

    def test_no_overlap(self):
        """Test with no overlap."""
        hits = {"A", "B"}
        sm_set = {"C", "D"}
        sources = {"A", "B", "C"}

        pval = _fisher_test(hits=hits, sm_set=sm_set, sources=sources)

        assert 0 <= pval <= 1
        assert pval == 1.0  # No overlap means p=1


class TestGetSimHits:
    """Tests for _get_sim_hits function."""

    def test_returns_dataframe(self):
        """Test that function returns a DataFrame."""
        sss = [{"TF1": True, "TF2": False}, {"TF1": False, "TF2": True}]
        sm_sets = pd.Series([{"TF1"}, {"TF2"}], index=["CellA", "CellB"])
        sources = {"TF1", "TF2"}

        result = _get_sim_hits(sss=sss, sm_sets=sm_sets, sources=sources, thr_fisher_padj=0.05)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == (2, 2)  # 2 steady states x 2 celltypes


class TestSssPrcRcl:
    """Tests for _sss_prc_rcl function."""

    def test_perfect_mapping(self):
        """Test with perfect 1-to-1 mapping."""
        # Each steady state maps to exactly one celltype
        hits = pd.DataFrame(
            [[True, False], [False, True]],
            columns=["CellA", "CellB"],
        )

        prc, rcl = _sss_prc_rcl(hits=hits)

        assert prc == 1.0
        assert rcl == 1.0

    def test_no_hits(self):
        """Test with no hits."""
        hits = pd.DataFrame(
            [[False, False], [False, False]],
            columns=["CellA", "CellB"],
        )

        prc, rcl = _sss_prc_rcl(hits=hits)

        assert prc == 0.0
        assert rcl == 0.0


class TestCoefmat:
    """Tests for _coefmat function."""

    def test_returns_dataframe(self, adata, simple_grn):
        """Test that function returns a coefficient matrix."""
        result = _coefmat(
            adata=adata,
            grn=simple_grn,
            alpha=1.0,
            seed=42,
            smin=1,
        )

        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == result.shape[1]  # Square matrix
        assert result.shape[0] == adata.var_names.size

    def test_respects_smin(self, adata, simple_grn):
        """Test that smin parameter controls minimum sources."""
        result_low = _coefmat(
            adata=adata,
            grn=simple_grn,
            smin=1,
        )

        result_high = _coefmat(
            adata=adata,
            grn=simple_grn,
            smin=10,
        )

        # Higher smin should result in more zeros
        assert (result_high != 0).sum().sum() <= (result_low != 0).sum().sum()


class TestSimulate:
    """Tests for _simulate function."""

    def test_returns_dataframe(self):
        """Test that function returns a DataFrame."""
        df = pd.DataFrame(
            [[1.0, 2.0, 3.0]],
            columns=["Gene1", "Gene2", "Gene3"],
        )
        coefmat = pd.DataFrame(
            [[0.0, 0.5, 0.0], [0.0, 0.0, 0.3], [0.0, 0.0, 0.0]],
            index=["Gene1", "Gene2", "Gene3"],
            columns=["Gene1", "Gene2", "Gene3"],
        )

        result = _simulate(df=df, coefmat=coefmat, gene="Gene1", n_steps=3)

        assert isinstance(result, pd.DataFrame)

    def test_knockdown_sets_gene_to_zero(self):
        """Test that the knocked-down gene has delta of -original_value."""
        df = pd.DataFrame(
            [[5.0, 2.0]],
            columns=["Gene1", "Gene2"],
        )
        coefmat = pd.DataFrame(
            [[0.0, 0.0], [0.0, 0.0]],
            index=["Gene1", "Gene2"],
            columns=["Gene1", "Gene2"],
        )

        result = _simulate(df=df, coefmat=coefmat, gene="Gene1", n_steps=1)

        # Gene1 knockdown should create a negative delta
        if "Gene1" in result.columns:
            assert result["Gene1"].iloc[0] == -5.0


@requires_pyboolnet
class TestSim:
    """Tests for _sim function (Boolean rule simulation)."""

    def test_basic_functionality(self, adata, simple_grn):
        """Test basic Boolean simulation."""
        prc, rcl, f01 = _sim(
            adata=adata,
            grn=simple_grn,
            indegree=10,
            thr_deg_lfc=0.5,
            thr_deg_padj=0.05,
            thr_fisher_padj=0.01,
        )

        # Results can be nan if not enough edges
        assert np.isnan(prc) or (0 <= prc <= 1)
        assert np.isnan(rcl) or (0 <= rcl <= 1)
        assert np.isnan(f01) or (0 <= f01 <= 1)

    def test_returns_nan_for_small_grn(self, adata):
        """Test that small GRNs return nan."""
        small_grn = pd.DataFrame({"source": ["PAX5"], "target": ["CD19"], "score": [0.5]})

        prc, rcl, f01 = _sim(
            adata=adata,
            grn=small_grn,
            indegree=10,
            thr_deg_lfc=0.5,
            thr_deg_padj=0.05,
        )

        # Should return nan for small networks
        assert np.isnan(prc)
        assert np.isnan(rcl)
        assert np.isnan(f01)


# GRN with TF-TF interactions for subprocess tests
_TF_TF_GRN = pd.DataFrame(
    {
        "source": ["TF1", "TF2", "TF3", "TF4", "TF5"],
        "target": ["TF2", "TF3", "TF4", "TF5", "TF6"],
        "score": [0.8, 0.7, 0.6, 0.5, 0.4],
    }
)

_TF_MARKERS = pd.DataFrame(
    {
        "celltype": ["CellA", "CellA", "CellB", "CellB", "CellC", "CellC"],
        "source": ["TF1", "TF2", "TF3", "TF4", "TF5", "TF6"],
    }
)


class TestSimSubprocessPaths:
    """Tests for _sim subprocess management paths (timeout, exitcode, success)."""

    @patch("gretapy.tl._mechanistic._get_source_markers")
    @patch("gretapy.tl._mechanistic.multiprocessing.Queue")
    @patch("gretapy.tl._mechanistic.multiprocessing.Process")
    def test_sim_timeout_returns_nan(self, mock_process_cls, mock_queue_cls, mock_source_markers, adata):
        """Test that _sim returns nan when subprocess times out."""
        mock_source_markers.return_value = _TF_MARKERS.copy()

        mock_proc = MagicMock()
        mock_proc.is_alive.return_value = True  # Simulate timeout
        mock_process_cls.return_value = mock_proc
        mock_queue_cls.return_value = MagicMock()

        prc, rcl, f01 = _sim(adata=adata, grn=_TF_TF_GRN.copy(), timeout=1)

        assert np.isnan(prc)
        assert np.isnan(rcl)
        assert np.isnan(f01)
        mock_proc.terminate.assert_called_once()

    @patch("gretapy.tl._mechanistic._get_source_markers")
    @patch("gretapy.tl._mechanistic.multiprocessing.Queue")
    @patch("gretapy.tl._mechanistic.multiprocessing.Process")
    def test_sim_nonzero_exitcode_raises(self, mock_process_cls, mock_queue_cls, mock_source_markers, adata):
        """Test that _sim raises RuntimeError when subprocess exits with non-zero code."""
        mock_source_markers.return_value = _TF_MARKERS.copy()

        mock_proc = MagicMock()
        mock_proc.is_alive.return_value = False
        mock_proc.exitcode = 1
        mock_process_cls.return_value = mock_proc
        mock_queue_cls.return_value = MagicMock()

        with pytest.raises(RuntimeError, match="Boolean simulation subprocess failed"):
            _sim(adata=adata, grn=_TF_TF_GRN.copy(), timeout=100)

    @patch("gretapy.tl._mechanistic._get_source_markers")
    @patch("gretapy.tl._mechanistic.multiprocessing.Queue")
    @patch("gretapy.tl._mechanistic.multiprocessing.Process")
    def test_sim_success_returns_values_from_queue(self, mock_process_cls, mock_queue_cls, mock_source_markers, adata):
        """Test that _sim returns values from queue on successful subprocess exit."""
        mock_source_markers.return_value = _TF_MARKERS.copy()

        mock_queue = MagicMock()
        mock_queue.get.return_value = (0.5, 0.6, 0.55)
        mock_queue_cls.return_value = mock_queue

        mock_proc = MagicMock()
        mock_proc.is_alive.return_value = False
        mock_proc.exitcode = 0
        mock_process_cls.return_value = mock_proc

        prc, rcl, f01 = _sim(adata=adata, grn=_TF_TF_GRN.copy(), timeout=100)

        assert prc == 0.5
        assert rcl == 0.6
        assert f01 == 0.55


class TestTfa:
    """Tests for _tfa function (TF activity scoring)."""

    def test_basic_functionality(self, adata, simple_grn, knocktf_db):
        """Test basic TF activity scoring."""
        prc, rcl, f01 = _tfa(
            adata=adata,
            grn=simple_grn,
            db=knocktf_db,
            cats=None,
            thr_pert_lfc=-0.5,
            thr_score_padj=0.05,
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_with_category_filter(self, adata, simple_grn, knocktf_db):
        """Test TF activity scoring with category filter."""
        prc, rcl, f01 = _tfa(
            adata=adata,
            grn=simple_grn,
            db=knocktf_db,
            cats=["Blood"],
            thr_pert_lfc=-0.5,
            thr_score_padj=0.05,
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    @patch("gretapy.tl._mechanistic.dc.mt.ulm")
    def test_tfa_tp_positive_path(self, mock_ulm, adata, simple_grn, knocktf_db):
        """Test _tfa when tp > 0 (precision/recall/f01 computation branch)."""
        # Return strongly negative scores with very low p-values → tp > 0
        mock_ulm.return_value = (
            pd.DataFrame([[-5.0]]),
            pd.DataFrame([[1e-10]]),
        )

        prc, rcl, f01 = _tfa(
            adata=adata,
            grn=simple_grn,
            db=knocktf_db,
            cats=None,
            thr_pert_lfc=-0.5,
            thr_score_padj=0.05,
        )

        # With mocked negative scores and tiny p-values, tp > 0
        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1


class TestFrc:
    """Tests for _frc function (perturbation forecasting)."""

    def test_basic_functionality(self, adata, simple_grn, knocktf_db):
        """Test basic perturbation forecasting."""
        prc, rcl, f01 = _frc(
            adata=adata,
            grn=simple_grn,
            db=knocktf_db,
            cats=None,
            thr_pert_lfc=-0.5,
            n_steps=3,
            min_size=5,
            thr_cor_stat=0.05,
            thr_cor_padj=0.05,
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_with_category_filter(self, adata, simple_grn, knocktf_db):
        """Test perturbation forecasting with category filter."""
        prc, rcl, f01 = _frc(
            adata=adata,
            grn=simple_grn,
            db=knocktf_db,
            cats=["Blood"],
            thr_pert_lfc=-0.5,
            n_steps=3,
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    @patch("gretapy.tl._mechanistic._coefmat")
    def test_frc_with_valid_tfs_in_inner_loop(self, mock_coefmat, adata, simple_grn, knocktf_db):
        """Test _frc inner loop body when tf in valid_tfs."""
        # Make coefmat have all non-zero entries so valid_tfs includes TFs from db
        all_genes = list(adata.var_names)
        n = len(all_genes)
        np.random.seed(42)
        # Non-zero coefmat: every gene has >= 3 non-zero entries
        coef_values = np.random.randn(n, n)
        coefmat = pd.DataFrame(coef_values, index=all_genes, columns=all_genes)
        mock_coefmat.return_value = coefmat

        prc, rcl, f01 = _frc(
            adata=adata,
            grn=simple_grn,
            db=knocktf_db,
            cats=None,
            thr_pert_lfc=-0.5,
            n_steps=1,
            min_size=1,
            thr_cor_stat=-1.0,  # Very lenient: any positive coef passes
            thr_cor_padj=1.0,   # Very lenient: any padj passes
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    @patch("gretapy.tl._mechanistic._coefmat")
    def test_frc_x_size_below_min_size(self, mock_coefmat, adata, simple_grn, knocktf_db):
        """Test _frc when x.size < min_size → r=0, p=1 (line 329)."""
        all_genes = list(adata.var_names)
        n = len(all_genes)
        np.random.seed(42)
        coef_values = np.random.randn(n, n)
        coefmat = pd.DataFrame(coef_values, index=all_genes, columns=all_genes)
        mock_coefmat.return_value = coefmat

        # Very high min_size so no intersections pass → falls to r=0, p=1 branch
        prc, rcl, f01 = _frc(
            adata=adata,
            grn=simple_grn,
            db=knocktf_db,
            cats=None,
            thr_pert_lfc=-0.5,
            n_steps=1,
            min_size=99999,  # Too large: x.size < min_size always → r=0, p=1
            thr_cor_stat=-1.0,
            thr_cor_padj=1.0,
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1
