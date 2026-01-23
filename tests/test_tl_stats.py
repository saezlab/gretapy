"""Tests for gretapy.tl._stats module."""

import numpy as np
import pandas as pd
import pytest

from gretapy.tl._stats import _get_grn_stats, _map_regions, _ocoeff, ocoeff, stats

# ============================================================================
# Fixtures specific to stats module tests
# ============================================================================


@pytest.fixture
def grn_a():
    """First GRN for overlap tests."""
    return pd.DataFrame(
        {
            "source": ["TF1", "TF1", "TF2", "TF2", "TF3"],
            "target": ["Gene1", "Gene2", "Gene3", "Gene4", "Gene5"],
            "cre": ["chr1-100-200", "chr1-300-400", "chr2-100-200", "chr2-300-400", "chr3-100-200"],
            "score": [0.8, 0.7, 0.9, 0.6, 0.5],
        }
    )


@pytest.fixture
def grn_b():
    """Second GRN for overlap tests (partial overlap with grn_a)."""
    return pd.DataFrame(
        {
            "source": ["TF1", "TF1", "TF2", "TF4", "TF4"],
            "target": ["Gene1", "Gene3", "Gene3", "Gene6", "Gene7"],
            "cre": ["chr1-100-200", "chr1-500-600", "chr2-100-200", "chr4-100-200", "chr4-300-400"],
            "score": [0.9, 0.6, 0.8, 0.7, 0.5],
        }
    )


@pytest.fixture
def grn_c():
    """Third GRN for overlap tests (no overlap with grn_a)."""
    return pd.DataFrame(
        {
            "source": ["TF5", "TF5", "TF6"],
            "target": ["Gene8", "Gene9", "Gene10"],
            "cre": ["chr5-100-200", "chr5-300-400", "chr6-100-200"],
            "score": [0.8, 0.7, 0.6],
        }
    )


@pytest.fixture
def grn_no_cre():
    """GRN without cre column."""
    return pd.DataFrame(
        {
            "source": ["TF1", "TF1", "TF2"],
            "target": ["Gene1", "Gene2", "Gene3"],
            "score": [0.8, 0.7, 0.9],
        }
    )


@pytest.fixture
def grn_with_duplicates():
    """GRN with duplicate edges."""
    return pd.DataFrame(
        {
            "source": ["TF1", "TF1", "TF1", "TF2"],
            "target": ["Gene1", "Gene1", "Gene2", "Gene3"],
            "cre": ["chr1-100-200", "chr1-100-200", "chr1-300-400", "chr2-100-200"],
            "score": [0.8, 0.8, 0.7, 0.9],
        }
    )


@pytest.fixture
def empty_grn():
    """Empty GRN DataFrame."""
    return pd.DataFrame(columns=["source", "target", "cre", "score"])


@pytest.fixture
def grn_overlapping_cres():
    """GRN with CREs that overlap with grn_a but are not identical."""
    return pd.DataFrame(
        {
            "source": ["TF1", "TF1", "TF2", "TF4"],
            "target": ["Gene1", "Gene2", "Gene3", "Gene6"],
            # chr1-150-250 overlaps with chr1-100-200
            # chr1-350-450 overlaps with chr1-300-400
            # chr2-150-250 overlaps with chr2-100-200
            # chr4-100-200 has no overlap with grn_a
            "cre": ["chr1-150-250", "chr1-350-450", "chr2-150-250", "chr4-100-200"],
            "score": [0.8, 0.7, 0.9, 0.6],
        }
    )


@pytest.fixture
def grn_non_overlapping_cres():
    """GRN with CREs that do not overlap with grn_a."""
    return pd.DataFrame(
        {
            "source": ["TF1", "TF2"],
            "target": ["Gene1", "Gene2"],
            # These CREs are on different chromosomes or non-overlapping positions
            "cre": ["chr10-100-200", "chr11-100-200"],
            "score": [0.8, 0.7],
        }
    )


# ============================================================================
# _map_regions helper function tests
# ============================================================================


class TestMapRegions:
    """Tests for _map_regions helper function."""

    def test_returns_dataframe(self):
        """Test that function returns a DataFrame."""
        regions_a = np.array(["chr1-100-200", "chr1-300-400"])
        regions_b = np.array(["chr1-150-250", "chr1-350-450"])

        result = _map_regions(regions_a, regions_b)

        assert isinstance(result, pd.DataFrame)

    def test_output_columns(self):
        """Test that output has correct columns."""
        regions_a = np.array(["chr1-100-200"])
        regions_b = np.array(["chr1-150-250"])

        result = _map_regions(regions_a, regions_b)

        assert set(result.columns) == {"region_a", "region_b"}

    def test_overlapping_regions(self):
        """Test that overlapping regions are correctly mapped."""
        regions_a = np.array(["chr1-100-200", "chr1-300-400"])
        regions_b = np.array(["chr1-150-250", "chr1-500-600"])

        result = _map_regions(regions_a, regions_b)

        # chr1-150-250 overlaps with chr1-100-200
        assert len(result) == 1
        assert result.iloc[0]["region_a"] == "chr1-100-200"
        assert result.iloc[0]["region_b"] == "chr1-150-250"

    def test_no_overlapping_regions(self):
        """Test with no overlapping regions."""
        regions_a = np.array(["chr1-100-200", "chr1-300-400"])
        regions_b = np.array(["chr2-100-200", "chr3-100-200"])

        result = _map_regions(regions_a, regions_b)

        assert result.empty

    def test_multiple_overlaps(self):
        """Test when one region overlaps with multiple regions."""
        regions_a = np.array(["chr1-100-300"])  # Large region
        regions_b = np.array(["chr1-50-150", "chr1-200-350"])  # Two overlapping regions

        result = _map_regions(regions_a, regions_b)

        # Both regions_b overlap with the single region_a
        assert len(result) == 2
        assert set(result["region_b"]) == {"chr1-50-150", "chr1-200-350"}
        assert all(result["region_a"] == "chr1-100-300")

    def test_identical_regions(self):
        """Test with identical regions."""
        regions_a = np.array(["chr1-100-200", "chr2-100-200"])
        regions_b = np.array(["chr1-100-200", "chr2-100-200"])

        result = _map_regions(regions_a, regions_b)

        assert len(result) == 2

    def test_empty_regions_a(self):
        """Test with empty regions_a."""
        regions_a = np.array([])
        regions_b = np.array(["chr1-100-200"])

        result = _map_regions(regions_a, regions_b)

        assert result.empty

    def test_empty_regions_b(self):
        """Test with empty regions_b."""
        regions_a = np.array(["chr1-100-200"])
        regions_b = np.array([])

        result = _map_regions(regions_a, regions_b)

        assert result.empty


# ============================================================================
# _ocoeff helper function tests
# ============================================================================


class TestOcoeffHelper:
    """Tests for _ocoeff helper function."""

    def test_returns_float(self, grn_a, grn_b):
        """Test that function returns a float."""
        result = _ocoeff(grn_a, grn_b, on=["source"])

        assert isinstance(result, float)

    def test_perfect_overlap(self, grn_a):
        """Test overlap coefficient with identical DataFrames."""
        result = _ocoeff(grn_a, grn_a.copy(), on=["source"])

        assert result == 1.0

    def test_no_overlap(self, grn_a, grn_c):
        """Test overlap coefficient with no overlap."""
        result = _ocoeff(grn_a, grn_c, on=["source"])

        assert result == 0.0

    def test_partial_overlap_source(self, grn_a, grn_b):
        """Test partial overlap on source column."""
        result = _ocoeff(grn_a, grn_b, on=["source"])

        # grn_a has TF1, TF2, TF3 (3 unique)
        # grn_b has TF1, TF2, TF4 (3 unique)
        # Intersection: TF1, TF2 (2)
        # ocoeff = 2 / min(3, 3) = 2/3
        assert result == pytest.approx(2 / 3)

    def test_partial_overlap_target(self, grn_a, grn_b):
        """Test partial overlap on target column."""
        result = _ocoeff(grn_a, grn_b, on=["target"])

        # grn_a has Gene1-5 (5 unique)
        # grn_b has Gene1, Gene3, Gene6, Gene7 (4 unique)
        # Intersection: Gene1, Gene3 (2)
        # ocoeff = 2 / min(5, 4) = 2/4 = 0.5
        assert result == pytest.approx(0.5)

    def test_partial_overlap_edge(self, grn_a, grn_b):
        """Test partial overlap on edge (source, target) columns."""
        result = _ocoeff(grn_a, grn_b, on=["source", "target"])

        # grn_a edges: (TF1,Gene1), (TF1,Gene2), (TF2,Gene3), (TF2,Gene4), (TF3,Gene5) - 5
        # grn_b edges: (TF1,Gene1), (TF1,Gene3), (TF2,Gene3), (TF4,Gene6), (TF4,Gene7) - 5
        # Intersection: (TF1,Gene1), (TF2,Gene3) - 2
        # ocoeff = 2 / min(5, 5) = 0.4
        assert result == pytest.approx(0.4)

    def test_empty_dataframe_returns_zero(self, grn_a, empty_grn):
        """Test that empty DataFrame returns 0.0."""
        result = _ocoeff(grn_a, empty_grn, on=["source"])

        assert result == 0.0

    def test_both_empty_returns_zero(self, empty_grn):
        """Test that two empty DataFrames return 0.0."""
        result = _ocoeff(empty_grn, empty_grn.copy(), on=["source"])

        assert result == 0.0

    def test_handles_duplicates(self, grn_with_duplicates, grn_a):
        """Test that duplicates are handled correctly (dropped)."""
        result = _ocoeff(grn_with_duplicates, grn_a, on=["source", "target"])

        # Should work without error and deduplicate
        assert 0.0 <= result <= 1.0

    def test_asymmetric_sizes(self):
        """Test overlap coefficient with asymmetric sizes."""
        df_a = pd.DataFrame({"source": ["TF1", "TF2"]})
        df_b = pd.DataFrame({"source": ["TF1", "TF2", "TF3", "TF4", "TF5"]})

        result = _ocoeff(df_a, df_b, on=["source"])

        # Intersection: TF1, TF2 (2)
        # ocoeff = 2 / min(2, 5) = 2/2 = 1.0
        assert result == 1.0

    def test_use_overlap_false_requires_exact_match(self, grn_a, grn_overlapping_cres):
        """Test that use_overlap=False requires exact CRE matches."""
        result = _ocoeff(grn_a, grn_overlapping_cres, on=["cre"], use_overlap=False)

        # No exact CRE matches between grn_a and grn_overlapping_cres
        assert result == 0.0

    def test_use_overlap_true_matches_overlapping_cres(self, grn_a, grn_overlapping_cres):
        """Test that use_overlap=True matches overlapping CREs."""
        result = _ocoeff(grn_a, grn_overlapping_cres, on=["cre"], use_overlap=True)

        # grn_a CREs: chr1-100-200, chr1-300-400, chr2-100-200, chr2-300-400, chr3-100-200 (5 unique)
        # grn_overlapping_cres: chr1-150-250, chr1-350-450, chr2-150-250, chr4-100-200 (4 unique)
        # Overlapping pairs:
        #   chr1-150-250 overlaps chr1-100-200
        #   chr1-350-450 overlaps chr1-300-400
        #   chr2-150-250 overlaps chr2-100-200
        # 3 CREs from smaller set (4) have overlaps
        # ocoeff = 3 / min(5, 4) = 3/4 = 0.75
        assert result == pytest.approx(0.75)

    def test_use_overlap_no_overlapping_cres(self, grn_a, grn_non_overlapping_cres):
        """Test use_overlap=True with non-overlapping CREs."""
        result = _ocoeff(grn_a, grn_non_overlapping_cres, on=["cre"], use_overlap=True)

        # No overlapping regions
        assert result == 0.0

    def test_use_overlap_identical_cres(self, grn_a):
        """Test use_overlap=True with identical CREs."""
        result = _ocoeff(grn_a, grn_a.copy(), on=["cre"], use_overlap=True)

        # Identical CREs should give perfect overlap
        assert result == 1.0

    def test_use_overlap_only_affects_cre_column(self, grn_a, grn_b):
        """Test that use_overlap only affects CRE comparison."""
        # For non-CRE columns, use_overlap should not change behavior
        result_with = _ocoeff(grn_a, grn_b, on=["source"], use_overlap=True)
        result_without = _ocoeff(grn_a, grn_b, on=["source"], use_overlap=False)

        assert result_with == result_without

    def test_use_overlap_with_edge_columns(self, grn_a, grn_b):
        """Test that use_overlap does not affect multi-column comparisons."""
        result_with = _ocoeff(grn_a, grn_b, on=["source", "target"], use_overlap=True)
        result_without = _ocoeff(grn_a, grn_b, on=["source", "target"], use_overlap=False)

        assert result_with == result_without

    def test_use_overlap_asymmetric_sizes(self):
        """Test use_overlap with asymmetric CRE set sizes."""
        df_a = pd.DataFrame({"cre": ["chr1-100-200", "chr1-300-400"]})  # 2 CREs
        df_b = pd.DataFrame(
            {"cre": ["chr1-150-250", "chr1-350-450", "chr2-100-200", "chr3-100-200", "chr4-100-200"]}
        )  # 5 CREs

        result = _ocoeff(df_a, df_b, on=["cre"], use_overlap=True)

        # df_a is smaller (2)
        # Both CREs in df_a overlap with CREs in df_b
        # ocoeff = 2 / min(2, 5) = 2/2 = 1.0
        assert result == 1.0

    def test_use_overlap_empty_dataframe(self, grn_a, empty_grn):
        """Test use_overlap with empty DataFrame."""
        result = _ocoeff(grn_a, empty_grn, on=["cre"], use_overlap=True)

        assert result == 0.0


# ============================================================================
# ocoeff function tests
# ============================================================================


class TestOcoeff:
    """Tests for ocoeff function."""

    def test_returns_dataframe(self, grn_a, grn_b):
        """Test that function returns a DataFrame."""
        result = ocoeff({"grn_a": grn_a, "grn_b": grn_b})

        assert isinstance(result, pd.DataFrame)

    def test_output_columns(self, grn_a, grn_b):
        """Test that output has correct columns."""
        result = ocoeff({"grn_a": grn_a, "grn_b": grn_b})

        expected_columns = ["grn_a", "grn_b", "source", "cre", "target", "edge"]
        assert list(result.columns) == expected_columns

    def test_single_pair(self, grn_a, grn_b):
        """Test with single pair of GRNs."""
        result = ocoeff({"A": grn_a, "B": grn_b})

        assert len(result) == 1
        assert result.iloc[0]["grn_a"] == "A"
        assert result.iloc[0]["grn_b"] == "B"

    def test_three_grns(self, grn_a, grn_b, grn_c):
        """Test with three GRNs (3 pairs)."""
        result = ocoeff({"A": grn_a, "B": grn_b, "C": grn_c})

        # 3 GRNs = 3 choose 2 = 3 pairs
        assert len(result) == 3

    def test_four_grns(self, grn_a, grn_b, grn_c):
        """Test with four GRNs (6 pairs)."""
        grn_d = grn_a.copy()
        result = ocoeff({"A": grn_a, "B": grn_b, "C": grn_c, "D": grn_d})

        # 4 GRNs = 4 choose 2 = 6 pairs
        assert len(result) == 6

    def test_less_than_two_grns_raises(self, grn_a):
        """Test that less than 2 GRNs raises ValueError."""
        with pytest.raises(ValueError, match="At least 2 GRNs are required"):
            ocoeff({"only_one": grn_a})

    def test_empty_dict_raises(self):
        """Test that empty dict raises ValueError."""
        with pytest.raises(ValueError, match="At least 2 GRNs are required"):
            ocoeff({})

    def test_values_in_range(self, grn_a, grn_b, grn_c):
        """Test that all coefficients are in [0, 1] range."""
        result = ocoeff({"A": grn_a, "B": grn_b, "C": grn_c})

        for col in ["source", "target", "edge"]:
            assert all(result[col] >= 0)
            assert all(result[col] <= 1)

    def test_identical_grns_have_perfect_overlap(self, grn_a):
        """Test that identical GRNs have overlap coefficient of 1.0."""
        result = ocoeff({"A": grn_a, "B": grn_a.copy()})

        assert result.iloc[0]["source"] == 1.0
        assert result.iloc[0]["cre"] == 1.0
        assert result.iloc[0]["target"] == 1.0
        assert result.iloc[0]["edge"] == 1.0

    def test_no_overlap_has_zero_coefficient(self, grn_a, grn_c):
        """Test that non-overlapping GRNs have coefficient of 0.0."""
        result = ocoeff({"A": grn_a, "C": grn_c})

        assert result.iloc[0]["source"] == 0.0
        assert result.iloc[0]["cre"] == 0.0
        assert result.iloc[0]["target"] == 0.0
        assert result.iloc[0]["edge"] == 0.0

    def test_missing_cre_column_returns_nan(self, grn_no_cre, grn_a):
        """Test that missing cre column returns NaN for cre coefficient."""
        result = ocoeff({"A": grn_no_cre, "B": grn_a})

        assert np.isnan(result.iloc[0]["cre"])

    def test_both_missing_cre_returns_nan(self, grn_no_cre):
        """Test that both missing cre columns returns NaN."""
        grn_no_cre_2 = grn_no_cre.copy()
        result = ocoeff({"A": grn_no_cre, "B": grn_no_cre_2})

        assert np.isnan(result.iloc[0]["cre"])

    def test_cre_overlap_uses_genomic_overlap(self, grn_a, grn_overlapping_cres):
        """Test that CRE comparison uses genomic overlap, not exact match."""
        result = ocoeff({"A": grn_a, "B": grn_overlapping_cres})

        # CREs are not identical but overlap genomically
        # Should have non-zero overlap coefficient
        assert result.iloc[0]["cre"] > 0.0
        assert result.iloc[0]["cre"] == pytest.approx(0.75)

    def test_cre_overlap_non_overlapping_regions(self, grn_a, grn_non_overlapping_cres):
        """Test CRE overlap with non-overlapping genomic regions."""
        result = ocoeff({"A": grn_a, "B": grn_non_overlapping_cres})

        # CREs do not overlap genomically
        assert result.iloc[0]["cre"] == 0.0

    def test_cre_overlap_vs_exact_match_different(self, grn_a, grn_overlapping_cres):
        """Test that genomic overlap gives different result than exact match would."""
        result = ocoeff({"A": grn_a, "B": grn_overlapping_cres})

        # With exact match, CRE overlap would be 0.0
        # With genomic overlap, it should be > 0.0
        assert result.iloc[0]["cre"] > 0.0

        # But source/target/edge still use exact match
        # so they may have different values
        assert result.iloc[0]["source"] >= 0.0
        assert result.iloc[0]["target"] >= 0.0
        assert result.iloc[0]["edge"] >= 0.0


# ============================================================================
# _get_grn_stats helper function tests
# ============================================================================


class TestGetGrnStats:
    """Tests for _get_grn_stats helper function."""

    def test_returns_tuple(self, grn_a):
        """Test that function returns a tuple."""
        result = _get_grn_stats(grn_a)

        assert isinstance(result, tuple)
        assert len(result) == 5

    def test_correct_n_sources(self, grn_a):
        """Test correct number of unique sources."""
        n_s, _, _, _, _ = _get_grn_stats(grn_a)

        # grn_a has TF1, TF2, TF3
        assert n_s == 3

    def test_correct_n_cres(self, grn_a):
        """Test correct number of unique CREs."""
        _, n_c, _, _, _ = _get_grn_stats(grn_a)

        # grn_a has 5 unique CREs
        assert n_c == 5

    def test_correct_n_targets(self, grn_a):
        """Test correct number of unique targets."""
        _, _, n_t, _, _ = _get_grn_stats(grn_a)

        # grn_a has Gene1-5
        assert n_t == 5

    def test_correct_n_edges(self, grn_a):
        """Test correct number of unique edges."""
        _, _, _, n_e, _ = _get_grn_stats(grn_a)

        # grn_a has 5 unique (source, target) pairs
        assert n_e == 5

    def test_correct_mean_regulon_size(self, grn_a):
        """Test correct mean regulon size."""
        _, _, _, _, n_r = _get_grn_stats(grn_a)

        # TF1: 2 targets, TF2: 2 targets, TF3: 1 target
        # Mean = (2 + 2 + 1) / 3 = 5/3
        assert n_r == pytest.approx(5 / 3)

    def test_no_cre_column_returns_zero(self, grn_no_cre):
        """Test that missing cre column returns 0 for n_cres."""
        _, n_c, _, _, _ = _get_grn_stats(grn_no_cre)

        assert n_c == 0

    def test_handles_duplicates(self, grn_with_duplicates):
        """Test that duplicates are handled correctly."""
        n_s, _, _, n_e, _ = _get_grn_stats(grn_with_duplicates)

        # grn_with_duplicates has TF1, TF2 (2 sources)
        assert n_s == 2
        # Unique edges: (TF1,Gene1), (TF1,Gene2), (TF2,Gene3) = 3
        assert n_e == 3

    def test_empty_grn(self, empty_grn):
        """Test with empty GRN."""
        n_s, n_c, n_t, n_e, n_r = _get_grn_stats(empty_grn)

        assert n_s == 0
        assert n_c == 0
        assert n_t == 0
        assert n_e == 0
        assert n_r == 0.0


# ============================================================================
# stats function tests
# ============================================================================


class TestStats:
    """Tests for stats function."""

    def test_returns_dataframe(self, grn_a):
        """Test that function returns a DataFrame."""
        result = stats(grn_a)

        assert isinstance(result, pd.DataFrame)

    def test_output_columns(self, grn_a):
        """Test that output has correct columns."""
        result = stats(grn_a)

        expected_columns = ["name", "n_sources", "n_cres", "n_targets", "n_edges", "mean_regulon_size"]
        assert list(result.columns) == expected_columns

    def test_single_dataframe_input(self, grn_a):
        """Test with single DataFrame input."""
        result = stats(grn_a)

        assert len(result) == 1
        assert result.iloc[0]["name"] == "grn"

    def test_dict_input(self, grn_a, grn_b):
        """Test with dictionary input."""
        result = stats({"GRN_A": grn_a, "GRN_B": grn_b})

        assert len(result) == 2
        assert set(result["name"]) == {"GRN_A", "GRN_B"}

    def test_correct_values_single_grn(self, grn_a):
        """Test correct statistics for single GRN."""
        result = stats(grn_a)

        row = result.iloc[0]
        assert row["n_sources"] == 3
        assert row["n_cres"] == 5
        assert row["n_targets"] == 5
        assert row["n_edges"] == 5
        assert row["mean_regulon_size"] == pytest.approx(5 / 3)

    def test_multiple_grns_correct_values(self, grn_a, grn_b, grn_c):
        """Test correct statistics for multiple GRNs."""
        result = stats({"A": grn_a, "B": grn_b, "C": grn_c})

        assert len(result) == 3

        # Check GRN A
        row_a = result[result["name"] == "A"].iloc[0]
        assert row_a["n_sources"] == 3
        assert row_a["n_edges"] == 5

        # Check GRN B
        row_b = result[result["name"] == "B"].iloc[0]
        assert row_b["n_sources"] == 3  # TF1, TF2, TF4
        assert row_b["n_edges"] == 5

        # Check GRN C
        row_c = result[result["name"] == "C"].iloc[0]
        assert row_c["n_sources"] == 2  # TF5, TF6
        assert row_c["n_edges"] == 3

    def test_grn_without_cre(self, grn_no_cre):
        """Test GRN without cre column."""
        result = stats(grn_no_cre)

        assert result.iloc[0]["n_cres"] == 0

    def test_empty_grn(self, empty_grn):
        """Test with empty GRN."""
        result = stats(empty_grn)

        row = result.iloc[0]
        assert row["n_sources"] == 0
        assert row["n_cres"] == 0
        assert row["n_targets"] == 0
        assert row["n_edges"] == 0
        assert row["mean_regulon_size"] == 0.0

    def test_grn_with_duplicates(self, grn_with_duplicates):
        """Test GRN with duplicate edges."""
        result = stats(grn_with_duplicates)

        row = result.iloc[0]
        # Unique edges after deduplication: 3
        assert row["n_edges"] == 3

    def test_preserves_dict_order(self):
        """Test that dictionary order is preserved in output."""
        grn1 = pd.DataFrame({"source": ["TF1"], "target": ["Gene1"]})
        grn2 = pd.DataFrame({"source": ["TF2"], "target": ["Gene2"]})
        grn3 = pd.DataFrame({"source": ["TF3"], "target": ["Gene3"]})

        result = stats({"First": grn1, "Second": grn2, "Third": grn3})

        assert list(result["name"]) == ["First", "Second", "Third"]
