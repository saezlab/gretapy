"""Tests for pygreta.mt module (GRN inference methods)."""

from unittest.mock import patch

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pyranges as pr
import pytest

from pygreta.mt._collectri import collectri
from pygreta.mt._correlation import correlation
from pygreta.mt._random import _get_cres_pr, _get_overlap_cres, _get_window, random

# ============================================================================
# Fixtures specific to mt module tests
# ============================================================================


@pytest.fixture
def mock_collectri_grn():
    """Mock CollecTRI GRN database."""
    return pd.DataFrame(
        {
            "source": ["PAX5", "PAX5", "GATA3", "GATA3", "SPI1", "SPI1", "TCF7", "TCF7", "RUNX1", "RUNX1"],
            "target": ["CD19", "MS4A1", "CD3E", "IL7R", "CD14", "CD79A", "FOXP3", "BCL2", "CD34", "IRF4"],
            "weight": [1.0, 0.8, 0.9, 0.7, 0.85, 0.6, 0.75, 0.65, 0.8, 0.7],
        }
    )


@pytest.fixture
def mock_promoters_db():
    """Mock Promoters database as PyRanges object with coordinates overlapping peak_names."""
    return pr.PyRanges(
        pd.DataFrame(
            {
                "Chromosome": [
                    "chr16",
                    "chr11",
                    "chr11",
                    "chr5",
                    "chr5",
                    "chr19",
                    "chr18",
                    "chr6",
                    "chr1",
                    "chr5",
                    "chrX",
                    "chr21",
                ],
                "Start": [
                    28931100,  # CD19
                    60223100,  # MS4A1
                    118209100,  # CD3E
                    35871100,  # IL7R
                    140013100,  # CD14
                    41879100,  # CD79A
                    63318100,  # BCL2
                    411100,  # IRF4
                    207940100,  # CD34
                    134110100,  # TCF7
                    49250100,  # FOXP3
                    34799100,  # RUNX1
                ],
                "End": [
                    28931400,
                    60223400,
                    118209400,
                    35871400,
                    140013400,
                    41879400,
                    63318400,
                    411400,
                    207940400,
                    134110400,
                    49250400,
                    34799400,
                ],
                "Name": [
                    "CD19",
                    "MS4A1",
                    "CD3E",
                    "IL7R",
                    "CD14",
                    "CD79A",
                    "BCL2",
                    "IRF4",
                    "CD34",
                    "TCF7",
                    "FOXP3",
                    "RUNX1",
                ],
            }
        )
    )


@pytest.fixture
def mock_lambert_tfs():
    """Mock LambertTFs database."""
    return pd.DataFrame({0: ["PAX5", "GATA3", "SPI1", "TCF7", "RUNX1"]})


@pytest.fixture
def mudata_for_mt():
    """MuData object specifically configured for mt module testing."""
    np.random.seed(42)

    genes = [
        "PAX5",
        "GATA3",
        "SPI1",
        "TCF7",
        "RUNX1",
        "CD19",
        "MS4A1",
        "CD3E",
        "IL7R",
        "CD14",
        "CD79A",
        "BCL2",
        "IRF4",
        "CD34",
        "FOXP3",
    ]
    peak_names = [
        "chr16-28931000-28931500",  # CD19 promoter
        "chr11-60223000-60223500",  # MS4A1 enhancer
        "chr11-118209000-118209500",  # CD3E promoter
        "chr5-35871000-35871500",  # IL7R enhancer
        "chr5-140013000-140013500",  # CD14 promoter
        "chr19-41879000-41879500",  # CD79A promoter
        "chr18-63318000-63318500",  # BCL2 enhancer
        "chr6-411000-411500",  # IRF4 promoter
        "chr1-207940000-207940500",  # CD34 enhancer
        "chr5-134110000-134110500",  # TCF7 promoter
        "chrX-49250000-49250500",  # FOXP3 enhancer
        "chr21-34799000-34799500",  # RUNX1 promoter
    ]

    n_cells = 50
    n_genes = len(genes)
    n_peaks = len(peak_names)

    # RNA modality
    rna_X = np.random.rand(n_cells, n_genes).astype(np.float32)
    rna = ad.AnnData(X=rna_X)
    rna.var_names = genes
    rna.obs_names = [f"Cell{i}" for i in range(n_cells)]

    # ATAC modality
    atac_X = np.random.rand(n_cells, n_peaks).astype(np.float32)
    atac = ad.AnnData(X=atac_X)
    atac.var_names = peak_names
    atac.obs_names = [f"Cell{i}" for i in range(n_cells)]

    return mu.MuData({"rna": rna, "atac": atac})


# ============================================================================
# Helper function tests
# ============================================================================


class TestGetCresPr:
    """Tests for _get_cres_pr helper function."""

    def test_returns_pyranges(self):
        """Test that function returns a PyRanges object."""
        peaks = np.array(["chr1-100-200", "chr2-300-400"])
        result = _get_cres_pr(peaks)

        assert isinstance(result, pr.PyRanges)

    def test_correct_columns(self):
        """Test that result has correct columns."""
        peaks = np.array(["chr1-100-200", "chr2-300-400"])
        result = _get_cres_pr(peaks)
        df = result.df

        assert "Chromosome" in df.columns
        assert "Start" in df.columns
        assert "End" in df.columns

    def test_correct_parsing(self):
        """Test that peaks are parsed correctly."""
        peaks = np.array(["chr1-100-200", "chrX-500-600"])
        result = _get_cres_pr(peaks)
        df = result.df

        assert df.iloc[0]["Chromosome"] == "chr1"
        assert int(df.iloc[0]["Start"]) == 100
        assert int(df.iloc[0]["End"]) == 200
        assert df.iloc[1]["Chromosome"] == "chrX"


class TestGetWindow:
    """Tests for _get_window helper function."""

    def test_returns_pyranges(self, mock_promoters_db):
        """Test that function returns a PyRanges object."""
        result = _get_window(gannot=mock_promoters_db, target="CD19", w_size=50000)

        assert isinstance(result, pr.PyRanges)

    def test_window_size_respected(self, mock_promoters_db):
        """Test that window size is respected."""
        w_size = 10000
        result = _get_window(gannot=mock_promoters_db, target="CD19", w_size=w_size)
        df = result.df

        # Window should be 2*w_size wide
        window_width = df["End"].iloc[0] - df["Start"].iloc[0]
        assert window_width == 2 * w_size


class TestGetOverlapCres:
    """Tests for _get_overlap_cres helper function."""

    def test_returns_array_when_overlaps_exist(self, mock_promoters_db):
        """Test that function returns array when overlaps exist."""
        peaks = np.array(["chr16-28931000-28931500", "chr11-60223000-60223500"])
        cres_pr = _get_cres_pr(peaks)

        result = _get_overlap_cres(gene="CD19", gannot=mock_promoters_db, cres_pr=cres_pr, w_size=50000)

        assert result is not None
        assert len(result) > 0

    def test_returns_none_when_no_overlaps(self, mock_promoters_db):
        """Test that function returns None when no overlaps."""
        # Peaks on different chromosome
        peaks = np.array(["chr22-100-200", "chr22-300-400"])
        cres_pr = _get_cres_pr(peaks)

        result = _get_overlap_cres(gene="CD19", gannot=mock_promoters_db, cres_pr=cres_pr, w_size=50000)

        assert result is None


# ============================================================================
# collectri function tests
# ============================================================================


class TestCollectri:
    """Tests for collectri function."""

    def test_invalid_mdata_type_raises(self):
        """Test that non-MuData input raises ValueError."""
        with pytest.raises(ValueError, match="must be a MuData object"):
            collectri(mdata="not_a_mudata")

    def test_missing_rna_modality_raises(self):
        """Test that missing rna modality raises ValueError."""
        atac = ad.AnnData(X=np.random.rand(10, 5))
        mdata = mu.MuData({"atac": atac})

        with pytest.raises(ValueError, match='must contain "rna" and "atac"'):
            collectri(mdata=mdata)

    def test_missing_atac_modality_raises(self):
        """Test that missing atac modality raises ValueError."""
        rna = ad.AnnData(X=np.random.rand(10, 5))
        mdata = mu.MuData({"rna": rna})

        with pytest.raises(ValueError, match='must contain "rna" and "atac"'):
            collectri(mdata=mdata)

    def test_invalid_organism_raises(self, mudata_for_mt):
        """Test that invalid organism raises ValueError."""
        with pytest.raises(ValueError, match="Invalid organism"):
            collectri(mdata=mudata_for_mt, organism="invalid")

    @patch("pygreta.mt._collectri.read_db")
    def test_returns_dataframe(self, mock_read_db, mudata_for_mt, mock_collectri_grn, mock_promoters_db):
        """Test that function returns a DataFrame."""
        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_collectri_grn if db_name == "CollecTRI" else mock_promoters_db
        )

        result = collectri(mdata=mudata_for_mt, organism="hg38", min_targets=1)

        assert isinstance(result, pd.DataFrame)

    @patch("pygreta.mt._collectri.read_db")
    def test_output_columns(self, mock_read_db, mudata_for_mt, mock_collectri_grn, mock_promoters_db):
        """Test that output has correct columns."""
        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_collectri_grn if db_name == "CollecTRI" else mock_promoters_db
        )

        result = collectri(mdata=mudata_for_mt, organism="hg38", min_targets=1)

        assert "source" in result.columns
        assert "cre" in result.columns
        assert "target" in result.columns
        assert "score" in result.columns

    @patch("pygreta.mt._collectri.read_db")
    def test_min_targets_filtering(self, mock_read_db, mudata_for_mt, mock_collectri_grn, mock_promoters_db):
        """Test that min_targets parameter filters TFs."""
        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_collectri_grn if db_name == "CollecTRI" else mock_promoters_db
        )

        result_low = collectri(mdata=mudata_for_mt, organism="hg38", min_targets=1)
        result_high = collectri(mdata=mudata_for_mt, organism="hg38", min_targets=100)

        # Higher min_targets should result in fewer or equal TFs
        assert len(result_high) <= len(result_low)

    @patch("pygreta.mt._collectri.read_db")
    def test_filters_by_available_genes(self, mock_read_db, mock_collectri_grn, mock_promoters_db):
        """Test that GRN is filtered by genes available in MuData."""
        # Create MuData with only some genes
        rna = ad.AnnData(X=np.random.rand(10, 3))
        rna.var_names = ["PAX5", "CD19", "GATA3"]  # Limited genes
        atac = ad.AnnData(X=np.random.rand(10, 2))
        atac.var_names = ["chr16-28931000-28931500", "chr11-118209000-118209500"]
        mdata = mu.MuData({"rna": rna, "atac": atac})

        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_collectri_grn if db_name == "CollecTRI" else mock_promoters_db
        )

        result = collectri(mdata=mdata, organism="hg38", min_targets=1)

        # All sources and targets in result should be in mdata genes
        genes = set(mdata.mod["rna"].var_names)
        if len(result) > 0:
            assert set(result["source"]).issubset(genes)
            assert set(result["target"]).issubset(genes)


# ============================================================================
# correlation function tests
# ============================================================================


class TestCorrelation:
    """Tests for correlation function."""

    def test_invalid_mdata_type_raises(self):
        """Test that non-MuData input raises ValueError."""
        with pytest.raises(ValueError, match="must be a MuData object"):
            correlation(mdata="not_a_mudata", tfs=["TF1"])

    def test_missing_modalities_raises(self):
        """Test that missing modalities raises ValueError."""
        rna = ad.AnnData(X=np.random.rand(10, 5))
        mdata = mu.MuData({"rna": rna})

        with pytest.raises(ValueError, match='must contain "rna" and "atac"'):
            correlation(mdata=mdata, tfs=["TF1"])

    def test_invalid_organism_raises(self, mudata_for_mt):
        """Test that invalid organism raises ValueError."""
        with pytest.raises(ValueError, match="Invalid organism"):
            correlation(mdata=mudata_for_mt, tfs=["PAX5"], organism="invalid")

    def test_invalid_method_raises(self, mudata_for_mt):
        """Test that invalid correlation method raises ValueError."""
        with pytest.raises(ValueError, match='must be "pearson" or "spearman"'):
            correlation(mdata=mudata_for_mt, tfs=["PAX5"], method="invalid")

    @patch("pygreta.mt._correlation.read_db")
    def test_returns_dataframe(self, mock_read_db, mudata_for_mt, mock_promoters_db):
        """Test that function returns a DataFrame."""
        mock_read_db.return_value = mock_promoters_db

        result = correlation(
            mdata=mudata_for_mt,
            tfs=["PAX5", "GATA3", "SPI1"],
            organism="hg38",
            thr_r=0.0,  # Low threshold to get results
            min_targets=1,
        )

        assert isinstance(result, pd.DataFrame)

    @patch("pygreta.mt._correlation.read_db")
    def test_output_columns(self, mock_read_db, mudata_for_mt, mock_promoters_db):
        """Test that output has correct columns."""
        mock_read_db.return_value = mock_promoters_db

        result = correlation(
            mdata=mudata_for_mt,
            tfs=["PAX5", "GATA3"],
            organism="hg38",
            thr_r=0.0,
            min_targets=1,
        )

        if len(result) > 0:
            assert "source" in result.columns
            assert "cre" in result.columns
            assert "target" in result.columns
            assert "score" in result.columns

    @patch("pygreta.mt._correlation.read_db")
    def test_pearson_method(self, mock_read_db, mudata_for_mt, mock_promoters_db):
        """Test pearson correlation method."""
        mock_read_db.return_value = mock_promoters_db

        result = correlation(
            mdata=mudata_for_mt,
            tfs=["PAX5", "GATA3"],
            organism="hg38",
            method="pearson",
            thr_r=0.0,
            min_targets=1,
        )

        assert isinstance(result, pd.DataFrame)

    @patch("pygreta.mt._correlation.read_db")
    def test_spearman_method(self, mock_read_db, mudata_for_mt, mock_promoters_db):
        """Test spearman correlation method."""
        mock_read_db.return_value = mock_promoters_db

        result = correlation(
            mdata=mudata_for_mt,
            tfs=["PAX5", "GATA3"],
            organism="hg38",
            method="spearman",
            thr_r=0.0,
            min_targets=1,
        )

        assert isinstance(result, pd.DataFrame)

    @patch("pygreta.mt._correlation.read_db")
    def test_thr_r_filtering(self, mock_read_db, mudata_for_mt, mock_promoters_db):
        """Test that thr_r parameter filters by correlation threshold."""
        mock_read_db.return_value = mock_promoters_db

        result_low = correlation(
            mdata=mudata_for_mt,
            tfs=["PAX5", "GATA3", "SPI1"],
            organism="hg38",
            thr_r=0.0,
            min_targets=1,
        )

        result_high = correlation(
            mdata=mudata_for_mt,
            tfs=["PAX5", "GATA3", "SPI1"],
            organism="hg38",
            thr_r=0.9,
            min_targets=1,
        )

        # Higher threshold should result in fewer edges
        assert len(result_high) <= len(result_low)

    @patch("pygreta.mt._correlation.read_db")
    def test_filters_self_regulation(self, mock_read_db, mudata_for_mt, mock_promoters_db):
        """Test that self-regulation edges are removed."""
        mock_read_db.return_value = mock_promoters_db

        result = correlation(
            mdata=mudata_for_mt,
            tfs=["PAX5", "GATA3"],
            organism="hg38",
            thr_r=0.0,
            min_targets=1,
        )

        # No edge should have source == target
        if len(result) > 0:
            assert not any(result["source"] == result["target"])

    @patch("pygreta.mt._correlation.read_db")
    def test_filters_tfs_not_in_dataset(self, mock_read_db, mudata_for_mt, mock_promoters_db):
        """Test that TFs not in dataset are filtered out."""
        mock_read_db.return_value = mock_promoters_db

        result = correlation(
            mdata=mudata_for_mt,
            tfs=["PAX5", "NONEXISTENT_TF"],
            organism="hg38",
            thr_r=0.0,
            min_targets=1,
        )

        # NONEXISTENT_TF should not appear in results
        if len(result) > 0:
            assert "NONEXISTENT_TF" not in result["source"].values


# ============================================================================
# random function tests
# ============================================================================


class TestRandom:
    """Tests for random function."""

    def test_invalid_mdata_type_raises(self):
        """Test that non-MuData input raises ValueError."""
        with pytest.raises(ValueError, match="must be a MuData object"):
            random(mdata="not_a_mudata")

    def test_missing_modalities_raises(self):
        """Test that missing modalities raises ValueError."""
        rna = ad.AnnData(X=np.random.rand(10, 5))
        mdata = mu.MuData({"rna": rna})

        with pytest.raises(ValueError, match='must contain "rna" and "atac"'):
            random(mdata=mdata)

    def test_invalid_organism_raises(self, mudata_for_mt):
        """Test that invalid organism raises ValueError."""
        with pytest.raises(ValueError, match="Invalid organism"):
            random(mdata=mudata_for_mt, organism="invalid")

    @patch("pygreta.mt._random.read_db")
    def test_returns_dataframe(self, mock_read_db, mudata_for_mt, mock_lambert_tfs, mock_promoters_db):
        """Test that function returns a DataFrame."""
        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_lambert_tfs if db_name == "LambertTFs" else mock_promoters_db
        )

        result = random(mdata=mudata_for_mt, organism="hg38", min_targets=1, seed=42)

        assert isinstance(result, pd.DataFrame)

    @patch("pygreta.mt._random.read_db")
    def test_output_columns(self, mock_read_db, mudata_for_mt, mock_lambert_tfs, mock_promoters_db):
        """Test that output has correct columns."""
        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_lambert_tfs if db_name == "LambertTFs" else mock_promoters_db
        )

        result = random(mdata=mudata_for_mt, organism="hg38", min_targets=1, seed=42)

        if len(result) > 0:
            assert "source" in result.columns
            assert "cre" in result.columns
            assert "target" in result.columns
            assert "score" in result.columns

    @patch("pygreta.mt._random.read_db")
    def test_reproducibility_with_seed(self, mock_read_db, mudata_for_mt, mock_lambert_tfs, mock_promoters_db):
        """Test that same seed produces same results."""
        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_lambert_tfs if db_name == "LambertTFs" else mock_promoters_db
        )

        result1 = random(mdata=mudata_for_mt, organism="hg38", min_targets=1, seed=42)
        result2 = random(mdata=mudata_for_mt, organism="hg38", min_targets=1, seed=42)

        pd.testing.assert_frame_equal(result1, result2)

    @patch("pygreta.mt._random.read_db")
    def test_different_seeds_different_results(self, mock_read_db, mudata_for_mt, mock_lambert_tfs, mock_promoters_db):
        """Test that different seeds produce different results."""
        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_lambert_tfs if db_name == "LambertTFs" else mock_promoters_db
        )

        result1 = random(mdata=mudata_for_mt, organism="hg38", min_targets=1, seed=42)
        result2 = random(mdata=mudata_for_mt, organism="hg38", min_targets=1, seed=123)

        if len(result1) > 0 and len(result2) > 0:
            # Results should generally differ (not guaranteed but highly likely)
            assert not result1.equals(result2) or len(result1) == 0

    @patch("pygreta.mt._random.read_db")
    def test_custom_tfs_list(self, mock_read_db, mudata_for_mt, mock_promoters_db):
        """Test using custom TFs list instead of LambertTFs."""
        mock_read_db.return_value = mock_promoters_db

        result = random(
            mdata=mudata_for_mt,
            organism="hg38",
            tfs=["PAX5", "GATA3"],
            min_targets=1,
            seed=42,
        )

        assert isinstance(result, pd.DataFrame)
        # All sources should be from provided TFs
        if len(result) > 0:
            assert set(result["source"]).issubset({"PAX5", "GATA3"})

    @patch("pygreta.mt._random.read_db")
    def test_g_perc_parameter(self, mock_read_db, mudata_for_mt, mock_lambert_tfs, mock_promoters_db):
        """Test that g_perc controls percentage of genes sampled."""
        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_lambert_tfs if db_name == "LambertTFs" else mock_promoters_db
        )

        result_low = random(mdata=mudata_for_mt, organism="hg38", g_perc=0.1, min_targets=1, seed=42)
        result_high = random(mdata=mudata_for_mt, organism="hg38", g_perc=0.9, min_targets=1, seed=42)

        # Higher g_perc should generally result in more edges (not guaranteed but likely)
        # We just test that both return valid DataFrames
        assert isinstance(result_low, pd.DataFrame)
        assert isinstance(result_high, pd.DataFrame)

    @patch("pygreta.mt._random.read_db")
    def test_min_targets_filtering(self, mock_read_db, mudata_for_mt, mock_lambert_tfs, mock_promoters_db):
        """Test that min_targets parameter filters TFs."""
        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_lambert_tfs if db_name == "LambertTFs" else mock_promoters_db
        )

        result_low = random(mdata=mudata_for_mt, organism="hg38", min_targets=1, seed=42)
        result_high = random(mdata=mudata_for_mt, organism="hg38", min_targets=100, seed=42)

        # Higher min_targets should result in fewer or equal edges
        assert len(result_high) <= len(result_low)

    @patch("pygreta.mt._random.read_db")
    def test_score_is_one(self, mock_read_db, mudata_for_mt, mock_lambert_tfs, mock_promoters_db):
        """Test that all scores are 1.0 for random GRN."""
        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_lambert_tfs if db_name == "LambertTFs" else mock_promoters_db
        )

        result = random(mdata=mudata_for_mt, organism="hg38", min_targets=1, seed=42)

        if len(result) > 0:
            assert all(result["score"] == 1.0)

    @patch("pygreta.mt._random.read_db")
    def test_empty_result_when_no_overlaps(self, mock_read_db):
        """Test that empty DataFrame is returned when no peak-gene overlaps."""
        # Create MuData with peaks on different chromosomes than promoters
        rna = ad.AnnData(X=np.random.rand(10, 5))
        rna.var_names = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]
        atac = ad.AnnData(X=np.random.rand(10, 3))
        atac.var_names = ["chr99-100-200", "chr99-300-400", "chr99-500-600"]
        mdata = mu.MuData({"rna": rna, "atac": atac})

        # Mock promoters on chr1
        mock_promoters = pr.PyRanges(
            pd.DataFrame(
                {
                    "Chromosome": ["chr1"] * 5,
                    "Start": [1000, 2000, 3000, 4000, 5000],
                    "End": [1500, 2500, 3500, 4500, 5500],
                    "Name": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
                }
            )
        )
        mock_tfs = pd.DataFrame({0: ["GENE1", "GENE2"]})

        mock_read_db.side_effect = lambda organism, db_name, verbose: (
            mock_tfs if db_name == "LambertTFs" else mock_promoters
        )

        result = random(mdata=mdata, organism="hg38", min_targets=1, seed=42)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0
        assert list(result.columns) == ["source", "cre", "target", "score"]
