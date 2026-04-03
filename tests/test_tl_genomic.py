"""Tests for gretapy.tl._genomic module."""

import numpy as np
import pandas as pd
import pyranges as pr

from gretapy.tl._genomic import _cre, _cre_column, _grn_to_pr, _peaks_to_pr


class TestGrnToPr:
    """Tests for _grn_to_pr function."""

    def test_basic_conversion(self, simple_grn):
        """Test basic GRN to PyRanges conversion."""
        result = _grn_to_pr(grn=simple_grn)

        assert isinstance(result, pr.PyRanges)
        assert "Chromosome" in result.df.columns
        assert "Start" in result.df.columns
        assert "End" in result.df.columns

    def test_with_column_parameter(self, simple_grn):
        """Test conversion with column parameter to extract Name."""
        result = _grn_to_pr(grn=simple_grn, column="source")

        assert "Name" in result.df.columns
        assert set(result.df["Name"]) == {"PAX5", "EBF1", "GATA3", "TCF7", "SPI1", "CEBPA"}

    def test_empty_grn(self):
        """Test with empty GRN returns None."""
        empty_grn = pd.DataFrame(columns=["source", "target", "cre"])
        result = _grn_to_pr(grn=empty_grn)

        assert result is None

    def test_parses_cre_format_correctly(self, simple_grn):
        """Test that CRE format chrX-start-end is parsed correctly."""
        result = _grn_to_pr(grn=simple_grn)

        # Check that coordinates are parsed correctly
        assert "chr16" in result.df["Chromosome"].values
        assert "chr11" in result.df["Chromosome"].values


class TestPeaksToPr:
    """Tests for _peaks_to_pr function."""

    def test_basic_conversion(self):
        """Test basic peaks to PyRanges conversion."""
        peaks = ["chr1-1000-2000", "chr2-3000-4000"]
        result = _peaks_to_pr(peaks=peaks)

        assert isinstance(result, pr.PyRanges)
        assert len(result) == 2

    def test_parses_coordinates(self):
        """Test that peak coordinates are parsed correctly."""
        peaks = ["chr16-28931000-28931500"]
        result = _peaks_to_pr(peaks=peaks)

        df = result.df
        assert df["Chromosome"].iloc[0] == "chr16"
        assert df["Start"].iloc[0] == 28931000
        assert df["End"].iloc[0] == 28931500

    def test_with_numpy_array(self):
        """Test that numpy array input works."""
        peaks = np.array(["chr1-100-200", "chr2-300-400"])
        result = _peaks_to_pr(peaks=peaks)

        assert isinstance(result, pr.PyRanges)
        assert len(result) == 2


class TestCreColumn:
    """Tests for _cre_column function."""

    def test_basic_functionality(self, simple_grn, pyranges_db):
        """Test basic CRE-target linkage evaluation."""
        genes = ["CD19", "MS4A1", "CD3E", "IL7R"]
        peaks = [
            "chr16-28931000-28931500",
            "chr11-60223000-60223500",
            "chr11-118209000-118209500",
            "chr5-35871000-35871500",
        ]

        prc, rcl, f01 = _cre_column(
            grn=simple_grn,
            db=pyranges_db,
            genes=genes,
            peaks=peaks,
            cats=None,
            column="target",
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_with_source_column(self, simple_grn, pyranges_db):
        """Test CRE evaluation using source column (TF binding)."""
        genes = ["PAX5", "GATA3", "CD19", "MS4A1", "CD3E", "IL7R"]
        peaks = [
            "chr16-28931000-28931500",
            "chr11-60223000-60223500",
            "chr11-118209000-118209500",
            "chr5-35871000-35871500",
        ]

        # Create a db that has TF names instead of gene names
        db = pr.PyRanges(
            pd.DataFrame(
                {
                    "Chromosome": ["chr16", "chr11"],
                    "Start": [28931100, 60223100],
                    "End": [28931400, 60223400],
                    "Name": ["PAX5", "PAX5"],
                    "Score": ["B cell", "T cell"],
                }
            )
        )

        prc, rcl, f01 = _cre_column(
            grn=simple_grn,
            db=db,
            genes=genes,
            peaks=peaks,
            cats=None,
            column="source",
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_with_category_filter(self, simple_grn, pyranges_db):
        """Test CRE evaluation with category filtering."""
        genes = ["CD19", "MS4A1", "CD3E", "IL7R"]
        peaks = [
            "chr16-28931000-28931500",
            "chr11-60223000-60223500",
            "chr11-118209000-118209500",
            "chr5-35871000-35871500",
        ]

        prc, rcl, f01 = _cre_column(
            grn=simple_grn,
            db=pyranges_db,
            genes=genes,
            peaks=peaks,
            cats=["B cell"],
            column="target",
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1


class TestCre:
    """Tests for _cre function."""

    def test_basic_functionality(self, simple_grn, pyranges_db):
        """Test basic CRE regions evaluation."""
        peaks = [
            "chr16-28931000-28931500",
            "chr11-60223000-60223500",
            "chr11-118209000-118209500",
            "chr5-35871000-35871500",
        ]

        prc, rcl, f01 = _cre(
            grn=simple_grn,
            db=pyranges_db,
            peaks=peaks,
            cats=None,
            reverse=False,
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_reverse_mode(self, simple_grn, pyranges_db):
        """Test CRE evaluation in reverse mode (e.g., for blacklist)."""
        peaks = [
            "chr16-28931000-28931500",
            "chr11-60223000-60223500",
            "chr11-118209000-118209500",
            "chr5-35871000-35871500",
        ]

        prc, rcl, f01 = _cre(
            grn=simple_grn,
            db=pyranges_db,
            peaks=peaks,
            cats=None,
            reverse=True,
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_with_category_filter(self, simple_grn, pyranges_db):
        """Test CRE evaluation with category filtering."""
        peaks = [
            "chr16-28931000-28931500",
            "chr11-60223000-60223500",
            "chr11-118209000-118209500",
            "chr5-35871000-35871500",
        ]

        prc, rcl, f01 = _cre(
            grn=simple_grn,
            db=pyranges_db,
            peaks=peaks,
            cats=["B cell", "T cell"],
            reverse=False,
        )

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_no_overlap(self):
        """Test with no overlap between GRN CREs and database."""
        grn = pd.DataFrame(
            {
                "source": ["TF1"],
                "target": ["Gene1"],
                "cre": ["chr1-1000-2000"],
                "score": [0.5],
            }
        )
        db = pr.PyRanges(
            pd.DataFrame(
                {
                    "Chromosome": ["chr99"],
                    "Start": [999999],
                    "End": [1000000],
                    "Name": ["Gene1"],
                    "Score": ["Type1"],
                }
            )
        )
        peaks = ["chr1-1000-2000"]

        prc, rcl, f01 = _cre(grn=grn, db=db, peaks=peaks, cats=None, reverse=False)

        # No overlap means no true positives
        assert prc == 0.0
        assert f01 == 0.0
