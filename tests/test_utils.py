"""Tests for gretapy._utils module (shared utility functions)."""

from unittest.mock import MagicMock, mock_open, patch

import pandas as pd
import pyranges as pr
import pytest

import gretapy as gt


class TestShowOrganisms:
    """Tests for show_organisms function."""

    def test_returns_list(self):
        """Test that show_organisms returns a list."""
        result = gt.show_organisms()
        assert isinstance(result, list)

    def test_contains_expected_organisms(self):
        """Test that result contains expected organisms."""
        result = gt.show_organisms()
        assert "hg38" in result
        assert "mm10" in result


class TestShowDatasets:
    """Tests for show_datasets function."""

    def test_returns_dataframe(self):
        """Test that show_datasets returns a DataFrame."""
        result = gt.show_datasets(organism="hg38")
        assert isinstance(result, pd.DataFrame)

    def test_dataframe_columns(self):
        """Test that DataFrame has expected columns."""
        result = gt.show_datasets(organism="hg38")
        expected_cols = {"name", "pubmed", "geo"}
        assert expected_cols == set(result.columns)

    def test_invalid_organism_raises(self):
        """Test that invalid organism raises assertion error."""
        with pytest.raises(AssertionError):
            gt.show_datasets(organism="invalid_organism")


class TestShowTerms:
    """Tests for show_terms function."""

    @patch("gretapy._utils.os.path.isfile", return_value=True)
    @patch("gretapy._utils.pd.read_csv")
    def test_returns_dataframe(self, mock_read_csv, mock_isfile):
        """Test that show_terms returns a DataFrame (file already cached)."""
        # Use side_effect so each call returns a fresh DataFrame (avoids "organism already exists")
        mock_read_csv.side_effect = lambda *a, **kw: pd.DataFrame({"db_name": ["HPA"], "term": ["B cell"]})
        result = gt.show_terms(organism="hg38")
        assert isinstance(result, pd.DataFrame)

    @patch("gretapy._utils.pd.read_csv")
    @patch("gretapy._utils.shutil.copyfileobj")
    @patch("gretapy._utils._download")
    @patch("gretapy._utils.os.path.isfile", return_value=False)
    def test_downloads_when_file_missing(self, mock_isfile, mock_download, mock_copy, mock_read_csv):
        """Test that show_terms downloads when file is not cached."""
        mock_data = MagicMock()
        mock_download.return_value = mock_data
        mock_read_csv.side_effect = lambda *a, **kw: pd.DataFrame({"db_name": ["HPA"], "term": ["B cell"]})
        with patch("builtins.open", mock_open()):
            result = gt.show_terms(organism="hg38")
        assert isinstance(result, pd.DataFrame)
        assert mock_download.called

    def test_invalid_organism_raises(self):
        """Test that invalid organism raises assertion error."""
        with pytest.raises(AssertionError):
            gt.show_terms(organism="invalid_organism")


class TestShowGenomeAnnotation:
    """Tests for show_genome_annotation function."""

    @patch("gretapy._utils.os.path.isfile", return_value=True)
    @patch("gretapy._utils.pr.read_bed")
    def test_returns_pyranges(self, mock_read_bed, mock_isfile):
        """Test that show_genome_annotation returns a PyRanges object (file cached)."""
        mock_read_bed.return_value = pr.PyRanges()
        result = gt.show_genome_annotation(organism="hg38")
        assert isinstance(result, pr.PyRanges)

    @patch("gretapy._utils.pr.read_bed")
    @patch("gretapy._utils.shutil.copyfileobj")
    @patch("gretapy._utils._download")
    @patch("gretapy._utils.os.path.isfile", return_value=False)
    def test_downloads_when_file_missing(self, mock_isfile, mock_download, mock_copy, mock_read_bed):
        """Test that show_genome_annotation downloads when file is not cached."""
        mock_data = MagicMock()
        mock_download.return_value = mock_data
        mock_read_bed.return_value = pr.PyRanges()
        with patch("builtins.open", mock_open()):
            result = gt.show_genome_annotation(organism="hg38")
        assert isinstance(result, pr.PyRanges)
        assert mock_download.called

    def test_invalid_organism_raises(self):
        """Test that invalid organism raises assertion error."""
        with pytest.raises(AssertionError):
            gt.show_genome_annotation(organism="invalid_organism")


class TestShowMetrics:
    """Tests for show_metrics function."""

    def test_returns_dataframe(self):
        """Test that show_metrics returns a DataFrame."""
        result = gt.show_metrics()
        assert isinstance(result, pd.DataFrame)

    def test_dataframe_columns(self):
        """Test that DataFrame has expected columns."""
        result = gt.show_metrics()
        expected_cols = {"organism", "category", "metric", "db"}
        assert expected_cols == set(result.columns)

    def test_filter_by_organism_hg38(self):
        """Test filtering by hg38 organism."""
        result = gt.show_metrics(organism="hg38")
        assert isinstance(result, pd.DataFrame)
        # When filtered, organism column is dropped
        assert "organism" not in result.columns
        assert len(result) > 0

    def test_filter_by_organism_mm10(self):
        """Test filtering by mm10 organism."""
        result = gt.show_metrics(organism="mm10")
        assert isinstance(result, pd.DataFrame)
        assert "organism" not in result.columns
        assert len(result) > 0

    def test_invalid_organism_raises(self):
        """Test that invalid organism raises assertion error."""
        with pytest.raises(AssertionError):
            gt.show_metrics(organism="invalid_organism")

    def test_contains_expected_categories(self):
        """Test that result contains expected metric categories."""
        result = gt.show_metrics()
        categories = result["category"].unique()
        expected_categories = {"Literature", "Genomic", "Predictive", "Mechanistic"}
        assert expected_categories.issubset(set(categories))
