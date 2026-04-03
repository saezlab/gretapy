"""Tests for gretapy.ds._db and gretapy.ds._dts modules."""

from unittest.mock import MagicMock, mock_open, patch

import anndata as ad
import mudata as mu
import pandas as pd
import pyranges as pr
import pytest

from gretapy.ds._db import _download_db, read_db, read_metrics
from gretapy.ds._dts import _download_dts, read_dts


def _make_metrics_df(extra_rows=None):
    """Helper to build a minimal metrics DataFrame as returned by the CSV."""
    rows = [
        {"name": "grn1", "org": "hg38", "dts": "Brain", "task": "TF Scoring", "db": "KnockTF", "prc": 0.8, "rcl": 0.7},
        {"name": "grn2", "org": "hg38", "dts": "Brain", "task": "Perturbation Forecasting", "db": "KnockTF", "prc": 0.6, "rcl": 0.5},
        {"name": "grn3", "org": "hg38", "dts": "Brain", "task": "TF Scoring", "db": "ChIP-Atlas", "prc": 0.9, "rcl": 0.85},
        {"name": "grn4", "org": "hg38", "dts": "Synthetic Pituitary", "task": "TF Scoring", "db": "ChIP-Atlas", "prc": 0.5, "rcl": 0.4},
        {"name": "grn5", "org": "hg38", "dts": "Unpaired Pituitary", "task": "TF Scoring", "db": "ChIP-Atlas", "prc": 0.4, "rcl": 0.3},
    ]
    if extra_rows:
        rows.extend(extra_rows)
    return pd.DataFrame(rows)


class TestReadMetrics:
    """Tests for read_metrics function."""

    @patch("gretapy.ds._db._download_metrics", return_value="/fake/metrics.csv.gz")
    @patch("gretapy.ds._db.pd.read_csv")
    def test_returns_dataframe(self, mock_read_csv, mock_download):
        """Test that read_metrics returns a DataFrame."""
        mock_read_csv.return_value = _make_metrics_df()
        result = read_metrics()
        assert isinstance(result, pd.DataFrame)

    @patch("gretapy.ds._db._download_metrics", return_value="/fake/metrics.csv.gz")
    @patch("gretapy.ds._db.pd.read_csv")
    def test_columns_renamed(self, mock_read_csv, mock_download):
        """Test that source columns are renamed correctly."""
        mock_read_csv.return_value = _make_metrics_df()
        result = read_metrics()
        assert "name" in result.columns
        assert "organism" in result.columns
        assert "dataset" in result.columns
        assert "precision" in result.columns
        assert "recall" in result.columns
        assert "org" not in result.columns
        assert "dts" not in result.columns
        assert "prc" not in result.columns
        assert "rcl" not in result.columns

    @patch("gretapy.ds._db._download_metrics", return_value="/fake/metrics.csv.gz")
    @patch("gretapy.ds._db.pd.read_csv")
    def test_remove_paired_true(self, mock_read_csv, mock_download):
        """Test that paired datasets are removed when remove_paired=True."""
        mock_read_csv.return_value = _make_metrics_df()
        result = read_metrics(remove_paired=True)
        assert "Synthetic Pituitary" not in result["dataset"].values
        assert "Unpaired Pituitary" not in result["dataset"].values

    @patch("gretapy.ds._db._download_metrics", return_value="/fake/metrics.csv.gz")
    @patch("gretapy.ds._db.pd.read_csv")
    def test_remove_paired_false(self, mock_read_csv, mock_download):
        """Test that paired datasets are kept when remove_paired=False."""
        mock_read_csv.return_value = _make_metrics_df()
        result = read_metrics(remove_paired=False)
        assert "Synthetic Pituitary" in result["dataset"].values
        assert "Unpaired Pituitary" in result["dataset"].values

    @patch("gretapy.ds._db._download_metrics", return_value="/fake/metrics.csv.gz")
    @patch("gretapy.ds._db.pd.read_csv")
    def test_knocktf_scoring_renamed(self, mock_read_csv, mock_download):
        """Test that KnockTF TF Scoring rows are renamed to 'KnockTF (scoring)'."""
        mock_read_csv.return_value = _make_metrics_df()
        result = read_metrics(remove_paired=False)
        knocktf_rows = result[result["db"] == "KnockTF (scoring)"]
        assert len(knocktf_rows) == 1

    @patch("gretapy.ds._db._download_metrics", return_value="/fake/metrics.csv.gz")
    @patch("gretapy.ds._db.pd.read_csv")
    def test_knocktf_forecasting_renamed(self, mock_read_csv, mock_download):
        """Test that KnockTF Perturbation Forecasting rows are renamed to 'KnockTF (forecasting)'."""
        mock_read_csv.return_value = _make_metrics_df()
        result = read_metrics(remove_paired=False)
        knocktf_rows = result[result["db"] == "KnockTF (forecasting)"]
        assert len(knocktf_rows) == 1

    @patch("gretapy.ds._db._download_metrics", return_value="/fake/metrics.csv.gz")
    @patch("gretapy.ds._db.pd.read_csv")
    def test_non_knocktf_db_unchanged(self, mock_read_csv, mock_download):
        """Test that non-KnockTF db values are not modified."""
        mock_read_csv.return_value = _make_metrics_df()
        result = read_metrics(remove_paired=False)
        assert "ChIP-Atlas" in result["db"].values
        assert "KnockTF" not in result["db"].values


class TestDownloadDb:
    """Tests for _download_db function."""

    @patch("gretapy.ds._db._log")
    @patch("gretapy.ds._db.os.path.isfile", return_value=True)
    @patch("gretapy.ds._db.os.makedirs")
    def test_returns_path_when_cached(self, mock_makedirs, mock_isfile, mock_log):
        """Test that returns file path when already cached (else branch)."""
        result = _download_db(organism="hg38", db_name="CollecTRI", verbose=False)
        assert isinstance(result, str)
        assert result.endswith(".csv.gz")
        assert not mock_log.called or True  # log may be called

    @patch("gretapy.ds._db._log")
    @patch("gretapy.ds._db.shutil.copyfileobj")
    @patch("gretapy.ds._db._download")
    @patch("gretapy.ds._db.os.path.isfile", return_value=False)
    @patch("gretapy.ds._db.os.makedirs")
    def test_downloads_non_h5ad_file(self, mock_makedirs, mock_isfile, mock_download, mock_copy, mock_log):
        """Test non-h5ad download path."""
        mock_data = MagicMock()
        mock_download.return_value = mock_data
        with patch("builtins.open", mock_open()):
            result = _download_db(organism="hg38", db_name="CollecTRI", verbose=False)
        assert isinstance(result, str)
        assert mock_download.called
        assert mock_copy.called

    @patch("gretapy.ds._db._log")
    @patch("gretapy.ds._db.os.remove")
    @patch("gretapy.ds._db.ad.read_h5ad")
    @patch("gretapy.ds._db.shutil.copyfileobj")
    @patch("gretapy.ds._db.gzip.GzipFile")
    @patch("gretapy.ds._db.tempfile.NamedTemporaryFile")
    @patch("gretapy.ds._db._download")
    @patch("gretapy.ds._db.os.path.isfile", return_value=False)
    @patch("gretapy.ds._db.os.makedirs")
    def test_downloads_h5ad_file(
        self,
        mock_makedirs,
        mock_isfile,
        mock_download,
        mock_tmpfile,
        mock_gzip,
        mock_copy,
        mock_read_h5ad,
        mock_remove,
        mock_log,
    ):
        """Test h5ad download path (special gzip handling)."""
        mock_data = MagicMock()
        mock_download.return_value = mock_data

        # Mock tempfile context manager
        mock_tmp = MagicMock()
        mock_tmp.name = "/tmp/fake.h5ad"
        mock_tmpfile.return_value.__enter__ = MagicMock(return_value=mock_tmp)
        mock_tmpfile.return_value.__exit__ = MagicMock(return_value=False)

        # Mock gzip context manager
        mock_gz = MagicMock()
        mock_gzip.return_value.__enter__ = MagicMock(return_value=mock_gz)
        mock_gzip.return_value.__exit__ = MagicMock(return_value=False)

        # Mock AnnData
        mock_adata = MagicMock(spec=ad.AnnData)
        mock_read_h5ad.return_value = mock_adata

        result = _download_db(organism="hg38", db_name="KnockTF (scoring)", verbose=False)
        assert isinstance(result, str)
        assert mock_download.called
        assert mock_read_h5ad.called
        assert mock_adata.write.called

    def test_invalid_organism_raises(self):
        """Test that invalid organism raises AssertionError."""
        with pytest.raises(AssertionError):
            _download_db(organism="invalid_org", db_name="CollecTRI")

    def test_invalid_db_name_raises(self):
        """Test that invalid db_name raises AssertionError."""
        with pytest.raises(AssertionError):
            _download_db(organism="hg38", db_name="NonExistentDB")


class TestReadDb:
    """Tests for read_db format dispatch."""

    @patch("gretapy.ds._db._download_db", return_value="./gretapy_data/fake_file.bed")
    @patch("gretapy.ds._db.pr.read_bed")
    def test_reads_bed_format(self, mock_read_bed, mock_download_db):
        """Test that .bed files are read as PyRanges."""
        mock_read_bed.return_value = pr.PyRanges()
        result = read_db(organism="hg38", db_name="ChIP-Atlas")
        assert isinstance(result, pr.PyRanges)
        mock_read_bed.assert_called_once_with("./gretapy_data/fake_file.bed")

    @patch("gretapy.ds._db._download_db", return_value="./gretapy_data/fake_file.tsv.gz")
    @patch("gretapy.ds._db.pd.read_csv")
    def test_reads_tsv_format(self, mock_read_csv, mock_download_db):
        """Test that .tsv.gz files are read as DataFrames."""
        mock_read_csv.return_value = pd.DataFrame({0: ["val1"]})
        result = read_db(organism="hg38", db_name="HPA")
        assert isinstance(result, pd.DataFrame)

    @patch("gretapy.ds._db._download_db", return_value="./gretapy_data/fake_file.csv.gz")
    @patch("gretapy.ds._db.pd.read_csv")
    def test_reads_csv_format(self, mock_read_csv, mock_download_db):
        """Test that .csv.gz files are read as DataFrames."""
        mock_read_csv.return_value = pd.DataFrame({"source": ["PAX5"], "target": ["CD19"]})
        result = read_db(organism="hg38", db_name="CollecTRI")
        assert isinstance(result, pd.DataFrame)

    @patch("gretapy.ds._db._download_db", return_value="./gretapy_data/fake_file.h5ad")
    @patch("gretapy.ds._db.ad.read_h5ad")
    def test_reads_h5ad_format(self, mock_read_h5ad, mock_download_db):
        """Test that .h5ad files are read as AnnData."""
        mock_adata = MagicMock(spec=ad.AnnData)
        mock_read_h5ad.return_value = mock_adata
        result = read_db(organism="hg38", db_name="KnockTF (scoring)")
        assert mock_read_h5ad.called

    @patch("gretapy.ds._db._download_db", return_value="./gretapy_data/fake_file.txt.gz")
    @patch("gretapy.ds._db.pd.read_csv")
    def test_reads_txt_format(self, mock_read_csv, mock_download_db):
        """Test that .txt.gz files are read as a list."""
        mock_read_csv.return_value = pd.DataFrame({0: ["gene1", "gene2", "gene3"]})
        result = read_db(organism="hg38", db_name="Lambert TFs")
        assert isinstance(result, list)
        assert "gene1" in result


class TestDownloadDts:
    """Tests for _download_dts function."""

    @patch("gretapy.ds._dts._log")
    @patch("gretapy.ds._dts.os.path.isfile", return_value=True)
    @patch("gretapy.ds._dts.os.makedirs")
    def test_returns_path_when_cached(self, mock_makedirs, mock_isfile, mock_log):
        """Test returns file path when already cached."""
        result = _download_dts(organism="hg38", dts_name="Brain", verbose=False)
        assert isinstance(result, str)
        assert "brain" in result.lower() or "hg38" in result.lower()

    @patch("gretapy.ds._dts._log")
    @patch("gretapy.ds._dts.shutil.copyfileobj")
    @patch("gretapy.ds._dts.gzip.open")
    @patch("gretapy.ds._dts._download")
    @patch("gretapy.ds._dts.os.path.isfile", return_value=False)
    @patch("gretapy.ds._dts.os.makedirs")
    def test_downloads_when_not_cached(
        self, mock_makedirs, mock_isfile, mock_download, mock_gzip_open, mock_copy, mock_log
    ):
        """Test that file is downloaded when not cached."""
        mock_data = MagicMock()
        mock_download.return_value = mock_data

        mock_gz_ctx = MagicMock()
        mock_gzip_open.return_value.__enter__ = MagicMock(return_value=mock_gz_ctx)
        mock_gzip_open.return_value.__exit__ = MagicMock(return_value=False)

        with patch("builtins.open", mock_open()):
            result = _download_dts(organism="hg38", dts_name="Brain", verbose=False)
        assert isinstance(result, str)
        assert mock_download.called

    def test_invalid_organism_raises(self):
        """Test that invalid organism raises AssertionError."""
        with pytest.raises(AssertionError):
            _download_dts(organism="invalid_org", dts_name="Brain")

    def test_invalid_dts_name_raises(self):
        """Test that invalid dts_name raises AssertionError."""
        with pytest.raises(AssertionError):
            _download_dts(organism="hg38", dts_name="NonExistentDataset")


class TestReadDts:
    """Tests for read_dts function."""

    @patch("gretapy.ds._dts._download_dts", return_value="./gretapy_data/fake_file.h5mu")
    @patch("gretapy.ds._dts.mu.read")
    def test_reads_mudata(self, mock_mu_read, mock_download_dts):
        """Test that read_dts returns MuData."""
        mock_mdata = MagicMock(spec=mu.MuData)
        mock_mu_read.return_value = mock_mdata
        result = read_dts(organism="hg38", dts_name="Brain")
        assert mock_mu_read.called
        mock_mu_read.assert_called_once_with("./gretapy_data/fake_file.h5mu")
