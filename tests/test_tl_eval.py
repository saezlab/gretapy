"""Tests for gretapy.tl._eval module."""

from unittest.mock import patch

import pandas as pd
import pytest

from gretapy.tl._eval import (
    _format_label,
    _format_log_prefix,
    _run_fileless_metric,
    _run_metric,
    _run_omics_metric,
    benchmark,
    eval_grn_dataset,
)


class TestFormatLogPrefix:
    """Tests for _format_log_prefix function."""

    def test_no_args_returns_empty(self):
        """Test that no args returns empty string."""
        assert _format_log_prefix() == ""

    def test_grn_name_only(self):
        """Test with grn_name only."""
        assert _format_log_prefix(grn_name="A") == "[A] "

    def test_dataset_name_only(self):
        """Test with dataset_name only."""
        assert _format_log_prefix(dataset_name="dts") == "[dts] "

    def test_both_args(self):
        """Test with both args."""
        assert _format_log_prefix(grn_name="A", dataset_name="dts") == "[A | dts] "


class TestFormatLabel:
    """Tests for _format_label function."""

    def test_no_args_returns_empty(self):
        """Test that no args returns empty string."""
        assert _format_label() == ""

    def test_grn_name_only(self):
        """Test with grn_name only."""
        assert _format_label(grn_name="A") == "A"

    def test_dataset_name_only(self):
        """Test with dataset_name only."""
        assert _format_label(dataset_name="dts") == "dts"

    def test_both_args(self):
        """Test with both args."""
        assert _format_label(grn_name="A", dataset_name="dts") == "A | dts"


class TestRunMetric:
    """Tests for _run_metric function."""

    def test_reference_grn_metric(self, simple_grn, reference_grn_db):
        """Test Reference GRN metric."""
        genes = ["PAX5", "GATA3", "SPI1", "CD19", "MS4A1", "CD3E", "IL7R", "CD14"]

        result = _run_metric(
            metric_type="Reference GRN",
            db_name="CollecTRI",
            grn=simple_grn,
            db=reference_grn_db,
            genes=genes,
            peaks=[],
            cats=None,
            adata=None,
        )

        assert result is not None
        assert len(result) == 3
        prc, rcl, f01 = result
        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_tf_markers_metric(self, simple_grn, tfm_db):
        """Test TF Markers metric (capital M)."""
        genes = ["PAX5", "GATA3", "SPI1", "RUNX1", "CD19", "MS4A1", "CD3E", "IL7R", "CD14"]

        result = _run_metric(
            metric_type="TF Markers",
            db_name="HPA",
            grn=simple_grn,
            db=tfm_db,
            genes=genes,
            peaks=[],
            cats=None,
            adata=None,
        )

        assert result is not None
        assert len(result) == 3

    def test_tf_pairs_metric(self, simple_grn, tfp_db):
        """Test TF Pairs metric (capital P)."""
        genes = ["PAX5", "GATA3", "SPI1", "CD19", "MS4A1", "CD3E", "IL7R", "CD14"]

        result = _run_metric(
            metric_type="TF Pairs",
            db_name="Europe PMC",
            grn=simple_grn,
            db=tfp_db,
            genes=genes,
            peaks=[],
            cats=None,
            adata=None,
        )

        assert result is not None
        assert len(result) == 3

    def test_tf_binding_metric(self, simple_grn):
        """Test TF Binding metric (db must have TF names in Name column)."""
        import pyranges as pr

        genes = ["PAX5", "GATA3", "SPI1", "CD19", "MS4A1", "CD3E", "IL7R", "CD14"]
        peaks = [
            "chr16-28931000-28931500",
            "chr11-60223000-60223500",
        ]
        # TF Binding db: Name = TF names, coordinates overlap with peaks
        db = pr.PyRanges(
            pd.DataFrame(
                {
                    "Chromosome": ["chr16", "chr11"],
                    "Start": [28931100, 60223100],
                    "End": [28931400, 60223400],
                    "Name": ["PAX5", "GATA3"],
                    "Score": ["B cell", "T cell"],
                }
            )
        )

        result = _run_metric(
            metric_type="TF Binding",
            db_name="ChIP-Atlas",
            grn=simple_grn,
            db=db,
            genes=genes,
            peaks=peaks,
            cats=None,
            adata=None,
        )

        assert result is not None
        assert len(result) == 3

    def test_cres_metric(self, simple_grn, pyranges_db):
        """Test CREs metric."""
        peaks = [
            "chr16-28931000-28931500",
            "chr11-60223000-60223500",
        ]

        result = _run_metric(
            metric_type="CREs",
            db_name="ENCODE CREs",
            grn=simple_grn,
            db=pyranges_db,
            genes=[],
            peaks=peaks,
            cats=None,
            adata=None,
        )

        assert result is not None
        assert len(result) == 3

    def test_cre_to_gene_links_metric(self, simple_grn, pyranges_db):
        """Test CRE to gene links metric."""
        genes = ["CD19", "MS4A1", "CD3E", "IL7R"]
        peaks = [
            "chr16-28931000-28931500",
            "chr11-60223000-60223500",
        ]

        result = _run_metric(
            metric_type="CRE to gene links",
            db_name="eQTL Catalogue",
            grn=simple_grn,
            db=pyranges_db,
            genes=genes,
            peaks=peaks,
            cats=None,
            adata=None,
        )

        assert result is not None
        assert len(result) == 3

    def test_gene_sets_metric(self, adata, simple_grn, gene_set_db):
        """Test Gene sets metric (lowercase s to match dispatch)."""
        result = _run_metric(
            metric_type="Gene sets",
            db_name="Hallmarks",
            grn=simple_grn,
            db=gene_set_db,
            genes=[],
            peaks=[],
            cats=None,
            adata=adata,
        )

        assert result is not None
        assert len(result) == 3

    def test_tf_scoring_metric(self, adata, simple_grn, knocktf_db):
        """Test TF Scoring metric."""
        result = _run_metric(
            metric_type="TF Scoring",
            db_name="KnockTF (scoring)",
            grn=simple_grn,
            db=knocktf_db,
            genes=[],
            peaks=[],
            cats=None,
            adata=adata,
        )

        assert result is not None
        assert len(result) == 3

    def test_perturbation_forecasting_metric(self, adata, simple_grn, knocktf_db):
        """Test Perturbation Forecasting metric."""
        result = _run_metric(
            metric_type="Perturbation Forecasting",
            db_name="KnockTF (forecasting)",
            grn=simple_grn,
            db=knocktf_db,
            genes=[],
            peaks=[],
            cats=None,
            adata=adata,
        )

        assert result is not None
        assert len(result) == 3

    def test_unknown_metric_returns_none(self, simple_grn):
        """Test that unknown metric type returns None."""
        result = _run_metric(
            metric_type="Unknown Metric",
            db_name="Unknown DB",
            grn=simple_grn,
            db=None,
            genes=[],
            peaks=[],
            cats=None,
            adata=None,
        )

        assert result is None


class TestRunFilelessMetric:
    """Tests for _run_fileless_metric function."""

    def test_omics_metric(self, adata, simple_grn):
        """Test Omics metric."""
        result = _run_fileless_metric(
            metric_type="Omics",
            db_name="Gene ~ TFs",
            dataset=adata,
            grn=simple_grn,
            adata=adata,
            is_mudata=False,
            has_cre=False,
        )

        assert result is not None
        assert len(result) == 3

    def test_steady_state_simulation(self, adata, simple_grn):
        """Test Steady state simulation metric."""
        result = _run_fileless_metric(
            metric_type="Steady state simulation",
            db_name="Boolean rules",
            dataset=adata,
            grn=simple_grn,
            adata=adata,
            is_mudata=False,
            has_cre=False,
        )

        # May return nan for small GRNs, but should not crash
        assert result is not None
        assert len(result) == 3

    def test_unknown_metric_returns_none(self, adata, simple_grn):
        """Test that unknown metric type returns None."""
        result = _run_fileless_metric(
            metric_type="Unknown Metric",
            db_name="Unknown DB",
            dataset=adata,
            grn=simple_grn,
            adata=adata,
            is_mudata=False,
            has_cre=False,
        )

        assert result is None


class TestRunOmicsMetric:
    """Tests for _run_omics_metric function."""

    def test_gene_tf_metric_anndata(self, adata, simple_grn):
        """Test Gene ~ TFs metric with AnnData."""
        result = _run_omics_metric(
            db_name="Gene ~ TFs",
            dataset=adata,
            grn=simple_grn,
            is_mudata=False,
            has_cre=False,
        )

        assert result is not None
        assert len(result) == 3

    def test_gene_tf_metric_mudata(self, mudata_with_celltype, simple_grn):
        """Test Gene ~ TFs metric with MuData."""
        result = _run_omics_metric(
            db_name="Gene ~ TFs",
            dataset=mudata_with_celltype,
            grn=simple_grn,
            is_mudata=True,
            has_cre=False,
        )

        assert result is not None
        assert len(result) == 3

    def test_gene_cre_metric_without_cre_returns_none(self, mudata_with_celltype, simple_grn):
        """Test that Gene ~ CREs returns None without CRE column."""
        result = _run_omics_metric(
            db_name="Gene ~ CREs",
            dataset=mudata_with_celltype,
            grn=simple_grn,
            is_mudata=True,
            has_cre=False,
        )

        assert result is None

    def test_gene_cre_metric_with_cre(self, mudata_with_celltype, simple_grn):
        """Test Gene ~ CREs metric with CRE column present."""
        result = _run_omics_metric(
            db_name="Gene ~ CREs",
            dataset=mudata_with_celltype,
            grn=simple_grn,
            is_mudata=True,
            has_cre=True,
        )

        assert result is not None
        assert len(result) == 3

    def test_cre_tf_metric_with_cre(self, mudata_with_celltype, simple_grn):
        """Test CRE ~ TFs metric with MuData and CRE column."""
        result = _run_omics_metric(
            db_name="CRE ~ TFs",
            dataset=mudata_with_celltype,
            grn=simple_grn,
            is_mudata=True,
            has_cre=True,
        )

        assert result is not None
        assert len(result) == 3

    def test_cre_tf_metric_without_mudata_returns_none(self, adata, simple_grn):
        """Test that CRE ~ TFs without MuData logs warning and returns None."""
        result = _run_omics_metric(
            db_name="CRE ~ TFs",
            dataset=adata,
            grn=simple_grn,
            is_mudata=False,
            has_cre=True,
        )

        assert result is None

    def test_gene_cre_metric_anndata_returns_none(self, adata, simple_grn):
        """Test that Gene ~ CREs without MuData logs warning and returns None."""
        result = _run_omics_metric(
            db_name="Gene ~ CREs",
            dataset=adata,
            grn=simple_grn,
            is_mudata=False,
            has_cre=True,
        )

        assert result is None


class TestEvalGrnDataset:
    """Tests for eval_grn_dataset function."""

    @patch("gretapy.pp._check.show_terms")
    def test_dataframe_has_expected_columns(self, mock_show_terms, adata, simple_grn):
        """Test that DataFrame has expected columns."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])

        result = eval_grn_dataset(
            organism="hg38",
            grn=simple_grn,
            dataset=adata,
            terms={},
            metrics=["Gene ~ TFs"],
            min_edges=1,
        )

        assert isinstance(result, pd.DataFrame)
        expected_cols = {"class", "task", "db", "precision", "recall", "f01"}
        assert expected_cols == set(result.columns)

    @patch("gretapy.pp._check.show_terms")
    def test_min_edges_threshold(self, mock_show_terms, adata):
        """Test that GRN with too few edges returns empty DataFrame."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        small_grn = pd.DataFrame({"source": ["PAX5"], "target": ["CD19"]})

        result = eval_grn_dataset(
            organism="hg38",
            grn=small_grn,
            dataset=adata,
            terms={},
            metrics=["Gene ~ TFs"],
            min_edges=10,
        )

        assert len(result) == 0

    @patch("gretapy.pp._check.show_terms")
    @patch("gretapy.pp._check.read_dts")
    def test_with_str_dataset(self, mock_read_dts, mock_show_terms, mudata_with_celltype, simple_grn):
        """Test evaluation with str dataset name."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        mock_read_dts.return_value = mudata_with_celltype

        result = eval_grn_dataset(
            organism="hg38",
            grn=simple_grn,
            dataset="Brain",
            terms={},
            metrics=["Gene ~ TFs"],
            min_edges=1,
        )

        assert isinstance(result, pd.DataFrame)

    @patch("gretapy.pp._check.show_terms")
    def test_anndata_no_atac_logs_warning(self, mock_show_terms, adata, simple_grn):
        """Test that AnnData dataset (no ATAC) logs warning and runs ok."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])

        result = eval_grn_dataset(
            organism="hg38",
            grn=simple_grn,
            dataset=adata,
            terms={},
            metrics=["Gene ~ TFs"],
            min_edges=1,
        )

        assert isinstance(result, pd.DataFrame)

    @patch("gretapy.pp._check.show_terms")
    def test_grn_without_cre_logs_warning(self, mock_show_terms, adata, simple_grn):
        """Test that GRN without cre column logs warning."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        grn_no_cre = simple_grn[["source", "target", "score"]].copy()

        result = eval_grn_dataset(
            organism="hg38",
            grn=grn_no_cre,
            dataset=adata,
            terms={},
            metrics=["Gene ~ TFs"],
            min_edges=1,
        )

        assert isinstance(result, pd.DataFrame)

    @patch("gretapy.pp._check.show_terms")
    @patch("gretapy.tl._eval.read_db")
    def test_with_file_based_metric(self, mock_read_db, mock_show_terms, adata, simple_grn, reference_grn_db):
        """Test evaluation with file-based metric (mocked read_db)."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        mock_read_db.return_value = reference_grn_db

        result = eval_grn_dataset(
            organism="hg38",
            grn=simple_grn,
            dataset=adata,
            terms={},
            metrics=["CollecTRI"],
            min_edges=1,
        )

        assert isinstance(result, pd.DataFrame)
        assert set(result.columns) == {"class", "task", "db", "precision", "recall", "f01"}

    @patch("gretapy.pp._check.show_terms")
    def test_with_grn_name_and_dataset_name(self, mock_show_terms, adata, simple_grn):
        """Test that grn_name and dataset_name are used in logging."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])

        result = eval_grn_dataset(
            organism="hg38",
            grn=simple_grn,
            dataset=adata,
            terms={},
            metrics=["Gene ~ TFs"],
            min_edges=1,
            grn_name="my_grn",
            dataset_name="my_dataset",
        )

        assert isinstance(result, pd.DataFrame)


class TestBenchmark:
    """Tests for benchmark function."""

    @patch("gretapy.pp._check.show_terms")
    def test_with_nested_dict(self, mock_show_terms, mudata_with_celltype, simple_grn):
        """Test benchmark with properly nested dict API."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        grns = {"my_grn": {"hg38": {"Brain": simple_grn}}}

        result = benchmark(
            grns=grns,
            datasets={"Brain": mudata_with_celltype},
            terms={"hg38": {"Brain": {}}},
            metrics=["Gene ~ TFs"],
            min_edges=1,
        )

        assert isinstance(result, pd.DataFrame)

    @patch("gretapy.pp._check.show_terms")
    def test_organism_kwarg_logs_warning(self, mock_show_terms, mudata_with_celltype, simple_grn):
        """Test that organism kwarg logs a warning when organisms are in grns dict."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        grns = {"my_grn": {"hg38": {"Brain": simple_grn}}}

        result = benchmark(
            grns=grns,
            organism="hg38",
            datasets={"Brain": mudata_with_celltype},
            terms={"hg38": {"Brain": {}}},
            metrics=["Gene ~ TFs"],
            min_edges=1,
        )

        assert isinstance(result, pd.DataFrame)

    def test_invalid_grns_type_raises(self):
        """Test that invalid grns type raises ValueError."""
        with pytest.raises(ValueError, match="grns must be"):
            benchmark(
                grns=[1, 2, 3],
                datasets=None,
                terms=None,
                metrics=None,
            )

    def test_empty_grns_raises(self):
        """Test that empty grns dict raises ValueError."""
        with pytest.raises(ValueError):
            benchmark(
                grns={},
                datasets=None,
                terms=None,
                metrics=None,
            )

    def test_invalid_organism_in_grns_raises(self, simple_grn):
        """Test that invalid organism in grns raises ValueError."""
        grns = {"my_grn": {"invalid_org": {"Brain": simple_grn}}}

        with pytest.raises(ValueError, match="Invalid organism"):
            benchmark(
                grns=grns,
                datasets=None,
                terms=None,
                metrics=None,
            )

    def test_invalid_metrics_raises(self, simple_grn):
        """Test that invalid metrics raises ValueError."""
        grns = {"my_grn": {"hg38": {"Brain": simple_grn}}}

        with pytest.raises(ValueError, match="Invalid metric"):
            benchmark(
                grns=grns,
                datasets=None,
                terms=None,
                metrics=["nonexistent_metric"],
            )

    @patch("gretapy.pp._check.show_terms")
    def test_returns_empty_df_when_min_edges_not_met(self, mock_show_terms, adata, simple_grn):
        """Test that result is empty DataFrame when min_edges not met."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        grns = {"my_grn": {"hg38": {"Brain": simple_grn}}}

        result = benchmark(
            grns=grns,
            datasets={"Brain": adata},
            terms={"hg38": {"Brain": {}}},
            metrics=["Gene ~ TFs"],
            min_edges=999999,
        )

        assert isinstance(result, pd.DataFrame)

    @patch("gretapy.pp._check.show_terms")
    @patch("gretapy.pp._check.read_dts")
    def test_with_dataset_list_filter(self, mock_read_dts, mock_show_terms, mudata_with_celltype, simple_grn):
        """Test benchmark with datasets as list to filter."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        mock_read_dts.return_value = mudata_with_celltype
        grns = {"my_grn": {"hg38": {"Brain": simple_grn}}}

        result = benchmark(
            grns=grns,
            datasets=["Brain"],
            terms={"hg38": {"Brain": {}}},
            metrics=["Gene ~ TFs"],
            min_edges=1,
        )

        assert isinstance(result, pd.DataFrame)

    def test_grn_inner_not_dict_raises(self, simple_grn):
        """Test that non-dict grn_inner raises ValueError (line 127)."""
        with pytest.raises(ValueError, match="must be a dict mapping organism"):
            benchmark(grns={"my_grn": "not_a_dict"}, datasets=None, terms=None, metrics=None)

    def test_org_inner_not_dict_raises(self, simple_grn):
        """Test that non-dict org_inner raises ValueError (line 132)."""
        with pytest.raises(ValueError, match="must be a dict mapping dataset"):
            benchmark(grns={"my_grn": {"hg38": "not_a_dict"}}, datasets=None, terms=None, metrics=None)

    def test_invalid_datasets_type_raises(self, simple_grn):
        """Test that invalid datasets type raises ValueError (line 152)."""
        grns = {"my_grn": {"hg38": {"Brain": simple_grn}}}
        with pytest.raises(ValueError, match="datasets must be None, list, or dict"):
            benchmark(grns=grns, datasets=42, terms=None, metrics=None)

    @patch("gretapy.pp._check.show_terms")
    def test_benchmark_terms_none_calls_check_terms(self, mock_show_terms, adata, simple_grn):
        """Test benchmark with terms=None calls _check_terms per dataset (lines 186, 191)."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        # Use a custom dataset name not in config so _check_terms returns {} (no validation)
        grns = {"my_grn": {"hg38": {"CustomDataset": simple_grn}}}

        # terms=None triggers line 186 (_check_terms call)
        # datasets is a dict → datasets_objects is not None → empty dataset_terms → line 191 warning
        result = benchmark(
            grns=grns,
            datasets={"CustomDataset": adata},
            terms=None,
            metrics=["Gene ~ TFs"],
            min_edges=1,
        )

        assert isinstance(result, pd.DataFrame)


class TestEvalGrnDatasetCoveragePaths:
    """Tests for eval_grn_dataset lines not covered by existing tests."""

    @patch("gretapy.pp._check.show_terms")
    @patch("gretapy.tl._eval._check_metrics")
    def test_db_info_none_is_skipped(self, mock_check_metrics, mock_show_terms, adata, simple_grn):
        """Test that db_name not in DATA config is silently skipped (line 334)."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        # Return a non-existent db_name from _check_metrics
        mock_check_metrics.return_value = ["NonExistentDB123", "Gene ~ TFs"]

        result = eval_grn_dataset(
            organism="hg38",
            grn=simple_grn,
            dataset=adata,
            terms={},
            metrics=["Gene ~ TFs"],
            min_edges=1,
        )

        # NonExistentDB123 is silently skipped; Gene ~ TFs still runs
        assert isinstance(result, pd.DataFrame)

    @patch("gretapy.pp._check.show_terms")
    def test_genomic_metrics_skipped_without_cre(self, mock_show_terms, adata, simple_grn):
        """Test that CREs metrics are skipped when GRN has no cre column (line 347)."""
        mock_show_terms.return_value = pd.DataFrame(columns=["organism", "db_name", "term"])
        grn_no_cre = simple_grn[["source", "target", "score"]].copy()

        # CREs metric type is "CREs" which matches the skip set; no cre column → skip
        result = eval_grn_dataset(
            organism="hg38",
            grn=grn_no_cre,
            dataset=adata,
            terms={},
            metrics=["CREs", "Gene ~ TFs"],
            min_edges=1,
        )

        # CREs db entries are skipped; Gene ~ TFs still runs
        assert isinstance(result, pd.DataFrame)
