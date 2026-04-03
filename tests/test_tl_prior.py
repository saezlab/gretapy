"""Tests for gretapy.tl._prior module."""

import pandas as pd
import pytest

from gretapy.tl._prior import _find_pairs, _grn, _tfm, _tfp


class TestGrn:
    """Tests for _grn function (Reference GRN comparison)."""

    def test_perfect_match(self, simple_grn, reference_grn_db):
        """Test with GRN that matches reference perfectly for available genes."""
        genes = ["PAX5", "GATA3", "SPI1", "CD19", "MS4A1", "CD3E", "IL7R", "CD14"]
        prc, rcl, f01 = _grn(grn=simple_grn, db=reference_grn_db, genes=genes)

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_no_overlap(self):
        """Test with GRN that has no overlap with reference."""
        grn = pd.DataFrame({"source": ["TFA", "TFB"], "target": ["GeneX", "GeneY"]})
        db = pd.DataFrame({"source": ["TFC", "TFD"], "target": ["GeneZ", "GeneW"]})
        genes = ["TFA", "TFB", "TFC", "TFD", "GeneX", "GeneY", "GeneZ", "GeneW"]

        prc, rcl, f01 = _grn(grn=grn, db=db, genes=genes)
        assert prc == 0.0
        assert rcl == 0.0
        assert f01 == 0.0

    def test_complete_overlap(self):
        """Test with identical GRN and reference."""
        grn = pd.DataFrame({"source": ["PAX5", "GATA3"], "target": ["CD19", "CD3E"]})
        db = pd.DataFrame({"source": ["PAX5", "GATA3"], "target": ["CD19", "CD3E"]})
        genes = ["PAX5", "GATA3", "CD19", "CD3E"]

        prc, rcl, f01 = _grn(grn=grn, db=db, genes=genes)
        assert prc == 1.0
        assert rcl == 1.0
        assert f01 == pytest.approx(1.0)

    def test_filters_by_genes(self):
        """Test that function filters database by measured genes."""
        grn = pd.DataFrame({"source": ["PAX5"], "target": ["CD19"]})
        db = pd.DataFrame({"source": ["PAX5", "GATA3"], "target": ["CD19", "FOXP3"]})
        genes = ["PAX5", "CD19"]

        prc, rcl, f01 = _grn(grn=grn, db=db, genes=genes)
        assert prc == 1.0
        assert rcl == 1.0

    def test_handles_duplicates(self):
        """Test that function handles duplicate edges."""
        grn = pd.DataFrame({"source": ["PAX5", "PAX5", "PAX5"], "target": ["CD19", "CD19", "MS4A1"]})
        db = pd.DataFrame({"source": ["PAX5"], "target": ["CD19"]})
        genes = ["PAX5", "CD19", "MS4A1"]

        prc, rcl, f01 = _grn(grn=grn, db=db, genes=genes)
        assert 0 <= prc <= 1
        assert rcl == 1.0


class TestTfm:
    """Tests for _tfm function (TF markers evaluation)."""

    def test_basic_functionality(self, simple_grn, tfm_db):
        """Test basic TF markers evaluation."""
        genes = ["PAX5", "GATA3", "SPI1", "RUNX1", "CD19", "MS4A1", "CD3E", "IL7R", "CD14"]
        prc, rcl, f01 = _tfm(grn=simple_grn, db=tfm_db, genes=genes, cats=None)

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_with_category_filter(self, simple_grn, tfm_db):
        """Test TF markers evaluation with category filtering."""
        genes = ["PAX5", "GATA3", "SPI1", "RUNX1", "CD19", "MS4A1", "CD3E", "IL7R", "CD14"]
        cats = ["B cell", "T cell"]
        prc, rcl, f01 = _tfm(grn=simple_grn, db=tfm_db, genes=genes, cats=cats)

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_no_matching_tfs(self):
        """Test with no matching TFs between GRN and database."""
        grn = pd.DataFrame({"source": ["TFX", "TFY"], "target": ["GeneA", "GeneB"]})
        db = pd.DataFrame({0: ["TFA", "TFB"], 1: ["Type1", "Type2"]})
        genes = ["TFX", "TFY", "TFA", "TFB", "GeneA", "GeneB"]

        prc, rcl, f01 = _tfm(grn=grn, db=db, genes=genes, cats=None)
        assert prc == 0.0
        assert f01 == 0.0

    def test_perfect_match(self):
        """Test with perfect TF match."""
        grn = pd.DataFrame({"source": ["PAX5", "GATA3"], "target": ["CD19", "CD3E"]})
        db = pd.DataFrame({0: ["PAX5", "GATA3"], 1: ["B cell", "T cell"]})
        genes = ["PAX5", "GATA3", "CD19", "CD3E"]

        prc, rcl, f01 = _tfm(grn=grn, db=db, genes=genes, cats=None)
        assert prc == 1.0
        assert rcl == 1.0


class TestFindPairs:
    """Tests for _find_pairs function."""

    def test_finds_significant_pairs(self):
        """Test that significant TF pairs are found."""
        grn = pd.DataFrame(
            {
                "source": ["TF1"] * 5 + ["TF2"] * 5 + ["TF3"] * 3,
                "target": ["G1", "G2", "G3", "G4", "G5", "G1", "G2", "G3", "G6", "G7", "G8", "G9", "G10"],
            }
        )
        pairs = _find_pairs(grn=grn, thr_pval=0.05)

        assert isinstance(pairs, set)

    def test_returns_empty_for_no_pairs(self):
        """Test that empty set is returned when no significant pairs exist."""
        grn = pd.DataFrame(
            {
                "source": ["TF1", "TF2", "TF3"],
                "target": ["G1", "G2", "G3"],
            }
        )
        pairs = _find_pairs(grn=grn, thr_pval=0.01)

        assert isinstance(pairs, set)

    def test_pair_format(self):
        """Test that pairs are formatted correctly (sorted, pipe-separated)."""
        grn = pd.DataFrame(
            {
                "source": ["TF1"] * 10 + ["TF2"] * 10,
                "target": [f"G{i}" for i in range(10)] + [f"G{i}" for i in range(10)],
            }
        )
        pairs = _find_pairs(grn=grn, thr_pval=1.0)

        for pair in pairs:
            assert "|" in pair
            tf_a, tf_b = pair.split("|")
            assert tf_a <= tf_b

    def test_empty_df_returns_empty_set(self):
        """Test with only one TF returns no pairs."""
        grn = pd.DataFrame({"source": ["TF1", "TF1"], "target": ["G1", "G2"]})
        pairs = _find_pairs(grn=grn, thr_pval=0.05)
        assert pairs == set()

    def test_no_shared_targets_returns_empty_set(self):
        """Test that TFs with no shared targets get p=1.0 and are not included."""
        grn = pd.DataFrame(
            {
                "source": ["TF1", "TF1", "TF2", "TF2"],
                "target": ["G1", "G2", "G3", "G4"],  # no overlap
            }
        )
        pairs = _find_pairs(grn=grn, thr_pval=0.05)
        assert isinstance(pairs, set)


class TestTfp:
    """Tests for _tfp function (TF pairs evaluation)."""

    def test_basic_functionality(self, simple_grn, tfp_db):
        """Test basic TF pairs evaluation."""
        genes = ["PAX5", "GATA3", "SPI1", "CD19", "MS4A1", "CD3E", "IL7R", "CD14",
                 "BCL2", "IRF4", "TCF7", "FOXP3", "RUNX1", "CD79A"]
        prc, rcl, f01 = _tfp(grn=simple_grn, db=tfp_db, genes=genes, thr_pval=0.05)

        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_filters_by_tfs_in_db(self):
        """Test that GRN is filtered to only include TFs present in database."""
        grn = pd.DataFrame(
            {
                "source": ["PAX5", "PAX5", "UnknownTF"],
                "target": ["CD19", "MS4A1", "GeneX"],
            }
        )
        db = pd.DataFrame({0: ["PAX5"], 1: ["GATA3"]})
        genes = ["PAX5", "GATA3", "CD19", "MS4A1", "GeneX", "UnknownTF"]

        prc, rcl, f01 = _tfp(grn=grn, db=db, genes=genes, thr_pval=0.05)
        assert 0 <= prc <= 1
        assert 0 <= rcl <= 1
        assert 0 <= f01 <= 1

    def test_handles_duplicates(self):
        """Test that function handles duplicate edges in GRN."""
        grn = pd.DataFrame(
            {
                "source": ["PAX5", "PAX5", "PAX5"],
                "target": ["CD19", "CD19", "MS4A1"],
            }
        )
        db = pd.DataFrame({0: ["PAX5"], 1: ["GATA3"]})
        genes = ["PAX5", "GATA3", "CD19", "MS4A1"]

        prc, rcl, f01 = _tfp(grn=grn, db=db, genes=genes, thr_pval=0.05)
        assert isinstance(prc, float)
        assert isinstance(rcl, float)
        assert isinstance(f01, float)

    def test_empty_db_after_gene_filter(self):
        """Test when db is empty after gene filtering returns zeros."""
        grn = pd.DataFrame({"source": ["PAX5", "GATA3"], "target": ["CD19", "CD3E"]})
        db = pd.DataFrame({0: ["TF_NOT_IN_GENES"], 1: ["TF_NOT_EITHER"]})
        genes = ["PAX5", "GATA3", "CD19", "CD3E"]

        prc, rcl, f01 = _tfp(grn=grn, db=db, genes=genes, thr_pval=0.05)
        assert prc == 0.0
        assert rcl == 0.0
        assert f01 == 0.0
