import os
import re

import yaml
from importlib.resources import files

PATH_DATA = os.path.join(".", "gretapy_data")
ID_ZENODO = 19600656
URL_STR = f"https://zenodo.org/records/{ID_ZENODO}/files/"
URL_END = "?download=1"


def _slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


DATA = {
    "hg38": {
        "gann": "hg38_gann.bed.gz",
        "terms": "hg38_terms.csv.gz",
        "dbs": {
            # Extra
            "Lambert TFs": {
                "fname": "hg38_tf_lambert.txt.gz",
                "metric": None,
            },
            "DoRothEA": {
                "fname": "hg38_gst_dorothea.csv.gz",
                "metric": None,
            },
            # Literature
            "HPA": {
                "fname": "hg38_tfm_hpa.tsv.gz",
                "metric": "TF Markers",
            },
            "TF-Marker": {
                "fname": "hg38_tfm_tfmdb.tsv.gz",
                "metric": "TF Markers",
            },
            "Europe PMC": {
                "fname": "hg38_tfp_europmc.tsv.gz",
                "metric": "TF Pairs",
            },
            "IntAct": {
                "fname": "hg38_tfp_intact.tsv.gz",
                "metric": "TF Pairs",
            },
            "CollecTRI": {
                "fname": "hg38_gst_collectri.csv.gz",
                "metric": "Reference GRN",
            },
            # Genomic
            "ChIP-Atlas": {
                "fname": "hg38_tfb_chipatlas.bed.gz",
                "metric": "TF Binding",
            },
            "ReMap 2022": {
                "fname": "hg38_tfb_remap2022.bed.gz",
                "metric": "TF Binding",
            },
            "UniBind": {
                "fname": "hg38_tfb_unibind.bed.gz",
                "metric": "TF Binding",
            },
            "ENCODE Blacklist": {
                "fname": "hg38_cre_blacklist.bed.gz",
                "metric": "CREs",
            },
            "ENCODE CREs": {
                "fname": "hg38_cre_encode.bed.gz",
                "metric": "CREs",
            },
            "GWAS Catalog": {
                "fname": "hg38_cre_gwascatalogue.bed.gz",
                "metric": "CREs",
            },
            "phastCons": {
                "fname": "hg38_cre_phastcons.bed.gz",
                "metric": "CREs",
            },
            "Promoters": {
                "fname": "hg38_cre_promoters.bed.gz",
                "metric": "CREs",
            },
            "Zhang21": {
                "fname": "hg38_cre_zhang21.bed.gz",
                "metric": "CREs",
            },
            "eQTL Catalogue": {
                "fname": "hg38_c2g_eqtlcatalogue.bed.gz",
                "metric": "CRE to Gene",
            },
            # Predictive
            "Hallmarks": {
                "fname": "hg38_gst_hall.csv.gz",
                "metric": "Gene Sets",
            },
            "KEGG": {
                "fname": "hg38_gst_kegg.csv.gz",
                "metric": "Gene Sets",
            },
            "Reactome": {
                "fname": "hg38_gst_reac.csv.gz",
                "metric": "Gene Sets",
            },
            "PROGENy": {
                "fname": "hg38_gst_prog.csv.gz",
                "metric": "Gene Sets",
            },
            "Gene ~ TFs": {
                "fname": None,
                "metric": "Omics",
            },
            "Gene ~ CREs": {
                "fname": None,
                "metric": "Omics",
            },
            "CRE ~ TFs": {
                "fname": None,
                "metric": "Omics",
            },
            # Mechanistic
            "KnockTF (scoring)": {
                "fname": "hg38_prt_knocktf.h5ad.gz",
                "metric": "TF Scoring",
            },
            "KnockTF (forecasting)": {
                "fname": "hg38_prt_knocktf.h5ad.gz",
                "metric": "Perturbation Forecasting",
            },
            "Boolean rules": {
                "fname": None,
                "metric": "Steady State Simulation",
            },
        },
        "dts": {
            "Brain": {
                "fname": "hg38_dts_brain_all.h5mu",
                "pubmed": "38491288",
                "geo": "GSE193688",
            },
            "Breast": {
                "fname": "hg38_dts_breast_all.h5mu",
                "pubmed": "39122969",
                "geo": None,
            },
            "Embryo": {
                "fname": "hg38_dts_embryo_all.h5mu",
                "pubmed": "37369347",
                "geo": "GSE218314",
            },
            "Eye": {
                "fname": "hg38_dts_eye_all.h5mu",
                "pubmed": "36277849",
                "geo": "GSE196235",
            },
            "Kidney": {
                "fname": "hg38_dts_kidney_all.h5mu",
                "pubmed": "37333123",
                "geo": None,
            },
            "Lung": {
                "fname": "hg38_dts_lung_all.h5mu",
                "pubmed": "39266564",
                "geo": "GSE241468",
            },
            "Heart": {
                "fname": "hg38_dts_heart_all.h5mu",
                "pubmed": "37438528",
                "geo": None,
            },
            "PBMC": {
                "fname": "hg38_dts_pbmc10k_all.h5mu",
                "pubmed": None,
                "geo": None,
            },
            "Pituitary": {
                "fname": "hg38_dts_pitupair_all.h5mu",
                "pubmed": "35263594",
                "geo": "GSE178454",
            },
            "Repr. Fibroblasts": {
                "fname": "hg38_dts_reprofibro_all.h5mu",
                "pubmed": "37873116",
                "geo": "GSE242419",
            },
            "Skin": {
                "fname": "hg38_dts_skin_all.h5mu",
                "pubmed": "36121124",
                "geo": "GSE207335",
            },
        },
    },
    "mm10": {
        "terms": "mm10_terms.csv.gz",
        "dbs": {
            # Extra
            "Lambert TFs": {
                "fname": "mm10_tf_lambert.txt.gz",
                "metric": None,
            },
            "DoRothEA": {
                "fname": "mm10_gst_dorothea.csv.gz",
                "metric": None,
            },
            # Literature
            "CollecTRI": {
                "fname": "mm10_gst_collectri.csv.gz",
                "metric": "Reference GRN",
            },
            # Genomic
            "ChIP-Atlas": {
                "fname": "mm10_tfb_chipatlas.bed.gz",
                "metric": "TF Binding",
            },
            "ReMap 2022": {
                "fname": "mm10_tfb_remap2022.bed.gz",
                "metric": "TF Binding",
            },
            "UniBind": {
                "fname": "mm10_tfb_unibind.bed.gz",
                "metric": "TF Binding",
            },
            "ENCODE Blacklist": {
                "fname": "mm10_cre_blacklist.bed.gz",
                "metric": "CREs",
            },
            "ENCODE CREs": {
                "fname": "mm10_cre_encode.bed.gz",
                "metric": "CREs",
            },
            "phastCons": {
                "fname": "mm10_cre_phastcons.bed.gz",
                "metric": "CREs",
            },
            "Promoters": {
                "fname": "mm10_cre_promoters.bed.gz",
                "metric": "CREs",
            },
            # Predictive
            "Hallmarks": {
                "fname": "mm10_gst_hall.csv.gz",
                "metric": "Gene Sets",
            },
            "Reactome": {
                "fname": "mm10_gst_reac.csv.gz",
                "metric": "Gene Sets",
            },
            "PROGENy": {
                "fname": "mm10_gst_prog.csv.gz",
                "metric": "Gene Sets",
            },
            "Gene ~ TFs": {
                "fname": None,
                "metric": "Omics",
            },
            "Gene ~ CREs": {
                "fname": None,
                "metric": "Omics",
            },
            "CRE ~ TFs": {
                "fname": None,
                "metric": "Omics",
            },
            # Mechanistic
            "KnockTF (scoring)": {
                "fname": "mm10_prt_knocktf.h5ad.gz",
                "metric": "TF Scoring",
            },
            "KnockTF (forecasting)": {
                "fname": "mm10_prt_knocktf.h5ad.gz",
                "metric": "Perturbation Forecasting",
            },
            "Boolean rules": {
                "fname": None,
                "metric": "Steady State Simulation",
            },
        },
        "dts": {
            "Palate": {
                "fname": "mm10_dts_epalate_all.h5mu",
                "pubmed": "38280850",
                "geo": "GSE218576",
            },
        },
    },
}

METRIC_CATS = {
    "TF Markers": "Literature",
    "TF Pairs": "Literature",
    "Reference GRN": "Literature",
    "TF Binding": "Genomic",
    "CREs": "Genomic",
    "CRE to Gene": "Genomic",
    "Gene Sets": "Predictive",
    "Omics": "Predictive",
    "TF Scoring": "Mechanistic",
    "Perturbation Forecasting": "Mechanistic",
    "Steady State Simulation": "Mechanistic",
}


def _load_all_terms() -> None:
    data_pkg = files("gretapy.data")
    for org, org_data in DATA.items():
        for dts_name, dts_data in org_data.get("dts", {}).items():
            fname = f"{org}_{_slug(dts_name)}_terms.yaml"
            with data_pkg.joinpath(fname).open("r") as f:
                dts_data["terms"] = yaml.safe_load(f)


_load_all_terms()
