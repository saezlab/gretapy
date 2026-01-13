import os

PATH_DATA = os.path.join(".", "pygreta_data")
ID_ZENODO = 17872739
URL_STR = f"https://zenodo.org/records/{ID_ZENODO}/files/"
URL_END = "?download=1"

DATA = {
    "hg38": {
        "terms": "hg38_terms.csv.gz",
        "dbs": {
            # Prior knowledge
            "Human Protein Atlas (HPA)": {
                "fname": "hg38_tfm_hpa.tsv.gz",
                "metric": "TF markers",
            },
            "TF-Marker": {
                "fname": "hg38_tfm_tfmdb.tsv.gz",
                "metric": "TF markers",
            },
            "Europe PMC": {
                "fname": "hg38_tfp_europmc.tsv.gz",
                "metric": "TF pairs",
            },
            "IntAct": {
                "fname": "hg38_tfp_intact.tsv.gz",
                "metric": "TF pairs",
            },
            "CollecTRI": {
                "fname": "hg38_gst_collectri.csv.gz",
                "metric": "Reference GRN",
            },
            # Genomic
            "ChIP-Atlas": {
                "fname": "hg38_tfb_chipatlas.bed.gz",
                "metric": "TF binding",
            },
            "ReMap 2022": {
                "fname": "hg38_tfb_remap2022.bed.gz",
                "metric": "TF binding",
            },
            "UniBind": {
                "fname": "hg38_tfb_unibind.bed.gz",
                "metric": "TF binding",
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
                "metric": "CRE to gene links",
            },
            # Predictive
            "Hallmarks": {
                "fname": "hg38_gst_hall.csv.gz",
                "metric": "Gene sets",
            },
            "KEGG": {
                "fname": "hg38_gst_kegg.csv.gz",
                "metric": "Gene sets",
            },
            "Reactome": {
                "fname": "hg38_gst_reac.csv.gz",
                "metric": "Gene sets",
            },
            "PROGENy": {
                "fname": "hg38_gst_prog.csv.gz",
                "metric": "Gene sets",
            },
            "gene ~ TFs": {
                "fname": None,
                "metric": "Omics",
            },
            "gene ~ CREs": {
                "fname": None,
                "metric": "Omics",
            },
            "CRE ~ TFs": {
                "fname": None,
                "metric": "Omics",
            },
            # Mechanistic
            "KnockTF (scoring)": {
                "fname": "hg38_prt_knocktf.h5ad",
                "metric": "TF scoring",
            },
            "KnockTF (forecasting)": {
                "fname": "hg38_prt_knocktf.h5ad",
                "metric": "Perturbation forecasting",
            },
            "Boolean rules": {
                "fname": None,
                "metric": "Steady state simulation",
            },
        },
        "dts": {
            "brain": {
                "fname": "hg38_dts_brain.h5mu",
                "pubmed": "",
                "geo": "",
            },
            "breast": {"fname": "hg38_dts_breast.h5mu", "pubmed": "", "geo": ""},
            "embryo": {"fname": "hg38_dts_embryo.h5mu", "pubmed": "", "geo": ""},
            "eye": {"fname": "hg38_dts_eye.h5mu", "pubmed": "", "geo": ""},
            "kidney": {"fname": "hg38_dts_kidney.h5mu", "pubmed": "", "geo": ""},
            "lung": {"fname": "hg38_dts_lung.h5mu", "pubmed": "", "geo": ""},
            "heart": {"fname": "hg38_dts_heart.h5mu", "pubmed": "", "geo": ""},
            "pbmc10k": {
                "fname": "hg38_dts_pbmc10k.h5mu",
                "pubmed": "",
                "geo": "",
                "terms": {
                    "TF-Marker": [
                        "B cell",
                        "Dendritic cell",
                        "Lymphoid cell",
                        "Macrophage cell",
                        "Macrophages cell",
                        "Macrophagocyte cell",
                        "Mononuclear cell",
                        "Natural killer cell",
                        "Peripheral blood cell",
                        "T cell",
                    ],
                    "Human Protein Atlas (HPA)": [
                        "B",
                        "Dendritic cells",
                        "Erythroid cells",
                        "Lymphoid tissue",
                        "Macrophages",
                        "Monocytes",
                        "NK",
                        "T",
                    ],
                    "ChIP-Atlas": [
                        "B-Cell",
                        "B-cell (CD19+)",
                        "Bone Marrow Cells",
                        "Bone marrow mononuclear cells",
                        "Bone marrow nuclear cells",
                        "CD20+ B cells",
                        "CD4+ T cells",
                        "CD8+ T cells",
                        "Dendritic Cells",
                        "Macrophages",
                        "Mast Cells",
                        "Memory B cells",
                        "Memory T cells",
                        "Monocytes",
                        "Monocytes-CD14+",
                        "Naive B cells",
                        "Natural Killer T-Cells",
                        "Natural killer cells",
                        "T cells",
                        "Th1 Cells",
                        "Th17 Cells",
                        "Th2 Cells",
                        "Treg",
                    ],
                    "ReMap 2022": [
                        "B-cell",
                        "CD4",
                        "CD4-pos",
                        "CD8",
                        "T-cell",
                        "Th1",
                        "Th17",
                        "macrophage",
                        "monocyte",
                        "peripheral-blood-mononuclear-cell",
                        "primary-B-cell",
                        "primary-monocyte",
                    ],
                    "UniBind": [
                        "B cells",
                        "CD4 CD25 CD45RA T cells",
                        "CD4 CD25 T cells",
                        "CD4 T cells",
                        "T cells",
                        "Th17 cells",
                        "Th2 cells",
                        "macrophage",
                        "macrophage hypo il",
                        "macrophage il4",
                        "macrophages",
                        "macrophages ifng",
                        "macrophages il4",
                        "macrophages tpp",
                        "monocyte derived macrophages",
                        "monocyte macro",
                        "monocytes blood",
                        "monocytes of peripheral blood",
                        "peripheral blood mononuclear cells",
                    ],
                    "GWAS Catalog": [
                        "chronic lymphocytic leukemia",
                        "Hodgkins lymphoma, multiple myeloma, chronic lymphocytic leukemia",
                        "acute myeloid leukemia",
                        "acute lymphoblastic leukemia",
                        "B-cell lymphoma/leukemia 10 measurement",
                        "induced myeloid leukemia cell differentiation protein Mcl-1 measurement",
                        "leukemia inhibitory factor receptor measurement",
                        "leukemia inhibitory factor measurement",
                        "lymphoid leukemia",
                        "B-cell acute lymphoblastic leukemia",
                        "chronic myelogenous leukemia",
                        "childhood acute lymphoblastic leukemia",
                    ],
                    "eQTL Catalogue": [
                        "B cell",
                        "CD16+ monocyte",
                        "CD4+ T cell",
                        "CD8+ T cell",
                        "NK cell",
                        "Tfh cell",
                        "Th1 cell",
                        "Th17 cell",
                        "Th2 cell",
                        "Treg memory",
                        "Treg naive",
                        "blood",
                        "macrophage",
                        "monocyte",
                        "neutrophil",
                    ],
                    "KnockTF (scoring)": [
                        "Blood",
                        "Haematopoietic_and_lymphoid_tissue",
                        "Haematopoietic_and_lymphoid_tissue_Blood",
                    ],
                    "KnockTF (forecasting)": [
                        "Blood",
                        "Haematopoietic_and_lymphoid_tissue",
                        "Haematopoietic_and_lymphoid_tissue_Blood",
                    ],
                },
            },
            "pitupair": {"fname": "hg38_dts_pitupair.h5mu", "pubmed": "", "geo": ""},
            "reprofibro": {"fname": "hg38_dts_reprofibro.h5mu", "pubmed": "", "geo": ""},
            "skin": {"fname": "hg38_dts_skin.h5mu", "pubmed": "", "geo": ""},
        },
    },
    "mm10": {
        "terms": "mm10_terms.csv.gz",
        "dbs": {
            # Prior knowledge
            "CollecTRI": {
                "fname": "mm10_gst_collectri.csv.gz",
                "metric": "Reference GRN",
            },
            # Genomic
            "ChIP-Atlas": {
                "fname": "mm10_tfb_chipatlas.bed.gz",
                "metric": "TF binding",
            },
            "ReMap 2022": {
                "fname": "mm10_tfb_remap2022.bed.gz",
                "metric": "TF binding",
            },
            "UniBind": {
                "fname": "mm10_tfb_unibind.bed.gz",
                "metric": "TF binding",
            },
            "ENCODE Blacklist": {
                "fname": "hg38_cre_blacklist.bed.gz",
                "metric": "CREs",
            },
            "ENCODE CREs": {
                "fname": "hg38_cre_encode.bed.gz",
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
            # Predictive
            "Hallmarks": {
                "fname": "mm10_gst_hall.csv.gz",
                "metric": "Gene sets",
            },
            "Reactome": {
                "fname": "mm10_gst_reac.csv.gz",
                "metric": "Gene sets",
            },
            "PROGENy": {
                "fname": "mm10_gst_prog.csv.gz",
                "metric": "Gene sets",
            },
            "gene ~ TFs": {
                "fname": None,
                "metric": "Omics",
            },
            "gene ~ CREs": {
                "fname": None,
                "metric": "Omics",
            },
            "CRE ~ TFs": {
                "fname": None,
                "metric": "Omics",
            },
            # Mechanistic
            "KnockTF (scoring)": {
                "fname": "m10_prt_knocktf.h5ad.gz",
                "metric": "TF scoring",
            },
            "KnockTF (forecasting)": {
                "fname": "m10_prt_knocktf.h5ad.gz",
                "metric": "Perturbation forecasting",
            },
            "Boolean rules": {
                "fname": None,
                "metric": "Steady state simulation",
            },
        },
        "dts": {
            "pbmc10k": {"fname": "mm10_dts_pbmc.h5mu", "pubmed": "", "geo": ""},
        },
    },
}

METRIC_CATS = {
    "TF markers": "Prior Knowledge",
    "TF pairs": "Prior Knowledge",
    "Reference GRN": "Prior Knowledge",
    "TF binding": "Genomic",
    "CREs": "Genomic",
    "CRE to gene links": "Genomic",
    "Gene sets": "Predictive",
    "Omics": "Predictive",
    "TF scoring": "Mechanistic",
    "Perturbation forecasting": "Mechanistic",
    "Steady state simulation": "Mechanistic",
}
