# Myeloid differentiation trajectories shape response to combined immuno-chemotherapy across multiple solid tumors

Analysis code for the study **"Myeloid differentiation trajectories shape
response to combined immuno-chemotherapy across multiple solid tumors."**

This repository contains the R / shell scripts used to integrate, annotate and
interrogate the **tumor-infiltrating myeloid (TIM)** compartment across four
solid malignancies — **NSCLC, TNBC, PC and ESCC** — together with the
downstream analyses (myeloid subclustering, cell-cell communication,
pseudotime trajectory, co-expression networks and statistical modelling of
clinical response).

**Cohort.** 71 tumor samples from 53 patients across four cohorts:
NSCLC (`GSE207422`), TNBC (`GSE169246`), PC (`GSE267814`) and ESCC
(`OMIX005710`). After quality control, **401,550 cells** were retained;
**74,085 myeloid and mast cells** were resolved into **12 functional subtypes**.

**Key findings.** Two conserved, clinically opposing macrophage subsets —
**Macro-APOE** (cholesterol-metabolic, matrix-remodelling, associated with
non-response) and **Macro-MS4A6A** (antigen-presenting, trafficking,
associated with favourable response) — arise from a shared monocyte-derived
differentiation trajectory. Co-expression network rewiring (hdWGCNA) further
links successful response to a retreat from hyper-metabolic programs toward
homeostatic, tissue-stabilising networks.

> **Note** — only **code and documentation** are tracked in this repository.
> Raw sequencing data, intermediate `Seurat`/`RDS` objects, figures and tables
> are intentionally **excluded** (see [Data availability](#data-availability)).
> Experimental / deprecated scripts are preserved under [`archive/`](archive/).

---

## Table of contents

- [Repository layout](#repository-layout)
- [Analysis workflow](#analysis-workflow)
- [Module descriptions](#module-descriptions)
- [Requirements](#requirements)
- [Data availability](#data-availability)
- [How to run](#how-to-run)
- [Archive](#archive)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Repository layout

```
.
├── README.md                # this file
├── LICENSE                  # MIT
├── CITATION.cff             # machine-readable citation
├── .gitignore               # whitelist: code + docs only
├── archive/                 # deprecated / test scripts (kept for provenance)
├── data/                    # data-import helper scripts (per cohort)
│   └── ESCC/  PC/  AML/  NSCLC/  TNBC1/  TNBC2/
└── WorkDir/                 # main analysis, one folder per step
    ├── preparing/           #  1. load / QC / merge per cohort
    ├── DoubletFinder/       #  2. doublet detection & removal
    ├── SCT_Integrate/       #  3. SCTransform + RPCA integration
    ├── Log_Integrate/       #  3b. LogNormalize-based integration
    ├── Annotation/          #  4. lineage annotation & validation
    ├── Clustering_subtype/  #  5. myeloid / T-NK subclustering
    ├── Description/         #  6. log-linear models & correspondence analysis
    ├── GeneSetEnrichment/   #  7. GSEA / over-representation (KEGG)
    ├── Cellchat/            #  8. cell-cell communication
    ├── Monocle/             #  9. pseudotime / trajectory (monocle3)
    ├── hdWGCNA/             # 10. co-expression networks
    ├── paper_writing/       # 11. main-figure generation
    ├── NMF/                 #  (additional) NMF programs
    ├── Prediction/          #  (additional) signature-based prediction
    └── Shiny/               #  (additional) interactive Shiny application
```

The four cohorts analysed in the manuscript are **NSCLC, TNBC, PC and ESCC**
(the manuscript's TNBC cohort corresponds to `TNBC1` / `GSE169246`).

The repository additionally retains several **exploratory datasets and modules
that are _not_ part of the final manuscript**, kept only for completeness and
provenance:

- the **AML** cohort (`data/AML/`, `WorkDir/**/AML/`, `GSE198052`);
- a second TNBC dataset, **TNBC2** (`GSE266919`), used for cross-checking;
- the melanoma dataset **GSE123813**, used in
  `WorkDir/Clustering_subtype/code/Myeloid/Map_Annotation.R`;
- the `NMF/`, `Prediction/` and `Shiny/` analysis modules.

Only the four manuscript cohorts and the modules from `preparing` through
`paper_writing` are described in the paper.

Most step folders follow the same internal convention:

```
<module>/
├── code/     # analysis scripts (+ sbatch*.sh for the SLURM cluster)
├── plot/     # figures   (NOT tracked)
└── data/     # objects   (NOT tracked; *.txt here are str() dumps for reference)
```

---

## Analysis workflow

The scripts are organised as a linear pipeline. Numbers below correspond to
the folder ordering and to the computational dependency between steps.

```
Raw data (GEO / NGDC)
        |
        v
[1] preparing            load, QC, per-cohort merge
        |
        v
[2] DoubletFinder        detect & remove doublets
        |
        v
[3] SCT_Integrate        SCTransform + RPCA integration (IntegrateLayers)
[3b] Log_Integrate       LogNormalize integration  (alternative strategy)
        |
        v
[4] Annotation           major cell-lineage annotation & validation
        |
        v
[5] Clustering_subtype   myeloid subclustering (T/NK via starCAT in parallel)
        |
        +-- [6] Description          log-linear models / correspondence analysis
        +-- [7] GeneSetEnrichment    GSEA / ORA (KEGG, clusterProfiler)
        +-- [8] Cellchat             ligand-receptor communication
        +-- [9] Monocle              monocle3 pseudotime trajectory
        +-- [10] hdWGCNA             co-expression modules
                    |
                    v
              [11] paper_writing   main figures
```

Scripts were executed on a SLURM cluster; the accompanying `sbatch*.sh`
files document resource requests and the run order.

---

## Module descriptions

| # | Folder | What it does |
|---|--------|--------------|
| 1 | `preparing/` | Per-cohort loading, quality control, metadata and sample renaming for ESCC, PC, NSCLC, TNBC. |
| 2 | `DoubletFinder/` | Doublet detection per sample (pN = 0.25, expected doublet rate = 0.06) and removal after integration. |
| 3 | `SCT_Integrate/` | SCTransform normalisation, RPCA (IntegrateLayers) integration and integration-quality evaluation. |
| 3b | `Log_Integrate/` | LogNormalize-based integration used as a parallel strategy. |
| 4 | `Annotation/` | Major cell-lineage annotation, marker inspection, and validation against original study annotations. |
| 5 | `Clustering_subtype/` | Myeloid subset re-integration, clustering sweeps, marker discovery and subtype annotation; parallel T/NK branch (`T_NK/`, starCAT annotation against the TCAT.V1 reference). |
| 6 | `Description/` | Linear mixed-effects models (lme4), log-linear interaction modelling (MASS) and correspondence analysis (`ca`) of subtypes across cancers. |
| 7 | `GeneSetEnrichment/` | GSEA / ORA on the KEGG database via `clusterProfiler`. |
| 8 | `Cellchat/` | CellChat cell-cell communication between myeloid subsets and T-cell states. |
| 9 | `Monocle/` | monocle3 pseudotime / trajectory, branch analysis and branch-specific GSEA. |
| 10 | `hdWGCNA/` | hdWGCNA co-expression networks for ESCC and PC (metacells, k = 25), module gene-set curation and DME analysis. |
| 11 | `paper_writing/` | Scripts reproducing the main result figures (`Result1`-`Result4`). |
| — | `NMF/`, `Prediction/`, `Shiny/` | Additional, exploratory analyses not included in the final manuscript. |

**TIM subtypes resolved (12).** Three dendritic-cell subsets (DC-HLA,
DC-LAMP3, DC-CPVL), one mast-cell population, two monocyte lineages
(Mono-FCN1, Mono-TIMP1), five macrophage subtypes (Macro-MS4A6A, Macro-APOE,
Macro-FOSB, Macro-CCL, Macro-MARCO) and one neutrophil population
(Neutro-FCGR3B).

---

## Requirements

Analyses were written in **R** (Seurat **v5.2.1**). Core packages:

**Single-cell** — `Seurat` (v5), `SeuratObject`, `sctransform`, `DoubletFinder`,
`SingleCellExperiment`, `SeuratDisk`, `hdf5r`, `scSHC`.

**Trajectory / networks** — `monocle3`, `hdWGCNA`, `WGCNA`, `igraph`, `ggraph`, `tidygraph`.

**Communication / enrichment** — `CellChat`, `clusterProfiler`, `UCell`; T-cell
annotation via **starCAT** / TCAT.V1 (Python, run through `reticulate`).

**Statistics** — `lme4`, `lmerTest`, `MASS` (log-linear models), `ca`
(correspondence analysis), `NMF`, `nnet`.

**Visualisation** — `ggplot2`, `cowplot`, `patchwork`, `ggrepel`, `ggpubr`,
`RColorBrewer`, `clustree`, `pheatmap`, `gridExtra`, `pdftools`.

Shell scripts assume a **SLURM** scheduler. A minimal installation helper is
provided at `WorkDir/hdWGCNA/code/install.R`.

---

## Data availability

This repository does **not** ship raw data or large intermediate objects.

The study used four publicly available single-cell datasets:

| Cancer type | Accession | Repository |
|-------------|-----------|------------|
| NSCLC | `GSE207422` | GEO |
| TNBC  | `GSE169246` | GEO |
| PC    | `GSE267814` | GEO |
| ESCC  | `OMIX005710` | NGDC |

The `*.txt` files under `WorkDir/*/data/**/str(...).txt` are lightweight
`str()` dumps of the `Seurat` objects; they are kept only to document object
structure (cell numbers, assays, metadata columns) and contain no expression
data.

---

## How to run

1. Download the four datasets above and place them under `data/<COHORT>/`
   following the paths hard-coded (or configured) in `WorkDir/preparing/code/`.
2. Run the modules **in order** (1 -> 11). Each module folder contains its own
   `sbatch*.sh`; submit with `sbatch <script>.sh`, or source the `.R` scripts
   interactively for smaller steps.
3. Paths are currently absolute/cluster-specific. Adjust the input/output
   paths at the top of each script to match your environment before running.

> The scripts are shared **for transparency and reproducibility**, not as a
> turn-key pipeline: they are the exact cluster scripts used for the study.

---

## Archive

`archive/` keeps scripts that were used during development but are **not part
of the final pipeline** (unit tests, ad-hoc experiments, failed builds, trial
runs). They are retained so that the analysis history remains traceable.
Examples: `WorkDir/Description/code/test.R`, `WorkDir/Shiny/build_shiny_failed.txt`,
`WorkDir/*/test*.R`.

---

## Citation

If you use this code, please cite the accompanying manuscript:

```
Wang G, Huang Y, Ding Y, Gao X, Zhai S, Zhang F.
Myeloid differentiation trajectories shape response to combined
immuno-chemotherapy across multiple solid tumors. (2025).
```

A machine-readable `CITATION.cff` is included and will be updated with the
journal reference and DOI once the paper is published.

---

## License

Released under the **MIT License** — see [`LICENSE`](LICENSE).

---

## Contact

**Guosheng Wang** — [@GoulsonWang](https://github.com/GoulsonWang)

Issues and pull requests are welcome.