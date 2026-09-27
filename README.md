# Myeloid Cell Atlas across Human Cancers

Analysis code for the pan-cancer **myeloid / macrophage cell atlas** study.
This repository contains the R / shell scripts used to build, annotate, and
interrogate the myeloid compartment across six single-cell RNA-seq cohorts
(**ESCC, PC, AML, NSCLC, TNBC1, TNBC2**), together with the downstream
analyses (subclustering, cell-cell communication, trajectory, co-expression
networks, NMF meta-programs, signature-based prediction and an interactive
Shiny explorer).

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
│   ├── ESCC/  PC/  AML/  NSCLC/  TNBC1/  TNBC2/
└── WorkDir/                 # main analysis, one folder per step
    ├── preparing/           #  1. load / QC / merge per cohort
    ├── DoubletFinder/       #  2. doublet detection & removal
    ├── SCT_Integrate/       #  3. SCTransform-based integration
    ├── Log_Integrate/       #  3b. LogNormalize-based integration
    ├── Annotation/          #  4. major-lineage annotation
    ├── Clustering_subtype/  #  5. myeloid / T-NK subclustering
    ├── Description/         #  6. subtype composition & correspondence
    ├── GeneSetEnrichment/   #  7. GSEA / over-representation analysis
    ├── Cellchat/            #  8. cell-cell communication
    ├── NMF/                 #  9. NMF meta-programs
    ├── Monocle/             # 10. pseudotime / trajectory (monocle3)
    ├── hdWGCNA/             # 11. co-expression networks
    ├── Prediction/          # 12. signature-based prediction
    ├── paper_writing/       # 13. main-figure generation
    └── Shiny/               # 14. interactive Shiny application
```

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
Raw data (GEO / in-house)
        |
        v
[1] preparing            load, QC, per-cohort merge
        |
        v
[2] DoubletFinder        detect & remove doublets
        |
        v
[3] SCT_Integrate        SCTransform integration
[3b] Log_Integrate       LogNormalize integration  (alternative strategies)
        |
        v
[4] Annotation           major-lineage annotation
        |
        v
[5] Clustering_subtype   myeloid subclustering (T/NK in parallel)
        |
        +-- [6] Description          composition / correspondence tests
        +-- [7] GeneSetEnrichment    GSEA / ORA
        +-- [8] Cellchat             ligand-receptor communication
        +-- [9] NMF                  meta-program decomposition
        +-- [10] Monocle             pseudotime trajectory
        +-- [11] hdWGCNA             co-expression modules
        +-- [12] Prediction          signature transfer to bulk / other data
                    |
                    v
              [13] paper_writing   main figures
                    |
                    v
              [14] Shiny           interactive explorer
```

Scripts were executed on a SLURM cluster; the accompanying `sbatch*.sh`
files document resource requests and the run order.

---

## Module descriptions

| # | Folder | What it does |
|---|--------|--------------|
| 1 | `preparing/` | Per-cohort loading, quality control, metadata and sample renaming for ESCC, PC, AML, NSCLC, TNBC1, TNBC2. |
| 2 | `DoubletFinder/` | Doublet detection per sample and removal after integration. |
| 3 | `SCT_Integrate/` | SCTransform normalisation, integration and integration-quality evaluation (per cohort + total data). |
| 3b | `Log_Integrate/` | LogNormalize-based integration used as a parallel/validation strategy. |
| 4 | `Annotation/` | Major cell-lineage annotation, marker inspection, and validation on held-out cohorts (TNBC1, NSCLC). |
| 5 | `Clustering_subtype/` | Myeloid subset re-integration, clustering sweeps, marker discovery and subtype annotation; parallel T/NK branch (`T_NK/`, SCTransform / starCAT annotation). |
| 6 | `Description/` | Distribution tests, high-dimensional contingency tables and correspondence analysis of subtypes across cancers. |
| 7 | `GeneSetEnrichment/` | GSEA / ORA input preparation, background-gene handling and radar plots. |
| 8 | `Cellchat/` | CellChat cell-cell communication analysis on the merged object. |
| 9 | `NMF/` | Non-negative matrix factorisation to extract meta-programs (with and without regression of covariates). |
| 10 | `Monocle/` | monocle3 pseudotime / trajectory, branch analysis and branch-specific GSEA. |
| 11 | `hdWGCNA/` | hdWGCNA co-expression networks for ESCC and PC, module gene-set curation and DME analysis. |
| 12 | `Prediction/` | Signature-based prediction of myeloid states in bulk / external datasets (e.g. TCGA ESCC, GSE197677). |
| 13 | `paper_writing/` | Scripts reproducing the main result figures (`Result1`-`Result4`). |
| 14 | `Shiny/` | `ShinyCell2`-based interactive application for exploring the atlas (`shinyApp/`). |

---

## Requirements

Analyses were written in **R** (>= 4.3). Core packages:

**Single-cell** — `Seurat` (v5), `SeuratObject`, `sctransform`, `DoubletFinder`,
`SingleCellExperiment`, `SeuratDisk`, `hdf5r`, `scSHC`.

**Trajectory / networks** — `monocle3`, `hdWGCNA`, `WGCNA`, `igraph`, `ggraph`, `tidygraph`.

**Communication / enrichment** — `CellChat`, `clusterProfiler`, `UCell`.

**Factorisation / stats** — `NMF`, `lme4`, `lmerTest`, `nnet`, `MASS`, `ca`.

**Visualisation** — `ggplot2`, `cowplot`, `patchwork`, `ggrepel`, `ggpubr`,
`RColorBrewer`, `clustree`, `pheatmap`, `gridExtra`, `pdftools`.

**App / IO** — `shiny`, `shinyhelper`, `DT`, `bslib`, `ShinyCell2`, `data.table`.

Shell scripts assume a **SLURM** scheduler. `Python` (via `reticulate`) is
used indirectly by a few tools.

A minimal installation helper is provided at `WorkDir/hdWGCNA/code/install.R`.

---

## Data availability

This repository does **not** ship raw data or large intermediate objects.

- Public datasets were downloaded from **GEO**; the corresponding accession
  numbers are referenced in the scripts (e.g. `GSE198052` for AML,
  `GSE266919` for TNBC2, `GSE197677` and TCGA for the prediction analyses).
- In-house cohorts (ESCC, PC, NSCLC, TNBC1) should be requested from the
  corresponding authors, subject to institutional policy.
- Processed objects are available from the authors on reasonable request,
  or via the data-accession statement of the accompanying paper.

The `*.txt` files under `WorkDir/*/data/**/str(...).txt` are lightweight
`str()` dumps of the `Seurat` objects; they are kept only to document object
structure (cell numbers, assays, metadata columns) and contain no expression
data.

---

## How to run

1. Obtain the raw / processed datasets and place them under `data/<COHORT>/`
   following the paths hard-coded (or configured) in `WorkDir/preparing/code/`.
2. Run the modules **in order** (1 -> 14). Each module folder contains its own
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
<AUTHORS>. <TITLE>. <JOURNAL> (<YEAR>). doi:<DOI>
```

A machine-readable `CITATION.cff` is included and will be updated once the
paper is published.

---

## License

Released under the **MIT License** — see [`LICENSE`](LICENSE).

---

## Contact

**Guosheng Wang** — [@GoulsonWang](https://github.com/GoulsonWang)

Issues and pull requests are welcome.