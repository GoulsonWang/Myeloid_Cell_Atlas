# Archive

This folder contains scripts that were used **during development** of the
myeloid cell atlas but are **not part of the final, published pipeline**.

They are kept here (rather than deleted) so that the analysis history stays
traceable and so that reviewers can see the experiments that were run.

## What is in here

| File | Reason for archiving |
|------|----------------------|
| `WorkDir/Description/code/test.R` | ad-hoc test script (superseded by `Description.R`) |
| `WorkDir/NMF/code/test_nmf.R` | NMF parameter trial (superseded by `NR_Myeloid.R`) |
| `WorkDir/SCT_Integrate/code/ESCC/test1.R` | integration trial |
| `WorkDir/SCT_Integrate/code/ESCC/test_Fun.R` | function prototyping |
| `WorkDir/SCT_Integrate/code/TNBC2/test1.R` | integration trial |
| `WorkDir/SCT_Integrate/code/PC/MT&HB_test.R` | QC threshold trial |
| `WorkDir/Log_Integrate/code/TotalData/IntegrateLayers_Test.R` | layer-integration trial |
| `WorkDir/preparing/code/TNBC1/ResultTest_Rename.R` | cell-renaming trial |
| `WorkDir/preparing/code/PC/load(single sample test).R` | single-sample load test |
| `WorkDir/Clustering_subtype/code/Myeloid/Test_Orig_nCells.R` | cell-count trial |
| `WorkDir/Clustering_subtype/code/T_NK/starcat_annotation_test.R` | annotation trial |
| `WorkDir/Shiny/install_test.R` | dependency install test |
| `WorkDir/hdWGCNA/code/learn_example/learn_example.R` | tutorial / learning script |
| `WorkDir/Shiny/build_shiny_failed.txt` | failed deploy log |
| `WorkDir/SCT_Integrate/data/ESCC/str(test1_SCT).txt` | object dump of a trial run |
| `WorkDir/SCT_Integrate/data/ESCC/str(test1_old).txt` | object dump of a trial run |

The original relative paths are preserved inside `archive/`, so any file
can be traced back to the folder it came from.