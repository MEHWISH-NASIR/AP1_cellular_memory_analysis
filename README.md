# AP-1 Cellular Memory: Kinase Analysis

## Overview

This repository contains a reproducible downstream analysis of kinase-associated
transcriptional changes during AP-1-associated cellular memory and
pharmacological interference.

The analysis addresses three main questions:

1. Which kinases respond during the initial dexamethasone (Dex) treatment?
2. Which kinase expression changes are detected in the post-Dex memory contrast?
3. Are these memory-associated kinase changes altered by JNK-IN-8 or T-5224 interference?

Differential-expression modelling is performed using `limma`, with kinase
annotation based on the KinHub/OpenKinome reference list.

---

## Study context

The analysis is based on the public data associated with the AP-1 cellular-memory study:

**Li J, Ravindran PT, O'Farrell A, et al.
AP-1 mediates cellular adaptation and memory formation.
Nature Communications (2026).**

The RNA-seq input files are obtained from:

```text
arjunrajlaboratory/cellularmemory
extractedData/burninrnaseq/202208/
```

This repository represents a downstream computational re-analysis and does not
replace the analyses or conclusions of the original publication.

---

## Reproducible workflow

The complete analysis can be executed using a single command:

```bash
Rscript run_all.R
```

The pipeline runs:

```text
00_fetch_data.R
        ↓
01_setup_and_data.R
        ↓
02_initial_Dex_response.R
        ↓
03_persistent_kinase_memory.R
        ↓
04_kinase_interference.R
        ↓
05_integrated_kinase_trajectory.R
        ↓
06_make_final_kinase_figure.R
```

`00_fetch_data.R` automatically creates the `data/` directory and downloads the
required public input files when they are absent.

---

## Fresh-clone reproduction

Clone the repository and enter the project directory:

```bash
git clone <repository-url>
cd AP1_cellular_memory_analysis
```

If `renv` is not already installed:

```bash
Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org")'
```

Restore the package environment recorded in `renv.lock`:

```bash
Rscript -e 'renv::restore()'
```

Then reproduce the complete analysis:

```bash
Rscript run_all.R
```

The data do not need to be downloaded manually. Missing external inputs are
retrieved automatically by `00_fetch_data.R`.

The project uses `renv` so that the package versions used for the analysis can
be restored from `renv.lock`.

---

## Input data

The workflow requires:

```text
data/log2CPM_matrix.txt
data/studydesign.tsv
data/dex_log2FC.txt
data/dextable.txt
data/KinHubKinaseList.csv
```

The AP-1 RNA-seq files are downloaded from the public cellular-memory repository.

The kinase reference table is downloaded from the OpenKinome/KinHub resource.

The large expression matrix is not required to be committed directly to this
repository because it can be reproduced through the data-fetching script.

---

## 01 — Data setup

Script:

```text
01_setup_and_data.R
```

This stage:

- loads the expression matrix;
- loads the experimental design;
- loads the KinHub kinase reference list;
- aligns metadata with expression-matrix sample columns;
- identifies `Sample_147` as present in the design but absent from the expression matrix;
- retains the 191 samples represented in both datasets;
- verifies that metadata and expression-column order are identical.

The analysis contains 476 measured kinase genes.

---

## 02 — Initial Dex response

Script:

```text
02_initial_Dex_response.R
```

Primary contrast:

```text
Dex - DMSO
```

Among the 476 measured kinases:

- 54 had nominal P < 0.05;
- 4 had kinase-level FDR < 0.05;
- TRIB3 remained significant after genome-wide FDR correction.

The four kinase-level FDR-significant initial responders were:

- TRIB3
- DDR2
- MAPK4
- SGK3

Output:

```text
results/Initial_Dex_kinase_response.csv
```

---

## 03 — Post-Dex memory-associated kinase changes

Script:

```text
03_persistent_kinase_memory.R
```

Primary memory contrast:

```text
pre_T_1 - T_1
```

Among 476 measured kinases:

- 71 had nominal P < 0.05;
- 12 had kinase-level FDR < 0.10;
- 9 had kinase-level FDR < 0.05;
- all 9 also had genome-wide FDR < 0.05.

The nine memory-associated kinase candidates were:

- ACVR1C
- PIM1
- TESK2
- TRIB3
- ERBB3
- SNRK
- CDK2
- STK32A
- PRKCE

Outputs:

```text
results/Memory_all_476_kinases.csv
results/Memory_9_significant_kinases.csv
```

---

## 04 — Pharmacological interference

Script:

```text
04_kinase_interference.R
```

Memory effects are evaluated under:

- no inhibitor;
- JNK-IN-8;
- T-5224.

Interaction contrasts test whether the memory effect changes under each
inhibitor relative to the no-inhibitor condition.

Among the nine memory-associated kinase candidates:

- four JNK-IN-8 interactions had nominal P < 0.05;
- one T-5224 interaction had nominal P < 0.05;
- none remained significant after candidate-level FDR correction.

These inhibitor-related signals are therefore treated as suggestive rather than
definitive.

Output:

```text
results/Memory_9_kinases_interference.csv
```

---

## 05 — Integrated kinase trajectories

Script:

```text
05_integrated_kinase_trajectory.R
```

This script explicitly reads the intermediate result CSV files generated by
scripts 02–04 and therefore does not require objects to already exist in the R
session.

The nine memory-associated kinases are classified as:

- 7 memory-emergent candidates;
- 1 initial nominal responder with the same direction of change: TESK2;
- 1 significant initial responder with direction reversal in memory: TRIB3.

Output:

```text
results/Integrated_9_kinase_trajectory.csv
```

---

## 06 — Final kinase figure

Script:

```text
06_make_final_kinase_figure.R
```

This script directly reads:

```text
results/Integrated_9_kinase_trajectory.csv
```

and can therefore run in a fresh R session.

The final figure compares:

- Initial Dex
- Memory
- Memory + JNK-IN-8
- Memory + T-5224

Outputs:

```text
figures/Final_Kinase_Memory_Interference.png
figures/Final_Kinase_Memory_Interference.pdf
```

---

## Main interpretation

Nine kinases show significant expression differences in the post-Dex memory
contrast.

Most are not significant initial Dex responders, indicating that the
post-treatment memory-associated transcriptional state is not simply a
continuation of the acute Dex response.

TRIB3 is distinct because it shows a strong initial Dex response followed by an
opposite-direction effect in the memory contrast.

Several nominal inhibitor-interaction signals are observed for JNK-IN-8 and
T-5224, but none survive candidate-level FDR correction. These effects should
therefore be interpreted cautiously.

---

## HTML analysis report

A rendered analysis report is included at:

```text
docs/index.html
```

The source report is:

```text
analysis_report.Rmd
```

To rebuild the HTML report:

```bash
Rscript render_report.R
```

The report contains the major kinase result tables, trajectory classifications,
interference results, final figure, and software-environment information.

R Markdown renders `.Rmd` documents to HTML using `rmarkdown::render()` and
Pandoc. The report renderer writes the final document into the `docs/`
directory.

---

## Software environment

The R environment is managed using:

```text
renv.lock
```

The lockfile records the package environment used by the project, including the
main analysis packages:

```text
limma
dplyr
ggplot2
```

and the report-generation packages such as:

```text
rmarkdown
knitr
```

Check environment consistency with:

```bash
Rscript -e 'renv::status()'
```

Restore it with:

```bash
Rscript -e 'renv::restore()'
```

---

## Repository structure

```text
AP1_cellular_memory_analysis/
│
├── 00_fetch_data.R
├── 01_setup_and_data.R
├── 02_initial_Dex_response.R
├── 03_persistent_kinase_memory.R
├── 04_kinase_interference.R
├── 05_integrated_kinase_trajectory.R
├── 06_make_final_kinase_figure.R
│
├── run_all.R
├── render_report.R
├── analysis_report.Rmd
│
├── renv.lock
├── .Rprofile
├── renv/
│
├── data/
├── results/
├── figures/
│
├── docs/
│   └── index.html
│
├── README.md
└── .gitignore
```

---

## Reproducibility status

The workflow has been tested by:

- automatically fetching the input data into an empty `data/` directory;
- running `01_setup_and_data.R` without pre-existing input data;
- running scripts 05 and 06 in fresh R sessions;
- executing the complete workflow with `Rscript run_all.R`;
- verifying the nine integrated kinase results after refactoring;
- confirming numerical agreement with the previously generated integrated table;
- running the complete pipeline under the `renv` project environment;
- rendering the complete HTML analysis report.

The project can therefore be reproduced from the documented input sources,
analysis scripts, package lockfile, and single pipeline entry point.
