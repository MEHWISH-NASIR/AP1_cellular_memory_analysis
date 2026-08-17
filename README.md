# AP-1 Cellular Memory: Kinase Analysis

## Overview

This repository contains an analysis of kinase-associated transcriptional changes during drug-induced cellular memory and pharmacological interference.

The analysis focuses on three questions:

1. Which kinases respond during the initial dexamethasone (Dex) treatment?
2. Which kinase expression changes are detected in the post-Dex memory contrast?
3. Are these memory-associated changes altered by JNK-IN-8 or T-5224 interference?

The analysis was performed using the study expression matrix and experimental design, with differential-expression modelling implemented using `limma`.

---

## Analysis workflow

The analysis is organized into six scripts:

### 01 — Data setup

`01_setup_and_data.R`

- Loads the expression matrix and experimental design.
- Loads the KinHub kinase reference list.
- Matches expression samples to the experimental design.
- Identifies one metadata sample (`Sample_147`) that is absent from the expression matrix.
- Restricts downstream analysis to the 191 samples represented in both datasets.
- Identifies 476 measured kinases.

### 02 — Initial Dex response

`02_initial_Dex_response.R`

Compares:

`Dex - DMSO`

Among the 476 measured kinases:

- 54 had nominal P < 0.05.
- 4 had kinase-level FDR < 0.05.
- TRIB3 remained significant after genome-wide FDR correction.

The four kinase-level FDR-significant initial responders were:

- TRIB3
- DDR2
- MAPK4
- SGK3

### 03 — Post-Dex memory-associated kinase changes

`03_persistent_kinase_memory.R`

The primary memory contrast is:

`pre_T_1 - T_1`

The six relevant treatment/interference conditions were modelled together to obtain a common empirical-Bayes variance estimate.

Among 476 measured kinases:

- 71 had nominal P < 0.05.
- 12 had kinase-level FDR < 0.10.
- 9 had kinase-level FDR < 0.05.
- All 9 also had genome-wide FDR < 0.05.

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

### 04 — Pharmacological interference

`04_kinase_interference.R`

Memory effects were evaluated under:

- no inhibitor
- JNK-IN-8
- T-5224

Interaction contrasts tested whether the memory effect changed under each inhibitor relative to the no-inhibitor condition.

For the nine memory-associated kinase candidates:

- Four JNK-IN-8 interactions had nominal P < 0.05.
- One T-5224 interaction had nominal P < 0.05.
- None of these interaction effects remained significant after FDR correction across the nine candidates.

These results are therefore treated as suggestive rather than definitive evidence of inhibitor-dependent modulation.

### 05 — Integrated kinase trajectories

`05_integrated_kinase_trajectory.R`

Initial Dex response, post-Dex memory effect, and interference results were combined into a single table.

The nine memory-associated kinases were classified as:

- 7 memory-emergent candidates
- 1 initial nominal responder with the same direction of change (TESK2)
- 1 significant initial responder showing direction reversal in the memory contrast (TRIB3)

### 06 — Final figure

`06_make_final_kinase_figure.R`

Generates the final effect-size profiles for the nine memory-associated kinases across:

- Initial Dex
- Memory
- Memory + JNK-IN-8
- Memory + T-5224

The figure is provided in PNG and PDF formats.

---

## Main interpretation

The analysis identifies a set of nine kinases with significant expression differences in the post-Dex memory contrast.

Most of these kinases were not significant initial Dex responders, suggesting that the post-treatment state is not simply a continuation of the acute Dex transcriptional response.

TRIB3 is distinct because it shows a strong initial Dex response followed by an opposite-direction effect in the memory contrast.

JNK-IN-8 and T-5224 produce several nominal interaction signals, but none survive candidate-level FDR correction. Therefore, inhibitor-dependent effects should be interpreted cautiously.

---

## Repository structure

```text
AP1_cellular_memory_analysis/
|
|-- 01_setup_and_data.R
|-- 02_initial_Dex_response.R
|-- 03_persistent_kinase_memory.R
|-- 04_kinase_interference.R
|-- 05_integrated_kinase_trajectory.R
|-- 06_make_final_kinase_figure.R
|
|-- data/
|
|-- results/
|
|-- figures/
|
|-- README.md
|-- .gitignore