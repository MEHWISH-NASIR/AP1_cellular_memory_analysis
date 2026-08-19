# ============================================================
# 01_setup_and_data.R
# Load and validate inputs for AP-1 cellular-memory analysis
# ============================================================


# ------------------------------------------------------------
# Ensure required input data are available
# ------------------------------------------------------------

required_inputs <- c(
  "data/log2CPM_matrix.txt",
  "data/studydesign.tsv",
  "data/KinHubKinaseList.csv"
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs) > 0) {

  message(
    "Required input file(s) missing: ",
    paste(basename(missing_inputs), collapse = ", ")
  )

  if (!file.exists("00_fetch_data.R")) {
    stop(
      "00_fetch_data.R was not found. ",
      "Cannot automatically retrieve required input data."
    )
  }

  message("Running 00_fetch_data.R ...")
  source("00_fetch_data.R")
}


# Final input validation

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs) > 0) {
  stop(
    "Required input file(s) still missing: ",
    paste(missing_inputs, collapse = ", ")
  )
}


# ------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------

library(limma)
library(dplyr)


# ------------------------------------------------------------
# Load expression data
# ------------------------------------------------------------

expr <- read.delim(
  "data/log2CPM_matrix.txt",
  check.names = FALSE
)

dim(expr)
head(expr[, 1:5])


# ------------------------------------------------------------
# Load experimental design
# ------------------------------------------------------------

design <- read.delim(
  "data/studydesign.tsv",
  stringsAsFactors = FALSE
)

dim(design)
head(design)


# ------------------------------------------------------------
# Load kinase reference list
# ------------------------------------------------------------

kinase_list <- read.csv(
  "data/KinHubKinaseList.csv",
  stringsAsFactors = FALSE
)

dim(kinase_list)
head(kinase_list)

dim(expr)
dim(design)
dim(kinase_list)


# ------------------------------------------------------------
# Validate sample alignment
# ------------------------------------------------------------

sample_cols <- colnames(expr)[
  colnames(expr) != "geneID"
]

# Identify samples present in design but absent from expression matrix
missing_samples <- setdiff(
  design$samples,
  sample_cols
)

cat(
  "Samples present in design but absent from expression matrix:",
  paste(missing_samples, collapse = ", "),
  "\n"
)


# Keep only samples represented in expression matrix

design_analysis <- design[
  design$samples %in% sample_cols,
]

# Match metadata order exactly to expression columns

design_analysis <- design_analysis[
  match(sample_cols, design_analysis$samples),
]


# ------------------------------------------------------------
# Final sample validation
# ------------------------------------------------------------

stopifnot(
  nrow(design_analysis) == length(sample_cols),
  identical(design_analysis$samples, sample_cols)
)

cat(
  "Expression samples:", length(sample_cols), "\n",
  "Analysis design samples:", nrow(design_analysis), "\n",
  "Sample order identical:",
  identical(design_analysis$samples, sample_cols),
  "\n"
)
