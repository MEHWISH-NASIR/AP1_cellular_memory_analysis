library(limma)
library(dplyr)
#Load expression data
expr <- read.delim(
  "data/log2CPM_matrix.txt",
  check.names = FALSE
)

dim(expr)
head(expr[, 1:5])
# load experimental design
design <- read.delim(
  "data/studydesign.tsv",
  stringsAsFactors = FALSE
)

dim(design)
head(design)
#Load kinase reference list
kinase_list <- read.csv(
  "data/KinHubKinaseList.csv",
  stringsAsFactors = FALSE
)

dim(kinase_list)
head(kinase_list)
dim(expr)
dim(design)
dim(kinase_list)
#Validate sample alignment
sample_cols <- colnames(expr)[
  colnames(expr) != "geneID"
]

# Identify samples missing from expression matrix
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

# Final validation
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