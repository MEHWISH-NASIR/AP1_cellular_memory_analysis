library(limma)
#initial-response samples
initial_conditions <- c(
  "DMSO_1",
  "Dex"
)

initial_design <- design_analysis[
  design_analysis$condition %in% initial_conditions,
]

initial_design <- initial_design[
  match(
    c(
      paste0("Sample_", 1:3),
      paste0("Sample_", 4:6)
    ),
    initial_design$samples
  ),
]
# Check selected samples
initial_design[
  ,
  c(
    "samples",
    "condition",
    "pretreat",
    "treatment"
  )
]
#Extract expression matrix
expr_initial <- expr[
  ,
  c(
    "geneID",
    initial_design$samples
  )
]

# Gene symbols become row names
expr_initial_mat <- as.matrix(
  expr_initial[, -1]
)
rownames(expr_initial_mat) <- expr_initial$geneID
# Validate dimensions
dim(expr_initial_mat)

stopifnot(
  ncol(expr_initial_mat) == nrow(initial_design),
  identical(
    colnames(expr_initial_mat),
    initial_design$samples
  )
)
#Construct limma design
group_initial <- factor(
  initial_design$condition,
  levels = c(
    "DMSO_1",
    "Dex"
  )
)
design_initial <- model.matrix(
  ~0 + group_initial
)
colnames(design_initial) <- c(
  "DMSO",
  "Dex"
)
design_initial
initial_design[
  ,
  c(
    "samples",
    "condition",
    "pretreat",
    "treatment"
  )
]

dim(expr_initial_mat)

design_initial
#Fit initial Dex response
fit_initial <- lmFit(
  expr_initial_mat,
  design_initial
)
contrast_initial <- makeContrasts(
  Dex_vs_DMSO = Dex - DMSO,
  levels = design_initial
)
fit_initial2 <- contrasts.fit(
  fit_initial,
  contrast_initial
)
fit_initial2 <- eBayes(
  fit_initial2
)
# Extract genome-wide results
dex_initial_all <- topTable(
  fit_initial2,
  coef = "Dex_vs_DMSO",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
dim(dex_initial_all)
head(dex_initial_all, 10)
#Define kinase symbols
kinase_symbols <- unique(
  kinase_list$HGNC.Name
)
kinase_symbols <- kinase_symbols[
  !is.na(kinase_symbols) &
    kinase_symbols != ""
]
length(kinase_symbols)
# Kinases represented in expression dataset
measured_kinases <- intersect(
  kinase_symbols,
  rownames(dex_initial_all)
)
length(measured_kinases)
#Extract initial Dex kinase response
dex_kinases_all <- dex_initial_all[
  rownames(dex_initial_all) %in% measured_kinases,
]
dex_kinases_all <- dex_kinases_all[
  order(dex_kinases_all$P.Value),
]
# Kinase-specific BH correction
dex_kinases_all$kinase_FDR <- p.adjust(
  dex_kinases_all$P.Value,
  method = "BH"
)
#Summarise initial response
initial_summary <- c(
  measured_kinases = nrow(dex_kinases_all),
  
  raw_P_lt_0.05 =
    sum(dex_kinases_all$P.Value < 0.05),
  
  kinase_FDR_lt_0.10 =
    sum(dex_kinases_all$kinase_FDR < 0.10),
  
  kinase_FDR_lt_0.05 =
    sum(dex_kinases_all$kinase_FDR < 0.05),
  
  genome_FDR_lt_0.05 =
    sum(dex_kinases_all$adj.P.Val < 0.05)
)
initial_summary
#Show strongest kinase responses
head(
  dex_kinases_all[
    ,
    c(
      "logFC",
      "AveExpr",
      "P.Value",
      "adj.P.Val",
      "kinase_FDR"
    )
  ],
  20
)


#optional 
initial_summary
head(
  dex_kinases_all[
    ,
    c("logFC", "P.Value", "adj.P.Val", "kinase_FDR")
  ],
  10
)

# save results
initial_kinase_results <- data.frame(
  kinase = rownames(dex_kinases_all),
  dex_kinases_all,
  row.names = NULL
)
write.csv(
  initial_kinase_results,
  "results/Initial_Dex_kinase_response.csv",
  row.names = FALSE
)
cat(
  "Saved:",
  "results/Initial_Dex_kinase_response.csv",
  "\n"
)
file.exists(
  "results/Initial_Dex_kinase_response.csv"
)