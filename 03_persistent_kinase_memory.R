library(limma)
#Select memory/interference samples
memory_conditions <- c(
  "T_1",
  "pre_T_1",
  "T_8_1",
  "pre_T_8_1",
  "T_5224_1",
  "pre_T_5224_1"
)
memory_design <- design_analysis[
  design_analysis$condition %in% memory_conditions,
]
memory_design[
  ,
  c(
    "samples",
    "condition",
    "pretreat",
    "treatment"
  )
]
table(memory_design$condition)
#Build expression matrix
expr_memory <- expr[
  ,
  c(
    "geneID",
    memory_design$samples
  )
]
expr_memory_mat <- as.matrix(
  expr_memory[, -1]
)
rownames(expr_memory_mat) <- expr_memory$geneID
# Validate sample allignment
stopifnot(
  ncol(expr_memory_mat) == nrow(memory_design),
  identical(
    colnames(expr_memory_mat),
    memory_design$samples
  )
)

dim(expr_memory_mat)
#Build six-condition design matrix
group_memory <- factor(
  memory_design$condition,
  levels = c(
    "T_1",
    "pre_T_1",
    "T_8_1",
    "pre_T_8_1",
    "T_5224_1",
    "pre_T_5224_1"
  )
)
design_memory <- model.matrix(
  ~0 + group_memory
)
colnames(design_memory) <- c(
  "T1",
  "preT1",
  "T8",
  "preT8",
  "T5224",
  "preT5224"
)
design_memory

#Fit joint memory model
fit_memory <- lmFit(
  expr_memory_mat,
  design_memory
)
contrast_memory <- makeContrasts(
  Memory_No_Inhibitor = preT1 - T1,
  levels = design_memory
)
contrast_memory
fit_memory2 <- contrasts.fit(
  fit_memory,
  contrast_memory
)
fit_memory2 <- eBayes(
  fit_memory2
)
# Check model statistics
fit_memory2$df.residual[1]
fit_memory2$df.prior
fit_memory2$s2.prior
#Extract genome-wide memory results
memory_all <- topTable(
  fit_memory2,
  coef = "Memory_No_Inhibitor",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
dim(memory_all)
head(memory_all, 10)
memory_kinases <- memory_all[
  rownames(memory_all) %in% measured_kinases,
]
memory_kinases <- memory_kinases[
  order(memory_kinases$P.Value),
]
# Kinase-panel FDR
memory_kinases$kinase_FDR <- p.adjust(
  memory_kinases$P.Value,
  method = "BH"
)
memory_summary <- c(
  kinase_total =
    nrow(memory_kinases),
  
  raw_P_lt_0.05 =
    sum(memory_kinases$P.Value < 0.05),
  
  kinase_FDR_lt_0.10 =
    sum(memory_kinases$kinase_FDR < 0.10),
  
  kinase_FDR_lt_0.05 =
    sum(memory_kinases$kinase_FDR < 0.05),
  
  genome_FDR_lt_0.05 =
    sum(memory_kinases$adj.P.Val < 0.05)
)

memory_summary
#Significant kinase-memory candidates
persistent_kinases_FDR05 <- memory_kinases[
  memory_kinases$kinase_FDR < 0.05,
]
persistent_kinases_FDR05[
  ,
  c(
    "logFC",
    "AveExpr",
    "P.Value",
    "adj.P.Val",
    "kinase_FDR"
  )
]
memory_kinase_results <- data.frame(
  kinase = rownames(memory_kinases),
  memory_kinases,
  row.names = NULL
)
write.csv(
  memory_kinase_results,
  "results/Memory_all_476_kinases.csv",
  row.names = FALSE
)
persistent_kinase_results <- data.frame(
  kinase = rownames(persistent_kinases_FDR05),
  persistent_kinases_FDR05,
  row.names = NULL
)
write.csv(
  persistent_kinase_results,
  "results/Memory_9_significant_kinases.csv",
  row.names = FALSE
)
file.exists(
  "results/Memory_all_476_kinases.csv"
)

file.exists(
  "results/Memory_9_significant_kinases.csv"
)