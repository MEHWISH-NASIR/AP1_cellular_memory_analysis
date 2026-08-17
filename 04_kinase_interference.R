#Test whether the kinase-memory response identified
library(limma)
#Define interference contrasts
contrast_interference <- makeContrasts(
  
  Memory_No_Inhibitor =
    preT1 - T1,
  
  Memory_JNKIN8 =
    preT8 - T8,
  
  Memory_T5224 =
    preT5224 - T5224,
  
  JNK_Interaction =
    (preT8 - T8) -
    (preT1 - T1),
  
  T5224_Interaction =
    (preT5224 - T5224) -
    (preT1 - T1),
  
  levels = design_memory
)
contrast_interference
#Apply interference contrasts
fit_interference <- contrasts.fit(
  fit_memory,
  contrast_interference
)

fit_interference <- eBayes(
  fit_interference
)
memory_no_inhib <- topTable(
  fit_interference,
  coef = "Memory_No_Inhibitor",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

memory_jnk <- topTable(
  fit_interference,
  coef = "Memory_JNKIN8",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)
memory_t5224 <- topTable(
  fit_interference,
  coef = "Memory_T5224",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)
jnk_interaction <- topTable(
  fit_interference,
  coef = "JNK_Interaction",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)
jnk_interaction <- topTable(
  fit_interference,
  coef = "JNK_Interaction",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)
t5224_interaction <- topTable(
  fit_interference,
  coef = "T5224_Interaction",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)
memory_no_inhib[
  rownames(persistent_kinases_FDR05),
  c("logFC", "P.Value", "adj.P.Val")
]
memory9 <- rownames(
  persistent_kinases_FDR05
)
memory9_interference <- data.frame(
  
  kinase = memory9,
  
  No_Inhib_logFC =
    memory_no_inhib[memory9, "logFC"],
  
  No_Inhib_P =
    memory_no_inhib[memory9, "P.Value"],
  
  JNKIN8_logFC =
    memory_jnk[memory9, "logFC"],
  
  JNKIN8_P =
    memory_jnk[memory9, "P.Value"],
  
  T5224_logFC =
    memory_t5224[memory9, "logFC"],
  T5224_P =
    memory_t5224[memory9, "P.Value"],
  
  JNK_interaction =
    jnk_interaction[memory9, "logFC"],
  
  JNK_interaction_P =
    jnk_interaction[memory9, "P.Value"],
  
  T5224_interaction =
    t5224_interaction[memory9, "logFC"],
  
  T5224_interaction_P =
    t5224_interaction[memory9, "P.Value"]
)
#Candidate-set FDR
memory9_interference$JNK_candidate_FDR <- p.adjust(
  memory9_interference$JNK_interaction_P,
  method = "BH"
)

memory9_interference$T5224_candidate_FDR <- p.adjust(
  memory9_interference$T5224_interaction_P,
  method = "BH"
)
memory9_interference
write.csv(
  memory9_interference,
  "results/Memory_9_kinases_interference.csv",
  row.names = FALSE
)

file.exists(
  "results/Memory_9_kinases_interference.csv"
)