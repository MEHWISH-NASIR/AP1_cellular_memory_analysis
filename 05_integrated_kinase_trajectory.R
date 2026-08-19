
required_result_files <- c(
  "results/Initial_Dex_kinase_response.csv",
  "results/Memory_9_significant_kinases.csv",
  "results/Memory_9_kinases_interference.csv"
)

missing_result_files <- required_result_files[
  !file.exists(required_result_files)
]

if (length(missing_result_files) > 0) {
  stop(
    "Required intermediate result file(s) missing: ",
    paste(missing_result_files, collapse = ", "),
    ". Run scripts 02-04 first or use run_all.R."
  )
}


initial_kinase_results <- read.csv(
  "results/Initial_Dex_kinase_response.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

persistent_kinase_results <- read.csv(
  "results/Memory_9_significant_kinases.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

interference_results <- read.csv(
  "results/Memory_9_kinases_interference.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)



# Validate required columns

required_initial_columns <- c(
  "kinase",
  "logFC",
  "P.Value",
  "adj.P.Val",
  "kinase_FDR"
)

required_memory_columns <- c(
  "kinase",
  "logFC",
  "P.Value",
  "adj.P.Val",
  "kinase_FDR"
)

required_interference_columns <- c(
  "kinase",
  "JNKIN8_logFC",
  "JNK_interaction",
  "JNK_interaction_P",
  "T5224_logFC",
  "T5224_interaction",
  "T5224_interaction_P"
)

stopifnot(
  all(required_initial_columns %in% colnames(initial_kinase_results)),
  all(required_memory_columns %in% colnames(persistent_kinase_results)),
  all(required_interference_columns %in% colnames(interference_results))
)

stopifnot(
  !anyDuplicated(initial_kinase_results$kinase),
  !anyDuplicated(persistent_kinase_results$kinase),
  !anyDuplicated(interference_results$kinase)
)
# Reconstruct objects required by the original script


dex_kinases_all <- initial_kinase_results
rownames(dex_kinases_all) <- dex_kinases_all$kinase
dex_kinases_all$kinase <- NULL


persistent_kinases_FDR05 <- persistent_kinase_results
rownames(persistent_kinases_FDR05) <-
  persistent_kinases_FDR05$kinase
persistent_kinases_FDR05$kinase <- NULL

memory_no_inhib <- persistent_kinases_FDR05


rownames(interference_results) <-
  interference_results$kinase


memory_jnk <- data.frame(
  logFC = interference_results$JNKIN8_logFC,
  row.names = rownames(interference_results)
)

memory_t5224 <- data.frame(
  logFC = interference_results$T5224_logFC,
  row.names = rownames(interference_results)
)

jnk_interaction <- data.frame(
  logFC = interference_results$JNK_interaction,
  P.Value = interference_results$JNK_interaction_P,
  row.names = rownames(interference_results)
)

t5224_interaction <- data.frame(
  logFC = interference_results$T5224_interaction,
  P.Value = interference_results$T5224_interaction_P,
  row.names = rownames(interference_results)
)

cat(
  "Loaded intermediate kinase results successfully.\n"
)


















memory9 <- rownames(
  persistent_kinases_FDR05
)

memory9
#Check candidates are present in all required result tables
stopifnot(
  all(memory9 %in% rownames(dex_kinases_all)),
  all(memory9 %in% rownames(memory_no_inhib)),
  all(memory9 %in% rownames(memory_jnk)),
  all(memory9 %in% rownames(memory_t5224)),
  all(memory9 %in% rownames(jnk_interaction)),
  all(memory9 %in% rownames(t5224_interaction))
)
#Build integrated table
integrated_kinase_table <- data.frame(
  
  kinase = memory9,
  
  # Initial Dex response
  Initial_Dex_logFC =
    dex_kinases_all[memory9, "logFC"],
  
  Initial_Dex_P =
    dex_kinases_all[memory9, "P.Value"],
  
  Initial_Dex_FDR =
    dex_kinases_all[memory9, "adj.P.Val"],
  
  Initial_Dex_kinaseFDR =
    dex_kinases_all[memory9, "kinase_FDR"],
  # Memory response
  Memory_logFC =
    memory_no_inhib[memory9, "logFC"],
  
  Memory_P =
    memory_no_inhib[memory9, "P.Value"],
  
  Memory_FDR =
    memory_no_inhib[memory9, "adj.P.Val"],
  
  Memory_kinaseFDR =
    persistent_kinases_FDR05[memory9, "kinase_FDR"],
  # JNK-IN-8
  JNKIN8_Memory_logFC =
    memory_jnk[memory9, "logFC"],
  
  JNK_Interaction =
    jnk_interaction[memory9, "logFC"],
  
  JNK_Interaction_P =
    jnk_interaction[memory9, "P.Value"],
  # T-5224
  T5224_Memory_logFC =
    memory_t5224[memory9, "logFC"],
  
  T5224_Interaction =
    t5224_interaction[memory9, "logFC"],
  
  T5224_Interaction_P =
    t5224_interaction[memory9, "P.Value"]
)
integrated_kinase_table$JNK_candidate_FDR <-
  p.adjust(
    integrated_kinase_table$JNK_Interaction_P,
    method = "BH"
  )

integrated_kinase_table$T5224_candidate_FDR <-
  p.adjust(
    integrated_kinase_table$T5224_Interaction_P,
    method = "BH"
  )
integrated_kinase_table
integrated_kinase_table$Initial_status <- ifelse(
  integrated_kinase_table$Initial_Dex_kinaseFDR < 0.05,
  "Significant",
  "Not significant"
)
integrated_kinase_table$Trajectory <- ifelse(
  
  integrated_kinase_table$Initial_Dex_kinaseFDR < 0.05 &
    sign(integrated_kinase_table$Initial_Dex_logFC) ==
    sign(integrated_kinase_table$Memory_logFC),
  
  "Persistent_same_direction",
  
  ifelse(
    
    integrated_kinase_table$Initial_Dex_kinaseFDR < 0.05 &
      sign(integrated_kinase_table$Initial_Dex_logFC) !=
      sign(integrated_kinase_table$Memory_logFC),
    
    "Persistent_direction_reversal",
    
    "Memory_emergent"
  )
)
integrated_kinase_table[
  ,
  c(
    "kinase",
    "Initial_Dex_logFC",
    "Initial_Dex_kinaseFDR",
    "Memory_logFC",
    "Memory_kinaseFDR",
    "Initial_status",
    "Trajectory"
  )
]
integrated_kinase_table$Trajectory <- ifelse(
  
  integrated_kinase_table$Initial_Dex_kinaseFDR < 0.05 &
    sign(integrated_kinase_table$Initial_Dex_logFC) ==
    sign(integrated_kinase_table$Memory_logFC),
  
  "Initial_responder_same_direction",
  
  ifelse(
    
    integrated_kinase_table$Initial_Dex_kinaseFDR < 0.05 &
      sign(integrated_kinase_table$Initial_Dex_logFC) !=
      sign(integrated_kinase_table$Memory_logFC),
    
    "Initial_responder_direction_reversal",
    
    ifelse(
      
      integrated_kinase_table$Initial_Dex_P < 0.05 &
        sign(integrated_kinase_table$Initial_Dex_logFC) ==
        sign(integrated_kinase_table$Memory_logFC),
      
      "Initial_nominal_same_direction",
      
      "Memory_emergent"
    )
  )
)
integrated_kinase_table[
  ,
  c(
    "kinase",
    "Initial_Dex_logFC",
    "Initial_Dex_P",
    "Initial_Dex_kinaseFDR",
    "Memory_logFC",
    "Memory_kinaseFDR",
    "Trajectory"
  )
]
integrated_kinase_table$JNK_interference_status <- ifelse(
  integrated_kinase_table$JNK_candidate_FDR < 0.05,
  "FDR_significant",
  ifelse(
    integrated_kinase_table$JNK_Interaction_P < 0.05,
    "Nominal_only",
    "Not_significant"
  )
)

integrated_kinase_table$T5224_interference_status <- ifelse(
  integrated_kinase_table$T5224_candidate_FDR < 0.05,
  "FDR_significant",
  ifelse(
    integrated_kinase_table$T5224_Interaction_P < 0.05,
    "Nominal_only",
    "Not_significant"
  )
)
write.csv(
  integrated_kinase_table,
  "results/Integrated_9_kinase_trajectory.csv",
  row.names = FALSE
)
table(integrated_kinase_table$Trajectory)

table(integrated_kinase_table$JNK_interference_status)

table(integrated_kinase_table$T5224_interference_status)

file.exists(
  "results/Integrated_9_kinase_trajectory.csv"
)