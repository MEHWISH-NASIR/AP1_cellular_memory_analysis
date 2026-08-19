
cat(
  "\n",
  "============================================================\n",
  " AP-1 cellular-memory kinase analysis\n",
  "============================================================\n\n",
  sep = ""
)
# Define analysis workflow

pipeline_scripts <- c(
  "00_fetch_data.R",
  "01_setup_and_data.R",
  "02_initial_Dex_response.R",
  "03_persistent_kinase_memory.R",
  "04_kinase_interference.R",
  "05_integrated_kinase_trajectory.R",
  "06_make_final_kinase_figure.R"
)

# Validate that all scripts are available


missing_scripts <- pipeline_scripts[
  !file.exists(pipeline_scripts)
]

if (length(missing_scripts) > 0) {
  stop(
    "Required pipeline script(s) missing: ",
    paste(missing_scripts, collapse = ", ")
  )
}
# Ensure output directories exist
dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)

# Run analysis sequentially

for (script in pipeline_scripts) {

  cat(
    "\n",
    "------------------------------------------------------------\n",
    "Running: ", script, "\n",
    "------------------------------------------------------------\n",
    sep = ""
  )

  start_time <- Sys.time()

  tryCatch(

    {
      source(
        script,
        echo = FALSE,
        chdir = FALSE
      )
    },

    error = function(e) {

      cat(
        "\nPIPELINE FAILED\n",
        "Script: ", script, "\n",
        "Error: ", conditionMessage(e), "\n",
        sep = ""
      )

      stop(
        "Analysis stopped while running ",
        script,
        call. = FALSE
      )
    }
  )

  elapsed <- difftime(
    Sys.time(),
    start_time,
    units = "secs"
  )

  cat(
    "Completed: ",
    script,
    " [",
    round(as.numeric(elapsed), 1),
    " sec]\n",
    sep = ""
  )
}

# Validate key final outputs

final_outputs <- c(
  "results/Integrated_9_kinase_trajectory.csv",
  "figures/Final_Kinase_Memory_Interference.png",
  "figures/Final_Kinase_Memory_Interference.pdf"
)

missing_outputs <- final_outputs[
  !file.exists(final_outputs)
]

if (length(missing_outputs) > 0) {
  stop(
    "Pipeline completed but expected output file(s) are missing: ",
    paste(missing_outputs, collapse = ", ")
  )
}

# Pipeline completion

cat(
  "\n",
  "============================================================\n",
  " PIPELINE COMPLETED SUCCESSFULLY\n",
  "============================================================\n",
  "Final outputs:\n",
  paste0("  - ", final_outputs, collapse = "\n"),
  "\n\n",
  sep = ""
)
