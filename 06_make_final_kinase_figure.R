library(ggplot2)

integrated_file <- "results/Integrated_9_kinase_trajectory.csv"

if (!file.exists(integrated_file)) {
  stop(
    "Required result file is missing: ",
    integrated_file,
    ". Run 05_integrated_kinase_trajectory.R first or use run_all.R."
  )
}

integrated_kinase_table <- read.csv(
  integrated_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_columns <- c(
  "kinase",
  "Initial_Dex_logFC",
  "Memory_logFC",
  "JNKIN8_Memory_logFC",
  "T5224_Memory_logFC",
  "Trajectory"
)

stopifnot(
  all(required_columns %in% colnames(integrated_kinase_table)),
  !anyDuplicated(integrated_kinase_table$kinase)
)

if (nrow(integrated_kinase_table) == 0) {
  stop("Integrated kinase table contains no rows.")
}

cat(
  "Loaded integrated kinase trajectory successfully.\n"
)

cat(
  "Kinases loaded:",
  nrow(integrated_kinase_table),
  "\n"
)

dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)

plot_data <- data.frame(

  kinase = rep(
    integrated_kinase_table$kinase,
    times = 4
  ),

  stage = rep(
    c(
      "Initial Dex",
      "Memory",
      "Memory + JNK-IN-8",
      "Memory + T-5224"
    ),
    each = nrow(integrated_kinase_table)
  ),

  logFC = c(
    integrated_kinase_table$Initial_Dex_logFC,
    integrated_kinase_table$Memory_logFC,
    integrated_kinase_table$JNKIN8_Memory_logFC,
    integrated_kinase_table$T5224_Memory_logFC
  )
)

plot_data$stage <- factor(
  plot_data$stage,
  levels = c(
    "Initial Dex",
    "Memory",
    "Memory + JNK-IN-8",
    "Memory + T-5224"
  )
)

plot_data$Trajectory <- integrated_kinase_table$Trajectory[
  match(
    plot_data$kinase,
    integrated_kinase_table$kinase
  )
]

stopifnot(
  nrow(plot_data) ==
    4 * nrow(integrated_kinase_table),

  !anyNA(plot_data$logFC),

  !anyNA(plot_data$Trajectory)
)

cat(
  "Plot data rows:",
  nrow(plot_data),
  "\n"
)

p_kinase <- ggplot(
  plot_data,
  aes(
    x = stage,
    y = logFC,
    group = kinase
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey50"
  ) +

  geom_line(
    linewidth = 0.8,
    colour = "grey55"
  ) +

  geom_point(
    size = 2.8,
    colour = "black"
  ) +

  facet_wrap(
    ~ kinase,
    ncol = 3
  ) +

  labs(
    title =
      "Kinase expression across initial response, memory and interference",

    subtitle =
      "Nine significant memory-associated kinase candidates",

    x = NULL,

    y = "log2 fold change"
  ) +

  theme_bw(
    base_size = 11
  ) +

  theme(

    strip.text = element_text(
      face = "bold",
      size = 11
    ),

    plot.title = element_text(
      face = "bold",
      size = 14
    ),

    plot.subtitle = element_text(
      size = 10
    ),

    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      size = 8
    ),

    panel.grid.minor = element_blank(),

    legend.position = "none"
  )

ggsave(
  filename =
    "figures/Final_Kinase_Memory_Interference.png",

  plot = p_kinase,

  width = 14,

  height = 8,

  dpi = 600,

  bg = "white"
)


ggsave(
  filename =
    "figures/Final_Kinase_Memory_Interference.pdf",

  plot = p_kinase,

  width = 14,

  height = 8
)


png_exists <- file.exists(
  "figures/Final_Kinase_Memory_Interference.png"
)

pdf_exists <- file.exists(
  "figures/Final_Kinase_Memory_Interference.pdf"
)

cat(
  "PNG created:",
  png_exists,
  "\n"
)

cat(
  "PDF created:",
  pdf_exists,
  "\n"
)

stopifnot(
  png_exists,
  pdf_exists
)

cat(
  "Final kinase figure generation complete.\n"
)
