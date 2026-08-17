library(ggplot2)
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
dim(plot_data)

head(plot_data, 15)

table(plot_data$stage)

table(plot_data$kinase)
# -----------------------------
# 5. Add trajectory information
# -----------------------------

plot_data$Trajectory <- integrated_kinase_table$Trajectory[
  match(
    plot_data$kinase,
    integrated_kinase_table$kinase
  )
]


# -----------------------------
# 6. Create final trajectory plot
# -----------------------------

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
    title = "Kinase expression across initial response, memory and interference",
    subtitle = "Nine significant memory-associated kinase candidates",
    x = NULL,
    y = "log2 fold change"
  ) +
  
  theme_bw(base_size = 11) +
  
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
p_kinase
plot_data$stage <- factor(
  plot_data$stage,
  levels = c(
    "Initial Dex",
    "Memory",
    "Memory + JNK-IN-8",
    "Memory + T-5224"
  ),
  labels = c(
    "Initial\nDex",
    "Memory",
    "+ JNK-IN-8",
    "+ T-5224"
  )
)
labs(
  title = "Kinase effect-size profiles across Dex memory and interference",
  subtitle = "Nine FDR-significant kinases in the post-Dex memory contrast",
  x = NULL,
  y = "log2 fold change"
)

ggsave(
  filename = "figures/Final_Kinase_Memory_Interference.png",
  plot = p_kinase,
  width = 14,
  height = 8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "figures/Final_Kinase_Memory_Interference.pdf",
  plot = p_kinase,
  width = 14,
  height = 8
)

file.exists(
  "figures/Final_Kinase_Memory_Interference.png"
)