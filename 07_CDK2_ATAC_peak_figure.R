library(ggplot2)

infile <- "results/ATAC_integration/CDK2_significant_peak_trajectory.tsv"

out_png <- "figures/CDK2_significant_ATAC_peak_trajectory.png"
out_pdf <- "figures/CDK2_significant_ATAC_peak_trajectory.pdf"

dir.create("figures", showWarnings = FALSE)

x <- read.delim(infile, check.names = FALSE)

condition_order <- c(
  "3d_DMSO",
  "3d_Dex",
  "Dex_to_DMSO",
  "DMSO_to_Tram",
  "Dex_to_Tram"
)

x$condition <- factor(x$condition, levels = condition_order)

replicates <- rbind(
  data.frame(
    condition = x$condition,
    replicate = "A",
    ATAC = x$replicate_A
  ),
  data.frame(
    condition = x$condition,
    replicate = "B",
    ATAC = x$replicate_B
  )
)

means <- data.frame(
  condition = x$condition,
  mean_ATAC = x$mean_ATAC
)

p <- ggplot() +
  geom_line(
    data = means,
    aes(x = condition, y = mean_ATAC, group = 1),
    linewidth = 0.8
  ) +
  geom_point(
    data = means,
    aes(x = condition, y = mean_ATAC),
    size = 3
  ) +
  geom_point(
    data = replicates,
    aes(x = condition, y = ATAC),
    position = position_jitter(width = 0.07, height = 0),
    size = 2.5,
    shape = 1,
    stroke = 1
  ) +
  annotate(
    "text",
    x = 4.5,
    y = max(replicates$ATAC) * 1.12,
    label = "Dex->Tram vs Tram-only:\nDESeq2 log2FC = +1.028, FDR = 0.00349",
    size = 3.6
  ) +
  labs(
    title = "CDK2 memory-associated ATAC peak",
    subtitle = "chr12:56358987–56360175 (hg19)",
    x = NULL,
    y = "ATAC accessibility (bigWig mean0)",
    caption = "Open circles: biological replicates; filled points: condition means."
  ) +
  scale_x_discrete(
    labels = c(
      "3d DMSO",
      "3d Dex",
      "Dex->DMSO",
      "Tram-only",
      "Dex->Tram"
    )
  ) +
  coord_cartesian(
    ylim = c(0, max(replicates$ATAC) * 1.28)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

ggsave(out_png, p, width = 7.5, height = 5.5, dpi = 300)
ggsave(out_pdf, p, width = 7.5, height = 5.5)

cat("PNG created:", file.exists(out_png), "\n")
cat("PDF created:", file.exists(out_pdf), "\n")
cat("Mean trajectory:\n")
print(means)
