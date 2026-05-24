# 4. Differential expression analysis: TNBC vs non-TNBC ----
#
# Differential expression analysis is performed between clinically defined TNBC
# and non-TNBC primary tumor samples.
#
# First, a global DEA is performed using permissive thresholds in order to retain
# all tested features for visualization in the volcano plot.
#
# Since TNBC is provided as the second condition, positive logFC values are
# interpreted as higher expression in TNBC compared with non-TNBC.
#
# DEGs are then defined using FDR < 0.01 and |log2FC| >= 1, corresponding to
# statistically significant genes with at least a two-fold expression change.

# 4.1. Global DEA keeping all analyzed features ----

dataDEA_all <- TCGAanalyze_DEA(
  mat1 = dataFilt[, samplesNonTNBC],
  mat2 = dataFilt[, samplesTNBC],
  Cond1type = "Non_TNBC",
  Cond2type = "TNBC",
  fdr.cut = 1,
  logFC.cut = 0,
  method = "glmLRT",
  pipeline = "edgeR"
)

# 4.2. Prepare DEA table for visualization ----

dea_all_tbl <- as.data.frame(dataDEA_all)
dea_all_tbl$feature_id <- rownames(dea_all_tbl)

# Avoid infinite values in the volcano plot when FDR = 0
dea_all_tbl$minus_log10_FDR <- -log10(
  pmax(dea_all_tbl$FDR, .Machine$double.xmin)
)

# Classify features according to DEA thresholds
dea_all_tbl$group <- "Not significant"

dea_all_tbl$group[
  dea_all_tbl$FDR < 0.01 & dea_all_tbl$logFC >= 1
] <- "Upregulated in TNBC"

dea_all_tbl$group[
  dea_all_tbl$FDR < 0.01 & dea_all_tbl$logFC <= -1
] <- "Downregulated in TNBC"

# 4.3. Volcano plot ----

volcano_plot <- ggplot(
  dea_all_tbl,
  aes(x = logFC, y = minus_log10_FDR, color = group)
) +
  geom_point(alpha = 0.6, size = 1) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(0.01),
    linetype = "dashed",
    color = "grey40"
  ) +
  scale_color_manual(
    values = c(
      "Downregulated in TNBC" = "#2C7BB6",
      "Not significant" = "grey75",
      "Upregulated in TNBC" = "#D7191C"
    )
  ) +
  labs(
    title = "Differential expression analysis: TNBC vs non-TNBC",
    subtitle = "TCGA-BRCA primary tumor samples",
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = ""
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

print(volcano_plot)

# 4.4. Filter DEGs from the global DEA ----

degs_tbl <- dea_all_tbl[
  dea_all_tbl$FDR < 0.01 & abs(dea_all_tbl$logFC) >= 1,
]

rownames(degs_tbl) <- degs_tbl$feature_id

dataDEGs <- degs_tbl

cat("Upregulated DEGs in TNBC:", sum(degs_tbl$logFC > 0), "\n")
cat("Downregulated DEGs in TNBC:", sum(degs_tbl$logFC < 0), "\n")
cat("Total DEGs:", nrow(degs_tbl), "\n")

deg_genes <- if ("gene_name" %in% colnames(degs_tbl)) {
  unique(degs_tbl$gene_name)
} else {
  unique(rownames(degs_tbl))
}

cat("Analyzed features:", nrow(dea_all_tbl), "\n")
cat("DEGs TNBC vs non-TNBC:", length(deg_genes), "\n")

# 4.5. Percentage plot of gene types by group ----
gene_type_summary <- degs_tbl %>%
  count(group, gene_type)

ggplot(gene_type_summary, aes(x = group, y = n, fill = gene_type)) +
  geom_col(position = "fill") +
  scale_y_continuous(
    labels = function(x) paste0(round(x * 100), "%")
  ) +
  labs(
    title = "Distribution of gene types among DEGs",
    subtitle = "TNBC vs non-TNBC",
    x = "",
    y = "Percentage of DEGs",
    fill = "Gene type"
  ) +
  theme_classic()
