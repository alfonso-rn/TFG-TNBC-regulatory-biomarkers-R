# 5. Survival analysis within TNBC ----
#
# This script evaluates the prognostic value of TNBC-associated DEGs.
#
# Survival analysis is performed only within clinically defined TNBC patients,
# using overall survival and gene expression data.

# 5.1. Prepare clinical survival data ----

clin <- GDCquery_clinic("TCGA-BRCA", "clinical")

clin_surv <- clin %>%
  dplyr::transmute(
    patient = bcr_patient_barcode,
    OS_time = dplyr::coalesce(
      as.numeric(days_to_death),
      as.numeric(days_to_last_follow_up)
    ),
    OS_event = as.numeric(vital_status == "Dead")
  ) %>%
  dplyr::filter(
    patient %in% tnbc_barcodes,
    !is.na(OS_time),
    !is.na(OS_event),
    OS_time > 0
  )

# 5.2. Match TNBC expression data with clinical data ----

expr_tnbc <- dataFilt[, samplesTNBC]

# Convert sample barcodes to patient barcodes
colnames(expr_tnbc) <- substr(colnames(expr_tnbc), 1, 12)

common_patients <- intersect(
  colnames(expr_tnbc),
  clin_surv$patient
)

expr_tnbc <- expr_tnbc[, common_patients]

clin_tnbc <- clin_surv[
  match(common_patients, clin_surv$patient),
]

cat("TNBC patients with expression and survival data:", nrow(clin_tnbc), "\n")

# 5.3. Univariate Cox analysis for each DEG ----

genes_for_survival <- intersect(
  rownames(degs_tbl),
  rownames(expr_tnbc)
)

cox_gene <- function(gene) {
  
  expr <- log2(as.numeric(expr_tnbc[gene, ]) + 1)
  
  # Exclude genes with low variability
  if (length(unique(expr[!is.na(expr)])) < 5) return(NULL)
  if (sd(expr, na.rm = TRUE) == 0) return(NULL)
  
  fit <- try(
    suppressWarnings(
      survival::coxph(
        survival::Surv(clin_tnbc$OS_time, clin_tnbc$OS_event) ~ scale(expr),
        control = survival::coxph.control(iter.max = 50)
      )
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) return(NULL)
  
  fit_summary <- summary(fit)
  
  gene_symbol <- if ("gene_name" %in% colnames(degs_tbl)) {
    as.character(degs_tbl[gene, "gene_name"])
  } else {
    gene
  }
  
  res <- data.frame(
    gene_id = gene,
    gene_name = gene_symbol,
    HR = fit_summary$coefficients[1, "exp(coef)"],
    lower95 = fit_summary$conf.int[1, "lower .95"],
    upper95 = fit_summary$conf.int[1, "upper .95"],
    pvalue = fit_summary$coefficients[1, "Pr(>|z|)"]
  )
  
  # Remove unstable Cox estimates
  if (!is.finite(res$HR) || res$HR > 100 || res$HR < 0.01) return(NULL)
  
  return(res)
}

cox_results <- dplyr::bind_rows(
  lapply(genes_for_survival, cox_gene)
) %>%
  dplyr::mutate(
    FDR = p.adjust(pvalue, method = "BH")
  ) %>%
  dplyr::arrange(FDR)

head(cox_results, 20)

cat("Genes tested in Cox analysis:", nrow(cox_results), "\n")

# 5.4. Kaplan-Meier plots for the top 4 prognostic genes ----

top4 <- cox_results %>%
  dplyr::slice_head(n = 4)

plot_km_gene <- function(gene_id, gene_name, FDR) {
  
  df_km <- data.frame(
    OS_time = clin_tnbc$OS_time,
    OS_event = clin_tnbc$OS_event,
    expr = log2(as.numeric(expr_tnbc[gene_id, ]) + 1)
  )
  
  df_km$group <- ifelse(
    df_km$expr > median(df_km$expr, na.rm = TRUE),
    "High expression",
    "Low expression"
  )
  
  fit_km <- survival::survfit(
    survival::Surv(OS_time, OS_event) ~ group,
    data = df_km
  )
  
  survminer::ggsurvplot(
    fit_km,
    data = df_km,
    pval = TRUE,
    risk.table = FALSE,
    title = paste0(gene_name, " (FDR = ", signif(FDR, 3), ")"),
  )
}

km_plots <- lapply(seq_len(nrow(top4)), function(i) {
  plot_km_gene(
    gene_id = top4$gene_id[i],
    gene_name = top4$gene_name[i],
    FDR = top4$FDR[i]
  )
})

survminer::arrange_ggsurvplots(
  km_plots,
  print = TRUE,
  ncol = 2,
  nrow = 2
)
