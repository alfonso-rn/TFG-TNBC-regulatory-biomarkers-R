# 8. Identification of transcription factors associated with affected CRMs ----
#
# This script retrieves transcription factors associated with prioritized CRMs
# affected by non-truncating mutations in differentially expressed genes.
#
# RBioGateway is used to query transcription factors associated with each CRM.
# Mutation-CRM-TF relationships are then integrated and summarized to prioritize
# candidate transcription factors.

# 8.1. Query transcription factors associated with CRMs ----

get_tfac_for_crm <- function(crm_id) {
  
  res <- tryCatch(
    RBioGateway::crm2tfac(crm_id),
    error = function(e) NULL
  )
  
  if (!is.data.frame(res) || nrow(res) == 0) {
    return(NULL)
  }
  
  res %>%
    dplyr::mutate(crm_name = crm_id)
}

crm_tfac_results <- dplyr::bind_rows(
  lapply(crm_to_info$crm_name, get_tfac_for_crm)
)

if (nrow(crm_tfac_results) == 0) {
  crm_tfac_results <- data.frame(
    crm_name = character(),
    tfac_name = character(),
    database = character(),
    articles = character(),
    evidence = character()
  )
}

cat("CRM-TF associations retrieved:", nrow(crm_tfac_results), "\n")

cat("Unique transcription factors:", dplyr::n_distinct(crm_tfac_results$tfac_name), "\n")

# 8.2. Integrate mutations, CRMs and transcription factors ----

tf_candidates <- crm_overlap_results %>%
  dplyr::filter(crm_name %in% crm_to_info$crm_name) %>%
  dplyr::left_join(
    crm_tfac_results,
    by = "crm_name",
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(!is.na(tfac_name)) %>%
  dplyr::left_join(
    crm_priority %>%
      dplyr::select(crm_name, n_mutations, n_patients, n_genes),
    by = "crm_name"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_patients),
    dplyr::desc(n_mutations),
    crm_name,
    tfac_name
  )

cat("Mutation-CRM-TF rows:", nrow(tf_candidates), "\n")

cat( "Mutations associated with TFs:", dplyr::n_distinct(tf_candidates$mutation_id), "\n")

cat("CRMs associated with TFs:", dplyr::n_distinct(tf_candidates$crm_name), "\n")

cat("Unique TFs:", dplyr::n_distinct(tf_candidates$tfac_name), "\n")

# 8.3. Summary by transcription factor ----

tf_summary_by_tf <- tf_candidates %>%
  dplyr::group_by(tfac_name) %>%
  dplyr::summarise(
    n_crms = dplyr::n_distinct(crm_name),
    n_mutations = dplyr::n_distinct(mutation_id),
    n_patients = dplyr::n_distinct(patient),
    n_genes = dplyr::n_distinct(Hugo_Symbol),
    associated_genes = paste(sort(unique(Hugo_Symbol)), collapse = "; "),
    databases = paste(sort(unique(stats::na.omit(database))), collapse = "; "),
    articles = paste(sort(unique(stats::na.omit(articles))), collapse = "; "),
    evidence = paste(sort(unique(stats::na.omit(evidence))), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_patients),
    dplyr::desc(n_mutations),
    dplyr::desc(n_crms)
  )

head(tf_summary_by_tf, 20)

# 8.4. Visualization of the main candidate transcription factors ----

if (nrow(tf_summary_by_tf) > 0) {
  
  top_tf_plot <- tf_summary_by_tf %>%
    dplyr::slice_head(n = 20)
  
  plot_tf <- ggplot2::ggplot(
    top_tf_plot,
    ggplot2::aes(
      x = reorder(tfac_name, n_mutations),
      y = n_mutations
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = associated_genes),
      hjust = -0.05,
      size = 3
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.45))
    ) +
    ggplot2::labs(
      title = "Main transcription factors associated with affected CRMs",
      x = "Transcription factor",
      y = "Number of associated mutations"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.margin = ggplot2::margin(5.5, 160, 5.5, 5.5)
    )
  
  print(plot_tf)
}
