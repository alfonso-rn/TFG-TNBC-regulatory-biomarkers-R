# 7. Identification of CRMs affected by non-truncating mutations in DEGs ----
#
# This script identifies cis-regulatory modules (CRMs) that overlap with
# non-truncating mutations located in differentially expressed genes (DEGs).
#
# RBioGateway is used to query CRMs overlapping each mutation and to retrieve
# genomic information for the prioritized CRMs.

# 7.1. Identify CRMs overlapping candidate mutations ----

get_crm_for_mutation <- function(i) {
  
  mut <- mutations_for_crm[i, ]
  
  res <- RBioGateway::getCRMs_by_overlap(
    chromosome = mut$chr_biogateway,
    start = mut$Start_Position,
    end = mut$End_Position
  )
  
  if (!is.data.frame(res) || nrow(res) == 0) {
    return(NULL)
  }
  
  res %>%
    dplyr::mutate(
      mutation_id = mut$mutation_id,
      Hugo_Symbol = mut$Hugo_Symbol,
      patient = mut$patient,
      Tumor_Sample_Barcode = mut$Tumor_Sample_Barcode,
      Chromosome = mut$Chromosome,
      chr_biogateway = mut$chr_biogateway,
      Start_Position_mut = mut$Start_Position,
      End_Position_mut = mut$End_Position,
      Variant_Classification = mut$Variant_Classification,
      Variant_Type = mut$Variant_Type,
      Reference_Allele = mut$Reference_Allele,
      Tumor_Seq_Allele2 = mut$Tumor_Seq_Allele2
    )
}

crm_overlap_results <- dplyr::bind_rows(
  lapply(seq_len(nrow(mutations_for_crm)), get_crm_for_mutation)
)

unique_crm_mutations <- crm_overlap_results %>%
  dplyr::distinct(
    mutation_id,
    Tumor_Sample_Barcode,
    Hugo_Symbol,
    Chromosome,
    Start_Position_mut,
    End_Position_mut,
    Reference_Allele,
    Tumor_Seq_Allele2
  )

# Results summary table
crm_summary_table <- tibble::tibble(
  Metric = c(
    "Candidate non-truncating mutations in DEGs",
    "CRM-mutation overlaps",
    "Unique mutations overlapping CRMs",
    "Tumor samples with CRM-overlapping mutations",
    "Unique affected CRMs",
    "Genes with non-truncating mutations in CRMs"
  ),
  Value = c(
    nrow(mutations_for_crm),
    nrow(crm_overlap_results),
    nrow(unique_crm_mutations),
    length(unique(crm_overlap_results$Tumor_Sample_Barcode)),
    length(unique(crm_overlap_results$crm_name)),
    length(unique(crm_overlap_results$Hugo_Symbol))
  )
)

crm_summary_table

# 7.2. Prioritize affected CRMs ----

crm_priority <- crm_overlap_results %>%
  dplyr::group_by(crm_name, start, end) %>%
  dplyr::summarise(
    n_mutations = dplyr::n_distinct(mutation_id),
    n_patients = dplyr::n_distinct(patient),
    n_genes = dplyr::n_distinct(Hugo_Symbol),
    associated_genes = paste(sort(unique(Hugo_Symbol)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_patients),
    dplyr::desc(n_mutations),
    dplyr::desc(n_genes)
  )

cat("Prioritized affected CRMs:", nrow(crm_priority), "\n")

print(n = 30, crm_priority)

# 7.3. Recurrent affected CRMs plot ----

top_n_crms <- 30

crm_plot_data <- crm_priority %>%
  arrange(desc(n_patients), desc(n_mutations)) %>%
  slice_head(n = top_n_crms) %>%
  mutate(
    main_gene = ifelse(
      grepl(";", associated_genes),
      paste0(sub(";.*", "", associated_genes), " + others"),
      associated_genes
    ),
    crm_id = sub("^crm/", "", crm_name),
    crm_id = factor(crm_id, levels = rev(crm_id))
  )

ggplot(crm_plot_data, aes(x = n_patients, y = crm_id)) +
  geom_segment(
    aes(x = 0, xend = n_patients, yend = crm_id, color = main_gene),
    linewidth = 0.8,
    alpha = 0.8
  ) +
  geom_point(aes(size = n_mutations, color = main_gene), alpha = 0.95) +
  labs(
    title = "Recurrent CRMs affected by non-truncating mutations in TNBC",
    subtitle = paste("Top", top_n_crms, "CRMs ranked by number of affected patients"),
    x = "Number of TNBC patients with overlapping mutations",
    y = "Affected CRM",
    color = "Associated gene",
    size = "Number of mutations"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

# 7.4. Retrieve genomic information for prioritized CRMs ----

n_crm_info <- min(100, nrow(crm_priority))

crm_to_info <- crm_priority %>%
  dplyr::slice_head(n = n_crm_info)

get_crm_info_simple <- function(crm_id) {
  
  res <- RBioGateway::getCRM_info(crm_id)
  
  if (!is.data.frame(res) || nrow(res) == 0) {
    return(NULL)
  }
  
  res %>%
    dplyr::mutate(crm_name = crm_id)
}

crm_info_results <- dplyr::bind_rows(
  lapply(crm_to_info$crm_name, get_crm_info_simple)
)

crm_info_results <- crm_info_results %>%
  dplyr::left_join(
    crm_priority,
    by = "crm_name"
  )

cat("CRMs annotated with getCRM_info():", nrow(crm_info_results), "\n")