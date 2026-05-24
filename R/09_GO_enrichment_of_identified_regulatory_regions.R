# 9. GO enrichment of identified regulatory regions ----

# 9.1. Define background and foreground CRMs ----

crm_go_background <- crm_priority %>%
  dplyr::slice_head(n = min(1000, nrow(crm_priority)))

crm_go_foreground <- crm_go_background %>%
  dplyr::filter(n_patients >= 5)

if (nrow(crm_go_foreground) < 30) {
  crm_go_foreground <- crm_go_background %>%
    dplyr::filter(n_patients >= 4)
}

cat("GO background CRMs:", nrow(crm_go_background), "\n")
cat("GO foreground CRMs:", nrow(crm_go_foreground), "\n")


# 9.2. Map CRMs to associated genes ----

get_genes_for_crm <- function(crm_id) {
  
  res <- tryCatch(
    RBioGateway::crm2gene(crm_id),
    error = function(e) NULL
  )
  
  if (!is.data.frame(res) || nrow(res) == 0 || !"gene_name" %in% colnames(res)) {
    return(tibble::tibble())
  }
  
  res %>%
    dplyr::mutate(crm_name = crm_id) %>%
    dplyr::distinct(crm_name, gene_name, .keep_all = TRUE)
}

crm_gene_background <- dplyr::bind_rows(
  lapply(unique(crm_go_background$crm_name), get_genes_for_crm)
)

if (nrow(crm_gene_background) == 0) {
  stop("No CRM-associated genes were retrieved.")
}

crm_gene_foreground <- crm_gene_background %>%
  dplyr::filter(crm_name %in% crm_go_foreground$crm_name)

genes_background <- unique(crm_gene_background$gene_name)
genes_foreground <- unique(crm_gene_foreground$gene_name)

cat("Background genes:", length(genes_background), "\n")
cat("Foreground genes:", length(genes_foreground), "\n")


# 9.3. Map genes to proteins ----

get_proteins_for_gene <- function(gene_symbol) {
  
  res <- tryCatch(
    RBioGateway::gene2prot(gene_symbol, taxon = "Homo sapiens"),
    error = function(e) NULL
  )
  
  if (is.null(res) || is.data.frame(res) || length(res) == 0) {
    return(tibble::tibble())
  }
  
  if (length(res) == 1 &&
      grepl("No data available|Incorrect", res, ignore.case = TRUE)) {
    return(tibble::tibble())
  }
  
  tibble::tibble(
    gene_name = gene_symbol,
    protein_name = unique(as.character(res))
  )
}

gene_protein_background <- dplyr::bind_rows(
  lapply(genes_background, get_proteins_for_gene)
) %>%
  dplyr::distinct(gene_name, protein_name)

if (nrow(gene_protein_background) == 0) {
  stop("No proteins were retrieved.")
}

gene_protein_foreground <- gene_protein_background %>%
  dplyr::filter(gene_name %in% genes_foreground)

proteins_background <- unique(gene_protein_background$protein_name)
proteins_foreground <- unique(gene_protein_foreground$protein_name)

cat("Background proteins:", length(proteins_background), "\n")
cat("Foreground proteins:", length(proteins_foreground), "\n")


# 9.4. Retrieve GO annotations for proteins ----

format_go <- function(res, protein_id, ontology) {
  
  if (!is.data.frame(res) || nrow(res) == 0) {
    return(tibble::tibble())
  }
  
  id_col <- intersect(
    switch(
      ontology,
      BP = c("bp_id", "go_id"),
      CC = c("cc_id", "go_id"),
      MF = c("mf_id", "go_id")
    ),
    colnames(res)
  )[1]
  
  term_col <- intersect(
    switch(
      ontology,
      BP = c("bp", "bp_label", "go_term"),
      CC = c("cc", "cc_label", "go_term"),
      MF = c("mf", "mf_label", "go_term")
    ),
    colnames(res)
  )[1]
  
  if (is.na(id_col) || is.na(term_col)) {
    return(tibble::tibble())
  }
  
  tibble::tibble(
    protein_name = protein_id,
    ontology = ontology,
    go_id = as.character(res[[id_col]]),
    go_term = as.character(res[[term_col]])
  ) %>%
    dplyr::filter(!is.na(go_id), !is.na(go_term)) %>%
    dplyr::distinct()
}

query_go <- function(protein_id, ontology) {
  
  fun <- switch(
    ontology,
    BP = RBioGateway::prot2bp,
    CC = RBioGateway::prot2cc,
    MF = RBioGateway::prot2mf
  )
  
  res <- tryCatch(
    suppressWarnings(fun(protein_id)),
    error = function(e) e
  )
  
  if (inherits(res, "error")) {
    data_out <- tibble::tibble()
    status <- "error"
  } else if (is.character(res)) {
    data_out <- tibble::tibble()
    status <- "no_data"
  } else {
    data_out <- format_go(res, protein_id, ontology)
    status <- ifelse(nrow(data_out) > 0, "ok", "no_data")
  }
  
  list(
    data = data_out,
    log = tibble::tibble(
      protein_name = protein_id,
      ontology = ontology,
      status = status
    )
  )
}

run_go_queries <- function(query_grid, pause = 0.15) {
  
  if (nrow(query_grid) == 0) {
    return(list(data = tibble::tibble(), log = tibble::tibble()))
  }
  
  out <- lapply(seq_len(nrow(query_grid)), function(i) {
    
    if (i %% 100 == 0) {
      cat("GO queries:", i, "/", nrow(query_grid), "\n")
    }
    
    Sys.sleep(pause)
    query_go(query_grid$protein_name[i], query_grid$ontology[i])
  })
  
  list(
    data = dplyr::bind_rows(lapply(out, `[[`, "data")),
    log = dplyr::bind_rows(lapply(out, `[[`, "log"))
  )
}

go_grid <- expand.grid(
  protein_name = proteins_background,
  ontology = c("BP", "CC", "MF"),
  stringsAsFactors = FALSE
)

go_initial <- run_go_queries(go_grid, pause = 0.15)

retry_grid <- go_initial$log %>%
  dplyr::filter(status == "error") %>%
  dplyr::select(protein_name, ontology) %>%
  dplyr::distinct()

cat("GO queries retried:", nrow(retry_grid), "\n")

go_retry <- run_go_queries(retry_grid, pause = 1)

go_query_log_final <- dplyr::bind_rows(
  go_initial$log %>%
    dplyr::anti_join(retry_grid, by = c("protein_name", "ontology")),
  go_retry$log
)

go_annotations_background <- dplyr::bind_rows(
  go_initial$data,
  go_retry$data
) %>%
  dplyr::distinct(protein_name, ontology, go_id, go_term)

if (nrow(go_annotations_background) == 0) {
  stop("No GO annotations were retrieved.")
}

go_annotations_foreground <- go_annotations_background %>%
  dplyr::filter(protein_name %in% proteins_foreground)

cat("Background GO annotations:", nrow(go_annotations_background), "\n")
cat("Foreground GO annotations:", nrow(go_annotations_foreground), "\n")

print(go_query_log_final %>% dplyr::count(ontology, status))


# 9.5. Perform GO enrichment analysis ----

enrich_go_ontology <- function(go_annotations, foreground_proteins, ontology_name) {
  
  ann <- go_annotations %>%
    dplyr::filter(ontology == ontology_name) %>%
    dplyr::distinct(protein_name, go_id, go_term)
  
  universe <- unique(ann$protein_name)
  foreground <- intersect(unique(foreground_proteins), universe)
  
  N <- length(universe)
  n <- length(foreground)
  
  if (N == 0 || n == 0) {
    return(tibble::tibble())
  }
  
  ann %>%
    dplyr::group_by(go_id, go_term) %>%
    dplyr::summarise(
      K = dplyr::n_distinct(protein_name),
      k = dplyr::n_distinct(protein_name[protein_name %in% foreground]),
      foreground_proteins = paste(sort(unique(protein_name[protein_name %in% foreground])), collapse = "; "),
      background_proteins = paste(sort(unique(protein_name)), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::filter(k > 0) %>%
    dplyr::mutate(
      ontology = ontology_name,
      N = N,
      n = n,
      pvalue = stats::phyper(k - 1, K, N - K, n, lower.tail = FALSE),
      FDR = stats::p.adjust(pvalue, method = "BH"),
      fold_enrichment = (k / n) / (K / N)
    )
}

go_enrichment_results <- dplyr::bind_rows(
  enrich_go_ontology(go_annotations_background, proteins_foreground, "BP"),
  enrich_go_ontology(go_annotations_background, proteins_foreground, "CC"),
  enrich_go_ontology(go_annotations_background, proteins_foreground, "MF")
) %>%
  dplyr::mutate(
    FDR_global = stats::p.adjust(pvalue, method = "BH")
  ) %>%
  dplyr::arrange(FDR, pvalue, dplyr::desc(fold_enrichment))

significant_go <- go_enrichment_results %>%
  dplyr::filter(FDR < 0.05)

significant_go_global <- go_enrichment_results %>%
  dplyr::filter(FDR_global < 0.05)

cat("GO terms tested:", nrow(go_enrichment_results), "\n")
cat("Significant GO terms, FDR by ontology < 0.05:", nrow(significant_go), "\n")
cat("Significant GO terms, global FDR < 0.05:", nrow(significant_go_global), "\n")

print(significant_go, n = Inf)


# 9.6. Summarise and visualise GO enrichment results ----

go_summary <- go_enrichment_results %>%
  dplyr::group_by(ontology) %>%
  dplyr::summarise(
    n_terms = dplyr::n(),
    n_significant_FDR_0_05 = sum(FDR < 0.05),
    n_significant_global_FDR_0_05 = sum(FDR_global < 0.05),
    best_term = go_term[which.min(FDR)],
    best_FDR = min(FDR),
    best_global_FDR = min(FDR_global),
    .groups = "drop"
  )

print(go_summary)

top_go_plot <- go_enrichment_results %>%
  dplyr::mutate(
    FDR_plot = pmax(FDR, .Machine$double.xmin),
    label = paste0(ontology, ": ", go_term)
  ) %>%
  dplyr::slice_head(n = 20)

plot_go <- ggplot2::ggplot(
  top_go_plot,
  ggplot2::aes(
    x = stats::reorder(label, -log10(FDR_plot)),
    y = -log10(FDR_plot)
  )
) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "GO enrichment of recurrently affected CRMs",
    x = "GO term",
    y = "-log10(FDR by ontology)"
  ) +
  ggplot2::theme_classic(base_size = 13)

print(plot_go)
