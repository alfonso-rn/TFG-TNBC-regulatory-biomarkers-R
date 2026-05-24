# 6. Somatic mutations and non-truncating variants in TNBC ----
#
# This script downloads somatic mutation data from TCGA-BRCA and keeps mutations
# detected in clinically defined TNBC patients.
#
# Truncating variants are removed in order to prioritize non-truncating mutations
# in differentially expressed genes for downstream regulatory-region analyses.

# 6.1. Query and download TCGA-BRCA somatic mutation data ----

query_mut <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

GDCdownload(
  query_mut,
  method = "api",
  files.per.chunk = 30
)

maf <- GDCprepare(query_mut)

# 6.2. Keep mutations from TNBC patients ----

maf_tnbc <- maf %>%
  dplyr::mutate(
    patient = substr(Tumor_Sample_Barcode, 1, 12)
  ) %>%
  dplyr::filter(
    patient %in% tnbc_barcodes
  )

# 6.3. Remove truncating variants ----

truncating_classes <- c(
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "Translation_Start_Site",
  "Splice_Site"
)

maf_tnbc_nontrunc <- maf_tnbc %>%
  dplyr::filter(
    !(Variant_Classification %in% truncating_classes)
  )

cat("TNBC somatic mutations:", nrow(maf_tnbc), "\n")
cat("TNBC non-truncating mutations:", nrow(maf_tnbc_nontrunc), "\n")

# 6.4. Convert DEG identifiers to HGNC symbols ----

deg_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(degs_tbl),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
) %>%
  as.character() %>%
  unique() %>%
  stats::na.omit()

# 6.5. Identify DEGs with non-truncating somatic mutations ----

mutated_deg_nontrunc <- maf_tnbc_nontrunc %>%
  dplyr::filter(
    Hugo_Symbol %in% deg_symbols
  )

cat( "DEGs with non-truncating somatic mutations:",
     dplyr::n_distinct(mutated_deg_nontrunc$Hugo_Symbol), "\n")

# 6.6. Oncoplot of the most frequently mutated candidate genes ----

maf_oncoplot <- mutated_deg_nontrunc %>%
  dplyr::mutate(
    Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 12)
  ) %>%
  dplyr::filter(
    !is.na(Hugo_Symbol),
    !is.na(Tumor_Sample_Barcode),
    !is.na(Variant_Classification)
  ) %>%
  dplyr::distinct(
    Hugo_Symbol,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    Variant_Classification,
    Variant_Type,
    Tumor_Sample_Barcode,
    .keep_all = TRUE
  )

top_genes <- maf_oncoplot %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(
    n_patients = dplyr::n_distinct(Tumor_Sample_Barcode),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_patients)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(Hugo_Symbol)

maf_obj <- maftools::read.maf(
  maf = maf_oncoplot %>%
    dplyr::filter(Hugo_Symbol %in% top_genes),
  verbose = FALSE
)

maftools::oncoplot(
  maf = maf_obj,
  genes = top_genes,
  removeNonMutated = TRUE,
  sortByMutation = TRUE,
  showTumorSampleBarcodes = FALSE,
  fontSize = 0.8,
  titleText = "DEGs with non-truncating mutations in TNBC patients"
)

# 6.7. Prepare mutation coordinates for RBioGateway ----
  
  mutations_for_crm <- mutated_deg_nontrunc %>%
  dplyr::select(
    Hugo_Symbol,
    patient,
    Tumor_Sample_Barcode,
    Chromosome,
    Start_Position,
    End_Position,
    Variant_Classification,
    Variant_Type,
    Reference_Allele,
    Tumor_Seq_Allele2
  ) %>%
  dplyr::mutate(
    chr_biogateway = paste0("chr-", gsub("^chr[-]?", "", Chromosome)),
    Start_Position = as.integer(Start_Position),
    End_Position = as.integer(End_Position),
    mutation_id = paste(
      Tumor_Sample_Barcode,
      Hugo_Symbol,
      Chromosome,
      Start_Position,
      End_Position,
      Reference_Allele,
      Tumor_Seq_Allele2,
      sep = "_"
    )
  ) %>%
  dplyr::filter(
    !is.na(chr_biogateway),
    !is.na(Start_Position),
    !is.na(End_Position)
  )

cat("Candidate non-truncating mutations for CRM search:", nrow(mutations_for_crm), "\n")