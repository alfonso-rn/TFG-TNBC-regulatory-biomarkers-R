# 1. Identification of TNBC patients within TCGA-BRCA ----
#
# This script downloads clinical data from the TCGA-BRCA project and identifies
# patients with clinically defined triple-negative breast cancer (TNBC).
#
# In this analysis, TNBC patients are defined using clinical immunohistochemistry
# (IHC) information. Specifically, patients are classified as TNBC when the three
# following biomarkers are reported as negative:
#
#   - Estrogen receptor (ER)
#   - Progesterone receptor (PR)
#   - Human epidermal growth factor receptor 2 (HER2)
#
# It is important to note that TNBC is sometimes approximated using the PAM50
# molecular classification, especially by selecting Basal-like tumors. However,
# Basal-like breast cancer and clinically defined TNBC are not exactly equivalent.
# For this reason, this script identifies TNBC patients based on clinical receptor
# status rather than PAM50 subtype.
#
# If PAM50 Basal-like classification is preferred, it could be obtained using:
#
#   subtypes <- TCGAquery_subtype(tumor = "BRCA")
#   basal_barcodes <- unique(
#     subtypes$patient[subtypes$BRCA_Subtype_PAM50 == "Basal"]
#   )

# 1.1. Query and download TCGA-BRCA clinical data ----

# Query clinical data from the TCGA-BRCA project.
# The selected format is "BCR Biotab", which contains structured clinical tables.
brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

# Download the queried clinical data from the GDC.
GDCdownload(brca)

# Prepare the downloaded data into an R object.
clinicalBCR <- GDCprepare(brca)

# 1.2. Extract patient-level clinical information ----

# Display the available clinical tables contained in the downloaded object.
names(clinicalBCR) 

# For TCGA-BRCA, the patient-level clinical information is stored in:
#
#   clinical_patient_brca
#
# Extract the patient-level clinical table.
clinical_patient <- clinicalBCR[["clinical_patient_brca"]][-c(1, 2), ]
# In the BCR Biotab format, the first two rows of this table are not patient records.

# Display the available clinical variables.
names(clinical_patient)

# 1.3. Identify clinically defined TNBC patients using IHC status ----

# Columns used to identify clinically triple-negative patients:
#
#   er_status_by_ihc   -> Estrogen receptor status
#   pr_status_by_ihc   -> Progesterone receptor status
#   her2_status_by_ihc -> HER2 receptor status
#
# The following compact summary table shows the number of samples assigned to
# each possible category for the three IHC biomarkers.
receptor_status_summary <- sapply(
  clinical_patient[, c(
    "er_status_by_ihc",
    "pr_status_by_ihc",
    "her2_status_by_ihc"
  )],
  table,
  useNA = "ifany"
)

receptor_status_summary

# Classify patients as TNBC when ER, PR, and HER2 are all "Negative".
tnbc <- clinical_patient %>%
  dplyr::filter(
    er_status_by_ihc == "Negative",
    pr_status_by_ihc == "Negative",
    her2_status_by_ihc == "Negative"
  )

cat("Number of clinically defined TNBC patients identified:", nrow(tnbc), "\n")

# Extract TCGA barcodes for the identified TNBC patients.
tnbc_barcodes <- tnbc %>%
  dplyr::pull(bcr_patient_barcode)
