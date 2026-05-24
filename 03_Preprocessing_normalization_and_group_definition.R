# 3. Preprocessing, normalization and group definition: TNBC vs non-TNBC ----

# 3.1. Preprocess RNA-seq expression data ----

dataProc <- TCGAanalyze_Preprocessing(
  object = BRCAprimarytumor,
  cor.cut = 0.6,
  filename = "Analisis_Calidad_BRCA_TNBC.png"
)

# 3.2. Normalize gene expression data ----

dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataProc,
  geneInfo = geneInfoHT,
  method = "gcContent"
)

# 3.3. Filter lowly expressed genes ----

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile",
  qnt.cut = 0.25
)

# 3.4. Selection of sample groups ----

# Keep one sample per patient
sample_patients <- substr(colnames(dataFilt), 1, 12)

dataFilt <- dataFilt[, !duplicated(sample_patients)]

sample_patients <- substr(colnames(dataFilt), 1, 12)

# Define sample groups: TNBC vs non-TNBC
samplesTNBC <- colnames(dataFilt)[
  sample_patients %in% tnbc_barcodes
]

samplesNonTNBC <- colnames(dataFilt)[
  !(sample_patients %in% tnbc_barcodes)
]

cat("TNBC samples:", length(samplesTNBC), "\n")
cat("Non-TNBC samples:", length(samplesNonTNBC), "\n")