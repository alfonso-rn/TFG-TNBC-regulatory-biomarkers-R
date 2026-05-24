# 2. Download gene expression data from primary TCGA-BRCA tumors ----
#
# This script downloads RNA-seq gene expression data from the TCGA-BRCA project.
#
# In this analysis, only primary tumor samples are selected. Therefore, downstream
# comparisons will be performed as TNBC vs non-TNBC rather than Tumor vs Normal.

# 2.1. Query TCGA-BRCA RNA-seq expression data ----

BRCAproy <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type= "STAR - Counts",
  experimental.strategy = "RNA-Seq",
  access = "open",
  sample.type = "Primary Tumor"
)


# 2.2. Download & Prepare RNA-seq expression data ----

# The API method is used because this query retrieves many files.
GDCdownload(BRCAproy, method = "api",
            files.per.chunk = 30,
)

BRCAprimarytumor <- GDCprepare(BRCAproy)
