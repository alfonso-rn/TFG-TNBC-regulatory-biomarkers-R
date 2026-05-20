# 0. Environment configuration ----

# 0.1. Group required packages ----

# CRAN packages ()
cran_pkgs <- c(
  "dplyr", "survival", "tibble", "survminer", "ggplot2"
)

# Bioconductor packages
bioc_pkgs <- c(
  "TCGAbiolinks", "TCGAbiolinksGUI.data", "SummarizedExperiment",
  "edgeR", "AnnotationDbi", "org.Hs.eg.db", "EDASeq", "maftools"
)

# CRAN + Bioconductor + RBioGateway
all_pkgs <- c(
  bioc_pkgs, cran_pkgs, "RBioGateway"
)

# 0.2. Installation of packages ----

# If CRAN is not installed: https://cran.r-project.org
install.packages(cran_pkgs)

if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")
BiocManager::install(bioc_pkgs)

install.packages("devtools")
devtools::install_github("tecnomod-um/RBioGateway")

# 0.3. Load all packages ----

lapply(all_pkgs, library, character.only = TRUE)

message("Environment successfully configured")
