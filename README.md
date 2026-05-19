# Description
R workflow for analyzing TCGA-BRCA data to identify differentially expressed genes, survival-associated biomarkers and affected regulatory regions in triple-negative breast cancer using TCGAbiolinks and RBioGateway.

# Project overview
This repository contains the computational workflow developed for an undergraduate thesis focused on the integrative characterization of triple-negative breast cancer (TNBC) within the TCGA-BRCA cohort. The study uses TCGAbiolinks to retrieve and preprocess clinical, transcriptomic, and somatic mutation data from the Genomic Data Commons, and RBioGateway to contextualize candidate non-truncating mutations within cis-regulatory modules (CRMs) and their associated transcription factor binding sites (TFBS). The main objective is to identify differentially expressed genes with prognostic relevance in TNBC and to extend the analysis beyond conventional transcriptome-level comparisons by incorporating somatic variation and regulatory interpretation. 

# Workflow files
The repository is organized as a sequential and reproducible set of numbered scripts, where each file represents one major analytical stage of the study. 

## 00_Environment_configuration.R 
Initializes the computational environment, installs and loads the required packages, and records the software versions used in the analysis. This file ensures that the workflow can be reproduced under a controlled R/Bioconductor setup before any data are queried or processed. 

## 02_Download_of_primary_BRCA_expression_data.R

Queries and downloads harmonized RNA-seq gene expression data from primary breast tumor samples in the TCGA-BRCA project using TCGAbiolinks. This file establishes the transcriptomic dataset used for downstream comparisons between TNBC and non-TNBC tumors.

## 03_Preprocessing_normalization_and_group_definition.R

Performs preprocessing of the RNA-seq expression dataset, including sample filtering, quality control, normalization, expression filtering, and removal of duplicated or non-informative entries when required. This script also assigns each primary tumor sample to the TNBC or non-TNBC analytical group, creating the final expression matrix and phenotype information used in the differential expression analysis.

## 04_Differential_expression_analysis.R

Carries out the differential expression analysis between TNBC and non-TNBC primary breast tumors. This script identifies genes whose expression differs significantly between both groups, applies statistical thresholds to classify relevant genes, and generates graphical summaries such as volcano plots. The resulting differentially expressed genes constitute the first set of candidates for downstream prognostic and regulatory analyses.

## 05_Survival_analysis_in_TNBC.R

Integrates gene expression data with clinical follow-up information from TNBC patients to evaluate the association between gene expression and survival. This script applies survival modelling approaches, including univariate Cox proportional hazards analysis and Kaplan-Meier visualization for prioritized genes. Its purpose is to identify differentially expressed genes with potential prognostic value in the TNBC cohort.

## 06_Somatic_mutation_analysis_and_filtering_of_truncating_variants.R

Downloads and processes TCGA-BRCA masked somatic mutation data, focusing on mutations present in TNBC samples. This script filters the mutation dataset to retain variants affecting the genes prioritized in the expression and survival analyses, while excluding truncating or obvious loss-of-function mutation classes. The aim is to focus on non-truncating variants that may contribute to altered gene regulation rather than directly disrupting protein coding capacity.

## 07_Identification_of_regulatory_regions_affected_by_non_truncating_mutations.R

Transforms the prioritized non-truncating mutation set into genomic coordinates compatible with RBioGateway and searches for overlap with cis-regulatory modules. This script introduces the regulatory interpretation layer of the workflow by linking candidate TNBC-associated mutations to putative regulatory DNA regions catalogued in BioGateway.

## 08_Identification_of_TFBS_affected_in_prioritized_regulatory_regions.R

Uses RBioGateway to retrieve transcription factors associated with the affected regulatory regions and summarizes mutation–CRM–transcription factor relationships. This stage prioritizes candidate regulatory events that may affect transcription factor binding and contribute to the dysregulation of survival-associated genes in TNBC.

## 09_GO_enrichment_of_identified_regulatory_regions.R
Performs Gene Ontology (GO) enrichment analysis on the regulatory regions prioritized in the previous steps. This script functionally interprets the affected CRMs and their associated genes or transcription factors by identifying overrepresented biological processes, molecular functions, and cellular components potentially related to TNBC biology. The aim of this final analytical stage is to determine whether the regulatory regions affected by candidate non-truncating mutations are enriched in functional categories relevant to tumor progression, gene regulation, and patient prognosis.

# Current scope
The current version of the workflow implements a complete integrative analysis from TCGA-BRCA cohort definition to functional interpretation of candidate regulatory events. The pipeline sequentially identifies clinically defined TNBC patients, retrieves and preprocesses transcriptomic and mutational data, detects differentially expressed genes, evaluates their prognostic relevance, filters non-truncating mutations, maps them to cis-regulatory modules, explores potentially affected TFBS-associated regions, and performs GO enrichment analysis to contextualize the biological relevance of the prioritized regulatory findings.
