## _________________________________________________
##
## Ataxia Level 2 (AtaxiaSubtype) Analysis - Cerebellum
##
## Aim: calculate the mis-splicing ratio of the annotated junction after
## pseudobulking by tissue and ataxia subtype (level 2).
##
## Author: Mr. Guillermo Rocamora Pérez
##
## Date Created: 29/05/2023
##
## Copyright (c) Guillermo Rocamora Pérez, 2023
##
## Email: guillermorocamora@gmail.com
## _________________________________________________
##
## Notes:
##
## This script is the third step in the Splicing Noise Analysis. We must first
## generate the junction files (see "download_and_extraction.R" script) and the
## main pipeline script (see "splicing_noise_analysis.R"). In this script, we
## focus on Cerebellum Level 2 studies.
##
## Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________


# Initial setup ----

## Required libraries
shhh <- suppressPackageStartupMessages
shhh(library(GenomicRanges))
shhh(library(here))
shhh(library(doParallel))
shhh(library(foreach))
shhh(library(tidyverse))
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")

## Source helper functions
source(here::here("R/hf_project_analysis.R"))
source(here::here("R/hf_graph_and_themes.R"))
source(here::here("R/hf_subsample.R"))
source(here::here("R/hf_additional.R"))

## Relevant Paths
if(!exists("results_path")) results_path <- here::here("results/")
project_path <- file.path(results_path, "Cerebellum_Level_2/")

dir.create(project_path, showWarnings = F, recursive = T)

## Script parameters
level = "AtaxiaSubtype"
ataxia_subtypes = c("KnownAtaxia", "UnknownAtaxia")

doParallel::registerDoParallel(4)

## Load the initial information
if(!exists("metadata", inherits = F)){
  metadata_path <- here::here("metadata/metadata.csv") 
  multiqc_path <- here::here("metadata/multiqc_rseqc_read_distribution.txt")
  metadata <- readr::read_delim(metadata_path, show_col_types = FALSE) %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(Individual_ID = stringr::str_split(ID_anon, "_", simplify = T)[1]) %>%
    dplyr::filter(!(Diagnosis %in% c("CANVAS", "AIFM1"))) %>%
    extractReadDepthMultiQC(multiqc_path)   # Add total reads information from MultiQC
}

# Metadata & subsampling weights ----
metadata_cerebellum <- metadata %>%  dplyr::filter(Region == "Cerebellum") 
variance_df <- getVarianceDf(metadata_project = metadata_cerebellum,
                             results_path = results_path,
                             covariates = c("RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex"),
                             response_var = "mapped_junctions")
metadata_project <- metadata_cerebellum 

# Pipeline for each possible diagnosis ----
foreach(i = seq_along(ataxia_subtypes)) %do%{
  subtype <- ataxia_subtypes[i]
  subtype_path <- paste0(project_path, subtype, "/")
  dir.create(subtype_path, recursive = T, showWarnings = F)
  
  clusters <- c(subtype, "Control")
  metadata_subtype <- metadata_project %>%
    dplyr::filter(!!sym(level) %in% clusters)
  
  ## Subsample using Gower distance ----
  metadata_subsample <- subsampleGowerDistance(metadata_project = metadata_subtype, 
                                               level = level,
                                               id_field = "Individual_ID",
                                               covariates = c("RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex"),
                                               clusters = clusters, 
                                               weights = variance_df)
  output_figure = paste0(subtype_path, "metadata_distributions.png")
  plotMetadataSubsample(metadata_subsample, level, output_file = output_figure, ratio = 1.2)
  
  ## Calculation of the Mis-splicing ratio ----
  projectAnalysis(metadata_project = metadata_subsample,
                  project_path = subtype_path,
                  level = level,
                  clusters = clusters,
                  rw_disk = T,
                  overwrite = F)
  
  ## Wilcox paired signed-rank test ----
  logger::log_info("\t Loading the common introns.")
  common_introns <- getCommonIntrons(subtype_path)
  common_novel <- getCommonNovel(subtype_path, common_introns)
  
  logger::log_info("\t Executing the Wilcoxon paired signed-rank test.")
  wilcox_test <- MSRanalysis(common_introns, subtype_path, clusters)
}