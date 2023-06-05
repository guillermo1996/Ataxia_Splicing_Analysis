## _________________________________________________
##
## Ataxia Level 3 (Diagnosis) Analysis - Cerebellum
##
## Aim: calculate the mis-splicing ratio of the annotated junction after
## pseudobulking by tissue and ataxia diagnosis (level 3).
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
## focus on Cerebellum Level 3 studies.
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

logger::log_info("Starting Analysis for Cerebellum Level 3.")

## Relevant Paths
if(!exists("results_path")) results_path <- here::here("results/")
project_path <- file.path(results_path, "Cerebellum_Level_3/")

dir.create(project_path, showWarnings = F, recursive = T)

## Script parameters
level = "Diagnosis"
all_diagnosis <- c("FRDA", "SCA1", "SCA2", "SCA6")

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
metadata_cerebellum <- metadata %>% dplyr::filter(Region == "Cerebellum") 
variance_df <- getVarianceDf(metadata_project = metadata_cerebellum,
                             results_path = results_path,
                             covariates = c("RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex"),
                             response_var = "mapped_junctions")
metadata_project <- metadata_cerebellum 

# Pipeline for each possible diagnosis ----
foreach(i = seq_along(all_diagnosis)) %do%{
  diagnosis <- all_diagnosis[i]
  diagnosis_path <- paste0(project_path, diagnosis, "/")
  dir.create(diagnosis_path, recursive = T, showWarnings = F)
  
  clusters <- c(diagnosis, "Control")
  metadata_diagnosis <- metadata_project %>%
    dplyr::filter(!!sym(level) %in% clusters)
  
  ## Subsample using Gower distance ----
  metadata_subsample <- subsampleGowerDistance(metadata_project = metadata_diagnosis, 
                                               level = level,
                                               id_field = "ID_anon",
                                               covariates = c("RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex"),
                                               clusters = clusters, 
                                               weights = variance_df)
  
  ### Plot and store the distributions in disk. Remove if not necessary. The plot
  ### function is wrapped in "suppressMessages" to ignore a warning from "ggplot".
  output_figure = paste0(project_path, "metadata_distributions.png")
  suppressMessages(print(plotMetadataSubsample(metadata_subsample, level, output_file = output_figure, ratio = 1.2)))
  
  ## Calculation of the Mis-splicing ratio ----
  projectAnalysis(metadata_project = metadata_subsample,
                  project_path = diagnosis_path,
                  level = level,
                  clusters = clusters,
                  rw_disk = T,
                  overwrite = F)
  
  ## Wilcox paired signed-rank test ----
  logger::log_info("\t Loading the common introns.")
  common_introns <- getCommonIntrons(diagnosis_path)
  common_novel <- getCommonNovel(diagnosis_path, common_introns)
  
  logger::log_info("\t Executing the Wilcoxon paired signed-rank test.")
  wilcox_test <- MSRanalysis(common_introns, diagnosis_path, clusters)
}