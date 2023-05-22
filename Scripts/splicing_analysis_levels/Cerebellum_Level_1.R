## _________________________________________________
##
## Script title
##
## Aim: 
##
## Author: Mr. Guillermo Rocamora Pérez
##
## Date Created: 
##
## Copyright (c) Guillermo Rocamora Pérez, year
##
## Email: guillermorocamora@gmail.com
## _________________________________________________
##
## Notes:
##
##
## Please contact guillermorocamora@gmail.com for further assitance.
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

## Relevant Paths
if(!exists("results_path")) results_path <- here::here("results/")
project_path <- file.path(results_path, "Cerebellum_Level_1/")

dir.create(project_path, showWarnings = F, recursive = T)

## Script parameters
level = "Type"
clusters = c("Case", "Control")

doParallel::registerDoParallel(4)

## Load the initial information
if(!exists("metadata")){
  metadata_path <- here::here("metadata/metadata.csv") 
  multiqc_path <- here::here("metadata/multiqc_rseqc_read_distribution.txt")
  metadata <- readr::read_delim(metadata_path, show_col_types = FALSE) %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(Individual_ID = stringr::str_split(ID_anon, "_", simplify = T)[1]) %>%
    dplyr::filter(!(Diagnosis %in% c("CANVAS", "AIFM1"))) %>%
    extractReadDepthMultiQC(multiqc_path)   # Add total reads information from MultiQC
}

# Metadata & clustering ----
metadata_project <- metadata %>% 
  dplyr::filter(Region == "Cerebellum") %>%
  dplyr::filter(!(Type == "Control" & RIN >= 7))

## Subsample using Gower distance ----
metadata_subsample <- subsampleGowerDistance(metadata_project = metadata_project,
                                             level = level,
                                             clusters = clusters)

output_figure = paste0(project_path, "metadata_distributions.png")
plotMetadataSubsample(metadata_subsample, level, output_file = output_figure, ratio = 1.2)

# Calculation of the Mis-splicing ratio ----
projectAnalysis(metadata_project = metadata_subsample,
                project_path = project_path,
                level = level,
                clusters = clusters,
                rw_disk = T,
                overwrite = T)

# Wilcox paired signed-rank test ----
logger::log_info("\t Loading the common introns.")
common_introns <- getCommonIntrons(project_path)

logger::log_info("\t Executing the Wilcoxon paired signed-rank test.")
wilcox_test <- MSRanalysis(common_introns, project_path, clusters = c("Case", "Control"))