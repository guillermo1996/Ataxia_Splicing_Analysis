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

## Relevant Paths
if(!exists("results_path")) results_path <- here::here("results/")
project_path <- file.path(results_path, "Cerebellum_Level1/")

dir.create(project_path, showWarnings = F, recursive = T)

## Script parameters
level = "Type"
clusters = c("Case", "Control")

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

metadata_project <- metadata %>% filter(Region == "Cerebellum")
