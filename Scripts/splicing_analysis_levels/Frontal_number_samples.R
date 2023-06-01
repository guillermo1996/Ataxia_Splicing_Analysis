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
project_path <- file.path(results_path, "Frontal_sample_number/")

dir.create(project_path, showWarnings = F, recursive = T)

## Script parameters
level = "Type"
clusters = c("Case", "Control")

doParallel::registerDoParallel(4)

## Load the initial information
metadata_path <- here::here("metadata/metadata.csv") 
multiqc_path <- here::here("metadata/multiqc_rseqc_read_distribution.txt")
metadata <- readr::read_delim(metadata_path, show_col_types = FALSE) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(Individual_ID = stringr::str_split(ID_anon, "_", simplify = T)[1]) %>%
  dplyr::filter(!(Diagnosis %in% c("CANVAS", "AIFM1"))) %>%
  extractReadDepthMultiQC(multiqc_path)

# Metadata & clustering ----
metadata_frontal_control <- metadata %>% dplyr::filter(Region == "Frontal", Type == "Control", RIN >= 4)  
variance_df <- getVarianceDf(metadata_frontal_control, 
                             results_path = results_path,
                             output_file = here::here("variables/variance_explained_frontal_control.rds"))
metadata_project <- metadata_frontal_control %>%
  dplyr::mutate(Type = ifelse(ID_anon %in% c("CO15B_Frontal", "CO11A_Frontal", "CO20B_Frontal"), "Case", "Control"))

# Subsample and MSR for every value of N
foreach(n = 1:6) %do%{
  logger::log_info("Starting Wilcoxon paired signed-rank test for N = ", n, ".")
  metadata_subsample <- subsampleGowerDistance(metadata_project = metadata_project, level = level, clusters = clusters, weights = variance_df, n = n)
  project_n_path <- file.path(project_path, paste0("N", n, "/"))
  dir.create(project_n_path, showWarnings = F)
  
  output_figure = paste0(project_n_path, "metadata_distributions.png")
  plotMetadataSubsample(metadata_subsample, level, output_file = output_figure, ratio = 1.2)
  
  projectAnalysis(metadata_project = metadata_subsample,
                  project_path = project_n_path,
                  level = level,
                  clusters = clusters,
                  rw_disk = T,
                  overwrite = F)
  
  ## Wilcoxon test
  logger::log_info("\t Loading the common introns.")
  common_introns <- getCommonIntrons(project_n_path)
  logger::log_info("\t Executing the Wilcoxon paired signed-rank test.")
  wilcox_test <- MSRanalysis(common_introns, project_n_path, clusters = c("Case", "Control"))
  
  ## Mean MSR & Other stats
  splice_sites = c("Donor", "Acceptor")
  zero_prop <- foreach(j = seq(splice_sites)) %do%{
    splice_site <- splice_sites[j]
    
    MSR_Table <- common_introns %>%
      dplyr::select(ref_junID, MSR_Donor, MSR_Acceptor, cluster) %>%
      tidyr::pivot_wider(id_cols = ref_junID, names_from = c("cluster"), values_from = c(paste0("MSR_", splice_site)))
    
    MSR_Case <- MSR_Table %>% dplyr::arrange(ref_junID) %>% dplyr::pull(Case)
    MSR_Control <- MSR_Table %>% dplyr::arrange(ref_junID) %>% dplyr::pull(Control)
    
    Case_zero_prop <- sum(MSR_Case == 0)/length(MSR_Case)
    Control_zero_prop <- sum(MSR_Control == 0)/length(MSR_Control)
    
    tibble(splice_site = splice_site,
           Case_zero_prop = Case_zero_prop,
           Control_zero_prop = Control_zero_prop,
           Case_mean_MSR = mean(MSR_Case),
           Control_mean_MSR = mean(MSR_Control),
           Case_nonZeroMean_MSR = mean(MSR_Case[MSR_Case != 0]),
           Control_nonZeroMean_MSR = mean(MSR_Control[MSR_Control != 0]))
  } %>% dplyr::bind_rows() %>% dplyr::mutate(N = n)
  
  zero_prop %>% saveRDS(paste0(project_n_path, "zero_prop.rds"))
  return(c(wilcox_test, zero_prop))
}
