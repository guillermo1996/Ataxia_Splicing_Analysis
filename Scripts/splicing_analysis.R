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
shhh(library(logger))
shhh(library(GenomicRanges))
shhh(library(here))
shhh(library(doParallel))
shhh(library(foreach))
shhh(library(tidyverse))
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")

## Logger options
log_file <- here::here("logs/Ataxia_splicing_analysis.log")
logger::log_appender(logger::appender_tee(log_file, append = T))
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)

## Source helper functions
source(here::here("R/hf_additional.R"))
source(here::here("R/hf_junction_extraction.R"))
source(here::here("R/hf_junction_annotation.R"))
source(here::here("R/hf_distance_calculation.R"))
source(here::here("R/hf_never_misspliced_junctions.R"))
source(here::here("R/hf_database_generation.R"))

## Relevant Paths
metadata_path <- here::here("metadata/metadata.csv") 
multiqc_path <- here::here("metadata/multiqc_rseqc_read_distribution.txt")
main_samples_path <- here::here("data/samples/")
results_path <- here::here("results/")
dir.create(results_path, showWarnings = F, recursive = T)

bedtools_path = "/home/grocamora/tools/bedtools/"
fasta_path = "/home/grocamora/RytenLab-Research/Additional_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fordownload_path = "/home/grocamora/tools/fordownload/"
gtf_path = "/home/grocamora/RytenLab-Research/Additional_files/Homo_sapiens.GRCh38.105.gtf"
blacklist_path = "/home/grocamora/RytenLab-Research/Additional_files/hg38-blacklist.v2.bed"
u12_introns_path = "/home/grocamora/RytenLab-Research/Additional_files/minor_introns_tidy.rds"
u2_introns_path = "/home/grocamora/RytenLab-Research/Additional_files/major_introns_tidy.rds"

## Script parameters
rw_disk = T       # Whether to write and read from disk
num_cores = 4     # Number of multprocessing cores to use
overwrite =  F    # Whether to overwrite previously generated results

## Load the initial information
metadata <- readr::read_delim(metadata_path, show_col_types = FALSE) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(Individual_ID = stringr::str_split(ID_anon, "_", simplify = T)[1]) %>%
  dplyr::filter(!(Diagnosis %in% c("CANVAS", "AIFM1"))) %>%
  extractReadDepthMultiQC(multiqc_path)   # Add total reads information from MultiQC

# Main Pipeline ----
logger::log_info("Starting splicing noise analysis for Ataxia Bulk samples.")

## Read junctions from disk ----
all_reads_combined <- junctionReading(metadata = metadata, 
                                      main_samples_path = main_samples_path,
                                      num_cores = num_cores,
                                      rw_disk = rw_disk,
                                      overwrite = overwrite,
                                      output_path = results_path)

## Junction annotation ----
# Only ENSEMBL gtf are supported as of 19/05/2023
annotated_SR_details <- junctionAnnotation(all_reads_combined = all_reads_combined, 
                                           main_samples_path = main_samples_path,
                                           blacklist_path = blacklist_path,
                                           gtf_path = gtf_path,
                                           bedtools_path = bedtools_path,
                                           fasta_path = fasta_path,
                                           fordownload_path = fordownload_path,
                                           rw_disk = rw_disk,
                                           overwrite = overwrite,
                                           output_path = results_path)

## Distance calculations ----
### Raw distances
all_distances_raw <- juctionPairing(metadata = metadata,
                                    main_samples_path = main_samples_path,
                                    annotated_SR_details = annotated_SR_details,
                                    all_reads_combined = all_reads_combined,
                                    num_cores = num_cores,
                                    rw_disk = rw_disk,
                                    overwrite = overwrite,
                                    output_path = results_path)

### Distance pruning 
all_distances_pruned <- removeAmbiguousPairing(all_distances_raw,
                                               main_samples_path,
                                               rw_disk = rw_disk,
                                               overwrite = overwrite,
                                               output_path = results_path)

# Project specific pipeline ----
projectAnalysis <- function(metadata_project, 
                            project_path, 
                            level, 
                            clusters,
                            rw_disk = T,
                            overwrite = F){
  purrr::map_df(clusters, function(cluster){
    if (rw_disk & !overwrite & 
        file.exists(paste0(project_path, cluster, "_db_introns.rds")) &
        file.exists(paste0(project_path, cluster, "_db_novel.rds"))){
      logger::log_info("\t Ignoring cluster ", cluster, ". Database already generated.")
      return()
    }
    
    logger::log_info("\t Starting cluster ", cluster, ".")
    
    # Load the metadata and samples for the cluster
    metadata_cluster <- metadata_project %>%
      dplyr::filter(!!sym(level) == cluster)
    cluster_samples <- metadata_cluster$Correct_sample %>% unique
    logger::log_info("\t\t Number of samples in cluster: ", length(cluster_samples), ".")
    
    # Extract the split reads, annotated junctions and distances for the cluster
    # samples
    cluster_split_reads <- all_reads_combined %>%
      dplyr::select(junID, all_of(cluster_samples)) %>%
      dplyr::filter(if_any(all_of(cluster_samples), ~ !is.na(.)))
    cluster_annotated_SR_details <- annotated_SR_details %>%
      dplyr::filter(junID %in% cluster_split_reads$junID)
    cluster_distances_raw <- all_distances_raw %>%
      dplyr::filter(sample %in% cluster_samples)
    
    # Extract the distances pruned and generate the tidy version (adding the
    # junction counts to the distances)
    cluster_distances_pruned <- all_distances_pruned %>%
      dplyr::filter(sample %in% cluster_samples) %>%
      dplyr::distinct(novel_junID, ref_junID, .keep_all = T) %>%
      dplyr::select(-sample, -ref_counts, -novel_counts) %>%
      dplyr::relocate(novel_junID, ref_junID)
    cluster_distances_tidy <- addCounts(cluster_samples,
                                        cluster_distances_pruned,
                                        cluster_split_reads)
    
    # Combine the never mis-spliced junctions with the tidy distances
    cluster_distances_tidy_all <- addNeverMissplicedJunction(cluster_samples,
                                                             cluster_split_reads,
                                                             cluster_distances_raw,
                                                             cluster_annotated_SR_details,
                                                             cluster_distances_tidy,
                                                             project_path,
                                                             cluster,
                                                             rw_disk = F,
                                                             overwrite = F)
    
    # Generate the Database
    generateDB(cluster_distances_tidy_all,
               cluster_annotated_SR_details,
               cluster_path = project_path,
               cluster_name = cluster,
               u12_introns_path = u12_introns_path,
               u2_introns_path = u2_introns_path,
               rw_disk = T,
               overwrite = T)
    
    # Clear some memory (should be done automatically)
    rm(metadata_cluster, cluster_samples, cluster_split_reads, cluster_annotated_SR_details,
       cluster_distances_raw, cluster_distances_pruned, cluster_distances_tidy, cluster_distances_tidy_all)
  })
}

## Cerebellum - Level 1 (Type)
project_path <- file.path(results_path, "Cerebellum_Level1/")
level = "Type"
clusters = c("Case", "Control")
metadata_project <- metadata %>% filter(Region == "Cerebellum")
dir.create(project_path, showWarnings = F, recursive = T)

projectAnalysis(metadata_project = metadata_project,
                project_path = project_path,
                level = "Type",
                clusters = c("Case", "Control"))



## Cerebellum - Level 2 (AtaxiaSubtype)

## Cerebellum - Level 3 (Diagnosis)

## Frontal Cortex - Level 1 (Type)
project_path <- file.path(results_path, "Frontal_Level1/")
level = "Type"
clusters = c("Case", "Control")
metadata_project <- metadata %>% filter(Region == "Frontal")
dir.create(project_path, showWarnings = F, recursive = T)

projectAnalysis(metadata_project = metadata_project,
                project_path = project_path,
                level = level,
                clusters = clusters)

## Frontal Cortex - Level 2 (AtaxiaSubtype)

## Frontal Cortex - Level 3 (Diagnosis)