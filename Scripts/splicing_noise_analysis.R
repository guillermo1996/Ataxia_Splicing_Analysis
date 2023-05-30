## _________________________________________________
##
## Splicing Noise Analysis
##
## Aim: main script to run the splicing noise analysis. It loads the junctions
## previously extracted from BAM files, annotates them and generates pair of
## annotated junction - novel junction.
##
## The output depends on the clustering of samples, but a procedure common to
## all samples is executed in this script.
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
## This script is the second step in the Splicing Noise Analysis. We must first
## generate the junction files (see "download_and_extraction.R" script). Then,
## different scripts are executed for different studies or clustering
## approaches. In this script, a total of 6 analyses (2 tissues in three
## different levels) are studied.
##
## Please contact guillermorocamora@gmail.com for further assistance.
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
samtools_path = "/home/grocamora/tools/samtools/bin/"
fordownload_path = "/home/grocamora/tools/fordownload/"
blacklist_path = "/home/grocamora/RytenLab-Research/Additional_files/hg38-blacklist.v2.bed"
u12_introns_path = "/home/grocamora/RytenLab-Research/Additional_files/minor_introns_tidy.rds"
u2_introns_path = "/home/grocamora/RytenLab-Research/Additional_files/major_introns_tidy.rds"

# fasta_path = "/home/grocamora/RytenLab-Research/Additional_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# gtf_path = "/home/grocamora/RytenLab-Research/Additional_files/Homo_sapiens.GRCh38.105.gtf"
gtf_path = "/home/grocamora/RytenLab-Research/Additional_files/GENCODE/gencode.v39.annotation.gtf"
fasta_path = "/home/grocamora/RytenLab-Research/Additional_files/GENCODE/GRCh38.primary_assembly.genome.fa"

## Script parameters
rw_disk = T       # Whether to write and read from disk
num_cores = 4     # Number of multiprocessing cores to use
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
                                           samtools_path = samtools_path,
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
## Level 1 (Type) ----
source(here::here("Scripts/splicing_analysis_levels/Cerebellum_Level_1.R"))
source(here::here("Scripts/splicing_analysis_levels/Frontal_Level_1.R"))
## Level 2 (AtaxiaSubtype) ----
source(here::here("Scripts/splicing_analysis_levels/Cerebellum_Level_2.R"))
source(here::here("Scripts/splicing_analysis_levels/Frontal_Level_2.R"))
## Level 3 (Diagnosis) ----
source(here::here("Scripts/splicing_analysis_levels/Cerebellum_Level_3.R"))
source(here::here("Scripts/splicing_analysis_levels/Frontal_Level_3.R"))