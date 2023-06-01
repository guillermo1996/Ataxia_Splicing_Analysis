## _________________________________________________
##
## Download and extraction of BAM files
##
## Aim: to download the BAM files of Ataxia RNAseq samples from AWS and to
## extract both the junctions within the files and the coverage of Ataxia
## specific genes.
##
## Author: Mr. Guillermo Rocamora Pérez
##
## Date Created: 2023-02-07
##
## Copyright (c) Guillermo Rocamora Pérez, 2023
##
## Email: guillermorocamora@gmail.com
## _________________________________________________
##
## Notes:
##
## To download the files from the AWS server, it is necessary to load the
## credentials. This script requires the user to run first the script
## "login_aws.sh" and follow the instructions. Some security parameters need to
## be input in the "generate_crendetials.R" script.
##
## Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

# Initial setup ----

## Required libraries
shhh <- suppressPackageStartupMessages
shhh(library(logger, warn.conflicts = F))
shhh(library(foreach))
shhh(library(tidyverse, warn.conflicts = F))
shhh(library(aws.s3))
shhh(library(GenomicAlignments))
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")

## Custom source functions
source(here::here("R/hf_download_bam_files.R"))
source(here::here("R/hf_extract_coverage.R"))

## Logger options
log_file <- here::here("logs/download_and_extraction.log")
logger::log_appender(logger::appender_tee(log_file, append = T))
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)

## Relevant Paths
samtools_path <- file.path("/home/grocamora/tools/samtools/bin/")
regtools_path <- file.path("/home/grocamora/tools/regtools/build/")
s3md5_path = file.path("/home/grocamora/tools/aws-s3-integrity-check/s3md5/")

coverage_metadata_path <- file.path("/home/grocamora/Core_Projects/Ataxia_Splicing_Analysis/metadata/All_patient_coverage.csv")
coverage_bam_files_path <- file.path("/home/grocamora/Core_Projects/Ataxia_Splicing_Analysis/data/Gene_filtered_BAM/")
main_samples_path <- file.path("/home/grocamora/Core_Projects/Ataxia_Splicing_Analysis/data/samples/")
backup_samples_path <- file.path("/home/grocamora/Core_Projects/Ataxia_Splicing_Analysis/data/samples_wrong_name/")

metadata_path <- file.path("/home/grocamora/Core_Projects/Ataxia_Splicing_Analysis/metadata/metadata.csv")
credentials_path <- file.path("/home/grocamora/Core_Projects/Ataxia_Splicing_Analysis/variables/credentials/role_credentials.json")

dir.create(coverage_bam_files_path, recursive = T, showWarnings = F)
dir.create(backup_samples_path, recursive = T, showWarnings = F)
dir.create(main_samples_path, recursive = T, showWarnings = F)

## Script parameters
samtools_threads <- 2
samtools_memory <- "20G"
overwrite = F
minimal_credential_minutes <- 60

## Load the initial information
metadata <- readr::read_delim(metadata_path, show_col_types = FALSE)
ataxia_genes <- loadAtaxiaGenes()
expiration_time <- loadCredentials(credentials_path)

# Pipeline ----

## Loop through every row in the metadata table to locate the sample file on AWS
## and execute the pipeline
r <- foreach(i = seq(nrow(metadata))) %do%{
  ## Load the sample's metadata
  row <- metadata[i, ]
  sample_name <- row$sample_name
  sample_correct_name <- row$Correct_sample
  sample_diagnosis <- row$Diagnosis
  sample_type <- row$Type
  
  logger::log_info("Starting download and extraction process for ", sample_type, " sample ", sample_name, " (", sample_correct_name, ").")
  
  ## First, we check if the time left to download AWS files is less than 60
  ## minutes. If it is, Stop the download and junction process rather than
  ## interrupting the BAM download process because of expired credentials.
  if(estimateTimeLeft(expiration_time, minimal_credential_minutes)){
    logger::log_warn("\t Credentials are nearly expired. Download process stopped.")
    return()
  }
  
  ## Generate the files path
  download_link <- getSampleDownloadLink(sample_name, sample_type)
  file_path <- generateSamplePath(sample_name, sample_type, backup_samples_path)
  sort_path <- paste0(file_path, ".sort")
  junc_path <- paste0(file_path, ".sort.s0.junc")
  junc_correct_name_path <- paste0(main_samples_path, tolower(sample_type), "/", sample_correct_name, ".bam.sort.s0.junc")
  
  ## Check if BAM file is already extracted
  if (!overwrite & file.exists(junc_correct_name_path)) {
    logger::log_info("\t Ignoring extraction and download. Junction file ", sample_correct_name, " already found!\n\n")
    return()
  }
  
  ## Previously, if the BAM file was from a FRDA sample, we would load it from
  ## previously downloaded folder.
  
  # ## Change the BAM path if sample is FRDA
  # if(sample_diagnosis == "FRDA"){
  #   file_path <- paste0("/home/zchen/ataxia-rnaseq/Input/BAM/FRDA", generateSamplePath(sample_name, "", ""))
  #   logger::log_info("\t FRDA sample found.")# File path location defaulted at ", file_path)
  # }
  
  ## Download the BAM file if needed ----
  logger::log_info("\t Starting the BAM download process.")
  
  if(!file.exists(file_path)){
    aws.s3::save_object(download_link, file = file_path, overwrite = overwrite)
  }else{
    logger::log_info("\t\t BAM file already existst. Starting an integrity check.")
    integrity_check <- checkIntegrity(download_link, file_path, s3md5_path)
    
    ## If the integrity does not check, I rather not doing anything as of now.
    ## In the future, it should remove the old file and update with the new
    ## (redownload).
    if(!integrity_check) return()
  }
  
  ## Zhongbo's Script to extract the coverage
  logger::log_info("\t Starting the coverage extraction process.")
  extractCoverage(file_path,
                  sample_correct_name,
                  sample_diagnosis,
                  coverage_metadata_path,
                  coverage_bam_files_path)
  
  ## If the .sort files exists, we want to remove them. It is better to start
  ## the process again.
  if(file.exists(sort_path)){
    system2(command = "rm", args = c(sort_path))
    system2(command = "rm", args = c(paste0(sort_path, ".bai")))
  }
  
  ## Junction extraction process ----
  logger::log_info("\t Starting the sorting process.")
  system2(command = paste0(samtools_path, "samtools"), args = c(
    "sort", file_path,
    "-o", sort_path,
    "--threads", samtools_threads,
    "-m", samtools_memory
  ))
  if(!file.exists(sort_path)) return("Sorting failed")
  
  logger::log_info("\t Starting the indexing process.")
  system2(command = paste0(samtools_path, "samtools"), args = c(
    "index", sort_path,
    "-@", floor(samtools_threads/2)
  ))
  
  logger::log_info("\t Starting the extraction process.")
  system2(command = paste0(regtools_path, "regtools"), args = c(
    "junctions extract", sort_path,
    "-m 25",
    "-M 1000000",
    "-s 0",
    "-o", junc_path
  ))
  
  ## Remove the files that are not necessary----
  if(!file.exists(junc_path)){
    logger::log_error("\t No junction file found after the process. An error must have ocurred.")
    logger::log_error("\t Removing all intermediary files. Please run again the analysis.")
    
    #system2(command = "rm", args = c(file_path))
    system2(command = "rm", args = c(sort_path))
    system2(command = "rm", args = c(paste0(sort_path, ".bai")))
  }else{
    logger::log_info("\t Successfully extracted the junctions. Removing all files (BAM included).")
    system2(command = "rm", args = c(file_path))
    system2(command = "rm", args = c(sort_path))
    system2(command = "rm", args = c(paste0(sort_path, ".bai")))
  }
  
  ## Copy the extracted junction and rename to match sample name ----
  system2(command = "cp", args = c(junc_path, junc_correct_name_path))
  if(file.exists(junc_correct_name_path)){
    logger::log_info("\t Junction file with correct name stored (", sample_correct_name , ")\n\n")
  }else{
    logger::log_error("\t Junction file with correct name not found! An error must have ocurred.\n\n")
  }
} 