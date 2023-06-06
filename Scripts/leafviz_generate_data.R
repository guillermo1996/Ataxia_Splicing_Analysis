## _________________________________________________
##
## Leafviz - Prepare Results
##
## Aim: convert the results from Leafcutter to an .RData file to use leafviz.
##
## Author: Mr. Guillermo Rocamora Pérez
##
## Date Created: 05/06/2023
##
## Copyright (c) Guillermo Rocamora Pérez, 2023
##
## Email: guillermorocamora@gmail.com
## _________________________________________________
##
## Notes:
##
## This scripts requires a very specific folder structure and naming convenction
## of the Leafcutter results to work. All files required can be found in
## /home/grocamora/Core_Projects/Ataxia-Splicing-Analysis/Leafviz.
##
## A specific docker image from user "lroser" was also employed:
## \href{link}{https://hub.docker.com/r/lroser/run_leafcutter}.
##
## This script generates the executable file to run on the docker container. As
## such, it later requires the user with docker usage rights to execute the
## following command:
##
## `docker run -it --rm -v /home/grocamora/Core_Projects/Ataxia_Splicing_Analysis/Leafviz:/Leafviz lroser/run_leafcutter bash -c "/Leafviz/prepare_results.sh"`
## 
## Results still require additional treatment to be used in leafviz. You can
## find a follow-up script in "leafviz_fix_data.R"
##
## Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

# Initial setup ----

## Required libraries
shhh <- suppressPackageStartupMessages
shhh(library(foreach))
shhh(library(tidyverse))
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")

## Helper functions

#' Extract "_cluster_significance.txt" path
#'
#' @param level Numeric, type of study. Required to find the proper data file.
#' @param tissue Character, tissue being studied. Required to find the proper
#'   data file.
#' @param comb_name Character, comparison being studied. Required to find the
#'   proper data file.
#'
#' @return
#' @export
getClusterSignificance <- function(level, tissue, comb_name){
  if(level == 1){
    return(paste0("/Leafviz/Results/Type", level, "/", comb_name, "_cluster_significance.txt"))
  }else{
    return(paste0("/Leafviz/Results/Type", level, "/", tissue, "/", comb_name, "_cluster_significance.txt"))
  }
}

#' Extract "_effect_sizes.txt" path
#'
#' @param level Numeric, type of study. Required to find the proper data file.
#' @param tissue Character, tissue being studied. Required to find the proper
#'   data file.
#' @param comb_name Character, comparison being studied. Required to find the
#'   proper data file.
#'
#' @return
#' @export
getEffectSize <- function(level, tissue, comb_name){
  if(level == 1){
    return(paste0("/Leafviz/Results/Type", level, "/", comb_name, "_effect_sizes.txt"))
  }else{
    return(paste0("/Leafviz/Results/Type", level, "/", tissue, "/", comb_name, "_effect_sizes.txt"))
  }
}


#' Extract "_group_file.txt" path
#'
#' @param level Numeric, type of study. Required to find the proper data file.
#' @param tissue Character, tissue being studied. Required to find the proper
#'   data file.
#' @param comb_name Character, comparison being studied. Required to find the
#'   proper data file.
#'
#' @return
#' @export
getMetadataPath <- function(level, tissue, comb_name){
  if(level == 1){
    return(paste0("/Leafviz/group_files/", tissue, "/Case_vs_Control_group_file.txt"))
  }else if(level == 2){
    return(paste0("/Leafviz/group_files/", tissue, "/Type", level, "/", comb_name, "_group_file.txt"))
  }else{
    return(paste0("/Leafviz/group_files/", tissue, "/Type", level, "/reordered/", comb_name, "_group_file.txt"))
  }
}


#' Generates .RData output path
#'
#' @param level Numeric, type of study. Required to find the proper data file.
#' @param tissue Character, tissue being studied. Required to find the proper
#'   data file.
#' @param comb_name Character, comparison being studied. Required to find the
#'   proper data file.
#'
#' @return
#' @export
getOutputPath <- function(level, tissue, comb_name){
  if(level == 1){
    return(paste0("/Leafviz/Leafviz_data/type", level, "_", comb_name, ".RData"))
  }else{
    return(paste0("/Leafviz/Leafviz_data/type", level, "_", tissue, "_", comb_name, ".RData"))
  }
}


#' Write the "prepare_results.R" command into file
#'
#' @param script_path Character, path to where to store the output executable.
#' @param level Numeric, type of study. Required to find the proper data file.
#' @param tissue Character, tissue being studied. Required to find the proper data file.
#' @param comb_name Character, comparison being studied. Required to find the proper data file.
#'
#' @return
#' @export
generateProgramLine <- function(script_path, level, tissue, comb_name){
  cat("/home/leafviz/prepare_results.R", 
      "/Leafviz/diseasegroups_perind_numers.counts.gz",
      getClusterSignificance(level, tissue, comb_name),
      getEffectSize(level, tissue, comb_name),
      "/Leafviz/annotation_codes/gencode_hg38/gencode_hg38",
      paste0("-m ", getMetadataPath(level, tissue, comb_name)),
      paste0("-o ", getOutputPath(level, tissue, comb_name)),
      "\n",
      file = script_path,
      append = T)
}

# Create script ----
leafviz_path <- here::here("Leafviz")
script_path <- file.path(leafviz_path, "prepare_results.sh")
tissues <- c("Cerebellum", "Frontal")

## By default, "less" is not installed in the docker image, but it is required
## for the "prepare_results.R" script to work properly. This command install it.
cat("apt install less \n", file = script_path)

## Type 1 ----
level = 1
foreach(i = seq_along(tissues)) %do%{
  tissue <- tissues[i]
  comb_name <- stringr::str_c(tissue, "_Case_vs_Control")
  generateProgramLine(script_path, level, tissue, comb_name)
}

## Type 2 ----
ataxia_subtypes <- c("KnownAtaxia", "UnknownAtaxia", "Control")
type2_comb <- combn(ataxia_subtypes, 2, simplify = T) %>% apply(2, function(x) stringr::str_c(x[1], "_vs_", x[2]))

level = 2
foreach(i = seq_along(tissues)) %do%{
  tissue <- tissues[i]
  foreach(j = seq_along(type2_comb)) %do%{
    comb_name <- type2_comb[j]
    generateProgramLine(script_path, level, tissue, comb_name)
  }
}

## Type 3 ----
ataxia_diagnoses <- c("FRDA", "SCA1", "SCA2", "SCA6")
type3_comb <- sapply(ataxia_diagnoses, function(x) stringr::str_c(x, "_vs_Control"))

level = 3
foreach(i = seq_along(tissues)) %do%{
  tissue <- tissues[i]
  foreach(j = seq_along(type3_comb)) %do%{
    comb_name <- type3_comb[j]
    generateProgramLine(script_path, level, tissue, comb_name)
  }
}

# Make executable ----
system2(command = "chmod", args = (c("+x", script_path)))

# To run the command:
# docker run -it --rm -v /home/grocamora/Core_Projects/Ataxia_Splicing_Analysis/Leafviz:/Leafviz --name grocamora_leafcutter lroser/run_leafcutter bash -c "/Leafviz/prepare_results.sh"