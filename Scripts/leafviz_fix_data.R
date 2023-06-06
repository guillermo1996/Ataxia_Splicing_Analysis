## _________________________________________________
##
## Fix for Leafviz .RData
##
## Aim: fix cluster strand in .RData inputs for Leafviz.
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
## The .RData files created with the script "leafviz_generate_data.R" need a
## further fixing because of the strand information for clusters. It might be
## possible to fix it using the "prepare_resuts.R" command itself, but did not
## manage to find a way.
##
## The script removes all "_NA" text from the cluster IDs and saves the .RData
## in a different location.
##
## Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

leafviz_main_path <- "/home/grocamora/Core_Projects/Ataxia_Splicing_Analysis/Leafviz/"

broken_data_path <- file.path(leafviz_main_path, "Leafviz_data")
fixed_data_path <- file.path(leafviz_main_path, "Leafviz_data_fixed")

RData_files <- list.files(broken_data_path)

for(RData_file in RData_files){
  load(file.path(broken_data_path, RData_file))
  
  introns$clusterID <- gsub("_NA", "", introns$clusterID)
  rownames(counts) <- gsub("_NA", "", rownames(counts))
  clusters$clusterID <- gsub("_NA", "", clusters$clusterID)
  
  save(cluster_summary, clusters, counts, exons_table, intron_summary, introns, 
       introns_to_plot, meta, pca, sample_table, annotation_code, cluster_ids, code, 
       file = file.path(fixed_data_path, RData_file), compression_level = 9)
}
