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