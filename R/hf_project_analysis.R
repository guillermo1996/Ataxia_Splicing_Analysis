#' Analysis of a specific project or study level
#'
#' The final step in the splicing noise analysis requires a specific clustering
#' of the samples. This function takes the metadata of the samples, the
#' field/column to split by clusters and which clusters to combine.
#'
#' It combines the junctions for the clustered samples and measure de
#' mis-splicing ratios. It output two different tables: one with the annotated
#' introns and one with the novel junctions.
#'
#' @param metadata_project Dataframe containing all the metadata.
#' @param project_path Character vector, path to where to store the results.
#' @param level Character vector, field/column in the metadata dataframe to
#'   split the samples in clusters.
#' @param clusters List of character vector, the clusters that are being
#'   studied. Must be named after the field specific in the "level" argument.
#' @param rw_disk Whether to store the results in disk. By default, TRUE.
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#'
#' @return
#' @export
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
    cluster_samples <- metadata_cluster$ID_anon %>% unique
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

#' Extracts the common annotated introns
#'
#' Given a path to DB results from \code{projectAnalysis}, it reads the
#' annotated reference introns and extract the common across all clusters.
#'
#' @param project_path Character vector, path to where to store the results.
#' @param prune Boolean, whether to remove unnecessary columns from the
#'   annotated introns. Defaults to TRUE.
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#' @param output_file Character vector, path to where to store the results.
#'
#' @return dataframe containing the common annotated reference introns.
#' @export
getCommonIntrons <- function(project_path, prune = T, overwrite = F, output_file = "common_introns.rds"){
  if(file.exists(paste0(project_path, output_file)) & !overwrite){
    common_introns <- readRDS(paste0(project_path, output_file))
    return(common_introns)
  }
  
  intron_files <- list.files(project_path, pattern = "_db_introns.rds$")
  global_introns <- foreach(i = seq_along(intron_files)) %do%{
    cluster_file <- intron_files[i]
    cluster_name <- stringr::str_split_fixed(cluster_file, "_", 2)[1]
    cluster_introns <- readRDS(paste0(project_path, cluster_file)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(cluster = cluster_name)
    
    return(cluster_introns)
  } %>% dplyr::bind_rows()
  
  if(prune){
    global_introns <- global_introns %>%
      dplyr::select(ref_junID, ref_type, MSR_Donor, MSR_Acceptor, ref_reads, ref_ss5score, ref_ss3score, cluster)
  }
  
  common_introns <- global_introns %>%
    dplyr::group_by(ref_junID) %>%
    dplyr::filter(n() == length(intron_files)) %>%
    dplyr::ungroup()
  
  if(output_file != ""){
    common_introns %>% saveRDS(paste0(project_path, output_file))
  }
  
  return(common_introns)
}

#' Extracts the common novel junctions
#'
#' Given a path to DB results from \code{projectAnalysis}, it reads the novel
#' junctions and extract the associated to common reference introns.
#'
#' @param project_path Character vector, path to where to store the results.
#' @common_introns Dataframe containing the common annotated reference introns
#'   from \code{getCommonIntrons}.
#' @param prune Boolean, whether to remove unnecessary columns from the
#'   annotated introns. Defaults to TRUE.
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#' @param output_file Character vector, path to where to store the results.
#'
#' @return dataframe containing the novel junctions associated to common
#'   reference introns.
#' @export
getCommonNovel <- function(project_path, 
                           common_introns, 
                           prune = T, 
                           overwrite = F, 
                           output_file = "common_novel.rds"){
  if(file.exists(paste0(project_path, output_file)) & !overwrite){
    common_novel <- readRDS(paste0(project_path, output_file))
    return(common_novel)
  }
  
  novel_files <- list.files(project_path, pattern = "_db_novel.rds$")
  global_novel <- foreach(i = seq_along(novel_files)) %do%{
    cluster_file <- novel_files[i]
    cluster_name <- stringr::str_split_fixed(cluster_file, "_", 2)[1]
    cluster_introns <- readRDS(paste0(project_path, cluster_file)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(cluster = cluster_name)
    
    return(cluster_introns)
  } %>% dplyr::bind_rows()

  
  if(prune){
    global_novel <- global_novel %>%
      dplyr::select(novel_junID, ref_junID, novel_type, novel_reads, distance, novel_ss5score, novel_ss3score, ref_ss5score, ref_ss3score, cluster) %>%
      dplyr::mutate(delta_ss5score = ref_ss5score - novel_ss3score,
                    delta_ss3score = ref_ss3score - novel_ss3score) %>%
      dplyr::select(-ref_ss5score, -ref_ss3score, -novel_ss5score, -novel_ss3score)
  }
  
  common_novel <- global_novel %>%
    dplyr::filter(ref_junID %in% unique(common_introns$ref_junID))
  
  if(output_file != ""){
    common_novel %>% saveRDS(paste0(project_path, output_file))
  }
  
  return(common_novel)
}

#' Mis-splicing ratio Wilcox test analysis
#'
#' Given the dataframe of common annotated introns, it executes the Wilcoxon
#' paired signed-rank test between the clusters' MSR. Results are reported in a
#' dataframe which specifies the splicing site, the p.value of the test, the
#' effect size and its magnitude.
#'
#' @param common_introns Dataframe containing the common annotated reference
#'   introns from \code{getCommonIntrons}.
#' @param project_path Character vector, path to where to store the results.
#' @param clusters List of character vectors, the clusters that are being
#'   studied.
#' @param splice_sites List of character vectors, the splice sites. Defaults to:
#'   c("Donor", "Acceptor).
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#' @param output_file Character vector, path to where to store the results.
#'
#' @return Dataframe containing the results of the Wilcoxon paired signed-rank
#'   test.
#' @export
MSRanalysis <- function(common_introns, 
                        project_path, 
                        clusters, 
                        splice_sites = c("Donor", "Acceptor"), 
                        overwrite = F, 
                        output_file = "wilcox_test_MSR.rds"){
  if(file.exists(paste0(project_path, output_file)) & !overwrite){
    wilcox_test_MSR <- readRDS(paste0(project_path, output_file))
    return(wilcox_test_MSR)
  }
  
  wilcox_test_MSR <- foreach(i = seq(splice_sites), .packages = "tidyverse") %dopar%{
    splice_site <- splice_sites[i]
    
    MSR_Table <- common_introns %>%
      dplyr::select(ref_junID, MSR_Donor, MSR_Acceptor, cluster) %>%
      tidyr::pivot_wider(id_cols = ref_junID, names_from = c("cluster"), values_from = c(paste0("MSR_", splice_site)))
    
    MSR_clusterA <- MSR_Table %>% arrange(ref_junID) %>% pull(clusters[1])
    MSR_clusterB <- MSR_Table %>% arrange(ref_junID) %>% pull(clusters[2])
    
    wilcox_test <- coin::wilcoxsign_test(MSR_clusterA ~ MSR_clusterB, paired = TRUE)
    wilcox_coin_pvalue <- show(wilcox_test)$p.value %>% as.numeric()
    wilcox_Z <- wilcox_test@statistic@teststatistic
    wilcox_effsize <- wilcox_Z/sqrt(length(MSR_clusterA))
    wilcox_magnitude <- rstatix:::get_wilcox_effsize_magnitude(wilcox_effsize)
    
    tibble(splice_site = splice_site,
           p.value = wilcox_coin_pvalue,
           effect_size = wilcox_effsize,
           magnitude = wilcox_magnitude)
  } %>% dplyr::bind_rows()
  
  if(output_file != ""){
    wilcox_test_MSR %>% saveRDS(paste0(project_path, output_file))
  }
  
  return(wilcox_test_MSR)
}
