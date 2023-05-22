subsampleGowerDistance <- function(metadata_project,
                                   level,
                                   clusters,
                                   ref_cluster = NULL,
                                   weights = NULL,
                                   n = 1){
  # Default weights for "RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex". RIN
  # is usually the only relevant one.
  if(is.null(weights)) weights <- c(1, 0, 0, 0, 0)
  
  # Extract the cluster with the least number of samples or the reference
  # cluster
  if(!is.null(ref_cluster)){
    ref_metadata_cluster = metadata_project %>%
      dplyr::filter(!!sym(level) == ref_cluster)
  }else{
    ref_metadata_cluster = metadata_project %>% 
      dplyr::filter(!!sym(level) %in% clusters) %>%
      dplyr::group_by(!!sym(level)) %>% 
      dplyr::mutate(n = n()) %>%
      dplyr::ungroup() %>%
      dplyr::slice_min(n)
    
    ref_cluster <- ref_metadata_cluster %>% 
      dplyr::pull(!!sym(level)) %>%
      unique()
  }
  
  # For every other cluster, extract the samples with the lowest Gower distance
  metadata_subsample <- foreach(i = seq_along(clusters)) %do%{
    cluster <- clusters[i]
    if(cluster == ref_cluster) return()
    
    cluster_samples <- metadata_project %>%
      dplyr::filter(!!sym(level) == cluster)
    
    # Loop for the number of samples to extract per reference sample
    cluster_subsamples <- foreach(j = 1:n) %do%{
      # Loop thorugh all the reference samples
      foreach(k = 1:nrow(ref_metadata_cluster)) %do%{
        # Prepare the metadata to measure the Gower Distance
        ref_sample <- ref_metadata_cluster[k, ] %>% prepareMetadata()
        cluster_samples_df <- cluster_samples %>% prepareMetadata()
        
        cluster_id <- StatMatch::gower.dist(ref_sample, cluster_samples_df, var.weights = weights) %>%
          `rownames<-`(rownames(ref_sample)) %>%
          `colnames<-`(rownames(cluster_samples_df)) %>% .[1, ] %>% which.min()
        
        # Extract the matched sample
        matched_sample <- cluster_samples[cluster_id, ]
        cluster_samples <- cluster_samples[-cluster_id, ]
        
        return(matched_sample)
      } %>% dplyr::bind_rows()
    } %>% dplyr::bind_rows()
  } %>% append(list(ref_metadata_cluster)) %>% dplyr::bind_rows() %>% dplyr::ungroup()
  
  return(metadata_subsample)
}

prepareMetadata <- function(metadata, 
                            columns = c("Individual_ID", "RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex")){
  metadata %>% 
    dplyr::select(all_of(columns)) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Individual_ID")
}