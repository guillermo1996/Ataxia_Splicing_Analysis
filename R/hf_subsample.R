#' Subsamples by Gower distance with weights
#'
#' Given the metadata of the samples, it subsample the majority class to match
#' in number of samples the minority class. The subsampling is done by matching
#' samples together in order to minimize the Gower distance between the two
#' clusters.
#'
#' @param metadata_project Dataframe containing all the metadata.
#' @param level Character vector, field/column in the metadata dataframe to
#'   split the samples in clusters.
#' @param clusters List of character vector, the clusters that are being
#'   studied. Must be named after the field specific in the "level" argument.
#' @param ref_cluster Character vector, cluster that will be used as reference
#'   to match samples against. If not provided, defaults to the minority class.
#' @param weights Numeric vector, weights associated to each covariate. By
#'   default, only RIN is accounted for. Weights must be specific in order:
#'   c(RIN, PMI, BrainBank, AgeOfDeath, Sex). It can also accepts a dataframe of
#'   percentage of covariance explained as the one output by
#'   \code{getVarianceDf}.
#' @param n Numeric, number of samples of the majority class to use per samples
#'   in the minority class. Recommended to leave as default. Defaults to 1.
#'
#' @return dataframe containing the metadata of the samples selected.
#' @export
subsampleGowerDistance <- function(metadata_project,
                                       level,
                                       clusters,
                                       ref_cluster = NULL,
                                       weights = NULL,
                                       n = 1){
  # Default weights for "RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex". RIN
  # is usually the only relevant one.
  if(any(class(weights) == "data.frame")){
    weights <- weights %>% deframe %>% .[c("RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex")] %>% unname
  }else  if(is.null(weights)){
    weights <- c(1, 0, 0, 0, 0)
  }
  
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
      dplyr::slice_min(n) %>%
      dplyr::select(-n)
    
    ref_cluster <- ref_metadata_cluster %>% 
      dplyr::pull(!!sym(level)) %>%
      unique()
  }
  
  metadata_subsample <- foreach(i = seq_along(clusters)) %do%{
    cluster <- clusters[i]
    if(cluster == ref_cluster) return()
    
    cluster_metadata_samples <- metadata_project %>% dplyr::filter(!!sym(level) == cluster)
    cluster_samples <- cluster_metadata_samples %>% prepareMetadata()
    
    foreach(j = 1:n) %do%{
      ref_samples <- ref_metadata_cluster %>% prepareMetadata()
      
      foreach(k = 1:nrow(ref_metadata_cluster)) %do%{
        distances_df <- StatMatch::gower.dist(ref_samples, cluster_samples, var.weights = weights) %>%
          `rownames<-`(rownames(ref_samples)) %>%
          `colnames<-`(rownames(cluster_samples))
        
        min_ref <- names(which.min(apply(distances_df, 1, min)))
        min_cluster <- names(which.min(apply(distances_df, 2, min)))
        
        matched_samples <- rbind(ref_metadata_cluster[ref_metadata_cluster$Individual_ID == min_ref, ], 
                                 cluster_metadata_samples[cluster_metadata_samples$Individual_ID == min_cluster, ])
        ref_samples <- ref_samples[!(rownames(ref_samples) == min_ref), ]
        cluster_samples <- cluster_samples[!(rownames(cluster_samples) == min_cluster), ]
        
        return(matched_samples)
      } %>% dplyr::bind_rows()
    } %>% dplyr::bind_rows()
  } %>% dplyr::bind_rows() %>% dplyr::distinct()
  
  return(metadata_subsample)
}

# subsampleGowerDistance_old <- function(metadata_project,
#                                    level,
#                                    clusters,
#                                    ref_cluster = NULL,
#                                    weights = NULL,
#                                    n = 1){
#   # Default weights for "RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex". RIN
#   # is usually the only relevant one.
#   if(any(class(weights) == "data.frame")){
#     weights <- weights %>% deframe %>% .[c("RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex")] %>% unname
#   }else  if(is.null(weights)){
#     weights <- c(1, 0, 0, 0, 0)
#   }
#   
#   # Extract the cluster with the least number of samples or the reference
#   # cluster
#   if(!is.null(ref_cluster)){
#     ref_metadata_cluster = metadata_project %>%
#       dplyr::filter(!!sym(level) == ref_cluster)
#   }else{
#     ref_metadata_cluster = metadata_project %>% 
#       dplyr::filter(!!sym(level) %in% clusters) %>%
#       dplyr::group_by(!!sym(level)) %>% 
#       dplyr::mutate(n = n()) %>%
#       dplyr::ungroup() %>%
#       dplyr::slice_min(n)
#     
#     ref_cluster <- ref_metadata_cluster %>% 
#       dplyr::pull(!!sym(level)) %>%
#       unique()
#   }
#   
#   # For every other cluster, extract the samples with the lowest Gower distance
#   metadata_subsample <- foreach(i = seq_along(clusters)) %do%{
#     cluster <- clusters[i]
#     if(cluster == ref_cluster) return()
#     
#     cluster_samples <- metadata_project %>%
#       dplyr::filter(!!sym(level) == cluster)
#     
#     # Loop for the number of samples to extract per reference sample
#     cluster_subsamples <- foreach(j = 1:n) %do%{
#       # Loop thorugh all the reference samples
#       foreach(k = 1:nrow(ref_metadata_cluster)) %do%{
#         # Prepare the metadata to measure the Gower Distance
#         ref_sample <- ref_metadata_cluster[k, ] %>% prepareMetadata()
#         cluster_samples_df <- cluster_samples %>% prepareMetadata()
#         
#         cluster_id <- StatMatch::gower.dist(ref_sample, cluster_samples_df, var.weights = weights) %>%
#           `rownames<-`(rownames(ref_sample)) %>%
#           `colnames<-`(rownames(cluster_samples_df)) %>% .[1, ] %>% which.min()
#         
#         # Extract the matched sample
#         matched_sample <- cluster_samples[cluster_id, ]
#         cluster_samples <- cluster_samples[-cluster_id, ]
#         
#         return(matched_sample)
#       } %>% dplyr::bind_rows()
#     } %>% dplyr::bind_rows()
#   } %>% append(list(ref_metadata_cluster)) %>% dplyr::bind_rows() %>% dplyr::ungroup()
#   
#   return(metadata_subsample)
# }

#' Prepares the metadata for the Gower distance calculation.
#'
#' @param metadata Dataframe containing all the metadata.
#' @param columns List of character vectors, columns and order to use for the
#'   Gower distance calculation.
#'
#' @return dataframe with the metadata preparared for Gower distance
#'   calculation.
#' @export
prepareMetadata <- function(metadata, 
                            columns = c("Individual_ID", "RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex")){
  metadata %>% 
    dplyr::select(all_of(columns)) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Individual_ID")
}

generateSampleJunctionInformation <- function(results_path, output_file = ""){
  if(file.exists(output_file)){
    junction_information <- readRDS(output_file)
    return(junction_information)
  }
  all_distances_pruned <- readRDS(file.path(results_path, "all_distances_pruned.rds"))
  all_reads_combined <- readRDS(file.path(results_path, "all_reads_combined.rds"))
  
  # Number of novel junctions by sample
  number_novel_junctions <- all_distances_pruned %>% dplyr::group_by(sample) %>% dplyr::count(name = "number_novel_junctions")
  
  # Number of junctions
  junction_reads <- all_reads_combined %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand) %>%
    tidyr::pivot_longer(cols = -junID, names_to = "sample", values_to = "reads", values_drop_na = TRUE) %>%
    dplyr::group_by(sample)
  ## Unique
  number_junctions <- junction_reads %>% dplyr::count(name = "number_junctions")
  ## Reads
  number_junction_reads <- junction_reads %>% dplyr::summarise(mapped_junctions = sum(reads))
  
  # Return combination
  junction_information <- number_novel_junctions %>%
    dplyr::left_join(number_junctions, by = "sample") %>%
    dplyr::left_join(number_junction_reads, by = "sample")
  
  if(output_file != ""){
    junction_information %>% saveRDS(output_file)
  }
  
  return(junction_information)
}

getVarianceExplained <- function(metadata_extra, y_var = "mapped_junctions"){
  metadata_cor <- metadata_extra %>%
    dplyr::select(!!sym(y_var), RIN, PMI, Brain.Bank, Age_at_death, Sex) %>%
    dplyr::mutate(Brain.Bank = ifelse(Brain.Bank == "QSBB", 1, 0),
                  Sex = ifelse(Sex == "M", 1, 0))
  
  fm <- as.formula(paste(y_var, "~ ."))
  fit <- lm(fm, data = metadata_cor)
  af <- anova(fit)[1:5, ]
  afss <- af$"Sum Sq"
  variance_df <- cbind(af, Variance_explained = afss/sum(afss)) %>% 
    tibble::rownames_to_column("Predictor") %>% 
    tibble::as_tibble() %>%
    dplyr::select(Predictor, Variance_explained)
  
  return(variance_df)
}

getVarianceDf <- function(metadata_project, results_path, output_file){
  junction_information <- generateSampleJunctionInformation(results_path = results_path,
                                                            output_file = output_file)
  
  metadata_extra <- metadata_project %>%
    dplyr::left_join(junction_information, by = c("Correct_sample" = "sample"))
  
  variance_df <- getVarianceExplained(metadata_extra)
  
  return(variance_df)
}
