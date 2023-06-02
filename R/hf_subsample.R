#' Subsamples by Gower distance with weights
#'
#' Given the metadata of the samples, it subsamples the majority class to match
#' in number of samples the minority class. The subsampling is done by matching
#' samples together in order to minimize the Gower distance between the two
#' clusters.
#'
#' @param metadata_project Dataframe containing all the metadata.
#' @param level Character, field/column in the metadata dataframe to split the
#'   samples in clusters.
#' @param id_field Character, name of the metadata field that will be used as
#'   rownames for the procedure. Must be unique to each sample.
#' @param covariates Character vector, covariates that will be used to get the
#'   percentage of variance of explained by each one. Must much the naming of
#'   columns from the metadata argument, and keep a consistent ordering with the
#'   "covariates" argument in the \code{getVarianceDf} function.
#' @param clusters Character vector, the clusters that are being studied. Must
#'   be named after the field specific in the "level" argument.
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
                                   id_field,
                                   covariates,
                                   clusters,
                                   ref_cluster = NULL,
                                   weights = NULL,
                                   n = 1){
  # Return if input metadata is already balanced
  unique_counts <- metadata_project %>% 
    dplyr::filter(!!sym(level) %in% clusters) %>% 
    dplyr::pull(!!sym(level)) %>%
    table() %>%
    unique() %>%
    length()
  if(unique_counts == 1){
    logger::log_info("\t Same number of samples in both classes. No subsample needed.")
    return(metadata_project)
  }
  
  # Default weights for "RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex". RIN
  # is usually the only relevant one.
  if(any(class(weights) == "data.frame")){
    weights <- weights %>% deframe %>% .[covariates] %>% unname
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
  
  ## ######### Subsample technique ############
  ##
  ## The subsample steps are:
  ##
  ## 1. Loop through every cluster. Ignore the reference cluster (cluster with
  ## the minimum number of samples).
  ##
  ## 2. Get the cluster metadata. Prepare the metadata to apply the Gower
  ## distance function by calling "prepareMetadata()". It is important to modify
  ## that function so that the order of covariates is the appropriate. In this
  ## case, the order is: RIN, PMI, Brain.Bank, Age_of_death and Sex. The column
  ## "Individual_ID" is employed as the rownames.
  ##
  ## 3. Loop for the variable "n", which represents how many samples from the
  ## majority class we want for each of the minority class. Defaults to n = 1.
  ##
  ## 4. Get the reference cluster metadata and execute the "prepareMetadata()"
  ## function too.
  ##
  ## 5. Loop through every sample in the reference cluster.
  ##
  ## 6. In each iteration, measure the weighted Gower distance between all
  ## reference samples and cluster samples. Pair the samples with the minimum
  ## distance between them and remove both of them from their respective
  ## metadata dataframe (so that they are not selected again in the next
  ## iteration). After the loop is executed, we should have a matched samples
  ## for every reference sample. Since we are looping through "n" (step 3), we
  ## will repeat the process, matching a new set of cluster samples for every
  ## reference sample without repeating the previously matched.
  ##
  ## 7. Each loop returns the metadata of the paired samples. In the end, we
  ## should have a dataframe with the metadata of all samples and their match.
  ## If "n" is different from 1, the reference samples will be repeated. We
  ## apply a "dplyr::distinct()" to remove them.
  
  ## Step 1. Loop through every cluster.
  metadata_subsample <- foreach(i = seq_along(clusters)) %do%{
    cluster <- clusters[i]
    if(cluster == ref_cluster) return()
    
    ## Step 2. Create the cluster metadata and prepare the dataframe.
    cluster_metadata_samples <- metadata_project %>% dplyr::filter(!!sym(level) == cluster)
    cluster_samples <- cluster_metadata_samples %>% prepareMetadata(id_field, covariates)
    
    ## Step 3. Loop through "n".
    foreach(j = 1:n) %do%{
      ## Step 4. Create the reference metadata and prepare the dataframe.
      ref_samples <- ref_metadata_cluster %>% prepareMetadata(id_field, covariates)
      
      ## Step 5. Loop through every reference sample.
      foreach(k = 1:nrow(ref_metadata_cluster)) %do%{
        
        ## Step 6. Extract the distance between all reference and cluster
        ## samples. Remove the pairs from their respective dataframes.
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


#' Prepares the metadata for the Gower distance calculation.
#'
#' The first element in the argument "columns" will be used as rownames of the
#' returned dataframe.
#'
#' @param metadata Dataframe containing all the metadata.
#' @param columns List of character vectors, columns and order to use for the
#'   Gower distance calculation.
#'
#' @return dataframe with the metadata preparared for Gower distance
#'   calculation.
#' @export
prepareMetadata <- function(metadata, 
                            id_field = "Individual_ID",
                            covariates = c("RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex")){
  columns <- c(id_field, covariates)
  metadata %>% 
    dplyr::select(all_of(columns)) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = id_field)
}


#' Generates junction information for every sample
#'
#' To estimate the weights in the subsampling by Gower distance minimization, we
#' require some junction information for every sample. This functions measures:
#'
#' 1. Number of novel junctions.
#'
#' 2. Number of total unique junctions (reference + novel).
#'
#' 3. Number of reads associated to all junctions.
#'
#' @param results_path Path to where the results files of the splicing analysis
#'   are found.
#' @param output_file Path to where to store the results of the function. Leave
#'   empty to not store in disk. If provided and the file exists, results will
#'   be loaded from disk.
#'
#' @return dataframe with the information extracted.
#' @export
generateSampleJunctionInformation <- function(results_path, 
                                              output_file = ""){
  # If fule already exists, load it from disk.
  if(file.exists(output_file)){
    junction_information <- readRDS(output_file)
    return(junction_information)
  }
  
  # Read variables generated during the analysis.
  all_distances_pruned <- readRDS(file.path(results_path, "all_distances_pruned.rds"))
  all_reads_combined <- readRDS(file.path(results_path, "all_reads_combined.rds"))
  
  # Number of novel junctions by sample
  number_novel_junctions <- all_distances_pruned %>% 
    dplyr::group_by(sample) %>% 
    dplyr::count(name = "number_novel_junctions")
  
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


#' Extracts the percentage of variance explained by each covariate
#'
#' Given a metadata dataframe and both the response variable and the covariates,
#' it measures the percentage of variance explained by each covariate using
#' "anova".
#'
#' @param metadata_extra Dataframe containing the relevant metadata to each
#'   sample. Must contain both the covariates and the response variable.
#' @param covariates Character vector, covariates that will be used to get the
#'   percentage of variance of explained by each one. Must much the naming of
#'   columns from the metadata argument, and keep a consistent ordering with the
#'   "covariates" argument in the \code{subsampleGowerDistance} function.
#' @param response_var Character, name of the response variable to employ in the
#'   weights calculation. Defaults to the total number of mapped reads
#'   associated to the studied junctions ("mapped_junctions").
#'
#' @return
#' @export
#'
#' @examples
getVarianceExplained <- function(metadata_extra, 
                                 covariates = c("RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex"),
                                 response_var = "mapped_junctions"){
  metadata_cor <- metadata_extra %>%
    dplyr::select(all_of(c(response_var, covariates))) %>%
    dplyr::mutate(across(where(is.character), ~ as.factor(.)))
  
  fm <- as.formula(paste(response_var, "~ ."))
  fit <- lm(fm, data = metadata_cor)
  af <- anova(fit)[1:length(covariates), ]
  afss <- af$"Sum Sq"
  variance_df <- cbind(af, Variance_explained = afss/sum(afss)) %>% 
    tibble::rownames_to_column("Predictor") %>% 
    tibble::as_tibble() %>%
    dplyr::select(Predictor, Variance_explained)
  
  return(variance_df)
}


#' Calculates the percentage of variance explained by each covariate
#'
#' To measure the Gower distance between different samples, we apply weights to
#' each of the variables. These weights are obtained based of the percentage of
#' variance explained by each variable to the response variable of our choice.
#'
#' This function is a combination of \code{generateSampleJunctionInformation}
#' and \code{getVarianceExplained}. More information about the specifics of each
#' function in their respective documentation.
#'
#' @param metadata_project Dataframe containing the relevant metadata to each
#'   sample.
#' @param results_path Path to where the results files of the splicing analysis
#'   are found.
#' @param covariates Character vector, covariates that will be used to get the
#'   percentage of variance of explained by each one. Must much the naming of
#'   columns from the metadata argument, and keep a consistent ordering with the
#'   "covariates" argument in the \code{subsampleGowerDistance} function.
#' @param response_var Character, name of the response variable to employ in the
#'   weights calculation. Defaults to the total number of mapped reads
#'   associated to the studied junctions ("mapped_junctions").
#'
#' @return Dataframe containing the percentage of variance explained by each
#'   covariate.
#' @export
getVarianceDf <- function(metadata_project,
                          results_path,
                          covariates = c("RIN", "PMI", "Brain.Bank", "Age_at_death", "Sex"),
                          response_var = "mapped_junctions"){
  # Extract relevant junction information for each sample that will be used as
  # the response variable of the linear model.
  junction_information <- generateSampleJunctionInformation(results_path = results_path,
                                                            output_file = here::here("variables/sample_junction_information.rds"))
  
  # Combine the junction information to the current metadata. Take care of the
  # columns to combine with.
  metadata_extra <- metadata_project %>%
    dplyr::left_join(junction_information, by = c("Correct_sample" = "sample"))
  
  variance_df <- getVarianceExplained(metadata_extra = metadata_extra,
                                      covariates = covariates,
                                      response_var = response_var) 
  
  return(variance_df)
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