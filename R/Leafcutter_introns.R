#' Removes the strand from the cluster ID
#'
#' @param df data.frame containing a cluster ID in the "input_name" column.
#' @param use_strand boolean to specify whether to use the strand in the
#'   clusters' ID.
#' @param input_name character vector with the name of the input cluster ID.
#' @param output_name character vector with the name of the output cluster ID.
#'
#' @return data.frame with the renamed cluster ID
#' @export
removeClusterStrand <- function(df, use_strand, input_name = "cluster"){
  if(!use_strand){
    df <- df %>%
      tidyr::separate(col = !!sym(input_name), into = c("prefix", "cluster_id", "strand"), sep = "_", remove = T) %>%
      tidyr::unite(col = "cluster_id", c("prefix", "cluster_id")) %>%
      dplyr::select(-strand)
  }else{
    df <- df %>%
      dplyr::rename(cluster_id = !!sym(input_name))
  }
  
  return(df)
}

#' Extract the junctions from significant comparisons and clusters
#'
#' @param annotated_gr GRanges object with annotated junctions.
#' @param leafcutter_result leafcutter results of the study.
#' @param use_strand boolean to specify whether to use the strand in the
#'   clusters' ID.
#' @param tissue character vector with the specific tissue to study. The tissue
#'   name must be contained in the comparison name.
#'
#' @return data.frame with significant junctions
#' @export
extractSignificantIntrons <- function(annotated_gr, leafcutter_result, use_strand = F, tissue = NULL){
  # Get the significant comparisons and clusters
  significant_cluster <- leafcutter_result$significant_clusters_0.05_filter %>%
    # If tissue is provided, filter the comparisons
    {if(!is.null(tissue)) dplyr::filter(., grepl(tissue, comparison)) else .} %>%
    dplyr::mutate(cluster_id = str_c(chr, ":", cluster)) %>%
    removeClusterStrand(use_strand, input_name = "cluster_id") %>%  
    dplyr::filter(abs(deltapsi) >= 0.1) %>% 
    dplyr::distinct(comparison, cluster_id)
  
  # Extract the junctions with significant comparison and cluster
  significant_annotated <- annotated_gr %>%
    tibble::as_tibble() %>%
    dplyr::inner_join(significant_cluster)
  
  return(significant_annotated)
}

#' Convert leafcutter introns to be annotated
#'
#' Inputs a leafcutter splicing analysis and converts it to a GRanges object
#' that can be annotated.
#' 
#' Source: https://github.com/RHReynolds/LBD-seq-bulk-analyses/blob/main/R/convert_leafcutter.R
#'
#' @param leafcutter_result leafcutter results of the study.
#' @param use_strand boolean to specify whether to use the strand in the
#'   clusters' ID.
#'
#' @return GRanges object from leafcutter results.
#' @export
convert_leafcutter <- function(leafcutter_results, use_strand = F){
  # Split the intron and cluster_id columns
  leafcutter_results <- leafcutter_results %>% 
    tidyr::separate(col = "intron", into = c("chr", "intron_start", "intron_end", "cluster_id"), sep = ":") %>%
    dplyr::mutate(cluster_id = str_c(chr, ":", cluster_id)) %>%
    tidyr::separate(col = "cluster_id", into = c("prefix", "cluster_id", "strand"), sep = "_", remove = F) %>%
    dplyr::mutate(chr = str_replace(chr, "chr", ""))
  
  # Whether to include the strand in the GRanges object
  if(use_strand){
    leafcutter_gr <- leafcutter_results %>% 
      tidyr::unite(col = "cluster_id", c("prefix", "cluster_id", "strand"), sep = "_", remove = F) %>%
      GenomicRanges::makeGRangesFromDataFrame(start.field = "start",
                                              end.field = "end",
                                              seqnames.field = "chr",
                                              strand.field = "strand",
                                              keep.extra.columns = T)
  }else{
    leafcutter_gr <- leafcutter_results %>% 
      dplyr::select(-strand) %>%
      tidyr::unite(col = "cluster_id", c("prefix", "cluster_id"), sep = "_", remove = T) %>%
      GenomicRanges::makeGRangesFromDataFrame(start.field = "start",
                                              end.field = "end",
                                              seqnames.field = "chr",
                                              keep.extra.columns = T)
  }
  
  # Leafcutter intron definition adds an extra +1 to intron ends. For best
  # matching to ref. annotation must remove 1 bp
  GenomicRanges::end(leafcutter_gr) <- GenomicRanges::end(leafcutter_gr) - 1
  return(leafcutter_gr)
}

#' Remove the junctions from ambiguous genes.
#'
#' @param input_SR_details Tibble (or data.frame) object with the relevant
#'   junctions.
#'
#' @return Junctions assigned to only one gene.
#' @export
removeAmbiguousGenes <- function(input_SR_details){
  output_SR_details <- input_SR_details[which(sapply(input_SR_details$gene_id_junction, length) < 2), ]
  
  return(output_SR_details)
}

#' Prints the proportion of significant junctions by type
#'
#' Given the list of annotated junctions and the list of leafcutter results, this
#' function extract a given index of the lists (corresponding to different
#' studies) and prints a table with the proportion of significantly
#' differentially spliced junction by type.
#'
#' @param annotated list of annotated junctions. Each element corresponds to a specific study.
#' @param leafcutter_list list of leafcutter results. Each element corresponds to a specific study.
#' @param level index or name of the specific study in both the "annotated" and "leafcutter_list" arguments.
#' @param use_strand boolean to specify whether to use the strand in the clusters' ID.
#' @param deltapsi_filter boolean to specify whether to use a deltaPSI filter of |dPSI| >= 0.1
#'
#' @export
printTypeTable <- function(annotated, leafcutter_list, level, use_strand = F, deltapsi_filter = F){
  # Extract the significant comparisons and clusters
  significant_clusters <- leafcutter_list[[level]]$significant_clusters_0.05_filter %>%
    dplyr::mutate(cluster_id = str_c(chr, ":", cluster)) %>% 
    removeClusterStrand(use_strand, input_name = "cluster_id") %>%
    {if(deltapsi_filter) dplyr::filter(., abs(deltapsi) >= 0.1) else .} %>%
    dplyr::distinct(comparison, cluster_id)
  
  # Conver the annotated GRanges object to a tibble.
  annotated_df <- annotated[[level]] %>%
    tibble::as_tibble()
  
  # Extract the number of junctions by type that are significant.
  significant_counts <- annotated_df %>%
    dplyr::inner_join(significant_clusters) %>%
    dplyr::group_by(type) %>% # (Optional) Add comparison/tissue here.
    dplyr::summarise(significant = n())
  
  # Prints the proportion of significant junctions by type.
  annotated_df %>%
    dplyr::group_by(type) %>% # (Optional) Add comparison/tissue here.
    dplyr::summarise(all = n()) %>% 
    dplyr::inner_join(significant_counts) %>%
    dplyr::mutate(proportion = (significant/all) %>% signif(2)) %>% 
    dplyr::arrange(-all) %>% 
    `colnames<-`(c("Junction Type", "Count Junctions", "Significant Junctions", "Proportion")) %>%
    kableExtra::kbl(booktabs = T, linesep = "") %>%
    kableExtra::kable_classic(full_width = T, "hover", "striped", html_font = "Cambria", font_size = 14) %>%
    kableExtra::row_spec(0, bold = T, font_size = 16)
}

#' Plots the proportion of junction categories
#'
#' Given the list of annotated junctions and the list of leafcutter results,
#' this function extract a given index of the lists (corresponding to different
#' studies) and plots the proportion of each category.
#'
#' Additional parameters include the option to use the cluster strand or to
#' divide by tissue.
#'
#' @param annotated list of annotated junctions. Each element corresponds to a
#'   specific study.
#' @param leafcutter_list list of leafcutter results. Each element corresponds
#'   to a specific study.
#' @param level index or name of the specific study in both the "annotated" and "leafcutter_list" arguments.
#' @param use_strand boolean to specify whether to use the strand in the
#'   clusters' ID.
#' @param tissue character vector with the specific tissue to study. The tissue
#'   name must be contained in the comparison name.
#'
#' @return
plotProportionOfAnnotation <- function(annotated, leafcutter_list, level, use_strand = F, tissue = NULL){
  # Get the successful comparisons and clusters
  success_cluster <- leafcutter_list[[level]]$cluster_significance %>%
    # If tissue is provided, it filter the comparisons by tissue
    {if(!is.null(tissue)) dplyr::filter(., grepl(tissue, comparison)) else .} %>%
    dplyr::filter(status == "Success") %>%
    dplyr::distinct(comparison, cluster) %>%
    removeClusterStrand(use_strand)
  
  # Extract the introns with successful comparison and cluster
  success <- annotated[[level]] %>%
    tibble::as_tibble() %>%
    dplyr::inner_join(success_cluster)
  
  # Get introns with significant comparison and cluster
  significant <- extractSignificantIntrons(annotated[[level]], leafcutter_list[[level]], use_strand, tissue)
  
  # Plots will be stored in a list 
  intron_list <- list(success, significant)
  plot_list <- vector(mode = "list", 2)
  
  for(j in 1:length(intron_list)){
    # Calculate the proportion of each category by comparison
    count <- intron_list[[j]] %>% 
      dplyr::filter(!type == "ambig_gene") %>% 
      dplyr::group_by(comparison, type) %>% 
      dplyr::summarise(n = n()) %>% 
      dplyr::mutate(prop = n/sum(n)) %>% 
      dplyr::ungroup()
    
    plot_list[[j]] <- count %>% 
      dplyr::mutate(type = type %>% factor(levels = intron_type_levels %>% rev()),
                    #comparison = comparison_labels[comparison],
                    comparison = as_factor(comparison)) %>% 
      ggplot(aes(x = comparison, y = prop, fill = type), colour = "black") +
      geom_col(color = "black") +
      labs(x = "", y = "Proportion") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
      scale_x_discrete(expand = expansion(add = 0.7)) + 
      scale_fill_manual(name = "Acceptor/donor annotation", 
                        values = rev(c("#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#7E6148", "grey"))) +
      custom_gg_theme + 
      guides(fill = guide_legend(reverse = T))
  }
  
  # Arrange both plots in two columns
  ggpubr::ggarrange(plotlist = plot_list,
                    labels = c("a", "b"),
                    common.legend = TRUE, legend = "top") %>%
    ## If tissue is provided, add a title (remove next line to ignore title)
    {if(!is.null(tissue)) ggpubr::annotate_figure(., top = ggpubr::text_grob(paste0("Tissue: ", tissue), size = 16)) else .}
}

#' Plot the number of overlapping significant junctions
#'
#' First, we obtain the annotated introns from significant clusters and
#' comparisons. Then, we remove the junction category "ambig_gene" and count the
#' number of entries by intron ID and category.
#'
#' @param annotated list of annotated junctions. Each element corresponds to a
#'   specific study.
#' @param leafcutter_list list of leafcutter results. Each element corresponds
#'   to a specific study.
#' @param level index or name of the specific study in both the "annotated" and "leafcutter_list" arguments.
#' @param use_strand boolean to specify whether to use the strand in the
#'   clusters' ID.
#' @param tissue character vector with the specific tissue to study. The tissue
#'   name must be contained in the comparison name.
#'
#' @return
#' @export
plotNumberOverlappingIntrons <- function(annotated, leafcutter_list, level, use_strand = T, tissue = NULL){
  # Get introns with significant comparison and cluster
  significant <- extractSignificantIntrons(annotated[[level]], leafcutter_list[[level]], use_strand, tissue)
  
  # Find number of overlapped introns between the comparisons
  number_ovelaps <- significant %>% 
    dplyr::filter(!type == "ambig_gene") %>% 
    dplyr::mutate(unique_id = str_c(cluster_id, ":", start, ":", end),
                  type = type %>% factor(levels = intron_type_levels)) %>% 
    dplyr::group_by(unique_id, type) %>% 
    dplyr::mutate(n_overlaps = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n_overlaps = factor(n_overlaps, levels = 1:length(unique(comparison)))) %>%
    dplyr::distinct(unique_id, type, n_overlaps)
  
  # Plot the results
  ggplot(number_ovelaps, aes(x = n_overlaps, fill = type)) +
    geom_bar(position = position_dodge2(preserve = "single", padding = 0), color = "black") +
    facet_wrap(~ type, labeller = labeller(type = intron_type_labels)) +
    labs(x = "Number of overlaps") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_fill_manual(name = "Acceptor/donor annotation", 
                      values = c("#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#7E6148", "grey"),
                      labels = intron_type_labels) +
    custom_gg_theme +
    theme(axis.text.x = ggplot2::element_text(color = "black", size = 8, angle = 0, hjust = 0.5)) +
    ## If tissue is provided, add a title (remove next line to ignore title)
    {if(!is.null(tissue)) ggtitle(paste0("Tissue: ", tissue))}
}

#' Plot the proportion of overlapping significant junctions by comparison
#'
#' First, we obtain the annotated introns from significant clusters and
#' comparisons. Then, we remove the junction category "ambig_gene" and count the
#' number of entries by intron ID and category. The proportion of each category
#' is calculated for each number of overlaps.
#'
#' @param annotated list of annotated junctions. Each element corresponds to a
#'   specific study.
#' @param leafcutter_list list of leafcutter results. Each element corresponds
#'   to a specific study.
#' @param level index or name of the specific study in both the "annotated" and "leafcutter_list" arguments.
#' @param use_strand boolean to specify whether to use the strand in the
#'   clusters' ID.
#' @param tissue character vector with the specific tissue to study. The tissue
#'   name must be contained in the comparison name.
#'
#' @return
#' @export
plotDistributionAnnotatedTypes <- function(annotated, leafcutter_list, level, use_strand = F, tissue = NULL){
  # Get introns with significant comparison and cluster
  significant <- extractSignificantIntrons(annotated[[level]], leafcutter_list[[level]], use_strand, tissue)
  
  # Find number of overlapped introns between the comparisons
  number_ovelaps <- significant %>% 
    dplyr::filter(!type == "ambig_gene") %>% 
    dplyr::mutate(unique_id = str_c(cluster_id, ":", start, ":", end),
                  type = type %>% factor(levels = intron_type_levels)) %>% 
    dplyr::group_by(unique_id, type) %>% 
    dplyr::mutate(n_overlaps = n() %>% as_factor()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n_overlaps = factor(n_overlaps, levels = 1:length(unique(comparison))))
  
  # Find the proportion of each category for a given number of overlaps
  number_overlaps_proportion <- number_ovelaps %>% 
    dplyr::distinct(comparison, unique_id, n_overlaps, type) %>% 
    dplyr::group_by(comparison, n_overlaps, type) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::group_by(comparison) %>% 
    dplyr::mutate(prop = n/sum(n))
  
  # Plot the results
  ggplot(number_overlaps_proportion,
         aes(x = n_overlaps, y = prop, fill = type)) +
    geom_col(position = position_dodge2(preserve = "single", padding = 0), color = "black") +
    facet_grid(type ~ comparison, scales = "free_y", labeller = labeller(type = intron_type_labels)) +
    labs(y = "Proportion of differentially spliced introns") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_fill_manual(name = "Acceptor/donor annotation", 
                      values = c("#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#7E6148", "grey"),
                      labels = intron_type_labels) +
    custom_gg_theme +
    theme(axis.text.x = ggplot2::element_text(color = "black", size = 8, angle = 0, hjust = 0.5)) +
    {if(!is.null(tissue)) ggtitle(paste0("Tissue: ", tissue))}
}

#' Plot the proportion significant clusters containing only annotated introns
#'
#' Obtain the introns from significant clusters and comparisons and removes the
#' junction category "ambig_gene". The proportion of annotated introns in each
#' cluster and comparison is calculated and represented.
#'
#' @param annotated list of annotated junctions. Each element corresponds to a
#'   specific study.
#' @param leafcutter_list list of leafcutter results. Each element corresponds
#'   to a specific study.
#' @param level index or name of the specific study in both the "annotated" and "leafcutter_list" arguments.
#' @param use_strand boolean to specify whether to use the strand in the
#'   clusters' ID.
#' @param tissue character vector with the specific tissue to study. The tissue
#'   name must be contained in the comparison name.
#'
#' @return
#' @export
plotDistributionOnlyAnnotatedTypes <- function(annotated, leafcutter_list, level, use_strand = F, tissue = NULL){
  # Get introns with significant comparison and cluster
  significant <- extractSignificantIntrons(annotated[[level]], leafcutter_list[[level]], use_strand, tissue)
  
  # Find number of overlapped introns between the comparisons
  number_ovelaps <- significant %>% 
    dplyr::filter(!type == "ambig_gene") %>% 
    dplyr::mutate(unique_id = str_c(cluster_id, ":", start, ":", end),
                  type = type %>% factor(levels = intron_type_levels)) %>% 
    dplyr::group_by(unique_id, type) %>% 
    dplyr::mutate(n_overlaps = n() %>% as_factor()) %>%
    dplyr::ungroup()
  
  # Count the proportion of clusters with all their introns categorized as
  # annotated
  cluster_annot_count <- number_ovelaps %>% 
    dplyr::group_by(comparison, cluster_id) %>% 
    dplyr::summarise(introns_all_annotated = all(in_ref)) %>% 
    dplyr::group_by(comparison, introns_all_annotated) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::mutate(prop = n/sum(n)) %>% 
    dplyr::ungroup()
  
  # Plot the results
  cluster_annot_count %>% 
    dplyr::mutate(comparison = as_factor(comparison)) %>% 
    ggplot(aes(x = comparison, y = prop, fill = introns_all_annotated)) +
    geom_col(color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_fill_manual(name = "All introns annotated", 
                      values = c("#E64B35", "#3C5488")) +
    labs(y = "Proportion of clusters") +
    custom_gg_theme +
    {if(!is.null(tissue)) ggtitle(paste0("Tissue: ", tissue))}
}