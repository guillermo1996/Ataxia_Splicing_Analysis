#' Extract read depth from MultiQC output file
#'
#' @param metadata dataframe, contains the relevant information of all samples.
#'   Required to have a "sample_name" field.
#' @param multiqc_path character vector, path to the
#'   multiqc_rseqc_read_distribution file.
#'
#' @return metadata dataframe with a new column with the total reads.
#' @export
extractReadDepthMultiQC <- function(metadata, multiqc_path){
  multi <- readr::read_delim(multiqc_path, show_col_types = F)
  metadata_reads <- metadata %>% 
    dplyr::left_join(multi %>% 
                       dplyr::select(Sample, total_reads) %>%
                       dplyr::mutate(Sample = gsub("merged_", "", Sample)),
                     by = c("sample_name" = "Sample"))
  
  return(metadata_reads)
}

#' Load splice junctions from first pass STAR alignment
#'
#' Source:
#' https://github.com/RHReynolds/LBD-seq-bulk-analyses/blob/main/R/load_sj_df.R
#'
#' @param metadata dataframe, contains the relevant information of all samples.
#'   For the Ataxia research, it requires to have "sample_name" and
#'   "Correct_sample" columns.
#' @param sj_tab_path character vector, path to the directory containing the
#'   splice junctions files (SJ.out.tab).
#'
#' @return dataframe with all the junctions from the SJ.out.tab files.
#' @export
loadSJ <- function(metadata, sj_tab_path){
  sj_df <- foreach(i = 1:nrow(metadata)) %do%{
    # Use the metadata dataframe to locate the files.
    row <- metadata[i, ]
    wrong_sample_name <- row %>% pull(sample_name)
    correct_sample_name <- row %>% pull(Correct_sample)
    junc_path <- paste0(sj_tab_path, "merged_", wrong_sample_name, "_mapped_post_merge.BAM_SJ.out.tab")
    
    # If file does not exist, return an empty tibble
    if(!file.exists(junc_path)) return(tibble())
    
    # Read the junction file into a tibble
    junc <- readr::read_delim(junc_path,
                              delim = "\t",
                              col_names = c("chr", "intron_start", "intron_end", "strand", "intron_motif", "intron_annotation", 
                                            "unique_reads_junction", "multimap_reads_junction", "max_splice_alignment_overhang"),
                              col_types = "cdddddddd") %>%
      dplyr::mutate(sample = correct_sample_name)
  } %>% dplyr::bind_rows()
  
  return(sj_df)
}

# referenceWithCHR <- function(TxDb_ref){
#   seqnames_ref <- intronsByTranscript(TxDb_ref) %>% 
#     unlist() %>% 
#     GenomicRanges::seqnames() %>%
#     unique()
#   
#   chr_prop <- sum(grepl("chr", seqnames_ref))/length(seqnames_ref)
#   if(chr_prop > 0.5){
#     return(T)
#   }else{
#     return(F)
#   }
# }

#' Convert output of STAR junction aligment to a RangedSummarisedExperiment
#' class
#'
#' Source:
#' https://github.com/RHReynolds/LBD-seq-bulk-analyses/blob/main/R/create_rse_jx.R
#'
#' @param sj_df dataframe, contains all the junctions extracted using
#'   \code{loadSJ}.
#' @param gtf_path character vector, path to the reference transcriptome.
#' @param sample_info dataframe, contains the relevant information of all
#'   samples. For the Ataxia research, it requires to have the column specified
#'   in the \code{sample_id_col} argument.
#' @param sample_id_col character vector, ID field in the \code{sample_info}
#'   dataframe argument.
#'
#' @return RangedSummarisedExperiment. Output of STAR junction alignment in
#'   RangedSummarisedExperiment class.
#' @export
createRSE <- function(sj_df, gtf_path, sample_info, sample_id_col){
  # Load reference transcriptome and other function variables
  TxDb_ref <- dasper:::ref_load(gtf_path)
  strand_map = c("*", "+", "-")
  valid_seqnames <- c(str_c(1:22), "X", "Y")
  
  # Some reference transcriptomes (i.e. Gencode) use the characters "chr" before
  # the seqnames. To bring consistency, the "seqlevelsStyle" is set to "NCBI".
  seqlevelsStyle(TxDb_ref) <- c("NCBI")
  
  # Annotate the introns
  row_ranges <- sj_df %>%
    dplyr::distinct(chr, intron_start, intron_end, strand) %>%
    dplyr::mutate(strand = strand_map[strand+1],
                  junction_id = stringr::str_c(chr, ":", intron_start, "-", intron_end, ":", strand)) %>%
    GenomicRanges::GRanges() %>%
    diffloop::rmchr() %>%
    GenomeInfoDb::keepSeqlevels(value = valid_seqnames, pruning.mode = "tidy") %>%
    dasper::junction_annot(ref = TxDb_ref)
  
  # Create the junction counts table
  counts <- sj_df %>%
    dplyr::select(-intron_motif, -intron_annotation, -multimap_reads_junction, -max_splice_alignment_overhang) %>%
    dplyr::filter(chr %in% str_c("chr", valid_seqnames)) %>%
    dplyr::mutate(strand = strand_map[strand+1],
                  junction_id = stringr::str_c(chr, ":", intron_start, "-", intron_end, ":", strand)) %>%
    dplyr::select(sample, junction_id, unique_reads_junction) %>%
    tidyr::pivot_wider(values_from = unique_reads_junction, names_from = sample) %>%
    replace(is.na(.), 0)
  
  # QC tests
  junc_match <- all(row_ranges$junction_id %in% counts$junction_id)
  junc_order <- all(row_ranges$junction_id == counts$junction_id)
  if(junc_match){
    if(!junc_order){
      ids <- match(row_ranges$junction_id, counts$junction_id)
      counts <- counts[ids, ]
    }else{
      print("Order of junctions in row_ranges matches order of junctions in counts table.")
    }
  }else{
    stop("Not all junctions in row_ranges are present in the counts table.")
  }
  
  # Transform counts to matrix
  counts_matrix <- counts %>% 
    tibble::column_to_rownames("junction_id") %>%
    as.matrix()
  
  # QC tests
  sample_match <- all(colnames(counts_matrix) %in% sample_info[, sample_id_col, drop = T])
  sample_order <- all(colnames(counts_matrix) == sample_info[, sample_id_col, drop = T])
  if(sample_match == "TRUE"){
    if(sample_order == "FALSE"){
      ids <- match(colnames(counts_matrix), sample_info[, sample_id_col, drop = T])
      sample_info <- sample_info[ids,]
    }else{
      print("Order of samples in counts matrix matches order of samples in sample info.")
    }
  }else{
    stop("Not all samples in counts matrix are present in the sample info.")
  }
  
  # Generate colData
  colData <- sample_info %>% 
    dplyr::filter(!!sym(sample_id_col) %in% colnames(counts_matrix)) %>%
    tibble::column_to_rownames(sample_id_col)
  
  # Create the RSE_jx object
  rse_jx <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts_matrix),
                                                       rowRanges = row_ranges,
                                                       colData = colData)
  
  return(rse_jx)
}

#' Loads the ENCODE blacklisted regions into memory
#'
#' The different versions can be downloaded
#' \href{https://github.com/Boyle-Lab/Blacklist/tree/master/lists}{here}.
#'
#' @param blacklist_path path to the ENCODE blacklisted regions .bed file.
#'
#' @return the ENCODE blacklisted region in GRanges object.
#' @export
loadEncodeBlacklist <- function(blacklist_path) {
  if(!exists("encode_blacklist_hg38")){
    encode_blacklist_hg38 <<- rtracklayer::import(blacklist_path) %>% diffloop::rmchr()
  }
  
  return(encode_blacklist_hg38)
}

#' Remove the junctions from the ENCODE blacklisted regions
#'
#' @param rse_jx RangedSummarisedExperiment, output of \code{createRSE}.
#' @param encode_blacklist_hg38 GRanges class object with the blacklisted
#'   regions.
#'
#' @return RangedSummarisedExperiment with only the junctions that do not
#'   overlap with the blacklisted regions.
#' @export
removeEncodeBlacklistRegionsRSE <- function(rse_jx, encode_blacklist_hg38) {
  ## Look fot the overlaps between rse_jx and the ENCODE blacklisted region
  overlaps <- GenomicRanges::findOverlaps(
    query = encode_blacklist_hg38,
    subject = rowRanges(rse_jx),
    ignore.strand = F,
    type = "any"
  )
  
  idxs <- S4Vectors::subjectHits(overlaps)
  
  ## If an overlap is found, remove the junctions
  if(length(idxs) > 0) rse_jx <- rse_jx[-idxs, ]
  
  return(rse_jx)
}


#' Remove the junctions shorter than 25bp.
#'
#' @param rse_jx RangedSummarisedExperiment, output of \code{createRSE}.
#'
#' @return RangedSummarisedExperiment with only the junctions bigger than 25 bp.
#' @export
removeShortJunctionsRSE <- function(rse_jx) {
  idxs <- which(rowRanges(rse_jx) %>% width() < 25)
  if(length(idxs) > 0) rse_jx <- rse_jx[-idxs, ]
  
  return(rse_jx)
}

#' Remove the junctions from ambiguous genes.
#'
#' @param rse_jx RangedSummarisedExperiment, output of \code{createRSE}.
#'
#' @return RangedSummarisedExperiment with only the junctions associated to one
#'   or less genes.
#' @export
removeAmbiguousGenesRSE <- function(rse_jx) {
  idxs <- which(sapply(rowRanges(rse_jx)$gene_id_junction, length) >= 2)
  if(length(idxs) > 0) rse_jx <- rse_jx[-idxs, ]
  
  return(rse_jx)
}

#' Plots the category percentages by sample
#'
#' @param annot_junc_prop_df dataframe, contains the proportion of junctions for
#'   each category and each sample.
#' @param level character vector, specifies the level of study. In this
#'   function, is the name of the column that will split the graph and the
#'   Wilcoxon tests. For the Ataxia RNAseq data, only "Region", "AtaxiaSubtype"
#'   and "Diagnosis" are valid inputs.
#' @param split_tissue boolean, whether to split the graph in facets by tissue.
#'
#' @return
#' @export
plotJunctionCategories <- function(annot_junc_prop_df, level, split_tissue = F, ref_group = NULL){
  p <- ggplot(annot_junc_prop_df,
              aes(x = type, y = prop_junc, fill = !!sym(level))) +
    geom_boxplot() +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.15)), limits = c(0, NA), breaks = seq(0, 1, 0.2)) +
    labs(x = "", y = "Proportion of junctions") +
    scale_fill_manual(values = pal_jco("default", alpha = 0.9)(10)[c(3, 1, 9, 2, 7)]) + 
    ggpubr::geom_pwc(aes(group = !!sym(level)), ref.group = ref_group, vjust = -0.5, hide.ns = T,
                     label = " p = {p.adj.format}", p.adjust.method = "bonferroni", p.adjust.by = "group") +
    custom_gg_theme
  if(split_tissue) p <- p + facet_wrap(vars(Region), ncol = 1)
  
  return(p)
}

#' Plots the proportion of counts for each category and sample
#'
#' @param annot_junc_prop_df dataframe, contains the proportion of counts for
#'   each category and sample.
#' @param level character vector, specifies the level of study. In this
#'   function, is the name of the column that will split the graph and the
#'   Wilcoxon tests. For the Ataxia RNAseq data, only "Region", "AtaxiaSubtype"
#'   and "Diagnosis" are valid inputs.
#' @param tissue character vector, tissue to use to graph. Leave as NULL if no
#'   distinction is used.
#' @param ref_group character vector, category in the \code{level} field that
#'   will be compare against for the Wilcoxon signed-rank test.
#'
#' @return
#' @export
plotJunctionCountCategories <- function(annot_junc_prop_df, level, tissue = NULL, ref_group = NULL){
  if(!is.null(tissue)) annot_junc_prop_df <- annot_junc_prop_df %>% dplyr::filter(Region == tissue)
  
  p <- ggplot(annot_junc_prop_df, aes(x = !!sym(level), y = prop_count, fill = !!sym(level))) +
    geom_boxplot() +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    labs(x = "", y = "Proportion of junctions") +
    scale_fill_manual(values = pal_jco("default", alpha = 0.9)(10)[c(1, 9)]) +
    facet_wrap(vars(type), ncol = 3, scales = "free") +
    ggpubr::geom_pwc(aes(group = !!sym(level)), ref.group = ref_group, vjust = -0.3, hide.ns = T,
                     label = " p = {p.adj.format}", p.adjust.method = "bonferroni", p.adjust.by = "group") +
    custom_gg_theme_subtitle# + theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"))
  
  if(!is.null(tissue)){
    p <- p +
      scale_fill_manual(values = pal_jco("default", alpha = 0.9)(10)[c(3, 1, 9, 2, 7)]) + 
      ggtitle("", subtitle = paste0("Tissue: ", tissue))
  }
  
  return(p)
}