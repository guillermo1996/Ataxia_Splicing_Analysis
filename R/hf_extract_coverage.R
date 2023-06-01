extractCoverage <- function(file_path,
                            sample_correct_name,
                            sample_diagnosis,
                            coverage_metadata_path = "/home/grocamora/RytenLab-Research/11-Ataxia-bulk/Zhongbo_Script/All_patient_coverage.csv",
                            coverage_bam_files_path = "/home/grocamora/RytenLab-Research/11-Ataxia-bulk/Zhongbo_Script/all_filtered_BAM/"){
  ## Import pre-existing metadata
  if(file.exists(coverage_metadata_path)){
    logger::log_info("\t\t Previous metadata already found.")
    all_cov_df <- readr::read_delim(coverage_metadata_path, show_col_types = F)
    
    ## If sample is already found in metadata, then we skip the process
    found_sample_genes <- all_cov_df[all_cov_df$sample_correct_name == sample_correct_name, ] %>% 
      dplyr::pull(gene_to_filter) %>% 
      sort()
    if(isTRUE(all.equal(found_sample_genes, ataxia_genes$gene_name %>% sort))){
      logger::log_info("\t\t All ataxia genes already found in the previous metadata. Skipping the process.")
      return()
    }
  }else{
    logger::log_info("\t\t No previous metadata found.")
    all_cov_df <- tibble::tibble()
  }
  
  ## Import the BAM file and calculates the coverage
  logger::log_info("\t\t Loading the BAM file.")
  bam_file <- rtracklayer::import(file_path)
  logger::log_info("\t\t Extracting the coverage.")
  coverage_total <- IRanges::coverage(bam_file) 
  coverage_all <- coverage_total %>% unlist() %>% sum()
  
  ## Loop through every ataxia gene
  logger::log_info("\t\t Looping through ataxia genes.")
  for(j in 1:nrow(ataxia_genes)){
    ## Gene and chr
    gene_to_filter <- ataxia_genes$gene_name[j]
    chr_to_filter <- as.character(ataxia_genes$seqid[j])
    logger::log_info("\t\t\t Extracting coverage for gene ", gene_to_filter, ".")
    
    ## If metadata for the current sample and gene is already found, exits the
    ## loop
    if(nrow(all_cov_df) > 0){
      if(nrow(merge(all_cov_df %>% select(sample_correct_name, gene_to_filter),
                    tibble(sample_correct_name, gene_to_filter))) > 0) next
    }
    
    ## Filters the BAM file and calculates the coverage
    filtered_BAM <- bam_file[seqnames(bam_file) == chr_to_filter & 
                               start(bam_file) >= ataxia_genes$start[j] & 
                               end(bam_file) <= ataxia_genes$end[j]]
    filtered_coverage <- IRanges::coverage(filtered_BAM) 
    filtered_coverage_all <- filtered_coverage %>% unlist() %>% sum() %>% as.numeric()
    
    ## Save in disk the filtered BAM file
    bam_file_name <- paste(sample_correct_name, sample_diagnosis, gene_to_filter, sep = "_")
    con <- file.path(coverage_bam_files_path, bam_file_name)
    logger::log_info("\t\t\t Storing filtered BAM file for gene ", gene_to_filter, ".")
    rtracklayer::export(filtered_BAM, con, format = "bam")
    
    ## Combine the sample/gene metadata with the previous samples
    all_cov <- cbind(sample_correct_name, 
                     sample_diagnosis, 
                     gene_to_filter, 
                     coverage_all, 
                     filtered_coverage_all) %>% 
      tibble::as_tibble()
    all_cov_df <- rbind(all_cov_df, all_cov)
  }
  
  ## Saves the metadata to disk
  logger::log_info("\t\t Writing metadata fo file.")
  write.table(all_cov_df, coverage_metadata_path, row.names = F, quote = FALSE, sep = ",")
  
  ## Free some memory
  rm(bam_file, coverage_total, coverage_all, all_cov_df, gene_to_filter, chr_to_filter,
     filtered_BAM, filtered_coverage, filtered_coverage_all)
  invisible(gc())
}

#' Load ataxia genes information from reference transcriptome.
#'
#' @return dataframe with ataxia genes information.
#' @export
loadAtaxiaGenes <- function(){
  gtf <- rtracklayer::readGFF("/home/grocamora/RytenLab-Research/Additional_files/GENCODE/gencode.v38.annotation.gtf.gz")
  
  # Mendelian genes for ataxia samples genomic location:  
  ataxia_genes <- gtf %>% 
    dplyr::filter(gene_name %in% c("FXN", "ATXN1", "ATXN2", "ATXN7", "CACNA1A")) %>% 
    dplyr::distinct(gene_name, .keep_all = TRUE) %>% 
    dplyr::mutate(gene_name = str_remove(gene_name, "\\.[^.]*$")) %>%
    dplyr::mutate(gene_id = str_remove(gene_id, "\\.[^.]*$")) 
  
  return(ataxia_genes)
}