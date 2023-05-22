#' Read JUNC files and extract the junctions
#'
#' From the .junc files generated using
#' \href{https://regtools.readthedocs.io/en/latest/commands/junctions-extract/}{regtools
#' junctions extract}, this functions process the columns into the desired
#' format and joins the information from junctions across all samples.
#'
#' One file is generated at the end, which contains information about all the
#' junctions found inside the samples' .junc files. It also contains the reads
#' per sample of each junction.
#'
#' @param metadata Dataframe containing all the metadata.
#' @param main_samples_path Path to where the JUNC files are stored.
#' @param num_cores Number of multiprocessing cores to use. Memory requirements
#'   significantly increase with the number of cores.
#' @param rw_disk Whether to store the results in disk. By default, TRUE.
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#'
#' @return Data.frame containing every junction found in all samples. The
#'   dataframe will contain information about the junction and the reads of that
#'   junction in each sample.
#' @export
junctionReading <- function(metadata,
                            main_samples_path,
                            num_cores = 4,
                            rw_disk = T,
                            overwrite = F,
                            output_path = ""){
  logger::log_info("\t Starting the junction reading process.")
  
  ## If rw_disk argument is set to TRUE and overwrite is set to FALSE, it looks
  ## for the results of the function in disk. If they have already been
  ## generated, the function is not executed and the files are returned.
  if (rw_disk & !overwrite & file.exists(paste0(output_path, "all_reads_combined.rds"))) {
    logger::log_info("\t\t Ignoring junction reading.") 
    logger::log_info("\t\t File ", output_path, "all_reads_combined.rds already exists.")
    logger::log_info("\t\t Loading file from disk.")
    return(readRDS(paste0(output_path, "all_reads_combined.rds")))
  }
  
  ############ Read the junctions ############
  ##
  ## 1. Loops through every .junc file located in control and cases.
  ## 2. Reads the columns and apply transformations
  ## 3. Join all the dataframes together
  
  ## Multiprocessing generation. Add argument "output.file" to
  ## parallel::makeCluster() if you want to output information about the
  ## parallel execution.
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  ## Multiprocessing loop
  logger::log_info("\t\t Reading all extracted .bam files (num_cores = ", num_cores, ").")
  all_junc <- foreach(i = 1:nrow(metadata), .packages = c("tidyverse")) %dopar% {
    ## Definition of the variables
    sample_id <- metadata[i, ] %>% pull(Correct_sample)
    cluster <- metadata[i, ] %>% pull(Type) %>% tolower()
    junc_path <- paste0(main_samples_path, cluster, "/", sample_id, ".bam.sort.s0.junc")
    
    if(!file.exists(junc_path)) return(tibble())
    
    ## Read the junction file into a tibble using readr::read_table()
    ## The "locale" argument is to read comma separated values (i.e. 25,12)
    junc <- readr::read_table(
      junc_path,
      col_names = F,
      show_col_types = F,
      progress = F,
      locale = readr::locale(grouping_mark = "")
    )
    
    ## Transformations of the junctions
    junc <- junc %>%
      dplyr::select(
        chr = X1,
        start = X2,
        stop = X3,
        junID = X4,
        reads = X5,
        strand = X6,
        blockSizes = X11
      ) %>%
      dplyr::mutate(strand = ifelse(strand == "?", "*", strand)) %>%
      dplyr::mutate(sampleID = sample_id) %>%
      tidyr::separate(col = blockSizes, sep = ",", c("blockSizesStart", "blockSizesEnd"), conver = T) %>%
      dplyr::mutate(start = start + blockSizesStart + 1, stop = stop - blockSizesEnd) %>%
      GenomicRanges::GRanges() %>%
      diffloop::rmchr() %>%
      GenomeInfoDb::keepSeqlevels(value = c(seq(1, 22) %>% as.character(), "X", "Y"), pruning.mode = "tidy") %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        seqnames = as.character(seqnames),
        strand = as.character(strand)
      ) %>%
      dplyr::select(-junID, -blockSizesStart, -blockSizesEnd)
    
    return(junc)
  }
  ## Stop the parallel cluster
  parallel::stopCluster(cl)
  
  ## Combine the resulted list into a single dataframe
  all_junc <- all_junc %>%
    dplyr::bind_rows() %>% 
    dplyr::distinct()
  invisible(gc())
  
  ############ Group and name the junctions ############
  ##
  ## 1. Name every junction by its seqname, start and end.
  ## 2. Generates a dataframe with the number of reads per junction and sample.
  ## 3. Store dataframe in disk.
  logger::log_info("\t\t Grouping the junctions across samples.")
  
  ## Generate the name of the junction by its seqname, start and end
  all_junc <- all_junc %>%
    dplyr::group_by(seqnames, start, end) %>%
    dplyr::mutate(junID = dplyr::cur_group_id(), .before = seqnames) %>%
    dplyr::ungroup()
  
  ## Generate the reads table
  all_reads_combined <- all_junc %>%
    tidyr::pivot_wider(values_from = reads, names_from = sampleID)
  
  ## Save to disk only if rw_disk is set to TRUE
  if (rw_disk) {
    logger::log_info("\t\t Saving output to ", output_path, "all_reads_combined.rds")
    all_reads_combined %>% saveRDS(paste0(output_path, "all_reads_combined.rds"))
  }
  
  logger::log_info("\t\t Reading process completed successfully.")
  return(all_reads_combined)
}