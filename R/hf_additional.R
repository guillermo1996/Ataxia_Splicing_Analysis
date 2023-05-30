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
                     by = c("sample_name" = "Sample")) %>%
    dplyr::ungroup()
  
  return(metadata_reads)
}

#' Loads the reference genome into memory (ENSEMBL)
#'
#' The different versions can be downloaded
#' \href{http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/}{here}.
#'
#' @param gtf_path path to the reference genome .gtf file.
#'
#' @return the connection to the reference genome DB.
#' @export
loadEdb <- function(gtf_path) {
  if (!exists("edb")) {
    logger::log_info("\t\t Loading the v105 reference genome.")
    tryCatch({
      edb <<- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
      edb <<- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
    }, error = function(x){
      logger::log_warn("\t\t Reference transcriptome is not from ENSEMBL. Some functionality might be limited.")
      edb <<- loadGTF(gtf_path)
    })
      
  } else {
    logger::log_info("\t\t Variable 'edb' already loaded!")
  }
  
  return(edb)
}


#' Loads the reference genome into memory (generic)
#'
#' @param gtf_path path to the reference genome .gtf file.
#'
#' @return TxDb-class object
#' @export
loadGTF <- function(gtf_path){
  ref <- GenomicFeatures::makeTxDbFromGFF(gtf_path)
  GenomeInfoDb::seqlevelsStyle(ref) <- "NCBI"
  
  return(ref)
}

isGeneSymbolAvailable <- function(edb){
  return("SYMBOL" %in% AnnotationDbi::columns(edb))
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
  if (!exists("encode_blacklist_hg38")) {
    logger::log_info("\t\t Loading the v2 ENCODE blacklisted regions.")
    encode_blacklist_hg38 <<- rtracklayer::import(blacklist_path) %>% diffloop::rmchr()
  } else {
    logger::log_info("\t\t Variable 'encode_blacklist_hg38' is already loaded!")
  }
  
  return(encode_blacklist_hg38)
}

isFastaSeqnamesStyleNCBI <- function(fasta_path,
                                     samtools_path){
  fai_path <- paste0(fasta_path, ".fai")
  if(!file.exists(fai_path)){
    ## Index the fasta file
    system2(
      command = paste0(samtools_path, "samtools"),
      args = c(
        "faidx", fasta_path,
        "-o", paste0(fasta_path, ".fai")
      )
    )
  }
  
  fasta_seqnames <- readr::read_delim(fai_path, show_col_types = F, col_names = F) %>% pull(X1)
  chr_seqnames <- str_c("chr", c(1:22, "X", "Y"))
  
  return(!(all(chr_seqnames %in% fasta_seqnames)))
}