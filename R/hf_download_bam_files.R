loadCredentials <- function(credentials_path){
  credentials <- jsonlite::fromJSON(txt = credentials_path)
  credentials <- credentials$roleCredentials
  
  current_time <- as.numeric(Sys.time())*1000
  expiration_time <- credentials$expiration
  
  Sys.setenv("AWS_ACCESS_KEY_ID" = credentials$accessKeyId,
             "AWS_SECRET_ACCESS_KEY" = credentials$secretAccessKey,
             "AWS_DEFAULT_REGION" = "eu-west-2",
             "AWS_SESSION_TOKEN" = credentials$sessionToken)
  
  if(expiration_time > current_time){
    logger::log_info(paste0("Credentials successfully imported until ", 
                            as_datetime(expiration_time/1000),
                            " (", floor((expiration_time - current_time)/1000/3600), " hours)."))
    return(expiration_time)
  }else{
    logger::log_error(paste0("Current crendentials are expired! Stoping the download process..."))
    logger::log_error("Please, launch the 'login_aws.sh' script and follow instructions to update the credentials.")
    return(expiration_time)
  }
}

loadAtaxiaGenes <- function(){
  gtf <- rtracklayer::readGFF("/home/zchen/gtf/gencode.v38.annotation.gtf.gz")
  
  # Mendelian genes for ataxia samples genomic location:  
  ataxia_genes <- gtf %>% 
    dplyr::filter(gene_name %in% c("FXN", "ATXN1", "ATXN2", "ATXN7", "CACNA1A")) %>% 
    dplyr::distinct(gene_name, .keep_all = TRUE) %>% 
    dplyr::mutate(gene_name = str_remove(gene_name, "\\.[^.]*$")) %>%
    dplyr::mutate(gene_id = str_remove(gene_id, "\\.[^.]*$")) 
  
  return(ataxia_genes)
}

generateSampleFileName <- function(sample_name){
  sample_file_name <- paste0("merged_", sample_name, "_mapped_post_merge.BAM_Aligned.sortedByCoord.out.bam")
  
  return(sample_file_name)
}

getSampleDownloadLink <- function(sample_name,
                                  sample_type){
  if(sample_type == "Case"){
    download_link <- paste0("s3://ataxia-bulk-rnaseq/nextflow_first_attemp/Star_2_pass_by_indv/output/STAR/align/BAM_files/", generateSampleFileName(sample_name))
  }else if(sample_type == "Control"){
    download_link <- paste0("s3://ataxia-bulk-rnaseq/nextflow_first_attemp/Star_2_pass_by_indv/output/STAR/align/BAM_files/controls/", generateSampleFileName(sample_name))
  }else{
    download_link <- ""
  }
  
  return(download_link)
}

generateSamplePath <- function(sample_name,
                               sample_type,
                               main_samples_path){
  sample_path <- paste0(main_samples_path, tolower(sample_type), "/", generateSampleFileName(sample_name))
  
  return(sample_path)
}

checkIntegrity <- function(download_link,
                           file_path,
                           s3md5_path = "/home/grocamora/tools/aws-s3-integrity-check/s3md5/"){
  aws_etag <- attr(aws.s3::head_object(download_link, parse_response = F), "etag") %>%
    stringr::str_replace_all(fixed("\""), "")
  local_etag <- system(command = paste0(s3md5_path, "s3md5 16 ", file_path),
                       intern = TRUE)
  
  if(aws_etag == local_etag){
    logger::log_info("\t\t Integrity check pass!")
  }else{
    logger::log_error("\t\t Integrity check fail! Object ETag is different!")
    logger::log_error("\t\t\t ETag AWS    = ", aws_etag)
    logger::log_error("\t\t\t ETag local = ", local_etag)
  }
  
  return(aws_etag == local_etag)
}

estimateTimeLeft <- function(expiration_time,
                             minutes){
  current_time <- as.numeric(Sys.time())*1000
  
  if((current_time + minutes*60*1000) > expiration_time){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
