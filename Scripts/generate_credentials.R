## _________________________________________________
##
## Generate AWS credentials
##
## Author: Mr. Guillermo Rocamora Pérez
##
## Date Created: 2023-01-17
##
## Copyright (c) Guillermo Rocamora Pérez, 2023
##
## Email: guillermorocamora@gmail.com
## _________________________________________________
##
## Notes:
##
## Credentials are only generated for 12 hours.
##
## Some privacy information was removed from this script, so it will not run by
## default. You require your AWS account ID and your json cache file, usually
## located at "~/.aws/sso/cache/".
##
## Please contact guillermorocamora@gmail.com for further assitance.
## _________________________________________________

# Initial setup ----

## Required libraries
shhh <- suppressPackageStartupMessages
shhh(library(jsonlite))
shhh(library(logger))
shhh(library(lubridate))

## Logger options
logger_layout <- layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)

## Script parameters

# Enter your own credentials in here. Emptied for security reasons
account_id <- ""
role_name <- "S3BucketAccess"
profile <- "grocamora"
credentials_path <- here::here("variables/credentials/role_credentials.json")
always_renew <- T

# Generating credentials ----
logger::log_info("Generating new credentials in ", credentials_path, ".")

if(file.exists(credentials_path)){
  old_credentials <- jsonlite::fromJSON(txt = credentials_path)
  old_expiration_date <- old_credentials$roleCredentials$expiration
}else{
  old_expiration_date <- 0
}

current_time <- as.numeric(Sys.time())*1000

if(current_time > old_expiration_date | always_renew){
  if(always_renew){
    logger::log_info("\t Variable to always renew the credentials is set to TRUE.")
    logger::log_info("\t Ignoring previous credentials.")
  }else{
    logger::log_info(paste0("\t Old credentials expired in ", as_datetime(old_expiration_date/1000)))
  }
  
  # You need to know your sso cache file. It should be located in "~/.aws/sso/cache/"
  sso_cache_file <- path.expand("~/.aws/sso/cache/________.json")
  json_data <- jsonlite::fromJSON(txt = sso_cache_file)
  
  role_credentials <- jsonlite::fromJSON(system(paste("aws sso get-role-credentials",
                                                      paste("--account-id", account_id),
                                                      paste("--role-name", role_name),
                                                      paste("--profile", profile),
                                                      paste("--region", json_data$region),
                                                      paste("--access-token", json_data$accessToken)),
                                                intern = TRUE))
  new_expiration_date <- role_credentials$roleCredentials$expiration
  write(jsonlite::toJSON(role_credentials, pretty = TRUE), credentials_path)
  
  logger::log_info(paste0("\t New credentials expire in ", 
                          as_datetime(new_expiration_date/1000), 
                          " (", round((new_expiration_date - current_time)/1000/3600), " hours)."))
  # expiration_UNIX <- role_credentials$roleCredentials$expiration/1000
  # expiration_date <- stringr::str_replace_all(as_datetime(expiration_UNIX), " ", "_")
}else{
  logger::log_info(paste0("\t Current credentials expire in ", 
                          as_datetime(old_expiration_date/1000), 
                          " (", round((new_expiration_date - current_time)/1000/3600), " hours)."))
}
