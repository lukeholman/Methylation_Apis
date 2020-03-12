library(dbplyr)
library(tidyverse)
library(RSQLite)

files <- list.files("data/database_tables", full.names = T)

location_of_database <- "data"
db_filepath <- file.path(location_of_database, "apis_db.sqlite3")
unlink(db_filepath)

# Create the empty database on disk
apis_db <- src_sqlite(db_filepath, create = TRUE)
apis_db <- dbConnect(RSQLite::SQLite(), db_filepath)

add_table <- function(file){
  temp <- read.csv(file, stringsAsFactors = FALSE)
  if("entrez_id" %in% names(temp)) temp$entrez_id <- as.character(temp$entrez_id)
  copy_to(apis_db, 
          temp,
          gsub("[.]csv", "", tail(strsplit(file, split = "/")[[1]], 1)), 
          temporary = FALSE)
}

invisible(lapply(files, add_table))