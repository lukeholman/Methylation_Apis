library(dbplyr)
library(tidyverse)
library(RSQLite)

# files <- list.files("data/database_tables", full.names = T)
files <- list.files("data/database_tables_new", full.names = T)
files <- files[files != "data/database_tables_new/uniprot_entrez.txt"]

location_of_database <- "data"
db_filepath <- file.path(location_of_database, "apis_db.sqlite3")
unlink(db_filepath)

# Create the empty database on disk
apis_db <- src_sqlite(db_filepath, create = TRUE)
apis_db <- dbConnect(RSQLite::SQLite(), db_filepath)

add_table <- function(file){
  if(grepl("csv", file)) temp <- read_csv(file)
  else temp <- read_tsv(file)

  if("entrez_id" %in% names(temp)) temp$entrez_id <- as.character(temp$entrez_id)
  
  table_name <- gsub("[.]csv", "", tail(strsplit(file, split = "/")[[1]], 1))
  table_name <- gsub("[.]txt", "", table_name)
  
  copy_to(apis_db, 
          temp,
          table_name, 
          temporary = FALSE)
}

invisible(lapply(files, add_table))
