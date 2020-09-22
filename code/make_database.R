library(dbplyr)
library(tidyverse)
library(RSQLite)
library(GenomicFeatures)

######################################################################
# List the files in "data/database_tables", not including tables that need indexes assigned to them
files <- list.files("data/database_tables", full.names = T)
files <- files[!(files %in% c("data/database_tables/methyl_diff_results.tsv"))]

location_of_database <- "data"
db_filepath <- file.path(location_of_database, "apis_db.sqlite3")
unlink(db_filepath)


######################################################################
# Create the empty database on disk, and add the tables in data/database_tables
apis_db <- src_sqlite(db_filepath, create = TRUE)
apis_db <- dbConnect(RSQLite::SQLite(), db_filepath)

add_table <- function(file){
  if(grepl("csv", file)) temp <- read_csv(file)
  else temp <- read_tsv(file)
  
  if("entrez_id" %in% names(temp)) temp$entrez_id <- as.character(temp$entrez_id)
  
  table_name <- gsub("[.]csv", "", tail(strsplit(file, split = "/")[[1]], 1))
  table_name <- gsub("[.]txt", "", table_name)
  table_name <- gsub("[.]tsv", "", table_name)
  
  # Add the database tables present in the data/database_tables directory
  copy_to(apis_db, 
          temp,
          table_name, 
          temporary = FALSE)
}

invisible(lapply(files, add_table))

######################################################################
# add the massive list of C and T counts for each sample/site combination, with some indexes
copy_to(apis_db,
        vroom::vroom("data/database_tables/methyl_diff_results.tsv"),
        "methylation_counts",
        indexes = list("site", "sample"),
        temporary = FALSE)

######################################################################
# Make a small table to decode the sample names, and add that too

sample_ids <- vroom::vroom("data/database_tables/methyl_diff_results.tsv", col_select = 2) %>%
  distinct() %>% 
  mutate(caste = substr(sample, 1, 1),
         time = as.numeric(str_extract(sample, "[0-9]")),
         rep = str_extract(sample, "[ABCD]")) %>%
  mutate(caste = replace(caste, caste == "q", "Queen"),
         caste = replace(caste, caste == "w", "Worker"),
         caste = replace(caste, caste == "w", "t0"))

copy_to(apis_db, 
        sample_ids,
        "sample_ids", 
        temporary = FALSE)

######################################################################
# Make table of annotations for each of the sites and add to the DB
# There is >1 row per site, since some sites have multiple annotated features

site_annotations <- read_tsv("data/database_tables/site_info.tsv") %>% 
  dplyr::select(site, seqnames, start) %>%
  add_gene_exon_columns() %>%
  dplyr::select(-seqnames, -start)

copy_to(apis_db, 
        site_annotations,
        "site_annotations", 
        temporary = FALSE)

######################################################################
# Add the significant differentially methylated 500bp windows from Claire's BWASP pipeline:

bwasp_caste_500bp <- list.files("data/BWASP_results/Caste specific results", full.names = TRUE) %>%
  map_df(~ read_tsv(.x)) %>% dplyr::select(-width, -strand)

bwasp_timeQ_500bp <- list.files("data/BWASP_results/Temporal_Variation_Q_results", full.names = TRUE) %>%
  map_df(~ read_tsv(.x)) %>% dplyr::select(-width, -strand)

bwasp_timeW_500bp <- list.files("data/BWASP_results/Temporal_Variation_W_results", full.names = TRUE) %>%
  map_df(~ read_tsv(.x)) %>% dplyr::select(-width, -strand)

copy_to(apis_db, 
        bwasp_caste_500bp,
        "bwasp_caste_500bp", 
        temporary = FALSE)

copy_to(apis_db, 
        bwasp_timeQ_500bp,
        "bwasp_timeQ_500bp", 
        temporary = FALSE)

copy_to(apis_db, 
        bwasp_timeW_500bp,
        "bwasp_timeW_500bp", 
        temporary = FALSE)



######################################################################
# Concatenate the numC and numT counts by gene

# Connect to the database we have been creating in this file (note, uses a different name)
db <- dbConnect(SQLite(), "data/apis_db.sqlite3") # DBI::dbListTables(db)

# Get the by-site counts, with gene info added
by_site <- tbl(db, "methylation_counts") %>%
  mutate(prop = numC / (numC + numT)) %>% 
  left_join(tbl(db, "sample_ids"), by = "sample") %>% 
  left_join(tbl(db, "site_annotations"), by = "site")

# Count up the Cs and Ts by gene. This includes all the sites that are between the
# stop and start codons of all the exons of the gene (i.e. CDS and UTRs, but not promoters)

fully_covered_sites <- tbl(db, "site_info") %>% # about 5 million sites
  filter(n_samples==36) %>%
  pull(site)

by_gene <- by_site %>%
  filter(!is.na(exon_in) & 
           site %in% fully_covered_sites) %>% 
  group_by(exon_in, sample, caste, time, rep) %>%
  summarise(numC = sum(numC),
            numT = sum(numT),
            prop = NA,
            num_sites = n()) %>% 
  mutate(prop = numC / (numC + numT)) %>%
  collect(n = Inf) %>% ungroup() %>%
  dplyr::rename(gene = exon_in) %>%
  dplyr::select(gene, sample, numC, numT, num_sites)



copy_to(apis_db, 
        by_gene,
        "by_gene", 
        temporary = FALSE)

######################################################################
# Do the same again, but this time only use sites that are methylated at 10x greater than the background rate
# that is dictated by bisulphite conversion failure. Use this for analysis of methylation levels of genes (e.g. for network analysis)

methylated_sites <- tbl(db, "site_info") %>%
  filter(mean_prop_meth > 0.027) %>%
  distinct() %>% pull(site)

by_site <- by_site %>%
  filter(site %in% methylated_sites)


by_gene_trimmed <- by_site %>%
  filter(!is.na(exon_in) & 
           site %in% fully_covered_sites) %>% 
  group_by(exon_in, sample, caste, time, rep) %>%
  summarise(numC = sum(numC),
            numT = sum(numT),
            prop = NA,
            num_sites = n()) %>% 
  mutate(prop = numC / (numC + numT)) %>%
  collect(n = Inf) %>% ungroup() %>%
  dplyr::rename(gene = exon_in) %>%
  dplyr::select(gene, sample, numC, numT, num_sites)



copy_to(apis_db, 
        by_gene_trimmed,
        "by_gene_trimmed", 
        temporary = FALSE)



