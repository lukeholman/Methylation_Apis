library(tidyverse)
library(AnnotationHub)
library(GO.db)

hub <- AnnotationHub()
# query(hub, c("mellifera","sqlite")) # Get the number for mellifera
# Useful database of Apis orthologs xx <- hub[["AH10452"]]
zz <- hub[["AH62534"]] # Apis mellifera
DB_connection <- dbconn(zz)
# DBI::dbListTables(DB_connection)

# mapping gene names to gene symbols to Entrez IDs. Also chromosome locations
gene_names <- tbl(DB_connection, "gene_info") %>% 
  left_join(tbl(DB_connection, "chromosomes")) %>%
  left_join(tbl(DB_connection, "genes")) %>% collect(n=Inf) %>%
  dplyr::rename(gene_name = GENENAME,
                chromosome = CHR,
                gene_symbol = SYMBOL, entrez_id = GID)

# import and parse table of new and old gene names (Entrez, old Beebase, and new Beebase)
entrez.tbl <- read.delim("data/apis_gene_comparisons/am.gene_info.txt", stringsAsFactors = FALSE)[,c(2,5,6)] 
names(entrez.tbl) <- c("entrez.id", "beebase1", "beebase2")
entrez.tbl <- entrez.tbl %>% 
  mutate(entrez.id = as.character(entrez.id), 
         beebase2 = gsub("BEEBASE:", "", beebase2))
entrez.tbl$beebase2[entrez.tbl$beebase2 == "-"] <- entrez.tbl$beebase1[entrez.tbl$beebase2 == "-"]
entrez.tbl$beebase2[grep("\\|", entrez.tbl$beebase2)] <- unname(unlist(sapply(entrez.tbl$beebase2[grep("\\|", entrez.tbl$beebase2)], function(x){
  namess <- strsplit(x, split = "\\|")[[1]]
  hits <- str_detect(namess, "GB")
  if(sum(hits) == 0) return(NA)
  return(namess[hits])
})))
entrez.tbl$beebase1[entrez.tbl$beebase1 == "-"] <- NA
entrez.tbl$beebase2[entrez.tbl$beebase2 == "-"] <- NA

gene_names <- gene_names %>% left_join(entrez.tbl, by = c("entrez_id" = "entrez.id"))

bee_GO <- tbl(DB_connection, "go") %>% left_join(tbl(DB_connection, "gene_info")) %>% 
  dplyr::select(SYMBOL, GO, ONTOLOGY) %>% collect()

#Get the meanings for each GO term ID from the GO.db database
go_meanings <- suppressMessages(
  AnnotationDbi::select(GO.db, 
                        bee_GO$GO, c("GOID", "ONTOLOGY", "TERM")))
names(go_meanings) <- c("GO", "ontology", "term")
go_meanings <- distinct(go_meanings)

# Get the KEGG annotations for bee genes (mapped to Entrez gene names)
bee_kegg_download <- clusterProfiler::download_KEGG("ame", keggType = "KEGG", keyType = "kegg")
kegg_meanings <- bee_kegg_download[[2]] %>% dplyr::rename(kegg = from, name = to)
bee_kegg <- bee_kegg_download[[1]] %>% dplyr::rename(kegg = from) %>%
  left_join(gene_names %>% dplyr::select(gene_symbol, entrez_id), by = c("to" = "entrez_id")) %>% 
  dplyr::select(kegg, gene_symbol)
rm(list = c("bee_kegg_download", "hub", "zz", "DB_connection"))

gene_names <- gene_names %>% dplyr::select(-`_id`)
invisible(lapply(c("bee_GO", "bee_kegg", "gene_names", "go_meanings", "kegg_meanings"),
       function(x){
         write_csv(get(x), path = paste("data/database_tables/", x, ".csv", sep = ""))
       }))
