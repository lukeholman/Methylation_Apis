library(tidyverse)
library(AnnotationHub)
library(GO.db)

hub <- AnnotationHub()
# query(hub, c("mellifera","sqlite")) # Get the number for mellifera

# incomplete code attempting to get table of Droso orthologies
# xx <- hub[["AH10452"]]
# DB_connection <- dbconn(xx)
# ortho_proteins <- tbl(DB_connection, "Drosophila_melanogaster") %>%
#   collect() 
# library(UniProt.ws)
# up <- UniProt.ws(taxId=7460)
# egs = keys(up, "ORTHODB")
# ortho_proteins%>%pull() %in%
# Amel_proteins <- ortho_proteins %>% filter(species == "A.mellifera") %>% pull(inp_id)
# select(up, keys=c("Q9VH97"), columns=c("ORTHODB", "UNIPROTKB"), keytype = "ORTHODB")
# Run to get the uniprot names needed, then paste them here to get the matching Entrez IDs:
# https://www.uniprot.org/uploadlists/
# tbl(DB_connection, "Drosophila_melanogaster") %>% pull(inp_id) %>% write_lines(path = "~/Downloads/orthos.tsv")


# Get gene annotations for Apis mellifera
zz <- hub[["AH76825"]] 
DB_connection <- dbconn(zz) # DBI::dbListTables(DB_connection)

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


# For CLaire
# library(GO.db)
# claire_GO <- read.delim("~/Downloads/Missing_meanings.txt", header=F, stringsAsFactors = FALSE)
# #Get the meanings for each GO term ID from the GO.db database
# go_meanings <- suppressMessages(
#   AnnotationDbi::select(GO.db, 
#                         claire_GO$V1, c("GOID", "ONTOLOGY", "TERM")))
# names(go_meanings) <- c("GO", "ontology", "term")
# go_meanings <- distinct(go_meanings)
# write.table(go_meanings, "~/Downloads/Missing_meanings_fixed.txt", row.names = FALSE)


read_csv("data/caste_results.csv") %>% 
  left_join(tbl(db, "gene_names") %>% 
              dplyr::select(gene_symbol, gene_name, entrez_id, beebase) %>% 
              collect(n=Inf), 
            by = c("Gene symbol" = "gene_symbol")) %>% 
  filter(is.na(entrez_id)) %>%
  select(`Gene symbol`, `Gene name`) %>% distinct() %>% pull(`Gene symbol`) %>%
  write_lines("mystery_beebase.txt")
