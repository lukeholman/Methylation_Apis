library(tidyverse)
library(AnnotationHub)
library(GO.db)

hub <- AnnotationHub()
# query(hub, c("mellifera","sqlite")) # Get the number for mellifera

# Get gene annotations for Apis mellifera - NB, the code for Apis often changes over time: 
# run "query" again to find the new one. I guess AnnotationHub don't like keeping things consistent...
zz <- hub[["AH81619"]] 
DB_connection <- dbconn(zz) # DBI::dbListTables(DB_connection)

# mapping gene names to gene symbols to Entrez IDs. Also chromosome locations
gene_names <- tbl(DB_connection, "gene_info") %>% 
  left_join(tbl(DB_connection, "chromosomes"), by = "_id") %>%
  left_join(tbl(DB_connection, "genes"), by = "_id") %>% collect(n=Inf) %>%
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

############################################################
# Get the various gene names:
gene_names <- gene_names %>% left_join(entrez.tbl, by = c("entrez_id" = "entrez.id"))

# Remove this one gene. Seems to be an error - there are 2 genes called Fabp listed in gene_names with symbol = "Fabp",
# one called "fatty acid binding protein", one called "FABP-like protein", that seem to be the same gene
gene_names <- filter(gene_names, gene_name != "FABP-like protein")


############################################################
# Get the GO annotations for A. mellifera
bee_GO <- tbl(DB_connection, "go") %>% 
  left_join(tbl(DB_connection, "gene_info"), by = "_id") %>% 
  dplyr::select(SYMBOL, GO, ONTOLOGY) %>% collect()

############################################################
# Get the KEGG annotations for bee genes (mapped to Entrez gene names)
# bee_kegg_download <- clusterProfiler::download_KEGG("ame", keggType = "KEGG", keyType = "kegg")
# kegg_meanings <- bee_kegg_download[[2]] %>% dplyr::rename(kegg = from, name = to)
# bee_kegg <- bee_kegg_download[[1]] %>% dplyr::rename(kegg = from) %>%
#   left_join(gene_names %>% dplyr::select(gene_symbol, entrez_id), by = c("to" = "entrez_id")) %>% 
#   dplyr::select(kegg, gene_symbol)
# rm(list = c("bee_kegg_download", "hub", "zz", "DB_connection"))


############################################################
# Get the GO annotations for Drosophila orthologs of our Amel genes, and make an expanded set of GO annotations that assumes the Drosophila orthologs have the same functions as the corresponding Amel genes. It also means that the conserved genes are better-annotated than the bee-specific genes.

# Get the Dmel - Amel orthologs created using code from Warner et al. 2019 Nature Comms
Amel_Dmel_orthologs <- read_csv("data/Amel_Dmel_orthologs.csv") %>% dplyr::select(-OGGend)


# Get the Dmel GO annotations from the org.db
# BiocManager::install("org.Dm.eg.db") 
library(org.Dm.eg.db)
dros_info <- select(org.Dm.eg.db, keys= keys(org.Dm.eg.db), columns = c("SYMBOL","FLYBASECG", "GO", "ONTOLOGY", "ENTREZID"))

# Merge the Dmel and Amel GO annotations, and save
# This dataframe lists all unique GO terms for the Apis gene *and* its ortholog(s) in Drosophila
dros_ortho_GO <- dros_info %>% 
  left_join(Amel_Dmel_orthologs,
            by = c("SYMBOL" = "gene_Dmel")) %>%
  as_tibble() %>%
  filter(!is.na(gene_Amel)) %>%
  dplyr::select(gene_Amel, GO, ONTOLOGY, FLYBASECG, ENTREZID) %>%
  dplyr::rename(SYMBOL = gene_Amel) %>%
  full_join(bee_GO, by = c("SYMBOL", "GO", "ONTOLOGY")) %>%
  distinct() 

# Add annotations to the three major sex determination genes in Apis, namely Csd, Fem, and Dsx
# See https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000222#abstract2 for evidence of their function in sex differentiation
# tra2 is transformer2, a major player in Drosophila sex determination
# LOC113219102 is named "sex determination protein fruitless-like"; fruitless is important in Drosophila sex det.
custom_ones <- tibble(SYMBOL = c("Csd", "Fem", "Dsx", "tra2", "LOC113219102"),
                      GO = "GO:0007530", # this is the GO number for "sex determination"
                      ONTOLOGY = "BP")

dros_ortho_GO <- dros_ortho_GO %>% 
  full_join(custom_ones, by = c("SYMBOL", "GO", "ONTOLOGY")) %>%
  arrange(SYMBOL, ONTOLOGY) 

bee_GO <- bee_GO %>% 
  full_join(custom_ones, by = c("SYMBOL", "GO", "ONTOLOGY")) %>%
  arrange(SYMBOL, ONTOLOGY) 




############################################################
#Get the meanings for each GO term ID from the GO.db database
go_meanings <- suppressMessages(
  AnnotationDbi::select(GO.db, 
                        unique(c(bee_GO$GO, dros_ortho_GO$GO)), 
                        c("GOID", "ONTOLOGY", "TERM")))
names(go_meanings) <- c("GO", "ontology", "term")
go_meanings <- distinct(go_meanings)

############################################################
# Trim the GO annotations to just include the genes in Expression_data/Genes_2020_version.txt
gene_to_keep <- rownames(as.matrix(read.table("data/Expression_data/Genes_2020_version.txt", row.names = 1, header = T)))
bee_GO <- bee_GO %>% filter(SYMBOL %in% gene_to_keep)
dros_ortho_GO <- dros_ortho_GO %>% filter(SYMBOL %in% gene_to_keep)



############################################################
# Save the spreadsheets, which will become tables in the database 
gene_names <- gene_names %>% dplyr::select(-`_id`)
invisible(lapply(c("bee_GO",  "gene_names", "go_meanings", "dros_ortho_GO"), # , "kegg_meanings" "bee_kegg",
       function(x){
         write_csv(get(x), path = paste("data/database_tables/", x, ".csv", sep = ""))
       }))





