# The input, df, should be a dataframe/tibble with a column called "statistic", 
# and another column called "gene", e.g. LOC413091 or Csd
# need to have the database db loaded too
# GO_list tells the function to use either the standard A. mellifera GO annotations (bee_GO),
# or the GO terms supplemented with annotations of our genes' Drosophila orthologs (dros_ortho_GO)
# This site explains the meaning of the enrichment score metric (ES and NES): https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Enrichment_Score_(ES)

GO_and_KEGG_gsea <- function(df, GO_list, min.size = 5, keep.all = FALSE){
  
  p <- 0.01; if(keep.all) p <- 1 # Set the significance threshold
  
  # Set up the geneList object in the form needed by the fgsea function - 
  # a named, ranked vector of the test statistic (e.g. logFC, diff in % methylation)
  df <- df %>% arrange(statistic)
  geneList <- df$statistic
  names(geneList) <- df$gene
  gene_universe <- names(geneList)
  
  # Internal function to run GO enrichment
  GO_enrichment <- function(geneList, ontol, GO_list){
    
    pathways <- tbl(db, GO_list) %>% 
      left_join(tbl(db, "go_meanings"), by = c("GO")) %>% 
      distinct() %>%
      collect(n = Inf) %>%
      filter(SYMBOL %in% gene_universe, 
             ONTOLOGY == ontol) %>%       
      dplyr::select(-ONTOLOGY) 
    pathways <- with(pathways, split(SYMBOL, GO))
    
    # Run the GO enrichment K-S GSEA test
    result <- fgsea::fgsea(pathways, geneList,
                           minSize = min.size, maxSize = 500)
    
    # Collapse redundant pathways
    collapse_pathways <- fgsea::collapsePathways(result, pathways, geneList, nperm = 1000)
    pathways_to_keep <- c(collapse_pathways[[1]], names(collapse_pathways[[2]]))
    
    result <- result %>% 
      filter(pathway %in% pathways_to_keep) %>% 
      mutate(padj = p.adjust(pval)) %>% 
      filter(pval <= p) 
    
    result <- result %>% 
      left_join(tbl(db, "go_meanings") %>% 
                  dplyr::select(-ontology) %>% collect(n=Inf), 
                by = c("pathway" = "GO")) 
    
    if(nrow(result) == 0) return(NULL)
    if(ontol == "BP") Test_type <- "GO: Biological process"
    if(ontol == "MF") Test_type <- "GO: Molecular function"
    if(ontol == "CC") Test_type <- "GO: Cellular component"
    data.frame(Test_type = Test_type, 
               result %>% arrange(pval), 
               stringsAsFactors = FALSE)
  } # end GO_enrichment
  
  # Internal function to run KEGG enrichment      
  kegg_enrichment <- function(geneList, GO_list){
    
    org_code <- "dme"
    if(GO_list == "bee_GO") org_code <- "ame"
    
    apis_kegg <- clusterProfiler::download_KEGG(org_code) 
    apis_kegg_names <- apis_kegg[[2]]
    
    if(GO_list == "bee_GO"){
      apis_kegg_focal <- apis_kegg[[1]] %>% 
        left_join(tbl(db, "gene_names") %>%
                    dplyr::select(entrez_id, gene_symbol) %>%
                    collect(n=Inf) %>%
                    mutate(entrez_id = as.character(entrez_id)), 
                  by = c("to" = "entrez_id")) %>%
        filter(gene_symbol %in% gene_universe)
      pathways <- with(apis_kegg_focal, split(gene_symbol, from))
    }
    else if(GO_list %in% c("dros_ortho_GO", "dros_ortho_GO_SLIM")){
      apis_kegg_focal <- apis_kegg[[1]] %>% 
        mutate(to = str_remove_all(to, "Dmel_")) %>%
        left_join(tbl(db, "dros_ortho_GO") %>%
                    filter(!is.na(FLYBASECG)) %>%
                    dplyr::select(FLYBASECG, SYMBOL) %>%
                    collect(n=Inf), 
                  by = c("to" = "FLYBASECG")) %>%
        filter(!is.na(SYMBOL)) %>%
        filter(SYMBOL %in% gene_universe)      
      pathways <- with(apis_kegg_focal, split(SYMBOL, from))
    } # end kegg_enrichment
    
    # Run the KEGG enrichment K-S GSEA test
    result <- fgsea::fgsea(pathways, geneList,
                           minSize = min.size, maxSize = 500) 
    
    collapse_pathways <- fgsea::collapsePathways(result, pathways, geneList, nperm = 1000)
    pathways_to_keep <- c(collapse_pathways[[1]], names(collapse_pathways[[2]]))
    
    result <- result %>% 
      filter(pathway %in% pathways_to_keep) %>% 
      mutate(padj = p.adjust(pval)) %>% 
      filter(pval <= p) 
    
    result <- result %>% 
      left_join(apis_kegg_names %>% dplyr::rename(term = to), by = c("pathway" = "from")) %>%
      mutate(pathway = str_replace_all(pathway, "ame", "KEGG:"))
    
    if(nrow(result) == 0) return(NULL)
    data.frame(Test_type = "KEGG", 
               result %>% arrange(pval), 
               stringsAsFactors = FALSE) 
  }
  
  rbind(GO_enrichment(geneList, "BP", GO_list=GO_list),
        GO_enrichment(geneList, "MF", GO_list=GO_list),
        GO_enrichment(geneList, "CC", GO_list=GO_list),
        kegg_enrichment(geneList, GO_list=GO_list)) %>%
    mutate(gene_names = map(leadingEdge, ~ {
      foc_genes <- unique(.x)
      tbl(db, "gene_names") %>% 
        filter(gene_symbol %in% foc_genes) %>%
        pull(gene_name) %>% unique()
    }),
    term = as.character(term))
}


