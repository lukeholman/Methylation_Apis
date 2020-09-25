get_bwasp_sites <- function(p_cutoff, q_cutoff, site_list = NULL){
  
  df <- list.files("data/BWASP_results/sites", 
                   full.names = TRUE, recursive = TRUE) %>%
    map_df(~ read_tsv(.x) %>% # filter by p and q
             filter(pvalue < p_cutoff & qvalue < q_cutoff)) %>% 
    left_join(tbl(db, "site_info") %>%
                dplyr::select(site, seqnames, start) %>% collect(), 
              by = c("seqnames", "start")) 
  
  # optionally, filter out sites not in site_list
  if(!is.null(site_list)) df <- df %>% filter(site %in% site_list)
  
  df %>%
    full_join(tbl(db, "site_annotations") %>% collect()) %>%
    arrange(qvalue) %>% as.data.frame() %>%
    left_join(tbl(db, "gene_names") %>% 
                dplyr::select(gene_symbol, gene_name) %>% collect(), 
              by = c("exon_in" = "gene_symbol")) %>%
    dplyr::rename(exon_gene_name = gene_name) %>%
    left_join(tbl(db, "gene_names") %>% 
                dplyr::select(gene_symbol, gene_name) %>% collect(), 
              by = c("intron_in" = "gene_symbol")) %>%
    dplyr::rename(intron_gene_name = gene_name) %>%
    left_join(tbl(db, "gene_names") %>% 
                dplyr::select(gene_symbol, gene_name) %>% collect(), 
              by = c("promoter_of" = "gene_symbol")) %>%
    dplyr::rename(promoter_gene_name = gene_name) %>%
    left_join(tbl(db, "gene_names") %>% 
                dplyr::select(gene_symbol, gene_name) %>% collect(), 
              by = c("fiveUTR_in" = "gene_symbol")) %>%
    dplyr::rename(fiveUTR_gene_name = gene_name) %>%
    left_join(tbl(db, "gene_names") %>% 
                dplyr::select(gene_symbol, gene_name) %>% collect(), 
              by = c("threeUTR_in" = "gene_symbol")) %>%
    dplyr::rename(threeUTR_gene_name = gene_name) %>%
    filter(!is.na(seqnames)) %>%
    dplyr::select(site, everything()) %>%
    as_tibble()
}