my_workspace <- ls()

# Download it from here: https://github.com/warnerm/devnetwork/blob/master/results/DEtests.RData
load("~/Downloads/DEtests.RData")

loaded_objects <- ls()[!(ls() %in% my_workspace)]

grab_table <- function(x, comparison, tissue){
  
  genes <- rownames(x)
  
  data.frame(
    gene = genes,
    comparison = comparison,
    tissue = tissue,
    x, stringsAsFactors = FALSE
  ) %>% as_tibble()
  
}

bind_rows(
  grab_table(bee_sexDE[[1]], "Male vs Mated queen", "Adult Head"),
  grab_table(bee_sexDE[[2]], "Male vs Mated queen", "Adult Mesosoma"),
  grab_table(bee_sexDE[[3]], "Male vs Mated queen", "Adult Gaster"),
  grab_table(bee_VM[[1]], "Mated vs Virgin queen", "Adult Head"),
  grab_table(bee_VM[[2]], "Mated vs Virgin queen", "Adult Mesosoma"),
  grab_table(bee_VM[[3]], "Mated vs Virgin queen", "Adult Gaster"),
  grab_table(beeSocial[[1]], "Nurse vs Forager", "Adult Head"),
  grab_table(beeSocial[[2]], "Nurse vs Forager", "Adult Mesosoma"),
  grab_table(beeSocial[[3]], "Nurse vs Forager", "Adult Gaster"),
  grab_table(beeTests[[1]], "Queen vs Worker", "L2 larvae"),
  grab_table(beeTests[[2]], "Queen vs Worker", "L3 larvae"),
  grab_table(beeTests[[3]], "Queen vs Worker", "L4 larvae"),
  grab_table(beeTests[[4]], "Queen vs Worker", "L5 larvae"),
  grab_table(beeTests[[5]], "Queen vs Worker", "Pupae"),
  grab_table(beeTests[[6]], "Queen vs Worker", "Adult Head"),
  grab_table(beeTests[[7]], "Queen vs Worker", "Adult Mesosoma"),
  grab_table(beeTests[[8]], "Queen vs Worker", "Adult Gaster")) %>%
  mutate(comparison_tissue = paste(comparison, tissue, sep = ": ")) %>%
  dplyr::select(gene, comparison_tissue, logFC) %>%
  mutate(logFC = logFC * -1) %>% # after this, positive logFC means queen-, female-, or mated-Q biased expression
  spread(comparison_tissue, logFC) %>%
  write_csv(path = "data/apis_gene_comparisons/Warner2018_apis_caste_logFC.csv")

rm(list = loaded_objects)
rm(loaded_objects)
