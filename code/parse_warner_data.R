library(tidyverse)

my_workspace <- ls()

# Download it from here: https://github.com/warnerm/devnetwork/blob/master/results/DEtests.RData
load("~/Downloads/DEtests.RData")

loaded_objects <- ls()[!(ls() %in% my_workspace)]

grab_table <- function(x, comparison){
  
  genes <- rownames(x)
  
  data.frame(
    gene = genes,
    comparison = comparison,
    x, stringsAsFactors = FALSE
  ) %>% as_tibble()
  
}

bind_rows(
  grab_table(bee_sexDE[[1]], "Expression in the head: Queen vs Male"),
  grab_table(bee_sexDE[[2]], "Expression in the mesosoma: Queen vs Male"),
  grab_table(bee_sexDE[[3]], "Expression in the gaster: Queen vs Male"),
  # grab_table(bee_VM[[1]], "Mated vs Virgin queen", "Adult Head"),
  # grab_table(bee_VM[[2]], "Mated vs Virgin queen", "Adult Mesosoma"),
  # grab_table(bee_VM[[3]], "Mated vs Virgin queen", "Adult Gaster"),
  grab_table(beeSocial[[1]], "Expression in the head: Nurse vs Forager"),
  grab_table(beeSocial[[2]], "Expression in the mesosoma: Nurse vs Forager"),
  grab_table(beeSocial[[3]], "Expression in the gaster: Nurse vs Forager"),
  grab_table(beeTests[[1]], "Expression in L2 larvae: Queen vs Worker"),
  grab_table(beeTests[[2]], "Expression in L3 larvae: Queen vs Worker"),
  grab_table(beeTests[[3]], "Expression in L4 larvae: Queen vs Worker"),
  grab_table(beeTests[[4]], "Expression in L5 larvae: Queen vs Worker"),
  grab_table(beeTests[[5]], "Expression in Pupae: Queen vs Worker"),
  grab_table(beeTests[[6]], "Expression in the head: Queen vs Worker"),
  grab_table(beeTests[[7]], "Expression in the mesosoma: Queen vs Worker"),
  grab_table(beeTests[[8]], "Expression in the gaster: Queen vs Worker")) %>%
  dplyr::select(gene, comparison, logFC) %>%
  mutate(logFC = logFC * -1) %>% # after this, positive logFC means queen-, female-, or mated-Q biased expression
  mutate(logFC = replace(logFC, str_detect(comparison, "Nurse vs Forager"), # after this, positive logFC means nurses > foragers
                         logFC[str_detect(comparison, "Nurse vs Forager")] * -1)) %>%
  spread(comparison, logFC) %>%
  write_csv(path = "data/apis_gene_comparisons/Warner2018_apis_caste_logFC.csv")

rm(list = loaded_objects)
rm(loaded_objects)
