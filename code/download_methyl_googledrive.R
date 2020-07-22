library(tidyverse)
library(googledrive)
library(glue)

search_hit <- drive_find(pattern = "_meth") 
comparisons <- substr(search_hit %>% pull(name), 1, 4)
file_ids <- search_hit %>% pull(id)

lapply(1:length(file_ids), function(i){
  focal <- comparisons[i]
  filename <- glue("data/methyl_diff_results/{focal}.txt")
  
  drive_download(
    as_id(file_ids[i]), path = filename
  )
})

list.files("data/methyl_diff_results", full.names = T) %>%
  map_df(read_tsv) %>% write_tsv("data/methyl_diff_results/methyl_diff_results.tsv")

if(file.exists("data/methyl_diff_results/methyl_diff_results.tsv")){
  unlink(list.files("data/methyl_diff_results", full.names = T, pattern = ".txt"))
}

# To open:
# vroom::vroom("data/methyl_diff_results/methyl_diff_results.tsv")














