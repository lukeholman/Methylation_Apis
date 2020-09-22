# Code to download tissue-specific expression data for fruit fly genes from "Flyatlas2"
# See this page: http://flyatlas.gla.ac.uk/FlyAtlas2/index.html?page=help

setwd("~/Rprojects/Methylation_Apis")
library(tidyverse)
library(glue)

done_already <- str_remove_all(list.files("data/flyatlas2"), "[.]rds")
fbids <- read_csv("data/database_tables/dros_ortho_GO.csv") %>% 
  filter(!is.na(FLYBASE) & !(FLYBASE %in% done_already)) %>% pull(FLYBASE) %>% unique()


get_one_gene <- function(FBgn_id){
  print(FBgn_id)
  try(foc <- read_lines(glue("http://flyatlas.gla.ac.uk/FA2Direct/index.html?fbgn={FBgn_id}&tableOut=gene")))
  Sys.sleep(1)
  if(!exists("foc")){
    print("nope, link fails"); return(NULL)
  }
  if(foc == "An error has occurred."){
    print("nope, error"); return(NULL)
  } 
  foc <- foc[which(str_detect(foc, "Tissue") & str_detect(foc, "FPKM"))[1]:length(foc)] %>% read_tsv()
  
  fpkm <- bind_rows(
    foc %>% dplyr::select(Tissue, FPKM, SD, Enrichment) %>% mutate(type = "Male"),
    foc %>% dplyr::select(Tissue, FPKM_1, SD_1, Enrichment_1) %>% mutate(type = "Female") %>% rename_all(~ str_remove_all(.x, "_1")),
    foc %>% dplyr::select(Tissue, FPKM_2, SD_2, Enrichment_2) %>% mutate(type = "Larva") %>% rename_all(~ str_remove_all(.x, "_2"))
  )
  sex_diff <- foc %>% dplyr::select(Tissue, `M/F`, `p value`)
  list(fpkm, sex_diff) %>% saveRDS(glue("data/flyatlas2/{FBgn_id}.rds"))
}

# Download one file per gene, parse a bit, and save as rds
lapply(fbids, get_one_gene)

# Now zip all the files into two .csv files:
done_files <- list.files("data/flyatlas2", full.names = TRUE, pattern = "rds")
fbids <- str_extract(done_files, "FBgn[:digit:]+")

parse_one <- function(i, type){
  if(type == "fpkm") return(read_rds(done_files[i])[[1]] %>% mutate(gene = fbids[i]))
  if(type == "sex_diff") return(read_rds(done_files[i])[[2]] %>% mutate(gene = fbids[i]))
}

fpkm <- map_df(1:length(done_files), parse_one, "fpkm")
sex_diff <- map_df(1:length(done_files), parse_one, "sex_diff")
write_csv(fpkm, "data/database_tables/flyatlas2_fpkm.csv")
write_csv(sex_diff, "data/database_tables/flyatlas2_sex_diff.csv")

# Use to delete all the individual files to save space
# unlink(list.files("data/flyatlas2", full.names = TRUE, pattern = "rds"))

