# get the pheromone sensitivity data from Holman et al. 2019 Nat. Comms.
# the .db object is available on my Github: https://github.com/lukeholman/queen-pheromone-RNAseq/tree/master/data
my_db <- src_sqlite("~/Rprojects/queen-pheromone-RNAseq/data/queen_pheromone.db")
tbl(my_db, "ebseq_gene_am") %>%
  dplyr::select(gene, PostFC) %>%
  collect(n=Inf) %>%
  mutate(`Pheromone sensitivity` = -log(PostFC)) %>%
  dplyr::select(-PostFC) %>% write_csv("data/apis_gene_comparisons/pheromone_sensitivity.csv")

# Get the data from Warner et al.
# read_tsv("https://raw.githubusercontent.com/warnerm/devnetwork/master/data/bees.counts.txt") # Warner's gene-by-sample expression matrix
load("~/Downloads/DEtests.RData") # file is at: https://github.com/warnerm/devnetwork/raw/master/results/DEtests.RData, the link does not work with load(), but the file can be downloaded then opened
xx <- beeRes_allstage[[1]] 
for(i in 2:ncol(xx)) xx[,i] <- xx[,i] * -1
xx %>% write_csv("data/apis_gene_comparisons/Warner2019_apis_caste_logFC.csv")
