# library(tidyverse)
# library(DBI)
# library(RSQLite)
source("code/parse_warner_data.R")

gene_names <- tbl(db, "gene_names") %>% collect() 
#x <- read_tsv("data/HaV3.1_gene_names.txt")

# beebase_converter <- read_csv("data/database_tables_new/beebase1_to_beebase2.csv") 


# Get the present study's expresison limma::voom results, and associate with their beebase IDs for matching the older studies
present_study_data <- read_csv("output/caste_results.csv") %>% 
  dplyr::filter(Time == 8) %>%              # USES THE 8-HOUR CASTE DIFFERENCE AND EXPRESSION LEVEL RESULTS
  left_join(tbl(db, "gene_names") %>% 
              dplyr::select(gene_symbol, gene_name, entrez_id, beebase2) %>% 
              collect(n = Inf), 
            by = c("Gene symbol" = "gene_symbol")) %>% 
  dplyr::select(beebase2, logFC, AveExpr, `Gene symbol`, entrez_id, gene_name) %>%
  dplyr::rename(`Caste difference in expression in 8h larvae (positive: Q > W)` = logFC) %>%
  mutate(entrez_id = replace(entrez_id, is.na(entrez_id), gsub("LOC", "", `Gene symbol`[is.na(entrez_id)]))) %>%
  dplyr::rename(beebase = beebase2)

# Merge with the present study's caste difference in methylation data 
# (estimated by Fisher exact test in BWASP, then averaged over sites within each gene)
# Note that a positive meth.diff value in comparison A.vs.B means that the B methylation percentage is higher than the A methylation percentage.
# So, a meth.diff value of +10 means that workers have a 10% higher value of %mC than queens

bwasp_gene_meth <- read_tsv("data/meth_network_input/Network_persite.txt") %>%
  gather("sample", "methylation", -Genes) %>%
  mutate(caste = case_when(
    str_detect(sample, "queen") ~ "queen",
    str_detect(sample, "worker") ~ "worker",
    TRUE ~ "t0"
  )) %>% filter(caste!="t0") %>%
  mutate(time = str_extract(sample, "[2468]"),
         ct = paste(caste, time, sep = "")) %>%
  group_by(Genes, ct) %>%
  summarise(mean_meth = mean(methylation)) %>%
  spread(ct, mean_meth) %>%
  mutate(diff_meth2 = queen2 - worker2,
         diff_meth4 = queen4 - worker4,
         diff_meth6 = queen6 - worker6,
         diff_meth8 = queen8 - worker8) %>%
  dplyr::select(Genes, starts_with("diff")) %>%
  gather("Time", "diff_meth", -Genes) %>%
  mutate(Time = as.numeric(str_extract(Time, "[2468]"))) %>%
  filter(Time == 8) %>% dplyr::select(Genes, diff_meth) %>%
  dplyr::rename(`Caste difference in % methylation in 8h larvae (positive: Q > W)` = diff_meth)

present_study_data <- present_study_data %>%
  left_join(bwasp_gene_meth,
            by = c("Gene symbol" = "Genes"))

# Merge with the present study's median % methylation of sites within each gene (from BWASP)
present_study_data <- present_study_data %>%
  left_join(read_tsv("data/meth_network_input/Network_persite.txt") %>%
              gather(key, value, -Genes) %>%
              filter(str_detect(key, "8")) %>%
              group_by(Genes) %>%
              summarise(mean_percent_meth_8h = mean(value) * 100),
            by = c("Gene symbol" = "Genes")
  ) 


# og version using the network methylation data:
# present_study_data <- present_study_data %>%
#   left_join(read_tsv("data/Methylation_data/input_data_for_meth_network.txt") %>%
#               gather(sample, percent_meth, -Genes) %>%
#               filter(str_detect(sample, "8")) %>%  # USES THE 8-HOUR CASTE DIFFERENCE RESULTS
#               group_by(Genes) %>%
#               summarise(percent_meth = median(percent_meth, na.rm = T), .groups = "drop") %>%
#               filter(!is.na(percent_meth)),
#             by = c("Gene symbol" = "Genes"))




# codon adaptation index measurements provided by Brendan Hunt
# merged_data <- read.delim("data/apis_gene_comparisons/Amel_AllData_012709.txt",
#                           header=T, stringsAsFactors = FALSE) %>%
#   # dplyr::select(ID, CAI) %>% # get gene ID and codon
#   # left_join(beebase_converter, by = c("ID" = "beebase1")) %>%
#   # dplyr::rename(beebase = beebase2) %>%
#   # left_join(gene_names, by = "beebase") %>%
#   dplyr::select(ID, CAI)

# Merge present study with:
# gamma data from Harpur et al 2014 PNAS
# connectivity data from the present study (transcription and methylation networks)
# Also add data on caste-specific gene expression from Warner et al 2018 Nature Comms
# And finally the pheromone sensitivity data for Apis from Holman et al. 2019 Nature Comms
merged_data <- present_study_data %>%
  
  # left_join(read.delim("data/apis_gene_comparisons/Amel_AllData_012709.txt",
  #                      header = TRUE, stringsAsFactors = FALSE) %>%
  #             dplyr::select(ID, CAI), by = c("beebase" = "ID"))

  # full_join(read_csv("data/apis_gene_comparisons/apis_gene_methyl_CG_OE.csv") %>% 
  #             dplyr::select(-NCBI_TranscriptID, -GC_OE) %>% distinct(gene, .keep_all = T),
  #           by = c("beebase" = "gene")) %>%
  
  full_join(read_delim("data/apis_gene_comparisons/harpur_etal_gamma.txt", delim = "\t"),
            by = c("beebase" = "Gene")) %>%
  
  full_join(read_csv("output/gene_connectivity.csv") %>% # gene co-expression network connectivity
              dplyr::select(Gene, kTotal) %>% dplyr::rename(expr_connectivity = kTotal),
            by = c("Gene symbol" = "Gene")) %>%
  full_join(read_csv("output/meth_connectivity.csv") %>% # gene co-methylation network connectivity
              dplyr::select(Gene, kTotal) %>% dplyr::rename(meth_connectivity = kTotal),
              by = c("Gene symbol" = "Gene")) %>%
  full_join(read_csv("data/apis_gene_comparisons/Warner2018_apis_caste_logFC.csv") %>%
              rename_all(~ paste(.x, " (Warner 2018)", sep = "")),
            by = c("Gene symbol" = "gene (Warner 2018)")) %>%
  full_join(read_csv("data/apis_gene_comparisons/pheromone_sensitivity.csv"),
            by = c("beebase" = "gene"))

# Calculate log2 CpG O/E ratio - change the sign, so that high values mean high methylation
# NB the expression level (AveExpr column) is already log2 transformed, see ?topTable
# merged_data <- merged_data %>% 
#   mutate(CG_OE = -log2(CG_OE),
#          CG_OE = replace(CG_OE, is.infinite(CG_OE), NA))       


# Add the data from Wojciechowski et al. 2018 Genome Biology
# There is gene-level data on the caste difference in 3 different histone modifications (from ChIP-seq),
# These data were created from Wojciechowski et al.'s raw data, as described by this script:
# https://github.com/lukeholman/queen-pheromone-RNAseq/blob/master/code/wojciechowski_histone_analysis.R
merged_data <- merged_data %>% 
  left_join(read_csv("data/apis_gene_comparisons/wojciechowski_histone_data/H3K4me3.csv") %>% 
              mutate(H3K4me3_caste = caste_difference) %>% 
              dplyr::select(gene, H3K4me3_caste), by =  c("beebase" = "gene")) %>%
  left_join(read_csv("data/apis_gene_comparisons/wojciechowski_histone_data/H3K27ac.csv") %>% 
              mutate(H3K27ac_caste = caste_difference) %>% 
              dplyr::select(gene, H3K27ac_caste), by =  c("beebase" = "gene")) %>%
  left_join(read_csv("data/apis_gene_comparisons/wojciechowski_histone_data/H3K36me3.csv") %>% 
              mutate(H3K36me3_caste = caste_difference) %>% 
              dplyr::select(gene, H3K36me3_caste), by =  c("beebase" = "gene"))


# Add data from Ma et al. 2019 BMC Genomics
# They measured the response of the transcriptome to two larva-produced pheromones. 
# One is called "brood pheromone" (it's 10 chemicals), one is called EBO: (E)-beta-ocimene.
# They also compared expression profiles of workers that chose to forage on pollen or nectar
# Positive number means higher expression after pheromone treatment
ma_brood_phero <- read_csv("data/apis_gene_comparisons/Ma_2019_brood_pheromone_RNAseq/BP_pheromone.csv") %>% 
  dplyr::select(gene, log2FoldChange) %>% 
  dplyr::rename(brood_pheromone = log2FoldChange) %>%
  left_join(read_csv("data/apis_gene_comparisons/Ma_2019_brood_pheromone_RNAseq/EBO_pheromone.csv") %>% 
              dplyr::select(gene, log2FoldChange) %>% 
              dplyr::rename(EBO_pheromone = log2FoldChange), by = "gene") %>%
  left_join(read_csv("data/apis_gene_comparisons/Ma_2019_brood_pheromone_RNAseq/pollen_nectar.csv") %>% 
              dplyr::select(gene, log2FoldChange) %>% 
              dplyr::rename(nectar_pollen = log2FoldChange), by = "gene")


merged_data <- merged_data %>% 
  left_join(ma_brood_phero %>% mutate(gene = as.character(gene)), 
            by = c("entrez_id" = "gene"))

m.rename <- function(merged, col, new) {
  if(!(col %in% names(merged))) print(col)
  names(merged)[names(merged) == col] <- new
  merged
}


merged_data <- merged_data %>%
#  m.rename("CAI", "Codon usage bias\n(CAI)") %>%
  m.rename("mean_percent_meth_8h", "% methylated cytosines in gene body") %>%
  m.rename("gamma", "Positive selection\n(Gamma)") %>%
  m.rename("AveExpr", "Log2 mean expression level") %>%
  m.rename("expr_connectivity", "Connectivity in the\ntranscriptome") %>%
  m.rename("meth_connectivity", "Connectivity in the\nmethylome") %>%
  m.rename("H3K4me3_caste", "Caste difference in H3K4me3 at 96h (positive: Q > W; Woj. 2018)") %>%
  m.rename("H3K27ac_caste", "Caste difference in H3K27ac at 96h (positive: Q > W; Woj. 2018)") %>%
  m.rename("H3K36me3_caste", "Caste difference in H3K36me3 at 96h (positive: Q > W; Woj. 2018)") %>%
  m.rename("Pheromone sensitivity", "Expression response to queen pheromone (Holman 2019)") %>%
  m.rename("brood_pheromone", "Expression response to brood pheromone (Ma 2019)") %>%
  m.rename("EBO_pheromone", "Expression response to EBO pheromone (Ma 2019)") %>%
  m.rename("nectar_pollen", "Expression difference between nectar and pollen foragers (Ma 2019)") %>%
  rename_all(~ str_replace_all(.x, "Warner_", "Expression up in "))

merged_data <- merged_data %>%
  dplyr::rename(gene = `Gene symbol`) %>%
  dplyr::select(gene, beebase, entrez_id, gene_name, everything()) %>%
  arrange(gene)

# Merge in the data from Flyatlas2, showing the sex difference in expression in whole bodies 
# for each gene in Drosophila (converted to bee orthologs using the ortholog relationships from Warner et al)
sex_diff <- read_csv("data/database_tables/flyatlas2_sex_diff.csv") %>%
  left_join(read_csv("data/database_tables/dros_ortho_GO.csv") %>% 
              dplyr::select(FLYBASE, SYMBOL) %>% distinct(), by = c("gene" = "FLYBASE")) %>%
  arrange(SYMBOL) %>% dplyr::select(-gene) %>%
  mutate(sex_bias = -1 * log(as.numeric(`M/F`))) %>%
  filter(Tissue == "Whole body" & !is.na(sex_bias)) %>%
  dplyr::rename(`Female bias in expression (Flyatlas2)` = sex_bias) %>%
  dplyr::select(SYMBOL, `Female bias in expression (Flyatlas2)`)

# Merge in the data from Flyatlas2, showing the enrichment of each gene in each tissue type, for female Drosophila
# (converted to bee orthologs using the ortholog relationships from Warner et al)
# Enrichment is expressed as FPKM_i / max(FPKM), where i is tissue i and max is the maximum FPKM recorded for all the tissues.
flyatlas <- read_csv("data/database_tables/flyatlas2_fpkm.csv") %>%
  left_join(read_csv("data/database_tables/dros_ortho_GO.csv") %>% 
              dplyr::select(FLYBASE, SYMBOL) %>% distinct(), by = c("gene" = "FLYBASE")) %>%
  arrange(SYMBOL) %>% dplyr::select(-gene) %>%
  mutate(FPKM = as.numeric(FPKM)) %>%
  filter(type == "Female" & !(Tissue %in% c("Whole body", "Carcass"))) %>%
  split(.$SYMBOL) %>%
  map_df(~ .x %>% mutate(FPKM = FPKM / max(FPKM, na.rm = TRUE))) %>%
  dplyr::select(SYMBOL, Tissue, FPKM) %>%
  spread(Tissue, FPKM) %>%
  select_if(~ sum(!is.na(.x)) > 0) %>%
  rename_all(~ paste(.x, "(Flyatlas2)"))


merged_data <- merged_data %>%
  left_join(sex_diff, by = c("gene" = "SYMBOL")) %>%
  left_join(flyatlas, by = c("gene" = "SYMBOL (Flyatlas2)")) 


# Remove this variable - too many missing values
# merged_data <- merged_data %>%
#   select(-`DNA methylation frequency\n(CpG depletion)`)

rm(list = c("present_study_data", 
            "ma_brood_phero",
            "gene_names", "m.rename"))

