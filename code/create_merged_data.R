# library(tidyverse)
# library(DBI)
# library(RSQLite)
source("code/parse_warner_data.R")

gene_names <- tbl(db, "gene_names") %>% collect() 
x <- read_tsv("data/HaV3.1_gene_names.txt")

# beebase_converter <- read_csv("data/database_tables_new/beebase1_to_beebase2.csv") 


# Get the present study's limma::voom results, and associate with their beebase IDs for matching the older studies
present_study_data <- read_csv("data/caste_results.csv") %>% 
  dplyr::filter(Time == 8) %>%              # USES THE 8-HOUR CASTE DIFFERENCE AND EXPRESSION LEVEL RESULTS
  left_join(tbl(db, "gene_names") %>% 
              dplyr::select(gene_symbol, gene_name, gene, entrez_id, beebase) %>% 
              collect(n=Inf), 
            by = c("Gene symbol" = "gene_symbol")) %>% 
  dplyr::select(beebase, logFC, AveExpr, `Gene symbol`, entrez_id, gene, gene_name) %>%
  dplyr::rename(`Upregulation in queen-destined 8h larvae` = logFC) 

# codon adaptation index measurements provided by Brendan Hunt
# merged_data <- read.delim("data/apis_gene_comparisons/Amel_AllData_012709.txt", 
#                           header=T, stringsAsFactors = FALSE) %>% 
#   # dplyr::select(ID, CAI) %>% # get gene ID and codon 
#   # left_join(beebase_converter, by = c("ID" = "beebase1")) %>%
#   # dplyr::rename(beebase = beebase2) %>%                              
#   # left_join(gene_names, by = "beebase") %>% 
#   dplyr::select(ID, CAI)

# Merge present study with CAI data, then
# Add in the methylation data provided by Soojin Yi and Xin Wu (the data are from Galbraith et al PNAS)
# Also add the gamma data from Harpur et al 2014 PNAS
# Also add the connectivity data from the present study
# Also add data on caste-specific gene expression from Warner et al 2018 Nature Comms
# And finally the pheromone sensitivity data for Apis from Holman et al. 2019 Nature Comms
merged_data <- present_study_data %>%

  full_join(read_csv("data/apis_gene_comparisons/apis_gene_methyl_CG_OE.csv") %>% 
              dplyr::select(-NCBI_TranscriptID, -GC_OE) %>% distinct(gene, .keep_all = T),
            by = c("beebase" = "gene")) %>%
  
  full_join(read_delim("data/apis_gene_comparisons/harpur_etal_gamma.txt", delim = "\t"),
            by = c("beebase" = "Gene")) %>%
  full_join(read_csv("data/gene_connectivity.csv") %>% dplyr::select(Gene, kTotal),
            by = c("Gene symbol" = "Gene")) %>%
  full_join(read_csv("data/apis_gene_comparisons/Warner2018_apis_caste_logFC.csv") %>%
              rename_all(~ paste("Warner", .x, sep = "_")),
            by = c("Gene symbol" = "Warner_gene")) %>%
  full_join(read_csv("data/apis_gene_comparisons/pheromone_sensitivity.csv"),
            by = c("beebase" = "gene"))

# Calculate log2 CpG O/E ratio - change the sign, so that high values mean high methylation
# NB the expression level (AveExpr column) is already log2 transformed, see ?topTable
merged_data <- merged_data %>% 
  mutate(CG_OE = -log2(CG_OE),
         CG_OE = replace(CG_OE, is.infinite(CG_OE), NA))       


# Add the data from Wojciechowski et al. 2018 Genome Biology
# There is gene-level data on the caste difference in 3 different histone modifications (from ChIP-seq),
# These data were created from Wojciechowski et al.'s raw data, as described by this script:
# https://github.com/lukeholman/queen-pheromone-RNAseq/blob/master/code/wojciechowski_histone_analysis.R
merged_data <- merged_data %>% 
  left_join(read_csv("data/apis_gene_comparisons/wojciechowski_histone_data/H3K4me3.csv") %>% 
              mutate(H3K4me3_caste = caste_difference) %>% dplyr::select(gene, H3K4me3_caste), by =  c("beebase" = "gene")) %>%
  left_join(read_csv("data/apis_gene_comparisons/wojciechowski_histone_data/H3K27ac.csv") %>% 
              mutate(H3K27ac_caste = caste_difference) %>% dplyr::select(gene, H3K27ac_caste), by =  c("beebase" = "gene")) %>%
  left_join(read_csv("data/apis_gene_comparisons/wojciechowski_histone_data/H3K36me3.csv") %>% 
              mutate(H3K36me3_caste = caste_difference) %>% dplyr::select(gene, H3K36me3_caste), by =  c("beebase" = "gene")) %>%
  dplyr::rename(`Gene name` = gene_name,
         `Beebase ID` = gene)


# Add data from Ma et al. 2019 BMC Genomics
# They measured the response of the transcriptome to two larva-produced pheromones. 
# One is called "brood pheromone" (it's 10 chemicals), one is called EBO: (E)-beta-ocimene.
# They also compared expression profiles of workers that chose to forage on pollen or nectar
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
  m.rename("CAI", "Codon usage bias\n(CAI)") %>%
  m.rename("CG_OE", "DNA methylation frequency\n(CpG depletion)") %>%
  m.rename("Gene_body_methylation", "DNA methylation frequency\n(BiS-seq)") %>%
  m.rename("gamma", "Positive selection\n(Gamma)") %>%
  m.rename("AveExpr", "Log2 mean expression level") %>%
  m.rename("kTotal", "Connectivity in the\ntranscriptome") %>%
  m.rename("H3K4me3_caste", "Caste difference in\nH3K4me3 at 96h") %>%
  m.rename("H3K27ac_caste", "Caste difference in\nH3K27ac at 96h") %>%
  m.rename("H3K36me3_caste", "Caste difference in\nH3K36me3 at 96h") %>%
  m.rename("Pheromone sensitivity", "Upregulation in workers exposed\nto queen pheromone") %>%
  m.rename("brood_pheromone", "Upregulation in workers exposed\nto brood pheromone") %>%
  m.rename("EBO_pheromone", "Upregulation in workers exposed\nto EBO pheromone") %>%
  m.rename("nectar_pollen", "Upregulation in nectar-foraging workers")

merged_data <- merged_data %>%
  # dplyr::select(-`Beebase ID`, -`Gene name`, -beebase, -NCBI_GeneID, -entrez_id) %>%
  dplyr::rename(gene = `Gene symbol`) %>%
  dplyr::select(gene, beebase, NCBI_GeneID, entrez_id, `Gene name`, everything())

# Remove this variable - too many missing values
# merged_data <- merged_data %>%
#   select(-`DNA methylation frequency\n(CpG depletion)`)

rm(list = c("present_study_data", 
            "ma_brood_phero",
            "gene_names", "m.rename"))

