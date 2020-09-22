library(tidyverse)
library(googledrive)
library(glue)
library(vroom)

# Search for the right files on Drive
search_hit <- drive_ls("Methylation_August_2020")
comparisons <- search_hit %>% pull(name) %>% str_remove_all("[.]txt")
file_ids <- search_hit %>% pull(id)

# Download and tidy up each of the CT count files on Google Drive
lapply(1:length(file_ids), function(i){
  print(i)
  
  filename <- glue("data/methyl_diff_results/{comparisons[i]}.txt")
  as_id(file_ids[i]) %>% 
    drive_download(path = filename, overwrite = TRUE)
  
  # Load up the file, and add a column with the sample name
  foc <- read_tsv(filename) %>% 
    mutate(sample = comparisons[i]) 
  
  # fix any weird names (sometime it is named e.g. numCs1 or numCs2, presumably a bug in Claire's code)
  names(foc)[str_detect(names(foc), "numC")] <- "numC" 
  names(foc)[str_detect(names(foc), "numT")] <- "numT"
  
  # Reorder columns and save to disk again
  foc %>% 
    dplyr::select(sample, seqnames, start, numC, numT) %>%
    write_tsv(filename)
})

  

# Get the list of sites (seqname, start) out of each file 
# **(not actually unique yet!!)**
unique_sites <- list.files("data/methyl_diff_results", 
                           full.names = T, pattern = "txt") %>%
  map_df(~ vroom(.x, col_select = c(seqnames, start)) %>%
           distinct() %>%
           mutate(pasted_site = paste(seqnames, start, sep = "~")) %>%
           select(pasted_site)
  )

# Count the number of samples for which each site was measured
site_sample_counts <- unique_sites %>%
  group_by(pasted_site) %>%
  summarise(n_samples = n())

# Save the counts of how many sites are represented by n samples, e.g. there are 5198696 sites covered by all 36 samples
table(site_sample_counts$n_samples) %>% 
  enframe("n_samples", "n_sites") %>%
  write_tsv("output/site_sample_counts.tsv")

# Get sites covered by at least 34/36 samples
sites_to_keep <- site_sample_counts %>%
  filter(site_sample_counts >= 34) %>%
  pull(pasted_site)


# make the sites unique and give each an ID number
unique_sites <- unique_sites %>% 
  distinct() %>%
  filter(pasted_site %in% sites_to_keep) %>%
  mutate(site = 1:n())


# Concatenate all the files, and replace the site details with id number to save disk space
# Also, filter to only keep sites covered by at least 34/36 samples
methyl_diff_results <- list.files("data/methyl_diff_results", full.names = T, pattern = "txt") %>%
  map_df(~ vroom::vroom(.x) %>%
           distinct() %>%
           mutate(pasted_site = paste(seqnames, start, sep = "~")) %>%
           left_join(unique_sites, by = "pasted_site") %>%
           filter(!is.na(site)) %>%
           select(site, sample, starts_with("num"))
         ) 

# Save the massive CT file
methyl_diff_results %>% 
  write_tsv("data/database_tables/methyl_diff_results.tsv")


# Calculate the % methylation per site (Across all samples), and find the SD and range in prop_meth
# This is useful for finding sites
site_info <- methyl_diff_results %>%
  mutate(prop = numC / (numC + numT)) %>%
  group_by(site) %>%
  summarise(
    total_C_reads = sum(numC),
    total_T_reads = sum(numT),
    n_samples = n(),
    mean_prop_meth = total_C_reads / (total_C_reads + total_T_reads),
    sd_prop_meth = sd(prop),
    range_prop_meth = max(prop) - min(prop)) 

# Save the site details too
unique_sites %>% 
  mutate(split = strsplit(pasted_site, split = "~"),
         seqnames = map_chr(split, ~ .x[1]),
         start = map_chr(split, ~ .x[2])) %>%
  select(site, seqnames, start) %>% 
  left_join(site_info, by = "site") %>%
  write_tsv("data/database_tables/site_info.tsv")




# Delete the constituent files that were downloaded from Drive
if(file.exists("data/database_tables/methyl_diff_results.tsv")){
  unlink(list.files("data/methyl_diff_results", full.names = T, pattern = "txt"))
}














