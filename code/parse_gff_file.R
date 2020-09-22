# Parse all the interesting parts of the Apis mellifera .gff file into more R-friendly formats (takes overnight, lazy coding)

# Get all the unique sites in the methylation analysis (produced in download_methyl_googledrive.R)
sites <- read_tsv("data/methyl_diff_results/sites.tsv") %>%
  mutate(is_promoter=F, is_3UTR=F, is_5UTR=F, is_exon=F, 
         is_CDS=F, gene=F, beebase=F, gene_biotype=F)

# Get promoter annotations and add to the sites file
promoters <- read_tsv("data/Genome_annotation/Amel.promoter.gff3", 
         col_names = c("seqnames", "platform", "is_promoter", "start", "end", "dot", "strand", "dot2", "details")) %>%
  select("seqnames", "start", "end", "is_promoter")

for(i in 1:nrow(promoters)){
  sites$is_promoter[sites$seqnames==promoters$seqnames[i] &
                      sites$start>= promoters$start[i] &
                      sites$end<= promoters$end[i]] <- TRUE
}

#######################################################################
# Get 3-UTR annotations and add to the sites file
x3UTR <- read_tsv("data/Genome_annotation/Amel.pcg-3pUTR.gff3", 
         col_names = c("seqnames", "dot", "is_3UTR", "start", "end", "dot2", "strand", "dot3", "details")) %>%
  select("seqnames", "start", "end", "is_3UTR")

for(i in 1:nrow(x3UTR)){
  sites$is_3UTR[sites$seqnames==x3UTR$seqnames[i] &
                      sites$start>= x3UTR$start[i] &
                      sites$end<= x3UTR$end[i]] <- TRUE
}

#######################################################################
# Get 5-UTR annotations and add to the sites file
x5UTR <- read_tsv("data/Genome_annotation/Amel.pcg-5pUTR.gff3", 
                  col_names = c("seqnames", "dot", "is_5UTR", "start", "end", "dot2", "strand", "dot3", "details")) %>%
  select("seqnames", "start", "end", "is_5UTR")

for(i in 1:nrow(x5UTR)){
  sites$is_5UTR[sites$seqnames==x5UTR$seqnames[i] &
                  sites$start>= x5UTR$start[i] &
                  sites$end<= x5UTR$end[i]] <- TRUE
}

#######################################################################
# Get exon annotations and add to the sites file
exon <- read_tsv("data/Genome_annotation/Amel.exon.gff3", 
                 col_names = c("seqnames", "platform", "is_exon", "start", "end", "dot2", "strand", "dot3", "details")) %>%
  select("seqnames", "start", "end", "is_exon")

for(i in 1:nrow(exon)){
  sites$is_exon[sites$seqnames==exon$seqnames[i] &
                  sites$start>= exon$start[i] &
                  sites$end<= exon$end[i]] <- TRUE
}

#######################################################################
# Get coding sequences (CDS) annotations and add to the sites file
CDS <- read_tsv("data/Genome_annotation/Amel.pcg-CDS.gff3", 
                 col_names = c("seqnames", "platform", "is_CDS", "start", "end", "dot2", "strand", "dot3", "details")) %>%
  select("seqnames", "start", "end", "is_CDS")

for(i in 1:nrow(CDS)){
  sites$is_CDS[sites$seqnames==CDS$seqnames[i] &
                  sites$start>= CDS$start[i] &
                  sites$end<= CDS$end[i]] <- TRUE
}

#######################################################################
# Get gene annotations and add to the sites file
genes <- read_tsv("data/Genome_annotation/Amel.gene.gff3", 
         col_names = c("seqnames", "platform", "gene", "start", "end", "dot2", "strand", "dot3", "details"))[1:1000,] %>%
  select("seqnames", "start", "end", "details") %>%
  mutate(details = str_split(details, ";"),
         gene = str_remove_all(map_chr(details, ~ .x[1]), "ID=gene-"),
         gene_biotype = str_remove_all(map_chr(details, ~ .x[str_detect(.x, "biotype")]), "gene_biotype=")) %>%
  select(-details)

for(i in 1:nrow(genes)){
  rows <- which(sites$seqnames == genes$seqnames[i] &
                 sites$start >= genes$start[i] &
                 sites$end <= genes$end[i])
  sites$gene[rows] <- genes$gene[i]
  sites$beebase[rows] <- genes$beebase[i]
  sites$gene_biotype[rows] <- genes$gene_biotype[i]
}

#######################################################################
# To save space, make another spreadsheet with all the gene info (two name types, and "gene biotype") 
# Create a gene_id column for later spreadsheet joins
gene_annots <- sites %>% select(gene, beebase, gene_biotype) %>% 
  distinct() %>% mutate(gene_id = 1:n()) %>% select(gene_id, everything())

# Add gene_id and then remove the 3 gene columns
sites <- sites %>% 
  left_join(gene_annots, by = c("gene", "beebase", "gene_biotype")) %>% select(-gene, -beebase, -gene_biotype)

# Save both to disk
write_tsv(sites, "data/site_annotations_from_gff.tsv")
write_tsv(gene_annots, "data/database_tables/gene_annotations_from_gff.tsv")

