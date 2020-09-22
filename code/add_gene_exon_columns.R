# code to use the gff file to add some annotations to a data frame of seqnames and start points
# Useful e.g. for annotating gene names to a collection of methylated sites or 500bp windows of the genome

annot_db <- makeTxDbFromGFF(file = "data/Genome_annotation/Amel.pcg.gff3")

# Make a genes-to-transcripts mapping table
all_transcribed_regions <- exons(annot_db, columns = c("GENEID", "EXONRANK"))
tr <- transcriptsBy(annot_db, by=c("gene")) %>% unlist() 
genes_transcripts <- tibble(gene = tr@ranges@NAMES,
                            transcript = tr@elementMetadata@listData$tx_name)

# Find the genomic ranges for the introns in each transcript (and list gene for each transcript)
introns <- intronsByTranscript(annot_db, use.names=TRUE) %>% unlist()
rown <- introns@ranges@NAMES
introns <- introns %>% 
  as.data.frame(row.names = NULL) %>% mutate(transcript = rown) %>% 
  left_join(genes_transcripts, by = "transcript") %>% as("GRanges")


# Find the genomic ranges for the 5-UTRs (and list gene for each one)
fiveUTRs <- fiveUTRsByTranscript(annot_db, use.names=TRUE) %>% unlist()
rown <- fiveUTRs@ranges@NAMES
fiveUTRs <- fiveUTRs %>% 
  as.data.frame(row.names = NULL) %>% mutate(transcript = rown) %>% 
  left_join(genes_transcripts, by = "transcript") %>% as("GRanges")


# Find the genomic ranges for the 3-UTRs (and list gene for each one)
threeUTRs <- threeUTRsByTranscript(annot_db, use.names=TRUE) %>% unlist()
rown <- threeUTRs@ranges@NAMES
threeUTRs <- threeUTRs %>% 
  as.data.frame(row.names = NULL) %>% mutate(transcript = rown) %>% 
  left_join(genes_transcripts, by = "transcript") %>% as("GRanges")

# function to get the annotations out of the gff-database for a set of seqnames-start-end data (end column is optional)
add_gene_exon_columns <- function(df){
  
  # First, set up a GRanges object of the focal sites in df
  
  if("end" %in% names(df)){ # Assuming there is an "end" column...
    
    unique_sites <- df %>% 
      dplyr::select(seqnames, start, end) %>% 
      distinct()
    
    my_ranges <- GRanges(
      seqnames = pull(unique_sites, seqnames),
      ranges = IRanges(start = pull(unique_sites, start),
                       end = pull(unique_sites, end)))
  } else { #Otherwise, start = end (i.e. for single nucleotides)
    
    unique_sites <- df %>% 
      dplyr::select(seqnames, start) %>% 
      distinct()
    
    my_ranges <- GRanges(
      seqnames = pull(unique_sites, seqnames),
      ranges = IRanges(start = pull(unique_sites, start),
                       end = pull(unique_sites, start))) 
  }
  
  
  # Get the gene+exon names for each transcribed region (i.e. the protein coding sequence plus 5' and 3'-UTRs)
  # Note: the unnesting step means that each site may have more than one gene/exon combination associated with it
  ex <- exons(annot_db, columns = c("GENEID", "EXONRANK"))
  overlaps <- findOverlaps(my_ranges, ex) 
  
  ex <- as.data.frame(ex[overlaps@to, ]) %>% 
    dplyr::select(GENEID, EXONRANK) %>% 
    dplyr::rename(exon_in = GENEID) %>%
    dplyr::rename(exon_rank = EXONRANK) %>%
    bind_cols(unique_sites[overlaps@from, ]) %>%
    unnest("exon_in") %>% unnest("exon_rank") 
  
  # Now do promoter regions
  pr <- promoters(annot_db, columns = c("GENEID"), upstream=2000, downstream=200)
  overlaps <- findOverlaps(my_ranges, pr) 
  
  pr <- as.data.frame(pr[overlaps@to, ], row.names = NULL) %>% 
    dplyr::select(GENEID) %>% dplyr::rename(promoter_of = GENEID) %>%
    bind_cols(unique_sites[overlaps@from, ]) %>%
    unnest("promoter_of")
  
  # Now introns:
  overlaps <- findOverlaps(my_ranges, introns) 
  
  intr <- as.data.frame(introns[overlaps@to, ], row.names = NULL) %>% 
    dplyr::select(gene) %>% dplyr::rename(intron_in = gene) %>%
    bind_cols(unique_sites[overlaps@from, ]) %>%
    unnest("intron_in")
  
  # Now 5-UTRs
  overlaps <- findOverlaps(my_ranges, fiveUTRs) 
  
  five <- as.data.frame(fiveUTRs[overlaps@to, ], row.names = NULL) %>% 
    dplyr::select(gene) %>% dplyr::rename(fiveUTR_in = gene) %>%
    bind_cols(unique_sites[overlaps@from, ]) %>%
    unnest("fiveUTR_in")
  
  # Now 3-UTRs
  overlaps <- findOverlaps(my_ranges, threeUTRs) 
  
  three <- as.data.frame(threeUTRs[overlaps@to, ], row.names = NULL) %>% 
    dplyr::select(gene) %>% dplyr::rename(threeUTR_in = gene) %>%
    bind_cols(unique_sites[overlaps@from, ]) %>%
    unnest("threeUTR_in")
  
  df %>% 
    full_join(ex) %>% full_join(pr) %>% 
    full_join(intr) %>% full_join(five) %>% full_join(three) %>%
    mutate(translated = !is.na(exon_in) & is.na(fiveUTR_in) & is.na(threeUTR_in))
}