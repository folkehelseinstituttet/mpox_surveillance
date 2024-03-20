library(tidyverse)

# Read the ncbi metadatafile
ncbi <- read_tsv("/home/jonr/Prosjekter/mpox/phylogenetic/data/metadata.tsv.gz")

# Read the consensus sequences using seqinr
fasta_files <- list.files("/home/jonr/Prosjekter/mpox/phylogenetic/data/", ".fa$", full.names = T)

for (i in seq_along(fasta_files)) {
  seq <- seqinr::read.fasta(fasta_files[i])
  if (!exists("fasta")) {
    fasta <- seq
  } else {
    fasta <- c(fasta, seq)
  }
}

# Simplify the fasta headers and write out as a single fasta
seqinr::write.fasta(fasta, names = names(fasta), file.out = "/home/jonr/Prosjekter/mpox/phylogenetic/data/norw_seqs.fa")
# This file will be concatenated with the ncbi fasta file

# Create the metadata with the same headers as the ncbi metadatafile
fields <- colnames(ncbi)

metadata <- matrix(nrow = length(fasta), ncol = length(fields))
colnames(metadata) <- fields
metadata <- as.data.frame(metadata)

# Add the fasta header as strain name
for (i in seq_along(fasta)) {
  metadata$strain[i] <- names(fasta[i])
  metadata$accession[i] <- names(fasta[i])
  metadata$country[i] <- "Norway"
  metadata$region[i] <- "Europe"
  metadata$host[i] <- "Homo sapiens"
  metadata$QC_rare_mutations <- "good"
}

metadata <- tibble(metadata)

# Add sampling dates
metadata <- metadata %>% 
  mutate(date = case_when(
    strain == "2023122020" ~ "2022-06-22",
    strain == "202401196" ~ "2023-10-05",
    strain == "202401262" ~ "2022-06-22",
    strain == "202401263" ~ "2022-06-27",
    strain == "202401264" ~ "2022-07-25",
    strain == "202401265" ~ "2022-07-29",
    strain == "202401266" ~ "2022-06-22",
    strain == "202401267" ~ "2023-11-29",
    .default = date
  )) %>%  
  mutate(date_submitted = case_when(
    strain == "2023122020" ~ "2022-06-22",
    strain == "202401196" ~ "2023-10-05",
    strain == "202401262" ~ "2022-06-22",
    strain == "202401263" ~ "2022-06-27",
    strain == "202401264" ~ "2022-07-25",
    strain == "202401265" ~ "2022-07-29",
    strain == "202401266" ~ "2022-06-22",
    strain == "202401267" ~ "2023-11-29",
    .default = date_submitted
  )) 

# Read the nextclade results from the Norwegian sequences
nc <- read_tsv("/home/jonr/Prosjekter/mpox/phylogenetic/data/nextclade.tsv") %>% 
  mutate(seqName = as.character(seqName)) %>% 
  # keep relevant fields
  select(seqName, lineage, outbreak, clade, coverage, isReverseComplement, frameShifts)

# merge nextclade data
metadata <- metadata %>% 
  left_join(nc, by = c("strain" = "seqName"))

# Create final data structure
metadata <- metadata %>% 
  select(accession,
         genbank_accession_rev, 
         strain, 
         date,
         region,
         country,
         division,
         location,
         host,
         date_submitted,
         sra_accession,
         abbr_authors,
         reverse,
         authors,
         institution,
         "clade" = clade.y, # Use clade from Nextclade
         "outbreak" = outbreak.y, # Use outbreak from Nextclade
         "lineage" = lineage.y, # Use lineage from Nextclade
         "coverage" = coverage.y, # Use coverage from Nextclade
         missing_data,
         divergence,
         nonACGTN,
         QC_missing_data,
         QC_mixed_sites,
         QC_rare_mutations,
         QC_frame_shifts,
         QC_stop_codons,
         "frame_shifts" = frameShifts, # Use frame shifts from Nextclade
         "is_reverse_complement" = isReverseComplement) # Use reverse complement from Nextclade

metadata <- metadata %>% 
  mutate(date_submitted = as.Date(date_submitted))

# Merge own data with ncbi and write file
# Test if column structures are the same
if (identical(colnames(metadata), colnames(ncbi))) {
 final_metadata <- bind_rows(metadata, ncbi) 
} else {
  print("Colnames are not the same")
}

if (exists("final_metadata")) {
  write_tsv(final_metadata, "/home/jonr/Prosjekter/mpox/phylogenetic/data/metadata.tsv")
}


final_metadata <- final_metadata %>% 
  mutate(date_submitted = case_when(
    strain == "2023122020" ~ "2022-06-22",
    strain == "202401196" ~ "2023-10-05",
    strain == "202401262" ~ "2022-06-22",
    strain == "202401263" ~ "2022-06-27",
    strain == "202401264" ~ "2022-07-25",
    strain == "202401265" ~ "2022-07-29",
    strain == "202401266" ~ "2022-06-22",
    strain == "202401267" ~ "2023-11-29",
    .default = date_submitted
  )) 
