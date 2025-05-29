library(Biostrings)

# Function to calculate overall GC content
calculate_all_gc <- function(sequences_matrix) {
  gc_count <- sum(sequences_matrix == "G" | sequences_matrix == "C")
  total_count <- length(sequences_matrix)
  gc_content <- gc_count / total_count
  return(gc_content)
}

fasta_files <- list.files(path = "./Tree/", pattern = "\\_aligned.fasta$", full.names = TRUE)

gc_results <- lapply(fasta_files, function(file) {
  dna <- readDNAStringSet(file)
  
  all_chars <- unlist(strsplit(as.character(unlist(dna)), split = ""))
  
  gc <- calculate_all_gc(all_chars)
  
  data.frame(filename = basename(file), gc_content = gc)
})

HCV_gc_df <- do.call(rbind, gc_results)

### Calculate similarity using the pegas package
library(ape)
library(pegas)

fasta_files <- list.files(path = "./Tree/", pattern = "\\_aligned.fasta$", full.names = TRUE)

pi_results <- lapply(fasta_files, function(file) {
  
  dna <- read.dna(file, format = "fasta")
  
  pi <- tryCatch({
    nuc.div(dna)
  }, error = function(e) NA)  # Handle errors (e.g., if file is too small or not aligned)
  
  data.frame(filename = basename(file), pi = pi)
})

# Combine into a single data frame
HCV_sim_df <- do.call(rbind, pi_results)