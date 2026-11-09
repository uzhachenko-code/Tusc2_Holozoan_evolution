# ==============================================================================
# TUSC2 C-TERMINAL PDZ-BINDING EVOLUTIONARY ANALYSIS & HIGH-SPEED PERMUTATIONS
# Consensus: ...[VLIFY].[ACVILF]$ | Fully Corrected Indexing & Shuffling Logic
# ==============================================================================
rm(list = ls())
graphics.off()

library(dplyr)
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
library(tidyr)
library(ggplot2)

fasta_path <- "C:/Users/.../datamonkey_fused_flat_protein.fasta"
taxon_path <- "C:/Users/.../Taxon_ID_categorized.txt"
output_dir <- "C:/Users/..."

# 1. Load Taxonomy Metadata Table & Standardize Group Levels
tax_meta <- read.delim(taxon_path, header = TRUE, sep = "\t", strip.white = TRUE)
tax_meta$Flat_ID <- trimws(as.character(tax_meta$Flat_ID))
tax_meta$Group   <- trimws(as.character(tax_meta$Group))
tax_meta$Is_Mammal <- if_else(tax_meta$Group == "Mammal", "Mammal", "Non-Mammal")

# 2. Parse FASTA Alignment Lines
fasta_lines <- readLines(fasta_path)
headers <- c(); sequences <- c(); current_seq <- ""
for (line in fasta_lines) {
  line <- trimws(line)
  if (startsWith(line, ">")) {
    if (current_seq != "") { sequences <- c(sequences, current_seq); current_seq <- "" }
    headers <- c(headers, substring(line, 2))
  } else { current_seq <- paste0(current_seq, line) }
}
if (current_seq != "") { sequences <- c(sequences, current_seq) }

seq_matrix <- do.call(rbind, strsplit(sequences, ""))
total_cols <- ncol(seq_matrix)

# 3. Robust, Accurate PDZ Evaluator Function
# Correctly maps: Position -2 -> second to last | Position 0 -> final residue
evaluate_pdz_accurate <- function(mat) {
  n_rows <- nrow(mat)
  pdz_valid <- logical(n_rows)
  
  pos_minus_2 <- c("V", "L", "I", "F", "Y")
  pos_0       <- c("A", "C", "V", "I", "L", "F")
  
  for (r in 1:n_rows) {
    # Extract real amino acids by removing alignment gaps
    true_seq <- mat[r, mat[r, ] != "-"]
    len <- length(true_seq)
    
    if (len >= 6) {
      # Target the specific positions dictated by the consensus pattern
      res_minus_2 <- true_seq[len - 2]
      res_0       <- true_seq[len]
      
      if (res_minus_2 %in% pos_minus_2 && res_0 %in% pos_0) {
        pdz_valid[r] <- TRUE
      }
    }
  }
  return(pdz_valid)
}

# 4. Map Biological Alignment via Exact Token Matching
extract_id <- function(header_str) {
  match_pos <- regexpr("Taxon_\\d+", header_str, ignore.case = TRUE)
  if (match_pos != -1) return(regmatches(header_str, match_pos))
  return(NA)
}
clean_fasta_ids <- sapply(headers, extract_id)
matched_indices <- match(clean_fasta_ids, tax_meta$Flat_ID)

master_db <- data.frame(
  Species   = tax_meta$Species[matched_indices],
  Group     = tax_meta$Group[matched_indices],
  Is_Mammal = tax_meta$Is_Mammal[matched_indices],
  stringsAsFactors = FALSE
)

# Compute the accurate biological baseline
master_db$PDZ_Valid <- evaluate_pdz_accurate(seq_matrix)
master_db           <- na.omit(as.data.frame(master_db))

# Isolate Mammalian rows for permutation testing
mammal_row_indices <- which(master_db$Is_Mammal == "Mammal")
obs_mammal_hits    <- sum(master_db$PDZ_Valid[mammal_row_indices])
mammal_submatrix   = seq_matrix[mammal_row_indices, , drop = FALSE]

# 5. Run 10,000 Permutations (Accurate dynamic tracking)
cat("Executing 10,000 permutations with strict terminal tracking...\n")
set.seed(42)
n_permutations <- 10000
null_distribution <- numeric(n_permutations)

pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)
for (p in 1:n_permutations) {
  # Shuffle alignment columns globally
  shuffled_cols <- mammal_submatrix[, sample(total_cols), drop = FALSE]
  
  # Accurately compute hits based on the newly formed shuffled C-termini
  perm_eval <- evaluate_pdz_accurate(shuffled_cols)
  null_distribution[p] <- sum(perm_eval)
  
  setTxtProgressBar(pb, p)
}
close(pb)

empirical_p_val <- sum(null_distribution >= obs_mammal_hits) / n_permutations

# 6. Compute Mammalia Specificity Hypothesis (Fisher's Exact 2x2)
mammal_factor <- factor(master_db$Is_Mammal, levels = c("Non-Mammal", "Mammal"))
pdz_factor    <- factor(master_db$PDZ_Valid, levels = c(FALSE, TRUE))
fisher_table  <- table(mammal_factor, pdz_factor)
fisher_p_val  <- format.pval(fisher.test(fisher_table, alternative = "greater")$p.value, digits = 5)

# 7. Print Master Publications Report
cat("\n========================================================\n")
cat("          TUSC2 C-TERMINAL PDZ-BINDING EVOLUTIONARY REPORT \n")
cat("========================================================\n\n")

rates_by_level <- master_db %>%
  group_by(Group) %>%
  summarise(Total = n(), Hits = sum(PDZ_Valid), Rate = round((Hits/Total)*100, 1), .groups='drop')

print(as.data.frame(rates_by_level))
cat("\n--- Contingency Profiling ---\n")
print(fisher_table)
cat("\nFisher's Exact Test for Mammalia Enrichment (P-value): ", fisher_p_val, "\n")
cat("Empirical Random Shuffling Permutation P-value (10k):   ", empirical_p_val, "\n\n")

# ==============================================================================
# EXPORT EXTENDED SUPPLEMENTARY TABLE FOR TUSC2 C-TERMINAL PDZ MOTIF SCREENING
# Captures literal sequence variations and consensus rules for all species
# ==============================================================================

# 1. Dynamically extract the literal final 6-AA sequence for each species row
extracted_fragments <- character(nrow(seq_matrix))

for (r in 1:nrow(seq_matrix)) {
  true_seq <- seq_matrix[r, seq_matrix[r, ] != "-"]
  len <- length(true_seq)
  if (len >= 6) {
    extracted_fragments[r] <- paste(true_seq[(len-5):len], collapse = "")
  } else {
    extracted_fragments[r] <- paste(true_seq, collapse = "") # Fallback if sequence is extremely short
  }
}

# 2. Append the sequence fragment string back to your master database tracking
master_db$C_Terminal_6AA <- extracted_fragments

# 3. Compile and structure the expanded publication-ready table
supp_table <- master_db %>%
  dplyr::select(
    `Species Name`               = Species,
    `Evolutionary Clade`         = Group,
    `C-Terminal Sequence (6-AA)` = C_Terminal_6AA,
    `PDZ Motif Status`           = PDZ_Valid
  ) %>%
  mutate(
    `Target Consensus Rule`     = "...[VLIFY].[ACVILF]$",
    `PDZ Motif Status`           = if_else(`PDZ Motif Status`, "TRUE (Detected)", "FALSE (Absent)")
  ) %>%
  dplyr::select(
    `Species Name`, 
    `Evolutionary Clade`, 
    `Target Consensus Rule`, 
    `C-Terminal Sequence (6-AA)`, 
    `PDZ Motif Status`
  ) %>%
  arrange(`Evolutionary Clade`, `Species Name`)

# 4. Save directly to your project space
write.csv(
  supp_table, 
  file = file.path(output_dir, "Supplementary_Table_S3_PDZ_Motif_Screening.csv"), 
  row.names = FALSE
)

cat("Extended supplementary table with sequence fragments successfully exported!\n")
