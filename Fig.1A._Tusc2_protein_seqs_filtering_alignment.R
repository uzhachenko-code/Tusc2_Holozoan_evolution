library(Biostrings)
library(msa)

# 1. Read FASTA file
aa_tusc2_sequences <- readAAStringSet("C:/Users/.../Protein_seq_Tusc2_unfiltered.fasta")

# 2. Filter out proteins over 200aa threshold
max_size <- 200 
filtered_aa_tusc2_seqs <- aa_tusc2_sequences[width(aa_tusc2_sequences) <= max_size]

# 3. Extract the very first amino acid of each sequence
first_amino_acid <- subseq(filtered_aa_tusc2_seqs, start = 1, end = 1)

# 4. Keep only sequences where the first amino acid is "M"
filtered_aa_tusc2_seqs <- filtered_aa_tusc2_seqs[as.character(first_amino_acid) == "M"]

# 5. Save the final alignment to a new FASTA file
writeXStringSet(filtered_aa_tusc2_seqs, file = "C:/Users/.../filtered_aa_tusc2_seqs.fasta")

# 6. Run your initial alignment (including the long N-terminal sequences)
alignment <- msa(filtered_aa_tusc2_seqs, method = "ClustalOmega")
aligned_matrix <- as.matrix(alignment)

# 7. Find the consensus/majority start column 
# Look for the first column where the majority of sequences have an 'M'
find_conserved_start <- function(matrix) {
  for (col in 1:ncol(matrix)) {
    column_data <- matrix[, col]
    m_count <- sum(column_data == "M")
    # If more than 50% of your 500+ sequences have an M here, this is the true start
    if (m_count > (nrow(matrix) * 0.50)) {
      return(col)
    }
  }
  return(1)
}

true_start_col <- find_conserved_start(aligned_matrix)

# 8. Trim everything before this column for ALL sequences in the alignment
# This chops off the unaligned N-terminal extensions cleanly
trimmed_matrix <- aligned_matrix[, true_start_col:ncol(aligned_matrix)]

# 9. Remove alignment gaps ("-") to recover the cleaned, unaligned sequences
# This allows you to re-align them properly without any bias
cleaned_list <- lapply(1:nrow(trimmed_matrix), function(i) {
  seq_vec <- trimmed_matrix[i, ]
  no_gaps <- seq_vec[seq_vec != "-"]
  paste(no_gaps, collapse = "")
})

# 10. Convert back to AAStringSet and name them correctly
cleaned_seqs <- AAStringSet(unlist(cleaned_list))
names(cleaned_seqs) <- rownames(aligned_matrix)

# 11. Run the final, gap-free Multiple Sequence Alignment
final_alignment <- msa(cleaned_seqs, method = "ClustalOmega")

# 12. Convert the alignment object to a standard AAStringSet
aligned_set <- as(final_alignment, "AAStringSet")

# 13. Save it as a FASTA file in your project directory
writeXStringSet(aligned_set, file = "C:/Users/.../Tusc2_final_alignment.fasta")

