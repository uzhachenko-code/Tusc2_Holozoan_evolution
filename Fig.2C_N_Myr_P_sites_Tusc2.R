# ==============================================================================
# FUS1/TUSC2 POSITION-SPECIFIC KINASE CENSUS & MONTE CARLO MASTER PIPELINE
# ==============================================================================
# Clear workspace memory and force immediate garbage collection to free RAM
rm(list = ls())
gc()

library(Biostrings)

# ------------------------------------------------------------------------------
# STEP 1: DEFINE SYSTEM DIRECTORIES AND FILE PATHS
# ------------------------------------------------------------------------------
base_dir   <- "C:/Users/..."
taxon_id_f <- paste0(base_dir, "Taxon_ID.txt")
fasta_in   <- paste0(base_dir, "datamonkey_fused_flat_protein.fasta")

# Final manuscript report destination path
master_pos_out <- paste0(base_dir, "position_specific_kinase_statistics_final.txt")

# ------------------------------------------------------------------------------
# STEP 2: LOAD SEQUENCES AND CLASSIFICATION METADATA
# ------------------------------------------------------------------------------
cat("[1/5] Loading primary protein dataset and Taxon_ID classification maps...\n")
prot_seqs <- readAAStringSet(fasta_in)
taxon_map <- read.table(taxon_id_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(taxon_map) <- c("Flat_ID", "Species", "Vertebrate")

n_initial <- length(prot_seqs)

# ------------------------------------------------------------------------------
# STEP 3: N-TERMINAL QUALITY FILTERING (PURGE FALSE ORFs & FRAGMENTS)
# ------------------------------------------------------------------------------
cat("[2/5] Filtering out truncated fragments and false open reading frames...\n")
valid_names <- c()
valid_seqs  <- c()
for (i in 1:n_initial) {
  # Clean out all alignment gap dashes to extract the true primary sequence length
  clean_str <- gsub("-", "", as.character(prot_seqs[i]))
  
  # CORE FILTER: Must start strictly with MG and maintain an un-gapped length of >= 20aa
  if (startsWith(clean_str, "MG") && nchar(clean_str) >= 20) {
    valid_names <- c(valid_names, names(prot_seqs)[i])
    valid_seqs  <- c(valid_seqs, clean_str)
  }
}
n_clean <- length(valid_seqs)
cat("-> Filtered functional sequences retained: ", n_clean, " out of ", n_initial, "\n")

# ------------------------------------------------------------------------------
# STEP 4: DEFINE KINASE REGEX PROFILES
# ------------------------------------------------------------------------------
kinase_regex_db <- c(
  ATM = "([ST]Q)", CKI = "([ST].[ST]|[ST]..[ST])", CKII = "([ST]..[DE])", CaM_II = "(R..[ST])",
  DNAPK = "([ST]Q)", EGFR = "([DE].*Y.*[LIVF]|[LIVF].*Y.*[DE])", GSK3 = "([ST]...[ST])",
  INSR = "([RK].{1,3}Y|Y.{1,3}[DE])", PKA = "(([RK]{2}.[ST])|([RK].[ST]))", PKB = "(R.R..[ST])",
  PKC = "(([RK].{0,2}[ST])|([ST].{0,2}[RK]))", PKG = "([RK]{2,3}.[ST])", RSK = "((R..R..[ST])|(RR.[ST]))",
  SRC = "([ED]GY.[FP])", cdc2 = "([ST]P.[RK])", cdk5 = "([ST]P)", p38MAPK = "([ST]P)"
)

# ------------------------------------------------------------------------------
# STEP 5: MAP BIOLOGICAL RESIDUE-SPECIFIC COHORT FREQUENCIES
# ------------------------------------------------------------------------------
cat("[3/5] Mapping exact residue-by-residue kinase profiles across positions 1 to 20...\n")
is_vert_vec  <- taxon_map$Vertebrate[match(valid_names, taxon_map$Flat_ID)]
n_vert_clean <- sum(is_vert_vec == 1)
n_inv_clean  <- sum(is_vert_vec == 0)

pos_stats_list <- list()

for (pos in 1:20) {
  cat(sprintf("\rProcessing position %d/20...", pos))
  
  for (k in seq_along(kinase_regex_db)) {
    kinase_name   <- names(kinase_regex_db)[k]
    regex_pattern <- kinase_regex_db[k]
    
    vert_hits <- 0
    inv_hits  <- 0
    
    for (i in 1:n_clean) {
      window_20 <- substr(valid_seqs[i], 1, 20)
      
      # Catch all motif matches within the focus window
      matches_list <- gregexpr(regex_pattern, window_20, perl = TRUE)
      matches      <- matches_list[[1]]
      
      if (matches[1] > 0) {
        # Cast lengths and positions explicitly as plain numerics to prevent type crashes
        match_lens <- as.numeric(attr(matches_list[[1]], "match.length"))
        matches    <- as.numeric(matches)
        
        is_hit_at_pos <- FALSE
        
        for (m_idx in seq_along(matches)) {
          m_start <- matches[m_idx]
          m_end   <- m_start + match_lens[m_idx] - 1
          
          if (pos >= m_start && pos <= m_end) {
            # Core validation rule: structural coordinate must match a phospho-target
            if (substr(window_20, pos, pos) %in% c("S", "T", "Y")) {
              is_hit_at_pos <- TRUE
              break
            }
          }
        }
        
        if (is_hit_at_pos) {
          if (is_vert_vec[i] == 1) vert_hits <- vert_hits + 1 else inv_hits <- inv_hits + 1
        }
      }
    }
    
    # Calculate percentages
    v_pct <- (vert_hits / n_vert_clean) * 100
    i_pct <- (inv_hits / n_inv_clean) * 100
    
    # Chi-Square statistical validation pass for taxonomic frequency shifts
    contingency_table <- matrix(c(vert_hits, n_vert_clean - vert_hits, inv_hits, n_inv_clean - inv_hits), nrow = 2)
    chisq_p <- tryCatch(chisq.test(contingency_table)$p.value, error = function(e) return(1.0))
    
    key_id <- paste0(pos, "_", kinase_name)
    pos_stats_list[[key_id]] <- data.frame(
      Position          = pos,
      Kinase            = kinase_name,
      Vertebrate_Hits   = vert_hits,
      Vertebrate_Pct    = round(v_pct, 2),
      Invertebrate_Hits = inv_hits,
      Invertebrate_Pct  = round(i_pct, 2),
      ChiSq_P_Value     = chisq_p,
      Observed_Global   = vert_hits + inv_hits,
      stringsAsFactors  = FALSE
    )
  }
}
cat("\n")
master_pos_df <- do.call(rbind, pos_stats_list)

# ------------------------------------------------------------------------------
# STEP 6: PURE VECTORIZED 10,000 PERMUTATION MONTE CARLO PROOF
# ------------------------------------------------------------------------------
cat("[4/5] Launching vectorized high-speed Monte Carlo permutation engine...\n")
set.seed(42)
n_permutations <- 10000

# Flush memory overhead to prepare raw RAM buffer for the big sampling pass
gc()

observed_totals <- master_pos_df$Observed_Global
perm_matrix     <- matrix(0, nrow = n_permutations, ncol = nrow(master_pos_df))
colnames(perm_matrix) <- rownames(master_pos_df)

k_names <- master_pos_df$Kinase
k_pos   <- master_pos_df$Position

# Unpack and un-gap the FULL sequences once to preserve exact background chemistry proportions
full_clean_sequences <- lapply(prot_seqs[valid_names], function(x) gsub("-", "", as.character(x)))

pb <- txtProgressBar(min = 1, max = n_permutations, style = 3)
for (p in 1:n_permutations) {
  
  # Scramble the letters natively using a single vector loop to keep memory compact
  scrambled_windows <- character(n_clean)
  for (i in 1:n_clean) {
    raw_chars <- strsplit(full_clean_sequences[[i]], "")[[1]]
    scrambled_windows[i] <- paste(sample(raw_chars, min(20, length(raw_chars)), replace = FALSE), collapse = "")
  }
  
  # Audit all position-kinase pairs using rapid vectorized row masking
  for (row_idx in 1:nrow(master_pos_df)) {
    pos       <- k_pos[row_idx]
    kinase_nm <- k_names[row_idx]
    reg_pat   <- kinase_regex_db[kinase_nm]
    
    valid_mask <- nchar(scrambled_windows) >= pos
    if (sum(valid_mask) > 0) {
      match_flags <- grepl(reg_pat, scrambled_windows[valid_mask], perl = TRUE) & 
        (substr(scrambled_windows[valid_mask], pos, pos) %in% c("S", "T", "Y"))
      
      perm_matrix[p, row_idx] <- sum(match_flags)
    }
  }
  
  # Update console screen ticks every 100 iterations to optimize thread performance
  if (p %% 100 == 0 || p == n_permutations) {
    setTxtProgressBar(pb, p)
  }
}
close(pb)

# ------------------------------------------------------------------------------
# STEP 7: DISK EXPORT COMPILATION PASS
# ------------------------------------------------------------------------------
cat("[5/5] Computing final empirical p-values and exporting dataset matrix...\n")

sim_means   <- colMeans(perm_matrix)
empirical_p <- numeric(nrow(master_pos_df))

for (row_idx in 1:nrow(master_pos_df)) {
  obs <- observed_totals[row_idx]
  empirical_p[row_idx] <- sum(perm_matrix[, row_idx] >= obs) / n_permutations
}

master_pos_df$Simulated_Mean_Random <- round(sim_means, 3)
master_pos_df$Empirical_P_Value     <- sprintf("%.4f", empirical_p)

# Filter out empty zero rows to keep the output text table clean and high-density
final_filtered_report <- master_pos_df[master_pos_df$Observed_Global > 0, ]

# Write final data sheet directly to your target drive directory folder
write.table(final_filtered_report, file = master_pos_out, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=========================================================\n")
cat("🎉 SUCCESS: STRUCTURAL SELECTION MATRIX SAVED DIRECTLY!\n")
cat("=========================================================\n")
cat("File written to: ", master_pos_out, "\n")
cat("Total active sequence-coordinate intersections logged: ", nrow(final_filtered_report), "\n")
cat("=========================================================\n")
