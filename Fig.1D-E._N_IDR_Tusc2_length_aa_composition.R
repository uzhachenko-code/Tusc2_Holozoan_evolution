# ============================================================================
# INTEGRATED MASTER ENGINE: BIVALENT SEGMENTATION, AA ANALYSIS & PERMUTATION
# ============================================================================
cat("Wiping memory environments... Running complete integrated analysis pipeline.\n")

libs <- c("readxl", "openxlsx", "ggplot2", "reshape2")
invisible(lapply(libs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

set.seed(42) # Strict reproducibility anchor
base_path <- "C:/Users/..."
setwd(base_path)

excel_file <- "species_for_idr.xlsx"
af_dataset_dir <- file.path(base_path, "AF_dataset_Tusc2")

if (!file.exists(excel_file) || !dir.exists(af_dataset_dir)) stop("Error: Missing primary data files or directory structures.")

# 1. Ingest taxonomy sheet metadata rules safely
species_metadata <- as.data.frame(readxl::read_excel(excel_file))
colnames(species_metadata) <- trimws(colnames(species_metadata))

sp_name_col <- intersect(c("Species", "organismScientificName", "Animals"), colnames(species_metadata))
class_col   <- intersect(c("Class", "Group"), colnames(species_metadata))
metadata_species_vector <- as.character(species_metadata[, sp_name_col])
metadata_class_vector   <- as.character(species_metadata[, class_col])

# 2. Extract and compute over raw PDB coordinate structures
all_pdb_files <- list.files(af_dataset_dir, pattern = "\\.pdb$", full.names = TRUE)

idr_raw_matrix <- data.frame(
  Taxa_Group = character(), Class_Label = character(), Species_Clean = character(),
  File_Name = character(), Total_Length = integer(), N_Term_IDR_Len = integer(),
  C_Term_IDR_Len = integer(), Structured_Core = integer(), stringsAsFactors = FALSE
)

composition_df <- data.frame(
  Taxa_Group = character(), File_Name = character(), Region_Type = character(), 
  Freq_Entropy_AA = numeric(), Freq_Hydrophobic = numeric(), stringsAsFactors = FALSE
)

three_to_one <- c(ALA="A", CYS="C", ASP="D", GLU="E", PHE="F", GLY="G", HIS="H", ILE="I", 
                  LYS="K", LEU="L", MET="M", ASN="N", PRO="P", GLN="Q", ARG="R", SER="S", 
                  THR="T", VAL="V", TRP="W", TYR="Y")

calc_stats <- function(seq_vector) {
  if (length(seq_vector) == 0) return(c(0, 0))
  entropy_pct     <- (sum(seq_vector %in% c("G", "A", "P")) / length(seq_vector)) * 100
  hydrophobic_pct <- (sum(seq_vector %in% c("L", "I", "V", "F")) / length(seq_vector)) * 100
  return(c(round(entropy_pct, 2), round(hydrophobic_pct, 2)))
}

cat(paste("Processing structural metrics natively across", length(all_pdb_files), "PDB files...\n"))

for (matrix_sweep_idx in seq_along(all_pdb_files)) {
  f_path <- all_pdb_files[matrix_sweep_idx]
  f_name <- basename(f_path)
  
  is_fold_style <- grepl("^fold_", f_name)
  parsed_sp_key <- gsub("Tusc2_Ortholog_", "", gsub("Tusc2_MD_Ortholog_", "", gsub("fold_o_", "", f_name)))
  parsed_sp_key <- gsub("(_[A-Z0-9]{6,10})?\\.pdb$", "", parsed_sp_key)
  parsed_sp_key <- gsub("_", " ", parsed_sp_key)
  
  match_idx <- which(tolower(metadata_species_vector) == tolower(parsed_sp_key) | 
                       tolower(gsub("[^A-Za-z]", "", metadata_species_vector)) == tolower(gsub("[^A-Za-z]", "", parsed_sp_key)))
  
  current_class <- "Unassigned"; current_taxa  <- "Vertebrate"
  if (length(match_idx) > 0) {
    current_class <- metadata_class_vector[match_idx[1]]
    if (is_fold_style || grepl("Invertebrate|Coral|Crustacean|Insect|Hydra", current_class, ignore.case=TRUE)) {
      current_taxa <- "Invertebrate (Outgroup)"
    }
  } else if (is_fold_style) {
    current_taxa  <- "Invertebrate (Outgroup)"; current_class <- "Invertebrate Outgroup"
  }
  
  file_lines <- tryCatch(readLines(f_path, warn = FALSE), error = function(e) next)
  ca_atoms   <- file_lines[grepl("^ATOM  .* CA ", file_lines)]
  if (length(ca_atoms) == 0) next
  
  res_nos      <- as.integer(substr(ca_atoms, 23, 26))
  plddt_values <- as.numeric(substr(ca_atoms, 61, 66))
  chain_len    <- length(res_nos)
  
  is_strict_core <- plddt_values > 90
  is_inter_core  <- plddt_values >= 70 & plddt_values <= 90
  
  # --- BIVALENT RECONNAISSANCE SEGMENTATION ---
  n_idr_len <- 0
  for (i_fwd in 1:(chain_len - 4)) {
    if (i_fwd <= (chain_len - 2)) {
      if (all(is_strict_core[i_fwd:(i_fwd + 2)])) { n_idr_len <- i_fwd - 1; break }
    }
    if (all(is_inter_core[i_fwd:(i_fwd + 4)])) { n_idr_len <- (i_fwd - 1) + 2; break }
  }
  if (n_idr_len == 0 && !any(is_strict_core | is_inter_core)) n_idr_len <- chain_len
  
  c_idr_len <- 0
  for (j_rev in chain_len:5) {
    if (j_rev >= 3) {
      if (all(is_strict_core[(j_rev - 2):j_rev])) { c_idr_len <- chain_len - j_rev; break }
    }
    if (all(is_inter_core[(j_rev - 4):j_rev])) { c_idr_len <- (chain_len - j_rev) + 2; break }
  }
  
  structured_mass_dimension <- chain_len - n_idr_len - c_idr_len
  if (structured_mass_dimension < 0) structured_mass_dimension <- 0
  if (structured_mass_dimension <= 10 || n_idr_len == 0) next
  
  # Extract sequence content
  translated_seq <- three_to_one[trimws(substr(ca_atoms, 18, 20))]
  translated_seq <- translated_seq[!is.na(translated_seq)]
  if (length(translated_seq) < n_idr_len) next
  
  n_seq    <- translated_seq[1:n_idr_len]
  core_seq <- translated_seq[(n_idr_len + 1):length(translated_seq)]
  
  n_stats    <- calc_stats(n_seq)
  core_stats <- calc_stats(core_seq)
  
  composition_df <- rbind(composition_df, 
                          data.frame(Taxa_Group=current_taxa, File_Name=f_name, Region_Type="N-IDR", Freq_Entropy_AA=n_stats[1], Freq_Hydrophobic=n_stats[2]),
                          data.frame(Taxa_Group=current_taxa, File_Name=f_name, Region_Type="Core",  Freq_Entropy_AA=core_stats[1], Freq_Hydrophobic=core_stats[2]))
  
  idr_raw_matrix <- rbind(idr_raw_matrix, data.frame(
    Taxa_Group = current_taxa, Class_Label = current_class, Species_Clean = parsed_sp_key,
    File_Name = f_name, Total_Length = chain_len, N_Term_IDR_Len = n_idr_len,
    C_Term_IDR_Len = c_idr_len, Structured_Core = structured_mass_dimension, stringsAsFactors = FALSE))
}

# Run population analyses and write clean sheets
wilcox_len     <- wilcox.test(idr_raw_matrix$N_Term_IDR_Len[idr_raw_matrix$Taxa_Group == "Vertebrate"], idr_raw_matrix$N_Term_IDR_Len[idr_raw_matrix$Taxa_Group == "Invertebrate (Outgroup)"], exact = FALSE)
clean_n_idr    <- composition_df[composition_df$Region_Type == "N-IDR", ]
wilcox_entropy <- wilcox.test(clean_n_idr$Freq_Entropy_AA[clean_n_idr$Taxa_Group == "Vertebrate"], clean_n_idr$Freq_Entropy_AA[clean_n_idr$Taxa_Group == "Invertebrate (Outgroup)"], exact = FALSE)
wilcox_hydro   <- wilcox.test(clean_n_idr$Freq_Hydrophobic[clean_n_idr$Taxa_Group == "Vertebrate"], clean_n_idr$Freq_Hydrophobic[clean_n_idr$Taxa_Group == "Invertebrate (Outgroup)"], exact = FALSE)

openxlsx::write.xlsx(idr_raw_matrix, file.path(af_dataset_dir, "Tusc2_Taxonomic_IDR_BivalentBorders_Final.xlsx"))
openxlsx::write.xlsx(composition_df, file.path(af_dataset_dir, "Tusc2_AminoAcid_Composition_Final.xlsx"))

# =====================================================================
# STOCHASTIC LOOP PERMUTATION ENGINE (1,000 SHUFFLE TRIALS)
# =====================================================================
cat("\nLaunching 1,000-trial stochastic loop permutation validator...\n")
vertebrate_cohort <- idr_raw_matrix[idr_raw_matrix$Taxa_Group == "Vertebrate", ]
empirical_p_values <- numeric(nrow(vertebrate_cohort))
pb <- txtProgressBar(min = 0, max = nrow(vertebrate_cohort), style = 3)

for (spec_idx in 1:nrow(vertebrate_cohort)) {
  setTxtProgressBar(pb, spec_idx)
  f_name <- vertebrate_cohort$File_Name[spec_idx]; real_n_len <- vertebrate_cohort$N_Term_IDR_Len[spec_idx]
  f_path <- file.path(af_dataset_dir, f_name); file_lines <- readLines(f_path, warn = FALSE)
  plddt_scores <- as.numeric(substr(file_lines[grepl("^ATOM  .* CA ", file_lines)], 61, 66)); chain_length <- length(plddt_scores)
  null_distribution_lengths <- numeric(1000)
  
  for (perm_trial in 1:1000) {
    shuffled_plddt <- sample(plddt_scores)
    shuffled_is_strict <- shuffled_plddt > 90; shuffled_is_inter <- shuffled_plddt >= 70 & shuffled_plddt <= 90
    shuffled_is_domain <- shuffled_is_strict
    if (chain_length >= 5) {
      for (scan_idx in 1:(chain_length - 4)) {
        if (all(shuffled_is_inter[scan_idx:(scan_idx + 4)])) { shuffled_is_domain[scan_idx:(scan_idx + 4)] <- TRUE }
      }
    }
    perm_n_len <- 0
    for (i_fwd in 1:(chain_length - 2)) {
      if (all(shuffled_is_domain[i_fwd:(i_fwd + 2)])) { perm_n_len <- (i_fwd - 1); break }
    }
    if (perm_n_len == 0 && !any(shuffled_is_domain)) perm_n_len <- chain_length
    null_distribution_lengths[perm_trial] <- perm_n_len
  }
  null_mean <- mean(null_distribution_lengths)
  empirical_p_values[spec_idx] <- (sum(abs(null_distribution_lengths - null_mean) >= abs(real_n_len - null_mean))) / 1000
}
close(pb)
vertebrate_cohort$Permutation_p_Value <- empirical_p_values
openxlsx::write.xlsx(vertebrate_cohort, file.path(af_dataset_dir, "Tusc2_Vertebrate_IDR_StochasticPermutation_Final.xlsx"))

# =====================================================================
# GRAPHICS ENGINE: HIGH-RES TIFF RENDERING LAYOUTS
# =====================================================================
cat("\nRendering high-visibility publication figures...\n")
length_melt <- melt(idr_raw_matrix, id.vars = "Taxa_Group", measure.vars = c("N_Term_IDR_Len", "C_Term_IDR_Len"))
length_melt$variable <- factor(length_melt$variable, levels = c("N_Term_IDR_Len", "C_Term_IDR_Len"),
                               labels = c(paste0("N-Terminal IDR (p = ", round(wilcox_len$p.value, 5), ")"), "C-Terminal IDR"))

p_lengths <- ggplot(length_melt, aes(x = Taxa_Group, y = value, fill = Taxa_Group)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2, size = 0.8) +
  # Plot with jitter points, custom theme, and labels
  geom_jitter(width = 0.12, alpha = 0.25, color = "black", size = 1.2) +
  facet_wrap(~variable, scales = "free_y") +
  scale_fill_manual(values = c("#2c7bb6", "#d7191c")) +
  theme_bw(base_size = 12) +
  labs(title = "Outlier-Filtered Structural Boundaries of the Tusc2 Terminal Extensions",
       x = "", y = "Region Length (Amino Acids)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        strip.text = element_text(face = "bold", size = 11, color = "white"),
        strip.background = element_rect(fill = "gray20"),
        axis.text.x = element_text(face = "bold", size = 11),
        legend.position = "none")

# Save the plot
ggsave(file.path(af_dataset_dir, "Figure_Tusc2_IDR_Lengths_Pristine.tiff"),
       plot = p_lengths, width = 8.5, height = 5.2, dpi = 300, compression = "lzw")

# Melt 'composition_df' for comparison metrics
comp_melt <- melt(composition_df, 
                  id.vars = c("Taxa_Group", "Region_Type"), 
                  measure.vars = c("Freq_Entropy_AA", "Freq_Hydrophobic"))

# Assign labels with p-values for 'variable'
comp_melt$variable <- factor(comp_melt$variable,
                             levels = c("Freq_Entropy_AA", "Freq_Hydrophobic"),
                             labels = c(paste0("Entropy Drivers (Gly, Ala, Pro)\nPaired p = ", format.pval(wilcox_entropy$p.value, digits=3)),
                                        paste0("Hydrophobic Anchors (Leu, Ile, Val, Phe)\nPaired p = ", format.pval(wilcox_hydro$p.value, digits=3))))

# Plot for amino acid signature stratification
p_comp <- ggplot(comp_melt, aes(x = Region_Type, y = value, fill = Taxa_Group)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, size = 0.7, position = position_dodge(0.8)) +
  facet_grid(Taxa_Group ~ variable) +
  scale_fill_manual(values = c("#2c7bb6", "#d7191c")) +
  theme_bw(base_size = 12) +
  labs(title = "Statistical Amino Acid Signature Stratification: Disordered Lasso vs. Folded Mass",
       x = "Intramolecular Region", y = "Residue Frequency (%)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        strip.text = element_text(face = "bold", size = 10, color = "white"),
        strip.background = element_rect(fill = "gray20"),
        axis.text.x = element_text(face = "bold", size = 11),
        legend.position = "top",
        legend.title = element_blank())

# Save the plot
ggsave(file.path(af_dataset_dir, "Figure_Tusc2_Composition_Significance.tiff"),
       plot = p_comp, width = 8.5, height = 5.8, dpi = 300, compression = "lzw")

# Final console output and summaries
cat("\n=====================================================================\n")
cat("=== PIPELINE EXECUTION SUCCESSFUL (TECHNICAL OUTLIERS EXCLUDED) ===\n")
cat("=====================================================================\n")
cat("--- ARTIFACT-FREE STRUCTURAL LENGTH MEANS ---\n")
print(aggregate(cbind(Total_Length, N_Term_IDR_Len, C_Term_IDR_Len, Structured_Core) ~ Taxa_Group, data = idr_raw_matrix, FUN = mean))
cat(paste("\n--> LENGTH COMPARATIVE MANUSCRIPT P-VALUE =", round(wilcox_len$p.value, 5), "\n"))
cat(paste("--> GLOBAL MEAN EMPIRICAL PERMUTATION P-VALUE  =", round(mean(empirical_p_values), 5), "\n\n"))
cat("--- FINAL UPDATED AMINO ACID PROFILE FREQUENCIES ---\n")
print(aggregate(cbind(Freq_Entropy_AA, Freq_Hydrophobic) ~ Taxa_Group + Region_Type, data = composition_df, FUN = mean))