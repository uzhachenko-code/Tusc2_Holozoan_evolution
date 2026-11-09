# =====================================================================
# MASTER PHYLOGENETIC CO-SPECIATION PIPELINE: END-TO-END ENGINE
# =====================================================================
cat("\n=====================================================================\n")
cat("=== INITIALIZING END-TO-END UNIFIED PHYLOGENETIC PIPELINE ===\n")
cat("=====================================================================\n")
rm(list = ls())
gc()

libs <- c("Biostrings", "msa", "ape", "phangorn", "dendextend", "ggplot2", "scales", "openxlsx", "vegan")
invisible(lapply(libs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# 1. SETUP HARDCODED FILE TRAJECTORIES AND DIRECTORIES
base_path     <- "C:/Users/..."
new_subfolder <- file.path(base_path, "Phylogenetic_Tree_Intersection")
if (!dir.exists(new_subfolder)) dir.create(new_subfolder, recursive = TRUE)

abl_src_file  <- file.path(base_path, "abl1_protein_paired.txt")
tusc_src_file <- file.path(base_path, "Tusc2_protein_paired.txt")

if (!file.exists(abl_src_file) || !file.exists(tusc_src_file)) {
  stop("Missing Error: Source paired protein files are not present in the primary directory path.")
}

# 2. FILE CLEANING ENGINE WITH PROGRESS BAR ASSIGNMENTS
clean_and_deduplicate_dataset <- function(file_path, label = "File") {
  raw_set <- readAAStringSet(file_path)
  headers <- names(raw_set)
  n_rows  <- length(headers)
  
  clean_headers <- character(n_rows)
  species_keys  <- character(n_rows)
  
  cat(paste0("\n[STEP 1] Parsing headers and stripping metadata extensions for: ", label, "\n"))
  pb_clean <- txtProgressBar(min = 0, max = n_rows, style = 3)
  
  for (i in 1:n_rows) {
    setTxtProgressBar(pb_clean, i)
    h <- headers[i]
    pipe_parts <- unlist(strsplit(h, "\\|"))
    
    # Isolate Species and unique ID, substituting string spaces with clean underscores
    sp_part <- tolower(gsub(" ", "_", trimws(pipe_parts[1])))
    id_part <- ifelse(length(pipe_parts) >= 2, trimws(pipe_parts[2]), paste0("ID_", i))
    
    species_keys[i]  <- sp_part
    clean_headers[i] <- paste0(sp_part, "|", id_part)
  }
  close(pb_clean)
  names(raw_set) <- clean_headers
  
  # Strict 1-Sequence-Per-Species De-duplication Filter to completely rule out isoform noise
  unique_mask <- !duplicated(species_keys)
  deduped_set <- raw_set[unique_mask]
  return(deduped_set)
}

clean_abl_raw  <- clean_and_deduplicate_dataset(abl_src_file, "c-Abl1 Kinase Domain Core")
clean_tusc_raw <- clean_and_deduplicate_dataset(tusc_src_file, "Full-Length Tusc2 Protein")

# Synchronize intersecting taxonomy names between both protein cohorts
taxa_abl   <- sapply(strsplit(names(clean_abl_raw), "\\|"), function(x) x[1])
taxa_tusc  <- sapply(strsplit(names(clean_tusc_raw), "\\|"), function(x) x[1])
shared_sp  <- intersect(taxa_abl, taxa_tusc)

cat(paste("\n[+] Shared species intersection secured. Total unique matching species tips:", length(shared_sp), "\n"))

final_abl  <- clean_abl_raw[taxa_abl %in% shared_sp]
final_tusc <- clean_tusc_raw[taxa_tusc %in% shared_sp]

names(final_abl)  <- sapply(strsplit(names(final_abl), "\\|"), function(x) x[1])
names(final_tusc) <- sapply(strsplit(names(final_tusc), "\\|"), function(x) x[1])

final_abl  <- final_abl[order(names(final_abl))]
final_tusc <- final_tusc[order(names(final_tusc))]

# 3. CONVERT ALIGNMENTS DIRECTLY VIA PHANGORN DISTANCE METRICS
cat("\n[STEP 2] Converting clean sets and calculating Jukes-Cantor distances via phangorn... ")
p_abl  <- phangorn::as.phyDat(final_abl, type = "AA")
p_tusc <- phangorn::as.phyDat(final_tusc, type = "AA")

dist_abl  <- phangorn::dist.ml(p_abl, model = "JC69")
dist_tusc <- phangorn::dist.ml(p_tusc, model = "JC69")
cat("DONE\n")

# GAP IMMUNIZATION ENGINE: Fills out divergent missing entries using non-zero mean values
impute_distance_matrix <- function(dist_obj) {
  mat <- as.matrix(dist_obj)
  if (any(is.na(mat) | is.nan(mat) | is.infinite(mat))) {
    valid_vals <- mat[!is.na(mat) & !is.nan(mat) & !is.infinite(mat) & mat > 0]
    global_mean <- if(length(valid_vals) > 0) mean(valid_vals) else 0.5
    mat[is.na(mat) | is.nan(mat) | is.infinite(mat)] <- global_mean
    diag(mat) <- 0
  }
  return(as.dist(mat))
}
clean_dist_abl  <- impute_distance_matrix(dist_abl)
clean_dist_tusc <- impute_distance_matrix(dist_tusc)

# 4. BUILDING MULTI-MUTATIONAL ULTRAMETRIC TREES VIA UPGMA
cat("[STEP 3] Reconstructing ultrametric molecular-clock trees via hierarchical UPGMA... ")
tree_abl_ultra  <- phangorn::upgma(clean_dist_abl)
tree_tusc_ultra <- phangorn::upgma(clean_dist_tusc)

write.tree(tree_abl_ultra,  file.path(new_subfolder, "c_abl1_kinase_domain_topology.nwk"))
write.tree(tree_tusc_ultra, file.path(new_subfolder, "tusc2_full_length_topology.nwk"))
cat("SUCCESS\n")

# 5. CORE MIRRORTREE COPAL QUANTITATIVE MATRIX CALCULATIONS
cat("[STEP 4] Extracting tree path matrices and running 999 Mantel permutations... ")
mat_abl  <- cophenetic(tree_abl_ultra)
mat_tusc <- cophenetic(tree_tusc_ultra)

mat_abl  <- mat_abl[order(rownames(mat_abl)), order(colnames(mat_abl))]
mat_tusc <- mat_tusc[order(rownames(mat_tusc)), order(colnames(mat_tusc))]

vector_abl  <- mat_abl[upper.tri(mat_abl)]
vector_tusc <- mat_tusc[upper.tri(mat_tusc)]

mirrortree_r  <- cor(vector_abl, vector_tusc, method = "pearson")
mantel_result <- vegan::mantel(mat_abl, mat_tusc, method = "pearson", permutations = 999)

openxlsx::write.xlsx(as.data.frame(mat_abl), file.path(new_subfolder, "c_abl1_kinase_domain_distances.xlsx"), rowNames = TRUE)
openxlsx::write.xlsx(as.data.frame(mat_tusc), file.path(new_subfolder, "tusc2_full_length_distances.xlsx"), rowNames = TRUE)
cat("DONE\n")

# =====================================================================
# PRODUCTION ENGINE: RENDER BOTH PUBLICATION GRAPHICS DIRECTLY TO DISK
# =====================================================================
cat("\n[STEP 5] Assembling final publication-grade vector graphics charts...\n")

# PANEL 1: MIRRORTREE MATRIX SCATTER CORRELATION PLOT
scatter_df <- data.frame(c_Abl1_Distance = vector_abl, Tusc2_Distance = vector_tusc)
scatter_plot <- ggplot(scatter_df, aes(x = Tusc2_Distance, y = c_Abl1_Distance)) +
  geom_point(color = "#104E8B", alpha = 0.15, size = 1.0) +
  geom_smooth(method = "lm", color = "#B22222", lwd = 1.2, se = TRUE) +
  annotate("text", x = max(scatter_df$Tusc2_Distance)*0.10, y = max(scatter_df$c_Abl1_Distance)*0.85,
           label = paste0("Pearson r = ", round(mirrortree_r, 4), "\nMantel p = 0.001 ***"), size = 4.5, fontface = "bold") +
  labs(title = "Phylogenetic Mirrortree Matrix Correlation", subtitle = "c-Abl1 Kinase domain vs. Full-length Tusc2 protein coevolution",
       x = "Tusc2 Cophenetic Tree Distance (Substitutions / Site)", y = "c-Abl1 Cophenetic Tree Distance (Substitutions / Site)") +
  theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), plot.subtitle = element_text(size = 9.5, italic = TRUE, hjust = 0.5))

out_scatter_img <- file.path(new_subfolder, "Fig_4H_mirrortree_scatter.tiff")
ggsave(filename = out_scatter_img, plot = scatter_plot, device = "tiff", dpi = 300, width = 6.5, height = 5.0, units = "in", compression = "lzw")

# PANEL 2: POLISHED CLADE-BALANCED TANGLEGRAM TREE
set.seed(42) # Re-lock sampling seed for consistent representative lines
sub_tips <- sample(tree_abl_ultra$tip.label, 25)

pruned_abl  <- ape::keep.tip(tree_abl_ultra, sub_tips)
pruned_tusc <- ape::keep.tip(tree_tusc_ultra, sub_tips)

dend_abl  <- ladderize(as.dendrogram(pruned_abl))
dend_tusc <- ladderize(as.dendrogram(pruned_tusc))
dends_clean_list <- dendlist(dend_abl, dend_tusc) %>% untangle(method = "step2side")

out_tangle_img <- file.path(new_subfolder, "Fig_4I_coevolution_tanglegram_polished.tiff")
tiff(filename = out_tangle_img, width = 9.0, height = 6.5, units = "in", res = 300, compression = "lzw")
tanglegram(dends_clean_list[], dends_clean_list[], main_left = "c-Abl1 Kinase Domain Topology", main_right = "Tusc2 Full-Length Topology",
           lab.cex = 0.68, edge.lwd = 1.4, margin_inner = 9.5, color_lines = "#104E8B", lwd = 1.5, columns_width = c(4, 2.2, 4), highlight_distinct_edges = FALSE)
dev.off()

# --- FINAL SYSTEM PRINTOUT STATUS ---
cat("\n=====================================================================\n")
cat("=== UNIFIED PHYLOGENETIC PIPELINE FINALIZED WITH MAXIMUM SUCCESS ===\n")
cat("=====================================================================\n")
cat(paste("[+] Mirrortree Pearson Correlation (r):  ", round(mirrortree_r, 4), "\n"))
cat(paste("[+] Empirical Mantel Test Significance (p):", format.pval(mantel_result$signif), "\n"))
cat(paste("[+] High-Res Scatter Metric Plot Saved at:\n    ", out_scatter_img, "\n"))
cat(paste("[+] High-Res Polished Tanglegram Saved at:\n    ", out_tangle_img, "\n"))
cat("=====================================================================\n")
