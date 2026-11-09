# =====================================================================
# SCRIPT 1: DIAGNOSTIC FUSED MSA COHORT GENERATOR & EXPORTER
# =====================================================================
cat("Wiping variable registers and initializing diagnostic alignment exporter...\n")
rm(list = ls())
gc()

if (!requireNamespace("Biostrings", quietly = TRUE)) install.packages("Biostrings")
if (!requireNamespace("msa", quietly = TRUE)) install.packages("msa")
library(Biostrings)
library(msa)

base_path        <- "C:/Users/..."
output_subfolder <- file.path(base_path, "Allosteric_Core_Matrix")
setwd(base_path)

if (!dir.exists(output_subfolder)) dir.create(output_subfolder, recursive = TRUE)

tusc2_src_file  <- "Tusc2_protein.txt"
abl1_src_file   <- "abl1_protein.txt"

if (!file.exists(tusc2_src_file) || !file.exists(abl1_src_file)) stop("Error: Core source fasta files missing.")

all_abl_seqs    <- readAAStringSet(abl1_src_file)
tusc2_seqs      <- readAAStringSet(tusc2_src_file)
abl1_candidates <- all_abl_seqs[!grepl("abl2", names(all_abl_seqs), ignore.case = TRUE)]

# Clean metadata headers to match species names character-for-character
clean_header_species <- function(headers) {
  clean_sp <- character(length(headers))
  for (i in seq_along(headers)) {
    h <- headers[i]
    if (grepl('"organism_name":"[^"]+"', h)) {
      clean_sp[i] <- tolower(gsub(" ", "_", gsub('.*"organism_name":"([^"]+)".*', "\\1", h)))
    } else if (grepl("\\[[A-Za-z]+ [A-Za-z]+\\]", h)) {
      clean_sp[i] <- tolower(gsub(" ", "_", gsub("^.*\\[|\\].*$", "", h)))
    } else {
      clean_sp[i] <- tolower(gsub("^([^| ]+_[^| ]+).*$", "\\1", h))
    }
  }
  return(clean_sp)
}

names(abl1_candidates) <- clean_header_species(names(abl1_candidates))
names(tusc2_seqs)      <- clean_header_species(names(tusc2_seqs))
common_species_v       <- intersect(names(abl1_candidates), names(tusc2_seqs))

# Filter to shared species pairs
paired_abl1  <- abl1_candidates[names(abl1_candidates) %in% common_species_v]
paired_tusc2 <- tusc2_seqs[names(tusc2_seqs) %in% common_species_v]

# Run ClustalOmega Multiple Sequence Alignments
cat("\nExecuting ClustalOmega alignments over tracking cohorts...\n")
abl_align  <- as(msa(paired_abl1, method="ClustalOmega", type="protein", verbose=FALSE), "AAStringSet")
tusc_align <- as(msa(paired_tusc2, method="ClustalOmega", type="protein", verbose=FALSE), "AAStringSet")

# Isolate human reference row parameters to map structural columns natively
h_abl_idx  <- grep("homo_sapiens", names(abl_align), ignore.case = TRUE)[1]
h_tusc_idx <- grep("homo_sapiens", names(tusc_align), ignore.case = TRUE)[1]

human_abl_vec  <- unlist(strsplit(as.character(abl_align[h_abl_idx]), ""))
human_tusc_vec <- unlist(strsplit(as.character(tusc_align[h_tusc_idx]), ""))

# Identify absolute matrix column boundaries where human sequence tracks c-Abl resno 232-503
abl_start_col  <- which(cumsum(human_abl_vec != "-") == 232)
abl_end_col    <- which(cumsum(human_abl_vec != "-") == 503)

# Convert alignments directly to matrix structures to preserve structural alignment gaps (hyphens)
abl_matrix_raw  <- as.matrix(abl_align)
tusc_matrix_raw <- as.matrix(tusc_align)

# Crop matrix columns cleanly across the calculated human template tracking boundaries
abl_trimmed_mat  <- abl_matrix_raw[, abl_start_col:abl_end_col, drop = FALSE]
tusc_trimmed_mat <- tusc_matrix_raw  # Keep full Tusc2 alignment row length natively

# Re-extract clean species strings from rownames
abl_species_labels  <- tolower(gsub("^([^| ]+_[^| ]+).*$", "\\1", rownames(abl_trimmed_mat)))
tusc_species_labels <- tolower(gsub("^([^| ]+_[^| ]+).*$", "\\1", rownames(tusc_trimmed_mat)))
shared_taxa         <- intersect(abl_species_labels, tusc_species_labels)

fused_sequences <- character()
fused_headers   <- character()

# Merge matching organism tracks into a continuous, aligned multi-FASTA matrix file
for (sp in shared_taxa) {
  sp_abl_row  <- which(abl_species_labels == sp)[1]
  sp_tusc_row <- which(tusc_species_labels == sp)[1]
  
  # Concatenate the characters while keeping all alignment gap markers fully intact
  abl_str  <- paste(abl_trimmed_mat[sp_abl_row, ], collapse = "")
  tusc_str <- paste(tusc_trimmed_mat[sp_tusc_row, ], collapse = "")
  
  fused_seq_row <- paste0(abl_str, tusc_str)
  clean_header  <- paste0(toupper(sp), " | Aligned_Fused_Ortholog | c-Abl1(232-503)_Fused_Tusc2")
  
  fused_sequences <- c(fused_sequences, fused_seq_row)
  fused_headers   <- c(fused_headers, clean_header)
}

fused_output_set <- AAStringSet(fused_sequences)
names(fused_output_set) <- fused_headers

fasta_out_path <- file.path(output_subfolder, "Abl1_Kinase_Tusc2_Fused_Taxonomy_Aligned.fasta")
writeXStringSet(fused_output_set, filepath = fasta_out_path, format = "fasta")

cat("\n=====================================================================\n=== TRUE ALIGNED FUSED FASTA FILE EXPORTED SUCCESSFULLY ===\n=====================================================================\n")
cat(paste("Pristine intermediate aligned sequence file written directly to disk at:\n   ", fasta_out_path, "\n\n"))



# =====================================================================
# SCRIPT 2 FINAL ACCURACY EXTENSION: HUMAN-LOCK MASK, MI-APC & DOCKING JOIN
# =====================================================================
cat("Wiping environment variables and initializing human-lock validation suite...\n")
rm(list = setdiff(ls(), c("base_path", "output_subfolder", "docking_file", "meta_sheet_file")))
gc()

library(Biostrings)
library(openxlsx)
library(readxl)

# --- 1. RE-ESTABLISH ENVIRONMENT DEFAULTS ---
base_path        <- "C:/Users/..."
output_subfolder <- file.path(base_path, "Allosteric_Core_Matrix")
docking_file     <- file.path(base_path, "ClusPro results_kinase domain_cAbl_Tusc2", "all_cluspro_modes_interface_surfaces.xlsx")
meta_sheet_file  <- file.path(base_path, "metazoa_classification.xlsx")

fused_fasta_file <- file.path(output_subfolder, "Abl1_Kinase_Tusc2_Fused_Taxonomy_Aligned.fasta")
if (!file.exists(fused_fasta_file)) fused_fasta_file <- file.path(output_subfolder, "Abl1_Kinase_Tusc2_Fused_Ground_Truth.fasta")

if (!file.exists(fused_fasta_file)) stop("Error: Aligned intermediate fused multi-FASTA file not found on disk drive.")
if (!file.exists(docking_file)) stop("Error: Master docking workbook missing from pathways.")
if (!file.exists(meta_sheet_file)) stop("Error: Taxonomic classification workbook missing.")

# --- 2. INGEST FUSED TRACKS AND ISOLATE HOMO reference GRID ---
cat("[1/4] Loading intermediate aligned matrix arrays...\n")
fused_set <- readAAStringSet(fused_fasta_file)
meta_df   <- as.data.frame(readxl::read_excel(meta_sheet_file, sheet = 1))
meta_df$Species <- tolower(trimws(as.character(meta_df$Species)))

# Programmatically track taxonomic rows relative to your fasta headers cleanly
fasta_species_labels <- tolower(gsub("^([^| ]+_[^| ]+).*$", "\\1", names(fused_set)))

human_row_index <- grep("homo_sapiens", fasta_species_labels, ignore.case = TRUE)[1]
if (is.na(human_row_index)) stop("Error: Strict human reference row context missing from intermediate alignment.")

human_seq_vector <- unlist(strsplit(as.character(fused_set[human_row_index]), ""))
total_alignment_columns <- length(human_seq_vector)

# =====================================================================
# THE RECOVERED RULE: STRICT GAPLESS REFERENCE MAPPING
# =====================================================================
# Dynamically extract strictly the columns where human residues exist (drop other gap loops)
human_active_columns <- which(human_seq_vector != "-")
n_filtered_nodes     <- length(human_active_columns)

# Track clean coordinate residue numbers natively to prevent offset drift
human_abl_raw  <- unlist(strsplit(as.character(readAAStringSet(file.path(base_path, "abl1_protein.txt"))[1]), ""))
human_tusc_raw <- unlist(strsplit(as.character(readAAStringSet(file.path(base_path, "Tusc2_protein.txt"))[1]), ""))

abl_human_gaps  <- which(human_abl_raw != "-")
tusc_human_gaps <- which(human_tusc_raw != "-")

# Track lengths to cleanly separate c-Abl from Tusc2 nodes within the continuous sequence string
fused_matrix_raw <- as.matrix(fused_set)
len_abl_trimmed_cols <- which(cumsum(human_seq_vector != "-") == (503 - 232 + 1))[1]

# Convert alignment matrix directly to a numeric array based strictly on human active columns
cat("[2/4] Transflating alignment metrics to numeric probability matrices...\n")
pristine_align_mat <- fused_matrix_raw[, human_active_columns, drop = FALSE]
n_sequences <- nrow(pristine_align_mat)

aa_chars <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")
numeric_matrix <- matrix(0, nrow = n_sequences, ncol = n_filtered_nodes)
for (a in 1:21) { numeric_matrix[pristine_align_mat == aa_chars[a]] <- a }

# --- 3. SHANNON ENTROPY MUTUAL INFORMATION LOOP WITH APC ---
cat("[3/4] Running Mutual Information with Average Product Correction (MI-APC)... \n")
mi_matrix <- matrix(0, nrow = n_filtered_nodes, ncol = n_filtered_nodes)
pb_matrix <- txtProgressBar(min = 1, max = n_filtered_nodes - 1, style = 3)
for (i in 1:(n_filtered_nodes - 1)) {
  setTxtProgressBar(pb_matrix, i)
  for (j in (i + 1):n_filtered_nodes) {
    p_ij_counts <- table(factor(numeric_matrix[, i], levels=1:20), factor(numeric_matrix[, j], levels=1:20))
    N_real_pairs <- sum(p_ij_counts)
    if (N_real_pairs > 5) {
      p_ij <- p_ij_counts / N_real_pairs; p_i <- rowSums(p_ij); p_j <- colSums(p_ij); p_i_p_j <- p_i %o% p_j
      valid_bins <- p_ij > 0 & p_i_p_j > 0
      if (any(valid_bins)) {
        mi_matrix[i, j] <- sum(p_ij[valid_bins] * log2(p_ij[valid_bins] / p_i_p_j[valid_bins]))
        mi_matrix[j, i] <- mi_matrix[i, j]
      }
    }
  }
}
close(pb_matrix)

mean_rows <- rowMeans(mi_matrix); mean_all <- mean(mi_matrix)
apc_matrix <- matrix(0, nrow = n_filtered_nodes, ncol = n_filtered_nodes)
for (i in 1:n_filtered_nodes) {
  for (j in 1:n_filtered_nodes) {
    if (i != j && mi_matrix[i, j] > 0) apc_matrix[i, j] <- mi_matrix[i, j] - ((mean_rows[i] * mean_rows[j]) / mean_all)
  }
}

full_network_df <- as.data.frame(as.table(apc_matrix))
colnames(full_network_df) <- c("Filtered_Idx_I", "Filtered_Idx_J", "MI_APC_Score")
full_network_df$Filtered_Idx_I <- as.integer(full_network_df$Filtered_Idx_I)
full_network_df$Filtered_Idx_J <- as.integer(full_network_df$Filtered_Idx_J)

# Isolate cross-protein interactions cleanly
len_abl_boundary <- sum(human_seq_vector[1:len_abl_trimmed_cols] != "-")
inter_protein_rows <- full_network_df[full_network_df$Filtered_Idx_I <= len_abl_boundary & full_network_df$Filtered_Idx_J > len_abl_boundary, ]

inter_mean <- mean(inter_protein_rows$MI_APC_Score[inter_protein_rows$MI_APC_Score > 0], na.rm=TRUE)
inter_sd   <- sd(inter_protein_rows$MI_APC_Score[inter_protein_rows$MI_APC_Score > 0], na.rm=TRUE)
inter_protein_rows$Z_Score <- (inter_protein_rows$MI_APC_Score - inter_mean) / inter_sd

# Map residue numbers direct 1:1 relative to gapless human sequence vectors
inter_protein_rows$Abl1_PDB_Num  <- inter_protein_rows$Filtered_Idx_I + 231
inter_protein_rows$Tusc2_PDB_Num <- inter_protein_rows$Filtered_Idx_J - len_abl_boundary

# Enforce clean character labels matching standard human coordinates
inter_protein_rows$Human_Abl1_Residue  <- sapply(inter_protein_rows$Abl1_PDB_Num,  function(n) paste0(human_abl_raw[abl_human_gaps[n]], n))
inter_protein_rows$Human_Tusc2_Residue <- sapply(inter_protein_rows$Tusc2_PDB_Num, function(n) paste0(human_tusc_raw[tusc_human_gaps[n]], n))

ordered_df <- inter_protein_rows[inter_protein_rows$MI_APC_Score > 0, ]
ordered_df <- ordered_df[order(-ordered_df$Z_Score), ]

coevo_df <- data.frame(
  Human_Abl1_Residue  = ordered_df$Human_Abl1_Residue,
  Human_Tusc2_Residue = ordered_df$Human_Tusc2_Residue,
  Abl1_Col_Key        = ordered_df$Abl1_PDB_Num,
  Tusc2_Col_Key       = ordered_df$Tusc2_PDB_Num,
  Abl1_PDB_Num        = ordered_df$Abl1_PDB_Num,
  Tusc2_PDB_Num       = ordered_df$Tusc2_PDB_Num,
  MI_APC_Score        = ordered_df$MI_APC_Score,
  Z_Score             = ordered_df$Z_Score,
  stringsAsFactors    = FALSE
)
openxlsx::write.xlsx(coevo_df, file.path(output_subfolder, "fused_human_footprint_coevolution.xlsx"))

# --- 4. SECURE SPATIAL DOCKING NUMERIC CROSS-JOIN ---
cat("[4/4] Executing numeric cross-join against 3D ClusPro contact structures...\n")
docking_df <- as.data.frame(readxl::read_excel(docking_file))

clean_to_numeric_id <- function(column_vector) { as.integer(gsub("[^0-9]", "", as.character(column_vector))) }
docking_df$Abl1_PDB_Num  <- clean_to_numeric_id(docking_df$Abl1_Residue)
docking_df$Tusc2_PDB_Num <- clean_to_numeric_id(docking_df$Tusc2_Residue)
docking_df$Atom_Distance <- as.numeric(docking_df$Atom_Distance)

matched_hits <- merge(docking_df, coevo_df, by = c("Abl1_PDB_Num", "Tusc2_PDB_Num"))

if (nrow(matched_hits) > 0) {
  # Sort entries by ascending distance first so the shortest interaction model variant wins row 1
  matched_hits     <- matched_hits[order(matched_hits$Abl1_PDB_Num, matched_hits$Tusc2_PDB_Num, matched_hits$Atom_Distance), ]
  unique_pair_mask <- !duplicated(data.frame(matched_hits$Abl1_PDB_Num, matched_hits$Tusc2_PDB_Num))
  optimized_hits   <- matched_hits[unique_pair_mask, ]
  
  # Filter strictly using your original thresholds (Z >= 1.45 and distance <= 4.5 Å)
  final_manuscript_output <- optimized_hits[optimized_hits$Z_Score >= 1.45 & optimized_hits$Atom_Distance <= 4.50, ]
  final_manuscript_output$Tusc2_Domain <- ifelse(final_manuscript_output$Tusc2_PDB_Num >= 78 & final_manuscript_output$Tusc2_PDB_Num <= 98, "78-98_STEERING_ANCHOR", "1-80_FLEXIBLE_TETHER")
  
  final_output <- final_manuscript_output[order(-final_manuscript_output$Z_Score), c(
    "Scoring_Mode", "Model_File", "Abl1_Residue", "Tusc2_Residue", 
    "Atom_Distance", "Abl1_Col_Key", "Tusc2_Col_Key", "Abl1_PDB_Num", 
    "Tusc2_PDB_Num", "MI_APC_Score", "Z_Score", "Tusc2_Domain"
  )]
  
  openxlsx::write.xlsx(final_output, file.path(output_subfolder, "uninflated_docking_validated_by_coevolution.xlsx"))
  cat("\n=====================================================================\n=== MASTER MATRIX COEVOLUTION RECONSTRUCTION SUCCESSFUL ===\n=====================================================================\n")
  print(head(final_output, 20))
} else {
  cat("\nNotice: No direct numeric overlaps found. Outputting raw aligned coevolution scores standalone.\n")
}

# --- 5. EXTRACT CONSOLIDATED METAZOA CLASSIFICATION REGISTRY ---
cat("\nExtracting final taxonomic classification tracking matrix lists...\n")
species_used <- unique(fasta_species_labels)
taxa_registry_output <- meta_df[meta_df$Species %in% species_used, c("Species", "Group", "Vertebrate_status", "Class")]
taxa_registry_output <- taxa_registry_output[order(taxa_registry_output$Vertebrate_status, taxa_registry_output$Species), ]
openxlsx::write.xlsx(taxa_registry_output, file.path(output_subfolder, "taxa_registry_cohort_used.xlsx"))

# =====================================================================
# SCRIPT 5: COEVOLUTIONARY LANDSPACE HEATMAP GENERATOR
# =====================================================================
cat("Wiping memory buffers and initializing heatmap rendering engine...\n")
rm(list = ls())

libs <- c("readxl", "ggplot2", "openxlsx", "scales")
invisible(lapply(libs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# 1. ESTABLISH GROUND-TRUTH SUBFOLDER DIRECTORIES
base_path   <- "C:/Users/..."
target_file <- file.path(base_path, "MI_Z_scores_abl1_tusc2_n_idr_for_heatmap.xlsx")
out_img_dir <- file.path(base_path, "Allosteric_Core_Matrix")

if (!file.exists(target_file)) stop(paste("Error: Target Excel workbook missing at pathway:", target_file))
if (!dir.exists(out_img_dir)) dir.create(out_img_dir, recursive = TRUE)

# 2. INGEST EXCEL WORKBOOK NATIVELY 
cat("Ingesting datasheet matrix from your master local drive file...\n")
data_matrix <- as.data.frame(readxl::read_excel(target_file, sheet = 1))

# Enforce clean string and numeric conversions across column names
colnames(data_matrix)[1:2] <- c("Abl1_Axis", "Tusc2_Axis")
data_matrix$Z_Score <- as.numeric(data_matrix$Z_Score)

# 3. RE-ORDER AXES BY NATIVE STRUCTURAL RESIDUE NUMBERS
# Programmatic sort keys prevent R from listing residues alphabetically (e.g. keeps 454 before 481)
extract_numeric_id <- function(strings) { as.integer(gsub("[^0-9]", "", as.character(strings))) }
data_matrix$Abl1_Order  <- extract_numeric_id(data_matrix$Abl1_Axis)
data_matrix$Tusc2_Order <- extract_numeric_id(data_matrix$Tusc2_Axis)

data_matrix$Abl1_Axis  <- reorder(data_matrix$Abl1_Axis, data_matrix$Abl1_Order)
data_matrix$Tusc2_Axis <- reorder(data_matrix$Tusc2_Axis, data_matrix$Tusc2_Order)

# =====================================================================
# THE RECOVERED PLOT SHAPE: SCIENTIFIC REPORTS STANDARD GRAPHIC MATRIX
# =====================================================================
cat("Assembling tile geometries and painting color scaling matrices...\n")
heatmap_plot <- ggplot(data_matrix, aes(x = Abl1_Axis, y = Tusc2_Axis, fill = Z_Score)) +
  geom_tile(color = "white", lwd = 0.25, linetype = 1) +
  # Continuous high-contrast color scheme isolating your top outliers (GLU481 hotspots)
  scale_fill_gradientn(
    colors = c("#104E8B", "#87CEFA", "#FFFFFF", "#FF7F24", "#B22222"),
    values = rescale(c(min(data_matrix$Z_Score), 0, 1.45, 2.50, max(data_matrix$Z_Score))),
    name = "Coevolution\nZ-Score"
  ) +
  labs(
    title = "Cross-Protein Coevolutionary Landscape Tail Profile",
    subtitle = "Lineage-Balanced metazoan matrix unmasks allosteric docking hubs",
    x = "Human c-Abl1 Kinase Domain Residues",
    y = "Human Tusc2 N-Terminal IDR Residues"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b=5)),
    plot.subtitle = element_text(size = 10, italic = TRUE, hjust = 0.5, margin = margin(b=15)),
    axis.text.x   = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", color = "black"),
    axis.text.y   = element_text(face = "bold", color = "black"),
    axis.title    = element_text(face = "bold", size = 11),
    panel.grid    = element_blank(),
    legend.title  = element_text(face = "bold", size = 10),
    legend.position = "right"
  )

# 4. SAVE HIGH-RESOLUTION GRAPHIC TIFF FOR MANUSCRIPT SUBMISSION PACKAGE
out_tiff_path <- file.path(out_img_dir, "Fig_4F_coevolution_heatmap.tiff")
ggsave(
  filename = out_tiff_path, 
  plot = heatmap_plot, 
  device = "tiff", 
  dpi = 300, 
  width = 7.5, 
  height = 5.5, 
  units = "in", 
  compression = "lzw"
)

cat("\n=====================================================================\n=== HEATMAP PLOT SUCCESSFUL: TIFF SECURED AT 300 DPI ===\n=====================================================================\n")
cat(paste("Pristine, publication-grade image written directly to folder disk at:\n   ", out_tiff_path, "\n\n"))
print(heatmap_plot)


# =====================================================================
# SCRIPT 16: STRICT 15-28 WINDOW MOLECULAR PERMUTATION ENGINE (for A/G stretch)
# =====================================================================
cat("Wiping variable registers and initializing position-isolated 10k engine...\n")
rm(list = ls())
gc()

libs <- c("ggplot2", "reshape2", "openxlsx", "readxl", "scales", "Biostrings")
invisible(lapply(libs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# --- 1. DYNAMIC FILE PATH INGESTION VIA POPUPS ---
cat("\n[POPUP 1/2] Select your master 'metazoa_classification.xlsx' file...\n")
meta_file <- file.choose()

base_path        <- dirname(meta_file)
output_subfolder <- file.path(base_path, "Allosteric_Core_Matrix")
if (!dir.exists(output_subfolder)) dir.create(output_subfolder, recursive = TRUE)

cat("\nScanning for aligned/unaligned Tusc2 tracking sequences...\n")
tusc_fasta_check <- file.path(base_path, "Allosteric_Core_Matrix", "Tusc2_protein_aligned.fasta")
tusc_txt_check   <- file.path(base_path, "Tusc2_protein.txt")

if (file.exists(tusc_fasta_check)) {
  tusc_src_file <- tusc_fasta_check
} else if (file.exists(tusc_txt_check)) {
  tusc_src_file <- tusc_txt_check
} else {
  cat("\n[POPUP 2/2] Target sequence file not found automatically. Select manually...\n")
  tusc_src_file <- file.choose()
}

# 2. INGEST DATASETS NATIVELY
meta_df   <- as.data.frame(readxl::read_excel(meta_file, sheet = 1))
tusc_seqs <- Biostrings::readAAStringSet(tusc_src_file)

meta_df$Species           <- tolower(trimws(as.character(meta_df$Species)))
meta_df$Vertebrate_status <- trimws(as.character(meta_df$Vertebrate_status))
meta_df$Class             <- trimws(as.character(meta_df$Class))

# Map lineage factor levels exactly to your requested layout (Actinopterygii unified)
meta_df$Class_Group <- "Invertebrates"
is_vert_row <- !tolower(meta_df$Vertebrate_status) %in% c("invertebrates", "invertebrate")
meta_df$Class_Group[is_vert_row & grepl("mammalia|mammals", meta_df$Class, ignore.case=TRUE)] <- "Mammals"
meta_df$Class_Group[is_vert_row & grepl("aves", meta_df$Class, ignore.case=TRUE)]            <- "Aves"
meta_df$Class_Group[is_vert_row & grepl("reptilia|reptiles", meta_df$Class, ignore.case=TRUE)]  <- "Reptilia"
meta_df$Class_Group[is_vert_row & grepl("amphibia", meta_df$Class, ignore.case=TRUE)]          <- "Amphibia"
meta_df$Class_Group[is_vert_row & grepl("fish|actinopterygii|chondrichthyes", meta_df$Class, ignore.case=TRUE)] <- "Actinopterygii"

fasta_headers     <- names(tusc_seqs)
clean_fasta_names <- character(length(tusc_seqs))

for (i in seq_along(fasta_headers)) {
  h <- fasta_headers[i]
  if (grepl('"organism_name":"[^"]+"', h)) {
    clean_fasta_names[i] <- tolower(gsub(" ", "_", gsub('.*"organism_name":"([^"]+)".*', "\\1", h)))
  } else if (grepl("\\[[A-Za-z]+ [A-Za-z]+\\]", h)) {
    clean_fasta_names[i] <- tolower(gsub(" ", "_", gsub("^.*\\[|\\].*$", "", h)))
  } else if (grepl("\\|", h)) {
    pipe_parts <- unlist(strsplit(h, "\\|"))
    genus_sp_match <- grep("^[A-Za-z]+ [A-Za-z]+", trimws(pipe_parts), value = TRUE)
    if(length(genus_sp_match) > 0) {
      clean_fasta_names[i] <- tolower(gsub(" ", "_", genus_sp_match))
    } else {
      clean_fasta_names[i] <- tolower(gsub("^([^| ]+_[^| ]+).*$", "\\1", gsub(">", "", h)))
    }
  } else {
    clean_fasta_names[i] <- tolower(gsub("^([^| ]+_[^| ]+).*$", "\\1", gsub(">", "", h)))
  }
}
clean_fasta_names <- gsub(">", "", trimws(clean_fasta_names))

# --- 3. RUN DE-DUPLICATION AND 15-28 ISOLATED PERMUTATION LOOP ---
cat("\nRunning transcript de-duplication...\n")
unique_species_pool <- intersect(unique(clean_fasta_names), meta_df$Species)
results_list        <- list()

set.seed(42) # Lock seed for reproducible runs
cat("Executing 10,000 iterations isolated strictly within the 15-28 window boundaries...\n")
pb <- txtProgressBar(min = 0, max = length(unique_species_pool), style = 3)

for (idx in seq_along(unique_species_pool)) {
  setTxtProgressBar(pb, idx)
  sp <- unique_species_pool[idx]
  
  meta_row <- meta_df[meta_df$Species == sp, ]
  if (nrow(meta_row) == 0) next
  cl_group <- meta_row$Class_Group
  
  sp_fasta_row <- which(clean_fasta_names == sp) 
  raw_full_str <- as.character(tusc_seqs[sp_fasta_row])
  
  char_vec_full  <- unlist(strsplit(raw_full_str, ""))
  char_vec_clean <- char_vec_full[char_vec_full != "-"]
  
  # Ensure the sequence natively reaches at least position 28
  if (length(char_vec_clean) < 28) next 
  
  # ISOLATED WINDOW CAPTURE (Wipes out 2-42 padding entirely from background loops)
  window_15_28 <- char_vec_clean[15:28] 
  
  # LOCKED FRACTION MATH: Denominator is 14 columns, tracking strictly Alanine and Glycine
  obs_ag_fraction <- sum(window_15_28 %in% c("A", "G")) / 14
  
  # 10,000 PERMUTATION MONTE CARLO COUNTERS
  sim_shuffled_scores     <- numeric(10000)
  sim_alphabet_scores     <- numeric(10000)
  
  for (sim in 1:10000) {
    # Test 1: Position-Isolated Shuffle (Shuffles only the 14 residues present inside the window)
    sim_shuffled_scores[sim] <- sum(sample(window_15_28) %in% c("A", "G")) / 14
    
    # Test 2: Standard 20-letter background protein distribution sampling
    alphabet_sample <- sample(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"), 14, replace = TRUE)
    sim_alphabet_scores[sim] <- sum(alphabet_sample %in% c("A", "G")) / 14
  }
  
  p_val_shuffle <- sum(sim_shuffled_scores >= obs_ag_fraction) / 10000
  p_val_alpha   <- sum(sim_alphabet_scores >= obs_ag_fraction) / 10000
  
  results_list[[idx]] <- data.frame(
    Species         = sp,
    Class_Group     = cl_group,
    Obs_AG_Fraction = obs_ag_fraction,
    P_Value_Shuffle = p_val_shuffle,
    P_Value_Alpha   = p_val_alpha,
    stringsAsFactors = FALSE
  )
}
close(pb)

final_plot_df <- do.call(rbind, results_list)
final_plot_df <- final_plot_df[!is.na(final_plot_df$Class_Group), ]
final_plot_df$Class_Group <- factor(final_plot_df$Class_Group, 
                                    levels = c("Invertebrates", "Actinopterygii", "Amphibia", "Reptilia", "Aves", "Mammals"))

# --- 4. CALCULATE ELEVATED Y COORDINATES FOR THE VISUAL LABELS ---
summary_stats    <- aggregate(cbind(Obs_AG_Fraction, P_Value_Shuffle, P_Value_Alpha) ~ Class_Group, data = final_plot_df, FUN = mean)
max_scatter_dots <- aggregate(Obs_AG_Fraction ~ Class_Group, data = final_plot_df, FUN = max)

# Elevate Y anchors cleanly above the true highest empirical scatter dot to clear whitespace
summary_stats$Label_Y <- max_scatter_dots$Obs_AG_Fraction + 0.05

summary_stats$Significance_Label <- sapply(1:nrow(summary_stats), function(i) {
  lbl <- ""
  if (summary_stats$P_Value_Alpha[i] < 0.05)   lbl <- paste0(lbl, "*")
  if (summary_stats$P_Value_Shuffle[i] < 0.05) lbl <- paste0(lbl, " #")
  if (lbl == "") lbl <- "n.s."
  return(trimws(lbl))
})

# --- 5. RENDER CHRONOLOGICAL MATRIX PROPORTION BAR GRAPH ---
enrichment_plot <- ggplot(final_plot_df, aes(x = Class_Group, y = Obs_AG_Fraction)) +
  stat_summary(fun = mean, geom = "bar", fill = "#104E8B", color = "black", width = 0.55, alpha = 0.85, size = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.5, color = "black") +
  geom_jitter(color = "#FF7F24", alpha = 0.4, width = 0.15, size = 1.3) +
  geom_text(data = summary_stats, aes(x = Class_Group, y = Label_Y, label = Significance_Label), size = 4.2, fontface = "bold", color = "darkred") +
  labs(
    title = "Chronological Metazoan Scale of Localized Tusc2 A/G Fluid Hinge Enrichment",
    subtitle = "Proportions across residues 15-28 isolated from native flanking sequences (* p < 0.05 alphabet, # p < 0.05 shuffle)",
    x = "Taxonomic Lineage Grouping Classes",
    y = "Localized A/G Composition Proportion (Fraction)"
  ) +
  scale_y_continuous(labels = percent, limits = c(0, 1.05), expand = c(0, 0)) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b=4)),
    plot.subtitle = element_text(size = 9.5, italic = TRUE, hjust = 0.5, margin = margin(b=15)),
    axis.text     = element_text(face = "bold", color = "black", size = 10.5),
    axis.title    = element_text(face = "bold", size = 11),
    panel.grid    = element_blank()
  )

# SAVE HIGH-RESOLUTION WORKBOOKS AND TIFF IMAGES DIRECTLY TO MANUSCRIPT TARGET DIRECTORY
out_img_path <- file.path(output_subfolder, "Fig_4G_tusc2_isolated_15-28_proportions.tiff")
ggsave(filename = out_img_path, plot = enrichment_plot, device = "tiff", dpi = 300, width = 8.2, height = 5.4, units = "in", compression = "lzw")
openxlsx::write.xlsx(final_plot_df, file.path(output_subfolder, "tusc2_10k_isolated_15-28_ag_analysis.xlsx"))

cat("\n=====================================================================\n=== PIPELINE SUCCESSFUL: WINDOW RESTRICTED STRICTLY TO 15-28 CORIDOR ===\n=====================================================================\n")
print(summary_stats[, c("Class_Group", "Obs_AG_Fraction", "P_Value_Shuffle", "P_Value_Alpha", "Significance_Label")])
