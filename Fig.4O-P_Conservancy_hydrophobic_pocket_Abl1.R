# =====================================================================
# SCRIPT 25: RIGOROUS HELIX POSITION CONSERVATION & ENTROPY ENGINE
# =====================================================================
cat("Wiping variable registers and initializing structural helix conservation calculator...\n")
rm(list = setdiff(ls(), "base_path"))
gc()

libs <- c("Biostrings", "openxlsx")
invisible(lapply(libs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

base_path       <- "C:/Users/..."
abl_fasta_path  <- file.path(base_path, "abl1_protein_aligned.fasta")

if (!file.exists(abl_fasta_path)) stop("Missing Error: 'abl1_protein_aligned.fasta' not found on specified path.")

# 1. INGEST ALIGNMENT MATRIX NATIVELY
abl_align  <- Biostrings::readAAStringSet(abl_fasta_path)
mat_abl    <- as.matrix(abl_align)
n_taxa     <- nrow(mat_abl)

# 2. ISOLATE HUMAN REFERENCE ROW AND GENERATE COORDINATE MAP
human_row_idx <- grep("homo_sapiens\\|ABL1", names(abl_align), ignore.case = TRUE)
if (length(human_row_idx) == 0 || is.na(human_row_idx)) {
  stop("Context Error: Could not locate the 'homo_sapiens|ABL1' reference sequence in your FASTA headers.")
}

human_seq_vector <- unlist(strsplit(as.character(abl_align[human_row_idx]), ""))

# Map gapless human residue positions directly to absolute alignment columns
human_to_alignment_col <- which(human_seq_vector != "-")

# Define target cavity positions explicitly grouped by their structural helix domains
helix_mapping <- data.frame(
  Helix = c(rep("\u03b1E", 3), rep("\u03b1F", 3), rep("\u03b1H", 4), rep("\u03b1I'", 3)),
  Human_Pos = c(356, 359, 360, 448, 452, 454, 481, 483, 484, 487, 521, 525, 529),
  Human_AA  = c("A", "L", "L", "L", "A", "Y", "D", "C", "P", "V", "I", "V", "L"),
  stringsAsFactors = FALSE
)

conservation_results <- list()
cat(paste("Successfully mapped human template matrix framework. Processing", n_taxa, "Holozoan sequences...\n"))

# 3. COMPUTE SHANNON ENTROPY AND PID OVER HUMAN-LOCKED TARGETS
for (i in 1:nrow(helix_mapping)) {
  helix_name <- helix_mapping$Helix[i]
  p          <- helix_mapping$Human_Pos[i]
  native_aa  <- helix_mapping$Human_AA[i]
  
  # Translate clean human residue coordinate directly into the alignment matrix column index
  align_col <- human_to_alignment_col[p]
  
  if (is.na(align_col) || align_col > ncol(mat_abl)) {
    cat(paste0("   [!] Warning: Human position ", native_aa, p, " falls outside alignment boundaries. Skipping.\n"))
    next
  }
  
  # Extract the isolated vertical column vector across all aligned organisms
  col_vector <- mat_abl[, align_col]
  
  # Clean out structural padding gaps to evaluate strictly real, substituted residue distributions
  col_vector <- col_vector[col_vector != "-"]
  
  if (length(col_vector) == 0) next
  
  # Compute absolute residue frequency counts
  counts <- table(col_vector)
  probs  <- counts / sum(counts)
  
  # CALCULATE SHANNON EVOLUTIONARY ENTROPY (H)
  shannon_entropy <- -sum(probs * log2(probs))
  
  # CALCULATE PERCENT SEQUENCE IDENTITY (PID) OF THE DOMINANT RESIDUE
  dominant_res  <- names(which.max(counts))
  percent_ident <- max(probs) * 100
  
  conservation_results[[i]] <- data.frame(
    Helix              = helix_name,
    Human_Residue      = paste0(native_aa, p),
    Alignment_Column   = align_col,
    Consensus_Residue  = dominant_res,
    Percent_Identity   = round(percent_ident, 2),
    Shannon_Entropy    = round(shannon_entropy, 4),
    stringsAsFactors   = FALSE
  )
}

final_conservation_df <- do.call(rbind, conservation_results)

# Export spreadsheet metrics directly to your drive archives
output_xlsx_path <- file.path(base_path, "Abl_cavities", "abl1_helix_specific_pocket_conservation.xlsx")
if (!dir.exists(dirname(output_xlsx_path))) dir.create(dirname(output_xlsx_path), recursive = TRUE)
openxlsx::write.xlsx(final_conservation_df, output_xlsx_path)

cat("\n=====================================================================\n=== HELIX POSITION CONSERVATION METRICS SECURED FOR REVIEWERS ===\n=====================================================================\n")
print(final_conservation_df, row.names = FALSE)
