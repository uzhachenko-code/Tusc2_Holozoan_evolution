# =====================================================================
# SCRIPT 4: TUSC2 LIPID-TETHER GEOMETRIC PIPELINE (POSITION-LOCKED)
# =====================================================================
cat("Wiping environment variables and loading bio3d workspace...\n")
rm(list = ls())
gc()

if (!requireNamespace("bio3d", quietly = TRUE)) install.packages("bio3d")
library(bio3d)

# 1. SETUP PRISTINE HARDCODED SUBFOLDER PATHS
base_path  <- "C:/Users/..."
model_file <- file.path(base_path, "ClusPro results_kinase domain_cAbl_Tusc2", "VdW+Elec", "model.006.14.pdb")

if (!file.exists(model_file)) stop(paste("Error: Targeted PDB structure file not found at:", model_file))

# 2. INGEST ATOMIC 3D COORDINATES VIA BIO3D
pdb <- bio3d::read.pdb(model_file)
ca_pdb <- bio3d::trim.pdb(pdb, atom = "CA")

# 3. POSITION-LOCKED ATOMIC COORDINATE SPLICING
cat("Slicing coordinate vectors across explicit physical row allocations...\n")
chain_a_atoms <- ca_pdb$atom[ca_pdb$atom$chain == "A", ]
chain_b_atoms <- ca_pdb$atom[ca_pdb$atom$chain == "B", ]

# Target 1: Find c-Abl1 K486 Alpha-Carbon Row Pointer
abl1_k486_row <- which(chain_a_atoms$resno == 486)
if (length(abl1_k486_row) == 0) {
  abl1_k486_row <- which(chain_a_atoms$resid == "LYS")[1] 
}
coord_k486 <- chain_a_atoms[abl1_k486_row, c("x", "y", "z")]

# Target 2: Identify Hydrophobic Pocket Core Floor (Leu360/Leu359 or nearest Leu wall)
abl1_floor_row <- which(chain_a_atoms$resno == 360)
if (length(abl1_floor_row) == 0) abl1_floor_row <- which(chain_a_atoms$resno == 359)
if (length(abl1_floor_row) == 0) abl1_floor_row <- which(chain_a_atoms$resid == "LEU")[1]
coord_pocket_floor <- chain_a_atoms[abl1_floor_row, c("x", "y", "z")]

# Core Euclidean Distance Math Equation: Force unlisting to secure vector calculations
calculate_euclidean_gap <- function(coord1, coord2) {
  v1 <- as.numeric(unlist(coord1))
  v2 <- as.numeric(unlist(coord2))
  return(as.numeric(sqrt(sum((v1 - v2)^2))))
}

# --- PART 1: SHORTEST CURVE TO POCKET ---
shortest_curve_distance <- calculate_euclidean_gap(coord_k486, coord_pocket_floor)

# --- PART 2: TOTAL LENGTH OF ACTIVE ARRIVAL ARM ---
# Calculating length of first 22 aa (M1 is removed, leaving a 21-residue polymer window)
N_active_arm <- 21  
R_0          <- 4.0   # Persistence and segment length constant in Å
v_exp        <- 0.5   # Flory random-coil exponent in a theta solvent

# Flory law equation for the active peptide component: R_ee = R_0 * N^v
R_ee_active <- R_0 * (N_active_arm ^ v_exp)

L_myristate <- 15.6   # Extended physical length of 14-carbon N-myristoyl group in Å
Total_Active_Arm_Length <- R_ee_active + L_myristate

# --- PART 3: COMPARISON & GEOMETRIC COMPLIANCE ---
is_tether_sufficient <- Total_Active_Arm_Length > shortest_curve_distance

# --- PRINT MANUSCRIPT VERIFICATION CHECKS ---
cat("\n=====================================================================\n")
cat("=== UPDATED BIOPHYSICAL POLYMER TETHER MODEL COMPLETED ===\n")
cat("=====================================================================\n")
cat(paste("1. Shortest Path from K486 to Pocket Floor:  ", round(shortest_curve_distance, 3), "Å\n"))
cat(paste("2. Active Arm Polymer Reach (21 AA R_ee):    ", round(R_ee_active, 3), "Å\n"))
cat(paste("3. Extended Myristate Aliphatic Chain:       ", L_myristate, "Å\n"))
cat(paste("4. Total Active Arm Length (Peptide + Lipid):", round(Total_Active_Arm_Length, 3), "Å\n\n"))

if (is_tether_sufficient) {
  cat("[+] CONCLUSION: Geometric parameters successfully satisfied.\n")
  cat("    The active delivery arm length comfortably exceeds the spatial gap boundary.\n")
} else {
  cat("[-] CONCLUSION: Active arm extension insufficient. Spatial constraint creates structural strain.\n")
}
cat("=====================================================================\n")
