#------------------------------
#EXTRACTION OF 3P21.3 ORTHOLOGS
#------------------------------

library(httr)
library(jsonlite)
library(dplyr)

# Standard Unified Ensembl API Gateway
BASE_URL <- "https://rest.ensembl.org"

# Target 8 core tumor suppressor genes
target_genes <- c("RASSF1", "TUSC2", "ZMYND10", "NPRL2", "HYAL2", "CACNA2D2", "CYB561D2", "TMEM115")
taxon_id <- "33208" # NCBI Taxonomy ID for Metazoa (All animals)

# Query both database divisions consecutively on the same server
compara_divisions <- c("vertebrates", "metazoa")

all_orthologs_list <- list()

print("Step 1: Gathering all unique ortholog IDs across both database scopes...")

for (gene_symbol in target_genes) {
  message(paste0("\nProcessing: ", gene_symbol, "..."))
  
  for (division in compara_divisions) {
    homology_ext <- paste0("/homology/symbol/human/", gene_symbol)
    
    # Configure the query matrix parameter safely
    query_params <- list(
      type = "orthologues",
      target_taxon = taxon_id
    )
    if (division == "metazoa") {
      query_params$compara <- "metazoa"
    }
    
    hom_resp <- GET(
      url = paste0(BASE_URL, homology_ext),
      add_headers("Accept" = "application/json"),
      query = query_params
    )
    
    if (status_code(hom_resp) == 200) {
      hom_data <- fromJSON(content(hom_resp, "text", encoding = "UTF-8"), simplifyVector = FALSE)
      
      # FIX: Dynamically isolate the 'homologies' block depending on division layout
      homologies_list <- NULL
      if (division == "vertebrates") {
        if (!is.null(hom_data[["data"]]) && length(hom_data[["data"]]) > 0) {
          homologies_list <- hom_data[["data"]][[1]][["homologies"]]
        }
      } else if (division == "metazoa") {
        if (!is.null(hom_data[["data"]]) && !is.null(hom_data[["data"]][["homologies"]])) {
          homologies_list <- hom_data[["data"]][["homologies"]]
        }
      }
      
      if (!is.null(homologies_list) && length(homologies_list) > 0) {
        # Safely unpack list cell-by-cell to guarantee matching dimensions
        parsed_div <- lapply(homologies_list, function(item) {
          target <- item[["target"]]
          if (!is.null(target)) {
            data.frame(
              human_gene_name = gene_symbol,
              target_species  = ifelse(is.null(target[["species"]]), NA, target[["species"]]),
              ortholog_id     = ifelse(is.null(target[["id"]]), NA, target[["id"]]),
              division_source = division,
              stringsAsFactors = FALSE
            )
          }
        })
        
        div_df <- bind_rows(parsed_div)
        all_orthologs_list[[length(all_orthologs_list) + 1]] <- div_df
        message(paste0("   -> Found ", nrow(div_df), " rows in '", division, "' database scope."))
      } else {
        message(paste0("   -> No structural data found for '", division, "'."))
      }
    } else {
      warning(paste("   -> Failed API request for:", division, "Code:", status_code(hom_resp)))
    }
    Sys.sleep(0.15) # Polite rate limiting compliance
  }
}

# Bind and deduplicate results globally
if (length(all_orthologs_list) > 0) {
  orthologs_master <- bind_rows(all_orthologs_list) %>% 
    distinct(human_gene_name, target_species, ortholog_id, .keep_all = TRUE)
  
  message(paste("\nTotal unique ortholog IDs accumulated globally:", nrow(orthologs_master)))
} else {
  stop("No records retrieved. Confirm that rest.ensembl.org can be accessed from your machine.")
}

print("\nStep 2: Resolving precise genomic coordinates via batch lookups...")

unique_ids <- unique(orthologs_master$ortholog_id)
unique_ids <- unique_ids[!is.na(unique_ids)]

chunk_size <- 500
coord_mappings <- list()

for (i in seq(1, length(unique_ids), by = chunk_size)) {
  chunk_ids <- unique_ids[i:min(i + chunk_size - 1, length(unique_ids))]
  chunk_ids <- chunk_ids[!is.na(chunk_ids)]
  
  message(paste("Processing batch lookup for IDs:", i, "to", min(i + chunk_size - 1, length(unique_ids))))
  
  req_body <- toJSON(list(ids = chunk_ids))
  
  lookup_resp <- POST(
    url = paste0(BASE_URL, "/lookup/id"),
    add_headers("Accept" = "application/json", "Content-Type" = "application/json"),
    body = req_body
  )
  
  if (status_code(lookup_resp) == 200) {
    lookup_data <- fromJSON(content(lookup_resp, "text", encoding = "UTF-8"), simplifyVector = FALSE)
    
    parsed_coords <- lapply(names(lookup_data), function(id) {
      gene_info <- lookup_data[[id]]
      if (!is.null(gene_info)) {
        data.frame(
          ortholog_id     = id,
          ortholog_chrom  = ifelse(is.null(gene_info$seq_region_name), NA, gene_info$seq_region_name),
          ortholog_start  = ifelse(is.null(gene_info$start), NA, gene_info$start),
          ortholog_end    = ifelse(is.null(gene_info$end), NA, gene_info$end),
          ortholog_strand = ifelse(is.null(gene_info$strand), NA, gene_info$strand),
          stringsAsFactors = FALSE
        )
      }
    })
    coord_mappings[[length(coord_mappings) + 1]] <- bind_rows(parsed_coords)
  }
  Sys.sleep(0.2)
}

# Merge physical coordinates with master file
coords_df <- bind_rows(coord_mappings)
final_dataset <- left_join(orthologs_master, coords_df, by = "ortholog_id")

# Save complete file to active directory
write.csv(final_dataset, "all_metazoan_3p21_3_coordinates.csv", row.names = FALSE)

print("--- Processing Success! ---")
print(paste("Final output contains", nrow(final_dataset), "rows mapping both server ecosystems."))
print(head(final_dataset[, c("human_gene_name", "target_species", "ortholog_chrom", "ortholog_start")]))



#------------------------------
#CALCULATING 3P21.3 CLUSTERS
#------------------------------
library(dplyr)
library(tidyr)
library(readxl)   # Required to read your local Excel file

# --- CONFIGURATION PATHWAYS ---
EXCEL_PATH <- "C:/Users/.../Species.xlsx"

# Ensure your previous final_dataset is active in memory. If loading from a file:
# final_dataset <- read.csv("all_metazoan_3p21_3_coordinates.csv", stringsAsFactors = FALSE)

# ==============================================================================
# LEVEL 1: LOAD METADATA DICTIONARY & DEFINE LEVELS
# ==============================================================================
print("Level 1: Loading species metadata annotations from Excel...")
if (!file.exists(EXCEL_PATH)) {
  stop(paste("Excel file not found at path:", EXCEL_PATH))
}

# Load the dictionary and map the exact Excel header names to our working columns
species_metadata <- read_excel(EXCEL_PATH) %>%
  rename(
    target_species = Species, 
    common_name    = `Common name`, 
    taxonomic_type = Type
  )

# Inject human as a standard reference baseline row since it won't be in your orthologs file
human_meta <- data.frame(
  target_species = "homo_sapiens",
  common_name    = "Human",
  taxonomic_type = "Mammalia",
  stringsAsFactors = FALSE
)
full_metadata_pool <- bind_rows(species_metadata, human_meta)

# ==============================================================================
# LEVEL 2: COMPUTE CHROMOSOMAL GROUPING & CLUSTER SPANS
# ==============================================================================
print("Level 2: Grouping by species and chromosomes to isolate cluster architectures...")

# Core logic function to divide tandem segmental duplications into separate sub-clusters
assign_sub_clusters <- function(gene_names) {
  if (length(gene_names) == 0) return(integer(0))
  sub_ids <- integer(length(gene_names))
  current_cluster <- 1
  seen_in_current_cluster <- c(gene_names[1])
  sub_ids[1] <- current_cluster
  
  if (length(gene_names) == 1) return(sub_ids)
  
  for (i in 2:length(gene_names)) {
    g <- gene_names[i]
    if (g %in% seen_in_current_cluster) {
      current_cluster <- current_cluster + 1
      seen_in_current_cluster <- c(g)
    } else {
      seen_in_current_cluster <- c(seen_in_current_cluster, g)
    }
    sub_ids[i] <- current_cluster
  }
  return(sub_ids)
}

processed_architectures <- final_dataset %>%
  filter(!is.na(ortholog_chrom), !is.na(ortholog_start), !is.na(ortholog_end)) %>%
  arrange(target_species, ortholog_chrom, ortholog_start) %>%
  group_by(target_species, ortholog_chrom) %>%
  mutate(
    sub_cluster_id = assign_sub_clusters(human_gene_name)
  ) %>%
  group_by(target_species, ortholog_chrom, sub_cluster_id) %>%
  summarise(
    max_genes_on_chrom    = n_distinct(human_gene_name),
    contains_tusc2         = ifelse("TUSC2" %in% human_gene_name, "YES", "NO"),
    cluster_start_bp       = min(ortholog_start),
    cluster_end_bp         = max(ortholog_end),
    cluster_span_kb        = (max(ortholog_end) - min(ortholog_start)) / 1000,
    left_border_gene       = human_gene_name[order(ortholog_start)][1],
    right_border_gene      = human_gene_name[order(ortholog_start)][n()],
    raw_arch_str           = paste(human_gene_name[order(ortholog_start)], collapse = " -> "),
    .groups = "drop"
  ) %>%
  mutate(
    clean_architecture = sapply(raw_arch_str, function(arch) {
      genes <- unlist(strsplit(arch, " -> "))
      forward_str <- paste(genes, collapse = " -> ")
      reverse_str <- paste(rev(genes), collapse = " -> ")
      ifelse(forward_str < reverse_str, forward_str, reverse_str)
    })
  )

# ==============================================================================
# LEVEL 3: TAXONOMIC ANNOTATION MATCHING & FACTOR LEVEL ORDERING
# ==============================================================================
print("Level 3: Merging taxonomic groups and standardizing levels...")

# Deduplicate the metadata dictionary pool to completely prevent the many-to-many warning
cleaned_metadata_pool <- full_metadata_pool %>%
  distinct(target_species, .keep_all = TRUE)

# Execute clean joining track using strict namespace paths
final_annotated_matrix <- processed_architectures %>%
  # Left join matches clean, unique rows natively
  left_join(cleaned_metadata_pool, by = "target_species") %>%
  
  # Set lineage type as a prioritized R factor mapping framework
  mutate(
    taxonomic_type = as.factor(taxonomic_type)
  ) %>%
  
  # FIX: Explicitly namespace select() using dplyr::select to override package hijacking
  dplyr::select(
    target_species, 
    common_name, 
    taxonomic_type, 
    ortholog_chrom, 
    sub_cluster_id, 
    max_genes_on_chrom, 
    cluster_span_kb, 
    left_border_gene, 
    right_border_gene, 
    clean_architecture,
    contains_tusc2
  ) %>%
  
  # Group rows logically: Lineage family, then completeness, then tightness
  arrange(
    taxonomic_type,
    desc(contains_tusc2),
    desc(max_genes_on_chrom),
    cluster_span_kb
  )

# Save the finalized, submission-ready spreadsheet
write.csv(final_annotated_matrix, "final_annotated_3p21_3_synteny.csv", row.names = FALSE)

print("SUCCESS! Lineage dataset compiled and saved cleanly.")


# Display a preview of the final annotated architecture table
print(head(final_annotated_matrix, 15))

#------------------------------------------------
#Direct Alignment with the Clean Architecture Map
#------------------------------------------------

# Clear the cached memory baseline variables to ensure clean calculations
if (exists("nn_bounds")) rm(nn_bounds)
if (exists("nn_distance_list")) rm(nn_distance_list)
if (exists("raw_full_coords")) rm(raw_full_coords)

library(dplyr)
library(tidyr)
library(stringr)

# ==============================================================================
# STEP 1: INTERACTIVE FILE SELECTORS (NO FIXED DIRECTORY STRINGS REQUIRED)
# ==============================================================================
print("Step 1: Please select your files in the pop-up windows...")

message("--> Click and select your summary file: 'final_annotated_3p21_3_synteny.csv'")
INPUT_CSV <- file.choose()

message("--> Click and select your coordinates file: 'all_metazoan_3p21_3_raw_with_duplicates.csv'")
RAW_COORDS <- file.choose()

OUTPUT_CSV <- file.path(dirname(INPUT_CSV), "final_annotated_3p21_3_synteny_augmented.csv")

annotated_clusters <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)
raw_coordinates     <- read.csv(RAW_COORDS, stringsAsFactors = FALSE)

# --- FIX: RECONSTRUCT TRUNCATED EXCEL/CSV FILE COLUMN HEADERS NATIVELY ---
if ("target_spe" %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(target_species = target_spe)
if ("common_r" %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(common_name = common_r)
if ("taxonomic." %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(taxonomic_type = taxonomic.)
if ("ortholog_c" %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(ortholog_chrom = ortholog_c)
if ("sub_cluste" %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(sub_cluster_id = sub_cluste)
if ("max_genes" %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(max_genes_on_chrom = max_genes)
if ("cluster_sp" %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(cluster_span_kb = cluster_sp)
if ("left_border" %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(left_border_gene = left_border)
if ("right_borde" %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(right_border_gene = right_borde)
if ("clean_arch" %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(clean_architecture = clean_arch)
if ("contains_t" %in% colnames(annotated_clusters)) annotated_clusters <- annotated_clusters %>% rename(contains_tusc2 = contains_t)

if ("target_spe" %in% colnames(raw_coordinates)) raw_coordinates <- raw_coordinates %>% rename(target_species = target_spe)
if ("ortholog_c" %in% colnames(raw_coordinates)) raw_coordinates <- raw_coordinates %>% rename(ortholog_chrom = ortholog_c)
if ("start" %in% colnames(raw_coordinates)) raw_coordinates <- raw_coordinates %>% rename(ortholog_start = start)
if ("end" %in% colnames(raw_coordinates)) raw_coordinates <- raw_coordinates %>% rename(ortholog_end = end)

annotated_clusters$ortholog_chrom <- as.character(annotated_clusters$ortholog_chrom)
raw_coordinates$ortholog_chrom     <- as.character(raw_coordinates$ortholog_chrom)

# ==============================================================================
# CORRECTED STEP 2: DYNAMIC BASELINE THRESHOLD CALIBRATION
# ==============================================================================
print("Step 2: Calibrating baseline standard deviation variance matrix...")

# Dynamic fallback: Identify the maximum available gene count present in your data
max_available_genes <- max(annotated_clusters$max_genes_on_chrom, na.rm = TRUE)
target_baseline_occupancy <- min(8, max_available_genes)

message(paste("--> Dynamic Baseline Target Selected:", target_baseline_occupancy, "genes"))

full_occupancy_targets <- annotated_clusters %>%
  filter(max_genes_on_chrom >= target_baseline_occupancy) %>%
  dplyr::select(target_species, ortholog_chrom)

raw_full_coords <- raw_coordinates %>%
  inner_join(full_occupancy_targets, by = c("target_species", "ortholog_chrom")) %>%
  filter(!is.na(ortholog_start), !is.na(ortholog_end))

unique_species <- unique(raw_full_coords$target_species)
nn_distance_list <- list()

for (sp in unique_species) {
  sp_df <- raw_full_coords %>% filter(target_species == sp)
  genes <- sp_df$human_gene_name
  coords <- (sp_df$ortholog_start + sp_df$ortholog_end) / 2 
  names(coords) <- genes
  
  # Measure distance using the dynamic highest available occupancy limit
  if (length(coords) < target_baseline_occupancy) next
  
  for (g1 in names(coords)) {
    other_coords <- coords[names(coords) != g1]
    if (length(other_coords) == 0) next
    min_dist_kb <- min(abs(coords[g1] - other_coords), na.rm = TRUE) / 1000
    
    nn_distance_list[[length(nn_distance_list) + 1]] <- data.frame(
      target_species = sp, 
      ortholog_chrom = sp_df$ortholog_chrom,
      gene_name      = g1, 
      nn_distance_kb = min_dist_kb, 
      stringsAsFactors = FALSE
    )
  }
}

# Compile bounds safely
nn_bounds <- bind_rows(nn_distance_list) %>%
  group_by(gene_name) %>%
  summarise(mean_kb = mean(nn_distance_kb, na.rm = TRUE), sd_kb = sd(nn_distance_kb, na.rm = TRUE), .groups = "drop") %>%
  mutate(threshold_2.5sigma = mean_kb + (2.5 * sd_kb))

print("--- Calibrated Thresholds Table ---")
print(nn_bounds)

# ==============================================================================
# STEP 3: DYNAMIC RECALCULATION LOOP FOR UPDATED CLUSTERS
# ==============================================================================
print("Step 3: Calculating and appending the updated cluster metadata...")

inflating_genes_vec           <- character(nrow(annotated_clusters))
max_genes_without_outliers    <- integer(nrow(annotated_clusters))
span_without_outliers_vec     <- numeric(nrow(annotated_clusters))
left_border_without_outliers  <- character(nrow(annotated_clusters))
right_border_without_outliers <- character(nrow(annotated_clusters))
architecture_without_outliers <- character(nrow(annotated_clusters))
outlier_metrics_summary_vec   <- character(nrow(annotated_clusters))

for (r in 1:nrow(annotated_clusters)) {
  current_sp   <- annotated_clusters$target_species[r]
  current_chrom <- annotated_clusters$ortholog_chrom[r]
  arch_string   <- annotated_clusters$clean_architecture[r]
  
  valid_genes <- unlist(strsplit(arch_string, " -> "))
  valid_genes <- gsub(" ", "", valid_genes) 
  
  sp_genes <- raw_coordinates %>%
    filter(target_species == current_sp, ortholog_chrom == current_chrom, human_gene_name %in% valid_genes, !is.na(ortholog_start)) %>%
    arrange(ortholog_start)
  
  if (nrow(sp_genes) < 2) {
    inflating_genes_vec[r]           <- "None"
    max_genes_without_outliers[r]    <- annotated_clusters$max_genes_on_chrom[r]
    span_without_outliers_vec[r]     <- annotated_clusters$cluster_span_kb[r]
    left_border_without_outliers[r]  <- annotated_clusters$left_border_gene[r]
    right_border_without_outliers[r] <- annotated_clusters$right_border_gene[r]
    architecture_without_outliers[r] <- annotated_clusters$clean_architecture[r]
    outlier_metrics_summary_vec[r]   <- "Incomplete Cluster"
    next
  }
  
  coords <- (sp_genes$ortholog_start + sp_genes$ortholog_end) / 2
  names(coords) <- sp_genes$human_gene_name
  
  flagged_outliers <- c()
  metrics_strings  <- c()
  
  for (g in names(coords)) {
    other_coords <- coords[names(coords) != g]
    if (length(other_coords) == 0) next
    obs_dist <- min(abs(coords[g] - other_coords), na.rm = TRUE) / 1000
    
    g_bound <- nn_bounds %>% filter(gene_name == g)
    if (nrow(g_bound) == 1) {
      if (obs_dist > g_bound$threshold_2.5sigma) {
        flagged_outliers <- c(flagged_outliers, g)
        metrics_strings  <- c(metrics_strings, paste0(g, "(Obs:", round(obs_dist, 1), "kb/Cut:", round(g_bound$threshold_2.5sigma, 1), "kb)"))
      }
    }
  }
  
  clean_genes_df <- sp_genes %>% filter(!(human_gene_name %in% flagged_outliers))
  
  if (nrow(clean_genes_df) >= 2) {
    new_span  <- (max(clean_genes_df$ortholog_end) - min(clean_genes_df$ortholog_start)) / 1000
    new_left  <- clean_genes_df$human_gene_name[which.min(clean_genes_df$ortholog_start)]
    new_right <- clean_genes_df$human_gene_name[which.max(clean_genes_df$ortholog_start)]
    
    raw_string  <- paste(clean_genes_df$human_gene_name[order(clean_genes_df$ortholog_start)], collapse = " -> ")
    genes_split <- unlist(strsplit(raw_string, " -> "))
    forward_str <- paste(genes_split, collapse = " -> ")
    reverse_str <- paste(rev(genes_split), collapse = " -> ")
    new_arch    <- ifelse(forward_str < reverse_str, forward_str, reverse_str)
    
    new_count   <- n_distinct(clean_genes_df$human_gene_name)
  } else {
    new_span  <- 0
    new_left  <- "None"
    new_right <- "None"
    new_arch  <- "Collapsed"
    new_count <- nrow(clean_genes_df)
  }
  
  inflating_genes_vec[r]           <- ifelse(length(flagged_outliers) == 0, "None", paste(flagged_outliers, collapse = ", "))
  max_genes_without_outliers[r]    <- ifelse(length(flagged_outliers) == 0, annotated_clusters$max_genes_on_chrom[r], new_count)
  span_without_outliers_vec[r]     <- round(ifelse(length(flagged_outliers) == 0, annotated_clusters$cluster_span_kb[r], new_span), 2)
  left_border_without_outliers[r]  <- ifelse(length(flagged_outliers) == 0, annotated_clusters$left_border_gene[r], new_left)
  right_border_without_outliers[r] <- ifelse(length(flagged_outliers) == 0, annotated_clusters$right_border_gene[r], new_right)
  architecture_without_outliers[r] <- ifelse(length(flagged_outliers) == 0, annotated_clusters$clean_architecture[r], new_arch)
  outlier_metrics_summary_vec[r]   <- ifelse(length(metrics_strings) == 0, "All Genes Within Baseline 2.5 Sigma", paste(metrics_strings, collapse = " | "))
}

print("SUCCESS! Recalculation loop complete.")


# --- STEP 4: AUGMENT ORIGINAL DATAFRAME MATRIX AND EXPORT ---
print("Step 4: Merging results and exporting your augmented CSV table...")

augmented_final_matrix <- annotated_clusters %>%
  mutate(
    inflating_genes               = inflating_genes_vec,
    updated_max_genes_on_chrom    = max_genes_without_outliers,
    updated_cluster_span_kb       = span_without_outliers_vec,
    updated_left_border_gene      = left_border_without_outliers,
    updated_right_border_gene     = right_border_without_outliers,
    updated_clean_architecture    = architecture_without_outliers,
    outlier_metrics_summary       = outlier_metrics_summary_vec
  )

write.csv(augmented_final_matrix, OUTPUT_CSV, row.names = FALSE)


#------------------------------------------------------------------
#NCBI ORTHOLOGS RETRIEVAL
#------------------------------------------------------------------
# Install the required offline taxonomy indexing library if not already present
if (!requireNamespace("taxizedb", quietly = TRUE)) install.packages("taxizedb")

library(dplyr)
library(tidyr)
library(taxizedb)

# --- FILE PATHWAY CONFIGURATIONS ---
TSV_PATH   <- "C:\\Users\\...\\3p21.3_cluster.tsv"
OUTPUT_CSV <- "C:\\Users\\...\\final_augmented_ncbi_synteny.csv"

# Core 8 functional tumor suppressor target anchors
TARGET_GENES <- c("RASSF1", "TUSC2", "ZMYND10", "NPRL2", "HYAL2", "CACNA2D2", "CYB561D2", "TMEM115")

print("Step 1: Loading raw NCBI multi-species TSV data matrix...")
if (!file.exists(TSV_PATH)) stop("Target file missing. Verify directory location.")

# Load the file, explicitly declaring tab-separation and mapping names row-by-row
raw_ncbi_data <- read.delim(TSV_PATH, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)

# Space-sensitive header mapping matching your exact file structure
colnames(raw_ncbi_data) <- c("GeneID", "Symbol", "Chromosome", "Accession", "Start", "Stop", "TaxID")

# ==============================================================================
# LEVEL 1: DISTANCE-BASED WINDOWING PIPELINE
# ==============================================================================
print("Level 1: Unifying genes within 1.5 Mb windows to fix fragmented clusters...")

cleaned_ncbi_pool <- raw_ncbi_data %>%
  mutate(
    Symbol = toupper(gsub(" ", "", Symbol)),
    Start  = as.numeric(Start),
    Stop   = as.numeric(Stop),
    TaxID  = as.character(TaxID)
  ) %>%
  filter(Symbol %in% TARGET_GENES, !is.na(Start), !is.na(Stop))

# Merges adjacent targets if they sit within 1.5 Mb of each other on the same scaffold
assign_clusters_by_distance <- function(starts, stops) {
  if (length(starts) == 0) return(integer(0))
  cluster_ids <- integer(length(starts))
  current_cluster <- 1
  cluster_ids <- current_cluster
  if (length(starts) == 1) return(cluster_ids)
  for (i in 2:length(starts)) {
    if ((starts[i] - stops[i-1]) > 1500000) {
      current_cluster <- current_cluster + 1
    }
    cluster_ids[i] <- current_cluster
  }
  return(cluster_ids)
}

processed_architectures <- cleaned_ncbi_pool %>%
  arrange(TaxID, Accession, Start) %>%
  group_by(TaxID, Accession) %>%
  mutate(sub_cluster_id = assign_clusters_by_distance(Start, Stop)) %>%
  group_by(TaxID, Accession, sub_cluster_id) %>%
  summarise(
    max_genes_on_chrom    = n_distinct(Symbol),
    contains_tusc2         = ifelse("TUSC2" %in% Symbol, "YES", "NO"),
    cluster_start_bp       = min(Start),
    cluster_end_bp         = max(Stop),
    cluster_span_kb        = (max(Stop) - min(Start)) / 1000,
    left_border_gene       = Symbol[which.min(Start)],
    right_border_gene      = Symbol[which.max(Start)],
    raw_arch_str           = paste(Symbol[order(Start)], collapse = " -> "),
    .groups = "drop"
  ) %>%
  mutate(
    clean_architecture = sapply(raw_arch_str, function(arch) {
      genes <- unlist(strsplit(arch, " -> "))
      forward_str <- paste(genes, collapse = " -> ")
      reverse_str <- paste(rev(genes), collapse = " -> ")
      ifelse(forward_str < reverse_str, forward_str, reverse_str)
    })
  )

# ==============================================================================
# LEVEL 2: COMPREHENSIVE LOCAL OFFLINE TAXONOMY RESOLUTION
# ==============================================================================
print("Level 2: Building local taxonomy database using high-speed local indexing...")

unique_taxids <- unique(processed_architectures$TaxID)
unique_taxids_num <- as.integer(unique_taxids[!is.na(unique_taxids)])

# Pull the master classification lineage nodes completely offline using taxizedb
print("   -> Fetching scientific classifications from local database file records...")
local_tax_data <- taxizedb::classification(unique_taxids_num, db = "ncbi")

taxonomy_dictionary <- list()

for (tid in names(local_tax_data)) {
  df <- local_tax_data[[tid]]
  if (!is.null(df) && is.data.frame(df)) {
    # Extract the true terminal scientific name
    scientific_name <- df$name[nrow(df)]
    
    # Trace through parent nodes to find the Class rank
    full_lineage <- paste(df$name, collapse = " | ")
    tax_class <- "Other Lineage"
    
    if (grepl("Mammalia", full_lineage, ignore.case=TRUE)) tax_class <- "Mammalia"
    else if (grepl("Amphibia", full_lineage, ignore.case=TRUE)) tax_class <- "Amphibia"
    else if (grepl("Lepidosauria|Reptilia|Testudines|Crocodylia|Squamata", full_lineage, ignore.case=TRUE)) tax_class <- "Reptilia"
    else if (grepl("Aves", full_lineage, ignore.case=TRUE)) tax_class <- "Aves"
    else if (grepl("Actinopterygii|Chondrichthyes|Teleostei", full_lineage, ignore.case=TRUE)) tax_class <- "Actinopterygii"
    else if (grepl("Insecta|Arthropoda", full_lineage, ignore.case=TRUE)) tax_class <- "Insecta (Invertebrate)"
    else if (grepl("Nematoda", full_lineage, ignore.case=TRUE)) tax_class <- "Nematoda (Invertebrate)"
    else if (grepl("Mollusca|Echinodermata|Cnidaria", full_lineage, ignore.case=TRUE)) tax_class <- "Invertebrate"
    
    taxonomy_dictionary[[tid]] <- data.frame(
      TaxID          = as.character(tid),
      target_species = as.character(scientific_name),
      common_name    = as.character(scientific_name),
      taxonomic_type = tax_class,
      stringsAsFactors = FALSE
    )
  }
}

taxonomy_df <- bind_rows(taxonomy_dictionary)

# ==============================================================================
# LEVEL 3: 2.5 SIGMA NEAREST-NEIGHBOR SPATIAL BASELINE
# ==============================================================================
print("Level 3: Calibrating baseline standard deviation variance matrix...")

max_observed_occupancy <- max(processed_architectures$max_genes_on_chrom)
full_clusters <- processed_architectures %>% 
  filter(max_genes_on_chrom >= min(7, max_observed_occupancy))

raw_full_coords <- cleaned_ncbi_pool %>%
  inner_join(full_clusters %>% dplyr::select(TaxID, Accession), by = c("TaxID", "Accession"))

nn_distance_list <- list()

for (t_id in unique(raw_full_coords$TaxID)) {
  sp_df <- raw_full_coords %>% filter(TaxID == t_id)
  genes <- sp_df$Symbol
  coords <- (sp_df$Start + sp_df$Stop) / 2 
  names(coords) <- genes
  
  if (length(coords) < 2) next
  
  for (g1 in names(coords)) {
    other_coords <- coords[names(coords) != g1]
    if (length(other_coords) == 0) next
    min_dist_kb <- min(abs(coords[g1] - other_coords), na.rm = TRUE) / 1000
    nn_distance_list[[length(nn_distance_list) + 1]] <- data.frame(
      TaxID = t_id, gene_name = g1, nn_distance_kb = min_dist_kb, stringsAsFactors = FALSE
    )
  }
}

nn_bounds <- bind_rows(nn_distance_list) %>%
  group_by(gene_name) %>%
  summarise(mean_kb = mean(nn_distance_kb), sd_kb = sd(nn_distance_kb), .groups = "drop") %>%
  mutate(threshold_2.5sigma = mean_kb + (2.5 * sd_kb))

# ==============================================================================
# ALGORITHM BLOCK: LEVEL 4 - RE-COMPILING CLUSTER TRACKS WITHOUT OUTLIERS
# ==============================================================================

# Initialize empty destination vectors matching row dimension of processed_architectures
inflating_genes_vec           <- character(length = nrow(processed_architectures))
updated_max_genes_vec         <- integer(length = nrow(processed_architectures))
updated_span_vec              <- numeric(length = nrow(processed_architectures))
updated_left_border_gene_vec  <- character(length = nrow(processed_architectures))
updated_right_border_gene_vec <- character(length = nrow(processed_architectures))
architecture_without_outliers <- character(length = nrow(processed_architectures))
outlier_metrics_summary_vec   <- character(length = nrow(processed_architectures))

# Begin physical row loop matrix scan
for (r in 1:nrow(processed_architectures)) {
  
  # Isolate active keys from base cluster tracking row
  t_id   <- processed_architectures$TaxID[r]
  acc    <- processed_architectures$Accession[r]
  arch   <- processed_architectures$clean_architecture[r]
  
  # Split text track arrow tokens into atomic family character strings
  valid_genes <- unlist(strsplit(arch, split = " -> "))
  
  # Filter raw coordinate table down exclusively to current locus track members
  sp_genes <- cleaned_ncbi_pool %>%
    filter(TaxID == t_id & Accession == acc & Symbol %in% valid_genes) %>%
    arrange(Start)
  
  # Guard condition handler for partial or single-gene loops
  if (nrow(sp_genes) < 2) {
    inflating_genes_vec[r]           <- "None"
    updated_max_genes_vec[r]         <- processed_architectures$max_genes_on_chrom[r]
    updated_span_vec[r]              <- processed_architectures$cluster_span_kb[r]
    updated_left_border_gene_vec[r]  <- processed_architectures$left_border_gene[r]
    updated_right_border_gene_vec[r] <- processed_architectures$right_border_gene[r]
    architecture_without_outliers[r] <- processed_architectures$clean_architecture[r]
    outlier_metrics_summary_vec[r]   <- "Incomplete"
    next # Jump instantly to subsequent row index block
  }
  
  # Force coordinates to numeric double vectors to prevent 2.14Gb integer overflow
  coords <- (as.numeric(sp_genes$Start) + as.numeric(sp_genes$Stop)) / 2
  names(coords) <- sp_genes$Symbol
  
  flagged_outliers <- character()
  metrics_strings  <- character()
  
  # Execute Nearest-Neighbor structural adjacency map sweep
  for (g in names(coords)) {
    other_coords <- coords[names(coords) != g]
    if (length(other_coords) == 0) next
    
    # Calculate coordinate physical base-pair span to the single closest neighbor element
    obs_dist <- min(abs(coords[g] - other_coords), na.rm = TRUE) / 1000
    
    # Handle NA or calculation exceptions natively
    if (is.na(obs_dist) | is.infinite(obs_dist)) next
    
    # Match against the pre-calibrated 2.5-sigma database background thresholds table
    g_bound <- nn_bounds %>% filter(gene_name == g)
    if (nrow(g_bound) == 1) {
      if (!is.na(g_bound$threshold_2.5sigma) & obs_dist > g_bound$threshold_2.5sigma) {
        flagged_outliers <- c(flagged_outliers, g)
        metrics_strings  <- c(metrics_strings, paste0(g, "(Obs:", round(obs_dist, 1), "kb/Cut:", round(g_bound$threshold_2.5sigma, 1), "kb)"))
      }
    }
  }
  
  # Subtract discovered outlying entries from local data matrix pool
  clean_genes_df <- sp_genes %>% filter(!(Symbol %in% flagged_outliers))
  
  # Recalculate whole cluster boundaries track metrics
  if (nrow(clean_genes_df) >= 2) {
    # Enforce numeric casting inside maximum evaluations to protect against overflow subtractions
    new_span  <- (max(as.numeric(clean_genes_df$Stop)) - min(as.numeric(clean_genes_df$Start))) / 1000
    new_left  <- clean_genes_df$Symbol[which.min(clean_genes_df$Start)]
    new_right <- clean_genes_df$Symbol[which.max(clean_genes_df$Start)]
    
    # Build direction-agnostic architecture text token string mapping layout
    raw_string  <- paste(clean_genes_df$Symbol[order(clean_genes_df$Start)], collapse = " -> ")
    genes_split <- unlist(strsplit(raw_string, split = " -> "))
    forward_str <- paste(genes_split, collapse = " -> ")
    reverse_str <- paste(rev(genes_split), collapse = " -> ")
    new_arch    <- ifelse(forward_str < reverse_str, forward_str, reverse_str)
    new_count   <- length(unique(clean_genes_df$Symbol))
  } else {
    # Fallback parameters if cluster components dissolve below boundary conditions
    new_span <- 0; new_left <- "None"; new_right <- "None"; new_arch <- "Collapsed"; new_count <- nrow(clean_genes_df)
  }
  
  # Inject processed loop values into vector storage locations matching cell index r
  inflating_genes_vec[r]           <- ifelse(length(flagged_outliers) == 0, "None", paste(flagged_outliers, collapse = ", "))
  updated_max_genes_vec[r]         <- ifelse(length(flagged_outliers) == 0, processed_architectures$max_genes_on_chrom[r], new_count)
  updated_span_vec[r]              <- round(ifelse(length(flagged_outliers) == 0, processed_architectures$cluster_span_kb[r], new_span), 2)
  updated_left_border_gene_vec[r]  <- ifelse(length(flagged_outliers) == 0, processed_architectures$left_border_gene[r], new_left)
  updated_right_border_gene_vec[r] <- ifelse(length(flagged_outliers) == 0, processed_architectures$right_border_gene[r], new_right)
  architecture_without_outliers[r] <- ifelse(length(flagged_outliers) == 0, processed_architectures$clean_architecture[r], new_arch)
  outlier_metrics_summary_vec[r]   <- ifelse(length(metrics_strings) == 0, "All Genes Within Baseline 2.5 Sigma", paste(metrics_strings, collapse = " | "))
}

# ==============================================================================
# ALGORITHM BLOCK: LEVEL 5 - MATRIX INTEGRATION AND FINAL CSV EXPORT
# ==============================================================================

print("Level 5: Binding tracking columns and executing final export...")

# Append Level 4 calculated tracking vector arrays directly onto the base data structure fields
augmented_matrix <- processed_architectures %>%
  mutate(
    inflating_genes               = inflating_genes_vec,
    updated_max_genes_on_chrom    = updated_max_genes_vec,
    updated_cluster_span_kb       = updated_span_vec,
    updated_left_border_gene      = updated_left_border_gene_vec,
    updated_right_border_gene     = updated_right_border_gene_vec,
    updated_clean_architecture    = architecture_without_outliers,
    outlier_metrics_summary       = outlier_metrics_summary_vec
  )

# Map dynamic local taxonomy lineage tables via TaxID key matching properties
if (exists("taxonomy_df") & nrow(taxonomy_df) > 0) {
  augmented_matrix <- augmented_matrix %>%
    left_join(taxonomy_df, by = "TaxID") %>%
    mutate(
      target_species = ifelse(is.na(target_species), paste0("NCBI_TaxID_", TaxID), target_species),
      common_name    = ifelse(is.na(common_name), paste0("Species_", TaxID), common_name),
      taxonomic_type = ifelse(is.na(taxonomic_type), "Other Lineage", taxonomic_type)
    )
} else {
  # Fallback configuration matrix blocks if local database indexing returns empty
  augmented_matrix <- augmented_matrix %>%
    mutate(
      target_species = paste0("TaxID_", TaxID), 
      common_name    = "Unknown Species", 
      taxonomic_type = "Unresolved"
    )
}

# Arrange fields systematically into publication formatting rules using explicit namespaces
final_output <- augmented_matrix %>%
  dplyr::select(
    target_species, common_name, taxonomic_type, ortholog_chrom = Accession, sub_cluster_id,
    max_genes_on_chrom, cluster_span_kb, left_border_gene, right_border_gene, clean_architecture, contains_tusc2,
    inflating_genes, updated_max_genes_on_chrom, updated_cluster_span_kb, updated_left_border_gene, 
    updated_right_border_gene, updated_clean_architecture, outlier_metrics_summary
  ) %>%
  # Re-order rows systematically, sorting first by lineage family and then by cluster completeness tightness
  arrange(taxonomic_type, updated_cluster_span_kb)

# Export the completed augmented matrix to file system directory
write.csv(final_output, file = OUTPUT_CSV, row.names = FALSE)

print("SUCCESS! Locus data augmented. Recalculated columns appended.")
print(table(final_output$taxonomic_type))

#-----------------------------------------------------------------
#Structural completeness of 3p21.3 cluster in Vertebrates (FIGURE)
#-----------------------------------------------------------------

# Install required visualization packages if not already present in the workspace
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

# --- FILE PATHWAY CONFIGURATIONS ---
INPUT_XLSX  <- "C:\\Users\\...\\3p21.3_cluster_mapping_vert.xlsx"
OUTPUT_PNG  <- "C:\\Users\\...\\Figure_3p21_locus_occupancy_completeness.png"

# ==============================================================================
# STEP 1: LOAD MASTER DATA FROM THE COMPREHENSIVE EXCEL WORKBOOK
# ==============================================================================
print("Step 1: Loading master data directly from your vertical mapping workbook...")
if (!file.exists(INPUT_XLSX)) stop("Master Excel workbook missing. Verify folder path.")

# Read the active sheet
master_data <- read_excel(INPUT_XLSX)

print("Step 2: Re-mapping column names out of formula evaluation boundaries...")

# FIX: Completely separate the header sniffer from case_when to prevent evaluation crashes
if ("raw_unitless_occupancy" %in% colnames(master_data)) {
  master_data$raw_unitless_occupancy <- as.numeric(master_data$raw_unitless_occupancy)
} else {
  # Locate whatever naming variant you used for the coding footprint sizes column
  detected_gene_len_col <- NA
  if ("total_gene_length_kb" %in% colnames(master_data)) {
    detected_gene_len_col <- "total_gene_length_kb"
  } else if ("total_gene_length" %in% colnames(master_data)) {
    detected_gene_len_col <- "total_gene_length"
  } else if ("gene_length_kb" %in% colnames(master_data)) {
    detected_gene_len_col <- "gene_length_kb"
  } else {
    matched <- colnames(master_data)[grepl("length|footprint|occupancy", colnames(master_data), ignore.case = TRUE)]
    if (length(matched) > 0) detected_gene_len_col <- matched[1]
  }
  
  if (!is.na(detected_gene_len_col)) {
    message(paste("--> Dynamically mapping coding sequence column:", detected_gene_len_col))
    # Calculate the ratio directly in the base dataframe safely
    master_data$raw_unitless_occupancy <- as.numeric(master_data[[detected_gene_len_col]]) / as.numeric(master_data$cluster_span_kb)
  } else {
    # If no length column exists, default safely to 1.0 to let the rest of the script run
    master_data$raw_unitless_occupancy <- 1.0
  }
}

clean_plot_data <- master_data %>%
  mutate(
    # Clean string trailing spaces from lineage categories
    taxonomic_type = trimws(as.character(taxonomic_type)),
    # Ensure maximum gene count factors handle integers perfectly
    updated_max_genes_on_chrom = as.factor(as.integer(updated_max_genes_on_chrom))
  ) %>%
  mutate(
    # Use case-insensitive partial string matching (grepl) to unify fish sub-labels
    taxonomic_type_clean = case_when(
      grepl("Actinopterygii", taxonomic_type, ignore.case = TRUE) ~ "Actinopterygii",
      grepl("Amphibia", taxonomic_type, ignore.case = TRUE)       ~ "Amphibia",
      grepl("Reptilia", taxonomic_type, ignore.case = TRUE)       ~ "Reptilia",
      grepl("Aves", taxonomic_type, ignore.case = TRUE)           ~ "Aves",
      grepl("Mammalia", taxonomic_type, ignore.case = TRUE)       ~ "Mammalia",
      TRUE ~ taxonomic_type
    )
  ) %>%
  # Filter down strictly to your 5 target publication classes
  filter(taxonomic_type_clean %in% c("Actinopterygii", "Amphibia", "Reptilia", "Aves", "Mammalia")) %>%
  mutate(
    raw_unitless_occupancy = pmin(1, pmax(0, as.numeric(raw_unitless_occupancy))),
    # Locked sequence timeline factor: Fish -> Amphibia -> Reptilia -> Aves -> Mammalia
    taxonomic_type_clean = factor(taxonomic_type_clean, levels = c("Actinopterygii", "Amphibia", "Reptilia", "Aves", "Mammalia"))
  ) %>%
  filter(!is.na(raw_unitless_occupancy))

# Nature/SciRep compliant palette configurations matching corrected factor sequences
lineage_colors <- c(
  "Actinopterygii" = "#2b8cbe", 
  "Amphibia"       = "#74c476", 
  "Reptilia"       = "#fe9929", 
  "Aves"           = "#df65b0", 
  "Mammalia"       = "#8856a7"
)

# ==============================================================================
# PANEL: LOCUS COMPLETENESS MATRIX (STACKED GENE FREQUENCY CHARTS)
# ==============================================================================
print("Panel B: Computing gene occupancy counts across target classes...")

plot_completeness <- ggplot(clean_plot_data, aes(x = taxonomic_type_clean, fill = updated_max_genes_on_chrom)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.2, width = 0.65) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_viridis_d(option = "D", direction = 1, name = "Core Genes\nIn Cluster") +
  labs(
    x = NULL,
    y = "Proportion of Documented Loci (%)",
    title = "B. Structural Completeness & Core Gene Retention Profiles"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 10)),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(face = "bold", size = 10),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# ==============================================================================
# WORKBOOK JOINING: PATCHWORK ASSEMBLY AND HIGH-RES EXPORT
# ==============================================================================
print("Workbook Factory: Integrating panel tracks and exporting high-resolution PNG...")

final_figure_panel <- plot_occupancy / plot_completeness + 
  plot_layout(heights = c(1.1, 1.0))

ggsave(
  filename = OUTPUT_PNG,
  plot = final_figure_panel,
  width = 8.5,
  height = 9.0,
  dpi = 400, 
  device = "png"
)

print("==============================================================================")
print("GRAPH SUCCESS! Multi-panel distribution chart compiled successfully.")
print("==============================================================================")
print(table(clean_plot_data$taxonomic_type_clean))


#-----------------------------------------------------------------------
#FIGURE FOR 3p21.3 CLUSTER IN VERTEBRATES
#-----------------------------------------------------------------------
# Install required graphics visualization libraries if missing
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# --- FILE PATHWAY CONFIGURATIONS ---
OUTPUT_PNG <- "C:\\Users\\...\\Figure_3p21_Locus_Synteny_Map.png"

# Load the comprehensive coordinate matrix with duplicated segments
wide_data <- data.frame(
  Gene = c(
    "RASSF1_a", "RASSF1_b", "TUSC2_a", "TUSC2_b", "ZMYND10_a", "NPRL2_a", 
    "HYAL2_a", "HYAL2_b", "CACNA2D2_a", "CACNA2D2_b", "CYB561D2_a", "TMEM115_a"
  ),
  Gene_Family = c(
    "RASSF1", "RASSF1", "TUSC2", "TUSC2", "ZMYND10", "NPRL2", 
    "HYAL2", "HYAL2", "CACNA2D2", "CACNA2D2", "CYB561D2", "TMEM115"
  ),
  
  Homo_sapiens_Start = c(50329788, NA, 50324909, NA, 50341112, 50347330, 50317808, NA, 50362613, NA, 50350862, 50354750),
  Homo_sapiens_End   = c(50340836, NA, 50328223, NA, 50345732, 50350775, 50322745, NA, 50504244, NA, 50354069, 50359521),
  
  Reptilia_Start     = c(45061815, NA, 45086402, NA, 45030628, 45015077, 45119727, NA, 44564409, NA, 44998299, 44997911),
  Reptilia_End       = c(45084051, NA, 45097078, NA, 45054735, 45030796, 45131388, NA, 44990191, NA, 45016750, 45007243),
  
  Aves_Start         = c(16970994, NA, 16975587, NA, 16962373, 16959916, 16982035, NA, 16834463, NA, 16958272, 16955260),
  Aves_End           = c(16975030, NA, 16979814, NA, 16970224, 16963505, 16984854, NA, 16954589, NA, 16959798, 16958136),
  
  Amphibia_Start     = c(127211953, NA, 127203033, NA, 127245488, 127912033, 127147727, NA, 127723234, NA, 127908596, 127901611),
  Amphibia_End       = c(127240683, NA, 127208433, NA, 127259149, 127925217, 127168488, NA, 127886445, NA, 127911998, 127905980),
  
  Pre_WGD_Start      = c(37863931, NA, 37873776, 13239813, 13218864, 13215960, 13243503, 37885268, 13367690, 37925000, 37858118, NA),
  Pre_WGD_End        = c(37873600, NA, 37876800, 13241681, 13229983, 13218817, 13247287, 37893303, 13544284, 37969855, 37860955, NA),
  
  Post_WGD_Start     = c(2197778,  2189555,  2195062,  NA,       772372,   769086,   815306,   NA,       850656,   NA,       2204003,  NA),
  Post_WGD_End       = c(2203897,  2194364,  2197439,  NA,       777860,   772232,   820534,   NA,       956943,   NA,       2206672,  NA),
  stringsAsFactors   = FALSE
)

print("Step 1: Reshaping data matrices...")
long_starts <- wide_data %>%
  dplyr::select(Gene, Gene_Family, contains("Start")) %>%
  pivot_longer(cols = -c(Gene, Gene_Family), names_to = "Lineage", values_to = "Start") %>%
  mutate(Lineage = gsub("_Start", "", Lineage))

long_ends <- wide_data %>%
  dplyr::select(Gene, Gene_Family, contains("End")) %>%
  pivot_longer(cols = -c(Gene, Gene_Family), names_to = "Lineage", values_to = "End") %>%
  mutate(Lineage = gsub("_End", "", Lineage))

synteny_matrix <- left_join(long_starts, long_ends, by = c("Gene", "Gene_Family", "Lineage")) %>%
  filter(!is.na(Start))

print("Step 2: Execating adaptive spacer contraction...")
standardized_tracks <- synteny_matrix %>%
  mutate(
    Cluster_Group = case_when(
      Lineage == "Pre_WGD" & Start < 20000000 ~ "C1",
      Lineage == "Pre_WGD" & Start >= 20000000 ~ "C2",
      Lineage == "Post_WGD" & Start < 1500000  ~ "C1",
      Lineage == "Post_WGD" & Start >= 1500000 ~ "C2",
      TRUE ~ "C1"
    )
  ) %>%
  group_by(Lineage, Cluster_Group) %>%
  mutate(
    sub_min = min(c(Start, End)),
    local_offset_start = (Start - sub_min) / 1000,
    local_offset_end   = (End - sub_min) / 1000
  ) %>%
  group_by(Lineage) %>%
  mutate(
    c1_max_kb = max(local_offset_end[Cluster_Group == "C1"]),
    Relative_Start_kb = ifelse(Cluster_Group == "C2", local_offset_start + c1_max_kb + 150, local_offset_start),
    Relative_End_kb   = ifelse(Cluster_Group == "C2", local_offset_end + c1_max_kb + 150, local_offset_end),
    Relative_Mid_kb   = (Relative_Start_kb + Relative_End_kb) / 2
  ) %>%
  ungroup() %>%
  mutate(
    Lineage_Label = case_when(
      Lineage == "Pre_WGD"       ~ "Actinopterygii (pre-WGD)",
      Lineage == "Post_WGD"      ~ "Actinopterygii (post-WGD)",
      Lineage == "Amphibia"      ~ "Amphibia",
      Lineage == "Reptilia"      ~ "Reptilia",
      Lineage == "Aves"          ~ "Aves",
      Lineage == "Homo_sapiens"  ~ "Homo sapiens"
    ),
    Lineage_Label = factor(Lineage_Label, levels = c(
      "Actinopterygii (pre-WGD)", "Actinopterygii (post-WGD)", 
      "Amphibia", "Reptilia", "Aves", "Homo sapiens"
    ))
  )

gene_palette <- c(
  "RASSF1"   = "#33a02c", "TUSC2"    = "#ff7f00", "ZMYND10"  = "#1f78b4", 
  "NPRL2"    = "#e31a1c", "HYAL2"    = "#6a3d9a", "CACNA2D2" = "#b15928", 
  "CYB561D2" = "#a6cee3", "TMEM115"  = "#b2df8a"
)

break_markers <- standardized_tracks %>%
  filter(Cluster_Group == "C2") %>%
  group_by(Lineage_Label) %>%
  summarise(Break_X = min(Relative_Start_kb) - 75, .groups = "drop")

# Build a discrete summary data frame of unique tracks to draw custom tight text labels
label_tracks_df <- standardized_tracks %>%
  distinct(Lineage_Label) %>%
  mutate(y_coord = as.numeric(Lineage_Label))

# ==============================================================================
# STEP 3: RENDERING THE ALIGNED MULTI-TRACK CANVAS
# ==============================================================================
print("Step 3: Compiling flush-aligned synteny graph panel...")

synteny_plot <- ggplot() +
  geom_line(data = standardized_tracks, 
            aes(x = Relative_Mid_kb, y = as.numeric(Lineage_Label), group = Gene_Family),
            color = "grey60", linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
  
  geom_segment(data = standardized_tracks %>% group_by(Lineage_Label) %>% 
                 summarise(Min_X = min(c(Relative_Start_kb, Relative_End_kb)), 
                           Max_X = max(c(Relative_Start_kb, Relative_End_kb)), .groups="drop"),
               aes(x = Min_X - 15, xend = Max_X + 15, y = as.numeric(Lineage_Label), yend = as.numeric(Lineage_Label)),
               color = "grey90", linewidth = 1.2) +
  
  geom_text(data = break_markers, aes(x = Break_X, y = as.numeric(Lineage_Label), label = "//"),
            color = "grey50", size = 4.0, fontface = "bold") +
  
  # FIX 1: Plot the text labels inside the chart canvas right next to the starting ticks.
  # Placed at -12 kb, right-justified (hjust = 1) to force them perfectly tight against the zero mark.
  geom_text(data = label_tracks_df, aes(x = -12, y = y_coord, label = Lineage_Label),
            hjust = 1, face = "bold", size = 3.4, fontface = "bold", color = "black") +
  
  geom_segment(data = standardized_tracks,
               aes(x = Relative_Start_kb, xend = Relative_End_kb, 
                   y = as.numeric(Lineage_Label), yend = as.numeric(Lineage_Label), 
                   color = Gene_Family),
               linewidth = 7.0, lineend = "butt") +
  
  # Force a slightly wider left horizontal boundary limit (limits = c(-190, 800)) 
  # to leave a clean, dedicated column track zone for the text labels without overlap clipping
  scale_x_continuous(labels = function(x) paste0(x, " kb"), breaks = seq(0, 800, by = 100), limits = c(-190, 800)) +
  scale_y_continuous(breaks = 1:6) +
  scale_fill_manual(values = gene_palette) +
  scale_color_manual(values = gene_palette) +
  labs(
    x = "Relative Genomic Scale",
    y = NULL,
    title = "3p21.3 TSG Cluster Multi-Species Synteny Layout (Normalized Physical Scale)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11, margin = margin(b = 15)),
    axis.text.x = element_text(color = "black"),
    # FIX 2: Completely strip out native text and line elements to wipe out the old whitespace gap
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 9),
    legend.box.margin = margin(t = 10)
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 5)))

ggsave(OUTPUT_PNG, plot = synteny_plot, width = 12.0, height = 5.5, dpi = 400, device = "png")

print("==============================================================================")
print("SUCCESS! Linear track map compiled. Class names are positioned directly flush.")
print("==============================================================================")

#-----------------------------------------------------------------------------
#HEATMAP FOR TOP-5 3p21.3 SUBCLUSTERS (figure)
#-----------------------------------------------------------------------------
# Install required visualization packages if missing
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

library(dplyr)
library(tidyr)
library(ggplot2)

# --- FILE PATHWAY CONFIGURATIONS ---
INPUT_TSV  = "C:\\Users\\...\\Table_3p21_Gene_Shuffling_Frequencies.txt"
OUTPUT_PNG = "C:\\Users\\...\\Figure_3p21_Top5_Frequency_Heatmap.png"

print("Step 1: Reading data matrix from your frequency table...")
if (!file.exists(INPUT_TSV)) stop("Shuffling frequency table missing. Run your step 2 aggregation first.")

raw_frequencies <- read.table(INPUT_TSV, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

print("Step 2: Processing and normalizing Top-5 architectures...")

target_classes <- c("Actinopterygii", "Amphibia", "Reptilia", "Aves", "Mammalia")
col_totals     <- colSums(raw_frequencies[, target_classes], na.rm = TRUE)

# Isolate ONLY the top 5 dominant architectures based on your Total_Species column
top_5_raw <- raw_frequencies %>%
  arrange(desc(Total_Species)) %>%
  slice(1:5)

# Convert counts into relative percentages per class column
heatmap_data <- top_5_raw %>%
  dplyr::select(Final_Architecture, all_of(target_classes)) %>%
  pivot_longer(cols = all_of(target_classes), names_to = "Class", values_to = "Raw_Count") %>%
  mutate(
    Class_Total = case_when(
      Class == "Actinopterygii" ~ col_totals["Actinopterygii"],
      Class == "Amphibia"       ~ col_totals["Amphibia"],
      Class == "Reptilia"       ~ col_totals["Reptilia"],
      Class == "Aves"           ~ col_totals["Aves"],
      Class == "Mammalia"       ~ col_totals["Mammalia"]
    ),
    # Calculate percentage out of the class's total population pool
    Relative_Frequency_Pct = (Raw_Count / Class_Total) * 100,
    
    # Enforce priority sequence: Fish -> Amphibia -> Reptilia -> Aves -> Mammalia
    Class = factor(Class, levels = c("Actinopterygii", "Amphibia", "Reptilia", "Aves", "Mammalia")),
    
    # Substitute clear formatting arrows for clean axis text line views
    Display_Architecture = gsub(" -> ", " → ", Final_Architecture)
  )

# Enforce a strict vertical order based on total global abundance 
arch_order <- unique(heatmap_data$Display_Architecture)
heatmap_data$Display_Architecture <- factor(heatmap_data$Display_Architecture, levels = rev(arch_order))

# ==============================================================================
# STEP 3: RENDERING THE VERTICAL LABEL BLUE-TO-RED MATRIX
# ==============================================================================
print("Step 3: Compiling clean vertical-axis blue-to-red heatmap...")

frequency_heatmap <- ggplot(heatmap_data, aes(x = Class, y = Display_Architecture, fill = Relative_Frequency_Pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  
  # FIX: Flipped the gradient. Blue (#08519c) is lowest (0%), Red (#ef3b2c) is highest (100%)
  scale_fill_gradient(
    low = "#08519c", 
    high = "#ef3b2c", 
    labels = function(x) paste0(x, "%"),
    name = "Relative\nFrequency\nIn Class"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Top 5 Syntenic 3p21.3 Architectures Across Vertebrata",
    subtitle = "Normalized by total column counts within each independent vertebrate lineage"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, color = "grey40", face = "italic", margin = margin(b = 15)),
    axis.text.x = element_text(face = "bold", color = "black", size = 11, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(face = "bold", color = "black", size = 8.5),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8)
  )

# Export the high-resolution file asset to disk
ggsave(OUTPUT_PNG, plot = frequency_heatmap, width = 9.0, height = 5.2, dpi = 400, device = "png")

print("==============================================================================")
print("SUCCESS! Clean blue-to-red frequency heatmap generated successfully:")
print(OUTPUT_PNG)
print("==============================================================================")


#-----------------------------------------------------------------------------
#HEATMAP FOR TOP-5 3p21.3 SUBCLUSTER ARCHITECTURES WITH >=3 GENES (figure)
#-----------------------------------------------------------------------------
# Install required visualization packages if missing
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- FILE PATHWAY CONFIGURATIONS ---
INPUT_XLSX <- "C:\\Users\\...\\3p21.3_cluster_mapping_vert.xlsx"
OUTPUT_PNG <- "C:\\Users\\...\\Figure_3p21_TUSC2_Architecture_MacroAverage_Barplot.png"

print("Step 1: Loading data from your vertical mapping workbook...")
if (!file.exists(INPUT_XLSX)) stop("Master Excel workbook missing. Verify folder path.")

master_data <- readxl::read_excel(INPUT_XLSX)

# Harmonize architecture columns
if ("updated_clean_architecture" %in% colnames(master_data)) {
  master_data$Final_Architecture <- master_data$updated_clean_architecture
} else {
  master_data$Final_Architecture <- master_data$clean_architecture
}

# Clean taxonomy headers
shuffling_df <- master_data %>%
  mutate(
    taxonomic_type = trimws(as.character(taxonomic_type)),
    Class = case_when(
      grepl("Actinopterygii", taxonomic_type, ignore.case = TRUE) ~ "Actinopterygii",
      grepl("Amphibia", taxonomic_type, ignore.case = TRUE)       ~ "Amphibia",
      grepl("Reptilia", taxonomic_type, ignore.case = TRUE)       ~ "Reptilia",
      grepl("Aves", taxonomic_type, ignore.case = TRUE)           ~ "Aves",
      grepl("Mammalia", taxonomic_type, ignore.case = TRUE)       ~ "Mammalia",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Class %in% c("Actinopterygii", "Amphibia", "Reptilia", "Aves", "Mammalia"),
         !is.na(Final_Architecture), Final_Architecture != "", Final_Architecture != "Collapsed")

# Calculate column denominators for independent class normalization totals
class_totals <- shuffling_df %>%
  group_by(Class) %>%
  tally(name = "Class_Total")

print("Step 2: Mining sliding window substrings (Length >= 3) containing TUSC2 exclusively...")
sub_block_list <- list()

for (i in 1:nrow(shuffling_df)) {
  current_class <- shuffling_df$Class[i]
  genes <- unlist(strsplit(shuffling_df$Final_Architecture[i], " -> "))
  if (length(genes) < 3) next # Enforce constraint strictly filtering out subclusters with 2 or fewer items
  
  # Generate all continuous internal combinations (sub-clusters) with a minimum length of 3 genes
  for (start_idx in 1:(length(genes)-2)) {
    for (end_idx in (start_idx+2):length(genes)) {
      sub_seq <- genes[start_idx:end_idx]
      
      # FIX: Restored valid comparison operator from %in= to standard %in%
      if (!("TUSC2" %in% sub_seq)) next
      
      sub_str <- paste(sub_seq, collapse = " → ")
      
      sub_block_list[[length(sub_block_list) + 1]] <- data.frame(
        Class       = current_class,
        Sub_Cluster = sub_str,
        stringsAsFactors = FALSE
      )
    }
  }
}

all_sub_clusters <- bind_rows(sub_block_list)

print("Step 3: Calculating column-wise frequencies and computing Macro-Averages...")
counts_per_class <- all_sub_clusters %>%
  group_by(Sub_Cluster, Class) %>%
  tally(name = "Sub_Count") %>%
  ungroup()

unique_tusc2_subclusters <- unique(all_sub_clusters$Sub_Cluster)

# Calculate the balanced scores to isolate the Top 5 overall configurations across all 5 classes
top_architectures_summary <- counts_per_class %>%
  complete(Sub_Cluster = unique_tusc2_subclusters, Class = c("Actinopterygii", "Amphibia", "Reptilia", "Aves", "Mammalia"), fill = list(Sub_Count = 0)) %>%
  left_join(class_totals, by = "Class") %>%
  mutate(Class_Percentage = (Sub_Count / Class_Total) * 100) %>%
  group_by(Sub_Cluster) %>%
  summarise(Avg_of_Averages_Pct = mean(Class_Percentage, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Avg_of_Averages_Pct)) %>%
  # Isolate top dominant configurations for plot representation
  slice(1:5) %>%
  mutate(
    # Invert factor levels to ensure the most dominant architecture plots at the top position
    Sub_Cluster = factor(Sub_Cluster, levels = rev(Sub_Cluster))
  )

# ==============================================================================
# STEP 4: RENDERING HIGH-RESOLUTION VERTICAL BARPLOT (ARCHITECTURES ON Y)
# ==============================================================================
print("Step 4: Compiling horizontal-bar architecture chart...")

barplot_panel <- ggplot(top_architectures_summary, aes(y = Sub_Cluster, x = Avg_of_Averages_Pct)) +
  # Draw crisp horizontal bars with a clean layout color fill profile
  geom_col(fill = "#3182bd", color = "grey20", linewidth = 0.5, width = 0.55, alpha = 0.85) +
  
  # Format percentages cleanly on the horizontal axis
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0), limits = c(0, max(top_architectures_summary$Avg_of_Averages_Pct) * 1.15)) +
  labs(
    x = "Macro-Average Frequency (Average of Class Percentages)",
    y = NULL,
    title = "Top 5 Conserved TUSC2-Core Subclusters Across Vertebrata",
    subtitle = "Length \u2265 3 containing TUSC2, ranked by unweighted macro-average index"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, color = "grey40", face = "italic", margin = margin(b = 15)),
    # Architecture text alignments read cleanly horizontally on the left border
    axis.text.y = element_text(face = "bold", color = "black", size = 9),
    axis.text.x = element_text(color = "black", size = 10),
    axis.title.x = element_text(face = "bold", size = 10, margin = margin(t = 12)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Export the high-resolution horizontal asset directly to disk
ggsave(OUTPUT_PNG, plot = barplot_panel, width = 9.2, height = 4.2, dpi = 400, device = "png")

print("==============================================================================")
print("SUCCESS! Standardized vertical-axis bar chart written safely to disk:")
print(OUTPUT_PNG)
print("==============================================================================")
print(top_architectures_summary)

#------------------------------------------------
#Chi-Square Test of Independence
#------------------------------------------------
# Ensure required data manipulation libraries are loaded
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")

library(readxl)
library(dplyr)
library(tidyr)

# --- FILE PATHWAY CONFIGURATIONS ---
INPUT_XLSX <- "C:\\Users\\...\\3p21.3_cluster_mapping_vert.xlsx"

print("Step 1: Loading comprehensive excel architecture tracks...")
if (!file.exists(INPUT_XLSX)) stop("Master Excel workbook missing. Verify folder path.")

master_data <- readxl::read_excel(INPUT_XLSX)

# Harmonize architecture columns natively based on your available fields
if ("updated_clean_architecture" %in% colnames(master_data)) {
  master_data$Final_Architecture <- master_data$updated_clean_architecture
} else {
  master_data$Final_Architecture <- master_data$clean_architecture
}

# Clean taxonomic group whitespace properties
shuffling_df <- master_data %>%
  mutate(
    taxonomic_type = trimws(as.character(taxonomic_type)),
    Class = case_when(
      grepl("Actinopterygii", taxonomic_type, ignore.case = TRUE) ~ "Actinopterygii",
      grepl("Amphibia", taxonomic_type, ignore.case = TRUE)       ~ "Amphibia",
      grepl("Reptilia", taxonomic_type, ignore.case = TRUE)       ~ "Reptilia",
      grepl("Aves", taxonomic_type, ignore.case = TRUE)           ~ "Aves",
      grepl("Mammalia", taxonomic_type, ignore.case = TRUE)       ~ "Mammalia",
      TRUE ~ "Other"
    )
  ) %>%
  # Filter strictly down to your 5 core publication vertebrate lineages
  filter(Class %in% c("Actinopterygii", "Amphibia", "Reptilia", "Aves", "Mammalia"),
         !is.na(Final_Architecture), Final_Architecture != "", Final_Architecture != "Collapsed")

# ==============================================================================
# STEP 2: CONSTURCT CON TINGENCY CROSS-TABULATION MATRIX
# ==============================================================================
print("Step 2: Cross-tabulating architectures by vertebrate lineage columns...")

contingency_table <- shuffling_df %>%
  # Ungroup completely to wipe out internal grouping metadata locks
  ungroup() %>%
  group_by(Final_Architecture, Class) %>%
  tally() %>%
  ungroup() %>%
  # Pivot into a clean wide matrix layout filled with structural zeros
  pivot_wider(names_from = Class, values_from = n, values_fill = 0)

# Extract ONLY the numeric count data streams, leaving text labels behind
numeric_matrix <- contingency_table %>%
  dplyr::select(Actinopterygii, Amphibia, Reptilia, Aves, Mammalia) %>%
  as.matrix()

# Map architectures back as matrix row markers for logging integrity
rownames(numeric_matrix) <- contingency_table$Final_Architecture

# ==============================================================================
# STEP 3: EXECUTE CHI-SQUARE INDEPENDENCE TEST WITH SIMULATED BASES
# ==============================================================================
print("Step 3: Calculating statistical independence parameters...")

# Run the native test using 2,000 Monte Carlo permutation replicates to handle zeroes safely
r_chi2 <- chisq.test(numeric_matrix, simulate.p.value = TRUE, B = 2000)

print("==============================================================================")
print("FINAL MANUSCRIPT STATISTICS: CHI-SQUARE TEST OF INDEPENDENCE")
print("==============================================================================")
print(r_chi2)
print("==============================================================================")

# Display your raw contingency matrix straight to the screen to inspect cell spreads
print("--- RAW OBSERVED CONTINGENCY MATRIX CODES ---")
print(numeric_matrix)

# Ensure required packages are loaded
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")

library(readxl)
library(dplyr)
library(tidyr)

# --- FILE PATHWAY CONFIGURATIONS ---
INPUT_XLSX <- "C:\\Users\\...\\3p21.3_cluster_mapping_vert.xlsx"

print("Step 1: Loading excel matrix...")
if (!file.exists(INPUT_XLSX)) stop("Master Excel workbook missing.")

master_data <- readxl::read_excel(INPUT_XLSX)

if ("updated_clean_architecture" %in% colnames(master_data)) {
  master_data$Final_Architecture <- master_data$updated_clean_architecture
} else {
  master_data$Final_Architecture <- master_data$clean_architecture
}

# Standardize lineage data lanes
shuffling_df <- master_data %>%
  mutate(
    taxonomic_type = trimws(as.character(taxonomic_type)),
    Class = case_when(
      grepl("Actinopterygii", taxonomic_type, ignore.case = TRUE) ~ "Actinopterygii",
      grepl("Amphibia", taxonomic_type, ignore.case = TRUE)       ~ "Amphibia",
      grepl("Reptilia", taxonomic_type, ignore.case = TRUE)       ~ "Reptilia",
      grepl("Aves", taxonomic_type, ignore.case = TRUE)           ~ "Aves",
      grepl("Mammalia", taxonomic_type, ignore.case = TRUE)       ~ "Mammalia",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Class %in% c("Actinopterygii", "Amphibia", "Reptilia", "Aves", "Mammalia"),
         !is.na(Final_Architecture), Final_Architecture != "", Final_Architecture != "Collapsed")

# ==============================================================================
# STEP 2: COMPILING OBSERVED CONTINGENCY MATRIX & TRUE CHI-SQ
# ==============================================================================
print("Step 2: Building contingency grid matrix...")

contingency_table <- shuffling_df %>%
  ungroup() %>%
  group_by(Final_Architecture, Class) %>%
  tally() %>%
  ungroup() %>%
  pivot_wider(names_from = Class, values_from = n, values_fill = 0)

numeric_matrix <- contingency_table %>%
  dplyr::select(Actinopterygii, Amphibia, Reptilia, Aves, Mammalia) %>%
  as.matrix()

# Calculate your true, real-world observed Chi-Square score baseline
observed_chisq <- as.numeric(chisq.test(numeric_matrix, correct = FALSE)$statistic)
print(paste("--> Real-World Observed Chi-Square Statistic Baseline:", round(observed_chisq, 2)))

# ==============================================================================
# STEP 3: HIGH-ITERATION MONTE CARLO PERMUTATION LOOP FOR P-VALUE
# ==============================================================================
print("Step 3: Simulating 5,000 randomized distributions to extract the empirical p-value...")

# Define the number of background permutations (higher B gives higher decimal precision)
B_PERMUTATIONS <- 5000
random_chisq_scores <- numeric(B_PERMUTATIONS)

set.seed(42) # Set a random seed to make the simulation perfectly reproducible

# Extract raw vector arrays for rapid scrambling loops
raw_arch  <- shuffling_df$Final_Architecture
raw_class <- shuffling_df$Class

for (b in 1:B_PERMUTATIONS) {
  # Randomly scramble the taxonomy column while keeping architecture fixed (Null Hypothesis)
  scrambled_class <- sample(raw_class)
  
  # Re-tabulate the randomized data frame in memory
  sim_matrix <- table(raw_arch, scrambled_class)
  
  # Run a fast Chi-Square metric evaluation on the scrambled grid
  random_chisq_scores[b] <- suppressWarnings(chisq.test(sim_matrix, correct = FALSE)$statistic)
}

# Calculate the exact empirical p-value: (Number of times random score >= observed score) / Total simulations
times_exceeded <- sum(random_chisq_scores >= observed_chisq)
empirical_p_value <- times_exceeded / B_PERMUTATIONS

# ==============================================================================
# STEP 4: PRESENTING FINAL MANUSCRIPT STATS
# ==============================================================================
print("==============================================================================")
print("DEDICATED P-VALUE EMPIRICAL PERMUTATION SUMMARY RESULTS")
print("==============================================================================")
print(paste("Observed Chi-Square Statistic  :", round(observed_chisq, 4)))
print(paste("Total Permutation Replicates    :", B_PERMUTATIONS))
print(paste("Times Scrambled Score Exceeded Real:", times_exceeded))

# Format the final p-value cleanly for your journal paper submission
if (empirical_p_value == 0) {
  print(paste("Empirical Significance P-Value  : p <", format(1 / B_PERMUTATIONS, scientific = FALSE)))
} else {
  print(paste("Empirical Significance P-Value  : p =", format(empirical_p_value, scientific = TRUE, digits = 4)))
}
print("==============================================================================")


#----------------------------------------------------------------------------
#RANDOMNESS TEST FOR CORE 3p21.3 SUBCLUSTERS
#----------------------------------------------------------------------------
# Install required visualization packages if missing
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")

library(readxl)
library(dplyr)
library(tidyr)

# --- FILE PATHWAY CONFIGURATIONS ---
INPUT_XLSX <- "C:\\Users\\...\\3p21.3_cluster_mapping_vert.xlsx"

print("Step 1: Loading data from your vertical mapping workbook...")
if (!file.exists(INPUT_XLSX)) stop("Master Excel workbook missing. Verify folder path.")

master_data <- readxl::read_excel(INPUT_XLSX)

# Harmonize architecture columns
if ("updated_clean_architecture" %in% colnames(master_data)) {
  master_data$Final_Architecture <- master_data$updated_clean_architecture
} else {
  master_data$Final_Architecture <- master_data$clean_architecture
}

# Filter out collapsed rows
valid_data <- master_data %>%
  filter(!is.na(Final_Architecture), Final_Architecture != "", Final_Architecture != "Collapsed")

# ==============================================================================
# STEP 2: MINING REAL-WORLD TARGET SUBCLUSTERS (OBSERVED BASELINE)
# ==============================================================================
print("Step 2: Extracting real-world TUSC2 subclusters...")
observed_list <- list()

for (i in 1:nrow(valid_data)) {
  genes <- unlist(strsplit(valid_data$Final_Architecture[i], " -> "))
  if (length(genes) < 3) next
  
  for (start_idx in 1:(length(genes)-2)) {
    for (end_idx in (start_idx+2):length(genes)) {
      sub_seq <- genes[start_idx:end_idx]
      if ("TUSC2" %in% sub_seq) {
        observed_list[[length(observed_list) + 1]] <- paste(sub_seq, collapse = " -> ")
      }
    }
  }
}

# Isolate Top 3 dominant subclusters to test
top_subclusters_to_test <- data.frame(Sub_Cluster = unlist(observed_list), stringsAsFactors=FALSE) %>%
  group_by(Sub_Cluster) %>%
  tally(name = "Observed_Count") %>%
  arrange(desc(Observed_Count)) %>%
  slice(1:3)

print("Target Subclusters to test:")
print(top_subclusters_to_test)

# ==============================================================================
# STEP 3: HIGH-SPEED VECTOR PERMUTATION WITH PROGRESS BAR
# ==============================================================================
print("Step 3: Launching high-speed vector scrambling permutation test...")

B_REPLICATES <- 2000
accidental_hits_matrix <- matrix(0, nrow = nrow(top_subclusters_to_test), ncol = B_REPLICATES)
rownames(accidental_hits_matrix) <- top_subclusters_to_test$Sub_Cluster

# Pre-split rows and targets into structural vector lists for massive performance boost
parsed_architectures <- lapply(valid_data$Final_Architecture, function(x) unlist(strsplit(x, " -> ")))
target_vectors       <- lapply(top_subclusters_to_test$Sub_Cluster, function(x) unlist(strsplit(x, " -> ")))

# Vector sequence matching function (Runs blazingly fast)
is_subvector <- function(sub, main) {
  n <- length(main)
  m <- length(sub)
  if (m > n) return(FALSE)
  for (i in 1:(n - m + 1)) {
    if (all(main[i:(i + m - 1)] == sub)) return(TRUE)
  }
  return(FALSE)
}

set.seed(123) # Re-lock random seeds

# INITIALIZE REAL-TIME PROGRESS BAR
pb <- txtProgressBar(min = 0, max = B_REPLICATES, style = 3)

for (b in 1:B_REPLICATES) {
  round_counts <- numeric(nrow(top_subclusters_to_test))
  names(round_counts) <- top_subclusters_to_test$Sub_Cluster
  
  for (i in seq_along(parsed_architectures)) {
    genes <- parsed_architectures[[i]]
    if (length(genes) < 3) next
    
    # Scramble the vector directly
    shuffled_genes <- sample(genes)
    
    # Check for presence using fast vector slicing instead of text-parsing regex strings
    for (t_idx in seq_along(target_vectors)) {
      if (is_subvector(target_vectors[[t_idx]], shuffled_genes)) {
        sub_name <- top_subclusters_to_test$Sub_Cluster[t_idx]
        round_counts[sub_name] <- round_counts[sub_name] + 1
      }
    }
  }
  
  accidental_hits_matrix[, b] <- round_counts
  
  # UPDATE THE PROGRESSION BAR VISUAL
  setTxtProgressBar(pb, b)
}
close(pb) # Safely close progress bar handle

# ==============================================================================
# STEP 4: PRESENT COMPACTED P-VALUE SIGNIFICANCE SUMMARY
# ==============================================================================
print("==============================================================================")
print("CORE SUBCLUSTER RANDOMNESS PERMUTATION TEST RESULTS")
print("==============================================================================")

results_summary <- top_subclusters_to_test %>%
  mutate(
    Max_Accidental_In_A_Single_Shuffle = apply(accidental_hits_matrix, 1, max),
    Times_Random_Beat_Observed         = sapply(Sub_Cluster, function(x) sum(accidental_hits_matrix[x, ] >= Observed_Count)),
    Empirical_P_Value                  = Times_Random_Beat_Observed / B_REPLICATES,
    Significance                       = ifelse(Empirical_P_Value < 0.05, "STATISTICALLY SIGNIFICANT (NON-RANDOM)", "NEUTRAL REARRANGEMENT"),
    Display_Sub_Cluster                = gsub(" -> ", " → ", Sub_Cluster)
  ) %>%
  dplyr::select(Display_Sub_Cluster, Observed_Count, Max_Accidental_In_A_Single_Shuffle, Empirical_P_Value, Significance)

print(results_summary)
print("==============================================================================")


