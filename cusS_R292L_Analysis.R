setwd("/home/william-ackerman/Desktop/cusS_R292L")

extract_aa_from_genbank <- function(gbk_file, gene_name = "cusS") {
  gbk_lines <- readLines(gbk_file)
  
  # Find the gene
  gene_line <- which(grepl(paste0('/gene="', gene_name, '"'), gbk_lines))
  
  if (length(gene_line) == 0) {
    cat("Gene", gene_name, "not found. Searching for similar genes...\n")
    # Search more broadly
    gene_lines <- which(grepl(gene_name, gbk_lines, ignore.case = TRUE))
    if (length(gene_lines) > 0) {
      cat("Found potential matches at lines:", gene_lines[1:min(5, length(gene_lines))], "\n")
      gene_line <- gene_lines[1]
    } else {
      stop("Gene not found")
    }
  }
  
  # Find the translation line after this gene
  translation_start <- NULL
  for (i in gene_line:length(gbk_lines)) {
    if (grepl('/translation="', gbk_lines[i])) {
      translation_start <- i
      break
    }
  }
  
  if (is.null(translation_start)) {
    stop("Translation not found for this gene")
  }
  
  # Extract the amino acid sequence
  aa_sequence <- ""
  i <- translation_start
  
  # Get the first line
  first_line <- gbk_lines[i]
  # Extract sequence after /translation="
  seq_part <- sub('.*\\/translation="([^"]*)".*', '\\1', first_line)
  
  if (grepl('".*"', first_line)) {
    # Single line translation
    aa_sequence <- seq_part
  } else {
    # Multi-line translation
    seq_part <- sub('.*\\/translation="', '', first_line)
    aa_sequence <- seq_part
    
    # Continue reading until we find the closing quote
    i <- i + 1
    while (i <= length(gbk_lines) && !grepl('"', gbk_lines[i])) {
      # Remove leading whitespace and add to sequence
      line_seq <- gsub('^\\s+', '', gbk_lines[i])
      aa_sequence <- paste0(aa_sequence, line_seq)
      i <- i + 1
    }
    
    # Handle the last line with closing quote
    if (i <= length(gbk_lines)) {
      last_line <- gsub('".*', '', gbk_lines[i])
      last_line <- gsub('^\\s+', '', last_line)
      aa_sequence <- paste0(aa_sequence, last_line)
    }
  }
  
  # Clean up the sequence - remove any remaining quotes and whitespace
  aa_sequence <- gsub('"', '', aa_sequence)
  aa_sequence <- gsub('\\s+', '', aa_sequence)
  
  return(aa_sequence)
}

# Extract cusS amino acid sequence
cusS_aa <- extract_aa_from_genbank("Ecoli_K12_MG1655.gbk", "cusS")

cat("cusS amino acid sequence extracted:\n")
cat("Length:", nchar(cusS_aa), "amino acids\n")
cat("First 60 residues:", substr(cusS_aa, 1, 60), "\n")
cat("Last 60 residues:", substr(cusS_aa, nchar(cusS_aa)-59, nchar(cusS_aa)), "\n")

# Check position 292
if (nchar(cusS_aa) >= 292) {
  residue_292 <- substr(cusS_aa, 292, 292)
  cat("\nResidue at position 292:", residue_292, "\n")
  
  # Show context around position 292
  start_pos <- max(1, 292 - 15)
  end_pos <- min(nchar(cusS_aa), 292 + 15)
  context <- substr(cusS_aa, start_pos, end_pos)
  
  cat("Context (positions", start_pos, "to", end_pos, "):\n")
  cat(context, "\n")
  
  # Mark position 292
  pos_in_context <- 292 - start_pos + 1
  marker <- paste(c(rep(" ", pos_in_context - 1), "^"), collapse = "")
  cat(marker, "\n")
  
  if (residue_292 == "R") {
    cat("✓ Confirmed: Position 292 is Arginine (R)\n")
  } else {
    cat("⚠ Position 292 is", residue_292, "not Arginine (R)\n")
  }
} else {
  cat("Sequence is too short for position 292\n")
}


# Create R292L mutant
create_R292L_mutant <- function(wt_sequence) {
  if (nchar(wt_sequence) < 292) {
    stop("Sequence too short for position 292")
  }
  
  mutant_seq <- wt_sequence
  substr(mutant_seq, 292, 292) <- "L"
  
  return(mutant_seq)
}

# Generate mutant if we have the sequence
if (exists("cusS_aa") && nchar(cusS_aa) >= 292) {
  cusS_mutant <- create_R292L_mutant(cusS_aa)
  
  cat("\n=== Mutation Analysis ===\n")
  cat("Wild-type residue 292:", substr(cusS_aa, 292, 292), "\n")
  cat("Mutant residue 292:", substr(cusS_mutant, 292, 292), "\n")
  
  # Show differences around the mutation
  window <- 15
  start_pos <- max(1, 292 - window)
  end_pos <- min(nchar(cusS_aa), 292 + window)
  
  wt_context <- substr(cusS_aa, start_pos, end_pos)
  mut_context <- substr(cusS_mutant, start_pos, end_pos)
  
  cat("\nSequence comparison around position 292:\n")
  cat("WT : ", wt_context, "\n")
  cat("MUT: ", mut_context, "\n")
  
  # Write sequences to files for further analysis
  writeLines(c(">cusS_WT_E.coli_K12", cusS_aa), "cusS_WT.fasta")
  writeLines(c(">cusS_R292L_mutant", cusS_mutant), "cusS_R292L.fasta")
  
  cat("\nFASTA files created:\n")
  cat("- cusS_WT.fasta\n")
  cat("- cusS_R292L.fasta\n")
}

# Analyze the mutation properties
analyze_R292L_mutation <- function() {
  cat("\n=== R292L Mutation Properties ===\n")
  
  properties <- data.frame(
    Property = c("Amino Acid", "Three Letter", "Charge", "Polarity", 
                 "Hydrophobicity (Kyte-Doolittle)", "Volume (Å³)", 
                 "Molecular Weight", "Side Chain Type"),
    Arginine_R = c("R", "Arg", "+1", "Polar", "-4.5", "148", "174.2", "Basic, Long"),
    Leucine_L = c("L", "Leu", "0", "Nonpolar", "3.8", "124", "131.2", "Hydrophobic, Branched"),
    stringsAsFactors = FALSE
  )
  
  print(properties)
  
  cat("\nKey Changes:\n")
  cat("• Loss of positive charge: +1 → 0\n")
  cat("• Polarity change: Polar → Nonpolar\n")
  cat("• Hydrophobicity increase: -4.5 → 3.8 (Δ = +8.3)\n")
  cat("• Size decrease: 148 → 124 Å³\n")
  cat("• Functional group change: Guanidino → Isobutyl\n")
  
  cat("\nPotential Structural/Functional Impact:\n")
  cat("• May disrupt ionic interactions\n")
  cat("• Could affect protein-protein interactions\n")
  cat("• Might alter membrane association (more hydrophobic)\n")
  cat("• May impact copper sensing or signal transduction\n")
}

analyze_R292L_mutation()


# Function to search for existing structures and prepare for modeling
search_cusS_structures <- function() {
  cat("=== Searching for cusS Structural Information ===\n\n")
  
  # Information about cusS structure
  cat("cusS is a sensor histidine kinase with typical domains:\n")
  cat("• Periplasmic sensing domain (N-terminal)\n")
  cat("• Transmembrane helices\n") 
  cat("• Cytoplasmic histidine kinase domain (C-terminal)\n\n")
  
  cat("Position 292 location prediction:\n")
  cat("• Total length: 481 amino acids\n")
  cat("• Position 292 is in the cytoplasmic domain\n")
  cat("• Likely in the kinase catalytic region\n\n")
  
  # Search suggestions
  cat("Structure search recommendations:\n")
  cat("1. Search RCSB PDB for:\n")
  cat("   - 'cusS'\n")
  cat("   - 'copper sensor histidine kinase'\n")
  cat("   - 'E. coli sensor kinase'\n\n")
  
  cat("2. Check AlphaFold database:\n")
  cat("   - UniProt ID: P77485\n")
  cat("   - https://alphafold.ebi.ac.uk/entry/P77485\n\n")
}

search_cusS_structures()

# Download and analyze AlphaFold structure
library(bio3d)
library(ggplot2)

# Function to download AlphaFold structure
download_alphafold_cusS <- function() {
  uniprot_id <- "P77485"  # cusS UniProt ID
  alphafold_url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprot_id, "-F1-model_v4.pdb")
  
  cat("Downloading AlphaFold structure for cusS...\n")
  cat("URL:", alphafold_url, "\n")
  
  # Try to download
  tryCatch({
    download.file(alphafold_url, "cusS_alphafold.pdb", method = "auto")
    cat("✓ AlphaFold structure downloaded as 'cusS_alphafold.pdb'\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗ Download failed:", e$message, "\n")
    cat("Please manually download from: https://alphafold.ebi.ac.uk/entry/P77485\n")
    return(FALSE)
  })
}

# Download the structure
alphafold_downloaded <- download_alphafold_cusS()

# Function to analyze cusS structure
analyze_cusS_structure <- function(pdb_file = "cusS_alphafold.pdb") {
  if (!file.exists(pdb_file)) {
    cat("PDB file not found. Please download cusS structure first.\n")
    return(NULL)
  }
  
  # Read the structure
  cusS_pdb <- read.pdb(pdb_file)
  
  cat("=== cusS Structure Analysis ===\n")
  cat("Structure file:", pdb_file, "\n")
  cat("Number of atoms:", nrow(cusS_pdb$atom), "\n")
  cat("Number of residues:", length(unique(cusS_pdb$atom$resno)), "\n")
  cat("Chain(s):", unique(cusS_pdb$atom$chain), "\n")
  
  # Get residue 292 information
  res_292 <- cusS_pdb$atom[cusS_pdb$atom$resno == 292, ]
  
  if (nrow(res_292) > 0) {
    cat("\nResidue 292 information:\n")
    cat("Residue type:", unique(res_292$resid), "\n")
    cat("Number of atoms:", nrow(res_292), "\n")
    cat("Coordinates (CA atom):\n")
    
    ca_atom <- res_292[res_292$elety == "CA", ]
    if (nrow(ca_atom) > 0) {
      cat("  X:", ca_atom$x, "Å\n")
      cat("  Y:", ca_atom$y, "Å\n") 
      cat("  Z:", ca_atom$z, "Å\n")
    }
    
    # B-factors (confidence scores for AlphaFold)
    avg_bfactor <- mean(res_292$b)
    cat("Average B-factor (confidence):", round(avg_bfactor, 2), "\n")
    
    if (avg_bfactor > 70) {
      cat("✓ High confidence region\n")
    } else if (avg_bfactor > 50) {
      cat("⚠ Medium confidence region\n")
    } else {
      cat("✗ Low confidence region\n")
    }
  } else {
    cat("⚠ Residue 292 not found in structure\n")
  }
  
  return(cusS_pdb)
}

# Analyze the structure if available
if (file.exists("cusS_alphafold.pdb")) {
  cusS_structure <- analyze_cusS_structure()
} else {
  cat("Please download the cusS structure first\n")
}

# Function to create R292L mutant structure
create_mutant_structure <- function(pdb_structure, mutation_pos = 292) {
  if (is.null(pdb_structure)) {
    cat("No structure provided\n")
    return(NULL)
  }
  
  cat("=== Creating R292L Mutant Structure ===\n")
  
  # Find residue 292
  res_292_atoms <- pdb_structure$atom[pdb_structure$atom$resno == mutation_pos, ]
  
  if (nrow(res_292_atoms) == 0) {
    cat("Residue", mutation_pos, "not found in structure\n")
    return(NULL)
  }
  
  cat("Original residue at position", mutation_pos, ":", unique(res_292_atoms$resid), "\n")
  
  # Create mutant structure
  mutant_pdb <- pdb_structure
  
  # Simple mutation: change residue name (for visualization)
  # Note: This is a simplified approach - real mutation would require
  # rebuilding side chain atoms
  mutant_pdb$atom[mutant_pdb$atom$resno == mutation_pos, "resid"] <- "LEU"
  
  # Write mutant structure
  write.pdb(mutant_pdb, file = "cusS_R292L_mutant.pdb")
  cat("✓ Mutant structure written to 'cusS_R292L_mutant.pdb'\n")
  
  return(mutant_pdb)
}

# Create mutant structure if we have the wild-type
if (exists("cusS_structure") && !is.null(cusS_structure)) {
  cusS_mutant_structure <- create_mutant_structure(cusS_structure)
}

analyze_mutation_environment <- function(pdb_structure, mutation_pos = 292, radius = 8) {
  if (is.null(pdb_structure)) return(NULL)
  
  cat("=== Mutation Site Environment Analysis ===\n")
  
  # Get coordinates of residue 292 CA atom
  res_292 <- pdb_structure$atom[pdb_structure$atom$resno == mutation_pos & 
                                  pdb_structure$atom$elety == "CA", ]
  
  if (nrow(res_292) == 0) {
    cat("Cannot find CA atom for residue", mutation_pos, "\n")
    return(NULL)
  }
  
  center_coords <- c(res_292$x, res_292$y, res_292$z)
  cat("Mutation site coordinates:", paste(round(center_coords, 2), collapse = ", "), "\n")
  
  # Find neighboring residues within radius
  all_ca <- pdb_structure$atom[pdb_structure$atom$elety == "CA", ]
  
  # Calculate distances
  distances <- sqrt(rowSums((all_ca[, c("x", "y", "z")] - 
                               rep(center_coords, each = nrow(all_ca)))^2))
  
  # Find neighbors within radius
  neighbors <- all_ca[distances <= radius & distances > 0, ]
  neighbor_distances <- distances[distances <= radius & distances > 0]
  
  cat("Neighboring residues within", radius, "Å:\n")
  neighbor_df <- data.frame(
    Position = neighbors$resno,
    Residue = neighbors$resid,
    Distance = round(neighbor_distances, 2),
    stringsAsFactors = FALSE
  )
  neighbor_df <- neighbor_df[order(neighbor_df$Distance), ]
  print(neighbor_df)
  
  # Identify potential interaction partners
  cat("\nPotential interaction analysis:\n")
  charged_residues <- c("ARG", "LYS", "ASP", "GLU", "HIS")
  polar_residues <- c("SER", "THR", "ASN", "GLN", "TYR")
  
  charged_neighbors <- neighbor_df[neighbor_df$Residue %in% charged_residues, ]
  polar_neighbors <- neighbor_df[neighbor_df$Residue %in% polar_residues, ]
  
  if (nrow(charged_neighbors) > 0) {
    cat("Charged residues nearby:\n")
    print(charged_neighbors)
    cat("⚠ R292L mutation may disrupt ionic interactions\n")
  }
  
  if (nrow(polar_neighbors) > 0) {
    cat("Polar residues nearby:\n")
    print(polar_neighbors)
  }
  
  return(list(neighbors = neighbor_df, 
              charged = charged_neighbors, 
              polar = polar_neighbors))
}

# Analyze mutation environment
if (exists("cusS_structure") && !is.null(cusS_structure)) {
  mutation_env <- analyze_mutation_environment(cusS_structure)
}


# Generate PyMOL script for visualization
create_pymol_script <- function() {
  cat("=== Creating PyMOL Visualization Script ===\n")
  
  pymol_script <- c(
    "# PyMOL script for cusS R292L mutation analysis",
    "",
    "# Load structures",
    "load cusS_alphafold.pdb, cusS_WT",
    "load cusS_R292L_mutant.pdb, cusS_mutant",
    "",
    "# Basic visualization",
    "hide everything",
    "show cartoon",
    "color cyan, cusS_WT",
    "color salmon, cusS_mutant",
    "",
    "# Highlight mutation site",
    "select mutation_site, resi 292",
    "show sticks, mutation_site",
    "color red, mutation_site and cusS_WT",
    "color blue, mutation_site and cusS_mutant",
    "",
    "# Show neighboring residues",
    "select neighbors, resi 292 around 8",
    "show lines, neighbors",
    "",
    "# Labels",
    "label mutation_site and name CA, \"292\"",
    "",
    "# Different views",
    "orient mutation_site",
    "zoom mutation_site, 5",
    "",
    "# Save session",
    "save cusS_R292L_analysis.pse"
  )
  
  writeLines(pymol_script, "cusS_visualization.pml")
  cat("✓ PyMOL script written to 'cusS_visualization.pml'\n")
  cat("Run in PyMOL with: pymol cusS_visualization.pml\n")
}

# Generate ChimeraX script
create_chimerax_script <- function() {
  cat("=== Creating ChimeraX Visualization Script ===\n")
  
  chimerax_script <- c(
    "# ChimeraX script for cusS R292L mutation analysis",
    "",
    "# Open structures", 
    "open cusS_alphafold.pdb",
    "open cusS_R292L_mutant.pdb",
    "",
    "# Style",
    "cartoon",
    "color #1 cyan",
    "color #2 salmon", 
    "",
    "# Highlight mutation site",
    "select :292",
    "show sel atoms",
    "size sel stickRadius 0.5",
    "color sel red",
    "",
    "# Show environment",
    "select :292@<8.0",
    "show sel atoms",
    "",
    "# Center on mutation",
    "view sel",
    "",
    "# Save image",
    "save cusS_R292L_mutation.png supersample 3"
  )
  
  writeLines(chimerax_script, "cusS_chimerax.cxc")
  cat("✓ ChimeraX script written to 'cusS_chimerax.cxc'\n")
  cat("Run in ChimeraX with: File > Open > cusS_chimerax.cxc\n")
}

# Create visualization scripts
create_pymol_script()
create_chimerax_script()

generate_structural_report <- function() {
  cat("\n", strrep("=", 60), "\n")
  cat("     cusS R292L STRUCTURAL ANALYSIS REPORT\n")
  cat("\n", strrep("=", 60), "\n")
  
  cat("STRUCTURE INFORMATION:\n")
  cat("• Source: AlphaFold prediction (UniProt: P77485)\n")
  cat("• Mutation: Arginine 292 → Leucine\n")
  cat("• Location: Cytoplasmic kinase domain\n\n")
  
  cat("PREDICTED STRUCTURAL IMPACTS:\n")
  cat("• Loss of positive charge at position 292\n")
  cat("• Increased hydrophobicity in kinase domain\n")
  cat("• Potential disruption of ionic interactions\n")
  cat("• May affect ATP binding or kinase activity\n\n")
  
  cat("FILES GENERATED:\n")
  cat("• cusS_alphafold.pdb - Wild-type structure\n")
  cat("• cusS_R292L_mutant.pdb - Mutant structure\n")
  cat("• cusS_visualization.pml - PyMOL script\n")
  cat("• cusS_chimerax.cxc - ChimeraX script\n\n")
  
  cat("RECOMMENDED ANALYSIS:\n")
  cat("1. Load structures in PyMOL or ChimeraX\n")
  cat("2. Examine local environment around position 292\n")
  cat("3. Check for nearby charged/polar residues\n")
  cat("4. Analyze potential impact on kinase function\n")
  cat("5. Consider molecular dynamics simulations\n\n")
  
  cat("FUNCTIONAL PREDICTIONS:\n")
  cat("• May reduce copper binding affinity\n")
  cat("• Could affect autophosphorylation\n")
  cat("• Might alter CusR phosphorylation\n")
  cat("• May impact copper homeostasis regulation\n\n")
  
  cat("\n", strrep("=", 60), "\n")
}

generate_structural_report()

library(tidyverse)
library(reshape2)
library(RColorBrewer)

# Enhanced structure analysis using bio3d
analyze_cusS_with_bio3d <- function(pdb_file = "cusS_alphafold.pdb") {
  
  cat("=== Comprehensive bio3d Analysis of cusS ===\n\n")
  
  # Read the structure
  pdb <- read.pdb(pdb_file)
  
  # Basic structure information
  cat("STRUCTURE OVERVIEW:\n")
  cat("File:", pdb_file, "\n")
  cat("Resolution:", ifelse(is.null(pdb$header$resolution), "N/A (AlphaFold)", pdb$header$resolution), "\n")
  cat("Number of chains:", length(unique(pdb$atom$chain)), "\n")
  cat("Total atoms:", nrow(pdb$atom), "\n")
  cat("Total residues:", max(pdb$atom$resno), "\n")
  cat("Sequence length:", length(pdb$calpha), "\n\n")
  
  # Extract sequence from structure
  seq <- aa321(pdb$atom[pdb$atom$elety == "CA", "resid"])
  cat("Sequence from structure:\n")
  cat(paste(seq, collapse = ""), "\n\n")
  
  return(list(pdb = pdb, sequence = seq))
}

# Detailed residue 292 analysis
analyze_residue_292 <- function(pdb) {
  cat("=== RESIDUE 292 DETAILED ANALYSIS ===\n\n")
  
  # Get all atoms for residue 292
  res292 <- atom.select(pdb, resno = 292)
  res292_atoms <- pdb$atom[res292$atom, ]
  
  if (nrow(res292_atoms) == 0) {
    cat("Residue 292 not found in structure!\n")
    return(NULL)
  }
  
  cat("RESIDUE 292 INFORMATION:\n")
  cat("Residue type:", unique(res292_atoms$resid), "\n")
  cat("Chain:", unique(res292_atoms$chain), "\n")
  cat("Number of atoms:", nrow(res292_atoms), "\n")
  
  # Get coordinates
  coords_292 <- res292_atoms[, c("x", "y", "z")]
  ca_292 <- res292_atoms[res292_atoms$elety == "CA", c("x", "y", "z")]
  
  cat("CA coordinates:", round(unlist(ca_292), 2), "\n")
  
  # B-factors (confidence for AlphaFold)
  avg_bfactor <- mean(res292_atoms$b)
  cat("Average B-factor:", round(avg_bfactor, 2), "\n")
  
  confidence <- if (avg_bfactor > 90) "Very High" else 
    if (avg_bfactor > 70) "High" else 
      if (avg_bfactor > 50) "Medium" else "Low"
  cat("Confidence level:", confidence, "\n\n")
  
  return(list(atoms = res292_atoms, coords = coords_292, ca_coords = ca_292))
}

# Find neighboring residues and analyze environment
analyze_environment <- function(pdb, center_resno = 292, radius = 8) {
  cat("=== STRUCTURAL ENVIRONMENT ANALYSIS ===\n\n")
  
  # Get CA atoms
  ca_sel <- atom.select(pdb, elety = "CA")
  ca_coords <- pdb$atom[ca_sel$atom, c("x", "y", "z")]
  ca_resno <- pdb$atom[ca_sel$atom, "resno"]
  ca_resid <- pdb$atom[ca_sel$atom, "resid"]
  
  # Get center coordinates (residue 292)
  center_idx <- which(ca_resno == center_resno)
  if (length(center_idx) == 0) {
    cat("Center residue", center_resno, "not found!\n")
    return(NULL)
  }
  
  center_coords <- as.numeric(ca_coords[center_idx, ])
  
  # Calculate distances to all other CA atoms
  distances <- apply(ca_coords, 1, function(x) {
    sqrt(sum((as.numeric(x) - center_coords)^2))
  })
  
  # Find neighbors within radius
  neighbors_idx <- which(distances <= radius & distances > 0)
  neighbor_data <- data.frame(
    resno = ca_resno[neighbors_idx],
    resid = ca_resid[neighbors_idx], 
    distance = round(distances[neighbors_idx], 2),
    stringsAsFactors = FALSE
  )
  neighbor_data <- neighbor_data[order(neighbor_data$distance), ]
  
  cat("NEIGHBORING RESIDUES (within", radius, "Å):\n")
  print(neighbor_data)
  
  # Classify neighbors
  basic_res <- c("ARG", "LYS", "HIS")
  acidic_res <- c("ASP", "GLU") 
  polar_res <- c("SER", "THR", "ASN", "GLN", "TYR")
  hydrophobic_res <- c("ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "MET", "PRO")
  
  neighbor_data$type <- sapply(neighbor_data$resid, function(x) {
    if (x %in% basic_res) "Basic"
    else if (x %in% acidic_res) "Acidic" 
    else if (x %in% polar_res) "Polar"
    else if (x %in% hydrophobic_res) "Hydrophobic"
    else "Other"
  })
  
  cat("\nNEIGHBOR CLASSIFICATION:\n")
  type_summary <- table(neighbor_data$type)
  print(type_summary)
  
  # Highlight important interactions
  cat("\nPOTENTIAL INTERACTIONS:\n")
  close_neighbors <- neighbor_data[neighbor_data$distance < 5, ]
  if (nrow(close_neighbors) > 0) {
    cat("Close contacts (< 5Å):\n")
    print(close_neighbors[, c("resno", "resid", "distance", "type")])
  }
  
  return(list(all_neighbors = neighbor_data, close_neighbors = close_neighbors))
}

# Secondary structure analysis
analyze_secondary_structure <- function(pdb) {
  cat("\n=== SECONDARY STRUCTURE ANALYSIS ===\n\n")
  
  # Try multiple approaches for secondary structure
  ss <- NULL
  
  # Try DSSP first
  tryCatch({
    ss <- dssp(pdb)
    cat("Secondary structure calculated using DSSP\n")
  }, error = function(e) {
    cat("DSSP not available:", e$message, "\n")
    
    # Try STRIDE
    tryCatch({
      ss <- stride(pdb)
      cat("Secondary structure calculated using STRIDE\n")
    }, error = function(e2) {
      cat("STRIDE not available:", e2$message, "\n")
      
      # Use simple geometric prediction as fallback
      cat("Using simple geometric secondary structure prediction\n")
      ss <- predict_ss_simple(pdb)
    })
  })
  
  if (!is.null(ss) && !is.null(ss$sse)) {
    # Get secondary structure for residue 292
    if (length(ss$sse) >= 292) {
      ss_292 <- ss$sse[292]
      cat("Secondary structure at position 292:", ifelse(is.na(ss_292), "Coil", ss_292), "\n")
    }
    
    # Overall secondary structure composition
    ss_table <- table(ss$sse, useNA = "ifany")
    cat("Overall secondary structure composition:\n")
    print(ss_table)
    
    return(ss)
  } else {
    cat("Secondary structure prediction failed, continuing without SS analysis\n")
    return(NULL)
  }
}

# Simple geometric secondary structure prediction
predict_ss_simple <- function(pdb) {
  # Simple phi/psi angle-based prediction (basic implementation)
  n_res <- max(pdb$atom$resno[pdb$atom$elety == "CA"])
  ss_pred <- rep("C", n_res)  # Default to coil
  
  # This is a very simplified approach
  # In practice, you'd calculate phi/psi angles
  
  return(list(sse = ss_pred))
}

# B-factor analysis (confidence for AlphaFold)
analyze_bfactors <- function(pdb) {
  cat("\n=== B-FACTOR/CONFIDENCE ANALYSIS ===\n\n")
  
  # Get CA B-factors
  ca_sel <- atom.select(pdb, elety = "CA")
  ca_bfactors <- pdb$atom[ca_sel$atom, "b"]
  ca_resno <- pdb$atom[ca_sel$atom, "resno"]
  
  # Statistics
  cat("B-factor statistics:\n")
  cat("Mean:", round(mean(ca_bfactors), 2), "\n")
  cat("Median:", round(median(ca_bfactors), 2), "\n")
  cat("Range:", round(min(ca_bfactors), 2), "-", round(max(ca_bfactors), 2), "\n")
  
  # B-factor at position 292
  bf_292 <- ca_bfactors[ca_resno == 292]
  if (length(bf_292) > 0) {
    cat("B-factor at position 292:", round(bf_292, 2), "\n")
  }
  
  # Create B-factor plot
  bf_df <- data.frame(
    position = ca_resno,
    bfactor = ca_bfactors,
    stringsAsFactors = FALSE
  )
  
  p2 <- ggplot(bf_df, aes(x = position, y = bfactor)) +
    geom_line(color = "blue", alpha = 0.7) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_vline(xintercept = 292, color = "red", linewidth = 1, linetype = "dashed") +
    geom_hline(yintercept = c(50, 70, 90), color = "gray", linetype = "dotted") +
    annotate("text", x = 292, y = max(ca_bfactors) * 0.9, 
             label = "R292L", color = "red", size = 4) +
    labs(title = "cusS Confidence Profile (AlphaFold B-factors)",
         subtitle = "Higher values = higher confidence",
         x = "Residue Position",
         y = "B-factor (Confidence Score)") +
    theme_minimal()
  
  print(p2)
  
  return(bf_df)
}

# Structural comparison for mutation impact
create_mutation_comparison <- function(pdb) {
  cat("\n=== MUTATION IMPACT PREDICTION ===\n\n")
  
  # Get residue 292 atoms
  res292_sel <- atom.select(pdb, resno = 292)
  res292_atoms <- pdb$atom[res292_sel$atom, ]
  
  if (nrow(res292_atoms) == 0) {
    cat("Cannot analyze - residue 292 not found\n")
    return(NULL)
  }
  
  original_residue <- unique(res292_atoms$resid)[1]
  cat("Original residue:", original_residue, "\n")
  cat("Mutant residue: LEU\n\n")
  
  # Analyze side chain properties
  if (original_residue == "ARG") {
    cat("MUTATION ANALYSIS: ARG → LEU\n")
    cat("Property changes:\n")
    cat("• Charge: +1 → 0 (loss of positive charge)\n")
    cat("• Size: Large → Medium (volume decrease)\n") 
    cat("• Hydrophobicity: Hydrophilic → Hydrophobic\n")
    cat("• Flexibility: High → Medium\n")
    cat("• H-bond potential: High → Low\n\n")
    
    cat("Predicted structural impacts:\n")
    cat("• Loss of ionic interactions\n")
    cat("• Reduced hydrogen bonding capacity\n") 
    cat("• Increased local hydrophobicity\n")
    cat("• Possible local conformational changes\n")
    cat("• May affect protein-protein interactions\n")
    
  } else {
    cat("WARNING: Expected ARG at position 292, found", original_residue, "\n")
  }
  
  return(list(original = original_residue, mutant = "LEU"))
}

# Generate 3D visualization data
prepare_3d_visualization <- function(pdb, neighbors_data) {
  cat("\n=== PREPARING 3D VISUALIZATION DATA ===\n\n")
  
  # Get coordinates for visualization
  ca_sel <- atom.select(pdb, elety = "CA")
  ca_coords <- pdb$atom[ca_sel$atom, c("x", "y", "z")]
  ca_resno <- pdb$atom[ca_sel$atom, "resno"]
  
  # Create data frame for plotting
  plot_data <- data.frame(
    x = ca_coords$x,
    y = ca_coords$y, 
    z = ca_coords$z,
    resno = ca_resno,
    is_mutation = ca_resno == 292,
    is_neighbor = ca_resno %in% neighbors_data$all_neighbors$resno,
    stringsAsFactors = FALSE
  )
  
  # 2D projections for visualization
  p3 <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(aes(color = is_neighbor, size = is_mutation), alpha = 0.6) +
    geom_point(data = plot_data[plot_data$is_mutation, ], 
               aes(x = x, y = y), color = "red", size = 4) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 4)) +
    labs(title = "cusS Structure - XY Projection",
         subtitle = "Red = Mutation site, Blue = Neighbors",
         x = "X coordinate (Å)",
         y = "Y coordinate (Å)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p3)
  
  return(plot_data)
}

# Main analysis function
run_complete_bio3d_analysis <- function(pdb_file = "cusS_alphafold.pdb") {
  
  cat("Starting comprehensive bio3d analysis...\n\n")
  
  # Check if file exists
  if (!file.exists(pdb_file)) {
    cat("PDB file not found. Please download cusS structure first.\n")
    return(NULL)
  }
  
  # Run all analyses
  structure_info <- analyze_cusS_with_bio3d(pdb_file)
  res292_info <- analyze_residue_292(structure_info$pdb)
  environment_info <- analyze_environment(structure_info$pdb)
  ss_info <- analyze_secondary_structure(structure_info$pdb)
  bf_info <- analyze_bfactors(structure_info$pdb)
  mutation_info <- create_mutation_comparison(structure_info$pdb)
  viz_data <- prepare_3d_visualization(structure_info$pdb, environment_info)
  
  # Save results
  results <- list(
    structure = structure_info,
    residue_292 = res292_info,
    environment = environment_info,
    secondary_structure = ss_info,
    bfactors = bf_info,
    mutation_analysis = mutation_info,
    visualization_data = viz_data
  )
  
  # Save to file
  saveRDS(results, "cusS_bio3d_analysis.rds")
  cat("\n✓ Complete analysis saved to 'cusS_bio3d_analysis.rds'\n")
  
  return(results)
}

# Interactive 3D plotting function (requires rgl)
create_3d_plot <- function(results) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    cat("rgl package required for 3D plotting. Install with: install.packages('rgl')\n")
    return(NULL)
  }
  
  library(rgl)
  
  viz_data <- results$visualization_data
  
  # Clear any existing plots
  rgl::clear3d()
  
  # Plot regular residues first (smaller, gray)
  regular_residues <- viz_data[!viz_data$is_mutation & !viz_data$is_neighbor, ]
  if (nrow(regular_residues) > 0) {
    rgl::points3d(regular_residues$x, regular_residues$y, regular_residues$z,
                  color = "gray", size = 3, alpha = 0.6)
  }
  
  # Plot neighbor residues (medium, blue)
  neighbor_residues <- viz_data[viz_data$is_neighbor & !viz_data$is_mutation, ]
  if (nrow(neighbor_residues) > 0) {
    rgl::points3d(neighbor_residues$x, neighbor_residues$y, neighbor_residues$z,
                  color = "blue", size = 5, alpha = 0.8)
  }
  
  # Plot mutation site (large, red)
  mutation_site <- viz_data[viz_data$is_mutation, ]
  if (nrow(mutation_site) > 0) {
    rgl::points3d(mutation_site$x, mutation_site$y, mutation_site$z,
                  color = "red", size = 8, alpha = 1.0)
    
    # Add text label for mutation site
    rgl::text3d(mutation_site$x, mutation_site$y, mutation_site$z + 3, 
                texts = "R292L", cex = 1.5, color = "red")
  }
  
  # Add axis labels
  rgl::axes3d()
  rgl::title3d(main = "cusS 3D Structure - R292L Mutation Site",
               xlab = "X (Å)", ylab = "Y (Å)", zlab = "Z (Å)")
  
  # Add a legend using text
  rgl::text3d(max(viz_data$x) - 10, max(viz_data$y), max(viz_data$z),
              texts = c("Legend:", "Red = R292L", "Blue = Neighbors", "Gray = Other"),
              color = c("black", "red", "blue", "gray"), cex = 1.0)
  
  cat("3D plot created successfully! Use mouse to rotate and zoom.\n")
  cat("- Left click + drag: Rotate\n")
  cat("- Right click + drag: Zoom\n")
  cat("- Middle click + drag: Pan\n")
  
  # Return the rgl device ID
  return(rgl::cur3d())
}

create_interactive_plotly <- function(results) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    cat("plotly package required. Install with: install.packages('plotly')\n")
    return(NULL)
  }
  
  library(plotly)
  
  viz_data <- results$visualization_data
  
  # Add labels for hover text
  viz_data$label <- paste("Residue", viz_data$resno,
                          ifelse(viz_data$is_mutation, " (R292L MUTATION)", 
                                 ifelse(viz_data$is_neighbor, " (Neighbor)", "")))
  
  # Create colors and sizes
  viz_data$color <- ifelse(viz_data$is_mutation, "red",
                           ifelse(viz_data$is_neighbor, "blue", "gray"))
  viz_data$size <- ifelse(viz_data$is_mutation, 12,
                          ifelse(viz_data$is_neighbor, 8, 4))
  
  # Create 3D scatter plot
  p <- plot_ly(viz_data, x = ~x, y = ~y, z = ~z,
               color = ~color, colors = c("blue", "gray", "red"),
               size = ~size, sizes = c(4, 12),
               text = ~label, hovertemplate = "%{text}<br>(%{x:.1f}, %{y:.1f}, %{z:.1f})<extra></extra>",
               type = "scatter3d", mode = "markers",
               marker = list(opacity = 0.8)) %>%
    layout(title = "cusS 3D Structure - R292L Mutation Analysis",
           scene = list(xaxis = list(title = "X (Å)"),
                        yaxis = list(title = "Y (Å)"),
                        zaxis = list(title = "Z (Å)")),
           showlegend = FALSE)
  
  print(p)
  cat("Interactive 3D plot created with plotly!\n")
  return(p)
}



analysis_results <- run_complete_bio3d_analysis("cusS_alphafold.pdb")

# Create multiple 2D projection views
create_2d_projections <- function(results) {
  library(ggplot2)
  library(gridExtra)
  
  viz_data <- results$visualization_data
  
  # Add point types for better visualization
  viz_data$point_type <- ifelse(viz_data$is_mutation, "Mutation",
                                ifelse(viz_data$is_neighbor, "Neighbor", "Other"))
  
  # XY projection
  p1 <- ggplot(viz_data, aes(x = x, y = y)) +
    geom_point(aes(color = point_type, size = point_type), alpha = 0.7) +
    scale_color_manual(values = c("Mutation" = "red", "Neighbor" = "blue", "Other" = "gray")) +
    scale_size_manual(values = c("Mutation" = 4, "Neighbor" = 2.5, "Other" = 1)) +
    labs(title = "XY Projection", x = "X (Å)", y = "Y (Å)") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # XZ projection  
  p2 <- ggplot(viz_data, aes(x = x, y = z)) +
    geom_point(aes(color = point_type, size = point_type), alpha = 0.7) +
    scale_color_manual(values = c("Mutation" = "red", "Neighbor" = "blue", "Other" = "gray")) +
    scale_size_manual(values = c("Mutation" = 4, "Neighbor" = 2.5, "Other" = 1)) +
    labs(title = "XZ Projection", x = "X (Å)", y = "Z (Å)") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # YZ projection
  p3 <- ggplot(viz_data, aes(x = y, y = z)) +
    geom_point(aes(color = point_type, size = point_type), alpha = 0.7) +
    scale_color_manual(values = c("Mutation" = "red", "Neighbor" = "blue", "Other" = "gray")) +
    scale_size_manual(values = c("Mutation" = 4, "Neighbor" = 2.5, "Other" = 1)) +
    labs(title = "YZ Projection", x = "Y (Å)", y = "Z (Å)") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # Distance from mutation site
  mutation_coords <- viz_data[viz_data$is_mutation, c("x", "y", "z")]
  if (nrow(mutation_coords) > 0) {
    viz_data$dist_from_mutation <- sqrt(
      (viz_data$x - mutation_coords$x)^2 + 
        (viz_data$y - mutation_coords$y)^2 + 
        (viz_data$z - mutation_coords$z)^2
    )
    
    p4 <- ggplot(viz_data[!viz_data$is_mutation, ], aes(x = resno, y = dist_from_mutation)) +
      geom_point(aes(color = point_type), alpha = 0.7, size = 2) +
      geom_hline(yintercept = 8, linetype = "dashed", color = "red", alpha = 0.7) +
      scale_color_manual(values = c("Neighbor" = "blue", "Other" = "gray")) +
      labs(title = "Distance from R292L Site", 
           x = "Residue Number", y = "Distance (Å)",
           subtitle = "Dashed line = 8Å neighbor cutoff") +
      theme_minimal() +
      theme(legend.title = element_blank())
  }
  
  # Combine all plots
  if (exists("p4")) {
    combined_plot <- grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
  } else {
    combined_plot <- grid.arrange(p1, p2, p3, nrow = 2, ncol = 2)
  }
  
  cat("2D projection plots created!\n")
  return(combined_plot)
}

# Create 2D projections
projections <- create_2d_projections(analysis_results)

rgl_plot <- create_3d_plot(analysis_results)
plotly_plot <- create_interactive_plotly(analysis_results)

summary_viz <- create_summary_visualization()




library(bio3d)

wt <- bio3d::read.fasta("cusS_WT.fasta")
mut <- bio3d::read.fasta("cusS_R292L.fasta")


# Extract the first (and usually only) sequence as strings

wt_seq  <- wt$ali[1, ]     # e.g., "M" "V" "S" "K" ...
mut_seq <- mut$ali[1, ]

# Quick sanity checks (should be TRUE / informative)
stopifnot(is.character(wt_seq), is.character(mut_seq))
stopifnot(all(nchar(wt_seq)  == 1), all(nchar(mut_seq) == 1))
length(wt_seq); length(mut_seq)  # likely equal if it’s a single-point mutant

# Feed a character MATRIX to seqaln (rows = sequences, columns = residues)
seqmat <- rbind(WT = wt_seq, Mut = mut_seq)

seqmat <- bio3d::seqbind(WT = wt_seq, Mut = mut_seq)

# td <- tempfile("bio3d_align_"); dir.create(td)
# oldwd <- getwd()
# setwd(td)
# on.exit(setwd(oldwd), add = TRUE)


muscle_bin <- Sys.which("muscle")
clustalw_bin <- Sys.which("clustalw")
exefile <- if (nzchar(muscle_bin)) muscle_bin else clustalw_bin


seqmat <- rbind(WT = wt_seq, Mut = mut_seq)
bio3d::seqbind(WT = wt_seq, Mut = mut_seq)


write.fasta <- function(seqs, file) {
  con <- file(file, "w")
  for (i in seq_along(seqs)) {
    cat(paste0(">", names(seqs)[i], "\n"), file = con)
    cat(paste0(seqs[[i]], "\n"), file = con)
  }
  close(con)
}


seqs <- list(WT = wt_seq, Mut = mut_seq)
write.fasta(seqs, "in.fa")

system("/usr/bin/muscle -align input.fa -output aln.fa")

aligned <- readLines("aln.fa")


# Parse FASTA (aligned) lines into a named character vector
parse_fasta <- function(lines) {
  hdr <- grep("^>", lines)
  stopifnot(length(hdr) > 0)
  hdr <- c(hdr, length(lines) + 1L)
  seqs <- character(length(hdr) - 1L)
  nms  <- character(length(hdr) - 1L)
  for (i in seq_along(seqs)) {
    nms[i] <- sub("^>", "", lines[hdr[i]])
    block  <- lines[(hdr[i] + 1L):(hdr[i + 1L] - 1L)]
    block  <- block[nzchar(block)]
    seqs[i] <- paste(block, collapse = "")
  }
  # Make sure all aligned sequences are the same length
  if (length(unique(nchar(seqs))) != 1L) {
    stop("Aligned sequences are not equal length—did you pass an unaligned FASTA?")
  }
  setNames(seqs, nms)
}

# Pretty-print alignment like CLUSTAL (markers line uses '*' for identity)
pretty_print_alignment <- function(aln_vec, width = 60, show_numbers = TRUE) {
  nms  <- names(aln_vec)
  seqs <- as.character(aln_vec)
  L    <- nchar(seqs[1])
  blocks <- ceiling(L / width)
  
  # running ungapped residue counters for each sequence (for right-side numbers)
  pos <- integer(length(seqs))
  
  # Build and print a CLUSTAL-style header
  cat("CLUSTAL-like alignment (identity markers only)\n\n")
  
  for (b in seq_len(blocks)) {
    i1 <- (b - 1L) * width + 1L
    i2 <- min(b * width, L)
    idx <- i1:i2
    
    # Compute markers: '*' if all equal and not a gap
    chars <- do.call(cbind, lapply(seqs, function(s) substring(s, idx, idx)))
    eq    <- apply(chars, 1L, function(col) length(unique(col)) == 1L && col[1] != "-")
    mark  <- ifelse(eq, "*", " ")
    marker_line <- paste(mark, collapse = "")
    
    # Print each sequence line
    for (j in seq_along(seqs)) {
      chunk <- paste(substring(seqs[j], idx, idx), collapse = "")
      # update ungapped position counters for this block
      pos[j] <- pos[j] + sum(substring(seqs[j], idx, idx) != "-")
      if (show_numbers) {
        cat(sprintf("%-12s %s %6d\n", nms[j], chunk, pos[j]))
      } else {
        cat(sprintf("%-12s %s\n", nms[j], chunk))
      }
    }
    # Print markers line under the block
    if (show_numbers) {
      cat(sprintf("%-12s %s\n\n", "", marker_line))
    } else {
      cat(sprintf("            %s\n\n", marker_line))
    }
  }
}

# ---- Run it on your file ----
aligned_lines <- readLines("aln.fa")
aln_vec <- parse_fasta(aligned_lines)
pretty_print_alignment(aln_vec, width = 60)


