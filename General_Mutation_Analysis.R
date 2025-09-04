setwd("/home/william-ackerman/Desktop/cusS_R292L")

# Load required libraries
library(bio3d)
library(tidyverse)
library(gridExtra)

extract_gene_from_genbank <- function(gbk_file, gene_name, mutation_type = "point") {
  gbk_lines <- readLines(gbk_file)
  
  cat("Searching for gene:", gene_name, "\n")
  
  # Find the gene - try multiple search strategies
  gene_patterns <- c(
    paste0('/gene="', gene_name, '"'),
    paste0('gene="', gene_name, '"'),
    paste0('/locus_tag=".*', gene_name),
    gene_name
  )
  
  gene_line <- NULL
  for (pattern in gene_patterns) {
    matches <- which(grepl(pattern, gbk_lines, ignore.case = TRUE))
    if (length(matches) > 0) {
      gene_line <- matches[1]
      cat("Found gene using pattern:", pattern, "at line", gene_line, "\n")
      break
    }
  }
  
  if (is.null(gene_line)) {
    # Broader search
    broad_matches <- which(grepl(gene_name, gbk_lines, ignore.case = TRUE))
    if (length(broad_matches) > 0) {
      cat("Gene not found with strict patterns. Potential matches at lines:\n")
      for (i in 1:min(10, length(broad_matches))) {
        cat("Line", broad_matches[i], ":", trimws(gbk_lines[broad_matches[i]]), "\n")
      }
      return(list(found = FALSE, matches = broad_matches))
    }
    stop("Gene not found")
  }
  
  # Extract gene information
  gene_info <- extract_gene_details(gbk_lines, gene_line)
  
  return(list(
    found = TRUE,
    gene_info = gene_info,
    gene_line = gene_line
  ))
}

# Extract detailed gene information
extract_gene_details <- function(gbk_lines, start_line) {
  info <- list()
  
  # Look backwards and forwards for CDS/gene features
  search_range <- max(1, start_line - 50):min(length(gbk_lines), start_line + 100)
  
  # Extract coordinates
  for (i in search_range) {
    line <- gbk_lines[i]
    
    # Look for CDS coordinates
    if (grepl("^\\s*CDS\\s+", line)) {
      coords <- extract_coordinates(line)
      if (!is.null(coords)) {
        info$coordinates <- coords
        info$strand <- coords$strand
      }
    }
    
    # Extract gene name
    if (grepl('/gene=', line)) {
      info$gene_name <- gsub('.*\\/gene="([^"]*)".*', '\\1', line)
    }
    
    # Extract product
    if (grepl('/product=', line)) {
      product_line <- line
      j <- i
      # Handle multi-line products
      while (j <= length(gbk_lines) && !grepl('".*"', product_line)) {
        j <- j + 1
        if (j <= length(gbk_lines)) {
          product_line <- paste(product_line, gbk_lines[j])
        }
      }
      info$product <- gsub('.*\\/product="([^"]*)".*', '\\1', product_line)
    }
    
    # Extract translation if available
    if (grepl('/translation=', line)) {
      info$translation <- extract_translation(gbk_lines, i)
    }
  }
  
  return(info)
}

# Extract coordinate information
extract_coordinates <- function(coord_line) {
  coords <- list()
  
  # Handle different coordinate formats
  if (grepl("complement", coord_line)) {
    coords$strand <- "reverse"
    coord_match <- regmatches(coord_line, regexpr("\\d+\\.\\.\\d+", coord_line))
  } else {
    coords$strand <- "forward"
    coord_match <- regmatches(coord_line, regexpr("\\d+\\.\\.\\d+", coord_line))
  }
  
  if (length(coord_match) > 0) {
    coord_parts <- strsplit(coord_match, "\\.\\.")[[1]]
    coords$start <- as.numeric(coord_parts[1])
    coords$end <- as.numeric(coord_parts[2])
    coords$length_bp <- coords$end - coords$start + 1
    coords$length_aa <- floor(coords$length_bp / 3)
  }
  
  return(coords)
}

# Extract translation sequence
extract_translation <- function(gbk_lines, start_line) {
  aa_sequence <- ""
  i <- start_line
  
  # Get the first line
  first_line <- gbk_lines[i]
  seq_part <- sub('.*\\/translation="([^"]*)".*', '\\1', first_line)
  
  if (grepl('".*"', first_line)) {
    # Single line translation
    aa_sequence <- seq_part
  } else {
    # Multi-line translation
    seq_part <- sub('.*\\/translation="', '', first_line)
    aa_sequence <- seq_part
    
    i <- i + 1
    while (i <= length(gbk_lines) && !grepl('"', gbk_lines[i])) {
      line_seq <- gsub('^\\s+', '', gbk_lines[i])
      aa_sequence <- paste0(aa_sequence, line_seq)
      i <- i + 1
    }
    
    if (i <= length(gbk_lines)) {
      last_line <- gsub('".*', '', gbk_lines[i])
      last_line <- gsub('^\\s+', '', last_line)
      aa_sequence <- paste0(aa_sequence, last_line)
    }
  }
  
  # Clean up the sequence
  aa_sequence <- gsub('"', '', aa_sequence)
  aa_sequence <- gsub('\\s+', '', aa_sequence)
  
  return(aa_sequence)
}

# MUTATION ANALYSIS FRAMEWORK 

# Create a mutation object
create_mutation <- function(gene_name, mutation_desc, mutation_type = "point") {
  mutation <- list(
    gene_name = gene_name,
    description = mutation_desc,
    type = mutation_type
  )
  
  # Parse different mutation types
  if (mutation_type == "point") {
    # e.g., "R292L"
    if (grepl("^[A-Z]\\d+[A-Z]$", mutation_desc)) {
      mutation$from_aa <- substr(mutation_desc, 1, 1)
      mutation$position <- as.numeric(gsub("[A-Z]", "", mutation_desc))
      mutation$to_aa <- substr(mutation_desc, nchar(mutation_desc), nchar(mutation_desc))
    }
  } else if (mutation_type == "deletion") {
    # e.g., "Δ8 bp (704‑711/1104 nt)"
    if (grepl("Δ.*bp", mutation_desc)) {
      bp_match <- regmatches(mutation_desc, regexpr("Δ\\d+", mutation_desc))
      mutation$deleted_bp <- as.numeric(gsub("Δ", "", bp_match))
      
      # Extract position range
      if (grepl("\\(\\d+‑\\d+", mutation_desc)) {
        range_match <- regmatches(mutation_desc, regexpr("\\d+‑\\d+", mutation_desc))
        range_parts <- strsplit(range_match, "‑")[[1]]
        mutation$start_pos <- as.numeric(range_parts[1])
        mutation$end_pos <- as.numeric(range_parts[2])
      }
      
      # Extract total length
      if (grepl("/\\d+\\s*nt", mutation_desc)) {
        total_match <- regmatches(mutation_desc, regexpr("\\d+(?=\\s*nt)", mutation_desc, perl = TRUE))
        mutation$total_length <- as.numeric(total_match)
      }
    }
  }
  
  class(mutation) <- c("mutation", mutation_type)
  return(mutation)
}

# Analyze any mutation
analyze_mutation <- function(mutation, gbk_file = "Ecoli_K12_MG1655.gbk") {
  cat("=== MUTATION ANALYSIS ===\n")
  cat("Gene:", mutation$gene_name, "\n")
  cat("Mutation:", mutation$description, "\n")
  cat("Type:", mutation$type, "\n\n")
  
  # Extract gene information
  gene_result <- extract_gene_from_genbank(gbk_file, mutation$gene_name)
  
  if (!gene_result$found) {
    cat("Gene not found in GenBank file\n")
    return(NULL)
  }
  
  gene_info <- gene_result$gene_info
  
  # Print gene information
  cat("=== GENE INFORMATION ===\n")
  if (!is.null(gene_info$product)) {
    cat("Product:", gene_info$product, "\n")
  }
  if (!is.null(gene_info$coordinates)) {
    cat("Location:", gene_info$coordinates$start, "-", gene_info$coordinates$end, 
        "(", gene_info$coordinates$strand, "strand )\n")
    cat("Length:", gene_info$coordinates$length_bp, "bp,", 
        gene_info$coordinates$length_aa, "aa\n")
  }
  
  # Type-specific analysis
  if (mutation$type == "point") {
    analyze_point_mutation(mutation, gene_info)
  } else if (mutation$type == "deletion") {
    analyze_deletion_mutation(mutation, gene_info)
  }
  
  return(list(mutation = mutation, gene_info = gene_info))
}

# Analyze point mutations
analyze_point_mutation <- function(mutation, gene_info) {
  cat("\n=== POINT MUTATION ANALYSIS ===\n")
  
  if (is.null(gene_info$translation)) {
    cat("No protein sequence available for analysis\n")
    return(NULL)
  }
  
  wt_seq <- gene_info$translation
  pos <- mutation$position
  from_aa <- mutation$from_aa
  to_aa <- mutation$to_aa
  
  if (nchar(wt_seq) < pos) {
    cat("Mutation position", pos, "is beyond protein length", nchar(wt_seq), "\n")
    return(NULL)
  }
  
  actual_aa <- substr(wt_seq, pos, pos)
  cat("Expected amino acid at position", pos, ":", from_aa, "\n")
  cat("Actual amino acid at position", pos, ":", actual_aa, "\n")
  
  if (actual_aa == from_aa) {
    cat("✓ Mutation notation matches sequence\n")
  } else {
    cat("⚠ Mutation notation does not match sequence!\n")
  }
  
  # Create mutant sequence
  mut_seq <- wt_seq
  substr(mut_seq, pos, pos) <- to_aa
  
  # Analyze amino acid properties
  analyze_aa_change(from_aa, to_aa, pos)
  
  # Write sequences
  write_sequences(gene_info$gene_name, wt_seq, mut_seq, paste0(from_aa, pos, to_aa))
  
  return(list(wt_seq = wt_seq, mut_seq = mut_seq))
}

# Analyze deletion mutations
analyze_deletion_mutation <- function(mutation, gene_info) {
  cat("\n=== DELETION MUTATION ANALYSIS ===\n")
  
  cat("Deleted region:", mutation$start_pos, "-", mutation$end_pos, "\n")
  cat("Deleted length:", mutation$deleted_bp, "bp\n")
  
  if (!is.null(mutation$total_length)) {
    cat("Total gene length:", mutation$total_length, "nt\n")
    cat("Percentage deleted:", round(100 * mutation$deleted_bp / mutation$total_length, 1), "%\n")
  }
  
  # Determine impact on reading frame
  if (mutation$deleted_bp %% 3 == 0) {
    cat("In-frame deletion (maintains reading frame)\n")
    deleted_aa <- mutation$deleted_bp / 3
    cat("Amino acids deleted:", deleted_aa, "\n")
  } else {
    cat("Frameshift deletion (disrupts reading frame)\n")
    cat("This will likely cause a premature stop codon\n")
  }
  
  # Analyze position context
  if (!is.null(gene_info$coordinates)) {
    rel_start <- (mutation$start_pos - gene_info$coordinates$start + 1) / gene_info$coordinates$length_bp
    rel_end <- (mutation$end_pos - gene_info$coordinates$start + 1) / gene_info$coordinates$length_bp
    
    cat("Relative position in gene:\n")
    cat("Start:", round(rel_start * 100, 1), "%\n")
    cat("End:", round(rel_end * 100, 1), "%\n")
    
    if (rel_start < 0.1) {
      cat("⚠ Deletion near start of gene - likely severe impact\n")
    } else if (rel_start > 0.9) {
      cat("Deletion near end of gene - may have less impact\n")
    } else {
      cat("Deletion in middle of gene\n")
    }
  }
  
  return(mutation)
}

# Analyze amino acid changes
analyze_aa_change <- function(from_aa, to_aa, position) {
  cat("\n=== AMINO ACID CHANGE ANALYSIS ===\n")
  
  # Amino acid properties
  aa_props <- data.frame(
    AA = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", 
           "F", "P", "S", "T", "W", "Y", "V"),
    Charge = c(0, 1, 0, -1, 0, 0, -1, 0, 0.5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    Hydrophobicity = c(1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5, 
                       3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2),
    Size = c("S", "L", "M", "M", "S", "M", "M", "S", "M", "M", "L", "L", "M", 
             "L", "M", "S", "M", "L", "L", "M"),
    Type = c("Nonpolar", "Basic", "Polar", "Acidic", "Special", "Polar", "Acidic",
             "Special", "Basic", "Nonpolar", "Nonpolar", "Basic", "Nonpolar",
             "Aromatic", "Special", "Polar", "Polar", "Aromatic", "Aromatic", "Nonpolar"),
    stringsAsFactors = FALSE
  )
  
  from_props <- aa_props[aa_props$AA == from_aa, ]
  to_props <- aa_props[aa_props$AA == to_aa, ]
  
  if (nrow(from_props) > 0 && nrow(to_props) > 0) {
    cat("Change:", from_aa, "→", to_aa, "at position", position, "\n")
    cat("Charge change:", from_props$Charge, "→", to_props$Charge, "\n")
    cat("Hydrophobicity change:", from_props$Hydrophobicity, "→", to_props$Hydrophobicity, "\n")
    cat("Type change:", from_props$Type, "→", to_props$Type, "\n")
    
    # Assess severity
    charge_change <- abs(from_props$Charge - to_props$Charge)
    hydro_change <- abs(from_props$Hydrophobicity - to_props$Hydrophobicity)
    
    severity_score <- charge_change * 2 + hydro_change
    
    if (severity_score < 2) {
      severity <- "Conservative"
    } else if (severity_score < 5) {
      severity <- "Moderate" 
    } else {
      severity <- "Radical"
    }
    
    cat("Change severity:", severity, "(score:", round(severity_score, 1), ")\n")
  }
}

# Write sequences to files
write_sequences <- function(gene_name, wt_seq, mut_seq = NULL, mut_desc = "") {
  # Write wild-type
  wt_filename <- paste0(gene_name, "_WT.fasta")
  writeLines(c(paste0(">", gene_name, "_WT_E.coli_K12"), wt_seq), wt_filename)
  cat("Wild-type sequence written to:", wt_filename, "\n")
  
  # Write mutant if provided
  if (!is.null(mut_seq)) {
    mut_filename <- paste0(gene_name, "_", mut_desc, ".fasta")
    writeLines(c(paste0(">", gene_name, "_", mut_desc, "_mutant"), mut_seq), mut_filename)
    cat("Mutant sequence written to:", mut_filename, "\n")
  }
}

# EXAMPLE ANALYSIS
cat("Analyzing ompC deletion mutation...\n\n")

# Create mutation object
ompC_mutation <- create_mutation(
  gene_name = "ompC",
  mutation_desc = "Δ8 bp (704‑711/1104 nt)",
  mutation_type = "deletion"
)

# Analyze the mutation
ompC_analysis <- analyze_mutation(ompC_mutation)

# You can also still analyze the original cusS mutation:
cat("\n\nAnalyzing original cusS mutation...\n\n")

cusS_mutation <- create_mutation(
  gene_name = "cusS", 
  mutation_desc = "R292L",
  mutation_type = "point"
)

cusS_analysis <- analyze_mutation(cusS_mutation)

# Function to create mutation summary report
create_mutation_report <- function(analyses) {
  cat("\n", strrep("=", 80), "\n")
  cat("                    MUTATION ANALYSIS SUMMARY REPORT\n")
  cat(strrep("=", 80), "\n\n")
  
  for (i in seq_along(analyses)) {
    analysis <- analyses[[i]]
    if (is.null(analysis)) next
    
    mut <- analysis$mutation
    gene <- analysis$gene_info
    
    cat("MUTATION", i, ":\n")
    cat("Gene:", mut$gene_name, "\n")
    cat("Type:", mut$type, "\n")
    cat("Description:", mut$description, "\n")
    
    if (!is.null(gene$product)) {
      cat("Product:", gene$product, "\n")
    }
    
    if (mut$type == "deletion") {
      cat("Impact: Likely", ifelse(mut$deleted_bp %% 3 == 0, "in-frame", "frameshift"), "deletion\n")
    }
    
    cat(strrep("-", 40), "\n\n")
  }
  
  cat("Files generated:\n")
  fasta_files <- list.files(pattern = "*.fasta")
  for (file in fasta_files) {
    cat("•", file, "\n")
  }
}

# Create summary report
create_mutation_report(list(ompC_analysis, cusS_analysis))

