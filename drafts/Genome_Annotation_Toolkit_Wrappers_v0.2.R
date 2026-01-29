#' Genome Annotation Toolkit
#' A collection of functions for working with genomic annotations and sequences
#' 
#' Required packages:
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(Rsamtools)
library(rentrez)
library(httr2)
library(xml2)
library(tidyverse)
library(cli)

# Setup and Data Loading Functions ====

#' Initialize Genome Resources
#'
#' @param gff_path Path to GFF3 file
#' @param fasta_path Path to FASTA file
#' @param auto_index Logical; automatically index FASTA if needed
#' @return List with gff (GRanges), fasta (DNAStringSet), and fa (FaFile) objects
#' @export
init_genome <- function(gff_path, fasta_path, auto_index = TRUE) {
  
  message("Loading GFF3 annotation...")
  gff <- import(gff_path)
  
  message("Loading FASTA sequence...")
  fasta <- readDNAStringSet(fasta_path)
  
  # Index FASTA if needed
  if (auto_index && !file.exists(paste0(fasta_path, ".fai"))) {
    message("Indexing FASTA file...")
    indexFa(fasta_path)
  }
  
  fa <- FaFile(fasta_path)
  open(fa)
  
  message("Genome resources loaded successfully!")
  
  list(
    gff = gff,
    fasta = fasta,
    fa = fa,
    seqnames = seqnames(seqinfo(fa))
  )
}

# Sequence Extraction Functions

#' Extract Sequences by Gene Name
#'
#' @param genome_obj Genome object from init_genome()
#' @param gene_pattern Pattern to match in Name field (regex supported)
#' @param translate Logical; translate to amino acid sequence
#' @param genetic_code Genetic code to use (default = "11" for bacteria)
#' @return DNAStringSet or AAStringSet with extracted sequences
#' @export
extract_sequences_by_name <- function(genome_obj, 
                                      gene_pattern, 
                                      translate = FALSE,
                                      genetic_code = "11") {
  
  # Find matching features
  features <- genome_obj$gff[grepl(gene_pattern, genome_obj$gff$Name)]
  
  if (length(features) == 0) {
    warning("No features found matching pattern: ", gene_pattern)
    return(NULL)
  }
  
  message(sprintf("Found %d feature(s) matching '%s'", length(features), gene_pattern))
  
  # Fix seqnames if there's only one chromosome
  fa_headers <- genome_obj$seqnames
  if (length(fa_headers) == 1) {
    seqlevels(features) <- fa_headers
    seqnames(features) <- Rle(fa_headers)
  }
  
  # Sanity check
  if (!all(as.character(seqnames(features)) %in% fa_headers)) {
    stop("Seqnames mismatch between GFF and FASTA!")
  }
  
  # Extract DNA sequences (strand-aware)
  dna_seqs <- getSeq(genome_obj$fa, features)
  names(dna_seqs) <- features$Name
  
  if (translate) {
    message("Translating to amino acids...")
    aa_seqs <- translate(dna_seqs, genetic.code = getGeneticCode(genetic_code))
    return(aa_seqs)
  }
  
  return(dna_seqs)
}


#' Extract Sequence from Genomic Coordinates
#'
#' @param genome_obj Genome object from init_genome()
#' @param seqname Chromosome/contig name
#' @param start Start coordinate
#' @param end End coordinate
#' @param strand Strand ("+" or "-"); default "+"
#' @return DNAStringSet with extracted sequence
#' @export
extract_sequence_by_coords <- function(genome_obj, seqname, start, end, strand = "+") {
  
  roi <- GRanges(
    seqnames = seqname,
    ranges = IRanges(start = start, end = end),
    strand = strand
  )
  
  seq <- getSeq(genome_obj$fa, roi)
  names(seq) <- sprintf("%s:%d-%d(%s)", seqname, start, end, strand)
  
  return(seq)
}


# Genomic Context and Region Analysis Functions

#' Get Genomic Context Around Features
#'
#' @param genome_obj Genome object from init_genome()
#' @param features GRanges object or gene name pattern
#' @param flank_size Flanking region size (bp) on each side
#' @param feature_filter Optional regex pattern to filter context features
#' @return GRanges with features in the flanking regions
#' @export
get_genomic_context <- function(genome_obj, 
                                features, 
                                flank_size = 20000,
                                feature_filter = NULL) {
  
  # If features is a string, find matching features
  if (is.character(features)) {
    features <- genome_obj$gff[grepl(features, genome_obj$gff$Name)]
  }
  
  if (length(features) == 0) {
    stop("No features provided or found")
  }
  
  # Expand features to include flanks
  ctx_region <- resize(features, 
                       width = width(features) + 2 * flank_size, 
                       fix = "center")
  
  # Find overlapping features
  ctx_feats <- subsetByOverlaps(genome_obj$gff, ctx_region)
  
  # Apply filter if provided
  if (!is.null(feature_filter)) {
    ctx_feats <- ctx_feats[grepl(feature_filter, ctx_feats$Name, ignore.case = TRUE)]
  }
  
  message(sprintf("Found %d features in %d bp flanking region(s)", 
                  length(ctx_feats), flank_size))
  
  return(ctx_feats)
}


#' Analyze Region of Interest
#'
#' @param genome_obj Genome object from init_genome()
#' @param seqname Chromosome/contig name
#' @param start Start coordinate
#' @param end End coordinate
#' @param flank Flanking region size (default = 0)
#' @param feature_type Feature type to extract (default = "CDS")
#' @param tidy Return tidy data frame (default = TRUE)
#' @return GRanges or tibble with features in ROI
#' @export
analyze_roi <- function(genome_obj, 
                        seqname, 
                        start, 
                        end, 
                        flank = 0,
                        feature_type = "CDS",
                        tidy = TRUE) {
  
  # Define ROI
  roi <- GRanges(
    seqnames = seqname,
    ranges = IRanges(start = start - flank, end = end + flank)
  )
  
  message(sprintf("Analyzing region: %s:%d-%d (±%d bp)", 
                  seqname, start, end, flank))
  
  # Get overlapping features
  roi_feats <- subsetByOverlaps(genome_obj$gff, roi)
  
  # Filter by feature type if specified
  if (!is.null(feature_type)) {
    roi_feats <- roi_feats[roi_feats$type == feature_type]
  }
  
  message(sprintf("Found %d %s feature(s) in ROI", 
                  length(roi_feats), feature_type))
  
  if (!tidy) {
    return(roi_feats)
  }
  
  # Convert to tidy format
  return(tidy_features(roi_feats))
}


#' Search Genome-wide for Features
#'
#' @param genome_obj Genome object from init_genome()
#' @param pattern Regex pattern to search in Name or Note fields
#' @param field Field to search ("Name", "Note", "Alias", or "any")
#' @param feature_type Feature type filter (default = "CDS")
#' @param tidy Return tidy data frame (default = TRUE)
#' @return GRanges or tibble with matching features
#' @export
#' 
search_features <- function(genome_obj, 
                            pattern, 
                            field = "Note",
                            feature_type = "CDS",
                            tidy = TRUE) {
  
  gff <- genome_obj$gff
  
  # Filter by feature type
  if (!is.null(feature_type)) {
    gff <- gff[gff$type == feature_type]
  }
  
  # Search in specified field(s)
  if (field == "any") {
    matches <- grepl(pattern, gff$Name) | 
      grepl(pattern, flatten_CharacterList(gff$Note)) |
      grepl(pattern, flatten_CharacterList(gff$Alias))
  } else if (field == "Note") {
    matches <- grepl(pattern, flatten_CharacterList(gff$Note))
  } else if (field == "Alias") {
    matches <- grepl(pattern, flatten_CharacterList(gff$Alias))
  } else if (field == "Name") {
    matches <- grepl(pattern, flatten_CharacterList(gff$Name))
  } else {
    matches <- grepl(pattern, gff[[field]])
  }
  
  result <- gff[matches]
  
  message(sprintf("Found %d feature(s) matching '%s' in %s field(s)", 
                  length(result), pattern, field))
  
  if (!tidy) {
    return(result)
  }
  
  return(tidy_features(result))
}

#' Extract FASTA Sequence from Region of Interest
#'
#' @param genome_obj Genome object from init_genome()
#' @param seqname Chromosome/contig name
#' @param start Start coordinate
#' @param end End coordinate
#' @param flank Flanking region size in bp (default = 0)
#' @param strand Strand to extract ("+" or "-" or "*" for both); default = "+"
#' @param reverse_complement Apply reverse complement for minus strand (default = TRUE)
#' @param output_file Optional file path to write FASTA (default = NULL)
#' @param seq_name Custom sequence name for FASTA header (default = auto-generated)
#' @return DNAStringSet with extracted sequence(s)
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic extraction
#' roi_seq <- get_roi_fasta(genome, "1", 1828501, 1978000)
#' 
#' # With flanking sequence
#' roi_seq <- get_roi_fasta(genome, "1", 1828501, 1978000, flank = 5000)
#' 
#' # Minus strand with reverse complement
#' roi_seq <- get_roi_fasta(genome, "1", 1828501, 1978000, strand = "-")
#' 
#' # Save to file
#' get_roi_fasta(genome, "1", 1828501, 1978000, 
#'               flank = 10000, 
#'               output_file = "roi_region.fasta")
#' }
get_roi_fasta <- function(genome_obj,
                          seqname,
                          start,
                          end,
                          flank = 0,
                          strand = "+",
                          reverse_complement = TRUE,
                          output_file = NULL,
                          seq_name = NULL) {
  
  # Input validation
  if (start > end) {
    stop("Start coordinate must be <= end coordinate")
  }
  
  if (flank < 0) {
    stop("Flank size must be >= 0")
  }
  
  if (!strand %in% c("+", "-", "*")) {
    stop("Strand must be '+', '-', or '*'")
  }
  
  # Calculate actual coordinates with flanking
  actual_start <- max(1, start - flank)
  actual_end <- end + flank
  
  # Get sequence length to validate end coordinate
  si <- seqinfo(genome_obj$fa)
  seq_lengths <- seqlengths(si)
  
  if (!seqname %in% names(seq_lengths)) {
    stop("Seqname '", seqname, "' not found in FASTA. Available: ", 
         paste(names(seq_lengths), collapse = ", "))
  }
  
  # Adjust end if it exceeds chromosome length
  max_length <- seq_lengths[seqname]
  if (actual_end > max_length) {
    warning(sprintf("End coordinate (%d) exceeds sequence length (%d). Adjusting to %d",
                    actual_end, max_length, max_length))
    actual_end <- max_length
  }
  
  # Generate sequence name if not provided
  if (is.null(seq_name)) {
    if (flank > 0) {
      seq_name <- sprintf("%s:%d-%d_flank%d(%s)", 
                          seqname, start, end, flank, strand)
    } else {
      seq_name <- sprintf("%s:%d-%d(%s)", 
                          seqname, start, end, strand)
    }
  }
  
  # Create GRanges object for extraction
  roi <- GRanges(
    seqnames = seqname,
    ranges = IRanges(start = actual_start, end = actual_end),
    strand = strand
  )
  
  # Extract sequence
  message(sprintf("Extracting sequence from %s:%d-%d (±%d bp, strand: %s)",
                  seqname, start, end, flank, strand))
  
  seq <- getSeq(genome_obj$fa, roi)
  
  # Apply reverse complement if on minus strand and requested
  if (strand == "-" && reverse_complement) {
    message("Applying reverse complement for minus strand")
    seq <- reverseComplement(seq)
  }
  
  # Set sequence name
  names(seq) <- seq_name
  
  # Report sequence info
  message(sprintf("Extracted %d bp (requested region: %d bp, with flank: %d bp)",
                  width(seq), end - start + 1, actual_end - actual_start + 1))
  
  # Write to file if requested
  if (!is.null(output_file)) {
    message("Writing sequence to: ", output_file)
    writeXStringSet(seq, filepath = output_file, format = "fasta")
  }
  
  return(seq)
}


#' Extract Multiple ROI Sequences as Batch
#'
#' @param genome_obj Genome object from init_genome()
#' @param roi_table Data frame with columns: seqname, start, end, and optionally strand, name
#' @param flank Flanking region size in bp (default = 0)
#' @param output_file Optional file path to write all sequences (default = NULL)
#' @return DNAStringSet with all extracted sequences
#' @export
#'
#' @examples
#' \dontrun{
#' # Define multiple regions
#' regions <- tibble(
#'   seqname = c("1", "1", "1"),
#'   start = c(100000, 200000, 300000),
#'   end = c(150000, 250000, 350000),
#'   name = c("region1", "region2", "region3")
#' )
#' 
#' # Extract all at once
#' seqs <- get_roi_fasta_batch(genome, regions, flank = 1000)
#' 
#' # Save to file
#' get_roi_fasta_batch(genome, regions, 
#'                     flank = 1000, 
#'                     output_file = "all_regions.fasta")
#' }
get_roi_fasta_batch <- function(genome_obj,
                                roi_table,
                                flank = 0,
                                output_file = NULL) {
  
  # Validate input table
  required_cols <- c("seqname", "start", "end")
  if (!all(required_cols %in% names(roi_table))) {
    stop("roi_table must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  message(sprintf("Extracting %d ROI sequences...", nrow(roi_table)))
  
  # Extract sequences for each ROI
  seq_list <- lapply(seq_len(nrow(roi_table)), function(i) {
    row <- roi_table[i, ]
    
    # Get strand if available
    strand <- if ("strand" %in% names(row)) row$strand else "+"
    
    # Get custom name if available
    seq_name <- if ("name" %in% names(row)) {
      if (flank > 0) {
        sprintf("%s_%s:%d-%d_flank%d", row$name, row$seqname, row$start, row$end, flank)
      } else {
        sprintf("%s_%s:%d-%d", row$name, row$seqname, row$start, row$end)
      }
    } else {
      NULL  # Will auto-generate in get_roi_fasta
    }
    
    # Extract sequence
    seq <- get_roi_fasta(
      genome_obj = genome_obj,
      seqname = row$seqname,
      start = row$start,
      end = row$end,
      flank = flank,
      strand = strand,
      seq_name = seq_name
    )
    
    return(seq)
  })
  
  # Combine all sequences
  all_seqs <- do.call(c, seq_list)
  
  message(sprintf("Successfully extracted %d sequences", length(all_seqs)))
  
  # Write to file if requested
  if (!is.null(output_file)) {
    message("Writing all sequences to: ", output_file)
    writeXStringSet(all_seqs, filepath = output_file, format = "fasta")
  }
  
  return(all_seqs)
}


#' Extract FASTA for Features with Flanking Regions
#'
#' @param genome_obj Genome object from init_genome()
#' @param features GRanges object or gene name pattern
#' @param flank_upstream Upstream flanking size (default = 0)
#' @param flank_downstream Downstream flanking size (default = 0)
#' @param include_feature Include the feature itself (default = TRUE)
#' @param output_file Optional file path to write sequences (default = NULL)
#' @return DNAStringSet with extracted sequences
#' @export
#'
#' @examples
#' \dontrun{
#' # Get acrB genes with 5kb upstream and 2kb downstream
#' seqs <- get_feature_fasta(genome, "acrB", 
#'                          flank_upstream = 5000, 
#'                          flank_downstream = 2000)
#' 
#' # Get only upstream promoter region (exclude gene)
#' promoters <- get_feature_fasta(genome, "acrB",
#'                               flank_upstream = 1000,
#'                               flank_downstream = 0,
#'                               include_feature = FALSE)
#' }
get_feature_fasta <- function(genome_obj,
                              features,
                              flank_upstream = 0,
                              flank_downstream = 0,
                              include_feature = TRUE,
                              output_file = NULL) {
  
  # If features is a string, find matching features
  if (is.character(features)) {
    pattern <- features
    features <- genome_obj$gff[grepl(pattern, genome_obj$gff$Name)]
    
    if (length(features) == 0) {
      stop("No features found matching pattern: ", pattern)
    }
    
    message(sprintf("Found %d feature(s) matching '%s'", length(features), pattern))
  }
  
  # Fix seqnames if needed
  fa_headers <- genome_obj$seqnames
  if (length(fa_headers) == 1) {
    seqlevels(features) <- fa_headers
    seqnames(features) <- Rle(fa_headers)
  }
  
  # Prepare extraction regions based on strand
  seq_list <- lapply(seq_along(features), function(i) {
    feat <- features[i]
    
    # Calculate coordinates accounting for strand
    if (as.character(strand(feat)) == "-") {
      # For minus strand, upstream is actually downstream in coordinates
      new_start <- start(feat) - flank_downstream
      new_end <- end(feat) + flank_upstream
    } else {
      # For plus strand
      new_start <- start(feat) - flank_upstream
      new_end <- end(feat) + flank_downstream
    }
    
    # If not including feature, adjust coordinates
    if (!include_feature) {
      if (as.character(strand(feat)) == "-") {
        new_start <- end(feat) + 1
        new_end <- end(feat) + flank_upstream
      } else {
        new_start <- start(feat) - flank_upstream
        new_end <- start(feat) - 1
      }
    }
    
    # Generate descriptive name
    feat_name <- if (!is.null(feat$Name)) {
      sprintf("%s_%s:%d-%d_up%d_down%d(%s)",
              feat$Name,
              as.character(seqnames(feat)),
              start(feat),
              end(feat),
              flank_upstream,
              flank_downstream,
              as.character(strand(feat)))
    } else {
      sprintf("%s:%d-%d_up%d_down%d(%s)",
              as.character(seqnames(feat)),
              start(feat),
              end(feat),
              flank_upstream,
              flank_downstream,
              as.character(strand(feat)))
    }
    
    # Extract sequence
    seq <- get_roi_fasta(
      genome_obj = genome_obj,
      seqname = as.character(seqnames(feat)),
      start = max(1, new_start),
      end = new_end,
      strand = as.character(strand(feat)),
      seq_name = feat_name
    )
    
    return(seq)
  })
  
  # Combine all sequences
  all_seqs <- do.call(c, seq_list)
  
  # Write to file if requested
  if (!is.null(output_file)) {
    message("Writing all feature sequences to: ", output_file)
    writeXStringSet(all_seqs, filepath = output_file, format = "fasta")
  }
  
  return(all_seqs)
}


#' Quick wrapper: Extract ROI with Annotations
#'
#' @param genome_obj Genome object from init_genome()
#' @param seqname Chromosome/contig name
#' @param start Start coordinate
#' @param end End coordinate
#' @param flank Flanking region size (default = 0)
#' @param get_sequence Extract FASTA sequence (default = TRUE)
#' @param get_features Get overlapping features (default = TRUE)
#' @param output_fasta Optional FASTA output file
#' @return List with sequence (if requested) and features (if requested)
#' @export
extract_roi_complete <- function(genome_obj,
                                 seqname,
                                 start,
                                 end,
                                 flank = 0,
                                 get_sequence = TRUE,
                                 get_features = TRUE,
                                 output_fasta = NULL) {
  
  result <- list()
  
  # Extract sequence if requested
  if (get_sequence) {
    message("\n=== Extracting sequence ===")
    result$sequence <- get_roi_fasta(
      genome_obj = genome_obj,
      seqname = seqname,
      start = start,
      end = end,
      flank = flank,
      output_file = output_fasta
    )
  }
  
  # Get features if requested
  if (get_features) {
    message("\n=== Finding overlapping features ===")
    result$features <- analyze_roi(
      genome_obj = genome_obj,
      seqname = seqname,
      start = start,
      end = end,
      flank = flank
    )
  }
  
  message("\n=== Extraction complete ===")
  return(result)
}

#' Submit Protein BLAST Search
#'
#' @param sequence AAString or character vector with protein sequence
#' @param database NCBI database (default = "nr")
#' @param wait_for_results Logical; wait for results (default = TRUE)
#' @param max_wait Maximum wait time in seconds (default = 300)
#' @return List with RID and results (if wait_for_results = TRUE)
#' @export
blast_protein <- function(sequence, 
                          database = "nr", 
                          wait_for_results = TRUE,
                          max_wait = 300) {
  
  blast_url <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
  
  # Convert sequence to character if needed
  if (is(sequence, "AAString") || is(sequence, "DNAString")) {
    sequence <- as.character(sequence)
  }
  
  message("Submitting BLAST search...")
  
  # Submit BLAST job
  resp <- request(blast_url) |>
    req_body_form(
      CMD = "Put",
      PROGRAM = "blastp",
      DATABASE = database,
      QUERY = sequence,
      FORMAT_TYPE = "XML"
    ) |>
    req_perform()
  
  text_resp <- resp_body_string(resp)
  rid <- str_match(text_resp, "RID = (\\S+)")[, 2]
  
  if (is.na(rid)) {
    stop("Failed to submit BLAST job")
  }
  
  message("BLAST job submitted. RID: ", rid)
  
  if (!wait_for_results) {
    return(list(rid = rid, results = NULL))
  }
  
  # Poll for results
  message("Waiting for BLAST results...")
  start_time <- Sys.time()
  
  while (difftime(Sys.time(), start_time, units = "secs") < max_wait) {
    Sys.sleep(10)
    
    status_resp <- request(blast_url) |>
      req_url_query(CMD = "Get", RID = rid, FORMAT_TYPE = "XML") |>
      req_perform()
    
    status_text <- resp_body_string(status_resp)
    
    if (grepl("Status=READY", status_text)) {
      message("BLAST results ready!")
      results <- parse_blast_xml(status_text)
      return(list(rid = rid, results = results))
    }
    
    if (grepl("Status=FAILED", status_text)) {
      stop("BLAST job failed")
    }
    
    message("Still waiting... (", round(difftime(Sys.time(), start_time, units = "secs")), "s)")
  }
  
  warning("BLAST search timed out. Retrieve results later with RID: ", rid)
  return(list(rid = rid, results = NULL))
}


#' Parse BLAST XML Results
#'
#' @param xml_text XML text from BLAST results
#' @param top_n Number of top hits to return (default = 10)
#' @return Tibble with top BLAST hits
#' @export
parse_blast_xml <- function(xml_text, top_n = 10) {
  
  doc <- read_xml(xml_text)
  hits <- xml_find_all(doc, "//Hit")
  
  if (length(hits) == 0) {
    message("No BLAST hits found")
    return(tibble())
  }
  
  result <- tibble(
    accession = xml_text(xml_find_all(hits, ".//Hit_accession")),
    title = xml_text(xml_find_all(hits, ".//Hit_def")),
    length = as.integer(xml_text(xml_find_all(hits, ".//Hit_len")))
  )
  
  if (!is.null(top_n)) {
    result <- head(result, top_n)
  }
  
  return(result)
}



#' Flatten CharacterList to Character Vector
#'
#' @param x CharacterList object
#' @return Character vector with comma-separated values
flatten_CharacterList <- function(x) {
  sapply(x, function(v) {
    if (length(v)) paste(v, collapse = ",") else NA_character_
  })
}


#' Convert GRanges to Tidy Data Frame
#'
#' @param gr GRanges object
#' @param sort_by Sort by coordinate (default = TRUE)
#' @return Tibble with key annotation fields
#' @export
tidy_features <- function(gr, sort_by = "start") {
  
  if (length(gr) == 0) {
    return(tibble())
  }
  
  # Sort if requested
  if (!is.null(sort_by) && sort_by == "start") {
    gr <- gr[order(start(gr))]
  }
  
  # Convert to data frame and tidy up
  df <- as.data.frame(gr)
  
  # Flatten CharacterList columns if they exist
  if ("Alias" %in% names(df)) {
    df$Alias <- flatten_CharacterList(gr$Alias)
  }
  if ("Note" %in% names(df)) {
    df$Note <- flatten_CharacterList(gr$Note)
  }
  
  # Select key columns (adjust based on your GFF structure)
  key_cols <- c("seqnames", "start", "end", "width", "strand", "type", 
                "Name", "Alias", "Note")
  available_cols <- intersect(key_cols, names(df))
  
  result <- df %>%
    select(all_of(available_cols)) %>%
    as_tibble()
  
  return(result)
}


#' Compare Protein Sequences
#'
#' @param seq1 First sequence (AAString or character)
#' @param seq2 Second sequence (AAString or character)
#' @param type Alignment type (default = "global")
#' @return List with percent identity and alignment object
#' @export
compare_sequences <- function(seq1, seq2, type = "global") {
  
  if (!requireNamespace("pwalign", quietly = TRUE)) {
    stop("Package 'pwalign' is required. Install with: BiocManager::install('pwalign')")
  }
  
  aln <- pwalign::pairwiseAlignment(seq1, seq2, type = type)
  pid <- pwalign::pid(aln)
  
  message(sprintf("Percent identity: %.2f%%", pid))
  
  list(
    percent_identity = pid,
    alignment = aln
  )
}

#' Complete Gene Analysis Workflow
#'
#' @param genome_obj Genome object from init_genome()
#' @param gene_pattern Gene name pattern
#' @param flank_size Flanking region size for context
#' @param blast Perform BLAST search (default = FALSE)
#' @return List with sequences, context, and optional BLAST results
#' @export
analyze_gene <- function(genome_obj, 
                         gene_pattern, 
                         flank_size = 20000,
                         blast = FALSE) {
  
  message("=== Gene Analysis Pipeline ===\n")
  
  # 1. Extract sequences
  dna_seqs <- extract_sequences_by_name(genome_obj, gene_pattern, translate = FALSE)
  aa_seqs <- extract_sequences_by_name(genome_obj, gene_pattern, translate = TRUE)
  
  message("\nSequence lengths:")
  print(width(aa_seqs))
  
  # 2. Get genomic context
  message("\n--- Genomic Context ---")
  context <- get_genomic_context(genome_obj, gene_pattern, flank_size)
  
  # 3. Compare sequences if multiple found
  comparisons <- NULL
  if (length(aa_seqs) > 1) {
    message("\n--- Sequence Comparisons ---")
    comparisons <- list()
    for (i in 1:(length(aa_seqs) - 1)) {
      for (j in (i + 1):length(aa_seqs)) {
        comp_name <- paste(names(aa_seqs)[i], "vs", names(aa_seqs)[j])
        comparisons[[comp_name]] <- compare_sequences(aa_seqs[i], aa_seqs[j])
      }
    }
  }
  
  # 4. Optional BLAST
  blast_results <- NULL
  if (blast && length(aa_seqs) > 0) {
    message("\n--- BLAST Analysis ---")
    blast_results <- lapply(seq_along(aa_seqs), function(i) {
      message(sprintf("\nBLASTing %s...", names(aa_seqs)[i]))
      blast_protein(aa_seqs[[i]])
    })
    names(blast_results) <- names(aa_seqs)
  }
  
  # Return comprehensive results
  list(
    dna_sequences = dna_seqs,
    protein_sequences = aa_seqs,
    genomic_context = tidy_features(context),
    sequence_comparisons = comparisons,
    blast_results = blast_results
  )
}

' Run BLASTP with explicit binary and DB paths; hermetically scoped '
#' @param query_faa Path to query FASTA/FAA file
#' @param db BLAST database prefix (e.g., "swissprot" or absolute "/path/to/swissprot")
#' @param dbdir Directory containing BLAST DBs; if db is relative, dbdir is prepended
#' @param blastp_bin Full path to blastp binary (preferred explicit control)
#' @param blastp_dir Directory containing blastp; used if blastp_bin is NULL
#' @param evalue E-value threshold
#' @param threads Integer number of threads
#' @param fields Vector of outfmt 6 columns
#' @param validate_db Logical; check DB via blastdbcmd before running
#' @param more_args Character vector of additional BLASTP args (e.g., c("-max_target_seqs","10"))
#' @return tibble of hits; aborts with rich diagnostics on failure
blastp_capture <- function(
    query_faa,
    db = "swissprot",
    dbdir = Sys.getenv("BLASTDB", ""),
    blastp_bin = Sys.getenv("BLASTP_BIN", ""),
    blastp_dir = Sys.getenv("BLASTP_PATH", ""),
    evalue = 1e-10,
    threads = parallel::detectCores(),
    fields = c("qseqid","sacc","stitle","pident","length","qcovs","bitscore","evalue","staxids"),
    validate_db = TRUE,
    more_args = character()
) {
  # Dependencies: avoid attachment
  if (!requireNamespace("cli", quietly = TRUE)) stop("Package 'cli' is required.")
  if (!requireNamespace("readr", quietly = TRUE)) stop("Package 'readr' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  
  # --- Resolve blastp binary deterministically ---
  resolve_blastp <- function(bin, dir) {
    if (nzchar(bin)) {
      if (file.exists(bin) && isTRUE(file.info(bin)$isdir)) file.path(bin, "blastp") else bin
    } else if (nzchar(dir)) {
      if (file.exists(dir) && isTRUE(file.info(dir)$isdir)) file.path(dir, "blastp") else dir
    } else {
      Sys.which("blastp")
    }
  }
  candidate_bin <- resolve_blastp(blastp_bin, blastp_dir)
  
  if (!nzchar(candidate_bin) || !file.exists(candidate_bin)) {
    cli::cli_abort(c(
      "Could not locate {.val blastp}.",
      i = "Set {.env BLASTP_BIN} to the full binary or {.env BLASTP_PATH} to its directory, ",
      i = "or ensure 'blastp' is on PATH."
    ))
  }
  if (isTRUE(file.info(candidate_bin)$isdir)) {
    cli::cli_abort(c(
      "Provided path is a directory, not a binary: {.path {candidate_bin}}",
      i = "The wrapper accepts a directory, but it must contain a {.file blastp} executable."
    ))
  }
  if (.Platform$OS.type == "unix" && file.access(candidate_bin, 1) != 0) {
    cli::cli_abort("blastp binary is not executable: {.path {candidate_bin}}")
  }
  blastp <- candidate_bin
  
  # Resolve blastdbcmd alongside blastp; fallback to PATH
  blastdbcmd <- file.path(dirname(blastp), "blastdbcmd")
  if (!file.exists(blastdbcmd)) {
    blastdbcmd <- Sys.which("blastdbcmd")
  }
  
  # --- Preflight inputs ---
  if (!file.exists(query_faa)) {
    cli::cli_abort("Query file not found: {.path {query_faa}}")
  }
  
  # Derive db_path: if db is relative and dbdir provided, join; otherwise absolute db is used verbatim.
  db_path <- if (nzchar(dbdir) && !grepl("^/", db)) file.path(dbdir, db) else db
  
  # --- Optional DB sanity check (quote the path; out-of-band shell may split otherwise) ---
  if (validate_db && nzchar(blastdbcmd)) {
    db_path_q <- shQuote(db_path, type = "sh")
    info <- tryCatch(
      system2(blastdbcmd, args = c("-db", db_path_q, "-info"), stdout = TRUE, stderr = TRUE),
      error = function(e) character()
    )
    if (length(info) == 0 || any(grepl("\\bError\\b", info))) {
      cli::cli_abort(c(
        "Could not open BLAST DB at {.path {db_path}}.",
        ">" = paste(info, collapse = "\n"),
        i = "Provide an absolute DB prefix or set {.env BLASTDB} correctly."
      ))
    }
  } else if (validate_db) {
    cli::cli_warn("blastdbcmd not found; skipping DB validation.")
  }
  
  # --- Compose outfmt 6 (must be quoted if a shell is involved) ---
  fmt <- paste("6", paste(fields, collapse = " "))
  fmt_q <- shQuote(fmt, type = "sh")
  
  # --- Build BLASTP args ---
  # Quote any path-like arguments that may contain spaces if we fall back to system2.
  query_q <- shQuote(query_faa, type = "sh")
  db_q    <- shQuote(db_path,   type = "sh")
  
  args <- c(
    "-query", query_faa,            # raw for processx
    "-db",    db_path,              # raw for processx
    "-evalue", as.character(evalue),
    "-seg",   "yes",
    "-comp_based_stats", "2",
    "-num_threads", as.integer(threads),
    "-outfmt", fmt
  )
  args_sh <- c(                      # quoted for system2 with stderr/stdout capture
    "-query", query_q,
    "-db",    db_q,
    "-evalue", as.character(evalue),
    "-seg",   "yes",
    "-comp_based_stats", "2",
    "-num_threads", as.integer(threads),
    "-outfmt", fmt_q
  )
  if (length(more_args)) {
    # preserve extra args; quote any that have spaces to be safe under a shell
    more_args_sh <- ifelse(grepl("\\s", more_args), shQuote(more_args, type = "sh"), more_args)
    args    <- c(args,    more_args)
    args_sh <- c(args_sh, more_args_sh)
  }
  
  cli::cli_alert_info("Running: {blastp} -db {db_path} -query {query_faa} -outfmt {fmt}")
  
  # --- Hermetically set BLASTDB just for the call ---
  old_BLASTDB <- Sys.getenv("BLASTDB", unset = NA_character_)
  on.exit({
    if (!is.na(old_BLASTDB)) Sys.setenv(BLASTDB = old_BLASTDB) else Sys.unsetenv("BLASTDB")
  }, add = TRUE)
  if (nzchar(dbdir)) Sys.setenv(BLASTDB = dbdir)
  
  # --- Execute and capture (prefer processx if available to avoid shell quoting entirely) ---
  use_px <- requireNamespace("processx", quietly = TRUE)
  if (isTRUE(use_px)) {
    res <- processx::run(blastp, args = args, error_on_status = FALSE, echo_cmd = FALSE)
    status <- res$status
    out    <- res$stdout
    err    <- res$stderr
  } else {
    out <- suppressWarnings(system2(command = blastp, args = args_sh, stdout = TRUE, stderr = TRUE))
    status <- attr(out, "status")
    err    <- attr(out, "stderr")
  }
  
  if (!is.null(status) && status != 0) {
    msg <- if (!is.null(err) && length(err)) err else out
    cli::cli_abort(c(
      "blastp exited with status {status}.",
      ">" = paste(msg, collapse = "\n"),
      i = "Check binary path, DB prefix (quoted if it contains spaces), permissions, and device mounts."
    ))
  }
  
  # --- Parse and rank ---
  df <- readr::read_tsv(I(out), col_names = fields, show_col_types = FALSE)
  df |>
    dplyr::mutate(dplyr::across(c(pident, qcovs, bitscore, evalue), suppressWarnings(as.numeric))) |>
    dplyr::arrange(evalue, dplyr::desc(bitscore), dplyr::desc(qcovs))
}

#' Run BLASTP on one or more query FAA files and return tidy hits
#' Inherits binary/db controls from blastp_capture; no reticulate involved.
blastp_roi <- function(
    faa_path,
    db = "swissprot",
    dbdir = Sys.getenv("BLASTDB", ""),
    blastp_bin = Sys.getenv("BLASTP_BIN", ""),
    blastp_dir = Sys.getenv("BLASTP_PATH", ""),
    evalue = 1e-10,
    threads = parallel::detectCores(),
    fields = c("qseqid","sacc","stitle","pident","length","qcovs","bitscore","evalue","staxids"),
    validate_db = TRUE,
    more_args = character()
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("cli", quietly = TRUE)) stop("Package 'cli' is required.")
  
  # Allow vector input; run sequentially unless you later want parallelization
  faa_path <- unique(faa_path)
  res_list <- vector("list", length(faa_path))
  
  for (i in seq_along(faa_path)) {
    q <- faa_path[i]
    cli::cli_alert_info("BLASTP ROI: {i}/{length(faa_path)} — {.path {q}}")
    
    res_list[[i]] <- blastp_capture(
      query_faa   = q,
      db          = db,
      dbdir       = dbdir,
      blastp_bin  = blastp_bin,
      blastp_dir  = blastp_dir,
      evalue      = evalue,
      threads     = threads,
      fields      = fields,
      validate_db = validate_db,
      more_args   = more_args
    )
  }
  
  # Bind rows; ensure required columns exist
  out <- dplyr::bind_rows(res_list)
  needed <- c("qseqid","sacc","stitle","pident","length","qcovs","bitscore","evalue")
  missing <- setdiff(needed, colnames(out))
  if (length(missing)) {
    stop("Missing columns in BLAST output: ", paste(missing, collapse = ", "))
  }
  
  out
}

#' Reduce BLAST hits with cov/identity thresholds; select best hit or top-N per query
reduce_hits <- function(
    hits,
    min_qcov = 40,
    min_pident = 25,
    besthit = TRUE,
    max_per_query = NULL
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!is.data.frame(hits) || nrow(hits) == 0) return(hits)
  
  h <- hits
  
  # Coerce numerics once, defensively; preserve NA semantics
  num_cols <- intersect(c("pident","qcovs","bitscore","evalue"), colnames(h))
  h <- dplyr::mutate(h, dplyr::across(dplyr::all_of(num_cols), suppressWarnings(as.numeric)))
  
  # Thresholds: allow NA to pass (as you did), which is reasonable for partial outputs
  if (!is.null(min_qcov) && "qcovs" %in% colnames(h)) {
    h <- dplyr::filter(h, is.na(qcovs) | qcovs >= min_qcov)
  }
  if (!is.null(min_pident) && "pident" %in% colnames(h)) {
    h <- dplyr::filter(h, is.na(pident) | pident >= min_pident)
  }
  
  h <- dplyr::group_by(h, qseqid)
  
  if (isTRUE(besthit)) {
    # Three-step tie-break: evalue → bitscore → qcovs
    h <- h %>%
      dplyr::slice_min(order_by = evalue, n = 1, with_ties = TRUE) %>%
      dplyr::slice_max(order_by = bitscore, n = 1, with_ties = TRUE) %>%
      dplyr::slice_max(order_by = qcovs,   n = 1, with_ties = FALSE)
  } else {
    h <- h %>%
      dplyr::arrange(evalue, dplyr::desc(bitscore), dplyr::desc(qcovs))
    if (!is.null(max_per_query) && is.finite(max_per_query)) {
      h <- dplyr::slice_head(h, n = max_per_query)
    }
  }
  
  dplyr::ungroup(h)
}

translate_cds <- function(dna, strand = "+", drop_terminal_stop = TRUE,
                          bact_gc = Biostrings::getGeneticCode("11")) {
  
  # specified bacterial/archaeal genetic code (11) 
  
  # Normalize orientation
  dna2 <- if (as.character(strand) == "-") Biostrings::reverseComplement(dna) else dna
  
  # Trim to nearest full codon to avoid fuzzy tail at ROI edges
  rem <- width(dna2) %% 3L
  if (rem != 0L) dna2 <- dna2[1:(width(dna2) - rem)]
  
  aa <- Biostrings::translate(dna2, genetic.code = bact_gc, if.fuzzy.codon = "X")
  aa <- as.character(aa)
  
  # Strip a terminal stop, keep internal stops as QC
  if (drop_terminal_stop) aa <- sub("\\*$", "", aa)
  aa
}

roi_cds_to_faa <- function(genome_obj, roi_tbl, 
                           bact_gc = Biostrings::getGeneticCode("11"),
                           outfile = tempfile(fileext = ".faa")) {
  
  # process data
  cds <- roi_tbl %>%
    filter(type == "CDS") %>%
    mutate(
      dna = pmap(list(seqnames, start, end, strand),
                 ~ extract_sequence_by_coords(genome_obj, ..1, ..2, ..3, ..4)),
      aa_raw = map2_chr(dna, strand, ~ as.character(Biostrings::translate(
        if (as.character(..2) == "-") Biostrings::reverseComplement(..1) else ..1,
        genetic.code = bact_gc, if.fuzzy.codon = "X"
      ))),
      aa = sub("\\*$", "", aa_raw),  # drop terminal stop only
      aa_len = nchar(aa),
      has_stop = str_detect(aa_raw, "\\*") & !str_detect(aa_raw, "\\*$")
    )
  
  # Header preference: Name, else Alias, else coord signature
  hdr <- cds$Name
  hdr[is.na(hdr) | hdr == ""] <- cds$Alias[is.na(hdr) | hdr == ""]
  fallback <- paste0(cds$seqnames, ":", cds$start, "-", cds$end, "(", cds$strand, ")")
  hdr[is.na(hdr) | hdr == ""] <- fallback[is.na(hdr) | hdr == ""]
  
  # Ensure uniqueness if the GFF repeats IDs
  if (any(duplicated(hdr))) {
    dup_idx <- which(duplicated(hdr) | duplicated(hdr, fromLast = TRUE))
    hdr[dup_idx] <- paste0(hdr[dup_idx], "_", seq_along(dup_idx))
  }
  
  lines <- unlist(map2(hdr, cds$aa, ~ c(paste0(">", .x), .y)))
  write_lines(lines, outfile)
  
  cds %>% mutate(qseqid = hdr, faa_path = outfile)
}

# Using the functions ----

setwd("/home/william-ackerman/Desktop/cnMOPS")

NIST_genome <- init_genome("reference.gff3", "reference.fasta")

extract_sequences_by_name(NIST_genome, "acrB", translate = TRUE)
extract_sequences_by_name(NIST_genome, "^bla", translate = TRUE)
extract_sequences_by_name(NIST_genome, "^fl", translate = TRUE)

NIST_genome$gff %>% as.data.frame() %>% filter(type == "CDS") %>% 
  filter(str_detect(Name, "^acr"))

NIST_genome$gff %>% as.data.frame() %>% filter(type == "CDS") %>% 
  filter(str_detect(Name, "^flh"))

NIST_genome$gff %>% as.data.frame() %>% filter(type == "CDS") %>% 
  filter(str_detect(Name, "^mot"))

NIST_genome$gff %>% as.data.frame() %>% filter(type == "CDS") %>% 
  filter(str_detect(Name, "^fli"))

NIST_genome$gff %>% as.data.frame() %>% filter(type == "CDS") %>% 
  filter(str_detect(Name, "^laf")) %>% select(Note)

NIST_genome$gff %>% as.data.frame() %>% filter(type == "CDS") %>% 
  filter(str_detect(Name, "^pur")) 

NIST_genome$gff %>% as.data.frame() %>% filter(type == "CDS") %>% 
  filter(str_detect(Name, "^pyr")) 

NIST_genome$gff %>% as.data.frame() %>% filter(type == "tRNA") %>% 
  filter(str_detect(Name, "^arg")) 

analyze_roi(NIST_genome, "1", 2369501, 2613000, flank = 10000) %>% 
  filter(str_detect(Name, "^flh"))


analyze_roi(NIST_genome, "1", 1706848, 1759613, flank = 5000) %>% 
  print(n=Inf) 

analyze_roi(NIST_genome, "1", 1706848, 1759613, flank = 50) %>% 
  print(n=Inf) -> prophage_tbl

save(prophage_tbl, file = "Prophage_tbl.RData")


prophage_cds <- roi_cds_to_faa(NIST_genome, prophage_tbl)


prophage_cds %>% print(n = Inf)

res <- blastp_capture(
  query_faa  = unique(prophage_cds$faa_path),
  db         = "swissprot",
  dbdir      = "/home/william-ackerman/Desktop/tmp_ncbi_blast",  # symlink target, no spaces
  blastp_dir = "/home/william-ackerman/Desktop/tmp_ncbi_blast/ncbi-blast-2.17.0+/bin",
  threads    = 32
)

res %>% 
  dplyr::group_by(qseqid) %>%
  dplyr::slice_min(order_by = evalue, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()


# Run on the ROI FASTA

hits   <- blastp_roi(unique(prophage_cds$faa_path), db = "swissprot",
                     dbdir      = "/home/william-ackerman/Desktop/tmp_ncbi_blast",  # symlink target, no spaces
                     blastp_dir = "/home/william-ackerman/Desktop/tmp_ncbi_blast/ncbi-blast-2.17.0+/bin")
best   <- reduce_hits(hits, min_qcov = 40, min_pident = 25, besthit = TRUE)

hits %>% print(n = Inf)

hits %>%
  dplyr::group_by(qseqid) %>%
  dplyr::slice_min(order_by = evalue, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()



setwd("/home/william-ackerman/Desktop/cnMOPS/Ecoli_MG1655")

MG1655_genome <- init_genome("reference.gff3", "reference.fasta")


NIST_genes <- 
  NIST_genome$gff %>% as.data.frame() %>% filter(seqnames == 1) %>% 
  filter(type == "CDS") %>% filter(!str_detect(Name, "JNMLBEGN_")) %>% 
  select(Name)

MG1655_genes <- 
  MG1655_genome$gff %>% as.data.frame() %>% filter(seqnames == "U00096") %>% 
  filter(type == "CDS") %>% filter(!str_detect(Name, "JNMLBEGN_")) %>% 
  select(Name)

NIST_short <- str_remove_all(NIST_genes$Name, "_1|_2|_3|_4|_5|_6") %>% unique()

head(NIST_genes)
head(MG1655_genes)

unique_NIST_genes <- setdiff(NIST_short, MG1655_genes$Name) %>% as.data.frame()

write.csv(unique_NIST_genes, file = "Genes_Only_in_NIST_not_MG1655.csv")

search_features(MG1655_genome , "transposase", field = "Note") %>% 
  as.data.frame()

search_features(MG1655_genome , "acr", field = "Name") %>% 
  as.data.frame()

analyze_roi(MG1655_genome, "U00096", 481254, 486408, flank = 10000) %>% 
  as.data.frame()

options(max.print = 100000) 
analyze_roi(NIST_genome, "1", 2369501, 2613000) %>% 
  as.data.frame() %>% write.csv(., file = "Second_region_genes.csv")


analyze_roi(NIST_genome, "1", 1828000, 1978000) %>% 
  as.data.frame() %>% write.csv(., file = "First_region_genes.csv")

# analyze_gene(NIST_genome, "acrB_1", blast = TRUE) # NEEDS WORK



analyze_roi(NIST_genome, "1", 1827608 , 1828588, flank = 0)

analyze_roi(NIST_genome, "1", 2369501, 2613000, flank = 0) %>% 
  arrange(start) %>% 
  as.data.frame() %>% filter(str_detect(Name, "amp"))

options(max.print = 3000)

IS621_candidate_1 <- 
  get_roi_fasta(
    genome_obj = NIST_genome,
    seqname = "1",
    start = 1827608,
    end = 1828588)

as.character(IS621_candidate_1) # confirmed US621 

roi_seq <- get_roi_fasta(
  genome_obj = NIST_genome,
  seqname = "1",
  start = 1,
  end = 4793125) %>% as.character()

fasta_path <- "NIST0056_chr1.fasta"
header     <- ">NIST0056_chr1"


wrap_seq <- function(x, width = 80) {
  n <- nchar(x)
  starts <- seq(1, n, by = width)
  vapply(starts, function(i) substr(x, i, min(i + width - 1, n)), "", USE.NAMES = FALSE)
}

writeLines(c(">NIST0056_chr1", wrap_seq(roi_seq, width = 80)), "NIST0056_chr1.fasta")
  
  
# makeblastdb -in NIST0056_chr1.fasta -dbtype nucl -out NIST0056_genome_db
# blastn -query IS26.fasta -db NIST0056_genome_db -out IS26_hits.txt -evalue 1e-10 -outfmt 6


# This is the region in A2-10 and C1-10 that had reduced copy numbe
analyze_roi(NIST_genome, "1", 1725534, 1759705, flank = 0) %>% 
  print(n=Inf)

extract_sequence_by_coords(NIST_genome, "1", 1725535, 1726107, strand = "-") %>% 
  reverse() %>% translate() %>% as.character()

analyze_roi(NIST_genome, "1", 1725534, 1759705, flank = 50000) 

extract_sequence_by_coords(NIST_genome, "1", 1725535 - 50000, 1726107 + 50000) 

