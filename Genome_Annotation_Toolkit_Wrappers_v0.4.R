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
library(GenomeInfoDb)
library(S4Vectors)

# Setup and Data Loading Functions ====

# A tiny infix helper to mirror rlang's %||%, avoiding hard dependency.
`%||%` <- function(x, y) if (is.null(x)) y else x
# Helpers

# Scrubber for common prefixes. Extend as your corpus dictates.
.scrub_prefixes <- function(x) {
  x <- as.character(x)
  x <- sub("^lcl\\|", "", x, perl = TRUE)
  x <- sub("^chr",    "", x, ignore.case = TRUE)
  x
}


# Return FASTA names and lengths from either DNAStringSet or FaFile
.get_fasta_index <- function(fa_obj) {
  if (inherits(fa_obj, "DNAStringSet")) {
    nm <- names(fa_obj) %||% character(length(fa_obj))
    len <- Biostrings::width(fa_obj)
    return(list(names = nm, lengths = setNames(as.integer(len), nm)))
  }
  if (inherits(fa_obj, "FaFile")) {
    if (!requireNamespace("Rsamtools", quietly = TRUE)) {
      stop("Rsamtools is required to inspect FaFile indices.")
    }
    idx <- Rsamtools::scanFaIndex(fa_obj)
    nm <- names(idx)
    len <- vapply(idx, function(rec) rec$seqlength, integer(1))
    return(list(names = nm, lengths = setNames(as.integer(len), nm)))
  }
  # Fallback if you only have headers
  list(names = character(), lengths = integer())
}

# Pre-filter the GFF to drop rows with missing start/end
clean_gff_for_import <- function(gff_path, 
                                 drop_invalid = TRUE, 
                                 verbose = TRUE) {
  if (!file.exists(gff_path)) stop("File not found: ", gff_path)
  
  msg <- function(...) if (verbose) message(...)
  
  msg("Reading GFF lines...")
  lines <- readLines(gff_path, warn = FALSE)
  
  # Remove comment/directive lines and blank lines
  is_comment <- grepl("^\\s*#", lines)
  is_blank   <- !nzchar(trimws(lines))
  lines <- lines[!(is_comment | is_blank)]
  
  if (length(lines) == 0) stop("No feature lines found after removing comments: ", gff_path)
  
  # Split into fields by tab
  parts <- strsplit(lines, "\t", fixed = TRUE)
  lens  <- lengths(parts)
  
  # Keep rows with at least 9 fields; truncate extras
  keep_len <- lens >= 9
  dropped_len <- sum(!keep_len)
  
  if (dropped_len > 0 && verbose) {
    msg("Dropping ", dropped_len, " lines with fewer than 9 tab-separated fields.")
  }
  
  parts <- parts[keep_len]
  # Build a 9-column matrix
  mat <- do.call(rbind, lapply(parts, function(x) x[1:9]))
  df  <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(df) <- c("seqid","source","type","start","end","score","strand","phase","attributes")
  
  # Coerce start/end to integers; treat "." as NA
  to_int <- function(x) {
    x[x %in% c("", ".", "NA", "NaN")] <- NA_character_
    suppressWarnings(as.integer(x))
  }
  df$start <- to_int(df$start)
  df$end   <- to_int(df$end)
  
  # Drop invalid rows if requested
  keep <- rep(TRUE, nrow(df))
  invalid_coord <- is.na(df$start) | is.na(df$end)
  inverted      <- !invalid_coord & (df$start > df$end)
  
  if (drop_invalid) {
    dropped <- sum(invalid_coord | inverted)
    if (dropped > 0 && verbose) {
      msg("Dropping ", dropped, " rows with missing or inverted coordinates.")
    }
    df <- df[!(invalid_coord | inverted), , drop = FALSE]
  } else if (any(invalid_coord | inverted)) {
    warning("There are rows with missing or inverted coordinates; rtracklayer::import() may fail.")
  }
  
  if (nrow(df) == 0) stop("No valid rows remain after filtering.")
  
  # Write to a temporary GFF3 with the required directive
  tmp <- tempfile(fileext = ".gff3")
  con <- file(tmp, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("##gff-version 3", con)
  
  # Write rows with tab separation, no quotes, no headers
  write.table(
    df, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  
  if (verbose) msg("Clean GFF written to: ", tmp)
  tmp
}


#' Harmonize GFF seqlevels to FASTA headers
#'
#' Attempts, seriatim:
#' 1) Identity: if all GFF seqlevels already match FASTA headers.
#' 2) Positional numeric mapping: if GFF levels look like "1","3","4",... and counts match FASTA.
#' 3) Region-guided mapping: infer mapping from top-level 'region' features, preserving their order.
#' 4) (Optional) Length-guided mapping: match region widths to FASTA contig lengths greedily.
#'
#' Returns a modified genome_obj with gff seqlevels renamed to match genome_obj$seqnames.
#' Emits informative cli messages; uses cli::warn on partial successes; cli::abort only when hopeless.


harmonize_gff_seqlevels <- function(genome_obj,
                                    use_length_fallback = FALSE) {
  
  
  # --- Required namespaces (fail fast, clear) ---
  stopifnot(requireNamespace("cli", quietly = TRUE))
  stopifnot(requireNamespace("rtracklayer", quietly = TRUE))
  stopifnot(requireNamespace("Biostrings", quietly = TRUE))
  stopifnot(requireNamespace("Rsamtools", quietly = TRUE))
  stopifnot(requireNamespace("GenomeInfoDb", quietly = TRUE))
  stopifnot(requireNamespace("GenomicRanges", quietly = TRUE))
  stopifnot(requireNamespace("S4Vectors", quietly = TRUE))
  
  # --- Basic object checks ---
  stopifnot(
    is.list(genome_obj),
    all(c("gff","fa","seqnames") %in% names(genome_obj))
  )
  
  
  gff <- genome_obj$gff
  fa_headers <- genome_obj$seqnames
  
  # ---- Compatibility shims ---------------------------------------------------
  .set_seqlevels <- function(gr, keep, pruning = "coarse") {
    # Try modern form with pruning.mode, else fall back
    ok <- TRUE
    res <- try({
      GenomeInfoDb::seqlevels(gr, pruning.mode = pruning) <- keep
      gr
    }, silent = TRUE)
    if (inherits(res, "try-error")) {
      ok <- FALSE
    } else {
      return(res)
    }
    # Fallback (older GenomeInfoDb)
    GenomeInfoDb::seqlevels(gr) <- keep
    gr
  }
  
  .rename_seqlevels <- function(gr, map_named_vector, pruning = "coarse") {
    # Try modern form with pruning.mode, else fall back
    res <- try({
      GenomeInfoDb::renameSeqlevels(gr, map_named_vector, pruning.mode = pruning)
    }, silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    # Fallback (older GenomeInfoDb)
    GenomeInfoDb::renameSeqlevels(gr, map_named_vector)
  }
  # ---------------------------------------------------------------------------
  
  cli::cli_h2("Reconciling GFF seqlevels with FASTA headers")
  cli::cli_li("FASTA headers: {.val {fa_headers}}")
  cli::cli_li("GFF seqlevels: {.val {GenomeInfoDb::seqlevels(gff)}}")
  
  gff_levels <- GenomeInfoDb::seqlevels(gff)
  
  # 1) Identity
  if (all(gff_levels %in% fa_headers)) {
    keep <- intersect(fa_headers, gff_levels)
    gff <- .set_seqlevels(gff, keep)
    cli::cli_alert_success("GFF is already harmonized with FASTA. Using existing labels.")
    genome_obj$gff <- gff
    return(genome_obj)
  }
  
  # helper for renaming + validation
  .try_map <- function(map_named_vector, label = NULL) {
    gff2 <- .rename_seqlevels(gff, map_named_vector)
    if (all(GenomeInfoDb::seqlevels(gff2) %in% fa_headers)) {
      if (!is.null(label)) {
        cli::cli_li("Applied mapping ({.emph {label}}): {.field {paste(names(map_named_vector), '->', unname(map_named_vector), collapse=', ')}}")
      } else {
        cli::cli_li("Applied mapping: {.field {paste(names(map_named_vector), '->', unname(map_named_vector), collapse=', ')}}")
      }
      return(gff2)
    }
    NULL
  }
  
  # 2) Positional numeric mapping when counts match
  looks_numeric <- all(grepl("^[0-9]+$", gff_levels))
  if (looks_numeric && length(gff_levels) == length(fa_headers)) {
    cli::cli_alert_info(
      "Attempting positional numeric mapping: {.val {gff_levels}} -> {.val {fa_headers}}"
    )
    # Map each *existing* GFF level to the FASTA header at the same position
    # e.g., c("1","3","4","5","7") -> c("Chromosome","Plasmid_1_(...)","Plasmid_2_(...)","Plasmid_3_(...)","Plasmid_4_(...)")
    map <- stats::setNames(fa_headers, gff_levels)
    gff2 <- .try_map(map, label = "positional")
    if (!is.null(gff2)) {
      cli::cli_alert_success("Positional numeric mapping succeeded.")
      genome_obj$gff <- gff2
      return(genome_obj)
    } else {
      cli::cli_warn("Positional numeric mapping failed. Will try region-guided mapping next.")
    }
  } else {
    cli::cli_alert_info("Skipping positional numeric mapping (levels not purely numeric or counts differ).")
  }
  
  # 3) Region-guided mapping (order of region rows -> order of FASTA headers)
  has_type <- "type" %in% names(S4Vectors::mcols(gff))
  if (has_type) {
    is_region <- S4Vectors::mcols(gff)$type == "region"
    if (any(is_region)) {
      regions <- gff[is_region]
      reg_levels <- unique(as.character(GenomicRanges::seqnames(regions)))
      cli::cli_alert_info("Found {.strong region} rows. Region seqnames: {.val {reg_levels}}")
      
      if (length(reg_levels) == length(fa_headers)) {
        map <- stats::setNames(fa_headers, reg_levels)  # e.g., "1"->"Chromosome"
        cli::cli_alert_info("Attempting region-guided mapping.")
        gff2 <- .try_map(map, label = "region-guided")
        if (!is.null(gff2)) {
          cli::cli_alert_success("Region-guided mapping succeeded.")
          genome_obj$gff <- gff2
          return(genome_obj)
        } else {
          cli::cli_warn("Region-guided mapping did not converge.")
        }
      } else {
        cli::cli_warn("Region count does not match FASTA count. Expected {.val {length(fa_headers)}}, got {.val {length(reg_levels)}}.")
      }
    } else {
      cli::cli_alert_info("No region features available in GFF; skipping region-guided mapping.")
    }
  } else {
    cli::cli_alert_info("No 'type' column present in GFF metadata; cannot use region-guided mapping.")
  }
  
  # 4) Optional length-guided greedy matching
  if (isTRUE(use_length_fallback)) {
    cli::cli_h3("Length-guided fallback mapping")
    
    # FASTA lengths
    if (!is.null(genome_obj$fasta)) {
      fa_lengths <- Biostrings::width(genome_obj$fasta)
      names(fa_lengths) <- genome_obj$seqnames
    } else {
      fai <- Rsamtools::scanFaIndex(genome_obj$fa)
      fa_lengths <- fai$seqlengths
      names(fa_lengths) <- names(fai)
    }
    
    if (has_type && any(S4Vectors::mcols(gff)$type == "region")) {
      regions <- gff[S4Vectors::mcols(gff)$type == "region"]
      gff_lengths <- GenomicRanges::width(regions)
      names(gff_lengths) <- as.character(GenomicRanges::seqnames(regions))
      
      fa_left <- fa_lengths
      map <- character(length(gff_lengths))
      names(map) <- names(gff_lengths)
      for (k in names(gff_lengths)) {
        best <- names(which.min(abs(fa_left - gff_lengths[[k]])))
        map[[k]] <- best
        fa_left <- fa_left[names(fa_left) != best]
      }
      gff2 <- .try_map(map, label = "length-guided")
      if (!is.null(gff2)) {
        cli::cli_alert_success("Length-guided mapping succeeded.")
        genome_obj$gff <- gff2
        return(genome_obj)
      } else {
        cli::cli_warn("Length-guided mapping failed to produce a consistent vocabulary.")
      }
    } else {
      cli::cli_warn("Cannot compute region lengths without region features; skipping length fallback.")
    }
  }
  
  cli::cli_abort(c(
    "x" = "Unable to harmonize GFF seqlevels to FASTA headers.",
    "i" = paste0("GFF levels: ", paste(gff_levels, collapse = ", ")),
    "i" = paste0("FASTA: ", paste(fa_headers, collapse = ", ")),
    "i" = "Inspect 'region' rows and FASTA order, or enable use_length_fallback=TRUE."
  ))
}

#' Initialize Genome Resources
#'
#' @param gff_path Path to GFF3 file
#' @param fasta_path Path to FASTA file
#' @param auto_index Logical; automatically index FASTA if needed
#' @return List with gff (GRanges), fasta (DNAStringSet), and fa (FaFile) objects
#' @export
#' 
#' 
init_genome <- function(
    gff_path,
    fasta_path,
    auto_index = TRUE,
    verbose = TRUE
) {
  # Dependencies
  stopifnot(requireNamespace("cli", quietly = TRUE))
  stopifnot(requireNamespace("rtracklayer", quietly = TRUE))
  stopifnot(requireNamespace("Biostrings", quietly = TRUE))
  stopifnot(requireNamespace("Rsamtools", quietly = TRUE))
  stopifnot(requireNamespace("GenomeInfoDb", quietly = TRUE))
  stopifnot(requireNamespace("GenomicRanges", quietly = TRUE))
  
  inform <- function(...) if (isTRUE(verbose)) cli::cli_inform(...)
  warn   <- function(...) cli::cli_warn(...)
  abort  <- function(...) cli::cli_abort(...)
  
  if (!file.exists(gff_path))   abort("GFF path not found: {gff_path}")
  if (!file.exists(fasta_path)) abort("FASTA path not found: {fasta_path}")
  
  # 1) Direct GFF import attempt
  inform("Loading GFF3 annotation (direct import) ...")
  gff_used <- gff_path
  import_errs <- list(direct = NULL, cleaned = NULL)
  
  gff <- try(rtracklayer::import(gff_path, format = "gff3"), silent = TRUE)
  if (inherits(gff, "try-error")) {
    import_errs$direct <- attr(gff, "condition")
    inform(c(
      "x" = "Direct GFF import failed.",
      ">" = conditionMessage(import_errs$direct)
    ))
    
    # 2) Clean and retry
    inform("Attempting pre-filter of GFF and re-import ...")
    tmp_gff <- try(clean_gff_for_import(gff_path, drop_invalid = TRUE, verbose = verbose),
                   silent = TRUE)
    if (inherits(tmp_gff, "try-error")) {
      import_errs$cleaned <- attr(tmp_gff, "condition")
      abort(c(
        "Unable to import GFF and cleaning also failed.",
        "• Direct import error:"  = conditionMessage(import_errs$direct),
        "• Cleaning failure:"     = conditionMessage(import_errs$cleaned)
      ))
    }
    
    gff_used <- tmp_gff
    gff2 <- try(rtracklayer::import(tmp_gff, format = "gff3"), silent = TRUE)
    if (inherits(gff2, "try-error")) {
      import_errs$cleaned <- attr(gff2, "condition")
      abort(c(
        "Unable to import GFF even after cleaning.",
        "• Direct import error:"  = conditionMessage(import_errs$direct),
        "• Cleaned import error:" = conditionMessage(import_errs$cleaned),
        "!" = "Consider inspecting the first 50 non-comment lines for malformed records."
      ))
    }
    gff <- gff2
  }
  
  # 3) FASTA load and index
  inform("Loading FASTA sequence ...")
  fasta <- Biostrings::readDNAStringSet(fasta_path)
  
  if (isTRUE(auto_index) && !file.exists(paste0(fasta_path, ".fai"))) {
    inform("Indexing FASTA (.fai not found) ...")
    Rsamtools::indexFa(fasta_path)
  }
  
  fa <- Rsamtools::FaFile(fasta_path)
  Rsamtools::open.FaFile(fa)
  
  # 4) Seqnames via FaFile index if possible
  seqs <- try(
    GenomicRanges::seqnames(GenomeInfoDb::seqinfo(fa)),
    silent = TRUE
  )
  if (inherits(seqs, "try-error")) {
    inform("Falling back to sequence names from FASTA object.")
    seqs <- names(fasta)
  }
  
  inform("Genome resources loaded successfully.")
  
  list(
    gff         = gff,
    fasta       = fasta,
    fa          = fa,
    seqnames    = as.character(seqs),
    gff_used    = gff_used,
    import_errs = import_errs
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
#' 
#' troublehooting
# genome_obj <- NIST_pgap_genome
# gene_pattern <- "acrB"
extract_sequences_by_name <- function(genome_obj, 
                                      gene_pattern, 
                                      translate = FALSE,
                                      genetic_code = "11",
                                      auto_harmonize = TRUE) {
  stopifnot(requireNamespace("cli", quietly = TRUE))
  stopifnot(requireNamespace("Biostrings", quietly = TRUE))
  stopifnot(requireNamespace("GenomeInfoDb", quietly = TRUE))
  stopifnot(requireNamespace("GenomicRanges", quietly = TRUE))
  stopifnot(requireNamespace("S4Vectors", quietly = TRUE))
  stopifnot(requireNamespace("Rsamtools", quietly = TRUE))
  
  gff <- genome_obj$gff
  fa_headers <- genome_obj$seqnames
  
  # ---------- helpers ----------
  .matches <- function(x, pat) {
    if (is.null(x)) return(rep(FALSE, length(GenomicRanges::seqnames(gff))))
    y <- tryCatch(as.character(x), error = function(e) rep(NA_character_, length(x)))
    !is.na(y) & grepl(pat, y, perl = TRUE)
  }
  .prefer_name_then_fallback <- function(gff, pat) {
    name_idx <- .matches(gff$Name, pat)
    if (any(name_idx, na.rm = TRUE)) return(gff[name_idx])
    # fallback only if Name yields nothing
    idx <- .matches(gff$gene, pat) |
      .matches(gff$locus_tag, pat) |
      .matches(gff$product, pat)
    idx[is.na(idx)] <- FALSE
    gff[idx]
  }
  .drop_extdb_dupes <- function(gr) {
    if (length(gr) == 0L) return(gr)
    is_extdb <- !is.na(gr$Name) & grepl("^extdb:", gr$Name)
    if (!any(is_extdb)) return(gr)
    # if a non-extdb feature exists with the same locus_tag (or gene), drop the extdb one
    key <- if ("locus_tag" %in% names(S4Vectors::mcols(gr))) gr$locus_tag else gr$gene
    if (is.null(key)) key <- gr$Name
    keep <- rep(TRUE, length(gr))
    if (!all(is.na(key))) {
      keys_ext <- key[is_extdb]
      keys_good <- key[!is_extdb]
      drop_ext <- is_extdb & key %in% keys_good
      keep <- !drop_ext
    } else {
      keep <- !is_extdb  # no keys to reconcile; just drop extdb
    }
    gr[keep]
  }
  .best_labels <- function(gr) {
    n <- length(gr)
    out <- rep(NA_character_, n)
    # first non-NA among gene, locus_tag, Name, product
    cand <- list(gr$gene, gr$locus_tag, gr$Name, gr$product)
    for (v in cand) {
      v <- tryCatch(as.character(v), error = function(e) rep(NA_character_, n))
      fill <- is.na(out) & !is.na(v)
      out[fill] <- v[fill]
    }
    # avoid pure extdb names if a better label exists; if all NA, synthesize
    all_extdb <- !is.na(out) & grepl("^extdb:", out)
    if (any(all_extdb)) {
      # try to replace with locus_tag or gene if present
      repl <- if ("locus_tag" %in% names(S4Vectors::mcols(gr))) gr$locus_tag else gr$gene
      repl <- tryCatch(as.character(repl), error = function(e) rep(NA_character_, n))
      swap <- all_extdb & !is.na(repl)
      out[swap] <- repl[swap]
    }
    out[is.na(out)] <- paste0("feature_", seq_len(n))[is.na(out)]
    make.unique(out, sep = "_")
  }
  # -
  
  # Select features honoring your original "Name-first" intent
  features <- .prefer_name_then_fallback(gff, gene_pattern)
  if (length(features) == 0L) {
    cli::cli_warn("No features found matching pattern: {.val {gene_pattern}}")
    return(NULL)
  }
  
  # Filter feature types depending on translation
  if (isTRUE(translate)) {
    # Prefer CDS for protein-coding extraction
    if ("type" %in% names(S4Vectors::mcols(features))) {
      cds_ok <- S4Vectors::mcols(features)$type %in% c("CDS")
      if (any(cds_ok)) features <- features[cds_ok]
    }
  } else {
    # Prefer gene features for DNA; fall back to CDS if no gene present
    if ("type" %in% names(S4Vectors::mcols(features))) {
      gene_ok <- S4Vectors::mcols(features)$type %in% c("gene")
      if (any(gene_ok)) {
        features <- features[gene_ok]
      } else {
        cds_ok <- S4Vectors::mcols(features)$type %in% c("CDS")
        if (any(cds_ok)) features <- features[cds_ok]
      }
    }
  }
  
  # Drop extdb duplicates when a better labeled counterpart exists
  features <- .drop_extdb_dupes(features)
  
  cli::cli_alert_success("Selected {length(features)} feature{?s} for extraction.")
  
  # Vocabulary check and auto-harmonize if needed
  fa_headers <- genome_obj$seqnames
  current_vocab_ok <- all(as.character(GenomicRanges::seqnames(features)) %in% fa_headers)
  if (!current_vocab_ok && isTRUE(auto_harmonize)) {
    cli::cli_alert_info("Seqname vocabulary mismatch detected. Invoking harmonizer...")
    genome_obj <- harmonize_gff_seqlevels(genome_obj)
    gff <- genome_obj$gff
    # Recompute features using the same predicate for determinism
    features <- .prefer_name_then_fallback(gff, gene_pattern)
    # Reapply type preference and dedup
    if (isTRUE(translate)) {
      if ("type" %in% names(S4Vectors::mcols(features))) {
        cds_ok <- S4Vectors::mcols(features)$type %in% c("CDS")
        if (any(cds_ok)) features <- features[cds_ok]
      }
    } else {
      if ("type" %in% names(S4Vectors::mcols(features))) {
        gene_ok <- S4Vectors::mcols(features)$type %in% c("gene")
        if (any(gene_ok)) {
          features <- features[gene_ok]
        } else {
          cds_ok <- S4Vectors::mcols(features)$type %in% c("CDS")
          if (any(cds_ok)) features <- features[cds_ok]
        }
      }
    }
    features <- .drop_extdb_dupes(features)
    # Recheck
    fa_headers <- genome_obj$seqnames
    current_vocab_ok <- all(as.character(GenomicRanges::seqnames(features)) %in% fa_headers)
    if (!current_vocab_ok) {
      bad <- setdiff(unique(as.character(GenomicRanges::seqnames(features))), fa_headers)
      cli::cli_abort(c(
        "x" = "Seqnames mismatch persists after harmonization.",
        "i" = paste0("Offending seqnames: ", paste(bad, collapse = ", ")),
        "i" = "Inspect GFF 'region' rows and FASTA headers for order or naming issues."
      ))
    } else {
      cli::cli_alert_success("Seqname vocabulary reconciled; proceeding.")
    }
  } else if (!current_vocab_ok) {
    bad <- setdiff(unique(as.character(GenomicRanges::seqnames(features))), fa_headers)
    cli::cli_abort(c(
      "x" = "Seqnames mismatch between GFF and FASTA.",
      "i" = paste0("Offending seqnames: ", paste(bad, collapse = ", ")),
      "i" = "Set auto_harmonize=TRUE or run harmonize_gff_seqlevels() beforehand."
    ))
  }
  
  # Extract
  dna_seqs <- Biostrings::getSeq(genome_obj$fa, features)
  names(dna_seqs) <- .best_labels(features)
  
  if (isTRUE(translate)) {
    cli::cli_alert_info("Translating to amino acids using genetic code {.val {genetic_code}}")
    aa_seqs <- Biostrings::translate(dna_seqs, genetic.code = Biostrings::getGeneticCode(genetic_code))
    return(aa_seqs)
  }
  dna_seqs
}

#' Extract Sequence from Genomic Coordinates
#'
#' @param genome_obj Genome object from init_genome() (expects $fa, $fasta (optional), $seqnames)
#' @param seqname Chromosome/contig name or numeric-like index (e.g., "1")
#' @param start Start coordinate (1-based, inclusive)
#' @param end End coordinate (1-based, inclusive)
#' @param strand Strand ("+", "-" or "*"); default "+"
#' @param auto_resolve If TRUE, resolve seqname to a valid FASTA header (default TRUE)
#' @param clamp If TRUE, clamp coordinates into contig bounds with a warning (default FALSE)
#' @return DNAStringSet of length 1 with the extracted sequence
#' @export
extract_sequence_by_coords <- function(genome_obj,
                                       seqname,
                                       start,
                                       end,
                                       strand = "+",
                                       auto_resolve = TRUE,
                                       clamp = FALSE) {
  # --- required namespaces ---
  stopifnot(requireNamespace("cli", quietly = TRUE))
  stopifnot(requireNamespace("GenomicRanges", quietly = TRUE))
  stopifnot(requireNamespace("IRanges", quietly = TRUE))
  stopifnot(requireNamespace("Biostrings", quietly = TRUE))
  stopifnot(requireNamespace("Rsamtools", quietly = TRUE))
  stopifnot(requireNamespace("S4Vectors", quietly = TRUE))
  
  # --- basic checks ---
  stopifnot(is.list(genome_obj), "fa" %in% names(genome_obj), "seqnames" %in% names(genome_obj))
  if (!strand %in% c("+", "-", "*")) {
    cli::cli_abort("Strand must be one of '+', '-', or '*'.")
  }
  s <- as.integer(start); e <- as.integer(end)
  if (is.na(s) || is.na(e) || s < 1L || e < 1L || s > e) {
    cli::cli_abort("Start/end must be positive integers with start <= end.")
  }
  
  fa_headers <- genome_obj$seqnames
  
  # --- resolver: map seqname -> a valid FASTA header ---
  .resolve_seqname <- function(snm) {
    q <- as.character(snm)
    
    # 1) exact
    if (q %in% fa_headers) return(q)
    
    # 2) case-insensitive exact
    m <- match(tolower(q), tolower(fa_headers))
    if (!is.na(m)) return(fa_headers[m])
    
    # 3) numeric-like -> positional header (e.g., "1" -> first FASTA header)
    if (grepl("^[0-9]+$", q)) {
      i <- as.integer(q)
      if (!is.na(i) && i >= 1 && i <= length(fa_headers)) {
        cli::cli_alert_info("Interpreting seqname {.val {q}} as positional index -> {.val {fa_headers[i]}}")
        return(fa_headers[i])
      }
    }
    
    # 4) prefix match (case-insensitive)
    starts <- which(startsWith(tolower(fa_headers), tolower(q)))
    if (length(starts) == 1L) return(fa_headers[starts])
    
    # 5) fuzzy match (small distance)
    close <- try(agrep(q, fa_headers, ignore.case = TRUE, max.distance = 0.1), silent = TRUE)
    if (!inherits(close, "try-error") && length(close) == 1L) return(fa_headers[close])
    
    # Fail with a crisp diagnostic
    cli::cli_abort(c(
      "x" = "Unknown seqname.",
      "i" = paste0("Requested: ", q),
      "i" = paste0("Available FASTA headers: ", paste(fa_headers, collapse = ", "))
    ))
  }
  
  target <- if (isTRUE(auto_resolve)) .resolve_seqname(seqname) else as.character(seqname)
  
  # --- determine contig length for bound checks (prefer in-memory DNAStringSet if present) ---
  lengths <- NULL
  if (!is.null(genome_obj$fasta) && methods::is(genome_obj$fasta, "Biostrings::DNAStringSet")) {
    lengths <- Biostrings::width(genome_obj$fasta)
    names(lengths) <- names(genome_obj$fasta)
  } else if (!is.null(genome_obj$fasta) && methods::is(genome_obj$fasta, "DNAStringSet")) {
    lengths <- Biostrings::width(genome_obj$fasta)
    names(lengths) <- names(genome_obj$fasta)
  } else {
    # fallback to FA index
    fai <- try(Rsamtools::scanFaIndex(genome_obj$fa), silent = TRUE)
    if (!inherits(fai, "try-error")) {
      # scanFaIndex returns an index with fields; names(fai) are contig names; use $seqlengths
      lengths <- fai$seqlengths
      names(lengths) <- names(fai)
    }
  }
  
  # --- check and optionally clamp coordinates ---
  if (!is.null(lengths) && target %in% names(lengths)) {
    L <- as.integer(lengths[[target]])
    new_s <- s; new_e <- e
    if (!is.na(L)) {
      if (s < 1L || e > L) {
        if (!clamp) {
          cli::cli_abort(c(
            "x" = "Coordinates out of bounds.",
            "i" = sprintf("Requested %s:%d-%d (strand %s); contig length = %d", target, s, e, strand, L)
          ))
        } else {
          new_s <- max(1L, s); new_e <- min(L, e)
          cli::cli_warn(sprintf("Clamped coordinates to %s:%d-%d (contig length %d).", target, new_s, new_e, L))
        }
      }
    }
    s <- new_s; e <- new_e
  }
  
  # --- build ROI and extract ---
  roi <- GenomicRanges::GRanges(
    seqnames = target,
    ranges = IRanges::IRanges(start = s, end = e),
    strand = strand
  )
  
  # getSeq will automatically reverse-complement if strand == "-"
  seq <- Biostrings::getSeq(genome_obj$fa, roi)
  names(seq) <- sprintf("%s:%d-%d(%s)", target, s, e, strand)
  seq
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
#' @param genome_obj Genome object from init_genome() (expects $gff and $seqnames)
#' @param seqname Chromosome/contig name (can be numeric-like, e.g., "1")
#' @param start Start coordinate (1-based)
#' @param end End coordinate (1-based)
#' @param flank Flanking region size (default = 0)
#' @param feature_type Feature type to extract (default = "CDS")
#' @param tidy Return tibble (TRUE) or GRanges (FALSE)
#' @param auto_resolve If TRUE, resolve seqname against GFF seqlevels (default TRUE)
#' @param drop_extdb If TRUE, prefer human labels over extdb placeholders (default TRUE)
#' @return GRanges or tibble with features in ROI
#' @export
analyze_roi <- function(genome_obj, 
                        seqname, 
                        start, 
                        end, 
                        flank = 0,
                        feature_type = "CDS",
                        tidy = TRUE,
                        auto_resolve = TRUE,
                        drop_extdb = TRUE) {
  # --- namespaces ---
  stopifnot(requireNamespace("cli", quietly = TRUE))
  stopifnot(requireNamespace("GenomicRanges", quietly = TRUE))
  stopifnot(requireNamespace("IRanges", quietly = TRUE))
  stopifnot(requireNamespace("S4Vectors", quietly = TRUE))
  stopifnot(requireNamespace("GenomeInfoDb", quietly = TRUE))
  stopifnot(requireNamespace("BiocGenerics", quietly = TRUE))
  stopifnot(requireNamespace("tibble", quietly = TRUE))
  
  gff <- genome_obj$gff
  
  # --------- resolver: map user seqname -> a level present in the GFF ---------
  .resolve_seqname_gff <- function(q, levels) {
    q <- as.character(q)
    # exact
    if (q %in% levels) return(q)
    # case-insensitive exact
    m <- match(tolower(q), tolower(levels))
    if (!is.na(m)) return(levels[m])
    # numeric-like -> positional level
    if (grepl("^[0-9]+$", q)) {
      i <- as.integer(q)
      if (!is.na(i) && i >= 1L && i <= length(levels)) {
        cli::cli_alert_info("Interpreting seqname {.val {q}} as positional index in GFF -> {.val {levels[i]}}")
        return(levels[i])
      }
    }
    # prefix
    starts <- which(startsWith(tolower(levels), tolower(q)))
    if (length(starts) == 1L) return(levels[starts])
    # fuzzy (light)
    close <- try(agrep(q, levels, ignore.case = TRUE, max.distance = 0.1), silent = TRUE)
    if (!inherits(close, "try-error") && length(close) == 1L) return(levels[close])
    cli::cli_abort(c(
      "x" = "Unknown seqname for GFF.",
      "i" = paste0("Requested: ", q),
      "i" = paste0("GFF seqlevels: ", paste(GenomeInfoDb::seqlevels(gff), collapse = ", "))
    ))
  }
  
  # --------- tidy normalization for PGAP fields ---------
  .tidy_features_pgappref <- function(gr, drop_extdb = TRUE) {
    n <- length(gr)
    # pull common columns safely
    mc <- S4Vectors::mcols(gr)
    getc <- function(x) if (x %in% names(mc)) as.character(mc[[x]]) else rep(NA_character_, n)
    
    gene       <- getc("gene")
    locus_tag  <- getc("locus_tag")
    name_raw   <- getc("Name")
    product    <- getc("product")
    type       <- getc("type")
    
    # Name priority: gene -> locus_tag -> Name (unless extdb) -> product -> ID
    name_out <- gene
    fill <- is.na(name_out) | name_out == ""
    name_out[fill] <- locus_tag[fill]
    fill <- is.na(name_out) | name_out == ""
    # use Name if not extdb:*, otherwise delay
    ok_name <- ifelse(!is.na(name_raw) & !grepl("^extdb:", name_raw), name_raw, NA_character_)
    name_out[fill] <- ok_name[fill]
    fill <- is.na(name_out) | name_out == ""
    name_out[fill] <- product[fill]
    fill <- is.na(name_out) | name_out == ""
    id <- getc("ID")
    name_out[fill] <- id[fill]
    
    # If still NA, synthesize
    name_out[is.na(name_out) | name_out == ""] <- paste0("feature_", which(is.na(name_out) | name_out == ""))
    
    # If drop_extdb, and we still have extdb in Name while locus_tag or gene exist, replace
    if (isTRUE(drop_extdb)) {
      is_ext <- grepl("^extdb:", name_out)
      replace_with <- ifelse(!is.na(gene) & gene != "", gene,
                             ifelse(!is.na(locus_tag) & locus_tag != "", locus_tag, name_out))
      name_out[is_ext] <- replace_with[is_ext]
    }
    
    tibble::tibble(
      seqnames = as.character(GenomicRanges::seqnames(gr)),
      start    = BiocGenerics::start(gr),
      end      = BiocGenerics::end(gr),
      width    = BiocGenerics::width(gr),
      strand   = as.character(GenomicRanges::strand(gr)),
      type     = as.character(type),
      Name     = name_out,
      Alias    = locus_tag,            # maps your desired Alias to locus_tag
      Note     = product               # maps your desired Note to product
    )
  }
  
  # ------------------- build ROI -------------------
  s <- as.integer(start); e <- as.integer(end); f <- as.integer(flank)
  if (is.na(s) || is.na(e) || s < 1L || e < 1L || s > e) {
    cli::cli_abort("Start/end must be positive integers with start <= end.")
  }
  
  # Resolve seqname against the GFF vocabulary
  target <- if (isTRUE(auto_resolve)) {
    .resolve_seqname_gff(seqname, GenomeInfoDb::seqlevels(gff))
  } else {
    as.character(seqname)
  }
  
  roi <- GenomicRanges::GRanges(
    seqnames = target,
    ranges   = IRanges::IRanges(start = max(1L, s - f), end = e + f)
  )
  
  cli::cli_inform(sprintf("Analyzing region: %s:%d-%d (\u00B1%d bp)",
                          as.character(seqname), s, e, f))
  
  # ------------------- overlap -------------------
  roi_feats <- IRanges::subsetByOverlaps(gff, roi)

  # Filter type if requested; drop region scaffolding
  if (!is.null(feature_type)) {
    mcols <- S4Vectors::mcols(roi_feats)
    if ("type" %in% names(mcols)) {
      roi_feats <- roi_feats[mcols$type == feature_type]
    }
  }
  
  cli::cli_inform(sprintf("Found %d %s feature(s) in ROI",
                          length(roi_feats), 
                          ifelse(is.null(feature_type), 
                                 "matching", 
                                 feature_type)))
  
  if (!isTRUE(tidy)) {
    return(roi_feats)
  }
  
  .tidy_features_pgappref(roi_feats, drop_extdb = drop_extdb)
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

# THIS SECTION MAY NEED TO BE DEPRECATED ====
#' Synchronize seqnames between a GRanges and FASTA headers/obj.
#' - DNAStringSet: FASTA gets renamed to scrubbed names; GRanges stays scrubbed.
#' - FaFile: GRanges gets renamed back to RAW file headers; file stays untouched.
#'
#' @param features   GRanges to align
#' @param fa_headers character vector of FASTA headers (raw as in file or DNAStringSet)
#' @param fa_obj     DNAStringSet or FaFile (optional but recommended)
#' @param set_lengths logical; if TRUE, set seqlengths(features) from FASTA and check bounds
#' @param verbose    logical
#' @return list(features=GRanges, fa_names_used=character, fa_obj=DNAStringSet/FaFile/or NULL)
synchronize_seqnames <- function(features,
                                 fa_headers,
                                 fa_obj = NULL,
                                 set_lengths = TRUE,
                                 verbose = TRUE) {
  stopifnot(inherits(features, "GRanges"))
  if (length(fa_headers) == 0) stop("fa_headers is empty.")
  
  # Set conservative pruning behavior during level changes
  o_old <- getOption("GenomeInfoDb::NA.pruning.mode")
  on.exit(options("GenomeInfoDb::NA.pruning.mode" = o_old), add = TRUE)
  options("GenomeInfoDb::NA.pruning.mode" = "coarse")
  
  # FASTA raw vs scrubbed
  fa_raw  <- as.character(fa_headers)
  fa_sync <- .scrub_prefixes(fa_raw)
  
  # Scrub features' seqlevels first
  GenomeInfoDb::seqlevels(features) <- .scrub_prefixes(GenomeInfoDb::seqlevels(features))
  
  # Keep only common (scrubbed) contigs
  keep_sync <- intersect(GenomeInfoDb::seqlevels(features), fa_sync)
  if (length(keep_sync) == 0) stop("After synchronization, no overlapping contigs between GFF and FASTA.")
  GenomeInfoDb::seqlevels(features, pruning.mode = "coarse") <- keep_sync
  if (isTRUE(verbose)) message("Synchronized seqlevels and pruned to intersection (scrubbed): ",
                               paste(keep_sync, collapse = ", "))
  
  # Inspect FASTA lengths when available
  fa_idx <- .get_fasta_index(fa_obj)
  
  if (inherits(fa_obj, "DNAStringSet")) {
    # Rename DNAStringSet to scrubbed names; keep GRanges scrubbed
    if (!is.null(names(fa_obj))) {
      names(fa_obj) <- .scrub_prefixes(names(fa_obj))
    } else if (length(fa_sync) == length(fa_obj)) {
      names(fa_obj) <- fa_sync
    }
    fa_names_used <- fa_sync
    if (isTRUE(verbose)) message("Synchronized DNAStringSet names to scrubbed scheme.")
    
  } else if (inherits(fa_obj, "FaFile")) {
    # Map GRanges scrubbed names -> RAW file headers, so getSeq(FaFile, ...) matches the file index
    map_scrub_to_raw <- setNames(fa_raw, fa_sync)
    present <- intersect(GenomeInfoDb::seqlevels(features), names(map_scrub_to_raw))
    if (length(present) == 0) stop("No mapping from scrubbed names to raw FaFile headers found.")
    
    # renameSeqlevels old->new (scrubbed -> raw)
    features <- GenomeInfoDb::renameSeqlevels(features,
                                              map = map_scrub_to_raw[present],
                                              pruning.mode = "coarse")
    fa_names_used <- fa_raw
    if (isTRUE(verbose)) message("Renamed GRanges seqlevels to RAW FaFile headers (e.g., 'lcl|1').")
    
  } else {
    # No FA object—default to scrubbed scheme on both sides
    fa_names_used <- fa_sync
    if (isTRUE(verbose)) message("No FASTA object provided; using scrubbed names for checks.")
  }
  
  # Optionally set seqlengths from FASTA and check for out-of-bounds
  if (isTRUE(set_lengths) && length(fa_idx$lengths) > 0) {
    # Work in the same namespace features now uses
    current_levels <- GenomeInfoDb::seqlevels(features)
    # Choose lengths in the same naming used by features
    lengths_in_use <- fa_idx$lengths[names(fa_idx$lengths) %in% current_levels]
    si <- GenomeInfoDb::Seqinfo(seqnames = names(lengths_in_use),
                                seqlengths = as.integer(lengths_in_use))
    GenomeInfoDb::seqinfo(features) <- GenomeInfoDb::merge(seqinfo(features), si)
    
    # Bounds check
    sn <- as.character(GenomicRanges::seqnames(features))
    en <- BiocGenerics::end(features)
    sl <- GenomeInfoDb::seqlengths(features)[sn]
    oob <- which(en > sl)
    if (length(oob) > 0) {
      offenders <- unique(paste0(sn[oob], ":", BiocGenerics::start(features)[oob], "-", en[oob],
                                 " > ", sl[oob]))
      stop("Out-of-bounds feature coordinates for FASTA:\n",
           paste(head(offenders, 10), collapse = "\n"),
           if (length(offenders) > 10) sprintf("\n... and %d more", length(offenders) - 10) else "")
    }
  }
  
  list(features = features, fa_names_used = fa_names_used, fa_obj = fa_obj)
}

# Attempt sync only if needed
try_synchronize_if_needed <- function(features,
                                      fa_headers,
                                      fa_obj = NULL,
                                      verbose = TRUE) {
  sn_vals <- as.character(GenomicRanges::seqnames(features))
  if (length(sn_vals) == 0) stop("Features carry no seqnames.")
  
  # Choose which set of names to test against initially:
  # If fa_obj is a FaFile we must compare against RAW headers; else against provided fa_headers.
  initial_headers <- if (inherits(fa_obj, "FaFile")) fa_headers else fa_headers
  
  if (all(sn_vals %in% initial_headers)) {
    return(list(features = features, fa_names_used = initial_headers, fa_obj = fa_obj))
  }
  
  if (isTRUE(verbose)) {
    bad <- unique(sn_vals[!(sn_vals %in% initial_headers)])
    message("Sanity check failed. Offending seqnames: ", paste(bad, collapse = ", "),
            "\nAttempting synchronization...")
  }
  
  out <- synchronize_seqnames(features, fa_headers, fa_obj = fa_obj, set_lengths = TRUE, verbose = verbose)
  
  # Re-check with the names actually in use post-sync
  sn_re <- as.character(GenomicRanges::seqnames(out$features))
  if (!all(sn_re %in% out$fa_names_used)) {
    bad2 <- unique(sn_re[!(sn_re %in% out$fa_names_used)])
    stop("Seqnames mismatch persists after synchronization.\n",
         "Offenders: ", paste(bad2, collapse = ", "), "\n",
         "GFF values: ", paste(sort(unique(sn_re)), collapse = ", "), "\n",
         "FASTA: ", paste(out$fa_names_used, collapse = ", "))
  }
  
  out
}


extract_features <- function(genome_obj,
                             gene_pattern,
                             translate = FALSE,
                             genetic_code = "11",
                             verbose = TRUE) {
  # 1) Find matching features
  gff <- genome_obj$gff
  if (!inherits(gff, "GRanges")) stop("genome_obj$gff must be a GRanges.")
  if (!("Name" %in% colnames(S4Vectors::mcols(gff)))) {
    stop("genome_obj$gff must have a 'Name' metadata column.")
  }
  
  features <- gff[grepl(gene_pattern, gff$Name)]
  if (length(features) == 0) {
    warning("No features found matching pattern: ", gene_pattern)
    return(NULL)
  }
  if (isTRUE(verbose)) message(sprintf("Found %d feature(s) matching '%s'",
                                       length(features), gene_pattern))
  
  # 2) First-pass sanity check, then conditional synchronization
  fa_headers <- genome_obj$seqnames
  sn_vals <- as.character(GenomicRanges::seqnames(features))
  
  if (!all(sn_vals %in% fa_headers)) {
    sync <- try_synchronize_if_needed(features, fa_headers, fa_obj = genome_obj$fa, verbose = verbose)
    features <- sync$features
    # Carry forward synchronized DNAStringSet names if applicable
    if (inherits(sync$fa_obj, "DNAStringSet")) {
      genome_obj$fa <- sync$fa_obj
    }
  }
  
  # 3) Extract sequences (strand-aware)
  dna_seqs <- Biostrings::getSeq(genome_obj$fa, features)
  names(dna_seqs) <- features$Name
  
  # 4) Optional translation
  if (isTRUE(translate)) {
    if (isTRUE(verbose)) message("Translating to amino acids...")
    aa_seqs <- Biostrings::translate(
      dna_seqs,
      genetic.code = Biostrings::getGeneticCode(genetic_code)
    )
    return(aa_seqs)
  }
  
  dna_seqs
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
#project_folder <- "/home/william-ackerman/Desktop/NIST0056_CNV_ANALYSIS_ORGANIZED"
#setwd(project_folder)
#setwd("./01_spine/data_clean/prokka_NIST_chr")

dir <- here::here()
setwd(dir)
NIST_pgap_genome <- init_genome("my_genome.gff3", "my_genome.fna")
extract_sequences_by_name(NIST_pgap_genome, "acrB", translate = TRUE)

extract_sequences_by_name(NIST_pgap_genome, "acrB", translate = TRUE)
extract_sequence_by_coords(NIST_pgap_genome, "1", 1834322, 1837471, "+")

analyze_roi(NIST_pgap_genome , "1", 2369501, 2613000, flank = 10000)

dir <- here::here("prokka_NIST_chr")
setwd(dir)
NIST_genome <- init_genome("reference.gff3", "reference.fasta")
extract_sequences_by_name(NIST_genome, "acrB", translate = TRUE)
extract_sequences_by_name(NIST_genome, "^bla", translate = TRUE)
extract_sequences_by_name(NIST_genome, "^fl", translate = TRUE)
  
NIST_genome$gff %>% as.data.frame() %>% filter(type == "CDS") %>% 
    filter(str_detect(Name, "^acr"))
analyze_roi(NIST_genome, "1", 2369501, 2613000, flank = 10000)

