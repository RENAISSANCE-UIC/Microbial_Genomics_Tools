# gbk_parser.R
# A dependency-light GenBank (.gb / .gbk) parser for R
# Handles: multi-record files, wrapped header fields, FEATURES block, multi-line qualifiers,
#          location strings with complement(), join(), order(), <, >, ^, and remote segments.

trim <- function(x) sub("\\s+$", "", sub("^\\s+", "", x))

# Split a multi-record GBK file by lines ending with "//"
split_records_by_slashes <- function(lines) {
  # Normalize and sanitize
  lines <- sub("\r$", "", lines)      # strip CR from CRLF files
  lines[is.na(lines)] <- ""           # replace NA with empty strings
  
  slash_idx <- which(grepl("^//\\s*$", lines))
  chunks <- list()
  
  if (length(slash_idx)) {
    prev <- 0L
    for (e in slash_idx) {
      s <- prev + 1L
      if (s <= e) {
        chunk <- lines[s:e]
        # drop the trailing '//' line if present
        if (length(chunk) && grepl("^//\\s*$", chunk[length(chunk)])) {
          chunk <- chunk[-length(chunk)]
        }
        # keep only if it has a LOCUS and any non-empty line
        if (length(chunk) &&
            any(grepl("^LOCUS\\b", chunk)) &&
            any(nzchar(trimws(chunk)))) {
          chunks[[length(chunks) + 1L]] <- chunk
        }
      }
      prev <- e
    }
    # If there is trailing content after the last //, include it only if it looks like a record
    if (prev < length(lines)) {
      tail_chunk <- lines[(prev + 1L):length(lines)]
      if (any(grepl("^LOCUS\\b", tail_chunk)) &&
          any(nzchar(trimws(tail_chunk)))) {
        chunks[[length(chunks) + 1L]] <- tail_chunk
      }
    }
  } else {
    # No slashes at all: treat whole file as one record if it has LOCUS
    if (any(grepl("^LOCUS\\b", lines))) {
      chunks <- list(lines)
    }
  }
  
  chunks
}


# Collects a block that begins with a tag (e.g., "DEFINITION") and may span multiple lines
collect_tag_block <- function(lines, start_idx, tag) {
  # GenBank conventions: tag is in col 1..12, content starts ~col 13 but actual files vary
  lab <- paste0("^", tag, "\\b")
  i <- start_idx
  if (!grepl(lab, lines[i])) stop("Expected tag '", tag, "' at line ", i)
  # First line content after the tag label
  first <- sub(lab, "", lines[i])
  first <- sub("^\\s*", "", first)
  i <- i + 1
  collected <- first
  while (i <= length(lines)) {
    # Continuation lines usually begin with spaces in the tag column region
    if (i <= length(lines) && grepl("^\\s{2,}\\S", lines[i]) && !grepl("^\\S", lines[i])) {
      collected <- c(collected, trim(lines[i]))
      i <- i + 1
    } else if (i <= length(lines) && grepl("^\\s*$", lines[i])) {
      # allow empty lines within a block
      collected <- c(collected, "")
      i <- i + 1
    } else {
      break
    }
  }
  list(value = trim(paste(collected, collapse = " ")), next_idx = i)
}

# Parse LOCUS line tokens (best-effort; field order varies across eras/divisions)
parse_locus_line <- function(line) {
  # Example:
  # LOCUS       SCU49845     5028 bp    DNA     circular PLN 21-JUN-1999
  x <- trim(sub("^LOCUS\\s+", "", line))
  toks <- strsplit(x, "\\s+")[[1]]
  res <- list(locus = NA, length_bp = NA, mol_type = NA, topology = NA,
              division = NA, date = NA)
  if (length(toks) >= 1) res$locus <- toks[1]
  # find a number followed by 'bp'
  bp_idx <- which(grepl("^[0-9]+$", toks))
  if (length(bp_idx)) {
    # length is first pure number; next token may be 'bp'
    len_idx <- bp_idx[1]
    res$length_bp <- suppressWarnings(as.integer(toks[len_idx]))
    # detect 'bp' token or attached 'bp'
    if (len_idx + 1 <= length(toks) && grepl("^bp$", toks[len_idx + 1], ignore.case = TRUE)) {
      # ok
    } else if (grepl("bp$", toks[len_idx], ignore.case = TRUE)) {
      res$length_bp <- suppressWarnings(as.integer(sub("bp$", "", toks[len_idx], ignore.case = TRUE)))
    }
  }
  # Heuristics for mol_type/topology/division/date
  # known mol types
  mol_types <- c("DNA","RNA","ss-DNA","ds-DNA","ss-RNA","ds-RNA","mRNA","tRNA","rRNA","uRNA","snRNA","snoRNA")
  mt_idx <- which(toks %in% mol_types)
  if (length(mt_idx)) res$mol_type <- toks[mt_idx[1]]
  topo_idx <- which(tolower(toks) %in% c("linear","circular"))
  if (length(topo_idx)) res$topology <- tolower(toks[topo_idx[1]])
  # Date is often last token in dd-MMM-yyyy
  if (length(toks) >= 1 && grepl("^[0-9]{2}-[A-Z]{3}-[0-9]{4}$", toks[length(toks)])) {
    res$date <- toks[length(toks)]
    # division often right before date
    if (length(toks) >= 2 && grepl("^[A-Z]{3}$", toks[length(toks)-1])) {
      res$division <- toks[length(toks)-1]
    }
  }
  res
}

# Parse ORGANISM block: organism line, then taxonomy lineage across wrapped lines
parse_source_block <- function(lines, start_idx) {
  # SOURCE line:
  if (!grepl("^SOURCE\\b", lines[start_idx])) stop("Expected SOURCE at ", start_idx)
  source_val <- trim(sub("^SOURCE\\s*", "", lines[start_idx]))
  i <- start_idx + 1
  org <- NA_character_
  lineage <- NA_character_
  if (i <= length(lines) && grepl("^\\s+ORGANISM\\b", lines[i])) {
    org <- trim(sub("^\\s+ORGANISM\\s*", "", lines[i]))
    i <- i + 1
    tax <- character()
    while (i <= length(lines) && grepl("^\\s{2,}\\S", lines[i]) && !grepl("^\\S", lines[i])) {
      tax <- c(tax, trim(lines[i]))
      i <- i + 1
    }
    lineage <- trim(paste(tax, collapse = " "))
  }
  list(source = source_val, organism = org, taxonomy = lineage, next_idx = i)
}

# -------- Location string parsing -------------------------------------------
# Split by commas at top level (paren-balanced)
split_top_level_commas <- function(s) {
  out <- character()
  buf <- ""
  depth <- 0L
  for (i in seq_len(nchar(s))) {
    ch <- substr(s, i, i)
    if (ch == "(") depth <- depth + 1L
    if (ch == ")") depth <- max(0L, depth - 1L)
    if (ch == "," && depth == 0L) {
      out <- c(out, buf)
      buf <- ""
    } else if (ch != ",") {
      buf <- paste0(buf, ch)
    }
  }
  c(out, buf)
}

parse_single_range_token <- function(tok) {
  # Handles forms: 123..456, <123..456, 123..>456, <123..>456, 123, 123^124, remote like ACC:123..456
  res <- list(
    kind = "range",        # "range" | "site" | "point"
    start = NA_integer_, end = NA_integer_,
    start_fuzzy = FALSE, end_fuzzy = FALSE,
    remote_acc = NA_character_
  )
  t <- gsub("\\s", "", tok)
  # remote segment ACC:coords
  if (grepl("^[A-Za-z0-9_.]+:", t)) {
    acc <- sub(":.*$", "", t)
    coords <- sub("^[A-Za-z0-9_.]+:", "", t)
    res$remote_acc <- acc
    t <- coords
  }
  if (grepl("^[<>]?[0-9]+\\^<?[0-9]+$", t)) {
    # site between bases, e.g., 123^124
    parts <- strsplit(t, "\\^")[[1]]
    s <- sub("^<", "", parts[1])
    e <- sub("^<", "", parts[2])
    res$kind <- "site"
    res$start <- as.integer(s)
    res$end   <- as.integer(e)
    res$start_fuzzy <- grepl("^<", parts[1])
    res$end_fuzzy   <- grepl("^<", parts[2])
    return(res)
  }
  if (grepl("^[<>]?[0-9]+\\.\\.<?[0-9]+$", t)) {
    # range
    parts <- strsplit(t, "\\.\\.")[[1]]
    s <- parts[1]; e <- parts[2]
    res$start_fuzzy <- grepl("^<", s)
    res$end_fuzzy   <- grepl("^<|^>", e)
    s <- sub("^<", "", s)
    e <- sub("^<|^>", "", e)
    res$start <- as.integer(s)
    res$end   <- as.integer(e)
    return(res)
  }
  if (grepl("^[<>]?[0-9]+$", t)) {
    # single base point
    v <- as.integer(sub("^<|^>", "", t))
    res$kind <- "point"
    res$start <- v
    res$end   <- v
    res$start_fuzzy <- grepl("^<|^>", t)
    res$end_fuzzy   <- res$start_fuzzy
    return(res)
  }
  # Unrecognized—return as NA; caller can fall back to location_string
  res$kind <- "unknown"
  res
}

parse_location_string <- function(loc) {
  # Returns list(strand, location_type, ranges=data.frame, ok=TRUE/FALSE)
  s <- gsub("\\s+", "", loc)
  strand <- 1L
  location_type <- "single"
  # complement
  if (grepl("^complement\\(", s)) {
    strand <- -1L
    s <- sub("^complement\\(", "", s)
    s <- sub("\\)$", "", s)
  }
  # join/order
  if (grepl("^join\\(", s)) {
    location_type <- "join"
    inner <- sub("^join\\(|\\)$", "", s)
    toks <- split_top_level_commas(inner)
  } else if (grepl("^order\\(", s)) {
    location_type <- "order"
    inner <- sub("^order\\(|\\)$", "", s)
    toks <- split_top_level_commas(inner)
  } else {
    toks <- c(s)
  }
  parsed <- lapply(toks, parse_single_range_token)
  df <- data.frame(
    kind        = vapply(parsed, `[[`, "", "kind"),
    start       = suppressWarnings(as.integer(vapply(parsed, `[[`, 0L, "start"))),
    end         = suppressWarnings(as.integer(vapply(parsed, `[[`, 0L, "end"))),
    start_fuzzy = as.logical(vapply(parsed, `[[`, FALSE, "start_fuzzy")),
    end_fuzzy   = as.logical(vapply(parsed, `[[`, FALSE, "end_fuzzy")),
    remote_acc  = vapply(parsed, `[[`, NA_character_, "remote_acc"),
    stringsAsFactors = FALSE
  )
  ok <- !all(is.na(df$start) & is.na(df$end))
  list(strand = strand, location_type = location_type, ranges = df, ok = ok)
}

# -------- FEATURES parsing ---------------------------------------------------

parse_features_block <- function(lines, start_idx, end_idx) {
  # FEATURES line at start_idx; ORIGIN or end at end_idx
  # Feature keys typically start at col 6; qualifiers at col 22 with '/'
  feats <- list()
  i <- start_idx + 1
  current <- NULL
  in_qual <- FALSE
  
  new_feature_line <- function(line) {
    grepl("^\\s{5}\\S", line) && !grepl("^\\s{21}/", line)
  }
  qualifier_line <- function(line) {
    grepl("^\\s{21}/", line)
  }
  continuation_line <- function(line) {
    grepl("^\\s{21}\\S", line) && !qualifier_line(line)
  }
  
  get_key_loc <- function(line) {
    # key in ~cols 6-20; location starts ~col 21
    key <- trim(substr(line, 6, 20))
    loc <- trim(substr(line, 21, nchar(line)))
    list(key = key, loc = loc)
  }
  
  while (i <= end_idx) {
    line <- lines[i]
    if (new_feature_line(line)) {
      # store previous feature
      if (!is.null(current)) feats[[length(feats) + 1]] <- current
      kl <- get_key_loc(line)
      current <- list(
        type = kl$key,
        location_string = kl$loc,
        qualifiers = list(),
        raw_qual_lines = character()
      )
      in_qual <- FALSE
      i <- i + 1
      next
    }
    if (!is.null(current) && continuation_line(line) && !in_qual) {
      # location continuation
      current$location_string <- paste0(current$location_string, " ", trim(line))
      i <- i + 1
      next
    }
    if (!is.null(current) && qualifier_line(line)) {
      in_qual <- TRUE
      current$raw_qual_lines <- c(current$raw_qual_lines, line)
      i <- i + 1
      next
    }
    if (!is.null(current) && in_qual && grepl("^\\s{21}\\S", line)) {
      # qualifier wrapped continuation (no leading '/')
      current$raw_qual_lines <- c(current$raw_qual_lines, line)
      i <- i + 1
      next
    }
    # otherwise, end of FEATURES region
    i <- i + 1
  }
  # append last
  if (!is.null(current)) feats[[length(feats) + 1]] <- current
  
  # Stitch qualifiers: combine wrapped lines, split into /key="value"
  stitch_qualifiers <- function(raw_lines) {
    if (length(raw_lines) == 0) return(list())
    # strip left padding (~21 spaces)
    left <- sub("^\\s{21}", "", raw_lines)
    # join lines, but preserve that wrapped lines without leading '/' belong to previous qualifier
    q <- list()
    buf <- ""
    for (ln in left) {
      if (startsWith(ln, "/")) {
        if (nzchar(buf)) q <- c(q, buf)
        buf <- ln
      } else {
        # continuation
        buf <- paste(buf, ln, sep = " ")
      }
    }
    if (nzchar(buf)) q <- c(q, buf)
    
    # parse /key=value; values may be quoted strings with spaces and quotes
    out <- list()
    for (entry in q) {
      # /key or /key=value
      m <- sub("^/", "", entry)
      kv <- strsplit(m, "=", fixed = TRUE)[[1]]
      key <- kv[1]
      val <- if (length(kv) >= 2) paste(kv[-1], collapse = "=") else NA_character_
      if (!is.na(val)) {
        # remove surrounding quotes if present; keep internal quotes
        val <- trim(val)
        if (startsWith(val, "\"") && endsWith(val, "\"")) {
          val <- substr(val, 2, nchar(val) - 1)
        }
        # de-wrap residue: collapse multiple spaces
        val <- gsub("\\s+", " ", val)
      }
      # collect multiple values under same key as character vector
      if (is.null(out[[key]])) out[[key]] <- val else out[[key]] <- c(out[[key]], val)
    }
    out
  }
  
  # Build final data.frame
  if (length(feats) == 0) {
    return(data.frame(
      type = character(), location_string = character(),
      strand = integer(), location_type = character(),
      start = integer(), end = integer(),
      ranges = I(list()), qualifiers = I(list()),
      gene = character(), locus_tag = character(), product = character(),
      protein_id = character(), translation = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  parsed <- lapply(feats, function(f) {
    quals <- stitch_qualifiers(f$raw_qual_lines)
    locp  <- parse_location_string(f$location_string)
    # summarize coords
    rng <- locp$ranges
    st  <- if (all(is.na(rng$start))) NA_integer_ else min(rng$start, na.rm = TRUE)
    en  <- if (all(is.na(rng$end)))   NA_integer_ else max(rng$end,   na.rm = TRUE)
    list(
      type = f$type,
      location_string = f$location_string,
      strand = locp$strand,
      location_type = locp$location_type,
      start = st, end = en,
      ranges = rng,
      qualifiers = quals,
      gene = if (!is.null(quals$gene)) quals$gene[[1]] else NA_character_,
      locus_tag = if (!is.null(quals$locus_tag)) quals$locus_tag[[1]] else NA_character_,
      product = if (!is.null(quals$product)) quals$product[[1]] else NA_character_,
      protein_id = if (!is.null(quals$protein_id)) quals$protein_id[[1]] else NA_character_,
      translation = if (!is.null(quals$translation)) paste(quals$translation, collapse = " ") else NA_character_
    )
  })
  
  df <- data.frame(
    type = vapply(parsed, `[[`, "", "type"),
    location_string = vapply(parsed, `[[`, "", "location_string"),
    strand = as.integer(vapply(parsed, `[[`, 0L, "strand")),
    location_type = vapply(parsed, `[[`, "", "location_type"),
    start = as.integer(vapply(parsed, `[[`, 0L, "start")),
    end   = as.integer(vapply(parsed, `[[`, 0L, "end")),
    stringsAsFactors = FALSE
  )
  df$ranges     <- I(lapply(parsed, `[[`, "ranges"))
  df$qualifiers <- I(lapply(parsed, `[[`, "qualifiers"))
  df$gene        <- vapply(parsed, `[[`, NA_character_, "gene")
  df$locus_tag   <- vapply(parsed, `[[`, NA_character_, "locus_tag")
  df$product     <- vapply(parsed, `[[`, NA_character_, "product")
  df$protein_id  <- vapply(parsed, `[[`, NA_character_, "protein_id")
  df$translation <- vapply(parsed, `[[`, NA_character_, "translation")
  df
}

# Extract the ORIGIN sequence block
parse_origin_sequence <- function(lines, origin_idx, end_idx) {
  if (is.na(origin_idx)) return(NA_character_)
  seq_lines <- lines[(origin_idx + 1):end_idx]
  seq_text <- paste(seq_lines, collapse = "")
  seq_text <- tolower(seq_text)
  seq_text <- gsub("[^acgtun]", "", seq_text)
  toupper(seq_text)
}

# Parse a single GBK record
parse_gbk_record <- function(lines) {
  # Identify major sections
  idx_locus     <- which(grepl("^LOCUS\\b", lines))[1]
  if (is.na(idx_locus)) stop("No LOCUS line found in record.")
  # Header tags we'll parse: DEFINITION, ACCESSION, VERSION, KEYWORDS, SOURCE/ORGANISM
  # FEATURES block
  idx_features  <- which(grepl("^FEATURES\\b", lines))[1]
  idx_origin    <- which(grepl("^ORIGIN\\b", lines))[1]
  idx_end       <- max(which(grepl("^//\\s*$", lines)), length(lines))
  
  # Metadata
  locus <- parse_locus_line(lines[idx_locus])
  
  # DEFINITION (multi-line)
  idx_def <- which(grepl("^DEFINITION\\b", lines))
  definition <- NA_character_
  next_i <- idx_locus + 1
  if (length(idx_def)) {
    defblk <- collect_tag_block(lines, idx_def[1], "DEFINITION")
    definition <- defblk$value
    next_i <- defblk$next_idx
  }
  
  # ACCESSION
  accession <- NA_character_
  idx_acc <- which(grepl("^ACCESSION\\b", lines))
  if (length(idx_acc)) {
    # ACCESSION may wrap; take first token as primary
    acc_blk <- collect_tag_block(lines, idx_acc[1], "ACCESSION")
    acc_str <- acc_blk$value
    accession <- strsplit(acc_str, "\\s+")[[1]][1]
  }
  
  # VERSION
  version <- NA_character_; gi <- NA_character_
  idx_ver <- which(grepl("^VERSION\\b", lines))
  if (length(idx_ver)) {
    ver_line <- trim(sub("^VERSION\\s*", "", lines[idx_ver[1]]))
    # e.g., U49845.1  GI:1293613
    parts <- strsplit(ver_line, "\\s+")[[1]]
    if (length(parts)) version <- parts[1]
    gi_idx <- which(grepl("^GI:", parts))
    if (length(gi_idx)) gi <- sub("^GI:", "", parts[gi_idx[1]])
  }
  
  # KEYWORDS (optional)
  keywords <- NA_character_
  idx_kw <- which(grepl("^KEYWORDS\\b", lines))
  if (length(idx_kw)) {
    kw_blk <- collect_tag_block(lines, idx_kw[1], "KEYWORDS")
    keywords <- kw_blk$value
    if (identical(keywords, ".")) keywords <- NA_character_
  }
  
  # SOURCE/ORGANISM/taxonomy
  idx_src <- which(grepl("^SOURCE\\b", lines))
  source_info <- list(source = NA_character_, organism = NA_character_,
                      taxonomy = NA_character_, next_idx = NA_integer_)
  if (length(idx_src)) {
    source_info <- parse_source_block(lines, idx_src[1])
  }
  
  # FEATURES
  features <- data.frame()
  if (!is.na(idx_features)) {
    # Determine end of features block
    feat_end <- if (!is.na(idx_origin)) idx_origin - 1 else idx_end
    features <- parse_features_block(lines, idx_features, feat_end)
  }
  
  # ORIGIN sequence
  seq <- if (!is.na(idx_origin)) parse_origin_sequence(lines, idx_origin, idx_end) else NA_character_
  
  # Replicon guess from source feature qualifiers (/plasmid, /chromosome)
  replicon_type <- NA_character_
  replicon_name <- NA_character_
  if (nrow(features)) {
    src_rows <- which(features$type == "source")
    if (length(src_rows)) {
      q <- features$qualifiers[[src_rows[1]]]
      if (!is.null(q$plasmid)) {
        replicon_type <- "plasmid"
        replicon_name <- q$plasmid[[1]]
      } else if (!is.null(q$chromosome)) {
        replicon_type <- "chromosome"
        replicon_name <- q$chromosome[[1]]
      }
    }
  }
  # Topology fallback from LOCUS if present
  if (is.na(replicon_type) && !is.null(locus$topology)) {
    # leave as NA (topology != replicon), we only set if derived from qualifiers
  }
  
  list(
    metadata = list(
      locus = locus$locus,
      length_bp = locus$length_bp,
      mol_type = locus$mol_type,
      topology = locus$topology,
      division = locus$division,
      date = locus$date,
      definition = definition,
      accession = accession,
      version = version,
      gi = gi,
      keywords = keywords,
      source = source_info$source,
      organism = source_info$organism,
      taxonomy = source_info$taxonomy,
      replicon_type = replicon_type,
      replicon_name = replicon_name
    ),
    features = features,
    sequence = seq
  )
}

# Public API: read a .gb/.gbk path (or vector of paths) and parse
read_gbk <- function(path) {
  stopifnot(length(path) == 1)
  lines <- readLines(path, warn = FALSE)
  # Normalize tabs to spaces
  lines <- gsub("\t", "    ", lines, fixed = TRUE)
  recs <- split_records_by_slashes(lines)
  lapply(recs, parse_gbk_record)
}

# Convenience: bind features across multi-record files with record metadata
bind_gbk_features <- function(gbk_list) {
  out <- list()
  for (i in seq_along(gbk_list)) {
    r <- gbk_list[[i]]
    if (!nrow(r$features)) next
    meta <- r$metadata
    df <- r$features
    # append record-level metadata as columns
    df$record_index <- i
    df$accession    <- meta$accession
    df$locus        <- meta$locus
    df$organism     <- meta$organism
    df$replicon_type <- meta$replicon_type
    df$replicon_name <- meta$replicon_name
    out[[length(out) + 1]] <- df
  }
  if (length(out)) do.call(rbind, out) else
    data.frame()
}

# Utility: write a quick CSV of CDS features
write_cds_table <- function(gbk_list, file) {
  feats <- bind_gbk_features(gbk_list)
  if (!nrow(feats)) {
    warning("No features to write.")
    return(invisible(NULL))
  }
  cds <- subset(feats, type == "CDS", select = c(
    accession, locus, organism, replicon_type, replicon_name,
    gene, locus_tag, product, protein_id, start, end, strand, location_string
  ))
  write.csv(cds, file, row.names = FALSE)
  invisible(file)
}

# Helper function for reverse complement
reverse_complement <- function(seq) {
  seq <- toupper(seq)
  complement_map <- c("A"="T", "T"="A", "G"="C", "C"="G", "N"="N")
  rev_seq <- rev(strsplit(seq, "")[[1]])
  comp_seq <- complement_map[rev_seq]
  paste(comp_seq, collapse="")
}


# ==== Using parser ====
#The specific mutation here is 2,312,039	Δ8 bp	50.8%	coding (704‑711/1104 nt)	ompC ←	outer membrane porin C lms4 

setwd("C:\\Users\\willi\\Downloads")
gb <- read_gbk("Ecoli_K12_MG1655.gbk")

features <- bind_gbk_features(gb)

# Find ompC gene (could be CDS or gene feature)
ompC_features <- features[grepl("ompC", features$gene, ignore.case = TRUE) | 
                            grepl("ompC", features$product, ignore.case = TRUE) |
                            grepl("ompC", features$locus_tag, ignore.case = TRUE), ]
mutation_pos <- 2312039
deletion_size <- 8
coding_pos_affected <- "704-711" 

ompC_features <- features[features$start <= mutation_pos & 
           features$end >= mutation_pos & 
           features$type == "CDS", ]

ompC_cds <- ompC_features[1, ]
sequence <- gb[[1]]$sequence
ompC_seq <- substr(sequence, ompC_cds$start, ompC_cds$end)
ompC_seq <- reverse_complement(ompC_seq)

deletion_size %% 3
del_start_codon_pos <- 704
# Show sequence around deletion site
del_nt_pos <- (del_start_codon_pos - 1) * 3 + 1  # Convert codon position to nucleotide

analyze_theoretical_frameshift <- function(codon_pos_start, deletion_size, total_codons) {
  cat("\n=== THEORETICAL FRAMESHIFT ANALYSIS ===\n")
  
  if(deletion_size %% 3 == 0) {
    cat("Deletion size (", deletion_size, " bp) is a multiple of 3\n")
    cat("This is an IN-FRAME deletion\n")
    cat("Effect: Removes", deletion_size/3, "amino acids starting at codon", ceiling(codon_pos_start/3), "\n")
    cat("Remaining protein length:", total_codons/3 - deletion_size/3, "amino acids\n")
    cat("Protein likely remains functional but with altered properties\n")
  } else {
    cat("Deletion size (", deletion_size, " bp) is NOT a multiple of 3\n")
    cat("This WILL cause a frameshift mutation\n")
    
    frame_shift <- deletion_size %% 3
    affected_codon <- ceiling(codon_pos_start / 3)
    
    cat("Frame shift:", frame_shift, "nucleotides\n")
    cat("Starting at codon ~", affected_codon, "\n")
    cat("Total codons in gene:", total_codons/3, "\n")
    
    # Calculate impact
    remaining_codons <- total_codons/3 - affected_codon + 1
    cat("Codons affected by frameshift:", remaining_codons, "\n")
    
    cat("\nPREDICTED IMPACT:\n")
    cat("- Frameshift will alter all amino acids downstream of position", affected_codon, "\n")
    cat("- High probability of premature stop codon\n")
    cat("- Likely results in truncated, non-functional protein\n")
    cat("- May trigger nonsense-mediated decay of mRNA\n")
    cat("- Expected to severely disrupt or eliminate ompC function\n")
  }
}

# Function to show translation impact (simplified)
show_translation_impact <- function(before_seq, deleted_seq, after_seq, del_pos) {
  # This is a simplified version - in practice you'd want more sophisticated translation
  cat("Note: Full translation analysis requires proper start codon identification\n")
  cat("and reading frame determination from the complete CDS sequence.\n")
}


if(del_nt_pos > 0 && del_nt_pos <= nchar(sequence)) {
  cat("\nSequence context (50bp before and after deletion site):\n")
  start_show <- max(1, del_nt_pos - 50)
  end_show <- min(nchar(sequence), del_nt_pos + deletion_size + 50)
  
  before_del <- substr(sequence, start_show, del_nt_pos - 1)
  deleted_seq <- substr(sequence, del_nt_pos, del_nt_pos + deletion_size - 1)
  after_del <- substr(sequence, del_nt_pos + deletion_size, end_show)
  
  cat("Before deletion:", before_del, "\n")
  cat("Deleted sequence: [", deleted_seq, "]\n")
  cat("After deletion:", after_del, "\n")
  
  # Translate before and after
  cat("\nTranslation impact:\n")
  show_translation_impact(before_del, deleted_seq, after_del, del_nt_pos)
}


deletion_context <- list(
  before = "TCACGCCGAACAAAAAGGCCAACACCTCGTCGATGGATTACTACCATCAG",
  deleted = "TTGCGTTA", # 8 bp deletion
  after = "TGCGGCGGAAAAATCGCGGCGTAAATTCCTCTATGACACCAACGTTGGGGC"
)

before_seq <- deletion_context$before
deleted_seq <- deletion_context$deleted
after_seq <- deletion_context$after

# Construct the search pattern (last 20bp of before + deleted + first 20bp of after)
search_pattern <- paste0(substr(before_seq, nchar(before_seq)-19, nchar(before_seq)), 
                         deleted_seq,
                         substr(after_seq, 1, 20))

original_seq <- ompC_seq
deletion_pos <- regexpr(search_pattern, original_seq, fixed = TRUE)

if(deletion_pos == -1) {
  # Try alternative approach - find the before sequence
  before_pattern <- substr(before_seq, nchar(before_seq)-19, nchar(before_seq))
  deletion_pos <- regexpr(before_pattern, original_seq, fixed = TRUE)
  deletion_pos <- deletion_pos + 20  # adjust for the search pattern
}

cat("Original ompC sequence length:", nchar(original_seq), "bp\n")
cat("Deletion: 8 bp (", deleted_seq, ")\n")
cat("Deletion position in sequence: ~", deletion_pos, "\n\n")

# Create mutant sequence
if(deletion_pos > 0) {
  # Find exact position by pattern matching
  full_context <- paste0(before_seq, deleted_seq, after_seq)
  context_pos <- regexpr(full_context, original_seq, fixed = TRUE)
  
  if(context_pos > 0) {
    before_del_pos <- context_pos + nchar(before_seq) - 1
    after_del_pos <- before_del_pos + nchar(deleted_seq) + 1
    
    mutant_seq <- paste0(
      substr(original_seq, 1, before_del_pos),
      substr(original_seq, after_del_pos, nchar(original_seq))
    )
  } else {
    # Fallback: create approximate mutant sequence
    cat("Using approximate deletion position\n")
    mutant_seq <- paste0(
      substr(original_seq, 1, deletion_pos - 1),
      substr(original_seq, deletion_pos + 8, nchar(original_seq))
    )
  }
} else {
  cat("Could not locate exact deletion site, creating approximate mutant sequence\n")
  # For demonstration, we'll create the mutant by removing 8bp from around position 704
  del_start <- 704
  mutant_seq <- paste0(
    substr(original_seq, 1, del_start - 1),
    substr(original_seq, del_start + 8, nchar(original_seq))
  )
}

cat("Mutant sequence length:", nchar(mutant_seq), "bp\n")
cat("Length difference:", nchar(original_seq) - nchar(mutant_seq), "bp\n\n")

# Translate both sequences
translate_sequence <- function(seq) {
  genetic_code <- c(
    "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L",
    "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
    "TAT"="Y", "TAC"="Y", "TAA"="*", "TAG"="*",
    "TGT"="C", "TGC"="C", "TGA"="*", "TGG"="W",
    "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
    "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
    "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
    "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R",
    "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M",
    "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
    "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K",
    "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
    "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
    "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
    "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E",
    "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"
  )
  
  protein <- ""
  for(i in seq(1, nchar(seq) - 2, 3)) {
    codon <- substr(seq, i, i + 2)
    if(nchar(codon) == 3) {
      aa <- genetic_code[[codon]]
      if(is.null(aa)) aa <- "X"  # unknown
      protein <- paste0(protein, aa)
      if(aa == "*") break  # stop at first stop codon
    }
  }
  protein
}

original_protein <- translate_sequence(original_seq)
mutant_protein <- translate_sequence(mutant_seq)

cat("TRANSLATION RESULTS:\n")
cat("Original protein length:", nchar(gsub("\\*", "", original_protein)), "amino acids\n")
cat("Mutant protein length:", nchar(gsub("\\*", "", mutant_protein)), "amino acids\n\n")

# Find where sequences diverge
min_len <- min(nchar(original_protein), nchar(mutant_protein))
divergence_pos <- 1

for(i in 1:min_len) {
  if(substr(original_protein, i, i) != substr(mutant_protein, i, i)) {
    divergence_pos <- i
    break
  }
}

cat("Sequences diverge at amino acid position:", divergence_pos, "\n")

# Show sequence alignment around divergence
context_start <- max(1, divergence_pos - 10)
context_end <- min(min_len, divergence_pos + 20)

cat("\nSEQUENCE ALIGNMENT (around frameshift site):\n")
cat("Original: ", substr(original_protein, context_start, context_end), "\n")
cat("Mutant  : ", substr(mutant_protein, context_start, context_end), "\n")
cat("Position: ", paste(context_start:context_end, collapse=""), "\n\n")

# Check for premature stop codons
original_stops <- gregexpr("\\*", original_protein)[[1]]
mutant_stops <- gregexpr("\\*", mutant_protein)[[1]]

cat("STOP CODON ANALYSIS:\n")
if(original_stops[1] > 0) {
  cat("Original protein stops at position:", original_stops[1], "\n")
} else {
  cat("Original protein: no stop codon found\n")
}

if(mutant_stops[1] > 0) {
  cat("Mutant protein stops at position:", mutant_stops[1], "\n")
  if(original_stops[1] > 0 && mutant_stops[1] < original_stops[1]) {
    cat("*** PREMATURE STOP CODON DETECTED ***\n")
    cat("Protein truncated by", original_stops[1] - mutant_stops[1], "amino acids\n")
  }
} else {
  cat("Mutant protein: no stop codon found in translated region\n")
}

cat("\nCONCLUSION:\n")
cat("The 8bp deletion creates a frameshift mutation that:\n")
cat("1. Changes the reading frame starting at position ~", divergence_pos, "\n")
cat("2. Alters all downstream amino acids\n")
if(mutant_stops[1] > 0 && (original_stops[1] < 0 || mutant_stops[1] < original_stops[1])) {
  cat("3. Introduces a premature stop codon\n")
  cat("4. Results in a truncated, likely non-functional protein\n")
} else {
  cat("3. May affect protein folding and function due to altered sequence\n")
}
cat("5. Expected to severely disrupt OmpC porin function\n")