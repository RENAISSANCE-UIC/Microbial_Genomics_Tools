# NEW GENBANK PARSER
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




# -------- FASTA export ====

# -------- FASTA export (revised) ---------------------------------------------

#' Write sequences from parsed GBK records to FASTA
#'
#' @param gbk_list List of parsed GBK records (output from read_gbk())
#' @param file Output .fna path
#' @param wrap_width Integer; wrap sequences at this width (default 80)
#' @param auto_label Logical; if TRUE (default), infer chromosome vs plasmids
#'   by sorting by length (longest = chromosome, rest = plasmids)
#' @param include_length Logical; append [length=X] to header (default TRUE)
#'
#' @return Invisibly returns output path
#' @export
write_gbk_fasta <- function(gbk_list, file, wrap_width = 80L,
                            auto_label = TRUE, include_length = TRUE) {
  if (length(gbk_list) == 0L) {
    cli::cli_warn("No records to write.")
    return(invisible(NULL))
  }
  
  # Helper: wrap sequence at specified width
  wrap_seq <- function(seq_str, width) {
    if (is.na(seq_str) || !nzchar(seq_str)) return("")
    if (width <= 0L) return(seq_str)
    pos <- seq(1L, nchar(seq_str), by = width)
    paste(substring(seq_str, pos, pos + width - 1L), collapse = "\n")
  }
  
  # Compute actual sequence lengths
  seq_lengths <- vapply(gbk_list, function(r) {
    seq <- r$sequence
    if (is.na(seq) || !nzchar(seq)) return(0L)
    nchar(seq)
  }, integer(1L))
  
  # Sort by length (descending) if auto_label enabled
  if (auto_label) {
    ord <- order(seq_lengths, decreasing = TRUE)
    gbk_list <- gbk_list[ord]
    seq_lengths <- seq_lengths[ord]
  }
  
  # Filter out records without sequence
  has_seq <- seq_lengths > 0L
  if (!any(has_seq)) {
    cli::cli_warn("No valid sequences to write.")
    return(invisible(NULL))
  }
  gbk_list <- gbk_list[has_seq]
  seq_lengths <- seq_lengths[has_seq]
  
  # Build FASTA entries
  entries <- character(length(gbk_list))
  plasmid_counter <- 1L
  
  for (i in seq_along(gbk_list)) {
    r <- gbk_list[[i]]
    meta <- r$metadata
    seq <- r$sequence
    
    # Determine replicon label
    if (auto_label) {
      if (i == 1L) {
        # First (longest) = chromosome
        # Check if definition mentions "chromosome"
        is_chr <- !is.na(meta$definition) && 
          grepl("chromosome", meta$definition, ignore.case = TRUE)
        replicon_label <- if (is_chr) "Chromosome" else "Chromosome"
      } else {
        # Subsequent = plasmids
        replicon_label <- sprintf("Plasmid_%d_(presumptive)", plasmid_counter)
        plasmid_counter <- plasmid_counter + 1L
      }
    } else {
      # Use existing metadata if present
      if (!is.na(meta$replicon_type) && !is.na(meta$replicon_name)) {
        replicon_label <- paste0(meta$replicon_type, "_", meta$replicon_name)
      } else if (!is.na(meta$accession)) {
        replicon_label <- meta$accession
      } else if (!is.na(meta$locus)) {
        replicon_label <- meta$locus
      } else {
        replicon_label <- sprintf("contig_%d", i)
      }
    }
    
    # Build header parts
    hdr_parts <- replicon_label
    if (!is.na(meta$organism)) {
      hdr_parts <- paste(hdr_parts, meta$organism)
    }
    if (!is.na(meta$topology)) {
      hdr_parts <- paste(hdr_parts, sprintf("[topology=%s]", meta$topology))
    }
    if (include_length) {
      hdr_parts <- paste(hdr_parts, sprintf("[length=%d]", seq_lengths[i]))
    }
    
    entries[i] <- paste0(">", hdr_parts, "\n", wrap_seq(seq, wrap_width))
  }
  
  writeLines(entries, file)
  cli::cli_inform("Wrote {length(entries)} sequence{?s} to {.path {file}}")
  invisible(file)
}


# -------- GFF3 export ====

#' Convert parsed GBK features to GFF3 format
#'
#' @param gbk_list List of parsed GBK records
#' @param file Output .gff3 path
#' @param include_sequence_region Logical; include ##sequence-region directives (default TRUE)
#' @param include_types Character vector; feature types to include (default: all)
#'
#' @details
#' GFF3 spec: coordinates are 1-based inclusive. Implements NCBI/INSDC conventions
#' for feature hierarchies (gene -> mRNA/CDS), ID/Parent relationships, and
#' standard attributes (gbkey, gene_biotype, etc.).
#'
#' @return Invisibly returns output path
#' @export
write_gbk_gff3 <- function(gbk_list, file, include_sequence_region = TRUE,
                           include_types = NULL) {
  if (length(gbk_list) == 0L) {
    cli::cli_warn("No records to convert.")
    return(invisible(NULL))
  }
  
  # Helper: URL-encode special chars in GFF3 attributes
  gff_escape <- function(x) {
    if (is.na(x) || !nzchar(x)) return("")
    x <- gsub("%", "%25", x, fixed = TRUE)
    x <- gsub(";", "%3B", x, fixed = TRUE)
    x <- gsub("=", "%3D", x, fixed = TRUE)
    x <- gsub("&", "%26", x, fixed = TRUE)
    x <- gsub(",", "%2C", x, fixed = TRUE)
    x <- gsub("\t", "%09", x, fixed = TRUE)
    x <- gsub("\n", "%0A", x, fixed = TRUE)
    x
  }
  
  # Helper: format GFF3 attributes (order matters for readability)
  format_attributes <- function(attrs_list) {
    # Standard order: ID, Name, Parent, Dbxref, Ontology_term, Note, then alphabetical
    priority <- c("ID", "Name", "Parent", "Dbxref", "Ontology_term", "Note")
    keys <- names(attrs_list)
    # Remove empty values
    attrs_list <- attrs_list[!vapply(attrs_list, function(x) is.na(x) || !nzchar(x), logical(1L))]
    if (length(attrs_list) == 0L) return(".")
    
    keys <- names(attrs_list)
    ordered_keys <- c(
      intersect(priority, keys),
      setdiff(keys, priority)
    )
    
    pairs <- vapply(ordered_keys, function(k) {
      v <- attrs_list[[k]]
      paste0(k, "=", gff_escape(v))
    }, character(1L))
    
    paste(pairs, collapse = ";")
  }
  
  # Helper: infer source column from qualifiers
  infer_source <- function(qualifiers) {
    # Check inference qualifier for method
    inf <- qualifiers[["inference"]]
    if (!is.null(inf) && !is.na(inf[[1]])) {
      # e.g., "COORDINATES: similar to AA sequence:RefSeq:NP_416675.1"
      if (grepl("Protein Homology", inf[[1]], ignore.case = TRUE)) return("Protein Homology")
      if (grepl("ab initio", inf[[1]], ignore.case = TRUE)) return("Ab initio prediction")
      if (grepl("RefSeq", inf[[1]], ignore.case = TRUE)) return("RefSeq")
    }
    "."  # unknown
  }
  
  # Helper: compute CDS phase
  compute_phase <- function(ranges_df, strand) {
    if (nrow(ranges_df) == 0L) return(".")
    # Phase for first segment is 0; subsequent depend on cumulative length
    phases <- rep("0", nrow(ranges_df))
    if (nrow(ranges_df) > 1L) {
      cumlen <- 0L
      for (i in seq_len(nrow(ranges_df))) {
        seg_len <- ranges_df$end[i] - ranges_df$start[i] + 1L
        phases[i] <- as.character(cumlen %% 3L)
        cumlen <- cumlen + seg_len
      }
    }
    phases
  }
  
  # Collect all output lines
  gff_lines <- list()
  seq_region_lines <- list()
  
  # Process each record
  for (rec_idx in seq_along(gbk_list)) {
    r <- gbk_list[[rec_idx]]
    meta <- r$metadata
    feats <- r$features
    
    # Seqid: use locus (the contig/chromosome identifier in LOCUS line)
    seqid <- if (!is.na(meta$locus) && nzchar(meta$locus)) {
      meta$locus
    } else if (!is.na(meta$accession)) {
      meta$accession
    } else {
      paste0("unknown_", rec_idx)
    }
    
    # Sequence length (from actual sequence)
    seq_len <- if (!is.na(r$sequence) && nzchar(r$sequence)) {
      nchar(r$sequence)
    } else {
      meta$length_bp %||% NA_integer_
    }
    
    # Add ##sequence-region directive
    if (include_sequence_region && !is.na(seq_len)) {
      seq_region_lines[[length(seq_region_lines) + 1L]] <- 
        sprintf("##sequence-region %s 1 %d", seqid, seq_len)
    }
    
    # Add top-level "region" feature (from "source" feature in GenBank)
    if (nrow(feats) > 0L) {
      src_idx <- which(feats$type == "source")
      if (length(src_idx) > 0L) {
        src_feat <- feats[src_idx[1], ]
        src_quals <- src_feat$qualifiers[[1]]
        
        # Extract taxon from db_xref
        taxon <- NA_character_
        if (!is.null(src_quals$db_xref)) {
          dbx <- src_quals$db_xref
          taxon_match <- grep("^taxon:", dbx, value = TRUE)
          if (length(taxon_match)) {
            taxon <- sub("^taxon:", "", taxon_match[1])
          }
        }
        
        # Region attributes
        region_attrs <- list(
          ID = sprintf("%s:1..%d", seqid, seq_len),
          Dbxref = if (!is.na(taxon)) paste0("taxon:", taxon) else NA_character_,
          Name = meta$definition %||% "ANONYMOUS",
          gbkey = "Src",
          genome = if (rec_idx == 1L) "chromosome" else "plasmid",
          mol_type = meta$mol_type %||% "genomic DNA"
        )
        
        region_line <- paste(
          seqid, "Local", "region", 1, seq_len, ".", "+", ".",
          format_attributes(region_attrs),
          sep = "\t"
        )
        gff_lines[[length(gff_lines) + 1L]] <- region_line
      }
    }
    
    # Filter features by type if requested
    if (!is.null(include_types)) {
      feats <- feats[feats$type %in% include_types, , drop = FALSE]
    }
    
    if (!nrow(feats)) next
    
    # Track gene IDs for Parent relationships
    gene_map <- list()  # locus_tag -> gene ID
    
    # First pass: create genes
    for (i in seq_len(nrow(feats))) {
      feat <- feats[i, ]
      if (feat$type != "gene") next
      
      quals <- feat$qualifiers[[1]]
      locus_tag <- feat$locus_tag
      if (is.na(locus_tag)) next
      
      gene_id <- paste0("gene-", locus_tag)
      gene_map[[locus_tag]] <- gene_id
      
      strand_char <- switch(as.character(feat$strand), "1" = "+", "-1" = "-", ".")
      
      gene_attrs <- list(
        ID = gene_id,
        Name = feat$gene %||% locus_tag,
        gbkey = "Gene",
        gene_biotype = "protein_coding",
        locus_tag = locus_tag
      )
      
      if (!is.na(feat$gene)) gene_attrs$gene <- feat$gene
      
      # Check for partial annotations
      if (!is.null(quals$partial)) {
        gene_attrs$partial <- "true"
      }
      
      gene_line <- paste(
        seqid, ".", "gene", feat$start, feat$end, ".", strand_char, ".",
        format_attributes(gene_attrs),
        sep = "\t"
      )
      gff_lines[[length(gff_lines) + 1L]] <- gene_line
    }
    
    # Second pass: create CDS and other features
    for (i in seq_len(nrow(feats))) {
      feat <- feats[i, ]
      feat_type <- feat$type
      
      # Skip source (already handled) and gene (handled above)
      if (feat_type %in% c("source", "gene")) next
      
      quals <- feat$qualifiers[[1]]
      ranges_df <- feat$ranges[[1]]
      locus_tag <- feat$locus_tag
      
      strand_char <- switch(as.character(feat$strand), "1" = "+", "-1" = "-", ".")
      source_col <- infer_source(quals)
      
      # Determine ID and Parent
      feat_id <- if (feat_type == "CDS" && !is.na(feat$protein_id)) {
        paste0("cds-", feat$protein_id)
      } else if (!is.na(locus_tag)) {
        paste0(tolower(feat_type), "-", locus_tag)
      } else {
        paste0(seqid, "_", feat_type, "_", i)
      }
      
      parent_id <- if (!is.na(locus_tag) && !is.null(gene_map[[locus_tag]])) {
        gene_map[[locus_tag]]
      } else {
        NA_character_
      }
      
      # Build attributes
      attrs <- list(
        ID = feat_id,
        Name = feat$protein_id %||% feat$product %||% feat$gene %||% NA_character_,
        Parent = parent_id,
        gbkey = toupper(feat_type)
      )
      
      # Add common qualifiers
      if (!is.na(locus_tag)) attrs$locus_tag <- locus_tag
      if (!is.na(feat$gene)) attrs$gene <- feat$gene
      if (!is.na(feat$product)) attrs$product <- feat$product
      if (!is.na(feat$protein_id)) attrs$protein_id <- feat$protein_id
      
      # Add inference
      if (!is.null(quals$inference)) {
        attrs$inference <- quals$inference[[1]]
      }
      
      # Parse and add GO terms
      if (!is.null(quals$note)) {
        note <- quals$note[[1]]
        # Extract GO terms from note
        go_funcs <- gregexpr("GO:[0-9]{7}", note, perl = TRUE)
        if (go_funcs[[1]][1] != -1L) {
          go_terms <- regmatches(note, go_funcs)[[1]]
          if (length(go_terms)) {
            attrs$Ontology_term <- paste(go_terms, collapse = ",")
          }
        }
      }
      
      # Add partial flag
      if (!is.null(quals$partial)) {
        attrs$partial <- "true"
      }
      
      # Add transl_table for CDS
      if (feat_type == "CDS" && !is.null(quals$transl_table)) {
        attrs$transl_table <- quals$transl_table[[1]]
      }
      
      # Compute phase for CDS
      phases <- if (feat_type == "CDS") {
        compute_phase(ranges_df, feat$strand)
      } else {
        rep(".", nrow(ranges_df))
      }
      
      # Write one line per segment (for join/order)
      for (seg_idx in seq_len(nrow(ranges_df))) {
        seg <- ranges_df[seg_idx, ]
        if (is.na(seg$start) || is.na(seg$end)) next
        
        phase_char <- phases[seg_idx]
        
        # Attributes only on first segment
        attrs_str <- if (seg_idx == 1L) {
          format_attributes(attrs)
        } else {
          format_attributes(list(ID = paste0(feat_id, ".seg", seg_idx), Parent = feat_id))
        }
        
        line <- paste(
          seqid, source_col, feat_type, seg$start, seg$end, ".", strand_char, phase_char, attrs_str,
          sep = "\t"
        )
        gff_lines[[length(gff_lines) + 1L]] <- line
      }
    }
  }
  
  if (length(gff_lines) == 0L) {
    cli::cli_warn("No features to write to GFF3.")
    return(invisible(NULL))
  }
  
  # Write output
  header <- "##gff-version 3"
  all_lines <- c(header, unlist(seq_region_lines), unlist(gff_lines))
  writeLines(all_lines, file)
  cli::cli_inform("Wrote {length(gff_lines)} feature line{?s} to {.path {file}}")
  invisible(file)
}

# -------- LOGGING ====
logger <- function(level = "INFO") {
  function(msg, ...) {
    if (Sys.getenv("LOG_LEVEL", "INFO") == "DEBUG" || level != "DEBUG") {
      cli::cli_inform(paste0("[", level, "] ", msg), ...)
    }
  }
}
# log_info <- logger("INFO")
# log_warn <- logger("WARN")
# log_debug <- logger("DEBUG")
# -------- Unit Tests ====
ex <- '
LOCUS       TEST0001       1200 bp    DNA     circular  BCT  01-JAN-2020
DEFINITION  Example plasmid with two CDS features and a joined gene.
ACCESSION   TEST0001
VERSION     TEST0001.1  GI:999999999
KEYWORDS    .
SOURCE      Synthetic construct
  ORGANISM  Artificial sequence
            other sequences; synthetic sequences.
REFERENCE   1  (bases 1 to 1200)
FEATURES             Location/Qualifiers
     source          1..1200
                     /organism="Artificial sequence"
                     /plasmid="pExample-1"
     gene            complement(join(100..200,300..350))
                     /gene="fooA"
                     /locus_tag="LTG_0001"
     CDS             complement(join(100..200,300..350))
                     /gene="fooA"
                     /product="fused protein alpha"
                     /protein_id="PROT_0001"
                     /translation="MSTNPKPQRKTK"
     CDS             700..>950
                     /gene="barB"
                     /product="hypothetical protein"
                     /protein_id="PROT_0002"
ORIGIN
        1 acgatcgatc gatcgatcga tcgatcgatc gatcgatcga tcgatcgatc gatcgatcga
       61 tcgatcgatc gatcgatcga tcgatcgatc gatcgatcga tcgatcgatc gatcgatcga
//
'
tf <- tempfile(fileext = ".gbk")
writeLines(ex, tf)
gb <- read_gbk(tf)

testthat::test_that("FASTA export works", {
  gbk <- read_gbk(tf)
  tmp <- tempfile(fileext = ".fna")
  write_gbk_fasta(gbk, tmp)
  lines <- readLines(tmp)
  testthat::expect_true(any(grepl("^>", lines)))  # has header
  testthat::expect_true(any(grepl("^[ACGTN]+$", lines)))  # has sequence
})

testthat::test_that("GFF3 export works", {
  gbk <- read_gbk(tf)
  tmp <- tempfile(fileext = ".gff3")
  write_gbk_gff3(gbk, tmp)
  lines <- readLines(tmp)
  testthat::expect_equal(lines[1], "##gff-version 3")
  testthat::expect_true(any(grepl("CDS", lines)))
})

# ==== Using parser ====
library(tidyverse)
here::here()
#setwd("/home/william-ackerman/Desktop/Genbank_GBK_FILES")
list.files(patter = "gbk")

gbk <- read_gbk("Ecoli_NIST0056_PGAP.gbk")

# Inspect
gbk[[1]]$metadata
head(gbk[[1]]$features[, c("type","start","end","strand","location_type","gene","product")])
gbk[[2]]$features$ranges[[1]]  # ranges for the joined/complement CDS
nchar(gbk[[2]]$sequence)     

map_dfr(gbk, ~ .x[["metadata"]])
gbk[[1]]$metadata$definition
map_dfr(gbk, ~ .x[["features"]][["ranges"]][[1]])

# 1. Export sequences to FASTA
write_gbk_fasta(gbk, "my_genome.fna")

# 2. Export features to GFF3
write_gbk_gff3(gbk, "my_genome.gff3")

# Optional: export only CDS features
write_gbk_gff3(gbk, "cds_only.gff3", include_types = "CDS")


ecoli_atcc <- read_gbk("Ecoli_NIST0056_PGAP.gbk")

all_feats <- bind_gbk_features(ecoli_atcc )

ecoli_atcc_gene_list <- all_feats %>% select(gene) %>% unique()

ecoli_atcc_gene_list %>%
  filter(str_detect(gene, "^cus") | str_detect(gene, "^sil"))
# gene
# 1 cusS
# 2 cusR
# 3 cusC
# 4 cusF
# 5 cusB
# 6 cusA

ecoli_nist <- read_gbk("Ecoli_NIST0056.gbk")


ecoli_nist[[1]]$sequence
ecoli_nist[[1]]$features %>% filter(gene == "ampC")
str_sub(ecoli_nist[[1]]$sequence, 2607270, 2608403)

library(Biostrings)

feats_nist <- bind_gbk_features(ecoli_nist  ) 
feats_nist %>% select(gene) %>% unique() %>%
   filter(str_detect(gene, "^cus") | str_detect(gene, "^sil"))

feats_nist %>%
  filter(str_detect(gene, "^ampC"))

Eco_nist_df <- feats_nist %>% select(gene, product)

Eco_nist_df %>% filter(str_detect(gene, "amp"))

Eco_atcc_df <- all_feats %>% select(gene, product)

write.csv(Eco_atcc_df, file = "Eco_atcc_df.csv")
write.table(Eco_atcc_df, file = "Eco_atcc_df.txt", sep = "\t", row.names = FALSE, quote = FALSE)




write.csv(Eco_nist_df, file = "Eco_nist_df.csv" )
write.table(Eco_nist_df, file = "Eco_nist_df.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# gene
# 1 cusS
# 2 cusR
# 3 cusC
# 4 cusF
# 5 cusB
# 6 cusA

ecoli_nist_p1 <- read_gbk("Ecoli_NIST0056 plasmid unnamed1.gbk")
feats2 <- bind_gbk_features(ecoli_nist_p1  ) 
feats2 %>% select(gene) %>% unique() %>%
  filter(str_detect(gene, "cus|sil|cop"))


feats2 %>% select(gene, product)

ecoli_nist_p2 <- read_gbk("Ecoli_NIST0056 plasmid unnamed2.gbk")
feats3 <- bind_gbk_features(ecoli_nist_p2  ) 
feats3 %>% select(gene) %>% unique() %>%
  filter(str_detect(gene, "^cus") | str_detect(gene, "^sil"))

ecoli_nist_p3 <- read_gbk("Ecoli_NIST0056 plasmid unnamed3.gbk")
feats4 <- bind_gbk_features(ecoli_nist_p3  ) 
feats4 %>% select(gene) %>% unique() %>%
  filter(str_detect(gene, "^cus") | str_detect(gene, "^sil"))

ecoli_nist_p4 <- read_gbk("Ecoli_NIST0056 plasmid unnamed4.gbk")
feats5 <- bind_gbk_features(ecoli_nist_p4  ) 
feats5 %>% select(gene) %>% unique() %>%
  filter(str_detect(gene, "^cus") | str_detect(gene, "^sil"))

colnames(feats4)
feats3 %>% filter(type == "CDS") %>% select(gene, product, locus)
feats4 %>% filter(type == "CDS") %>% select(gene, product, locus)
feats5 %>% filter(type == "CDS") %>% select(gene, product, locus)

feats5 %>% filter(product == "N-6 DNA methylase")

Enteroc_ATCC13047 <- read_gbk("Enteroc_ATCC13047.gbk" )
Entfeats <- bind_gbk_features(Enteroc_ATCC13047) 
Entfeats %>% select(gene) %>% unique() %>% 
  filter(str_detect(gene, "^cus") | str_detect(gene, "^sil"))

# gene
# 1 cusS
# 2 cusR
# 3 cusC
# 4 cusF
# 5 cusB
# 6 cusA
# 7 silE

Enteroc_ATCC13047_p1 <- read_gbk("Enteroc_ATCC13047_plasmid_pECL_A.gbk" )
Entfeatsp1 <- bind_gbk_features(Enteroc_ATCC13047_p1) 
Entfeatsp1 %>% select(gene) %>% unique() %>% 
  filter(str_detect(gene, "^cus") | str_detect(gene, "^sil"))
# gene
# 1  silP
# 2  silA
# 3  silB
# 4  cusF
# 5  silC
# 6  silR
# 7  silS
# 8  silE
# 9  cusB
# 10 cusC
# 11 cusR
# 12 cusS

Enteroc_ATCC13047_p2 <- read_gbk("Enteroc_ATCC13047_plasmid_pECL_B.gbk" )
Entfeatsp2 <- bind_gbk_features(Enteroc_ATCC13047_p2) 
Entfeatsp2 %>% select(gene) %>% unique() %>% 
  filter(str_detect(gene, "^cus") | str_detect(gene, "^sil"))






list.files(pattern = "gbk")
Kpneumo_ATCC13883 <- read_gbk("Kpneumo_ATCC13883.gbk")
Kpnfeats <- bind_gbk_features(Kpneumo_ATCC13883)

Kpnfeats %>% select(product) %>%
  filter(str_detect(product, "efflux")) 
Kpnfeats %>% select(product) %>%
  filter(str_detect(product, "Cus") | str_detect(product, "Sil")) 
# product
# 1     CusA/CzcA family heavy metal efflux RND transporter
# 2 copper/silver efflux system outer membrane protein CusC
# 3                                      sensor kinase CusS        

Kpneumo_ATCC13883_p1 <- read_gbk("Kpneumo_ATCC13883_plasmid_pKPHS1.gbk")
Kpnfeatsp1 <- bind_gbk_features(Kpneumo_ATCC13883_p1)
Kpnfeatsp1 %>% colnames()

Kpnfeatsp1 %>% select(gene, product) %>% # phage and porphyrin # phage-related genes
  filter(str_detect(product, "Sil") | str_detect(product, "pump") | str_detect(product, "protein"))  

Kpneumo_ATCC13883_p2 <- read_gbk("Kpneumo_ATCC13883_plasmid_pKPHS2.gbk") # F-like conjugative plasmid (tra/trb genes).
Kpnfeatsp2 <- bind_gbk_features(Kpneumo_ATCC13883_p2)
Kpnfeatsp2 %>% select(gene, product) %>% 
  filter(str_detect(product, "Sil") | str_detect(product, "pump") | str_detect(product, "protein")) 

Kpneumo_ATCC13883_p3 <- read_gbk("Kpneumo_ATCC13883_plasmid_pKPHS3.gbk")
Kpnfeatsp3 <- bind_gbk_features(Kpneumo_ATCC13883_p3)
Kpnfeatsp3 %>% select(gene, product) %>% 
  filter(str_detect(product, "Sil|sil")) 

Kpneumo_ATCC13883_p4 <- read_gbk("Kpneumo_ATCC13883_plasmid_pKPHS4.gbk")
Kpnfeatsp4 <- bind_gbk_features(Kpneumo_ATCC13883_p4)
Kpnfeatsp4 %>% select(gene, product) %>% 
  filter(str_detect(product, "Sil") | str_detect(product, "pump") | str_detect(product, "protein")) 

Kpneumo_ATCC13883_p5 <- read_gbk("Kpneumo_ATCC13883_plasmid_pKPHS5.gbk")
Kpnfeatsp5 <- bind_gbk_features(Kpneumo_ATCC13883_p5)
Kpnfeatsp5 %>% select(gene, product) %>% 
  filter(str_detect(product, "Sil") | str_detect(product, "pump") | str_detect(product, "protein")) %>% 
  filter(str_detect(product, "Sil|sil")) 

Kpneumo_ATCC13883_p6 <- read_gbk("Kpneumo_ATCC13883_plasmid_pKPHS6.gbk")
Kpnfeatsp6 <- bind_gbk_features(Kpneumo_ATCC13883_p6)
Kpnfeatsp6 %>% select(gene, product) %>% 
  filter(str_detect(product, "Sil") | str_detect(product, "pump") | str_detect(product, "protein")) %>% 
  filter(str_detect(product, "Sil|sil")) 





list.files(pattern = "gbk")
Saureus_ATCC25923 <- read_gbk("Saureus_ATCC25923.gbk")
Saufeats <- bind_gbk_features(Saureus_ATCC25923)

Saufeats %>% select(product) %>% 
  filter(str_detect(product, "pump")) 
Saufeats %>% select(product) %>%
  filter(str_detect(product, "Cop|copper|Cso"))  
# product
# 1 CsoR family transcriptional regulator
# 2                 copper chaperone CopZ

Saureus_ATCC25923_p1 <- read_gbk("Saureus_ATCC25923_plasmidS1945.gbk") 
Saufeatsp1 <- bind_gbk_features(Saureus_ATCC25923_p1)
Saufeatsp1 %>% select(product) %>% 
  filter(str_detect(product, "pump")) 
Saufeatsp1 %>% select(product) %>% 
  filter(str_detect(product, "Cop") | str_detect(product, "cop")) 
