
library(readr)
library(dplyr)
library(cli)
library(reticulate)

reticulate::use_condaenv("base")
Sys.setenv(BLASTDB = "/home/william-ackerman/Desktop/Link to Desktop/tmp_ncbi_blast/")

blastp_capture <- function(query_faa,
                           db = "swissprot",
                           dbdir = Sys.getenv("BLASTDB"),
                           evalue = 1e-10,
                           threads = parallel::detectCores(),
                           fields = c("qseqid","sacc","stitle","pident","length",
                                      "qcovs","bitscore","evalue","staxids"),
                           validate_db = TRUE) {
  # 1) Locate blastp
  blastp <- Sys.which("blastp")
  if (blastp == "") cli::cli_abort("blastp not found on PATH. Install BLAST+ or adjust PATH.")
  
  # 2) Resolve DB path
  db_path <- if (dbdir != "" && !grepl("^/", db)) file.path(dbdir, db) else db
  
  # 3) Optional DB sanity check (no temp files; capture to memory)
  if (validate_db) {
    info <- system2("blastdbcmd", args = c("-db", db, "-info"),
                    stdout = TRUE, stderr = "")
    if (length(info) == 0 || any(grepl("\\bError\\b", info))) {
      cli::cli_abort(c("Could not open BLAST DB.",
                       ">" = paste(info, collapse = "\n"),
                       i = "Set Sys.setenv(BLASTDB = '~/blastdb') or pass an absolute -db path."))
    }
  }
  
  # 4) Compose outfmt 6 (tabular) string
  fmt <- paste("6", paste(fields, collapse = " "))
  fmt_q <- shQuote(fmt) 
  
  # 5) Build args (keep SEG + composition-based stats)
  args <- c("-query", query_faa,
            "-db", db,
            "-evalue", format(evalue, scientific = TRUE),
            "-seg", "yes",
            "-comp_based_stats", "2",
            "-num_threads", threads,
            "-outfmt", fmt_q)
  
  cli::cli_alert_info("Running: {blastp} {paste(args, collapse=' ')}")
  
  # 6) Capture stdout in memory; print stderr to console (no files)
  out <- system2(command = blastp, args = args, stdout = TRUE, stderr = "")
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cli::cli_abort(c("blastp exited with status {status}.",
                     i = "Inspect console stderr for details."))
  }
  
  # 7) Parse in-line (no temp files)
  df <- readr::read_tsv(I(out), col_names = fields, show_col_types = FALSE)
  
  # 8) Tidy and rank
  df |>
    mutate(across(c(pident, qcovs, bitscore, evalue), suppressWarnings(as.numeric))) |>
    arrange(evalue, desc(bitscore), desc(qcovs))
}


query_faa <- tempfile(fileext = ".faa")
writeLines(c(
  ">query_seq",
  aa_seq), query_faa)

hits <- blastp_capture(query_faa)
print(hits, n = 10)








