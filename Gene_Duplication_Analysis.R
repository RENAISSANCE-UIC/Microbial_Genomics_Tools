setwd("/home/william-ackerman/Desktop/cnMOPS/star")

library(tidyverse)
library(data.table)

#' Read depth file
#' 
#' @param filepath Path to depth file (tab-separated: chrom, pos, depth)
#' @return data.table with depth information
read_depth_file <- function(filepath) {
  fread(filepath, col.names = c("chrom", "pos", "depth"))
}

#' Read and prepare gene reference file
#' 
#' @param filepath Path to BED file
#' @return tibble with chrom, start, end, gene columns
read_gene_reference <- function(filepath) {
  fread(filepath) %>%
    select(1, 2, 3, 10) %>%
    set_names(c("chrom", "start", "end", "gene_info"))
}

#' Parse gene attributes from info string
#' 
#' @param info_string String containing gene attributes (e.g., "Name=ampC;ID=123;...")
#' @return tibble with parsed attributes
parse_gene_info <- function(info_string) {
  tibble(gene_info = info_string) %>%
    mutate(
      name = str_extract(gene_info, "Name=([^;]+)", group = 1),
      id = str_extract(gene_info, "ID=([^;]+)", group = 1),
      alias = str_extract(gene_info, "Alias=([^;]+)", group = 1),
      note = str_extract(gene_info, "Note=([^;]+)", group = 1)
    )
}

#' Calculate mean coverage for a genomic region
#' 
#' @param depth_dt data.table with depth information
#' @param chrom_name Chromosome name
#' @param start_pos Start position
#' @param end_pos End position
#' @return Mean coverage value
calculate_region_coverage <- function(depth_dt, chrom_name, start_pos, end_pos) {
  depth_dt %>%
    filter(chrom == chrom_name, pos >= start_pos, pos <= end_pos) %>%
    pull(depth) %>%
    mean(na.rm = TRUE)
}

#' Calculate coverage for all genes
#' 
#' @param genes_df Dataframe with gene coordinates
#' @param depth_dt data.table with depth information
#' @return tibble with mean coverage for each gene
calculate_gene_coverage <- function(genes_df, depth_dt) {
  genes_df %>%
    rowwise() %>%
    mutate(
      mean_cov = calculate_region_coverage(depth_dt, chrom, start, end)
    ) %>%
    ungroup()
}

#' Normalize coverage and identify duplications
#' 
#' @param coverage_df Dataframe with mean_cov column
#' @param threshold Threshold for duplication detection (default: 1.8)
#' @return tibble with normalized coverage and duplication status
normalize_and_flag_duplications <- function(coverage_df, threshold = 1.8) {
  median_cov <- median(coverage_df$mean_cov, na.rm = TRUE)
  
  coverage_df %>%
    mutate(
      median_cov = median_cov,
      norm_cov = mean_cov / median_cov,
      is_duplicated = norm_cov >= threshold  # Fixed: >= instead of >
    )
}

#' Search for genes by name pattern
#' 
#' @param results_df Results dataframe
#' @param pattern Search pattern (regex)
#' @param field Field to search in (default: "name")
#' @return Filtered tibble
search_genes <- function(results_df, pattern, field = "name") {
  results_df %>%
    filter(str_detect(.data[[field]], regex(pattern, ignore_case = TRUE)))
}

#' Get duplicated genes
#' 
#' @param results_df Results dataframe with is_duplicated column
#' @return tibble of duplicated genes sorted by normalized coverage
get_duplicated_genes <- function(results_df) {
  results_df %>%
    filter(is_duplicated) %>%
    arrange(desc(norm_cov))
}

#' Main pipeline to analyze gene coverage
#' 
#' @param depth_file Path to depth file
#' @param reference_file Path to reference BED file
#' @param duplication_threshold Threshold for calling duplications (default: 1.8)
#' @return List containing results and summary statistics
analyze_gene_coverage <- function(depth_file, reference_file, duplication_threshold = 1.8) {
  # Read data
  message("Reading depth file...")
  depth_dt <- read_depth_file(depth_file)
  
  message("Reading reference file...")
  genes_df <- read_gene_reference(reference_file)
  
  # Calculate coverage
  message("Calculating gene coverage...")
  coverage_df <- calculate_gene_coverage(genes_df, depth_dt)
  
  # Normalize and flag
  message("Normalizing coverage...")
  results <- normalize_and_flag_duplications(coverage_df, duplication_threshold)
  
  # Parse gene information
  message("Parsing gene information...")
  results <- results %>%
    bind_cols(parse_gene_info(.$gene_info))
  
  # Summary statistics
  n_duplicated <- sum(results$is_duplicated, na.rm = TRUE)
  message(sprintf("Found %d duplicated genes (>= %.1fx)", n_duplicated, duplication_threshold))
  
  list(
    results = results,
    n_duplicated = n_duplicated,
    median_coverage = unique(results$median_cov)
  )
}

#' Batch process all depth files in a directory
#' 
#' @param depth_dir Directory containing depth files
#' @param reference_file Path to reference BED file
#' @param output_dir Directory for output files (default: same as depth_dir)
#' @param duplication_threshold Threshold for calling duplications (default: 1.8)
#' @param pattern Pattern to match depth files (default: "_depth.txt")
#' @param save_format Format to save results: "rdata", "csv", or "both" (default: "rdata")
#' @return Invisible list of all results
batch_process_coverage <- function(depth_dir, 
                                   reference_file,
                                   output_dir = depth_dir,
                                   duplication_threshold = 1.8,
                                   pattern = "_depth\\.txt$",
                                   save_format = "rdata") {
  
  # Get all depth files
  depth_files <- list.files(depth_dir, pattern = pattern, full.names = TRUE)
  
  if (length(depth_files) == 0) {
    stop("No depth files found matching pattern: ", pattern)
  }
  
  message(sprintf("Found %d depth files to process", length(depth_files)))
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  
  # Process each file
  all_results <- list()
  
  for (i in seq_along(depth_files)) {
    depth_file <- depth_files[i]
    
    # Extract sample prefix
    sample_prefix <- basename(depth_file) %>%
      str_remove(pattern)
    
    message(sprintf("\n[%d/%d] Processing: %s", i, length(depth_files), sample_prefix))
    
    # Try to process the file
    tryCatch({
      # Run analysis
      results <- analyze_gene_coverage(
        depth_file = depth_file,
        reference_file = reference_file,
        duplication_threshold = duplication_threshold
      )
      
      # Add sample name to results
      results$sample_name <- sample_prefix
      results$results <- results$results %>%
        mutate(sample = sample_prefix, .before = 1)
      
      # Save results
      if (save_format %in% c("rdata", "both")) {
        output_file <- file.path(output_dir, paste0(sample_prefix, ".RData"))
        save(results, file = output_file)
        message(sprintf("  ✓ Saved RData: %s", basename(output_file)))
      }
      
      if (save_format %in% c("csv", "both")) {
        csv_file <- file.path(output_dir, paste0(sample_prefix, "_results.csv"))
        write_csv(results$results, csv_file)
        message(sprintf("  ✓ Saved CSV: %s", basename(csv_file)))
        
        # Also save duplicated genes separately
        dup_file <- file.path(output_dir, paste0(sample_prefix, "_duplicated.csv"))
        get_duplicated_genes(results$results) %>%
          write_csv(dup_file)
        message(sprintf("  ✓ Saved duplications: %s", basename(dup_file)))
      }
      
      # Store in list
      all_results[[sample_prefix]] <- results
      
    }, error = function(e) {
      message(sprintf("  ✗ ERROR processing %s: %s", sample_prefix, e$message))
      all_results[[sample_prefix]] <- NULL
    })
  }
  
  message(sprintf("\n✓ Batch processing complete! Processed %d/%d files successfully", 
                  length(all_results), length(depth_files)))
  
  invisible(all_results)
}

#' Create summary table across all samples
#' 
#' @param results_list List of results from batch processing
#' @return tibble with summary statistics per sample
create_batch_summary <- function(results_list) {
  map_df(results_list, function(res) {
    tibble(
      sample = res$sample_name,
      n_genes = nrow(res$results),
      n_duplicated = res$n_duplicated,
      pct_duplicated = round(100 * n_duplicated / n_genes, 2),
      median_coverage = res$median_coverage,
      mean_coverage = mean(res$results$mean_cov, na.rm = TRUE),
      max_norm_cov = max(res$results$norm_cov, na.rm = TRUE)
    )
  })
}

#' Create combined duplication matrix across samples
#' 
#' @param results_list List of results from batch processing
#' @return Wide-format tibble with duplication status for each gene x sample
create_duplication_matrix <- function(results_list) {
  map_df(results_list, function(res) {
    res$results %>%
      select(chrom, start, end, name, norm_cov, is_duplicated)
  }, .id = "sample") %>%
    pivot_wider(
      id_cols = c(chrom, start, end, name),
      names_from = sample,
      values_from = norm_cov
    )
}

# ============================================================================
# USAGE EXAMPLES
# ============================================================================

setwd("/home/william-ackerman/Desktop/cnMOPS/read_depth")

list.files(pattern = "_depth.txt")
list.files(pattern = ".bed")

# Example 1: Full pipeline
all_results <- batch_process_coverage(
  depth_dir = ".",
  reference_file = "no_repeats.bed",
  output_dir = "./results_no_repeats",
  duplication_threshold = 1.8,
  save_format = "both"
)

# Create summary table
summary_table <- create_batch_summary(all_results)
summary_table %>%
  arrange(desc(n_duplicated)) %>%
  print(n = Inf)

# Create duplication matrix (genes x samples)
dup_matrix <- create_duplication_matrix(all_results)

write.csv(dup_matrix, file="Duplication_Matrix_Samtools_Depth_1.8X_REPEAT_MASKED.csv")


# Heatmap
dup_matrix$name <- make.unique(dup_matrix$name)

heatmap_data <- dup_matrix %>%
  select(name, `3-0_S1`:`A3-D10_S7`) %>%
  filter(!is.na(name)) %>% 
  column_to_rownames("name")

heatmap_matrix <- as.matrix(heatmap_data)

heatmap_matrix_clean <- heatmap_matrix[complete.cases(heatmap_matrix), ]

library(pheatmap)
pheatmap(heatmap_matrix_clean,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,  # Hide row names if too many
         fontsize_col = 8,
         color = colorRampPalette(c("blue", "white", "red"))(100))




# View all results sorted by coverage
results$results %>% 
  arrange(desc(norm_cov)) %>% 
  select(chrom, start, end, name, mean_cov, norm_cov, is_duplicated) %>%
  head(20)

# Example 2: Get duplicated genes
duplicated_genes <- get_duplicated_genes(results$results)
duplicated_genes %>%
  select(chrom, start, end, name, norm_cov, note) %>%
  print(n = Inf)

# Example 3: Search for specific gene
ampC_results <- search_genes(results$results, "ampC", field = "name")
ampC_results %>%
  select(chrom, start, end, name, norm_cov, is_duplicated)

# Example 4: Process multiple samples
process_multiple_samples <- function(depth_files, reference_file, threshold = 1.8) {
  map(depth_files, function(depth_file) {
    sample_name <- str_remove(basename(depth_file), "_depth\\.txt")
    message(sprintf("\nProcessing %s...", sample_name))
    
    results <- analyze_gene_coverage(depth_file, reference_file, threshold)
    
    results$results %>%
      mutate(sample = sample_name) %>%
      select(sample, everything())
  }) %>%
    bind_rows()
}

# Use with multiple files
# depth_files <- list.files(pattern = "*_depth.txt", full.names = TRUE)
# all_results <- process_multiple_samples(depth_files, "reference.bed")

# Example 5: Export results
export_results <- function(results_list, output_prefix) {
  # All results
  write_csv(
    results_list$results %>% select(chrom, start, end, name, id, mean_cov, norm_cov, is_duplicated, note),
    paste0(output_prefix, "_all_genes.csv")
  )
  
  # Duplicated genes only
  write_csv(
    get_duplicated_genes(results_list$results) %>% 
      select(chrom, start, end, name, id, norm_cov, note),
    paste0(output_prefix, "_duplicated_genes.csv")
  )
  
  message(sprintf("Results exported with prefix: %s", output_prefix))
}

# export_results(results, "A3-D10_S7")