
setwd("D:\\Desktop_Laptop_Xfer\\Nano-Enabled_Antimicrobials_Working_Group\\Stabryla_Bacterial_Evolution\\AMR_Motility\\AMR_Motility_breseq_Results_WEA")

#library(xml2)

library(rvest)
library(tidyverse)
library(janitor)

# Look for html file in working directory ----

files <- list.files(pattern=".html")


# Make S3 list object converting mutations in html to tibbles ----

mutation_table_list <- list()

for(i in 1:length(files)){
  
  entry_name <- str_remove(files[i], ".html") 
  
  table_list <- read_html(files[i]) %>% 
    html_nodes("table") %>% .[2] %>% 
    html_table(fill = TRUE) 
  
  mut_table_entry <- table_list[[1]] %>% 
    row_to_names(row_number = 1)
  
  mutation_table_list[[entry_name]] <- mut_table_entry  

}

# Functions for data wrangling ----
convert_freq_to_numeric <- function(df) {
  df %>% mutate(freq_num = str_remove_all(freq, "%") %>% as.numeric())
}

filter_mut_freq_x_pct <- function(df, x_pct = 20) {
  
  # check that any supplied user input for x_pct is valid
  if (!is.numeric(x_pct) || x_pct < 1 || x_pct > 100) {
    stop("Error: x_pct must be a numeric value between 1 and 100.")
  }
  
  # Make sure the data frame contains freq_num column
  # If absent, call convert_freq_to_numeric() function
  if (!"freq_num" %in% colnames(df)) {
    df <- convert_freq_to_numeric(df)
  }
  
    df %>% filter(freq_num >= x_pct )
}

make_mut_uid <- function(df){
  df %>%  mutate(mutation_uid = paste0(position, " ", mutation))
}

list_all_mutations <- function(df){

  df_with_uids <- lapply(working_data, make_mut_uid)
  
  all_mutation_uids <- unique(unlist(lapply(df_with_uids, 
                                            function(df) df$mutation_uid)))
  
  return(all_mutation_uids)
}

# Are there fixed mutations present in the entire dataset?
working_data <-
  lapply(mutation_table_list, filter_mut_freq_x_pct, x_pct = 20)

all_mutations <- list_all_mutations(working_data)

# Find the intersection of mutation_uids across all dataframes

df_with_uids <- lapply(working_data, make_mut_uid)

names(df_with_uids)

# Subtract the ancestral in 3-D0 and 4-D0 as well as U1-5-D10
ancestral_mutation_uids <- 
  Reduce(intersect, lapply(df_with_uids[c(1,8)], 
                           function(df) df$mutation_uid))

mutations_wo_sel_press_uids <- 
  Reduce(intersect, lapply(df_with_uids[26:30], 
                           function(df) df$mutation_uid))

common_mutation_uids <- 
  unique(c(ancestral_mutation_uids, 
           mutations_wo_sel_press_uids ))

# Filter common mutations
res <- 
  lapply(df_with_uids[16:30], function(df) 
  { df %>% filter(!df$mutation_uid %in% common_mutation_uids)} ) 

# Collate and export results ----
# Function to modify each data frame
modify_df <- function(df, name) {
  df <- df %>%
    add_column(name = name, .before = 1) %>%
    select(-freq_num, -mutation_uid)
  return(df)
}

modified_res <- lapply(names(res), function(name) modify_df(res[[name]], name))

# Bind the modified data frames into a master data frame
master_df <- bind_rows(modified_res)

# Export to spreadsheet
library(openxlsx)
write.xlsx(master_df, 
          file = "Spreadsheet_of_breseq_results-WEA-REDO.xlsx",
          fileEncoding = "UTF-8")
      



