
setwd("D:\\Desktop_Laptop_Xfer\\Nano-Enabled_Antimicrobials_Working_Group\\Stabryla_Bacterial_Evolution\\AMR_Motility\\AMR_Motility_breseq_Results_WEA")

#library(xml2)

pacman::p_load(
  rvest,
  tidyverse,
  janitor
)

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


names(working_data)
time_course_samples_low <- 
  c(
    "3-0_S1_mutations", "3-1_S2_mutations", "3-3_S3_mutations", 
    "3-4_S4_mutations", "3-7_S5_mutations", "3-8_S6_mutations",
    "3-10_S7_mutations")

time_course_samples_high <- 
  c(
    "4-0_S8_mutations", "4-1_S9_mutations", "4-3_S10_mutations",
    "4-4_S11_mutations", "4-5_S12_mutations", "4-6_S13_mutations",
    "4-8_S14_mutations", "4-10_S15_mutations"                          
  )

timecourse_data_low <- working_data[time_course_samples_low]

timecourse_data_high <- working_data[time_course_samples_high]


all_mutations <- list_all_mutations(timecourse_data)

# Find the intersection of mutation_uids across all dataframes

df_with_uids_low <- lapply(timecourse_data_low, make_mut_uid)

df_with_uids_high <- lapply(timecourse_data_high, make_mut_uid)

# Subtract the ancestral in 3-D0 and 4-DO
ancestral_mutation_uids_low <- 
  Reduce(intersect, 
         lapply(df_with_uids_low[1], 
                           function(df) df$mutation_uid))

ancestral_mutation_uids_high <- 
  Reduce(intersect, 
         lapply(df_with_uids_high[1], 
                function(df) df$mutation_uid))


# Filter common mutations
res_low <- 
  lapply(df_with_uids_low, function(df) 
  { df %>% filter(!df$mutation_uid %in% ancestral_mutation_uids_low)} ) 

res_high <- 
  lapply(df_with_uids_high, function(df) 
  { df %>% filter(!df$mutation_uid %in% ancestral_mutation_uids_high)} ) 


# Collate and export results ----
# Function to modify each data frame
modify_df <- function(df, name) {
  df <- df %>%
    add_column(name = name, .before = 1) %>%
    select(-freq_num, -mutation_uid)
  return(df)
}

modified_res_low <- lapply(names(res_low), 
                           function(name) modify_df(res_low[[name]], name))

modified_res_high <- lapply(names(res_high), 
                           function(name) modify_df(res_high[[name]], name))

# Bind the modified data frames into a master data frame
master_df_low <- bind_rows(modified_res_low)

master_df_high <- bind_rows(modified_res_high)

# Export to spreadsheet
library(openxlsx)
write.xlsx(master_df_low, 
          file = "Spreadsheet_of_breseq_results-WEA-low_AmpR.xlsx",
          fileEncoding = "UTF-8")

write.xlsx(master_df_high, 
           file = "Spreadsheet_of_breseq_results-WEA-high_AmpR.xlsx",
           fileEncoding = "UTF-8")      



