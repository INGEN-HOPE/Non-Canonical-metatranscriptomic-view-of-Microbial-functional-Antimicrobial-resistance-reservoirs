library(dplyr)
library(tidyr)
library(readr)
library(stringr)
host_files = list.files("~/Downloads/VARIANTS/final/read_ids_amr_host_sm/extracted_fasta_new/tsv_format_blast_result_per_read/blast_results_output_raw_with_NA where NA is no similarity found/all",pattern = ".tsv$", full.names = TRUE)

host_df = lapply(host_files, \(x) 
                 read_tsv(x))

sample_names <- sub("_all\\.tsv$", "", basename(host_files))
names(host_df) <- sample_names

head(host_files)


process_host_sample <- function(df) {
  df_clean <- df %>%
    mutate(QueryCover = as.numeric(gsub("%", "", QueryCover))) %>%
    distinct(ReadName, ScientificName, TotalScore, QueryCover, Evalue, PercIdent, .keep_all = TRUE)

  df_missing <- df_clean %>%
    group_by(ReadName) %>%
    filter(all(is.na(TotalScore))) %>%
    reframe(
      host_name = first(ScientificName),
      n_matches = NA_integer_,
      filter_status = "No hits found",
      across(c(TotalScore, QueryCover, PercIdent, Evalue, Description, Accession), \(x) first(x))
    )
  
  #Weighted Selection
  all_reads_processed <- df_clean %>%
    filter(!ReadName %in% df_missing$ReadName) %>%
    # Create a quality metric to balance Coverage and Identity
    mutate(quality_metric = (QueryCover * PercIdent) / 100) %>%
    group_by(ReadName) %>%
    #Highest quality balance first
    slice_max(order_by = quality_metric, n = 1, with_ties = TRUE) %>%
    # Lowest Evalue as tie breaker
    slice_min(order_by = Evalue, n = 1, with_ties = TRUE) %>%
    #Highest Total Score as final tie breaker
    slice_max(order_by = TotalScore, n = 1, with_ties = TRUE) %>%
    reframe(
      host_name = paste(unique(na.omit(ScientificName)), collapse = ", "),
      n_matches = n(),
      across(c(TotalScore, QueryCover, PercIdent, Evalue, Description, 
               CommonName, Taxid, Accession, MaxScore, AccLen), \(x) first(x))
    ) %>%
    mutate(filter_status = "Unfiltered Top Hit")
  
  #Re-combine
  overall_final <- bind_rows(all_reads_processed, df_missing)
  
  #Final Filter
  df_result <- all_reads_processed %>%
    filter(Evalue < 1e-05, PercIdent >= 90, QueryCover >= 50) %>%
    mutate(filter_status = "Passed quality filters")
  
  return(list(overall = overall_final, filtered = df_result))
}


# Apply to all samples
host_results_list_new <- lapply(host_df, process_host_sample)
names(host_results_list_new) <- sample_names

# Create separate lists for overall and filtered
host_overall_list_new <- lapply(host_results_list_new, function(x) x$overall)
host_filtered_list_new <- lapply(host_results_list_new, function(x) x$filtered)

# Create dataframes
host_overall_df_new <- bind_rows(host_overall_list_new, .id = "sample_id")
host_filtered_df_new <- bind_rows(host_filtered_list_new, .id = "sample_id")

dim(host_overall_df_new)
dim(host_filtered_df_new)

write.csv(host_overall_df_new, "all_host_reads_blast_output_new.csv", row.names = FALSE)
write.csv(host_filtered_df_new, "filter_host_reads_blast_output_new.csv", row.names = FALSE)



# #Sanity Check
# # Create key columns in both
# overall_keys <- host_overall_df_new %>% 
#   select(sample_id, ReadName, host_name) %>%
#   distinct()
# 
# filtered_keys <- host_filtered_df_new %>%
#   select(sample_id, ReadName, host_name) %>%
#   distinct()
# 
# #Rows in filtered but NOT in overall
# in_filtered_not_overall <- anti_join(filtered_keys, overall_keys, 
#                                      by = c("sample_id", "ReadName", "host_name"))
# cat("In filtered but NOT in overall:", nrow(in_filtered_not_overall), "\n")
# 
# #Rows in overall but NOT in filtered
# in_overall_not_filtered <- anti_join(overall_keys, filtered_keys, 
#                                      by = c("sample_id", "ReadName", "host_name"))
# cat("In overall but NOT in filtered:", nrow(in_overall_not_filtered), "\n")
# 
# #Check if same ReadName has DIFFERENT host_name between the two dfs
# merged_check_n <- inner_join(overall_keys, filtered_keys, 
#                              by = c("sample_id", "ReadName"),
#                              suffix = c("_overall", "_filtered")) %>%
#   filter(host_name_overall != host_name_filtered)
# 
# cat("ReadNames with DIFFERENT host_name between overall and filtered: should ideally be zero", nrow(merged_check_n), "\n")


# Filtering only abundant microbes with >0.1% RA
abund <- read.csv("~/Downloads/VARIANTS/FOR_FILTERING ABUNDANT TAMS.csv")
meta <- read.csv("~/Downloads/VARIANTS/meta.csv")

# Get the lists of allowed hosts for each condition
negative_hosts <- abund$Negative
nonsevere_hosts <- abund$Nonsevere
severe_hosts <- abund$Severe

# Create a mapping of sample to condition
sample_condition_map <- setNames(meta$condition, meta$sample_name)

# Function to filter host names based on allowed list
filter_host_names <- function(host_name, allowed_hosts) {
  if (is.na(host_name)) {
    return(NA_character_)
  }
  
  # Split by comma
  hosts <- str_split(host_name, ",\\s*")[[1]]
  
  # Keep only hosts that are in the allowed list
  filtered <- hosts[hosts %in% allowed_hosts]
  
  # Return comma-separated or NA if none match
  if (length(filtered) == 0) {
    return(NA_character_)
  } else {
    return(paste(filtered, collapse = ", "))
  }
}

# Add condition to host_filtered_df
host_filtered_df_n <- host_filtered_df_new %>%
  mutate(
    condition = sample_condition_map[sample_id]
  )

# Apply filtering based on condition
host_filtered_df_n <- host_filtered_df_n %>%
  rowwise() %>%
  mutate(
    host_name_to_assign = case_when(
      condition == "Negative" ~ filter_host_names(host_name, negative_hosts),
      condition == "Non-Severe" ~ filter_host_names(host_name, nonsevere_hosts),
      condition == "Severe" ~ filter_host_names(host_name, severe_hosts),
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# Reorder columns
host_filtered_df_n <- host_filtered_df_n %>%
  select(sample_id, ReadName, host_name, n_matches, host_name_to_assign, MaxScore, everything(), condition, -filter_status)

write.csv(host_filtered_df_n, "filtered_host_abundant_TAMs.csv", row.names = FALSE)

