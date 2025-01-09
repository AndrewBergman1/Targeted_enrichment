library(tidyverse)
library(GGally)
library(DESeq2)
library(pheatmap)
library(edgeR)

###############################################################
### LOAD DATA
###############################################################

# Function to read sample files and assign identifiers
read_and_assign <- function(file_names, sample_type, folder, col_names) {
  map2_dfr(file_names, 1:length(file_names), ~ {
    df <- read_tsv(file.path(folder, .x), col_names = col_names, col_types = cols())
    
    # Convert all columns to character to ensure consistency
    df <- df %>% mutate(across(everything(), as.character))
    
    df %>%
      mutate(
        sample_id = .y,
        sample_type = sample_type
      )
  })
}

# Define folder and column names
folder <- "C:/Users/abergm/Downloads/combined_blast_files/"
standard_col_names <- paste0("V", 1:17)

# File names for fecal samples
fecal_file_names <- c(
  "6323_1x_S89_trinity_combined_megablast.txt",
  "6323_1x_S96_trinity_combined_megablast.txt",
  "6323_10e-1_S90_trinity_combined_megablast.txt",
  "6323_10e-2_S91_trinity_combined_megablast.txt",
  "6323_10e-3_S92_trinity_combined_megablast.txt",
  "6323_10e-4_S93_trinity_combined_megablast.txt",
  "6323_10e-5_S94_trinity_combined_megablast.txt",
  "6323_10e-6_S95_trinity_combined_megablast.txt"
)

# File names for spike-in samples
spike_file_names <- c(
  "25ng_unenriched_S1_trinity_combined_megablast.txt",
  "25ng_VERO_E6_enriched_S8_trinity_combined_megablast.txt",
  "10ng_enriched_S2_trinity_combined_megablast.txt",
  "1ng_enriched_S3_trinity_combined_megablast.txt",
  "10e-1ng_enriched_S4_trinity_combined_megablast.txt",
  "5e-2ng_enriched_S6_trinity_combined_megablast.txt",
  "2.5e-2ng_enriched_S7_trinity_combined_megablast.txt",
  "10e-2ng_enriched_S5_trinity_combined_megablast.txt"
)

# Read and assign sample types
fecal_samples <- read_and_assign(fecal_file_names, "fecal", folder, standard_col_names)
spike_samples <- read_and_assign(spike_file_names, "spike_in", folder, standard_col_names)

# Combine all samples into one data frame
all_samples <- bind_rows(fecal_samples, spike_samples)

###############################################################
### SAMPLE METADATA
###############################################################

# Sample metadata for fecal samples
fecal_sample_data <- tibble(
  ID = 1:8,
  Sample_Name = c(
    "6323_1x(1)",
    "6323_1x(2)",
    "10e-1ng_enriched",
    "10e-2ng_enriched",
    "10e-3ng_enriched",
    "10e-4ng_enriched",
    "10e-5ng_enriched",
    "10e-6ng_enriched"
  ),
  sample_type = "fecal"
)

# Sample metadata for spike-in samples
spike_in_sample_data <- tibble(
  ID = 1:8,
  Sample_Name = c(
    "25ng_control",
    "25ng",
    "10ng",
    "1ng",
    "0.1ng",
    "0.05ng",
    "0.025ng",
    "0.01ng"
  ),
  sample_type = "spike_in"
)

# Combine metadata
all_sample_data <- bind_rows(fecal_sample_data, spike_in_sample_data)

###############################################################
### DEFINE FUNCTIONS
###############################################################

# Function to calculate total hits per sample
calculate_hits_per_sample <- function(data, sample_metadata) {
  data %>%
    group_by(sample_id, sample_type) %>%
    summarise(hits = n(), .groups = "drop") %>%
    left_join(sample_metadata, by = c("sample_id" = "ID", "sample_type" = "sample_type")) %>%
    select(sample_id, Sample_Name, sample_type, hits)
}

# Function to calculate unique hits per sample
calculate_unique_hits <- function(data, sample_metadata) {
  data %>%
    group_by(sample_id, V17, sample_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(sample_id, sample_type) %>% 
    summarise(unique_hits = n(), .groups = "drop") %>%
    left_join(sample_metadata, by = c("sample_id" = "ID", "sample_type" = "sample_type")) %>%
    select(sample_id, Sample_Name, sample_type, unique_hits)
}

calculate_plasmid_hit_prevalence <- function(data, sample_metadata) {
  result <- data %>%
    # Filter rows where V17 contains "plasmid" or "vector" (case-insensitive)
    filter(str_detect(V17, regex("plasmid|vector", ignore_case = TRUE))) %>%
    
    # Group by sample_id and sample_type
    group_by(sample_id, sample_type, V17) %>%
    
    # Count occurrences of each V17 within the groups
    summarise(n = n(), .groups = "drop") %>%
    
    # Rename sample_id to sample for consistency
    rename(sample = sample_id) %>%
    
    # Join with the sample metadata to get Sample_Name and other details
    left_join(sample_metadata, by = c("sample" = "ID", "sample_type" = "sample_type")) %>%
    
    # Select and arrange the desired columns
    select(V17, sample, Sample_Name, n)
  
  return(result)
}

calculate_non_plasmid_hit_prevalence <- function(data, sample_metadata) {
  result <- data %>%
    # Filter rows where V17 does NOT contain "plasmid" or "vector" (case-insensitive)
    filter(!str_detect(V17, regex("plasmid|vector|res|resistance|antibiotic|arg|ARG", ignore_case = TRUE))) %>%
    
    # Group by sample_id, sample_type, and V17
    group_by(sample_id, sample_type, V17) %>%
    
    # Count the number of occurrences for each V17 within the groups
    summarise(n = n(), .groups = "drop") %>%
    
    # Rename sample_id to sample to match sample_metadata
    rename(sample = sample_id) %>%
    
    # Join with the sample metadata to get Sample_Name and other details
    left_join(sample_metadata, by = c("sample" = "ID", "sample_type" = "sample_type")) %>%
    
    # Select and arrange the desired columns
    select(V17, sample, Sample_Name, n)
  
  return(result)
}

calculate_fraction_hits <- function(hits_df, plasmid_df, non_plasmid_df, arg_hits_df, sample_metadata) {
  hits_df %>%
    # Join plasmid_hits by matching sample_id in hits_df to sample in plasmid_df
    left_join(
      plasmid_df %>% 
        group_by(sample) %>% 
        summarise(plasmid_hits = sum(n, na.rm = TRUE), .groups = "drop"),
      by = c("sample_id" = "sample")
    ) %>%
    # Join non_plasmid_hits similarly
    left_join(
      non_plasmid_df %>% 
        group_by(sample) %>% 
        summarise(non_plasmid_hits = sum(n, na.rm = TRUE), .groups = "drop"),
      by = c("sample_id" = "sample")
    ) %>%
    # Join arg_hits similarly
    left_join(
      arg_hits_df %>% 
        group_by(sample) %>% 
        summarise(arg_hits = sum(n, na.rm = TRUE), .groups = "drop"),
      by = c("sample_id" = "sample")
    ) %>% 
    # Join with sample_metadata to get Sample_Name and sample_type
    left_join(
      sample_metadata, 
      by = c("sample_id" = "ID", "sample_type" = "sample_type")
    ) %>%
    # Calculate fractions
    mutate(
      frac_plasmids = plasmid_hits / hits,
      frac_non_plasmids = non_plasmid_hits / hits,
      frac_args = arg_hits / hits
    ) %>%
    # Select and arrange desired columns
    select(sample_id, Sample_Name.x, sample_type, frac_plasmids, frac_non_plasmids, frac_args, arg_hits, plasmid_hits, non_plasmid_hits) %>% 
    rename(Sample_Name=Sample_Name.x)
}
# Function to calculate isoform statistics
calculate_isoform_statistics <- function(data) {
  data %>%
    mutate(
      Isoform = str_extract(V1, "_i(\\d+)") %>% 
        str_replace("_i", "") %>% 
        as.integer()
    ) %>%
    group_by(V17) %>%
    summarise(n_isoforms = n_distinct(Isoform), .groups = "drop") %>%
    summarise(
      mean_isoform_count = mean(n_isoforms, na.rm = TRUE),
      min_isoform_count = min(n_isoforms, na.rm = TRUE),
      max_isoform_count = max(n_isoforms, na.rm = TRUE),
      median_isoform_count = median(n_isoforms, na.rm = TRUE),
      sd_isoform_count = sd(n_isoforms, na.rm = TRUE),
      IQR_isoform_count = IQR(n_isoforms, na.rm = TRUE)
    )
}

plot_plasmid_fraction <- function(frac_df, fecal = F) {
  
  if (fecal == T) {
    desired_order = c(
      "6323_1x(1)",
      "6323_1x(2)",
      "10e-1ng_enriched",
      "10e-2ng_enriched",
      "10e-3ng_enriched",
      "10e-4ng_enriched",
      "10e-5ng_enriched",
      "10e-6ng_enriched"
    )
  }
  else {
    desired_order = c("25ng",
                      "10ng",
                      "1ng",
                      "0.1ng",
                      "0.05ng",
                      "0.025ng",
                      "0.01ng"
    )
  }
  
  frac_df <- frac_df %>% 
    mutate(Sample_Name = factor(Sample_Name, levels = desired_order),
           frac_plasmids = frac_plasmids*post_enr_c,
           frac_non_plasmids = frac_non_plasmids*post_enr_c,
           frac_arg = frac_args)
  
  
  frac_df_long <- frac_df %>%
    select(Sample_Name , frac_plasmids, frac_non_plasmids, frac_args, post_enr_c) %>%
    pivot_longer(
      cols = c(starts_with("frac_")),
      names_to = "hit_type",
      values_to = "fraction"
    ) %>%
    mutate(hit_type = recode(hit_type,
                             frac_plasmids = "Plasmid/Vector",
                             frac_non_plasmids = "Non-Plasmid/Vector/ARG",
                             frac_args = "ARGs")) %>% 
    mutate(hit_type = factor(hit_type,
                             levels = c("Non-Plasmid/Vector/ARG", "ARGs", "Plasmid/Vector")))
  
  
  
  p <- ggplot(frac_df_long, aes(x = Sample_Name, y = fraction, fill = hit_type)) +
    geom_col(position = "stack") +
    theme_minimal(base_size = 14) +
    scale_fill_manual(values = c("Plasmid/Vector" = "steelblue",
                                 "ARGs" = "green4",
                                 "Non-Plasmid/Vector/ARG" = "orange"
                                 )) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15),
      legend.title = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold")
    ) +
    labs(
      title = "Fraction of ARGs/Plasmids per sample",
      x = "Sample",
      y = "ng/ul",
      fill = "Hit Type"
    )
  
  print(p)
}

add_dna_data <- function(frac_df, values) {
  frac_df <- frac_df %>% 
    mutate(post_enr_c = values)
}


calculate_plasmid_coverage <- function(plasmid_hit_prevalence,
                                       non_plasmid_hit_prevalence,
                                       arg_hit_prevalence,
                                       frac_df,
                                       sample_data,
                                       hits_per_sample,
                                       fecal = FALSE,
                                       custom_order = NULL) {
  
  # Define sample order based on 'fecal' flag or use 'custom_order' if provided
  if (fecal) {
    default_order <- c(
      "6323_1x(1)", "6323_1x(2)", "10e-1ng_enriched", 
      "10e-2ng_enriched", "10e-3ng_enriched", "10e-4ng_enriched", 
      "10e-5ng_enriched", "10e-6ng_enriched"
    )
  } else {
    default_order <- c("25ng", "10ng", "1ng", "0.1ng", "0.05ng", "0.025ng", "0.01ng")
    # sample_data <- sample_data %>% filter(Sample_Name != "25ng_control") # if needed
  }
  
  # Use custom_order if provided, else use default_order
  order <- if (!is.null(custom_order)) custom_order else default_order
  
  # Replace NA values in frac_df with 0
  frac_df[is.na(frac_df)] <- 0
  
  # Calculate total unique plasmids
  total_unique_plasmids <- plasmid_hit_prevalence %>%
    distinct(V17) %>%
    nrow()
  
  # Calculate total unique non-plasmids
  total_unique_non_plasmids <- non_plasmid_hit_prevalence %>%
    distinct(V17) %>%
    nrow()
  
  # Calculate total unique ARGs
  total_unique_args <- arg_hit_prevalence %>%
    distinct(V17) %>%
    nrow()
  
  # Summarize unique plasmid hits per sample
  n_unique_plasmid_hits <- plasmid_hit_prevalence %>%
    group_by(Sample_Name) %>%
    summarize(n_unique_plasmids = n_distinct(V17), .groups = 'drop')
  
  # Summarize unique non-plasmid hits per sample
  n_unique_non_plasmid_hits <- non_plasmid_hit_prevalence %>%
    group_by(Sample_Name) %>%
    summarize(n_unique_non_plasmids = n_distinct(V17), .groups = 'drop')
  
  # Summarize unique ARG hits per sample
  n_unique_arg_hits <- arg_hit_prevalence %>%
    group_by(Sample_Name) %>%
    summarize(n_unique_args = n_distinct(V17), .groups = 'drop')
  
  # Merge all metrics into a unified data frame
  unique_df <- frac_df %>%
    left_join(hits_per_sample, by = "Sample_Name") %>%
    left_join(n_unique_plasmid_hits, by = "Sample_Name") %>%
    left_join(n_unique_non_plasmid_hits, by = "Sample_Name") %>%
    left_join(n_unique_arg_hits, by = "Sample_Name") %>% 
    # Ensure all samples in 'order' are included
    complete(
      Sample_Name = order, 
      fill = list(
        plasmid_hits = 0, 
        non_plasmid_hits = 0, 
        arg_hits = 0, 
        hits = 0,
        n_unique_args = 0, 
        n_unique_plasmids = 0, 
        n_unique_non_plasmids = 0
      )
    ) %>%
    # Rename for consistency
    rename(sample = Sample_Name) %>%
    # Convert 'sample' to a factor with levels defined by 'order'
    mutate(sample = factor(sample, levels = order)) %>%
    # Arrange based on the factor levels
    arrange(sample) %>%
    # Calculate fraction plasmid coverage and fraction ARG coverage
    mutate(
      frac_plasmid_coverage = n_unique_plasmids / total_unique_plasmids,
      frac_arg_coverage = n_unique_args / total_unique_args
    )
  
  # ----- Transform data to long format -----
  # In this example, we compare fraction plasmid coverage vs fraction ARG coverage side by side
  coverage_long <- unique_df %>%
    select(sample, frac_plasmid_coverage, frac_arg_coverage) %>%
    tidyr::pivot_longer(
      cols = c(frac_plasmid_coverage, frac_arg_coverage),
      names_to = "coverage_type",
      values_to = "coverage_value"
    )
  
  # ----- Grouped Bar Plot: Plasmid vs. ARG coverage -----
  p_coverage <- ggplot(coverage_long, aes(x = sample, y = coverage_value, fill = coverage_type)) +
    geom_col(position = position_dodge(width = 0.9)) +
    scale_fill_manual(
      values = c("frac_plasmid_coverage" = "steelblue", 
                 "frac_arg_coverage" = "green4"),
      labels = c("Plasmid Coverage", "ARG Coverage")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_blank()
    ) +
    labs(
      title = "Fraction of Unique Plasmids and ARGs Covered per Sample",
      x = "Sample",
      y = "Fraction Covered"
    )
  
  # Print the plot
  print(p_coverage)
  
  # Return the updated data frame
  return(unique_df)
}


calculate_arg_hit_prevalence <- function(data, sample_metadata) {
  result <- data %>%
    # Filter rows where V17 does NOT contain "plasmid" or "vector" (case-insensitive)
    filter(str_detect(V17, regex("res|resistance|antibiotic|arg", ignore_case = TRUE))) %>%
    
    # Group by sample_id, sample_type, and V17
    group_by(sample_id, sample_type, V17) %>%
    
    # Count the number of occurrences for each V17 within the groups
    summarise(n = n(), .groups = "drop") %>%
    
    # Rename sample_id to sample to match sample_metadata
    rename(sample = sample_id) %>%
    
    # Join with the sample metadata to get Sample_Name and other details
    left_join(sample_metadata, by = c("sample" = "ID", "sample_type" = "sample_type")) %>%
    
    # Select and arrange the desired columns
    select(V17, sample, Sample_Name, n)
  
  return(result)
}

#------------#
# FECAL DATA
#------------#

fecal_data <- all_samples %>% filter(sample_type == "fecal")
fecal_metadata <- all_sample_data %>% filter(sample_type == "fecal")
fecal_dna_values <- c(1.34, 1.36, 1.58, 1.22, 0.846, 0.542, 0.486, 0.418)

hits_per_fecal_sample <- calculate_hits_per_sample(fecal_data, fecal_metadata)
unique_hits_fecal <- calculate_unique_hits(fecal_data, fecal_metadata)
plasmid_hit_prevalence <- calculate_plasmid_hit_prevalence(fecal_data, fecal_metadata)

arg_hit_prevalence <- calculate_arg_hit_prevalence(fecal_data, fecal_metadata) 

non_plasmid_hit_prevalence <- calculate_non_plasmid_hit_prevalence(fecal_data, fecal_metadata)

frac_df_fecal <- calculate_fraction_hits(hits_per_fecal_sample, plasmid_hit_prevalence, non_plasmid_hit_prevalence, arg_hit_prevalence, fecal_metadata)
isoform_stats_fecal <- calculate_isoform_statistics(fecal_data)
frac_df_fecal <- add_dna_data(frac_df_fecal, fecal_dna_values)
plot_plasmid_fraction(frac_df_fecal, fecal = T)


unique_fecal <- calculate_plasmid_coverage(plasmid_hit_prevalence, non_plasmid_hit_prevalence, arg_hit_prevalence, frac_df_fecal, fecal_sample_data, hits_per_fecal_sample, fecal = T) 



#----------------#
# Spike-In Data
#----------------#
spike_in_data <- all_samples %>% filter(sample_type == "spike_in")
spike_in_metadata <- all_sample_data %>% filter(sample_type == "spike_in")
spike_in_dna_values <- c(2.48, 5.22, 1.44, 1.2, 0.8, 0.61, 0.376)

hits_per_spike_sample <- calculate_hits_per_sample(spike_in_data, spike_in_metadata)
unique_hits_spike <- calculate_unique_hits(spike_in_data, spike_in_metadata)
plasmid_hit_prevalence_spike <- calculate_plasmid_hit_prevalence(spike_in_data, spike_in_metadata)

arg_hit_prevalence <- calculate_arg_hit_prevalence(spike_in_data, spike_in_metadata) 

non_plasmid_hit_prevalence_spike <- calculate_non_plasmid_hit_prevalence(spike_in_data, spike_in_metadata)
frac_df_spike_in <- calculate_fraction_hits(hits_per_spike_sample, plasmid_hit_prevalence_spike, non_plasmid_hit_prevalence_spike, arg_hit_prevalence, spike_in_metadata)
isoform_stats_spike <- calculate_isoform_statistics(spike_in_data)
frac_df_spike_in <- add_dna_data(frac_df_spike_in, spike_in_dna_values)
plot_plasmid_fraction(frac_df_spike_in)


unique_spike <- calculate_plasmid_coverage(plasmid_hit_prevalence_spike, non_plasmid_hit_prevalence_spike, arg_hit_prevalence, frac_df_spike_in, spike_in_sample_data, hits_per_spike_sample, fecal = F) 




