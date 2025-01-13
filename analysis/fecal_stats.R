library(tidyverse)
library(ComplexUpset)
library(RColorBrewer)
# Set WD
downloads <- "C:/Users/abergm/Downloads/fecal/"

# File set 1
temp <- c("25ng_unenriched_S1_kraken2_fastp_0.1.report",
             "25ng_VERO_E6_enriched_S8_kraken2_fastp_0.1.report",
             "10ng_enriched_S2_kraken2_fastp_0.1.report",
             "1ng_enriched_S3_kraken2_fastp_0.1.report",
             "10e-1ng_enriched_S4_kraken2_fastp_0.1.report",
             "5e-2ng_enriched_S6_kraken2_fastp_0.1.report",
             "2.5e-2ng_enriched_S7_kraken2_fastp_0.1.report",
             "10e-2ng_enriched_S5_kraken2_fastp_0.1.report")

# File set 2 
samples <- c("6323_1x_S89_kraken2_fastp_0.1.report",
             "6323_1x_S96_kraken2_fastp_0.1.report",
             "6323_10e-1_S90_kraken2_fastp_0.1.report",
             "6323_10e-2_S91_kraken2_fastp_0.1.report",
             "6323_10e-3_S92_kraken2_fastp_0.1.report",
             "6323_10e-4_S93_kraken2_fastp_0.1.report",
             "6323_10e-5_S94_kraken2_fastp_0.1.report",
             "6323_10e-6_S95_kraken2_fastp_0.1.report")

# Sample names
sample_names <- c("1x",
                  "1x(2)",
                  "10e-1",
                  "10e-2",
                  "10e-3",
                  "10e-4",
                  "10e-5",
                  "10e-6")

# Names file set 2
temp <- c("25ng_control",
                  "25ng",
                  "10ng",
                  "1ng",
                  "0.1ng",
                  "0.05ng",
                  "0.025ng",
                  "0.01ng")

# Read sample files
read_and_process <- function(file, sample_name) {
  df <- read.table(paste0(downloads, file),
                   sep = "\t",
                   header = FALSE,
                   quote = "",
                   stringsAsFactors = FALSE)
  colnames(df) <- c("percent_reads", "reads_clade", "reads_direct",
                    "rank_code", "tax_id", "name")
  
  df <- df %>%
    mutate(sample = sample_name)
  
  return(df)
}


# Plots all levels of classification
plot <- function(all_data, plot_dir) {
  
  ranks_to_plot <- c("D","K","P","C","O","F","G","S")
  
  # Create a directory for plots
  dir.create(plot_dir, showWarnings = FALSE)
  
  # order of samples
  desired_order <- c("1x", "1x(2)", "10e-1", "10e-2", 
                     "10e-3", "10e-4", "10e-5", "10e-6")
  
  # Color palette for plotting phyla
  phylum_num_colors <- 11  # Top 10 phyla + Others
  phylum_palette_name <- "Paired"
  
  # Get the colors from "Paired"
  phylum_colors <- brewer.pal(n = phylum_num_colors, name = phylum_palette_name)
  
  # Assign the first 10 colors to the top 10 phyla and the last color to "Others"
  others_color <- phylum_colors[length(phylum_colors)]
  phylum_colors <- phylum_colors[1:(phylum_num_colors - 1)]
  
  # Loop through each rank
  for (rank in ranks_to_plot) {
    # Filter data for the current rank
    rank_df <- all_data %>%
      filter(rank_code == rank)
    
    if (nrow(rank_df) == 0) {
      # If no data for this rank, skip it
      next
    }
    
    # *** Begin Modification: Grouping Top 10 Phyla and Others ***
    if (rank == "P") {
      # Determine the top 10 phyla
      top10_phyla <- rank_df %>%
        group_by(name) %>%
        summarize(total_reads = sum(reads_clade), .groups = 'drop') %>%
        arrange(desc(total_reads)) %>%
        slice_head(n = 10) %>%
        pull(name)
      
      # Replace phyla not in top 10 with "Others"
      rank_df <- rank_df %>%
        mutate(name = ifelse(name %in% top10_phyla, name, "Others"))
      
      # Reorder 'name' factor to have top 10 phyla first, then "Others"
      rank_df$name <- factor(rank_df$name, levels = c(top10_phyla, "Others"))
      
      # Assign color
      fill_scale <- scale_fill_manual(values = c(phylum_colors, others_color))
      
    } else {
      num_categories <- length(unique(rank_df$name))
      
      # Choose a palette based on the number of categories
      other_palette_name <- "Set3"
      max_palette_colors <- brewer.pal.info[other_palette_name, "maxcolors"]
      
      if (num_categories <= max_palette_colors) {
        other_colors <- brewer.pal(n = num_categories, name = other_palette_name)
      } else {
        # Generate additional colors if needed
        other_colors <- colorRampPalette(brewer.pal(max_palette_colors, other_palette_name))(num_categories)
      }
      
      fill_scale <- scale_fill_manual(values = other_colors)
    }
    
    # Normalize by the sum of reads_clade at this rank
    rank_totals <- rank_df %>%
      group_by(sample) %>%
      summarize(total_rank_reads = sum(reads_clade), .groups = 'drop')
    
    rank_df <- rank_df %>%
      left_join(rank_totals, by = "sample") %>%
      mutate(fraction = reads_clade / total_rank_reads)
    
    # Convert 'sample' to a factor with the specified order
    rank_df$sample <- factor(rank_df$sample, levels = desired_order)
    
    p <- ggplot(rank_df, aes(x = sample, y = fraction, fill = name)) +
      geom_col(position = "stack") +
      theme_minimal(base_size = 14) +  
      fill_scale +  
      guides(fill = guide_legend(ncol = 2, title = ifelse(rank == "P", "Phylum", rank))) + 
      theme(
        legend.text = element_text(size = 10),           
        legend.title = element_text(size = 12, face = "bold"),  
        legend.key.size = unit(0.8, "lines"),            
        legend.spacing.x = unit(0.4, "cm"),               
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      ) +
      labs(
        title = ifelse(rank == "P", "Distribution of Top 10 Phyla and Others", paste("Distribution of", rank, "Ranks")),
        x = "Sample",
        y = "Fraction of Reads"
      )
    
    # Display the plot
    print(p)
    
    # Save the plot
    outfile <- file.path(plot_dir, paste0("plot_", rank, ".png"))
    ggsave(outfile, p, width = 12, height = 8, dpi = 300)
  }
}

# Convert sample names into numeric dilution factors.
parse_dilution <- function(x) {
  # Remove parentheses if present
  x <- gsub("\\(.*\\)", "", x)
  
  # If it contains '10e-', we can interpret that as scientific notation
  # "10e-1" means 10^-1 = 0.1
  if (grepl("10e-", x)) {
    # Extract the exponent
    exp_part <- sub("10e-", "", x)
    as.numeric(paste0("1e-", exp_part))
  } else if (x == "1x") {
    1
  } else {
    # Default to 1 if something unexpected
    1
  }
}

# Plot the microbial abundance fold change in phyla
plot_enrichment <- function(all_data, plot_dir) {
  
  all_data <- all_data %>%
    mutate(name = trimws(name)) %>%
    mutate(dilution = sapply(sample, parse_dilution))
  
  # Get phylum-level data
  phylum_data <- all_data %>%
    filter(rank_code == "P")
  
  # Compute total phylum reads per sample and fractions
  phylum_totals <- phylum_data %>%
    group_by(sample) %>%
    summarize(total_phylum_reads = sum(reads_clade))
  
  phylum_data <- phylum_data %>%
    left_join(phylum_totals, by = "sample") %>%
    mutate(fraction = reads_clade / total_phylum_reads)
  
  # Identify unique phyla
  phyla <- unique(phylum_data$name)
  
  # Extract fractions for the baseline samples
  baseline_sample <- 1  # Adjust based on your dataset logic
  one_x_data <- phylum_data %>% filter(dilution == baseline_sample)
  
  # Compute mean baseline fractions, handling missing phyla with a default value
  one_x_fractions <- one_x_data %>%
    group_by(name) %>%
    summarize(mean_1x_fraction = mean(fraction, na.rm = TRUE)) %>%
    ungroup()
  
  default_baseline_value <- 1e-3  # Define a small default value for missing taxa
  
  # Join these baseline fractions back to the full dataset and compute fold change
  phylum_data <- phylum_data %>%
    left_join(one_x_fractions, by = "name") %>%
    mutate(
      mean_1x_fraction = ifelse(is.na(mean_1x_fraction), default_baseline_value, mean_1x_fraction),
      fold_change = ifelse(dilution != baseline_sample, fraction / mean_1x_fraction, NA)
    )
  
  # Create a table showing fold changes of each phylum at each dilution
  fold_change_table <- phylum_data %>%
    filter(dilution != baseline_sample) %>%
    select(name, sample, dilution, fraction, mean_1x_fraction, fold_change) %>%
    arrange(name, dilution) %>% 
    mutate(dilution = factor(dilution))
  
  print(fold_change_table)
  
  # Arrange table 
  fold_change_table <- fold_change_table %>%
    arrange(desc(fold_change))
  
  dilutions <- unique(fold_change_table$dilution)
  plots <- list()

  # loop through dilutions and plot the diffs compared to the control / 1x
  for (this_dilution in dilutions) {
    # Filter and arrange data for the current dilution
    sample_data <- fold_change_table %>% 
      filter(dilution == this_dilution)
    
    # Create the plot
    plot_title <- paste("Fold Change", this_dilution, "vs Baseline")
    
    p <- ggplot(sample_data, aes(x = reorder(name, fold_change), y = fold_change)) +
      geom_col(aes(fill = mean_1x_fraction == default_baseline_value), show.legend = TRUE) +
      coord_flip() +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      scale_fill_manual(
        values = c("TRUE" = "orange", "FALSE" = "gray"),
        labels = c("TRUE" = "Absent in Baseline", "FALSE" = "Present in Baseline"),
        name = "Baseline Presence"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1.5, size = 12),
        axis.text.y = element_text(vjust = 1, size = 7)
      ) +
      ggtitle(plot_title)
    
    # Append the plot to the list
    plots[[length(plots) + 1]] <- p
    
    ggsave(filename = paste0(plot_dir, "plot_dilution_", this_dilution, ".png"), plot = p, width = 8, height = 6)
  }
  
  return(fold_change_table)
}

# Create a table visualising common enrichments throuhgout the samples
find_common_enrichments <- function(fold_change_table, plot_dir) {

  # Calculate average enrichment per sample
  average_enrichment <- fold_change_table %>%
    group_by(sample) %>%
    summarize(mean_enrichment = mean(fold_change, na.rm = TRUE)) %>%
    arrange(desc(mean_enrichment))  # Sort by highest enrichment
  
  # Print the average enrichment table
  print(average_enrichment)
  
  # Filter for enriched microorganisms (fold change > 1)
  enriched_data <- fold_change_table %>%
    filter(fold_change > 1) %>%
    select(name, sample)  # Retain microorganism names and samples
  
  print(enriched_data)
  
  # Calculate overlaps between microorganisms
  overlap_table <- enriched_data %>%
    group_by(name) %>%
    summarize(num_samples = n(), samples = paste(unique(sample), collapse = ", ")) %>%
    arrange(desc(num_samples))  # Sort by the number of samples
  
  print(overlap_table)
  
  # Bar Plot of Overlaps
  ggplot(overlap_table, aes(x = reorder(name, -num_samples), y = num_samples)) +
    geom_col() +
    labs(
      title = "Overlapping Enriched Microorganisms Across Samples",
      x = "Microorganism",
      y = "Number of Samples"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12))
  
  ggsave(filename= paste0(plot_dir, "overlap_plot.png"), width = 10, height = 6)

  return(overlap_table)
}

# Plot the abundance of select phyla throughout the samples
assess_diffs_between_common_phyla <- function(fold_change_table) {
  interesting_phyla <- fold_change_table %>% 
    filter(name %in% c("Chordata")) %>% 
    mutate(name = as.factor(name))
  
  # Create the plot
  p <- ggplot(interesting_phyla, aes(x = sample, y = fold_change, fill = name)) +
    geom_col(position = "stack") +
    facet_wrap(~name, scales = "free") +
    scale_fill_brewer(palette = "Dark2") +  # Use a ColorBrewer palette for better aesthetics
    labs(title = "Fold Change Between Samples for Selected Phyla",
         x = "Sample",
         y = "Fold Change",
         fill = "Phylum") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
  
  print(p)
  ggsave(filename= paste0(plot_dir, "interesting_phyla.png"), width = 10, height = 6)
  
  return(interesting_phyla)
}

# Plot directory 
plot_dir <- "C:/Users/abergm/Downloads/hierarchy_plots/fecal/"
all_data <- map2_dfr(samples, sample_names, read_and_process)

# Plot all classification levels
plot(all_data, plot_dir)
# Calculate fold change
fold_change_table <- plot_enrichment(all_data, plot_dir)
fold_change_table <- fold_change_table %>% 
  arrange(desc(fold_change))

# Find common phyla throughout the dilution series
common_enrichment <- find_common_enrichments(fold_change_table, plot_dir)

# Find "interesting" phyla (not used)
interesting_phyla <- assess_diffs_between_common_phyla(fold_change_table)

