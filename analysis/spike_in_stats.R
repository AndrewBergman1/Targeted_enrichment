library(tidyverse)
library(RColorBrewer)

downloads <- "C:/Users/abergm/Downloads/"

samples <- c("25ng_unenriched_S1_kraken2_fastp_0.1.report",
          "25ng_VERO_E6_enriched_S8_kraken2_fastp_0.1.report",
          "10ng_enriched_S2_kraken2_fastp_0.1.report",
          "1ng_enriched_S3_kraken2_fastp_0.1.report",
          "10e-1ng_enriched_S4_kraken2_fastp_0.1.report",
          "5e-2ng_enriched_S6_kraken2_fastp_0.1.report",
          "2.5e-2ng_enriched_S7_kraken2_fastp_0.1.report",
          "10e-2ng_enriched_S5_kraken2_fastp_0.1.report")

fecal <- c("6323_1x_S89_kraken2_fastp_0.1.report",
             "6323_1x_S96_kraken2_fastp_0.1.report",
             "6323_10e-1_S90_kraken2_fastp_0.1.report",
             "6323_10e-2_S91_kraken2_fastp_0.1.report",
             "6323_10e-3_S92_kraken2_fastp_0.1.report",
             "6323_10e-4_S93_kraken2_fastp_0.1.report",
             "6323_10e-5_S94_kraken2_fastp_0.1.report",
             "6323_10e-6_S95_kraken2_fastp_0.1.report")

temp <- c("1x",
                  "1x(2)",
                  "10e-1",
                  "10e-2",
                  "10e-3",
                  "10e-4",
                  "10e-5",
                  "10e-6")

sample_names <- c("25ng",
          "25ng_control",
          "10ng",
          "1ng",
          "0.1ng",
          "0.05ng",
          "0.025ng",
          "0.01ng")

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
}# Define the plotting function
plot <- function(all_data, plot_dir) {
  
  # Define the ranks you want to plot:
  # D = domain, K = kingdom, P = phylum, C = class, O = order, F = family, G = genus, S = species
  ranks_to_plot <- c("D","K","P","C","O","F","G","S")
  
  # Create a directory for plots if it doesn't exist
  dir.create(plot_dir, showWarnings = FALSE)
  
  # Define the desired order of samples
  desired_order <- c("25ng_control", "25ng", "10ng", "1ng", 
                     "0.1ng", "0.05ng", "0.025ng", "0.01ng")
  
  # Define a color palette for publication-ready colors
  # Option 1: Using RColorBrewer's "Paired" palette (12 colors)
  # Option 2: Using RColorBrewer's "Set3" palette (12 colors)
  # Option 3: Using viridis package (requires installation)
  
  # Here, we'll use "Paired" for its distinct colors for phylum rank
  phylum_num_colors <- 11  # Top 10 phyla + Others
  phylum_palette_name <- "Paired"
  
  # Check if the palette has enough colors
  if (phylum_num_colors > brewer.pal.info[phylum_palette_name, "maxcolors"]) {
    stop(paste("Palette", phylum_palette_name, "does not have enough colors."))
  }
  
  # Get the colors from the palette for phylum rank
  phylum_colors <- brewer.pal(n = phylum_num_colors, name = phylum_palette_name)
  
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
      # Determine the top 10 phyla based on total reads across all samples
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
      
      # Assign colors for phyla
      fill_scale <- scale_fill_manual(values = phylum_colors)
      
    } else {
      # For other ranks, use a different color palette
      # Option 1: Automatically assign colors
      # Option 2: Use a predefined palette with enough colors
      
      # Determine the number of unique categories in 'name'
      num_categories <- length(unique(rank_df$name))
      
      # Choose a palette based on the number of categories
      # Here, we'll use "Set3" palette which has up to 12 colors
      other_palette_name <- "Set3"
      
      # If number of categories exceeds the palette, use a larger palette or repeat colors
      if (num_categories <= brewer.pal.info[other_palette_name, "maxcolors"]) {
        other_colors <- brewer.pal(n = num_categories, name = other_palette_name)
      } else {
        # Use a different palette or generate colors programmatically
        other_colors <- colorRampPalette(brewer.pal(12, other_palette_name))(num_categories)
      }
      
      # Assign colors
      fill_scale <- scale_fill_manual(values = other_colors)
      
      # Alternatively, you can use a continuous palette or another discrete palette
      # Example using viridis:
      # library(viridis)
      # fill_scale <- scale_fill_viridis_d(option = "D", begin = 0.1, end = 0.9)
    }
    # *** End Modification ***
    
    # Compute total reads for this rank per sample
    # Normalize by the sum of reads_clade at this rank
    rank_totals <- rank_df %>%
      group_by(sample) %>%
      summarize(total_rank_reads = sum(reads_clade), .groups = 'drop')
    
    rank_df <- rank_df %>%
      left_join(rank_totals, by = "sample") %>%
      mutate(fraction = reads_clade / total_rank_reads)
    
    # Convert 'sample' to a factor with the specified order
    rank_df$sample <- factor(rank_df$sample, levels = desired_order)
    
    # Create the ggplot with publication-ready colors
    p <- ggplot(rank_df, aes(x = sample, y = fraction, fill = name)) +
      geom_col(position = "stack") +
      theme_minimal(base_size = 14) +  # Increase base font size for readability
      fill_scale +  # Apply the defined color palette based on rank
      guides(fill = guide_legend(ncol = 2, title = ifelse(rank == "P", "Phylum", rank))) +  # Customize legend
      theme(
        legend.text = element_text(size = 10),            # Adjust legend text size
        legend.title = element_text(size = 12, face = "bold"),  # Legend title styling
        legend.key.size = unit(0.8, "lines"),             # Increase legend key size
        legend.spacing.x = unit(0.4, "cm"),               # Adjust horizontal spacing in legend
        axis.text.x = element_text(angle = 45, hjust = 1),# Tilt x-axis labels for readability
        axis.title = element_text(size = 14, face = "bold"),# Axis titles styling
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5) # Plot title styling
      ) +
      labs(
        title = ifelse(rank == "P", "Distribution of Top 10 Phyla and Others", paste("Distribution of", rank, "Ranks")),
        x = "Sample",
        y = "Fraction of Reads"
      )
    
    # Display the plot
    print(p)
    
    # Save the plot as a high-resolution PNG
    outfile <- file.path(plot_dir, paste0("plot_", rank, ".png"))
    ggsave(outfile, p, width = 12, height = 8, dpi = 300)
  }
}
# Convert sample names into numeric dilution factors.
# For example: "1x" -> 1; "10e-1" -> 0.1, "10e-2" -> 0.01, etc.
# We can write a small function to parse these:
parse_dilution <- function(x) {
  # Remove "_control" suffix if present
  x <- gsub("_control", "", x)
  
  # Convert units like "25ng" to numeric values (25, 10, 1, etc.)
  as.numeric(gsub("ng", "", x))
}


plot_enrichment <- function(all_data, plot_dir) {
  
  all_data <- all_data %>%
    mutate(name = trimws(name)) %>%
    mutate(dilution = parse_dilution(sample))  # Use the updated parser
  
  # Focus on phylum-level data
  phylum_data <- all_data %>%
    filter(rank_code == "P")
  
  # Compute total phylum reads per sample and fractions
  phylum_totals <- phylum_data %>%
    group_by(sample) %>%
    summarize(total_phylum_reads = sum(reads_clade))
  
  phylum_data <- phylum_data %>%
    left_join(phylum_totals, by = "sample") %>%
    mutate(fraction = reads_clade / total_phylum_reads)
  
  # Identify baseline (control) sample fractions
  baseline_sample <- "25ng_control"  # Adjust this if your control sample changes
  one_x_data <- phylum_data %>% filter(sample == baseline_sample)
  
  # Calculate mean fractions for the baseline sample
  one_x_fractions <- one_x_data %>%
    group_by(name) %>%
    summarize(mean_1x_fraction = mean(fraction, na.rm = TRUE)) %>%
    ungroup()
  
  default_baseline_value <- 1e-3
  
  phylum_data <- phylum_data %>%
    left_join(one_x_fractions, by = "name") %>%
    mutate(
      # Replace NA in mean_1x_fraction with a small value
      mean_1x_fraction = ifelse(is.na(mean_1x_fraction), default_baseline_value, mean_1x_fraction),
      # Calculate fold change
      fold_change = ifelse(
        sample != baseline_sample,
        fraction / mean_1x_fraction,
        NA
      )
    )
  
  # Create fold change table
  fold_change_table <- phylum_data %>%
    filter(sample != baseline_sample) %>%
    select(name, sample, dilution, fraction, mean_1x_fraction, fold_change) %>%
    arrange(dilution, desc(fold_change))
  
  # Print the fold change table
  print(fold_change_table)
  
  fold_change_table <- fold_change_table %>%
    arrange(desc(fold_change))
  
  dilutions <- unique(fold_change_table$dilution)
  plots <- list()
  
  for (this_dilution in dilutions) {
    # Filter and arrange data for the current dilution
    sample_data <- fold_change_table %>% 
      filter(dilution == this_dilution)
    
    # Create the plot with color-coded bars
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
        axis.text.y = element_text(vjust = 1, size = 12)
      ) +
      ggtitle(plot_title)
    
    # Append the plot to the list
    plots[[length(plots) + 1]] <- p
    
    ggsave(filename = paste0(plot_dir, "plot_dilution_", this_dilution, ".png"), plot = p, width = 8, height = 6)
  }
  
  return(fold_change_table)
}

find_common_enrichments <- function(fold_change_table, plot_dir) {
  
  # Step 1: Calculate average enrichment per sample
  average_enrichment <- fold_change_table %>%
    group_by(sample) %>%
    summarize(mean_enrichment = mean(fold_change, na.rm = TRUE)) %>%
    arrange(desc(mean_enrichment))  # Sort by highest enrichment
  
  # Print the average enrichment table
  print(average_enrichment)
  
  # Step 2: Filter for enriched microorganisms (fold change > 1)
  enriched_data <- fold_change_table %>%
    filter(fold_change > 1) %>%
    select(name, sample)  # Retain microorganism names and samples
  
  # Print enriched data
  print(enriched_data)
  
  # Step 3: Calculate overlaps between microorganisms
  overlap_table <- enriched_data %>%
    group_by(name) %>%
    summarize(num_samples = n(), samples = paste(unique(sample), collapse = ", ")) %>%
    arrange(desc(num_samples))  # Sort by the number of samples
  
  # Print overlap table
  print(overlap_table)
  
  # Step 4: Optional Visualization (Bar Plot of Overlaps)
  ggplot(overlap_table, aes(x = reorder(name, -num_samples), y = num_samples)) +
    geom_col() +
    labs(
      title = "Overlapping Enriched Microorganisms Across Samples",
      x = "Microorganism",
      y = "Number of Samples"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12))
  
  ggsave(filename= paste0(plot_dir, "microbial_overlap.png"), width = 10, height = 6)
  
  return(overlap_table)
}
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
  
  # Display the plot
  print(p)
  ggsave(filename= paste0(plot_dir, "interesting_phyla.png"), width = 10, height = 6)
  
  return(interesting_phyla)
}

# Define the output directory
plot_dir <- "C:/Users/abergm/Downloads/hierarchy_plots/"

# Process the data and create the plots
all_data <- map2_dfr(samples, sample_names, read_and_process)

plot(all_data, plot_dir)

# Generate fold change plots
fold_change_table <- plot_enrichment(all_data, plot_dir)

common_enrichment <- find_common_enrichments(fold_change_table, plot_dir)
interesting_phyla <- assess_diffs_between_common_phyla(fold_change_table)

sorted <- fold_change_table %>% 
  arrange(desc(fold_change))

