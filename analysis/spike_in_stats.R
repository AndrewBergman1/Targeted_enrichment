library(tidyverse)
library(RColorBrewer)

# WD
downloads <- "C:/Users/abergm/Downloads/"

# Sample names
samples <- c("25ng_unenriched_S1_kraken2_fastp_0.1.report",
          "25ng_VERO_E6_enriched_S8_kraken2_fastp_0.1.report",
          "10ng_enriched_S2_kraken2_fastp_0.1.report",
          "1ng_enriched_S3_kraken2_fastp_0.1.report",
          "10e-1ng_enriched_S4_kraken2_fastp_0.1.report",
          "5e-2ng_enriched_S6_kraken2_fastp_0.1.report",
          "2.5e-2ng_enriched_S7_kraken2_fastp_0.1.report",
          "10e-2ng_enriched_S5_kraken2_fastp_0.1.report")

# sample names
sample_names <- c("25ng",
          "25ng_control",
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
# Define the plotting function
plot <- function(all_data, plot_dir) {
  
  # Define the ranks for plotting
  ranks_to_plot <- c("D","K","P","C","O","F","G","S")
  
  # Create a directory for plots
  dir.create(plot_dir, showWarnings = FALSE)
  
  # Sample order
  desired_order <- c("25ng_control", "25ng", "10ng", "1ng", 
                     "0.1ng", "0.05ng", "0.025ng", "0.01ng")
  

  #define colors
  phylum_num_colors <- 11  # Top 10 phyla + Others
  phylum_palette_name <- "Paired"
  
  # Get the colors from the palette for phylum rank
  phylum_colors <- brewer.pal(n = phylum_num_colors, name = phylum_palette_name)

  # plot all levels of classification
  for (rank in ranks_to_plot) {
    # Filter data for the current rank
    rank_df <- all_data %>%
      filter(rank_code == rank)
    
    if (nrow(rank_df) == 0) {
      # If no data for this rank, skip it
      next
    }
    
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
      
      rank_df$name <- factor(rank_df$name, levels = c(top10_phyla, "Others"))
      
      # Assign colors for phyla
      fill_scale <- scale_fill_manual(values = phylum_colors)
      
    } else {
      num_categories <- length(unique(rank_df$name))
      
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
    }
    
    # Calculate total reads for this rank per sample
    rank_totals <- rank_df %>%
      group_by(sample) %>%
      summarize(total_rank_reads = sum(reads_clade), .groups = 'drop')
    
    rank_df <- rank_df %>%
      left_join(rank_totals, by = "sample") %>%
      mutate(fraction = reads_clade / total_rank_reads)
    
    # Convert 'sample' to a factor with the specified order
    rank_df$sample <- factor(rank_df$sample, levels = desired_order)
    
    # Plot classifications
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
    print(p)
    
    # Save the plot
    outfile <- file.path(plot_dir, paste0("plot_", rank, ".png"))
    ggsave(outfile, p, width = 12, height = 8, dpi = 300)
  }
}
# Convert sample names into numeric dilution factors.
parse_dilution <- function(x) {
  # Remove "_control" suffix if present
  x <- gsub("_control", "", x)
  
  as.numeric(gsub("ng", "", x))
}

# plot the fold change of microbes throughotu the samples
plot_enrichment <- function(all_data, plot_dir) {

  # remove whitespaces
  all_data <- all_data %>%
    mutate(name = trimws(name)) %>%
    mutate(dilution = parse_dilution(sample))  # Use the updated parser
  
  # Focus on phylum-level 
  phylum_data <- all_data %>%
    filter(rank_code == "P")
  
  # Compute total phylum reads per sampel
  phylum_totals <- phylum_data %>%
    group_by(sample) %>%
    summarize(total_phylum_reads = sum(reads_clade))

  # get fracs of each phyla
  phylum_data <- phylum_data %>%
    left_join(phylum_totals, by = "sample") %>%
    mutate(fraction = reads_clade / total_phylum_reads)
  
  # Identify baseline
  baseline_sample <- "25ng_control"  
  one_x_data <- phylum_data %>% filter(sample == baseline_sample)
  
  # Calculate mean fractions for the baseline sample
  one_x_fractions <- one_x_data %>%
    group_by(name) %>%
    summarize(mean_1x_fraction = mean(fraction, na.rm = TRUE)) %>%
    ungroup()
  
  default_baseline_value <- 1e-3 # to calculate fold change of phyla absent in ctrl

  # calculate fodl change
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
  
  print(fold_change_table)
  
  fold_change_table <- fold_change_table %>%
    arrange(desc(fold_change))
  
  dilutions <- unique(fold_change_table$dilution)
  plots <- list()

  # plot fold change for all samples
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
    
    plots[[length(plots) + 1]] <- p

    # save file
    ggsave(filename = paste0(plot_dir, "plot_dilution_", this_dilution, ".png"), plot = p, width = 8, height = 6)
  }
  
  return(fold_change_table)
}

# find common enrihcments throughout the sampels
find_common_enrichments <- function(fold_change_table, plot_dir) {

  # Get average enrichment per sample
  average_enrichment <- fold_change_table %>%
    group_by(sample) %>%
    summarize(mean_enrichment = mean(fold_change, na.rm = TRUE)) %>%
    arrange(desc(mean_enrichment))  # Sort by highest enrichment
  
  print(average_enrichment)

  # get fold change > 1 
  enriched_data <- fold_change_table %>%
    filter(fold_change > 1) %>%
    select(name, sample)  # Retain microorganism names and samples
  
  print(enriched_data)

  # Fnid overlapping microbes
  overlap_table <- enriched_data %>%
    group_by(name) %>%
    summarize(num_samples = n(), samples = paste(unique(sample), collapse = ", ")) %>%
    arrange(desc(num_samples))  # Sort by the number of samples
  
  print(overlap_table)
  
  # plot overlaps
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
# See how common phyla differs among the samples
assess_diffs_between_common_phyla <- function(fold_change_table) {
  interesting_phyla <- fold_change_table %>% 
    filter(name %in% c("Chordata")) %>% 
    mutate(name = as.factor(name))

  # plot "interesting" phyla
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

# Define the output directory
plot_dir <- "C:/Users/abergm/Downloads/hierarchy_plots/"

# Process the data
all_data <- map2_dfr(samples, sample_names, read_and_process)

# plot data
plot(all_data, plot_dir)

# Generate fold change plots
fold_change_table <- plot_enrichment(all_data, plot_dir)

common_enrichment <- find_common_enrichments(fold_change_table, plot_dir)
interesting_phyla <- assess_diffs_between_common_phyla(fold_change_table)

sorted <- fold_change_table %>% 
  arrange(desc(fold_change))

