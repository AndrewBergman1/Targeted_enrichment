library(Biostrings)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

original_probes_file <- args[1]
filtered_probes_file <- args[2]
semifiltered_probes_file <- args[3]  # Ensure you have this argument
output_plot <- args[4]

# Read the DNA sequences
original_probes <- readDNAStringSet(original_probes_file)
filtered_probes <- readDNAStringSet(filtered_probes_file)
semifiltered_probes <- readDNAStringSet(semifiltered_probes_file)  # Read the semifiltered probes

# Calculate nucleotide frequencies
original_nt <- as.data.frame(alphabetFrequency(original_probes))
filtered_nt <- as.data.frame(alphabetFrequency(filtered_probes))
semifiltered_nt <- as.data.frame(alphabetFrequency(semifiltered_probes))  # Calculate frequencies for semifiltered

# Create data frames for original, filtered, and semi-filtered probes
original_df <- data.frame(source = "original_probes", original_nt)
filtered_df <- data.frame(source = "filtered_probes", filtered_nt)
semifiltered_df <- data.frame(source = "semifiltered_probes", semifiltered_nt)

# Calculate Tm and GC for original probes
original_df <- original_df %>%
  mutate(Tm = 64.9 + 41 * (G + C - 16.4) / (A + T + G + C)) %>%
  mutate(GC = (G + C) / (A + T + G + C)) %>%
  select(c("A", "T", "C", "G", "GC", "Tm")) %>%
  mutate(source = "Original")

# Calculate Tm and GC for filtered probes
filtered_df <- filtered_df %>%
  mutate(Tm = 64.9 + 41 * (G + C - 16.4) / (A + T + G + C)) %>%
  mutate(GC = (G + C) / (A + T + G + C)) %>%
  select(c("A", "T", "C", "G", "GC", "Tm")) %>%
  mutate(source = "Filtered")

# Calculate Tm and GC for semi-filtered probes
semifiltered_df <- semifiltered_df %>%
  mutate(Tm = 64.9 + 41 * (G + C - 16.4) / (A + T + G + C)) %>%
  mutate(GC = (G + C) / (A + T + G + C)) %>%  # Ensure consistent calculation
  select(c("A", "T", "C", "G", "GC", "Tm")) %>%
  mutate(source = "Semi-filtered")

# Combine data frames
df <- rbind(original_df, filtered_df, semifiltered_df)

# Transform data frame to long format
df_long <- pivot_longer(df, cols = c("GC", "Tm"), values_to = "Value", names_to = "Metric")

# Create boxplot
q <- ggplot(df_long, aes(x = source, y = Value)) +
  geom_boxplot() +
  facet_wrap(~Metric, scales = "free")

# Save the plot
ggsave(output_plot, plot = q)

# Print summary statistics
cat("Number of Original probes:", length(original_probes), "\n",
    "GC% IQR of original probes:", IQR(original_df$GC), "\n",
    "Tm IQR of original probes:", IQR(original_df$Tm), "\n",
    "\n",
    "Number of filtered probes:", length(filtered_probes), "\n",
    "GC% IQR of filtered probes:", IQR(filtered_df$GC), "\n",
    "Tm IQR of filtered probes:", IQR(filtered_df$Tm), "\n",
    "\n",
    "Number of semi-filtered probes:", length(semifiltered_probes), "\n",
    "GC% IQR of semi-filtered probes:", IQR(semifiltered_df$GC), "\n",
    "Tm IQR of semi-filtered probes:", IQR(semifiltered_df$Tm), "\n"
)
