# Load necessary libraries
library(dplyr)
library(stringr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
output <- args[2]
plot <- args[3]

# Read the file
data <- readLines(input)

# Create a vector of indices for IDs and sequences
id_indices <- seq(1, length(data), by = 3)
sequence_indices <- id_indices + 1
structure_indices <- id_indices + 2

# Extract IDs, sequences, structures, and free energy
ids <- str_extract(data[id_indices], "\\d+")
sequences <- data[sequence_indices]
structures <- data[structure_indices]
free_energy <- as.numeric(str_extract(data[structure_indices], "-?\\d+\\.\\d+"))

# Create a data frame
rna_data <- data.frame(
  ID = ids,
  Sequence = sequences,
  Structure = structures,
  Free_Energy = free_energy,
  stringsAsFactors = FALSE
)

# Display the data frame
print(rna_data)

# Plotting free energy distribution
q <- ggplot(rna_data, aes(x = Free_Energy)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = -45, color = "red") +
  labs(title = "Free Energy Distribution", x = "Free Energy", y = "Count")

# Filter data based on free energy
filtered_data <- rna_data %>%
  filter(Free_Energy > -45)

# Write filtered data to CSV
write.csv(filtered_data, output, row.names = FALSE)

ggsave(plot, q)
