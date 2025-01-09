library(UpSetR)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
table <- args[1]

start_time <- Sys.time()
print(paste("Start time:", start_time))

# List all .txt files in the directory
datasets <- list.files("/home/abergm/bait_gen/outputs/kmer_freq/datasets/", pattern = "\\.txt$", full.names = TRUE)
syotti <- "/home/abergm/bait_gen/outputs/kmer_freq/syotti/syotti_kmers.txt"


# Combine syotti with datasets
all_files <- c(datasets, syotti)

# Initialize an empty list to store the datasets
db_list <- list()

# Loop through each file, read it, and add the file name as a new column
for (file in all_files) {
  name <- gsub("\\.txt$", "", basename(file))
  
  data <- read.delim(file, header = FALSE, col.names = c("kmer", "freq"), sep = "\t")

  if (nrow(data) == 0) {
    warning(paste("File is empty or has no valid data:", file))
    next
  }

  data$kmer <- toupper(data$kmer)
  data$source_file <- name

  db_list[[name]] <- data
}

load_time <- Sys.time()
print(paste("Data loaded: ", load_time))

# Combine all datasets into one
combined_data <- bind_rows(db_list)

combined_time <- Sys.time()
print(paste("Data combined: ", combined_time))

# Pivot the data to wide format
wide_data <- combined_data %>%
  pivot_wider(names_from = source_file, values_from = freq, values_fill = 0)

# Create a binary matrix
binary_matrix <- as.data.frame(ifelse(wide_data[,-1] > 0, 1, 0))

# Extract syotti_kmers data
syotti_kmers <- unique(db_list[[gsub("\\.txt$", "", basename(syotti))]]$kmer)

# Calculate overlaps
overlap_counts <- sapply(datasets, function(dataset) {
  dataset_kmers <- unique(db_list[[gsub("\\.txt$", "", basename(dataset))]]$kmer)
  length(intersect(syotti_kmers, dataset_kmers))
})

# Create a data frame for overlaps
overlap_table <- data.frame(
  Dataset = basename(datasets),
  Overlap_Count = overlap_counts,
  stringsAsFactors = FALSE
)

# Save the table to a TSV file
write.table(overlap_table, file = table, sep = "\t", row.names = FALSE, quote = FALSE)

# Print the overlap table to the console (optional)
print(overlap_table)

# Filter to keep only rows where "syotti_kmers" is involved
binary_matrix_filtered <- binary_matrix[binary_matrix$syotti_kmers == 1, ]

# Determine which datasets have overlaps with syotti_kmers
overlapping_datasets <- colnames(binary_matrix)[colSums(binary_matrix) > 0]

# Create the binary matrix only for the datasets that we want to include
filtered_binary_matrix <- binary_matrix_filtered[, overlapping_datasets]

# Combine the filtered binary matrix with the kmer column
upset_data <- cbind(kmer = wide_data$kmer[binary_matrix$syotti_kmers == 1], filtered_binary_matrix)

png("/home/abergm/bait_gen/outputs/plots/syotti_upset_plot.png", width = 1920, height = 1080)

# Create the UpSet plot for kmer overlap based on presence/absence
upset(
  upset_data,
  sets = overlapping_datasets,
  keep.order = TRUE,
  order.by = c("degree", "freq"),
  nsets = NA,
  nintersects = NA,
  text.scale = c(2, 2, 1.5, 1.5, 1.5, 1.5)
)

# Print end time for tracking duration
end_time <- Sys.time()
print(paste("End time:", end_time))

dev.off()
