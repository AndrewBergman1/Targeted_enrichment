library(UpSetR)
library(tidyverse)

start_time <- Sys.time()
print(paste("Start time:", start_time))

# List all .txt files in the directory
datasets <- list.files("/home/abergm/bait_gen/outputs/kmer_freq/datasets/", pattern = "\\.txt$", full.names = TRUE)
syotti <- "/home/abergm/bait_gen/outputs/kmer_freq/syotti/syotti_kmers.txt"
catch <- "/home/abergm/bait_gen/outputs/kmer_freq/catch/catch_kmers.txt"

# Combine syotti and catch with datasets
all_files <- c(datasets, syotti, catch)

# Initialize an empty list to store the datasets
db_list <- list()

# Loop through each file, read it, and add the file name as a new column
for (file in all_files) {
  # Extract the base file name (without directory and extension)
  name <- gsub("\\.txt$", "", basename(file))

  # Read the file (assuming it's a tabular dataset, like CSV or TSV)
  data <- read.delim(file, header = FALSE, col.names = c("kmer", "freq"), sep = "\t")  # Adjust separator as needed

  # Check if the file is empty or has no rows
  if (nrow(data) == 0) {
    warning(paste("File is empty or has no valid data:", file))
    next  # Skip to the next file
  }

  # Convert kmers to uppercase to ensure consistency in matching
  data$kmer <- toupper(data$kmer)

  # Add a new column with the file name
  data$source_file <- name

  # Store the modified dataset in the list, using the base name of the file
  db_list[[name]] <- data
}

load_time <- Sys.time()
print(paste("Data loaded: ", load_time))

# Combine all datasets into one
combined_data <- bind_rows(db_list)

combined_time <- Sys.time()
print(paste("Data combined: ", combined_time))

# Pivot the data to wide format: kmers as rows, datasets as columns
wide_data <- combined_data %>%
  pivot_wider(names_from = source_file, values_from = freq, values_fill = 0)

# Create a binary matrix: 1 if the kmer is present in the dataset (freq > 0), 0 otherwise
binary_matrix <- as.data.frame(ifelse(wide_data[,-1] > 0, 1, 0))  # Exclude kmer column and convert to binary presence/absence

# Filter out catch_kmers and ensure syotti_kmers is included
datasets_to_include <- setdiff(colnames(binary_matrix), "syotti_kmers")

# Filter the binary matrix to keep only rows where "syotti_kmers" is involved (i.e., set to 1)
binary_matrix_filtered <- binary_matrix %>% filter(catch_kmers == 1)

# Create the binary matrix only for the datasets that we want to include
filtered_binary_matrix <- binary_matrix_filtered[, datasets_to_include]

# Combine the filtered binary matrix with the kmer column
upset_data <- cbind(kmer = wide_data$kmer[binary_matrix$catch_kmers == 1], filtered_binary_matrix)

png("/home/abergm/bait_gen/outputs/plots/catch_upset_plot.png", width = 1920, height = 1080)

# Create the UpSet plot for kmer overlap based on presence/absence
upset(
  upset_data,
  sets = datasets_to_include,  # Use the filtered datasets (including syotti_kmers)
  keep.order = FALSE,
  order.by = c("degree", "freq"),
  nsets = NA,
  nintersects = NA,
  text.scale = c(2, 2, 1.5, 1.5, 1.5, 1.5)
)

# Print end time for tracking duration
end_time <- Sys.time()
print(paste("End time:", end_time))


dev.off()
