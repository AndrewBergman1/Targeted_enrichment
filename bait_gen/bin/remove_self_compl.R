library(dplyr)
library(tidyverse)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

probes <- args[1]
good_probes <- args[2]
output <- args[3]

probes <- readDNAStringSet(probes)
good_probes <- read.csv(good_probes)

# Extract the IDs of problematic probes
good_ids <- good_probes$ID

# Create a logical index to keep sequences that are not in problematic IDs
keep_indices <- (names(probes) %in% good_ids)

# Subset the Biostrings object to remove problematic probes
filtered_probes <- probes[keep_indices]

# Print or save the filtered Biostrings object
print(filtered_probes)


writeXStringSet(filtered_probes, output)
