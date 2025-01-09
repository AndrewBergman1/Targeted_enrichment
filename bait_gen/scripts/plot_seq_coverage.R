library(tidyverse)
library(Biostrings)

start_time <- Sys.time()
print(paste("Start time:", start_time))

# Get the command line arguments (file paths)
args <- commandArgs(trailingOnly = TRUE)

fasta_file <- args[1]
#catch_probes_file <- args[2]
syotti_probes_file <- args[2]

# Check if files exist
if (!file.exists(fasta_file)) stop("FASTA file not found: ", fasta_file)
if (!file.exists(syotti_probes_file)) stop("Syotti probes file not found: ", syotti_probes_file)
#if (!file.exists(catch_probes_file)) stop("Catch probes file not found: ", catch_probes_file)

# Read the files using readDNAStringSet
fasta <- readDNAStringSet(fasta_file, format = "fasta")
syotti_probes <- readDNAStringSet(syotti_probes_file, format = "fasta")
#catch_probes <- readDNAStringSet(catch_probes_file, format = "fasta")

# Print sample data to verify
print(paste("Number of sequences in FASTA file:", length(fasta)))
print(paste("Number of Syotti probes:", length(syotti_probes)))
#print(paste("Number of Catch probes:", length(catch_probes)))

# For bar plot on number of baits
sizes <- data.frame("syotti" = length(syotti_probes))
sizes_long <- pivot_longer(sizes, cols = c("syotti"), values_to = "Value", names_to = "method")

#catch_matches <- numeric(length(fasta))
syotti_matches <- numeric(length(fasta))

# Count the matches for Catch probes
#for (i in seq_len(length(catch_probes))) {
#  catch_match <- vcountPattern(as.character(catch_probes[i]), fasta)
#  catch_matches <- catch_matches + catch_match
#}

#catch_time <- Sys.time()
#print(paste("Catch finished at: ", catch_time))

# Count the matches for Syotti probes
for (i in seq_len(length(syotti_probes))) {
  syotti_match <- vcountPattern(as.character(syotti_probes[i]), fasta)
  syotti_matches <- syotti_matches + syotti_match
}

syotti_time <- Sys.time()
print(paste("Syotti finished at: ", syotti_time))

# Prepare coverage data for plotting
coverage <- data.frame(
  "sequence_number" = seq_len(length(fasta)),
  #"catch" = catch_matches,
  "syotti" = syotti_matches
)
coverage_long <- pivot_longer(coverage, cols = c("syotti"), values_to = "coverage", names_to = "method")

# Plot coverage
seq_cov <- ggplot(coverage_long, aes(x = sequence_number, y = log(coverage), color = method)) +
  stat_summary(fun = mean, geom = "point", size = 0.2) +  # Plots the mean of 'values'
  facet_wrap(~method) +
  ggtitle("Average Database Coverage: Syotti")

# Save the plot
ggsave(filename = "probe_sequence_coverage.png", plot = seq_cov, path = "outputs/plots/", create.dir = TRUE)

