library(Biostrings)
library(tibble)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
data <- args[1]
mm <- as.character(args[2])

fasta_file <- readDNAStringSet(data, format = "fasta")

fasta_file <- as.data.frame(fasta_file)
fasta_file <- rownames_to_column(fasta_file, var = "probe_name")

# Assuming fasta_file is your data frame
fasta_file <- fasta_file %>%
  mutate(
    A = str_count(x, "A"),
    T = str_count(x, "T"),
    C = str_count(x, "C"),
    G = str_count(x, "G"),
    Total = A + T + C + G,
    GC = (C + G) / Total,
    Tm = 64.9 + 41 * (G + C - 16.4) / Total,
    source = "Unfiltered_data"
  )

# Filter probes based on GC content
good_probes <- fasta_file %>%
  filter(GC >= 0.35 & GC <= 0.75)

good_probes <- good_probes %>% 
  mutate(source = "Filtered_data")


# Count out-of-bounds probes
probes_out_of_bounds <- nrow(fasta_file) - nrow(good_probes)

plot_df <- rbind(fasta_file, good_probes)

long_plot_df <- pivot_longer(plot_df, cols = c("GC", "Tm"), names_to = "Metric", values_to = "Value")

sample <- long_plot_df %>% sample_n(100)

q <- ggplot(long_plot_df, aes(x = source, y = Value, fill = source)) +
  geom_boxplot() +
  labs(title = "Boxplot of GC and Tm", x = "Metric", y = "Value") +
  facet_wrap(~Metric, scales = "free") +
  scale_fill_manual(values = c("Unfiltered_data" = "lightblue", "Filtered_data" = "lightgreen"))

file_path <- paste("/home/abergm/bait_gen/outputs/plots/", mm, "_GC_Tm_filter_plots.png", sep = "")

print(file_path)

ggsave(file_path, q, create.dir = TRUE)

print(paste(probes_out_of_bounds, "probes were filtered away from the original", 
            nrow(fasta_file), "probes, leaving", nrow(good_probes), "probes after filtering."))

print("The probes have been filtered based on GC%. GC% is between 35% and 75%.")

output_text <- paste(
  "Number of probes before filtering:", nrow(fasta_file), "\n",
  "Number of probes after filtering:", nrow(good_probes), "\n",
  "Number of probes removed:", probes_out_of_bounds, "\n",
  "Mean GC before filtering:", mean(fasta_file$GC), "\n",
  "Mean GC after filtering:", mean(good_probes$GC), "\n",
  "Mean Tm before filtering:", mean(fasta_file$Tm), "\n",
  "Mean Tm after filtering:", mean(good_probes$Tm), "\n",
  "IQR GC% before filtering:", IQR(fasta_file$GC), "\n",
  "IQR GC% after filtering:", IQR(good_probes$GC), "\n",
  "IQR Tm before filtering:", IQR(fasta_file$Tm), "\n",
  "IQR Tm after filtering:", IQR(good_probes$Tm), "\n"
)

# Write the output to a text file
cat(output_text, file = "/home/abergm/bait_gen/outputs/plots/GC_Tm_README.txt")

probes_to_save <- DNAStringSet(good_probes$x)
names(probes_to_save) <- good_probes$probe_name

writeXStringSet(probes_to_save, "/home/abergm/bait_gen/outputs/syotti/COMPLETED_0_SYOTTI_PROBES.fasta", format = "fasta")


