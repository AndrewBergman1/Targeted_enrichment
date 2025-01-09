library(tidyverse)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

output_file <- args[2]
syotti_probes <- args[1]

syotti_probes <- readDNAStringSet(syotti_probes, format = "fasta")
#catch_probes <- readDNAStringSet(catch_probes, format = "fasta")

syotti_nt <- data.frame(alphabetFrequency(syotti_probes))
#catch_nt <- data.frame(alphabetFrequency(catch_probes))

syotti_nt["Source"] <- rep("syotti_probes", nrow(syotti_nt))
#catch_nt["Source"] <- rep("catch_probes", nrow(catch_nt))

#total_freqs <- rbind(syotti_nt, catch_nt)
total_freqs <- syotti_nt

total_freqs <- total_freqs[, c("A", "T", "C", "G", "Source")]

total_freqs$GC <- (total_freqs$G + total_freqs$C) / rowSums(total_freqs[, c("A", "T", "C", "G")])

total_freqs$Tm <- 64.9 + 41 * (total_freqs$G + total_freqs$C - 16.4) / rowSums(total_freqs[, c("A", "T", "C", "G")])

long_freqs <- pivot_longer(total_freqs, cols = -Source, names_to = "Metric", values_to = "Value")

Tm_gc <- long_freqs %>% filter(Metric %in% c("GC", "Tm"))

# Corrected summarization
sum <- Tm_gc %>%
  group_by(Source, Metric) %>%
  summarize(
    mean = mean(Value),
    min = min(Value),
    max = max(Value),
    .groups = "drop"  # Ensures that the data is ungrouped after summarizing
  )

sum <- sum %>%
  mutate(y_pos = max + (max - min) * 0.05)

# Updated plotting code
p <- ggplot(Tm_gc, aes(x = Source, y = Value)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
  facet_wrap(~Metric, scales = "free") +
  geom_text(
    data = sum,
    aes(
      x = Source,
      y = y_pos,
      label = paste0(
        "Mean: ", round(mean, 2), "\n",
        "Min: ", round(min, 2), "\n",
        "Max: ", round(max, 2)
      )
    ),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3,
    color = "blue"
  )

ggsave(filename = output_file, plot = p, dpi = 300, width = 10, height = 8)

