library(tidyverse)


args <- commandArgs(trailingOnly = TRUE)

catch_probes <- args[1]


# Read the data
catch <- read_delim(catch_probes, delim = "\t",  # Adjust delimiter as needed
                    col_names = c("Query", "Target", "taxlineage", "Taxa_Lineage", "bits"))

# Summarize data
sum <- catch %>% 
  mutate(Genus = word(Taxa, 1, sep = " ")) %>% 
  group_by(Genus) %>% 
  summarise(n_probes_per_genus = n(), .groups = 'drop')

ggplot(sum, aes(x = log(n_probes_per_genus))) +
  geom_histogram(col = "white", fill = "steelblue", binwidth = 1)


sum_sum <- sum %>% 
  summarise(mean_probes_per_genus = mean(n_probes_per_genus), min_probes_per_genus = min(n_probes_per_genus), max_probes_per_genus = max(n_probes_per_genus))
