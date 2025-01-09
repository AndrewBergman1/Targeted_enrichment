library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

data_path <- args[1]
plot_path <- args[2]

# Read the data
data <- read.table(data_path, header = TRUE)  # Use header = TRUE to directly use the first row as column names
data$Overlap_Count <- as.numeric(data$Overlap_Count)

# Sort by Overlap_Count in descending order
data <- data %>% arrange(desc(Overlap_Count))

# Create the plot
q <- ggplot(data, aes(x = reorder(Dataset, -Overlap_Count), y = log(Overlap_Count))) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Save the plot
ggsave(plot_path, q)