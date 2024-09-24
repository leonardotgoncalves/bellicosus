# Set working directory
setwd("C:/")

############### Processing pixy output ###############

# Turning off scientific notation
# To turn it back on, change the value to zero
options(scipen = 999)

# ATTENTION: Must configure and rerun this for each species separately

#Important variables
# Species name
spp <- "PASC"
# Output from pixy, filename must be spp code + "_pixy_pi_10kb.txt"
data <- na.omit(read.delim(paste0(spp, "_pixy_pi_10kb.txt")))
# Window size (bp)
win_size <- 10000
# Output filename
filename <- paste0(spp, "_pixy_ROH_10kb.txt")

# Keep only lines that analyzed more than half of the window size
data <- subset(data, no_sites > win_size*0.5)

# Extract unique chromosome names and sort them alphabetically
unique_strings <- sort(unique(data$chromosome))

# Create a numbered index for each unique string
numbered_strings <- setNames(seq_along(unique_strings), unique_strings)

# Replace the strings in the original dataframe with their corresponding numbers
# ATTENTION: This will overwrite the previous pixy file loaded
data$chromosome <- numbered_strings[data$chromosome]
write.table(data, file = paste0(spp, "_pixy_pi_10kb.txt"), quote = F, sep = "\t")

# Calculating mean heterozygosity for species loaded
meanhet <- sum(data$count_diffs, na.rm=T)/sum(data$no_sites, na.rm=T)
cat("Mean heterozygosity for", spp, ":", round(meanhet*1000, digits = 2), "per kb")

# To run the following code, chromosome names must be a number (and only a number)
# Previous section of the code must've fixed that
# Iterate through each chromossome to search for runs of homozygosity

for (k in 1:length(unique(data$chromosome))){
  # Create a logical vector indicating values lower than 0.0005
  lower_than_threshold <- data$avg_pi[data$chromosome == k] < 0.0005

  # Compute lengths of consecutive runs of TRUE values
  run_lengths <- rle(lower_than_threshold)

  # Filter out runs that are TRUE and get their lengths
  consecutive_runs <- run_lengths$lengths[run_lengths$values]*win_size

  # Print the number of consecutive runs and their sizes
  cat("Chromosome:", k, "\n")
  cat("Number of consecutive runs:", length(consecutive_runs), "\n")
  cat("Sizes of consecutive runs:", consecutive_runs, "\n")
  cat("Percentage of windows with homozygosis:", length(lower_than_threshold)*100/length(data$window_pos_1), "\n")

  # If it was detected a ROH in the chromosome...
  if (length(consecutive_runs) > 0) {
  # Create a data frame with win_size repeated and consecutive_runs
    output_data <- data.frame(species = spp,
                              chromosome = k,
                              run_length = consecutive_runs)
  
  # Append the data frame to the file
  write.table(output_data, file = filename, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
}
}

# Now, outside R, I have concatenated the ROH files generated into a tab-separated file (ROH_count_all.txt)

############### Plotting ROH ###############

library(ggplot2)
library(dplyr)

# Read plot data
# Must be a tab-separated file with columns named species and count
plot_data <- read.delim("ROH_count_all.txt", header = T)

# Plot ROH of pascuorum
ROH.PASC <- filter(plot_data, species == "PASC") %>%
  ggplot(aes(x=count/1000, fill = species)) +
  scale_y_continuous(trans="sqrt", breaks = seq(0,600, by = 150)) +
  geom_bar() + theme_bw() +
  xlim(0,670) +
  scale_fill_manual(values="#DC6464") +
  labs(x="ROH length (kb)", y = "Frequency") +
  theme(legend.position = "none")

# Plot ROH of terrestris
ROH.TERR <- filter(plot_data, species == "TERR") %>%
  ggplot(aes(x=count/1000, fill = species)) + geom_bar() +
  scale_y_continuous(trans="sqrt", breaks = seq(0,600, by = 150)) +
  xlim(0,670) +
  theme_bw() +
  scale_fill_manual(values="#B778B3") +
  labs(x="ROH length (kb)", y = "Frequency") +
  theme(legend.position = "none")

# Plot ROH of bellicosus
ROH.BELI4 <- filter(plot_data, species == "BELI4") %>%
  ggplot(aes(x=count/1000, fill = species)) + geom_bar() +
  scale_y_continuous(trans="sqrt", breaks = seq(0,600, by = 150)) +
  xlim(0,670) +
  theme_bw() +
  scale_fill_manual(values="#5496CE") +
  labs(x="ROH length (kb)", y = "Frequency") +
  theme(legend.position = "none")


############### Plotting heterozygosity ###############
## Must run all species together because I've "reutilized"
## some of the variables

# Plot heterozygosity of terrestris (mean het. = 2.90 per kb)
het <- na.omit(read.delim("TERR_pixy_pi_10kb.txt"))
het <- subset(het, no_sites > 5000)

label_text <- het %>%
  group_by(chromosome) %>%
  summarise(
    lab = unique(chromosome),
    x = mean(window_pos_1),
    y = -1.5
  )

ggplot(het, aes(y=(avg_pi*1000), x = window_pos_1)) + 
  geom_point(size = 0.1, alpha = 0.5, col ="#B778B3") +
  scale_fill_manual(name=NULL)  +
  facet_wrap(~chromosome, 
             scales = "free_x", 
             strip.position = "bottom", 
             ncol = length(unique(het$chromosome))) +
  labs(x="Chromosome",
       y="Heterozygosity per kb",
       title = "Buff-tailed bumblebee",
       subtitle = "H = 2.90 per kb") +
  geom_text(data = label_text, mapping = aes(x = x, y = y, label = lab), vjust = 0.5, hjust = 0.5) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.spacing = unit(0, "mm")) -> het.TERR

# Plot heterozygosity of bellicosus (mean het. = 2.72 per kb)
het <- na.omit(read.delim("BELI4_pixy_pi_10kb.txt"))
het <- subset(het, no_sites > 5000)

label_text <- het %>%
  group_by(chromosome) %>%
  summarise(
    lab = unique(chromosome),
    x = mean(window_pos_1),
    y = -1.5
  )

ggplot(het, aes(y=(avg_pi*1000), x = window_pos_1)) + 
  geom_point(size = 0.1, alpha = 0.5, col ="#5496CE") +
  scale_fill_manual(name=NULL)  +
  facet_wrap(~chromosome, 
             scales = "free_x", 
             strip.position = "bottom", 
             ncol = length(unique(het$chromosome))) +
  labs(x="Chromosome",
       y="Heterozygosity per kb",
       title = "Bellicose bumblebee",
       subtitle = "H = 2.72 per kb") +
  geom_text(data = label_text, mapping = aes(x = x, y = y, label = lab), vjust = 0.5, hjust = 0.5) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.spacing = unit(0, "mm")) -> het.BELI4

# Plot heterozygosity of pascuorum (mean het. = 3.88 per kb)
het <- na.omit(read.delim("PASC_pixy_pi_10kb.txt"))
het <- subset(het, no_sites > 5000)

label_text <- het %>%
  group_by(chromosome) %>%
  summarise(
    lab = unique(chromosome),
    x = mean(window_pos_1),
    y = -1.5
  )

ggplot(het, aes(y=(avg_pi*1000), x = window_pos_1)) + 
  geom_point(size = 0.1, alpha = 0.5, col ="#DC6464") +
  scale_fill_manual(name=NULL)  +
  facet_wrap(~chromosome, 
             scales = "free_x", 
             strip.position = "bottom", 
             ncol = length(unique(het$chromosome))) +
  labs(x="Chromosome",
       y="Heterozygosity per kb",
       title = "Common carder bee",
       subtitle = "H = 3.88 per kb") +
  geom_text(data = label_text, mapping = aes(x = x, y = y, label = lab), vjust = 0.5, hjust = 0.5) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.spacing = unit(0, "mm")) -> het.PASC

############### Make a panel with the ROH and heterozygosity plots (Figure 2) ###############

library(patchwork)
library(svglite)

svglite(filename = "panel_chr_ROH.svg", width = 8.7, height = 6.5)

((het.BELI4 | ROH.BELI4) + plot_layout(widths = c(6,1))) / 
((het.PASC | ROH.PASC) + plot_layout(widths = c(6,1))) /
((het.TERR | ROH.TERR) + plot_layout(widths = c(6,1))) + 
plot_annotation(tag_levels = list(c("a", "d", "b", "e", "c", "f"))) + plot_layout(axis_titles = "collect")

dev.off()

# Must edit the output in a separate software to fix the chromosome labels
# (couldn't fix it in R)
