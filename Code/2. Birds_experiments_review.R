library(readxl)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(lubridate)
library(janitor)
library(vegan)
library(factoextra)
library(tibble)
library(MASS)
library(svglite)
library(pracma)
library(zoo)
library(ggthemes)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Average curve per species and subsequent AUC calculation ---------------------
# Read data
PNAE_birds <- read_excel("Data/Birds_data/Data/PNAE_birds.xlsx") %>% 
  dplyr::select(latin_name, old_latin_name, family, order)
viremia_experiments <- read_excel("Data/Birds_data/Data/Review_experimental_infections_birds_PNAE.xlsx") %>%
  clean_names()

# Editing
viremia_experiments$day <- as.numeric(viremia_experiments$day)
viremia_experiments <- viremia_experiments %>%
  filter(!is.na(day)) %>%
  group_by(curve) %>%
  mutate(
    individuals_tested = max(sample_size_original),
    survival = min(sample_size_survival)/max(sample_size_original)) %>%
  dplyr::select(species, host_family, host_order, country, wnv_strain, virus_genotype, virus_dose, 
                curve, individuals_tested, survival, day, titer) %>%
  arrange(curve, day)

# Add day 0
viremia_experiments <- viremia_experiments %>%
  group_by(curve) %>%
  arrange(curve, day) %>%
  slice(1) %>% # Take the first row of each group (first day)
  mutate(day = 0, titer = 0) %>% # Change the day to 0 and the titer to 0
  bind_rows(viremia_experiments) %>% # Combine with the original dataframe
  arrange(curve, day)

# NAs interpolation
viremia_experiments <- viremia_experiments %>%
  mutate(titer_interpolated = na.approx(titer, day, na.rm = FALSE)) %>%
  filter(!is.na(titer_interpolated)) %>%
  dplyr::select(curve, species, host_family, host_order, day, titer, titer_interpolated)

# Function to calculate AUC for titer values above 5 PFU
calculate_auc <- function(df, threshold = 5) {
  
  # Adjust all titer values: if below threshold, set to threshold
  df_adjusted <- df %>% mutate(titer_adjusted = pmax(titer_interpolated, threshold))
  
  # Calculate the AUC using adjusted titer values
  return(trapz(df_adjusted$day, df_adjusted$titer_adjusted - threshold))
}

# Average curves per species/family/order and day + NAs interpolation
curves_species <- viremia_experiments %>%
  dplyr::select(curve, species, host_family, host_order, day, titer) %>%
  group_by(species, host_family, host_order, day) %>%
  summarise(titer = if_else(all(is.na(titer)), NA_real_, mean(titer, na.rm = TRUE))) %>%
  mutate(titer_interpolated = na.approx(titer, day, na.rm = FALSE)) %>%
  filter(!is.na(titer_interpolated))

curves_family <- viremia_experiments %>%
  dplyr::select(curve, species, host_family, host_order, day, titer) %>%
  group_by(host_family, host_order, day) %>%
  summarise(titer = if_else(all(is.na(titer)), NA_real_, mean(titer, na.rm = TRUE))) %>%
  mutate(titer_interpolated = na.approx(titer, day, na.rm = FALSE)) %>%
  filter(!is.na(titer_interpolated))

curves_order <- viremia_experiments %>%
  dplyr::select(curve, species, host_family, host_order, day, titer) %>%
  group_by(host_order, day) %>%
  summarise(titer = if_else(all(is.na(titer)), NA_real_, mean(titer, na.rm = TRUE))) %>%
  mutate(titer_interpolated = na.approx(titer, day, na.rm = FALSE)) %>%
  filter(!is.na(titer_interpolated))

# Apply function to species/family/order curves
auc_species <- curves_species %>%
  group_by(species, host_family, host_order) %>%
  summarise(host_competence = calculate_auc(cur_data()), .groups = 'drop') 

auc_family <- curves_family %>%
  group_by(host_family, host_order) %>%
  summarise(host_competence_family = calculate_auc(cur_data()), .groups = 'drop') 

auc_order <- curves_order %>%
  group_by(host_order) %>%
  summarise(host_competence_order = calculate_auc(cur_data()), .groups = 'drop')

# Merge with PNAE list
auc_species$species_2 <- auc_species$species
PNAE_birds_WNV <- merge(PNAE_birds, auc_species, 
                        by.x = c("latin_name", "old_latin_name", "family", "order"),
                        by.y = c("species", "species_2", "host_family", "host_order"), 
                        all.x = TRUE)
PNAE_birds_WNV <- merge(PNAE_birds_WNV, auc_family, 
                        by.x = c("family", "order"), by.y = c("host_family", "host_order"), 
                        all.x = TRUE)
PNAE_birds_WNV <- merge(PNAE_birds_WNV, auc_order, 
                        by.x = c("order"), by.y = c("host_order"), 
                        all.x = TRUE)

PNAE_birds_WNV <- PNAE_birds_WNV %>% 
  dplyr::select(latin_name, old_latin_name, family, order, 
                host_competence, host_competence_family, host_competence_order)
writexl::write_xlsx(PNAE_birds_WNV, "Data/Birds_data/Outputs/PNAE_birds_WNV.xlsx")

# Plots ------------------------------------------------------------------------
# Plot 1: mean curves by family
# Calculate the maximum value of titer_interpolated for each host_family
family_max_titer <- curves_family %>%
  group_by(host_family) %>%
  summarise(max_titer = max(titer_interpolated, na.rm = TRUE)) %>%
  arrange(desc(max_titer))

# Reorder levels of host_family factor based on max_titer values
curves_family <- curves_family %>%
  mutate(host_family = factor(host_family, levels = family_max_titer$host_family))

# Create the plot
mean_curves_family_plot <- curves_family %>%
  filter(!day %in% c("8", "9", "10")) %>%
  ggplot(aes(x = day, y = titer_interpolated, color = host_family)) +
  geom_line() + 
  geom_point(size = 1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black", size = 0.5) +  # Horizontal line at threshold (y = 5)
  scale_x_continuous(breaks = seq(0, max(curves_family$day, na.rm = TRUE), by = 1)) +
  labs(title = "Viremia Levels by Family",
       x = "Day postinfection",
       y = "Log PFU/ml serum",
       color = "Family") +
  theme_minimal() +
  theme(legend.position = "right") # Adjust legend position if necessary

mean_curves_family_plot
filename <- "Plots/Birds/mean_curves_family.jpg"
ggsave(filename, dpi = 300, width = 10, height = 6, units = "in", type = "jpg", quality = 100)

rm(curves_family, family_max_titer, mean_curves_family_plot)

# Plot 2: mean curves by order
# Calculate the maximum value of titer_interpolated for each host_order
order_max_titer <- curves_order %>%
  group_by(host_order) %>%
  summarise(max_titer = max(titer_interpolated, na.rm = TRUE)) %>%
  arrange(desc(max_titer))

# Reorder levels of host_order factor based on max_titer values
curves_order <- curves_order %>%
  mutate(host_order = factor(host_order, levels = order_max_titer$host_order))

# Create the plot
mean_curves_order_plot <- curves_order %>%
  filter(!day %in% c("8", "9", "10")) %>%
  ggplot(aes(x = day, y = titer_interpolated, color = host_order)) +
  geom_line() + 
  geom_point(size = 1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black", size = 0.5) +  # Horizontal line at threshold (y = 5)
  scale_x_continuous(breaks = seq(0, max(curves_order$day, na.rm = TRUE), by = 1)) +
  labs(title = "Viremia Levels by Order",
       x = "Day postinfection",
       y = "Log PFU/ml serum",
       color = "Order") +
  theme_minimal() +
  theme(legend.position = "right") # Adjust legend position if necessary

mean_curves_order_plot
filename <- "Plots/Birds/mean_curves_order.jpg"
ggsave(filename, dpi = 300, width = 10, height = 6, units = "in", type = "jpg", quality = 100)

rm(order_max_titer, mean_curves_order_plot, curves_order)

# Plot 3: experiments + average curves per species
# Preparing dataframe
curves_species$curve <- "mean"
curves_species <- curves_species %>%
  dplyr::select(curve, species, host_family, host_order, day, titer, titer_interpolated)
viremia_experiments$curve <- as.character(viremia_experiments$curve)
viremia_experiments <- rbind(viremia_experiments, curves_species)

# Plot
species_curves_plot <- viremia_experiments %>%
  ggplot(aes(x = day, y = titer_interpolated, group = curve, color = ifelse(curve == "mean", "mean", "experiment"))) +
  geom_line() +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
  scale_color_manual(values = c("mean" = "black", "experiment" = "darkolivegreen3")) +
  facet_wrap(~species, scales = "free") +
  labs(
    # title = "Experiments and Average Viremia Curves per Species",
    x = "Day postinfection",
    y = "Log PFU/ml serum") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic", size = 10),
        axis.text = element_text(size = 10),  
        axis.title = element_text(size = 11),
        panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.grid.minor = element_blank())  

species_curves_plot
filename <- "Plots/Birds/species_curves.jpg"
ggsave(filename, dpi = 300, width = 10, height = 6, units = "in", type = "jpg", quality = 100)

rm(curves_species, species_curves_plot, viremia_experiments, filename, calculate_auc)

