################################################################################
# Experimental infections
################################################################################
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

# Function to calculate AUC for titer values above 4 PFU
calculate_auc <- function(df, threshold = 4) {
  
  # Adjust all titer values: if below threshold, set to threshold
  df_adjusted <- df %>% mutate(titer_adjusted = pmax(titer_interpolated, threshold))
  
  # Calculate the AUC using adjusted titer values
  return(trapz(df_adjusted$day, df_adjusted$titer_adjusted - threshold))
}

# Average curves per species and day + NAs interpolation
curves_species <- viremia_experiments %>%
  dplyr::select(curve, species, host_family, host_order, day, titer) %>%
  group_by(species, host_family, host_order, day) %>%
  summarise(titer = if_else(all(is.na(titer)), NA_real_, mean(titer, na.rm = TRUE))) %>%
  mutate(titer_interpolated = na.approx(titer, day, na.rm = FALSE)) %>%
  filter(!is.na(titer_interpolated))

# Apply function to species curves
auc_species <- curves_species %>%
  group_by(species, host_family, host_order) %>%
  summarise(host_competence = calculate_auc(cur_data()), .groups = 'drop')

# Merge with PNAE list
auc_species$species_2 <- auc_species$species
PNAE_birds_WNV <- merge(PNAE_birds, auc_species, 
                        by.x = c("latin_name", "old_latin_name", "family", "order"),
                        by.y = c("species", "species_2", "host_family", "host_order"), 
                        all.x = TRUE)

writexl::write_xlsx(PNAE_birds_WNV, "Data/Birds_data/Outputs/PNAE_birds_WNV_reanalysis.xlsx")

# Viremia curves plot
curves_species$curve <- "mean"
curves_species <- curves_species %>%
  dplyr::select(curve, species, host_family, host_order, day, titer, titer_interpolated)
viremia_experiments$curve <- as.character(viremia_experiments$curve)
viremia_experiments <- rbind(viremia_experiments, curves_species)

species_curves_plot <- viremia_experiments %>%
  ggplot(aes(x = day, y = titer_interpolated, group = curve, color = ifelse(curve == "mean", "mean", "experiment"))) +
  geom_line() +
  geom_hline(yintercept = 4, linetype = "dashed", color = "red", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
  scale_color_manual(values = c("mean" = "black", "experiment" = "darkolivegreen3")) +
  facet_wrap(~species, scales = "free") +
  labs(
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
filename <- "Plots/Birds/species_curves_reanalysis.jpg"
ggsave(filename, dpi = 300, width = 10, height = 6, units = "in", type = "jpg", quality = 100)

################################################################################
# Birds analysis
################################################################################
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
library(FactoMineR)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Working with the Atlas data --------------------------------------------------
# Read data
atlas_data <- read_delim("Data/Birds_data/Data/Abundance_data_atlas.csv")
colnames(atlas_data)[colnames(atlas_data) == "latin_name"] <- "atlas_name"
study_area <- read_excel("Data/Birds_data/Data/Study_area.xlsx")
PNAE_birds_WNV <- read_excel("Data/Birds_data/Outputs/PNAE_birds_WNV_reanalysis.xlsx")
correspondencia <- PNAE_birds_WNV %>% dplyr::select(latin_name, old_latin_name)

# Filter by study site (PNAE)
atlas_data <- merge(atlas_data, study_area, by.x = "utm_1x1", by.y = "UTM1X1") %>%
  clean_names() %>%
  dplyr::select(atlas_name, utm_1x1, zona, minim, maxim) # filter by area

# Fix names
atlas_names <- unique(atlas_data$atlas_name)
atlas_names_df <- as.data.frame(atlas_names)
latin_names <- unique(PNAE_birds_WNV$latin_name)
old_latin_names <- unique(PNAE_birds_WNV$old_latin_name)

result_list <- list()

for(i in 1:nrow(atlas_names_df)) {
  atlas_name <- atlas_names_df$atlas_name[i]
  
  # Search for a match in 'old_latin_name'
  match_old <- correspondencia[correspondencia$old_latin_name == atlas_name, ]
  
  # Search for a match in 'latin_name'
  match_new <- correspondencia[correspondencia$latin_name == atlas_name, ]
  
  # If there is a match in 'old_latin_name'
  if(nrow(match_old) > 0) {
    result_list[[i]] <- data.frame(
      atlas_name = atlas_name,
      old_latin_name = match_old$old_latin_name,
      latin_name = match_old$latin_name
    )
  }
  # If there is a match in 'latin_name'
  else if(nrow(match_new) > 0) {
    result_list[[i]] <- data.frame(
      atlas_name = atlas_name,
      old_latin_name = match_new$old_latin_name,
      latin_name = match_new$latin_name
    )
  }
}

correspondencia <- do.call(rbind, result_list)
rm(atlas_names_df, match_new, match_old, result_list, atlas_name, atlas_names,
   latin_names, old_latin_names, i)

PNAE_birds_WNV_atlas <- merge(PNAE_birds_WNV, correspondencia, by = c("latin_name", "old_latin_name")) %>%
  dplyr::select(atlas_name, family, order, host_competence)

# Add viremia
atlas_data <- merge(atlas_data, PNAE_birds_WNV_atlas, by = c("atlas_name"), all.x = TRUE) %>%
  dplyr::select(atlas_name, family, order, utm_1x1, zona, minim, maxim, host_competence)

# Reservoirs and non-reservoirs abundance
atlas_data_reservoir_species <- atlas_data %>% filter(host_competence > 0)
atlas_data_non_reservoir_species <- atlas_data %>% filter(host_competence == 0)

# Cell calculations ------------------------------------------------------------
cell_pv <- atlas_data %>% 
  group_by(utm_1x1) %>% 
  summarise(predicted_abundance = sum(maxim))

# Cell calculations (reservoir species) 
cell_pv_rs <- atlas_data_reservoir_species %>%
  group_by(utm_1x1) %>%
  summarise(predicted_abundance_rs = sum(maxim),
            host_capacity_rs = sum(maxim*host_competence))

cell_pv <- merge(cell_pv, cell_pv_rs, by = "utm_1x1")

# Cell calculations (non-reservoir species)
cell_pv_nrs <- atlas_data_non_reservoir_species %>%
  group_by(utm_1x1) %>%
  summarise(predicted_abundance_nrs = sum(maxim))

cell_pv <- merge(cell_pv, cell_pv_nrs, by = "utm_1x1")

# The proportion of reservoir species determine amplification and dilution effects
cell_pv$rs_proportion <- (cell_pv$predicted_abundance_rs)/
  (cell_pv$predicted_abundance_nrs + cell_pv$predicted_abundance_rs)

write.csv(cell_pv, "Data/Birds_data/Outputs/birds_UTM1x1_reanalysis.csv")
rm(study_area, atlas_data_non_reservoir_species, atlas_data_reservoir_species, 
   cell_pv, cell_pv_nrs, cell_pv_rs, correspondencia)

# Host capacity PCA ------------------------------------------------------------
# Filter data and calculate host capacity
host_capacity <- atlas_data %>%
  filter(!is.na(host_competence)) %>%
  group_by(atlas_name, host_competence) %>%
  summarise(abundance = sum(maxim), .groups = "drop") %>%
  mutate(host_capacity = host_competence * abundance) %>%
  ungroup()

# Add species
new_rows <- data.frame(
  atlas_name = c("Bubulcus ibis", "Nycticorax nycticorax"),
  host_competence = c(0.000000, 7.2735000),
  abundance = c(100, 100),
  host_capacity = c(0.0000, 727.35000)
)

host_capacity <- rbind(host_capacity, new_rows)

# Create a dataframe for PCA (excluding categorical variables)
pca_df <- host_capacity %>%
  dplyr::select(-atlas_name, -host_capacity)

# Perform PCA
pca_result <- PCA(pca_df, scale.unit = TRUE, graph = FALSE)
summary(pca_result) 

# Visualize principal components
fviz_eig(pca_result, addlabels = TRUE)

# K-means clustering on PCA space
num_components <- 2  # Number of principal components to use
pca_data <- as.data.frame(pca_result$ind$coord[, 1:num_components]) 

# Determine the optimal number of clusters using the elbow method
fviz_nbclust(pca_data, kmeans, method = "wss")

# Apply K-means clustering
set.seed(123)
num_clusters <- 3
kmeans_model <- kmeans(pca_data, centers = num_clusters)

# Assign clusters to the original data
host_capacity$cluster <- as.factor(kmeans_model$cluster)

# Adjusted cluster numbers
original_clusters <- host_capacity$cluster

new_clusters <- rep(NA, length(original_clusters))  # Vector vacío para asignar nuevos valores
new_clusters[original_clusters == 1] <- 2  # Donde era 1, ahora es 2
new_clusters[original_clusters == 2] <- 1  # Donde era 2, ahora es 1
new_clusters[original_clusters == 3] <- 3  # Mantener el 3 sin cambios

host_capacity$cluster <- as.factor(new_clusters)

# Log transformed values
host_capacity <- host_capacity %>%
  mutate(log_abundance = log1p(abundance),  # log1p(x) = log(1 + x), evita log(0)
         log_host_competence = log1p(host_competence))

# K-means plot
kmeans_plot <- ggplot(host_capacity, aes(x = log_abundance, y = log_host_competence, label = atlas_name)) +
  geom_point(aes(color = factor(cluster)), size = 8, alpha = 0.8) +  
  geom_label_repel(aes(label = atlas_name), 
                   size = 6, max.overlaps = 20, force = 20, 
                   box.padding = 2, point.padding = 1, min.segment.length = 0.1,
                   color = ifelse(host_capacity$atlas_name %in% 
                                    c("Anas platyrhynchos", "Columba livia", "Streptopelia decaocto", "Coturnix coturnix"), 
                                  "firebrick", "black"),
                   fontface = "italic") +
  labs(x = "Log(Abundance)", y = "Log(Host Competence)") +
  theme_minimal() +
  guides(color = guide_legend(title = "Cluster", override.aes = list(size = 6))) + 
  theme_minimal(base_size = 14) + 
  theme(
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    axis.text.x = element_text(size = 14, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black", size = 0.5),  
    panel.grid.major = element_line(color = "gray90", size = 0.5),  
    panel.grid.minor = element_blank(), 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.7, "cm"),
    legend.spacing.y = unit(0.2, "cm"),  
    legend.background = element_rect(fill = "white")
  )

kmeans_plot
filename <- "Plots/Birds/kmeans_reanalysis.jpg"
ggsave(filename, dpi = 300, width = 15, height = 7, units = "in", type = "jpg", quality = 100)

################################################################################
# Risk analysis
################################################################################
library(readxl)
library(tidyverse)
library(ggplot2)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Read data
birds_UTM1x1 <- read.csv("Data/Birds_data/Outputs/birds_UTM1x1_reanalysis.csv")
culex_UTM1x1 <- read.csv("Prediction/Culex/culex_UTM1x1.csv")

risk <- merge(birds_UTM1x1, culex_UTM1x1, by.x = "utm_1x1", by.y = "UTM1X1") %>%
  dplyr::select(utm_1x1, rs_proportion, counts) %>%
  group_by(utm_1x1, rs_proportion, counts) %>%
  summarise(risk = rs_proportion*counts)

risk$risk <- (risk$risk - min(risk$risk))/(max(risk$risk) - min(risk$risk))

write.csv(risk, "Data/Risk_data/risk_UTM1x1_reanalysis.csv")

