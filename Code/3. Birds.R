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

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Working with the Atlas data --------------------------------------------------
# Read data
atlas_data <- read_delim("Data/Birds_data/Data/Abundance_data_atlas.csv")
colnames(atlas_data)[colnames(atlas_data) == "latin_name"] <- "atlas_name"
study_area <- read_excel("Data/Birds_data/Data/Study_area.xlsx")
PNAE_birds_WNV <- read_excel("Data/Birds_data/Outputs/PNAE_birds_WNV.xlsx")
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
  dplyr::select(atlas_name, family, order, host_competence, host_competence_family, host_competence_order)

# Add viremia
atlas_data <- merge(atlas_data, PNAE_birds_WNV_atlas, by = c("atlas_name"), all.x = TRUE) %>%
  dplyr::select(atlas_name, family, order, utm_1x1, zona, minim, maxim, 
                host_competence, host_competence_family, host_competence_order)

abundance_ranking <- atlas_data %>%
  group_by(atlas_name, family, order) %>%
  summarise(abundancia = sum(maxim)) %>%
  arrange(-abundancia)

writexl::write_xlsx(abundance_ranking, "Data/Birds_data/Outputs/abundance_ranking.xlsx")

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

write.csv(cell_pv, "Data/Birds_data/Outputs/birds_UTM1x1.csv")

# Host capacity per species ----------------------------------------------------
host_capacity <- atlas_data %>%
  filter(!is.na(host_competence)) %>%
  group_by(atlas_name, host_competence) %>%
  summarise(abundance = sum(maxim)) %>%
  group_by(atlas_name, host_competence, abundance) %>%
  summarise(host_capacity = host_competence*abundance) %>%
  arrange(desc(host_capacity),desc(abundance))

rm(atlas_data_non_reservoir_species, atlas_data_reservoir_species, cell_pv_nrs, cell_pv_rs)

