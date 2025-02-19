library(readxl)
library(writexl)
library(tidyverse)
library(lubridate)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(janitor)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Loading data
loc.data <- paste0(getwd(), "/Data/")
mosquito_data <- read_excel(paste0(loc.data, 'Mosquito_data/Mosquito_Surveillance_SCM/Surveilance2001_2021_DataSCM.xlsx'),
                            skip = 3, col_names = TRUE , col_types = c("guess", "guess","guess",
                                                                       "guess", "guess", "guess", "guess","guess", "guess",
                                                                       "guess", "guess", "guess", "guess", "guess", "guess",
                                                                       "date", "guess", "guess", "date","guess", "guess",
                                                                       "guess", "guess", "guess", "guess", "guess","text",
                                                                       "guess", "text", "text", "guess")) %>% # hay que cambiar los números a texto 
  clean_names() %>%
  dplyr::select(trap_name, start_date, end_date, latitude, longitude, species, females, comment_trap) %>% 
  filter(females != "Not available") %>% 
  transform( # volvemos a transformar el texto a número
    females = as.numeric(females),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  ) %>%
  pivot_wider(names_from = species, values_from = females, values_fill = 0) %>% # Creamos las columnas para todas las especies y añadimos 0 en el caso de NULL values
  clean_names() %>%
  dplyr::select(trap_name, start_date, end_date, latitude, longitude, comment_trap, 
                culex_sp_pipiens, culex_theileri, culex_modestus, culex_sp, 
                anopheles_atroparvus, anopheles_maculipennis, aedes_albopictus) %>%
  rowwise()

# Anti-mosquito treatment
mosquito_data$comment_trap <- ifelse(grepl("treatment", mosquito_data$comment_trap), 1, 0)
colnames(mosquito_data)[colnames(mosquito_data) == "comment_trap"] ="AMT"

# Pivot longer
mosquito_data <- mosquito_data %>%
  pivot_longer(!c(trap_name, start_date, end_date, latitude, longitude, AMT), 
               names_to = "species", values_to = "females")

# Trapping effort
mosquito_data$trapping_effort <- mosquito_data$end_date - mosquito_data$start_date
mosquito_data$trapping_effort <- as.numeric(mosquito_data$trapping_effort)
mosquito_data=mosquito_data[-which(mosquito_data$trapping_effort>8),]

# Filter rows - SPECIES
culex_df = mosquito_data %>% 
  filter(species %in% c("culex_sp_pipiens", "culex_modestus", "culex_theileri", "culex_sp"))

# Land covers in 250m buffer
buffer_250 <- read_csv("Data/Land_covers/landcovers_buffer250.csv")
buffer_250 <- buffer_250 %>%
  group_by(trap_name, DN) %>%
  summarise(area = sum(area)) %>%
  pivot_wider(names_from = "DN", values_from = "area", values_fill = 0) %>%
  dplyr::select(trap_name, `2`, `4`, `5`,`6`,`7`,`9`,`10`,`11`,`14`,`15`,`16`,`17`,`19`,
                `20`,`21`,`22`,`23`,`24`)

buffer_250$'Urban area' <- buffer_250$`4` + buffer_250$`5` + buffer_250$`6` + buffer_250$`7`
buffer_250$'Wetland' <- buffer_250$`2` + buffer_250$`10`
buffer_250$'Herbaceous crops' <- buffer_250$`19` + buffer_250$`20`
buffer_250$'Rice fields' <- buffer_250$`21`

buffer_250 <- buffer_250 %>%
  dplyr::select('Urban area', 'Wetland', 'Herbaceous crops', 'Rice fields',
                `9`,`11`,`14`,`15`,`16`,`17`,`22`,`23`,`24`) %>%
  pivot_longer(!c(trap_name), names_to = "land_cover", values_to = "area") %>%
  arrange(trap_name, desc(area))

buffer_250 <- buffer_250 %>%
  group_by(trap_name) %>%
  mutate(has_rice_field = any(land_cover == "Rice fields" & area > 0)) %>%
  filter(ifelse(has_rice_field, land_cover == "Rice fields" & area > 0, area == max(area))) %>%
  ungroup() %>%
  dplyr::select(-has_rice_field, -area) 

traps <- read_excel("Data/Mosquito_data/traps.xlsx")
traps <- merge(traps, buffer_250, by = "trap_name") %>%
  dplyr::select(trap_name, trap_type, land_cover)

# Add trap_type and land_use variables
culex_df <- merge(culex_df, traps, by = "trap_name")
names(culex_df)[names(culex_df) == "land_cover"] <- "land_use"

# Counts/day per species and land use
counts_day <- culex_df %>%
  group_by(species, land_use) %>%
  summarise(trapping_effort = sum(trapping_effort),
            females = sum(females),
            counts_day = females / trapping_effort) %>%
  dplyr::select(species, land_use, counts_day) %>%
  mutate(counts_day = round(counts_day, digits = 2)) %>% 
  pivot_wider(names_from = land_use, values_from = counts_day)

# Species contribution to the total Culex abundance by land use
contribution <- culex_df %>%
  group_by(land_use) %>%
  mutate(total_culex_land_use = sum(females)) %>%  
  group_by(species, land_use) %>% 
  summarise(females = sum(females),  
            contribution = females*100 / first(total_culex_land_use)) %>% 
  dplyr::select(species, land_use, contribution) %>%
  mutate(contribution = round(contribution, digits = 2)) %>% 
  pivot_wider(names_from = land_use, values_from = contribution) %>%
  mutate(species = dplyr::recode(species,
                          "culex_sp_pipiens" = "Cx. pipiens",
                          "culex_sp" = "Culex sp.",
                          "culex_modestus" = "Cx. modestus",
                          "culex_theileri" = "Cx. theileri"))

contribution_long <- contribution %>%
  pivot_longer(cols = starts_with("Herbaceous crops"):starts_with("Wetland"), 
               names_to = "land_use", 
               values_to = "contribution") %>%
  group_by(land_use) %>%
  arrange(land_use, desc(contribution)) %>%
  mutate(species = factor(species, levels = unique(species))) %>%
  ungroup()

# Bar graph
colors <- c("snow3", brewer.pal(n = 4, name = "Set2"))
contribution_plot <- ggplot(
  contribution_long, aes(x = land_use, y = contribution, fill = species)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  labs(title = "Species contribution to the total Culex abundance by habitat",
       fill = "Species",
       x = "Land use",
       y = "Contribution (%)") +
  scale_fill_manual(values = colors) + 
  theme_minimal() +
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        axis.text.x = element_text(size = 14, color = "black"),  
        axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

contribution_plot

filename <- "Plots/Culex/contribution_plot.jpg"
ggsave(filename, dpi = 300, width = 12, height = 7, units = "in", type = "jpg", quality = 100)

