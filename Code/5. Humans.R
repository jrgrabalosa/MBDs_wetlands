library(readxl)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(janitor)
library(mgcv)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Població estacional ----------------------------------------------------------
# Data from https://www.idescat.cat/pub/?id=epe&n=9523&geo=com:02

stational_population <- read.csv("Data/Human_data/Població_estacional_Alt_Empordà.csv")

# Rename variables
colnames(stational_population)[colnames(stational_population) == "any"] ="year"
colnames(stational_population)[colnames(stational_population) == "trimestre"] ="trimester"
colnames(stational_population)[colnames(stational_population) == "població_estacional"] ="stational_population"

stational_population <- stational_population %>%
  dplyr::select(year, T1, T2, T3, T4) %>%
  pivot_longer(!c(year),
               names_to = "trimester", values_to = "stational_population")

stational_population$year <- as.factor(stational_population$year)
stational_population$trimester <- as.factor(stational_population$trimester)

stational_population %>%
  ggplot(aes(x = trimester, y = stational_population, fill = year, group = year)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2", name = "Year") +
  labs(x = "Trimester", y = "Stational population", fill = "Year")

# Mean stational population
mean_population <- stational_population %>%
  filter(year %in% c("2015", "2016", "2017", "2018")) %>%
  group_by(trimester) %>%
  summarise(mean_population = mean(stational_population))

mean_population <- mean_population %>%
  mutate(
    initial_week = c(1, 14, 27, 40),
    final_week = c(13, 26, 39, 52)
  )

# Aedes albopictus prediction
pred <- read.csv("Prediction/Aedes_albopictus/pred.csv")

agg_newdata <- pred %>%
  dplyr::select(year, week, pp) %>% 
  group_by(year, week) %>% 
  summarise(pp = mean(pp))

agg_newdata$year <- as.factor(agg_newdata$year)
agg_newdata$week <- as.integer(as.character(agg_newdata$week))

mean_counts <- agg_newdata %>%
  dplyr::select(week, pp) %>% 
  group_by(week) %>%
  summarise(pp = mean(pp))

mean_counts <- mean_counts %>%
  mutate(pp_scaled = pp * 10000)

# Plot
albopictus_humans_plot <- ggplot() +
  geom_rect(data = mean_population,
    aes(xmin = initial_week, xmax = final_week, ymin = 0, ymax = mean_population, fill = "Human population"),
    alpha = 0.6, color = "black", size = 0.3) +  # Black border and adjusted opacity for the bars
  geom_line(
    data = mean_counts, 
    aes(x = week, y = pp_scaled, color = "Aedes albopictus"), 
    size = 1) +
  scale_y_continuous(
    name = "FTE seasonal human population",
    sec.axis = sec_axis(~ . / 10000, name = "Average predicted mosquito counts per trap")) +  # Adjusted right axis
  scale_x_continuous(
    breaks = seq(1, max(agg_newdata$week), by = 5)) +  # Define custom X-axis tick marks
  labs(x = "Week", fill = NULL, color = NULL) +
  scale_fill_manual(values = c("Human population" = "grey")) +
  scale_color_manual(values = c("Aedes albopictus" = "#d62728")) +
  theme_minimal(base_size = 14) +  
  theme(
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),  
    axis.title.y.right = element_text(size = 16, margin = margin(l = 10)), 
    axis.text.x = element_text(size = 14, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),  
    axis.line = element_line(color = "black", size = 0.5),  
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.89, 0.55),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.7, "cm"),
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_rect(fill = "white"),
    legend.box = "vertical") +  # Group items in a single box
  scale_color_manual(
    labels = expression(italic("Aedes albopictus")),
    values = c("Aedes albopictus" = "#d62728")
  )

albopictus_humans_plot

filename <- "Plots/Humans/albopictus_humans.jpg"
ggsave(filename, dpi = 300, width = 15, height = 8, units = "in", type = "jpg", quality = 100)

# Human cases ------------------------------------------------------------------
# Delete all variables.
rm(albopictus_humans_plot, agg_newdata, mean_population, stational_population, pred)

# Girona
girona_human_cases <- read_excel("Data/Human_data/Girona_human_cases.xlsx") %>%
  clean_names() %>%
  dplyr::select(-paludisme, -febre_del_nil_occidental)
girona_human_cases <- girona_human_cases[-nrow(girona_human_cases),]
girona_human_cases$any_de_notificacio <- as.integer(girona_human_cases$any_de_notificacio)
girona_human_cases$setmana_epidemiologica <- as.integer(girona_human_cases$setmana_epidemiologica)

# Create complete week sequence
complete_weeks <- expand_grid(
  any_de_notificacio = 2015:2023,
  setmana_epidemiologica = 1:53
)

girona_human_cases <- complete_weeks %>%
  left_join(girona_human_cases, by = c("any_de_notificacio", "setmana_epidemiologica")) %>%
  mutate(across(c(dengue, zika, chikungunya), ~replace_na(., 0)))

girona_human_cases$aedes_borne <- girona_human_cases$dengue + girona_human_cases$chikungunya + girona_human_cases$zika

girona_human_cases_mean <- girona_human_cases %>%
  group_by(setmana_epidemiologica) %>%
  summarise(aedes_borne_girona = mean(aedes_borne))

# Barcelona
barcelona_human_cases <- read_excel("Data/Human_data/Barcelona_human_cases.xlsx") %>%
  clean_names() %>%
  dplyr::select(-paludisme, -febre_del_nil_occidental)
barcelona_human_cases <- barcelona_human_cases[-nrow(barcelona_human_cases),]
barcelona_human_cases$any_de_notificacio <- as.integer(barcelona_human_cases$any_de_notificacio)
barcelona_human_cases$setmana_epidemiologica <- as.integer(barcelona_human_cases$setmana_epidemiologica)

barcelona_human_cases <- complete_weeks %>%
  left_join(barcelona_human_cases, by = c("any_de_notificacio", "setmana_epidemiologica")) %>%
  mutate(across(c(dengue, zika, chikungunya), ~replace_na(., 0)))

barcelona_human_cases$aedes_borne <- barcelona_human_cases$dengue + barcelona_human_cases$chikungunya + barcelona_human_cases$zika

barcelona_human_cases_mean <- barcelona_human_cases %>%
  group_by(setmana_epidemiologica) %>%
  summarise(aedes_borne_barcelona = mean(aedes_borne))

# Catalunya
catalunya_human_cases <- read_excel("Data/Human_data/Catalunya_human_cases.xlsx") %>%
  clean_names() %>%
  dplyr::select(-paludisme, -febre_del_nil_occidental)
catalunya_human_cases <- catalunya_human_cases[-nrow(catalunya_human_cases),]
catalunya_human_cases$any_de_notificacio <- as.integer(catalunya_human_cases$any_de_notificacio)
catalunya_human_cases$setmana_epidemiologica <- as.integer(catalunya_human_cases$setmana_epidemiologica)

catalunya_human_cases <- complete_weeks %>%
  left_join(catalunya_human_cases, by = c("any_de_notificacio", "setmana_epidemiologica")) %>%
  mutate(across(c(dengue, zika, chikungunya), ~replace_na(., 0)))

catalunya_human_cases$aedes_borne <- catalunya_human_cases$dengue + catalunya_human_cases$chikungunya + catalunya_human_cases$zika

catalunya_human_cases_mean <- catalunya_human_cases %>%
  group_by(setmana_epidemiologica) %>%
  summarise(aedes_borne_catalunya = mean(aedes_borne))

human_cases_mean <- merge(catalunya_human_cases_mean, barcelona_human_cases_mean, by = "setmana_epidemiologica")
human_cases_mean <- merge(human_cases_mean, girona_human_cases_mean, by = "setmana_epidemiologica")

rm(barcelona_human_cases, catalunya_human_cases, girona_human_cases,
   barcelona_human_cases_mean, girona_human_cases_mean, complete_weeks)

# Combine human data + albopictus
mean_counts <- mean_counts %>%
  mutate(pp_scaled = pp * 2)

combined_data <- merge(human_cases_mean, mean_counts, by.x = "setmana_epidemiologica", by.y = "week") %>%
  pivot_longer(
    !c(setmana_epidemiologica), 
    names_to = "name", values_to = "value"
  ) %>%
  mutate(name = factor(name, levels = c("aedes_borne_catalunya", "aedes_borne_barcelona", "aedes_borne_girona", "pp_scaled"))) %>%
  arrange(name, setmana_epidemiologica) %>%
  filter(!(name == "pp"))

# Palette
my_palette <- c(
  "aedes_borne_catalunya" = "dodgerblue4",
  "aedes_borne_barcelona" = "seashell2",
  "aedes_borne_girona" = "sandybrown",
  "pp_scaled" = "#D62728"
)

# Plot
aedes_plot <- ggplot(data = combined_data, 
                     aes(x = setmana_epidemiologica, y = value, group = name)) +
  geom_bar(
    data = subset(combined_data, name %in% c("aedes_borne_catalunya", "aedes_borne_barcelona", "aedes_borne_girona")), 
    aes(fill = name), 
    stat = "identity", 
    position = position_dodge(width = 0.8),
    color = "black", 
    size = 0.2) +
  geom_line(
    data = subset(combined_data, name == "pp_scaled"), 
    aes(color = name), 
    size = 1) +
  scale_fill_manual(
    values = my_palette,
    labels = c(
      "aedes_borne_catalunya" = "Catalonia",
      "aedes_borne_barcelona" = "Barcelona",
      "aedes_borne_girona" = "Girona")) +
  scale_color_manual(
    values = my_palette, # Usa la misma paleta para la línea
    labels = c("pp_scaled" = expression(italic("Aedes albopictus")))) +
  labs(
    x = "Week", 
    y = "Average predicted mosquito counts per trap", 
    fill = "Region", 
    color = "Name") +
  scale_y_continuous(
    name = "Average number of imported human cases",
    sec.axis = sec_axis(~ . / 2, name = "Average predicted mosquito counts per trap")) +
  scale_x_continuous(
    breaks = seq(1, max(human_cases_mean$setmana_epidemiologica), by = 5)) +
  # Tema personalizado
  theme_minimal(base_size = 14) +  
  theme(
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),  
    axis.title.y.right = element_text(size = 16, margin = margin(l = 10)), 
    axis.text.x = element_text(size = 14, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),  
    axis.line = element_line(color = "black", size = 0.5),  
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    # legend.position = "none",
    legend.position = c(0.9, 0.6),
    legend.text = element_text(size = 14, hjust = 0),
    legend.key.size = unit(0.7, "cm"),
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_rect(fill = "white"),
    legend.box = "vertical"
  )

aedes_plot
filename <- "Plots/Humans/albopictus_human_cases.jpg"
ggsave(filename, dpi = 300, width = 15, height = 8, units = "in", type = "jpg", quality = 100)

rm(aedes_plot, filename, my_palette)

