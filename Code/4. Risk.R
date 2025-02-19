library(readxl)
library(tidyverse)
library(ggplot2)
library(vegan)
library(factoextra)
library(betareg)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Read data
birds_UTM1x1 <- read.csv("Data/Birds_data/Outputs/birds_UTM1x1.csv")
culex_UTM1x1 <- read.csv("Prediction/Culex/culex_UTM1x1.csv")

risk <- merge(birds_UTM1x1, culex_UTM1x1, by.x = "utm_1x1", by.y = "UTM1X1") %>%
  dplyr::select(utm_1x1, rs_proportion, counts) %>%
  group_by(utm_1x1, rs_proportion, counts) %>%
  summarise(risk = rs_proportion*counts)

# Scatter plot with regression line
ggplot(risk, aes(x = rs_proportion, y = counts)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Gráfico de correlación", x = "Birds ratio", y = "Culex spp. mean predicted counts")

risk$risk <- (risk$risk - min(risk$risk))/(max(risk$risk) - min(risk$risk))

write.csv(risk, "Data/Risk_data/risk_UTM1x1.csv")
rm(birds_UTM1x1, culex_UTM1x1)

# Land uses --------------------------------------------------------------------
land_covers <- read_csv("Data/Land_covers/landcovers_utm1x1.csv") %>%
  group_by(UTM1X1, DN) %>%
  summarise(area = sum(area)) %>%
  pivot_wider(names_from = DN, values_from = area, values_fill = 0)

land_covers$area_total <- rowSums(land_covers[, -1])
land_covers <- land_covers %>%
  mutate(across(1:20, ~ . / area_total)) %>%
  dplyr::select(-area_total)

land_covers$urban_area <- land_covers$`4` + land_covers$`5` + land_covers$`6`
land_covers$wetlands <- land_covers$`10` 
land_covers$rice_fields <- land_covers$`21`
land_covers$herbaceous_crops <- land_covers$`19` + land_covers$`20`

land_covers <- land_covers %>%
  dplyr::select(UTM1X1, urban_area, herbaceous_crops, wetlands, rice_fields)

# Merge ratio
risk <- merge(risk, land_covers, by.x = "utm_1x1", by.y = "UTM1X1") %>%
  dplyr::select(utm_1x1, risk, urban_area, herbaceous_crops, wetlands, rice_fields)

# If there are 0 or 1 values in the 'risk' column, transform them to the range (0, 1)
epsilon <- 0.0001  # A small value
risk$risk <- pmax(epsilon, pmin(1 - epsilon, risk$risk))

# Fit the beta model with land use proportions as explanatory variables
model <- betareg(risk ~ urban_area + herbaceous_crops + wetlands + rice_fields, data = risk)
summary(model)

# Predictions (optional): generate predicted values
risk$predicted_risk <- predict(model, type = "response")

# Visualize predicted vs observed values
ggplot(risk, aes(x = risk, y = predicted_risk)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed Risk", y = "Predicted Risk", title = "Observed vs. Predicted WNV Risk") +
  theme_minimal()


