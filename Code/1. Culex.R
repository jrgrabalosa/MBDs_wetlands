################################################################################
#DATAFRAME
################################################################################
library(readxl)
library(writexl)
library(tidyverse)
library(lubridate)
library(dplyr)
library(ggplot2)
library(purrr)
library(zoo)
library(parallelly)
library(parallel)
library(janitor)
library(pollen)

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

# Sum all the species
culex_df <- culex_df %>%
  group_by(trap_name, start_date, end_date, latitude, longitude, AMT, trapping_effort) %>%
  summarise(species="culex",females=sum(females))

# Add year, month and week columns
culex_df = culex_df %>% mutate(year=year(start_date))
culex_df$year <- as.factor(culex_df$year)
culex_df = culex_df %>% mutate(month=month(start_date))
culex_df$month <- as.factor(culex_df$month)
culex_df = culex_df %>% mutate(week=week(start_date))
culex_df$week <- as.factor(culex_df$week)

# Raw data plot
culex_plot <- culex_df %>%
  group_by(week, year, trap_name) %>%
  summarise(femalessum = sum(females)) %>%
  drop_na(femalessum) %>%
  group_by(trap_name) %>% 
  filter(any(femalessum != 0)) %>%
  ggplot(aes(x = week, y = femalessum,
             color = trap_name, group = trap_name)) +
  geom_point(size = 0.25) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~year, scales = "free") +
  ylab("Counts per trap") +
  xlab("Time in weeks") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(color = "Trap name")

filename <- "Plots/Culex/raw_data_plot.jpg"
ggsave(filename, dpi = 300, width = 8, height = 6, units = "in", type = "jpg", quality = 100)

# Meteo data
w1 <- read.csv("Data/Meteo_data/Meteo_data_W1.csv", sep = ";") %>%
  dplyr::select(DATA, TM, TX, TN, HRM, PPT) %>%
  mutate(
    DATA = as.POSIXct(DATA, tz = "GMT", format = "%d/%m/%Y"))

u2 <- read.csv("Data/Meteo_data/Meteo_data_U2.csv", sep = ";") %>%
  dplyr::select(DATA, TM, TX, TN, HRM, PPT) %>%
  mutate(
    DATA = as.POSIXct(DATA, tz = "GMT", format = "%d/%m/%Y"))

meteo_data <- w1 %>%
  left_join(u2, by = "DATA") %>%
  mutate(
    TM = ifelse(is.na(TM.x), TM.y, TM.x),
    TX = ifelse(is.na(TX.x), TX.y, TX.x),
    TN = ifelse(is.na(TN.x), TN.y, TN.x),
    HRM = ifelse(is.na(HRM.x), HRM.y, HRM.x), 
    PPT = ifelse(is.na(PPT.x), PPT.y, PPT.x), ) %>%
  dplyr::select(DATA, TM, TX, TN, HRM, PPT)

colnames(meteo_data)[colnames(meteo_data) == "DATA"] ="start_date"

# Mosquito weather index
meteo_data = meteo_data %>% mutate(
  FHme = case_when(HRM < 40~0, HRM  >95~0,
                   (HRM  >=40 & HRM <= 95)~
                     ((HRM/55)-(40/55)) ),
  FTme = case_when(TM<=15~0, TM>30~0, 
                   (TM>15 & TM <=20)~
                     (.2*TM)-3,
                   (TM>20 & TM<=25)~1,
                   (TM>25 & TM <= 30)~
                     (-.2*TM)+6),
  mwime = FHme*FTme)

# 7, 14 and 21-day cumulative weather variables
meteo_data <- meteo_data %>% 
  mutate(TM7 = rollmean(TM, k=7, align = "right", fill = NA),
         TM14 = rollmean(TM, k=14, align = "right", fill = NA),
         TM21 = rollmean(TM, k=21, align = "right", fill = NA),
         TX7 = rollmean(TX, k=7, align = "right", fill = NA),
         TX14 = rollmean(TX, k=14, align = "right", fill = NA),
         TX21 = rollmean(TX, k=21, align = "right", fill = NA),
         TN7 = rollmean(TN, k=7, align = "right", fill = NA),
         TN14 = rollmean(TN, k=14, align = "right", fill = NA),
         TN21 = rollmean(TN, k=21, align = "right", fill = NA),
         HRM7 = rollmean(HRM, k=7, align = "right", fill = NA),
         HRM14 = rollmean(HRM, k=14, align = "right", fill = NA),
         HRM21 = rollmean(HRM, k=21, align = "right", fill = NA),
         PPT7 = rollapplyr(PPT, 7, FUN=sum, align = "right", fill = NA),
         PPT14 = rollapplyr(PPT, 14, FUN=sum, align = "right", fill = NA),
         PPT21 = rollapplyr(PPT, 21, FUN=sum, align = "right", fill = NA), 
         mwi7 = rollmean(mwime, k = 7, align = "right", fill = NA),
         mwi14 = rollmean(mwime, k = 14, align = "right", fill = NA),
         mwi21 = rollmean(mwime, k=21, align = "right", fill = NA))

# Loop for averaging weather variables
for(i in 1:nrow(culex_df)){ 
  x0 <- culex_df$start_date[i]
  x1 <- culex_df$end_date[i]
  if(is.na(x0) || is.na(x1)){
    next
  }
  int <- seq(x0, x1, by= "day")
  
  TM <- mean(meteo_data$TM[meteo_data$start_date %in% int])
  TX <- mean(meteo_data$TX[meteo_data$start_date %in% int])
  TN <- mean(meteo_data$TN[meteo_data$start_date %in% int])
  HRM <- mean(meteo_data$HRM[meteo_data$start_date %in% int])
  PPT <- sum(meteo_data$PPT[meteo_data$start_date %in% int])
  mwime <- mean(meteo_data$mwime[meteo_data$start_date %in% int])
  TM7 <- meteo_data$TM7[meteo_data$start_date %in% x0]
  TM14 <- meteo_data$TM14[meteo_data$start_date %in% x0]
  TM21 <- meteo_data$TM21[meteo_data$start_date %in% x0]
  TX7 <- meteo_data$TX7[meteo_data$start_date %in% x0]
  TX14 <- meteo_data$TX14[meteo_data$start_date %in% x0]
  TX21 <- meteo_data$TX21[meteo_data$start_date %in% x0]
  TN7 <- meteo_data$TN7[meteo_data$start_date %in% x0]
  TN14 <- meteo_data$TN14[meteo_data$start_date %in% x0]
  TN21 <- meteo_data$TN21[meteo_data$start_date %in% x0]
  HRM7 <- meteo_data$HRM7[meteo_data$start_date %in% x0]
  HRM14 <- meteo_data$HRM14[meteo_data$start_date %in% x0]
  HRM21 <- meteo_data$HRM21[meteo_data$start_date %in% x0]
  PPT7 <- meteo_data$PPT7[meteo_data$start_date %in% x0]
  PPT14 <- meteo_data$PPT14[meteo_data$start_date %in% x0]
  PPT21 <-meteo_data$PPT21[meteo_data$start_date %in% x0]
  mwi7 <- meteo_data$mwi7[meteo_data$start_date %in% x0]
  mwi14 <- meteo_data$mwi14[meteo_data$start_date %in% x0]
  mwi21 <- meteo_data$mwi21[meteo_data$start_date %in% x0]
  
  culex_df$TM[i] <- TM
  culex_df$TX[i] <- TX
  culex_df$TN[i] <- TN
  culex_df$HRM[i] <- HRM
  culex_df$PPT[i] <- PPT
  culex_df$mwime[i] <- mwime
  culex_df$TM7[i] <- TM7
  culex_df$TM14[i] <- TM14
  culex_df$TM21[i] <- TM21
  culex_df$TX7[i] <- TX7
  culex_df$TX14[i] <- TX14
  culex_df$TX21[i] <- TX21
  culex_df$TN7[i] <- TN7
  culex_df$TN14[i] <- TN14
  culex_df$TN21[i] <- TN21
  culex_df$HRM7[i] <- HRM7
  culex_df$HRM14[i] <- HRM14
  culex_df$HRM21[i] <- HRM21
  culex_df$PPT7[i] <- PPT7
  culex_df$PPT14[i] <- PPT14
  culex_df$PPT21[i] <- PPT21
  culex_df$mwi7[i] <- mwi7
  culex_df$mwi14[i] <- mwi14
  culex_df$mwi21[i] <- mwi21
}

# Growing degree days
culex_df <-  bind_rows(mclapply(1:nrow(as.data.frame(culex_df)), function(i){
  # Verbose
  if (i %% 50 == 0) {print(i)}
  
  x = culex_df[i,]
  this_start_date = x$start_date
  this_end_date = x$end_date
  
  # Exposure week
  wth_v <- meteo_data %>%
    filter(start_date >= this_start_date & start_date <= this_end_date)
  
  l0gdd = gdd(tmax = wth_v$TX,
              tmin = wth_v$TN,
              tbase = 10,
              tbase_max = 30)
  
  weather_info <- x %>%
    mutate(
      l0gdd = tail(l0gdd, 1) 
    )
  
  # lag 7 days
  wth_v <- meteo_data %>%
    filter(start_date >= this_start_date - days(7) & start_date <= this_start_date) 
  
  l7gdd = gdd(tmax = wth_v$TX,
              tmin = wth_v$TN,
              tbase = 10,
              tbase_max = 30)
  
  weather_info <- weather_info %>%
    mutate(
      l7gdd = tail(l7gdd, 1)
    )
  
  # lag 14 days
  wth_v <- meteo_data %>%
    filter(start_date >= this_start_date - days(14) & start_date <= this_start_date) 
  
  l14gdd = gdd(tmax = wth_v$TX,
               tmin = wth_v$TN,
               tbase = 10,
               tbase_max = 30)
  weather_info <- weather_info %>%
    mutate(
      l14gdd = tail(l14gdd, 1)
    )
  
  # lag 21 days
  wth_v <- meteo_data %>%
    filter(start_date >= this_start_date - days(21) & start_date <= this_start_date) 
  
  l21gdd = gdd(tmax = wth_v$TX,
               tmin = wth_v$TN,
               tbase = 10,
               tbase_max = 30)
  weather_info <- weather_info %>%
    mutate(
      l21gdd = tail(l21gdd, 1)
    )
}, mc.cores = 4
))

# Land covers in 250m buffer
buffer_250 <- read_csv("Data/Land_covers/landcovers_buffer250.csv")
buffer_250 <- buffer_250 %>%
  group_by(trap_name, DN) %>%
  summarise(area = sum(area)) %>%
  pivot_wider(names_from = "DN", values_from = "area", values_fill = 0) %>%
  dplyr::select(trap_name, `2`, `4`, `5`,`6`,`7`,`9`,`10`,`11`,`14`,`15`,`16`,`17`,`19`,
                `20`,`21`,`22`,`23`,`24`)

buffer_250$urban_area <- buffer_250$`4` + buffer_250$`5` + buffer_250$`6` + buffer_250$`7`
buffer_250$wetland <- buffer_250$`2` + buffer_250$`10`
buffer_250$herbaceous_crops <- buffer_250$`19` + buffer_250$`20`
buffer_250$rice_field <- buffer_250$`21`

buffer_250 <- buffer_250 %>%
  dplyr::select(urban_area, wetland, herbaceous_crops, rice_field,
                `9`,`11`,`14`,`15`,`16`,`17`,`22`,`23`,`24`) %>%
  pivot_longer(!c(trap_name), names_to = "land_cover", values_to = "area") %>%
  arrange(trap_name, desc(area))

buffer_250 <- buffer_250 %>%
  group_by(trap_name) %>%
  mutate(has_rice_field = any(land_cover == "rice_field" & area > 0)) %>%
  filter(ifelse(has_rice_field, land_cover == "rice_field" & area > 0, area == max(area))) %>%
  ungroup() %>%
  dplyr::select(-has_rice_field, -area) 

traps <- read_excel("Data/Mosquito_data/traps.xlsx")
traps <- merge(traps, buffer_250, by = "trap_name") %>%
  dplyr::select(trap_name, trap_type, land_cover)
  
# Add trap_type and land_use variables
culex_df <- merge(culex_df, traps, by = "trap_name")
names(culex_df)[names(culex_df) == "land_cover"] <- "land_use"

# Save dataframe
write_csv(culex_df, "Data/Mosquito_data/Culex/culex_df.csv")
writexl::write_xlsx(culex_df, "Data/Mosquito_data/Culex/culex_df.xlsx")

################################################################################
#MODELS
################################################################################
library(readxl)
library(tidyverse)
library(lubridate)
library(dplyr)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(performance)
library(car)
library(MuMIn)
library(corrplot)
library(cmdstanr)
library(brms)
library(rstanarm)
library(loo)
library(tidybayes)
library(RColorBrewer)

# Delete all variables
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Read excel
culex_df <- read_excel("Data/Mosquito_data/Culex/culex_df.xlsx")

# Convert variables
culex_df$trap_name <- as.factor(culex_df$trap_name)
culex_df$AMT <- as.factor(culex_df$AMT)
culex_df$year <- factor(culex_df$year)
culex_df$month <- factor(culex_df$month, ordered = TRUE, levels = as.character(1:12))
culex_df$week <- factor(culex_df$week, ordered = TRUE, levels = as.character(1:53))
culex_df$trap_type <- as.factor(culex_df$trap_type)
culex_df$land_use <- as.factor(culex_df$land_use)

# Boxplots ---------------------------------------------------------------------
# Land use boxplot
boxplot_landuse <- ggplot(culex_df, aes(land_use, females)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Land use") +
  ylab("Mosquito counts") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "plain"))

# Trap name boxplot
boxplot_trap_name <- ggplot(culex_df, aes(trap_name, females)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Trap name") +
  ylab("Mosquito counts") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "plain"))

# Trap type boxplot
boxplot_trap_type <- ggplot(culex_df, aes(trap_type, females)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Trap type") +
  ylab("Mosquito counts") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "plain"))

# Generalized linear models ----------------------------------------------------
# Function collinearity
max.r <- function(x){
  corm <- cov2cor(vcov(x))
  corm <- as.matrix(corm)
  if (length(corm)==1){
    corm <- 0
    max(abs(corm))
  } else if (length(corm)==4){
    cormf <- corm[2:nrow(corm),2:ncol(corm)]
    cormf <- 0
    max(abs(cormf))
  } else {
    cormf <- corm[2:nrow(corm),2:ncol(corm)]
    diag(cormf) <- 0
    max(abs(cormf))
  }
}

# Full model
full_model <- glmer.nb(females ~ scale(TM) + scale(TM7) + scale(TM14) + scale(TM21) + 
                         scale(TX) + scale(TX7) + scale(TX14) + scale(TX21) + 
                         scale(TN) + scale(TN7) + scale(TN14) + scale(TN21) + 
                         scale(HRM) + scale(HRM7) + scale(HRM14) + scale(HRM21) + 
                         scale(PPT) + scale(PPT7) + scale(PPT14) + scale(PPT21) +
                         scale(l0gdd) + scale(l7gdd) + scale(l14gdd) + scale(l21gdd) +
                         log(trapping_effort) +
                         (1|year:trap_name) +
                         (1|trap_type) +
                         (1|land_use),
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                       data = culex_df)

summary(full_model)
car::vif(full_model) #vif <5 means low collinearity

options(na.action="na.fail")

full_model_dredge <- dredge(full_model, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 3), extra= c(max.r))
NCM <- get.models(full_model_dredge, subset = max.r<=0.6) # Retrieve non-collinear models (max.r <=0.6)
NCMDF <- model.sel(NCM) #Ranking and model selection
model.avg(NCMDF) #average
sw(NCMDF) #variable weights and models containing that variable. >0.8 vs 0.5

# test correlation
V <- culex_df %>%
  dplyr::select(TM, HRM7, PPT21)
M <-cor(V)
res1 <- cor.mtest(V, conf.level = .95)

corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"), p.mat = res1$p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05))

rm(V, M, res1)

# Save workspace
save.image("Models/Culex/glm_models/full_model.RData")

# Final model
culex_model <- glmer.nb(females ~ scale(TM) + scale(HRM7) + scale(PPT21) +
                          (1|year:trap_name) +
                          (1|trap_type) +
                          (1|land_use),
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                        data = culex_df)
save(culex_model, file= "Models/Culex/glm_models/culex_model.RData")
# load("Models/Culex/glm_models/culex_model.RData")

summary(culex_model)
car::vif(culex_model)

# Model and residuals check ----------------------------------------------------
simulationOutput <- simulateResiduals(fittedModel = culex_model, plot = F)
plot(simulationOutput)

testUniformity(simulationOutput)
testOutliers(simulationOutput)
testDispersion(simulationOutput) #no over dispersion
testQuantiles(simulationOutput)
testCategorical(simulationOutput, catPred = culex_df$land_use)
testZeroInflation(simulationOutput) #ratio >1 and p-value <0,05 means zero-inflation

# Zero inflation tests
countOnes <- function(x) sum(x == 0)  # testing for number of 0s
testGeneric(simulationOutput, summary = countOnes, alternative = "greater") # 0-inflation

fittedModel <- glmmTMB(females ~ scale(TM) + scale(HRM7) + scale(PPT21) +
                         (1|year:trap_name) +
                         (1|trap_type) +
                         (1|land_use),
                       ziformula = ~1,
                       family = nbinom1(link = "log"),
                       data = culex_df)

summary(fittedModel)

# Performance
check_model(culex_model)

# Bayesian models with BRMS ----------------------------------------------------
# zinb model
culex_model_zinb <- brm(
  bf(formula = females ~ scale(TM) + scale(HRM7) + scale(PPT21) +
       (1|year:trap_name) +
       (1|land_use)
  ),
  data = culex_df,
  family = zero_inflated_negbinomial(link = "log"),
  prior = set_prior("cauchy(0,2.5)", class="b"), 
  iter = 7000,
  chains = 4, 
  cores = 4,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.999))

save(culex_model_zinb, file= "Models/Culex/brms_models/culex_model_zinb.RData")
# load("Models/Culex/brms_models/culex_model_zinb.RData")

summary(culex_model_zinb)
density_plots <- plot(culex_model_zinb)
marginal_effects_plots <- marginal_effects(culex_model_zinb)
pp_check(culex_model_zinb)

# Save density plots
for (i in 1:length(density_plots)) {
  # Create the plot
  file_name <- density_plots[[i]]

  # Create a unique name for the file
  density_plot <- paste0("density_plots", i, ".png")

  # Save the plot
  ggsave(density_plot, file_name, width = 16, height = 8, dpi = 300, type = "jpg", limitsize = FALSE)
}

# Model checking (leave one out cross-validation)
loo <- loo(culex_model_zinb)
loo

problem_obs <- pareto_k_ids(loo, threshold = 1)
problematic_rows <- culex_df[problem_obs, ]
rm(problematic_rows)

# New model removing one problematic row
culex_df_cleaned <- culex_df[-problem_obs, ]

culex_model_zinb_2 <- brm(
  bf(formula = females ~ scale(TM) + scale(HRM7) + scale(PPT21) +
       (1|year:trap_name) +
       (1|land_use)
  ),
  data = culex_df_cleaned,
  family = zero_inflated_negbinomial(link = "log"),
  prior = set_prior("cauchy(0,2.5)", class="b"), 
  iter = 7000,
  chains = 4, 
  cores = 4,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.999))

save(culex_model_zinb_2, file = "Models/Culex/brms_models/culex_model_zinb_2.RData")
load("Models/Culex/brms_models/culex_model_zinb_2.RData")

# Model checking (leave one out cross-validation)
loo <- loo(culex_model_zinb)
loo

# Variability in mosquito counts across land use categories
# Extract the random effects from the model
random_effects_df <- as.data.frame(ranef(culex_model_zinb)$land_use)
random_effects_df$land_use <- rownames(random_effects_df)

# Land use variability plot
ggplot(random_effects_df, aes(x = land_use, y = Estimate.Intercept, 
                              ymin = Q2.5.Intercept, ymax = Q97.5.Intercept)) +
  geom_point(size = 4, color = "blue") +  # Points for the mean estimate
  geom_errorbar(width = 0.2, color = "black") +  # Error bars for confidence intervals
  labs(x = "Land Use", y = "Random Effect Estimate", 
       title = "Variability in Mosquito Counts Across Land Use Categories") +
  scale_x_discrete(labels = c("herbaceous_crops" = "Herbaceous crops", 
                              "rice_field" = "Rice fields", 
                              "urban_area" = "Urban area",
                              "wetland" = "Wetland")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate labels for readability
    axis.text.y = element_text(size = 12)
  )

################################################################################
#PREDICT
################################################################################
library(readxl)
library(tidyverse)
library(lubridate)
library(dplyr)
library(ggplot2)
library(purrr)
library(zoo)
library(parallelly)
library(parallel)
library(janitor)
library(pollen)
library(sf)
library(sp)
library(gstat)
library(raster)
library(RColorBrewer)
library(viridis)
library(cmdstanr)
library(brms)
library(rstanarm)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Dataframe preparation --------------------------------------------------------
# Add study area: generate grid and points
cat_st <- st_read(paste0("Data/Land_covers/study_area_shp.shp")) %>%
  st_transform(32631) 

cell_res <- 150 #distance between points (grid dimension)

new_points <- st_make_grid(cat_st, cellsize = c(cell_res, cell_res), # 150 meters (UTM)
                           what = "centers", square = TRUE) %>%
  st_sf() %>% 
  st_join(cat_st, join = st_intersects, left=FALSE) %>%
  dplyr::select(UTM1X1)%>%
  mutate(
    long = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2],
    trap_name = c(paste0("new_trap_", 1:nrow(.))))

# Erasing hypothetical traps in the sea
erase_traps <- read_excel("Data/Mosquito_data/erase_traps.xlsx")
erase_traps <- rename(erase_traps, trap_name = trap_nm)
new_points <- anti_join(new_points, erase_traps, by = "trap_name")

# Add land
url <- "Data/Land_covers/UsosCobertes_2017/MUCSC_2017_30_m_v_4.tif"

us <- raster(url)

new_points <- new_points %>%
  mutate(
    land_use = raster::extract(us, st_coordinates(new_points))) %>%
  drop_na(land_use) # extracting land uses

new_points <- new_points %>%
  mutate(
    land_use = ifelse(land_use %in% c(4, 5, 6, 7) , "urban_area",
                      ifelse(land_use %in% c(19, 20), "herbaceous_crops",
                             ifelse(land_use %in% c(2, 10), "wetland",
                                    ifelse(land_use == 21, "rice_field", "others"))))) %>%
  dplyr::select(UTM1X1, trap_name, long, lat, land_use) %>%
  mutate(trap_type = "new_trap_type")

st_geometry(new_points) <- NULL

# Add weather variables
# Loading data from Meteocat. The files are saved in local.
w1 <- read.csv("Data/Meteo_data/Meteo_data_W1.csv", sep = ";") %>%
  dplyr::select(DATA, TM, TX, TN, HRM, PPT) %>%
  mutate(
    DATA = as.POSIXct(DATA, tz = "GMT", format = "%d/%m/%Y"))

u2 <- read.csv("Data/Meteo_data/Meteo_data_U2.csv", sep = ";") %>%
  dplyr::select(DATA, TM, TX, TN, HRM, PPT) %>%
  mutate(
    DATA = as.POSIXct(DATA, tz = "GMT", format = "%d/%m/%Y"))

meteo_data <- w1 %>%
  left_join(u2, by = "DATA") %>%
  mutate(
    TM = ifelse(is.na(TM.x), TM.y, TM.x),
    TX = ifelse(is.na(TX.x), TX.y, TX.x),
    TN = ifelse(is.na(TN.x), TN.y, TN.x),
    HRM = ifelse(is.na(HRM.x), HRM.y, HRM.x), 
    PPT = ifelse(is.na(PPT.x), PPT.y, PPT.x), ) %>%
  dplyr::select(DATA, TM, TX, TN, HRM, PPT) 

colnames(meteo_data)[colnames(meteo_data) == "DATA"] ="start_date"

# Mosquito weather index
meteo_data = meteo_data %>% mutate(
  FHme = case_when(HRM < 40~0, HRM  >95~0,
                   (HRM  >=40 & HRM <= 95)~
                     ((HRM/55)-(40/55)) ),
  FTme = case_when(TM<=15~0, TM>30~0, 
                   (TM>15 & TM <=20)~
                     (.2*TM)-3,
                   (TM>20 & TM<=25)~1,
                   (TM>25 & TM <= 30)~
                     (-.2*TM)+6),
  mwime = FHme*FTme)

# 7, 14 and 21-day cumulative weather variables
meteo_data <- meteo_data %>% 
  mutate(
    start_date = as_date(start_date),
    TM7 = rollmean(TM, k=7, align = "right", fill = NA),
    TM14 = rollmean(TM, k=14, align = "right", fill = NA),
    TM21 = rollmean(TM, k=21, align = "right", fill = NA),
    TX7 = rollmean(TX, k=7, align = "right", fill = NA),
    TX14 = rollmean(TX, k=14, align = "right", fill = NA),
    TX21 = rollmean(TX, k=21, align = "right", fill = NA),
    TN7 = rollmean(TN, k=7, align = "right", fill = NA),
    TN14 = rollmean(TN, k=14, align = "right", fill = NA),
    TN21 = rollmean(TN, k=21, align = "right", fill = NA),
    HRM7 = rollmean(HRM, k=7, align = "right", fill = NA),
    HRM14 = rollmean(HRM, k=14, align = "right", fill = NA),
    HRM21 = rollmean(HRM, k=21, align = "right", fill = NA),
    PPT7 = rollapplyr(PPT, 7, FUN=sum, align = "right", fill = NA),
    PPT14 = rollapplyr(PPT, 14, FUN=sum, align = "right", fill = NA),
    PPT21 = rollapplyr(PPT, 21, FUN=sum, align = "right", fill = NA), 
    mwi7 = rollmean(mwime, k = 7, align = "right", fill = NA),
    mwi14 = rollmean(mwime, k = 14, align = "right", fill = NA),
    mwi21 = rollmean(mwime, k=21, align = "right", fill = NA))

# Loop for averaging weather variables
new_wth <- function(all_dates){
  
  wth_data <-data.frame()
  
  for(date in all_dates){ 
    
    x0 <- as_date(date) - 1
    x1 <- as_date(date)
    
    int <- seq(x0, x1, by= "day")
    
    TM <- mean(meteo_data$TM[meteo_data$start_date %in% int]) 
    TX <- mean(meteo_data$TX[meteo_data$start_date %in% int])
    TN <- mean(meteo_data$TN[meteo_data$start_date %in% int])
    HRM <- mean(meteo_data$HRM[meteo_data$start_date %in% int])
    PPT <- sum(meteo_data$PPT[meteo_data$start_date %in% int])
    mwime <- mean(meteo_data$mwime[meteo_data$start_date %in% int])
    TM7 <- meteo_data$TM7[meteo_data$start_date %in% x0]
    TM14 <- meteo_data$TM14[meteo_data$start_date %in% x0]
    TM21 <- meteo_data$TM21[meteo_data$start_date %in% x0]
    TX7 <- meteo_data$TX7[meteo_data$start_date %in% x0]
    TX14 <- meteo_data$TX14[meteo_data$start_date %in% x0]
    TX21 <- meteo_data$TX21[meteo_data$start_date %in% x0]
    TN7 <- meteo_data$TN7[meteo_data$start_date %in% x0]
    TN14 <- meteo_data$TN14[meteo_data$start_date %in% x0]
    TN21 <- meteo_data$TN21[meteo_data$start_date %in% x0]
    HRM7 <- meteo_data$HRM7[meteo_data$start_date %in% x0]
    HRM14 <- meteo_data$HRM14[meteo_data$start_date %in% x0]
    HRM21 <- meteo_data$HRM21[meteo_data$start_date %in% x0]
    PPT7 <- meteo_data$PPT7[meteo_data$start_date %in% x0]
    PPT14 <- meteo_data$PPT14[meteo_data$start_date %in% x0]
    PPT21 <-meteo_data$PPT21[meteo_data$start_date %in% x0]
    mwi7 <- meteo_data$mwi7[meteo_data$start_date %in% x0]
    mwi14 <- meteo_data$mwi14[meteo_data$start_date %in% x0]
    mwi21 <- meteo_data$mwi21[meteo_data$start_date %in% x0]
    
    wth <- data.frame(
      start_date = x1,
      TM = TM,
      TX = TX,
      TN = TN,
      HRM = HRM,
      PPT = PPT,
      mwime = mwime,
      TM7 = TM7,
      TM14 = TM14,
      TM21 = TM21,
      TX7 = TX7,
      TX14 = TX14,
      TX21 = TX21,
      TN7 = TN7,
      TN14 = TN14,
      TN21 = TN21,
      HRM7 = HRM7,
      HRM14 = HRM14,
      HRM21 = HRM21,
      PPT7 = PPT7,
      PPT14 = PPT14,
      PPT21 = PPT21,
      mwi7 = mwi7,
      mwi14 = mwi14,
      mwi21 = mwi21)
    
    wth_data <- rbind(wth_data, wth)} 
  
  return(wth_data)}

all_dates <- seq.Date(from = as_date("2015-01-01"), to = as_date("2018-12-31"), by = "day")

wth_pred <- new_wth(all_dates)

# Growing degree days
wth_pred <-  bind_rows(mclapply(1:nrow(as.data.frame(wth_pred)), function(i){
  # Verbose
  if (i %% 50 == 0) {print(i)}
  
  x = wth_pred[i,]
  this_start_date = x$start_date
  this_end_date = x$start_date
  
  # Exposure week
  wth_v <- meteo_data %>%
    filter(start_date >= this_start_date & start_date <= this_end_date)
  
  l0gdd = gdd(tmax = wth_v$TX,
              tmin = wth_v$TN,
              tbase = 10,
              tbase_max = 30)
  
  weather_info <- x %>%
    mutate(
      l0gdd = tail(l0gdd, 1) 
    )
  
  # lag 7 days
  wth_v <- meteo_data %>%
    filter(start_date >= this_start_date - days(7) & start_date <= this_start_date) 
  
  l7gdd = gdd(tmax = wth_v$TX,
              tmin = wth_v$TN,
              tbase = 10,
              tbase_max = 30)
  
  weather_info <- weather_info %>%
    mutate(
      l7gdd = tail(l7gdd, 1)
    )
  
  # lag 14 days
  wth_v <- meteo_data %>%
    filter(start_date >= this_start_date - days(14) & start_date <= this_start_date) 
  
  l14gdd = gdd(tmax = wth_v$TX,
               tmin = wth_v$TN,
               tbase = 10,
               tbase_max = 30)
  weather_info <- weather_info %>%
    mutate(
      l14gdd = tail(l14gdd, 1)
    )
  
  # lag 21 days
  wth_v <- meteo_data %>%
    filter(start_date >= this_start_date - days(21) & start_date <= this_start_date) 
  
  l21gdd = gdd(tmax = wth_v$TX,
               tmin = wth_v$TN,
               tbase = 10,
               tbase_max = 30)
  weather_info <- weather_info %>%
    mutate(
      l21gdd = tail(l21gdd, 1)
    )
}, mc.cores = 4
))

# Merge both databases
newdata <- merge(new_points, wth_pred)

# Remove data frames
rm(erase_traps,cat_st,cell_res,url,us,w1,u2,meteo_data,new_wth,all_dates,new_points,wth_pred)

# Add year, month and week columns
newdata <- newdata %>% mutate(year=year(start_date))
newdata$year <- as.factor(newdata$year)
newdata <- newdata %>% mutate(month=month(start_date))
newdata$month <- as.factor(newdata$month)
newdata <- newdata %>% mutate(week=week(start_date))
newdata$week <- as.factor(newdata$week)

# Add trapping effort
newdata$trapping_effort <- 1

# Save workspace
save.image("Prediction/Culex/prediction_prep.RData")

# Predictions ------------------------------------------------------------------
# Delete all variables.
rm(list = ls())

# Load prediction dataframe
load("Prediction/Culex/prediction_prep.RData")

# Convert variables
newdata$UTM1X1 <- as.factor(newdata$UTM1X1)
newdata$trap_name <- as.factor(newdata$trap_name)
newdata$land_use <- as.factor(newdata$land_use)
newdata$trap_type <- as.factor(newdata$trap_type)

# Load model
load("Models/Culex/brms_models/culex_model_zinb_2.RData")
load("Models/Culex/glm_models/culex_model.RData")

# # Prediction using glm model
# pred <- newdata %>%
#   dplyr::mutate(
#     pp = predict(culex_model, newdata = newdata, type = c("response"), allow.new.levels = TRUE))

# Prediction using Bayesian model
nrow_these_pred_points = nrow(newdata)
max_chunksize = 300000
ncores = 1
chunksize = min(as.integer((nrow_these_pred_points/ncores)), max_chunksize)

pred <- bind_rows(mclapply(seq(1, nrow(newdata), chunksize), function(i){
  print(i)

  data_chunk = newdata[i:min(nrow(newdata), (i+(chunksize-1))), ]
  flush.console()

  pp <- apply(posterior_predict(culex_model_zinb,
                                newdata = data_chunk,
                                allow_new_levels = TRUE,
                                re_formula = NULL,
                                ndraws = 1000),
              2, function(x) mean(x)) %>%
    as.data.frame()

  data_chunk <- data_chunk %>% dplyr::select(UTM1X1, long, lat,
                                             trap_name, land_use,
                                             start_date,
                                             year, month, week)

  pp <- bind_cols(data_chunk, pp)

  return(pp)

}, mc.cores = ncores))

colnames(pred)[colnames(pred) == "."] ="pp"

write.csv(pred, "Prediction/Culex/pred.csv")

# Mapping the predictions ------------------------------------------------------
# Delete all variables.
rm(list = ls())

# Read prediction
pred <- read.csv("Prediction/Culex/pred_2.csv")

# Create shapefile
newdata_sf <- pred %>%
  dplyr::select(trap_name, long, lat, pp) %>%
  group_by(trap_name, long, lat) %>%
  summarise(mean_counts = mean(pp)) %>%
  st_as_sf(., coords = c("long", "lat"), crs = 32631, remove = FALSE)

st_write(newdata_sf, "Prediction/Culex/culex_spatial_prediction.shp", append=FALSE) #save shapefile

# Average map
map_0 <- ggplot() +
  geom_tile(data = newdata_sf, aes(x = long, y = lat, fill = mean_counts), alpha = 0.8) +
  viridis::scale_fill_viridis() +
  theme_bw() +
  labs(title = "Average Prediction (2015-2018)", fill = "Mean counts", x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(size = 16, face = "bold"), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))

map_0
filename <- "Plots/Culex/culex_map_0.jpg"
ggsave(filename, dpi = 300, width = 7, height = 6, units = "in", type = "jpg", quality = 100)

# Temporal variability plot ----------------------------------------------------
agg_newdata <- pred %>%
  dplyr::select(trap_name, start_date, year, week, pp) %>%
  group_by(start_date, year, week) %>%
  summarise(pp = mean(pp)) %>%
  ungroup() %>%
  group_by(year, week) %>%
  summarise(pp = sum(pp)) # traps average + weekly sum

agg_newdata$year <- as.factor(agg_newdata$year)
agg_newdata$week <- as.integer(as.character(agg_newdata$week))

mean_counts <- agg_newdata %>%
  dplyr::select(week, pp) %>% 
  group_by(week) %>%
  summarise(pp = mean(pp))

agg_newdata <- rbind(agg_newdata, mean_counts %>% mutate(year = "Mean"))

# Define the color palette for the plot
my_palette <- c("red", brewer.pal(n = 4, name = "Set2"))

# Plot
temporal_variability_plot <- ggplot(data = agg_newdata, aes(x = week, y = pp, color = year, group = year)) +
  geom_point(size = 0.3) +  # Small points for visibility
  geom_line(size = 1, aes(linetype = year)) +  # Thicker lines for clarity
  scale_color_manual(values = my_palette, name = "Year", 
                     limits = c("Mean", "2015", "2016", "2017", "2018")) + 
  scale_linetype_manual(values = c("solid", rep("dashed", 4)), name = "Year", 
                        limits = c("Mean", "2015", "2016", "2017", "2018")) +
  labs(x = "Week", y = "Average predicted mosquito counts per trap", color = "Year", linetype = "Year") +
  scale_x_continuous(breaks = seq(1, max(agg_newdata$week), by = 5)) +  # Set breaks for the X-axis
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400)) + # Set breaks for the Y-axis
  theme_minimal(base_size = 14) + 
  theme(
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    axis.text.x = element_text(size = 14, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black", size = 0.5),  # Add black axis lines for better framing
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Light gray for major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = c(0.96, 0.5),  # Position the legend in the top-right corner
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.7, "cm"),
    legend.spacing.y = unit(0.2, "cm"),  # Reduce vertical spacing between legend items
    legend.background = element_rect(fill = "white")  # Set a white background for the legend
  )

temporal_variability_plot
filename <- "Plots/Culex/culex_temporal_variability.jpg"
ggsave(filename, dpi = 300, width = 15, height = 7, units = "in", type = "jpg", quality = 100)

rm(mean_counts)

# Raw data and prediction plot -------------------------------------------------
# Read Culex data
culex_df <- read_excel("Data/Mosquito_data/Culex/culex_df.xlsx")

# Convert variables
culex_df$trap_name <- as.factor(culex_df$trap_name)
culex_df$AMT <- as.factor(culex_df$AMT)
culex_df$year <- factor(culex_df$year)
culex_df$month <- factor(culex_df$month, ordered = TRUE, levels = as.character(1:12))
culex_df$week <- factor(culex_df$week, ordered = TRUE, levels = as.character(1:53))
culex_df$trap_type <- as.factor(culex_df$trap_type)
culex_df$land_use <- as.factor(culex_df$land_use)

# Calculate the mean number of female mosquitoes per week across years
culex_df_mean <- culex_df %>%
  filter(year %in% c("2015", "2016", "2017", "2018")) %>%
  group_by(week, year) %>%
  summarise(females = mean(females)) %>%
  ungroup() %>%
  group_by(week) %>%
  summarise(females = mean(females))

culex_df_mean$week <- as.integer(culex_df_mean$week)

# Filter predicted data (average year)
agg_newdata_mean <- pred %>%
  group_by(week, year) %>%
  summarise(pp = mean(pp)) %>%
  ungroup() %>%
  group_by(week) %>%
  summarise(pp = mean(pp))

# Define a scaling factor to make real and predicted values comparable
scale_factor <- max(agg_newdata_mean$pp) / max(culex_df_mean$females)

# Define manual breaks for the y-axis
y_breaks_pred <- c(0, 16, 32, 48)
y_breaks_real <- c(0, 7, 14, 21)  
y_breaks_scaled <- y_breaks_real * scale_factor

# Plot
observed_predicted_plot <- ggplot(data = agg_newdata_mean, aes(x = week)) +
  # Scaled observed data
  geom_line(data = culex_df_mean, aes(x = week, y = females * scale_factor, color = "Observed"), 
            size = 1) +
  # Predicted values
  geom_line(aes(y = pp, color = "Predicted"), size = 1) +  
  # Color and legend settings
  scale_color_manual(name = "Mosquito counts", 
                     values = c("Predicted" = "red", "Observed" = "black")) +
  # Axis settings with specific breaks
  scale_x_continuous(breaks = seq(1, max(agg_newdata_mean$week), by = 5)) +
  scale_y_continuous(
    name = "Average Predicted Mosquito Counts per Trap",  
    breaks = y_breaks_pred,
    sec.axis = sec_axis(~ . / scale_factor, 
                        name = "Average Observed Female Mosquito Counts",
                        breaks = y_breaks_real)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.grid.major.x = element_line(color = "gray90", size = 0.5),
    panel.grid.major.y = element_line(color = "gray90", size = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = c(0.9, 0.8),  # Adjust legend position
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.7, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.background = element_rect(fill = "white")
  )

observed_predicted_plot

# Abundance by UTM1x1 ----------------------------------------------------------
culex_UTM1x1 <- pred %>%
  dplyr::select(UTM1X1, trap_name, long, lat, pp) %>%
  group_by(UTM1X1) %>%
  summarise(counts = mean(pp)) 

write_csv(culex_UTM1x1, "Prediction/Culex/culex_UTM1x1.csv")
