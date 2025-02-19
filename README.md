# Arboviral disease risk in a Mediterranean wetland area

Here, we assess the potential risk of arboviral transmission in a high-risk, non-endemic Mediterranean region of northern Catalonia, analyzing mosquito vector populations (Aedes albopictus, Culex pipiens, Culex modestus and Culex theileri) and West Nile virus (WNV) avian hosts.

# Installation

Make sure to have R 4.3.2 (or later) and, optionally, RStudio installed for an enhanced development environment.

1. Clone the repository:
git clone https://github.com/username/project.git
cd project

2. Install required packages: Open R or RStudio and run the following code to install the required packages:

packages <- c("readxl", "writexl", "tidyverse", "lubridate", "dplyr", "ggplot2", "purrr", "zoo", "parallelly", "parallel", "janitor", "pollen", "lme4", "DHARMa", "glmmTMB", "performance", "car", "MuMIn", "corrplot", "cmdstanr", "brms", "rstanarm", "loo", "tidybayes", "RColorBrewer", "sf", "raster", "viridis", "ggrepel", "vegan", "factoextra", "tibble", "MASS", "svglite", "pracma", "ggthemes")

packages_to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(packages_to_install)) install.packages(packages_to_install) #Install packages that are not already installed

lapply(packages, library, character.only = TRUE) #Load packages

3. Run the project: Open the main script file (e.g., main.R) in R or RStudio and execute it:
source("main.R")

# Project Directory Structure

The project directory PNAE is organized as follows:

```markdown

PNAE/
├── Code/                  		    
│   │
│   ├── 1_Aedes_albopictus.R    	    # Dataframe preparation, models, and predictions for Aedes albopictus data
│   ├── 1_Culex.R               	    # Dataframe preparation, models, and predictions for Culex spp. data
│   ├── 1_Culex_species.R                   # Exploring Culex spp. (pipiens, modestus and theileri)
│   ├── 2_Birds_experiments_review.R        # Host competence calculations using experimental infection data
│   ├── 3_Birds.R               	    # WNV host spatial analysis
│   ├── 4_Risk.R                	    # Overlapping mosquitoes and birds
│   ├── 5_Humans.R              	    # Human population dynamics
│   ├── 6_Reanalysis.R    	            # Birds + risk reanalysis using an alternative infection threshold (4 PFU/ml)
│   │	
├── Data/
│   │ 
│   ├── Birds_data/       		
│   │   ├── Data/
│   │   │   ├── Abundance_data_atlas.csv                     	# Abundance data from the "Atlas of Nesting Birds of Catalonia"
│   │   │   ├── PNAE_birds.xlsx                             	# List of bird species
│   │   │   └── Review_experimental_infections_birds_PNAE.xlsx  # WNV experimental infections global review
│   │   └── Outputs/                                         	# Output folder
│   │	
│   ├── Human_data/
│   │   ├── Barcelona_human_cases.xlsx                       	# Imported cases in Barcelona
│   │   ├── Catalunya_human_cases.xlsx                       	# Imported cases in Catalunya
│   │   ├── Girona_human_cases.xlsx                          	# Imported cases in Girona
│   │   └── Població_estacional_Alt_Empordà.csv              	# FTE seasonal population (Alt Empordà region)
│   │	
│   ├── Land_covers/      		
│   │   ├── UsosCobertes_2017/                               	# Land cover map            
│   │   ├── landcovers_buffer250.csv                         	# Land covers in a 250 m radius buffer around each mosquito trap
│   │   ├── landcovers_utm1x1.csv                            	# Land covers in each 1x1 km grid cell
│   │   └── study_area_shp.shp                               	# Study area shapefile (to generate predictions)
│   │	
│   ├── Meteo_data/       		
│   │   ├── Meteo_data_U2.csv                       	     	# Secondary meteo station data
│   │   └── Meteo_data_W1.csv             		     	# Main meteo station data
│   │	
│   ├── Mosquito_data/    		
│   │   ├── Aedes_albopictus/                                	# Aedes albopictus output folder
│   │   ├── Culex/                                           	# Culex spp. output folder
│   │   └── Mosquito_Surveillance_SCM/                       	# 20-year mosquito surveillance database
│   │	
│   └── Risk_data/        		                     	# Risk assessment output folder
│
├── Models/               		
│   │
│   ├── Aedes_albopictus/
│   │   ├── brms_models/  		# Bayesian models for Aedes albopictus
│   │   └── glm_models/   		# GLM models for Aedes albopictus
│   └── Culex/
│       ├── brms_models/  		# Bayesian models for Culex spp.
│       └── glm_models/   		# GLM models for Culex spp.
│
├── Plots/                		# Plots and visualizations generated from data
│   │
│   ├── Aedes_albopictus/
│   ├── Birds/
│   └── Culex/
│
└── Predictions/          		# Contains space-time prediction files
    │
    ├── Aedes_albopictus/
    └── Culex/



