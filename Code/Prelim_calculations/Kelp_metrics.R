## Barkley Sound kelp forest data processing
## Author: Claire Attridge
## Origin date: August 2022
## Updated by Em Lim Dec 2024

# Loading base packages
library(tidyverse)

#### Kelp density, height, & biomass cleaning ----

####### Reading in our kelp density data
dens <- read_csv("Data/Team_kelp/kelp_density_2022.csv") %>%
  mutate(Macro = (Macro_5m2/5), Nereo=(Nereo_5m2/5)) %>% # Changing units to /m2 area
  rowwise() %>% # To sum across rows
  mutate(Kelp = sum(Macro,Nereo)) %>% # Sum macro and nereo to get total kelp dens / transect
  ungroup() %>%# Stop rowwise 
  filter(!str_detect(SiteName, 'REDO')) # Removing the site 'redos' for current purposes

####### Reading in our kelp height and biomass data
kelp <- read_csv("Data/Team_kelp/kelp_morphology_2022.csv") %>%
  as.data.frame() %>%
  filter(!str_detect(SiteName, 'REDO')) %>% # Removing the site 'redos' for current purposes
  rowwise() %>%
  mutate(Sub_diam_cm = (Sub_circ_cm/pi)) %>% # Converting sub-bulb circumference to diameter
  ungroup()

kelp2 <- kelp[-266,] # Removing the outlier point at Second Beach South!


# Equation for converting from sub-bulb diamater (cm) to biomass (g) as per our Nereo sub-bulb model (using equation for all 3 sample sites combined)
# See 'nereo_biomass.R' for the origin relationship
formula <- function(x){
  (150.7966*(x)^2 -216.2721*(x) + 315.0124)
}

# options(scipen=999) # Turning off scientific notation

# Statement for applying sub-bulb equation when biomass is absent (i.e., when it needs to happen)
kelp2 <- kelp %>%
  mutate(Biomass_g = ifelse(Species == "Nereo", formula(Sub_diam_cm), Biomass_g))


# Averaging to transect from individual kelp samples (i.e., averages of the 5 sampled individuals / transect for each species)
kelptrans <- kelp2 %>%
  mutate(SiteName = as.factor(SiteName), 
         Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName, Transect, Species) %>% # Averaging to transect and species
  summarise(HeightT = mean(Height_m, na.rm=TRUE),
            BiomassTind_g = mean(Biomass_g, na.rm=TRUE)) %>% # Ave individual biomass (g) / transect
  ungroup()

# replace nan with 0
kelptrans$BiomassTind_g[is.nan(kelptrans$BiomassTind_g)] <- 0
kelptrans$HeightT[is.nan(kelptrans$HeightT)] <- 0

# make df with column for macro biomass and column for nereo biomass
nereo <- kelptrans %>%
  filter(Species == "Nereo") %>%
  transmute(SiteName = SiteName,
            Transect = Transect,
            nereo_biomass_ind = BiomassTind_g,
            nereo_height_ind = HeightT)

macro <- kelptrans %>%
  filter(Species == "Macro") %>%
  transmute(SiteName = SiteName,
            Transect = Transect,
            macro_biomass_ind = BiomassTind_g,
            macro_height_ind = HeightT)

nereo_macro <- merge(nereo, macro, by = c("SiteName", "Transect"), all=TRUE)


# Bringing in the kelp density (transect level) data
kelpjoin <- merge(nereo_macro, dens, by = c("SiteName", "Transect"), all=TRUE)

# Getting to the average biomass / m2 area for each transect
kelptog <- kelpjoin %>%
  rowwise() %>% # To sum across rows
  mutate(macro_biomass_trans = Macro*macro_biomass_ind/1000, # multiply avg indiv kelp biomass by transect kelp density/m2
         nereo_biomass_trans = Nereo*nereo_biomass_ind/1000,
         biomass_trans_mean = sum(macro_biomass_trans, nereo_biomass_trans, na.rm=TRUE), # add nereo and macro density
         height_mean = mean(c(macro_height_ind, nereo_height_ind), na.rm=TRUE)) %>% # could do height better I think
  ungroup() # Stop rowwise 


# make df to use in kelp_pee_analysis
kelp_metrics <- read_csv("Data/Team_kelp/Output_data/site_names.csv") %>% # site codes and names
  left_join(read_csv("Data/Team_kelp/Output_data/forest_area.csv"), by = "SiteName") %>% # forest area
  left_join(kelptog, by = "SiteName") %>% # now join in biomass, density, height
  select(SiteName, site_code, Transect, Date, Time_start, Depth_m, RLS_dist, Macro_5m2, Nereo_5m2, Macro, Nereo, Kelp, macro_biomass_trans, nereo_biomass_trans, biomass_trans_mean, height_mean, Area_m2) %>%
  ungroup() # just to be safe

kelp_metrics$height_mean[is.nan(kelp_metrics$height_mean)] <- 0

# write_csv(kelp_metrics, "Data/Team_kelp/Output_data/kelp_metrics2024.csv")

  


# Average to site level right away -----
# I could also go straight to site level averages but I don't think we want to
# Averaging to site from individual kelp samples (i.e., averages of the 5 sampled individuals / transect)
kelptrans3 <- kelp2 %>%
  mutate(SiteName = as.factor(SiteName), 
         Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName, Species) %>% # Averaging to transect 
  summarise(HeightT = mean(Height_m, na.rm=TRUE),
            BiomassTind_g = mean(Biomass_g, na.rm=TRUE)) %>% # Ave individual biomass (g) / transect
  ungroup() %>%
  mutate(HeightT == ifelse(SiteName == "Less Dangerous Bay", 0, HeightT),
         BiomassTind_g == ifelse(SiteName == "Less Dangerous Bay", 0, BiomassTind_g))


# make df with column for macro biomass and column for nereo biomass
nereo2 <- kelptrans3 %>%
  filter(Species == "Nereo") %>%
  transmute(SiteName = SiteName,
            nereo_biomass_ind = BiomassTind_g,
            nereo_height_ind = HeightT)

macro2 <- kelptrans3 %>%
  filter(Species == "Macro") %>%
  transmute(SiteName = SiteName,
            macro_biomass_ind = BiomassTind_g,
            macro_height_ind = HeightT)

ind <- merge(nereo2, macro2, by = c("SiteName"), all=TRUE)

# Bringing in the kelp density (transect level) data my way
kelpjoin3 <- dens %>%
  group_by(SiteName) %>%
  summarise(macro_den_sum = sum(Macro_5m2),
            nereo_den_sum = sum(Nereo_5m2)) %>%
  left_join(ind, by = "SiteName") %>%
  rowwise() %>%
  mutate(total_macro_biomass = macro_den_sum*macro_biomass_ind,
         total_nereo_biomass = nereo_den_sum*nereo_biomass_ind,
         total_biomass_kg = sum(total_macro_biomass, total_nereo_biomass, na.rm=TRUE)/1000,
         Biomass_kieran = total_biomass_kg/20)

# ok now compare those estimates
em2 <- kelp_metrics %>%
  group_by(SiteName) %>%
  summarise(mean_bio_recalc = mean(biomass_trans_mean)) %>%
  mutate(mean_bio_recalc = round(mean_bio_recalc, digits = 5))

kieran <- kelpjoin3 %>%
  select(SiteName, Biomass_kieran)%>%
  mutate(Biomass_kieran = round(Biomass_kieran, digits = 5))

check2 <- em2 %>%
  left_join(kieran, by = "SiteName")



# Grouping/averaging from transect to site level old way -----
kelpgrp <- kelptog %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  group_by(SiteName, Depth_datum_m) %>% # Averaging to site & keeping depth (m)
  summarise(HeightM = weighted.mean(HeightT, na.rm=T), HeightSD = sd(HeightT, na.rm=T), # Ave height (m)
            BiomassM = mean(BiomassTkg, na.rm=T), BiomassSD = sd(BiomassTkg, na.rm=T), # Ave biomass (kg / m2)
            DensityM = mean(Kelp), DensitySD = sd(Kelp, na.rm=T), # Ave density any kelp (m2)
            MacroM = mean(Macro), MacroSD = sd(Macro, na.rm=T), # Ave density Macro (m2)
            NereoM = mean(Nereo), NereoSD = sd(Nereo, na.rm=T)) # Ave density Nereo (m2)


# Adding composition identity column for each site
kelpgrp <- kelpgrp %>%
  mutate(Composition = case_when(NereoM != 0 & MacroM != 0 ~ "Mixed",
                                 MacroM != 0 ~ "Macro",
                                 NereoM != 0 ~ "Nereo",
                                 TRUE ~ as.character("None")))
