## Barkley Sound kelp forest data processing
## Author: Claire Attridge
## Origin date: August 2022

## Updated by Em Lim Dec 2024

# Loading base packages
library(tidyverse)
#library(ggpubr)
library(gridExtra)

#### Kelp density, height, & biomass cleaning ----

####### Reading in our kelp density data

dens <- read_csv("Data/Team_kelp/kelp_density_2022.csv") %>%
  mutate(Macro = (Macro_5m2/5), Nereo=(Nereo_5m2/5)) %>% # Changing units to /m2 area
  rowwise() %>% # To sum across rows
  mutate(Kelp = sum(Macro,Nereo)) %>% # Sum macro and nereo to get total kelp dens / transect
  ungroup() # Stop rowwise 

# Removing the site 'redos' for current purposes
dens2 <- dens %>%
  filter(!str_detect(SiteName, 'REDO'))

# Dataframe for site-specific density averages 
densgrp <- dens2 %>%
  group_by(SiteName, Depth_datum_m) %>% # Averaging to site
  summarise(KelpM = mean(Kelp), KelpSD = sd(Kelp))



####### Reading in our kelp height and biomass data

kelp <- read_csv("Data/Team_kelp/kelp_morphology_2022.csv") %>%
  as.data.frame()

kelp2 <- kelp[-266,] # Removing the outlier point at Second Beach South!

# Removing the site 'redos' for current purposes
kelp3 <- kelp2 %>%
  filter(!str_detect(SiteName, 'REDO'))


# Converting sub-bulb circumference to diameter
kelp4 <- kelp3 %>%
  rowwise() %>%
  mutate(Sub_diam_cm = (Sub_circ_cm/pi)) %>%
  ungroup()

# Equation for converting from sub-bulb diamater (cm) to biomass (g) as per our Nereo sub-bulb model (using equation for all 3 sample sites combined)
# See 'nereo_biomass.R' for the origin relationship
formula <- function(x){
  (150.7966*(x)^2 -216.2721*(x) + 315.0124)
}

# options(scipen=999) # Turning off scientific notation

# Statement for applying sub-bulb equation when biomass is absent (i.e., when it needs to happen)
kelp5 <- kelp4 %>%
  mutate(Biomass_g = ifelse(Species == "Nereo", formula(Sub_diam_cm), Biomass_g))


# Averaging to transect from individual kelp samples (i.e., averages of the 5 sampled individuals / transect)
kelptrans <- kelp5 %>%
  mutate(SiteName = as.factor(SiteName), 
         Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName, Transect) %>% # Averaging to transect 
  summarise(HeightT = mean(Height_m, na.rm=TRUE),
            BiomassTind_g = mean(Biomass_g, na.rm=TRUE)) %>% # Ave individual biomass (g) / transect
  ungroup()

# do it again
kelptrans2 <- kelp5 %>%
  mutate(SiteName = as.factor(SiteName), 
         Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName, Transect, Species) %>% # Averaging to transect 
  summarise(HeightT = mean(Height_m, na.rm=TRUE),
            BiomassTind_g = mean(Biomass_g, na.rm=TRUE)) %>% # Ave individual biomass (g) / transect
  ungroup()

# replace nan with 0
is.nan.data.frame <- function(x) {
  do.call(cbind, lapply(x, is.nan))
}

kelptrans$BiomassTind_g[is.nan(kelptrans$BiomassTind_g)] <- 0
kelptrans$HeightT[is.nan(kelptrans$HeightT)] <- 0

kelptrans2$BiomassTind_g[is.nan(kelptrans2$BiomassTind_g)] <- 0
kelptrans2$HeightT[is.nan(kelptrans2$HeightT)] <- 0

# make df with column for macro biomass and column for nereo biomass
nereo <- kelptrans2 %>%
  filter(Species == "Nereo") %>%
  transmute(SiteName = SiteName,
            Transect = Transect,
            nereo_biomass_ind = BiomassTind_g,
            nereo_height_ind = HeightT)

macro <- kelptrans2 %>%
  filter(Species == "Macro") %>%
  transmute(SiteName = SiteName,
            Transect = Transect,
            macro_biomass_ind = BiomassTind_g,
            macro_height_ind = HeightT)

ind <- merge(nereo, macro, by = c("SiteName", "Transect"), all=TRUE)

# Bringing in the kelp density (transect level) data
kelpjoin <- merge(kelptrans, dens2, by = c("SiteName", "Transect"), all=TRUE)

# Bringing in the kelp density (transect level) data my way
kelpjoin2 <- merge(ind, dens2, by = c("SiteName", "Transect"), all=TRUE)

# Getting to the average biomass / m2 area for each transect
kelptog <- kelpjoin %>%
  rowwise() %>% # To sum across rows
  mutate(BiomassTind_kg = (BiomassTind_g/1000), # Convert ave individual kelp biomass from g to kg 
         BiomassTkg = (BiomassTind_kg*Kelp)) %>% # Multiplying this by the transect kelp densities / m2
  ungroup() # Stop rowwise 


# try to get more accurate biomass per transect
kelptog2 <- kelpjoin2 %>%
  rowwise() %>% # To sum across rows
  mutate(macro_biomass_trans = Macro*macro_biomass_ind/1000,
         nereo_biomass_trans = Nereo*nereo_biomass_ind/1000,
         biomass_trans_mean = sum(macro_biomass_trans, nereo_biomass_trans, na.rm=TRUE),
         height_mean = mean(c(macro_height_ind, nereo_height_ind), na.rm=TRUE))

# make df to use in kelp_pee_analysis
# height mean isn't how I'd personally do it, but leave it for now
kelp_metrics <- kelptog2 %>%
  select(SiteName, Transect, Date, Time_start, Depth_m, RLS_dist, Macro_5m2, Nereo_5m2, Macro, Nereo, Kelp, macro_biomass_trans, nereo_biomass_trans, biomass_trans_mean, height_mean)

# write_csv(kelp_metrics, "Data/Team_kelp/Output_data/kelp_metrics2024.csv")


# ok now compare those estimates
em <- kelptog2 %>%
  transmute(SiteName = SiteName,
            Transect = Transect,
            Biomass_recalc = round(total_biomass_trans_kg, digits = 5))

claire <- kelptog %>%
  transmute(SiteName = SiteName,
            Transect = Transect,
            Biomass_original = round(BiomassTkg, digits = 5))

check <- merge(claire, em, by = c("SiteName", "Transect"), all=TRUE) %>%
  mutate(compare = ifelse(Biomass_recalc == Biomass_original, "same", "diff"))

# save metrics for future code

  

# Average to site level right away -----
# Averaging to site from individual kelp samples (i.e., averages of the 5 sampled individuals / transect)
kelptrans3 <- kelp5 %>%
  mutate(SiteName = as.factor(SiteName), 
         Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName, Species) %>% # Averaging to transect 
  summarise(HeightT = mean(Height_m, na.rm=TRUE),
            BiomassTind_g = mean(Biomass_g, na.rm=TRUE)) %>% # Ave individual biomass (g) / transect
  ungroup() %>%
  mutate(HeightT == ifelse(SiteName == "Less Dangerous Bay", 0, HeightT),
         BiomassTind_g == ifelse(SiteName == "Less Dangerous Bay", 0, BiomassTind_g))


# make df with column for macro biomass and column for nereo biomass
nereo <- kelptrans3 %>%
  filter(Species == "Nereo") %>%
  transmute(SiteName = SiteName,
            nereo_biomass_ind = BiomassTind_g,
            nereo_height_ind = HeightT)

macro <- kelptrans3 %>%
  filter(Species == "Macro") %>%
  transmute(SiteName = SiteName,
            macro_biomass_ind = BiomassTind_g,
            macro_height_ind = HeightT)

ind <- merge(nereo, macro, by = c("SiteName"), all=TRUE)

# Bringing in the kelp density (transect level) data my way
kelpjoin3 <- dens2 %>%
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
em2 <- em %>%
  group_by(SiteName) %>%
  summarise(mean_bio_recalc = mean(Biomass_recalc)) %>%
  mutate(mean_bio_recalc = round(mean_bio_recalc, digits = 5))

claire2 <- claire %>%
  group_by(SiteName) %>%
  summarise(mean_bio_og = mean(Biomass_original))%>%
  mutate(mean_bio_og = round(mean_bio_og, digits = 5))

kieran <- kelpjoin3 %>%
  select(SiteName, Biomass_kieran)%>%
  mutate(Biomass_kieran = round(Biomass_kieran, digits = 5))

check2 <- claire2 %>%
  left_join(em2, by = "SiteName") %>%
  left_join(kieran, by = "SiteName")

#  mutate(claire_em = mean_bio_recalc == mean_bio_og,
#         em_kieran = mean_bio_recalc == Biomass_kieran)

# write to csv
# write_csv(check2, "Output/Output_data/recalc_biomass.csv")

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
