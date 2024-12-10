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
  mutate(BiomassTind_kg = (BiomassTind_g/1000)) %>% # Convert ave individual kelp biomass from g to kg 
  mutate(BiomassTkg = (BiomassTind_kg*Kelp)) %>% # Multiplying this by the transect kelp densities / m2
  ungroup() # Stop rowwise 
kelptog[sapply(kelptog, is.nan)] <- 0 # NaNs to 0s for working with


# try to get more accurate biomass per transect
kelptog2 <- kelpjoin2 %>%
  rowwise() %>% # To sum across rows
  mutate(macro_biomass_trans = Macro*macro_biomass_ind,
         nereo_biomass_trans = Nereo*nereo_biomass_ind,
         total_biomass_trans = sum(macro_biomass_trans, nereo_biomass_trans, na.rm=TRUE),
         total_biomass_trans_kg = total_biomass_trans/1000)


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
write_csv(check2, "Output/Output_data/recalc_biomass.csv")

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

#### Kelp area cleaning ----

# Loading packages for working with sf and mapping objects
library(scales)
library(sf)
library(ggsn)
library(maptools)
library(rgeos)

######## Prepping the backing map (in case you want to plot out any of the .shp files)

## setting projections at outset 
proj <- st_crs(3005) # ESPG for BC/Albers 
latlong <- st_crs(4326) # for baseline/existing WGS84 ("+proj=longlat +datum=WGS84 +no_defs")

## setting map extent for Barkley Sound 
ymax <- 48.922
ymin <- 48.80
xmax <- -125.05
xmin <- -125.26

## making corners for the area of interest 
corners <- st_multipoint(matrix(c(xmax,xmin,ymax,ymin),ncol=2)) %>% 
  st_sfc(crs=latlong) %>% # Origin crs as WGS84
  st_sf() %>%
  st_transform(proj) # Projecting into working crs (BC/ALBERS)
plot(corners)


######## Loading in the kelp forest shape files

# Get all files with the .shp extension from the correct folder (i.e., the QGIS cleaned versions)
path <- "./MSc_data/Data_new/Kelp_area/CleanGPX" # set path to the folder

n <- list.files(path, pattern = "*shp", full.names=FALSE) # Pull only the relevant names

shps <- list.files(path, pattern = "*shp", full.names=TRUE) # Pull full file location names
lshps <- as.list(shps) # Turn character vector into a list

names(lshps) <- sub(" area.shp", "", n) # Assign relevant names (cleaned) to elements of the list

f <- function(x){
  st_read(x, crs=latlong) # Base crs of the GPS used to create the shp files is WGS84 
}

allshps <- lapply(lshps, f) # Apply function to read in shape files as sf objects with crs WGS84

f2 <- function(x){
  st_transform(x, crs=proj) # TO reproject to project crs with units in metres (BC/Albers) 
}

allshps <- lapply(allshps, f2) # Apply function to translate shape files to BC/Albers


list2env(allshps, envir=.GlobalEnv) # Bring each unique sf object to the working environment for plotting             


####### Calculating area of each kelp forest shape file

# I know this is not the most efficient code solution,
# but it does work to produce the areas from each .shp file
a <- cbind(st_area(allshps[[1]]), as.character(names(allshps[1])))
b <- cbind(st_area(allshps[[2]]), as.character(names(allshps[2])))
c <- cbind(st_area(allshps[[3]]), as.character(names(allshps[3])))
d <- cbind(st_area(allshps[[4]]), as.character(names(allshps[4])))
e <- cbind(st_area(allshps[[5]]), as.character(names(allshps[5])))
f <- cbind(st_area(allshps[[6]]), as.character(names(allshps[6])))
g <- cbind(st_area(allshps[[7]]), as.character(names(allshps[7])))
h <- cbind(st_area(allshps[[8]]), as.character(names(allshps[8])))
i <- cbind(st_area(allshps[[9]]), as.character(names(allshps[9])))
j <- cbind(st_area(allshps[[10]]), as.character(names(allshps[10])))
k <- cbind(st_area(allshps[[11]]), as.character(names(allshps[11])))
l <- cbind(st_area(allshps[[12]]), as.character(names(allshps[12])))
m <- cbind(st_area(allshps[[13]]), as.character(names(allshps[13])))
n <- cbind(st_area(allshps[[14]]), as.character(names(allshps[14])))
o <- cbind(st_area(allshps[[15]]), as.character(names(allshps[15])))
p <- cbind(st_area(allshps[[16]]), as.character(names(allshps[16])))
q <- cbind(st_area(allshps[[17]]), as.character(names(allshps[17])))
r <- cbind(st_area(allshps[[18]]), as.character(names(allshps[18])))
s <- cbind(st_area(allshps[[19]]), as.character(names(allshps[19])))
t <- cbind(st_area(allshps[[20]]), as.character(names(allshps[20])))
u <- cbind(st_area(allshps[[21]]), as.character(names(allshps[21])))
v <- cbind(st_area(allshps[[22]]), as.character(names(allshps[22])))
w <- cbind(st_area(allshps[[23]]), as.character(names(allshps[23])))


# Generating data frame of kelp forest areas
areagrp <- data.frame(rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w)) %>%
  mutate(Area_m2 = as.numeric(X1), SiteName = as.factor(X2)) %>%
  dplyr::select(-c(X1, X2))


# Removing the site 'redos' for current purposes
areagrp <- areagrp %>%
  filter(!str_detect(SiteName, 'REDO'))


#### Final data joining ----

### Merging density, canopy height, biomass, and area data 
kelpdat <- merge(kelpgrp, areagrp, by = "SiteName", all=TRUE)

# Filtering out the no kelp sites (Less Dangerous Bay & Wizard I North)
kelpdat <- kelpdat %>%
  filter(SiteName != "Less Dangerous Bay") %>%
  filter(SiteName != "Wizard Islet North")

# Filtering out the 'redos' of site sampling in the fall
kelpdat <- kelpdat %>%
  filter(SiteName != "Ross Islet 2 REDO") %>%
  filter(SiteName != "Second Beach South REDO") %>%
  filter(SiteName != "Less Dangerous Bay REDO")


# saving a .csv file of the kelp metrics by site
write_csv(kelpdat, "Data/Team_kelp/Output_data/kelp_metrics_2022_final.csv", row.names=F)