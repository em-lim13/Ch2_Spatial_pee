# Play with abalone data
# Em Lim
# May 2024

# I cleaned up the abalone data from the RLS data in the "RLS_data_scrub.R" file and saved the clean csv "abalone_RLS.csv"

# this script is the beginning of an exploration of that data

# run this to install the packages!
renv::restore()

# Load packages and data -----
# Load packages
library(readr)
library(tidyr)
library(dplyr) # for general data wrangling
library(lubridate) # for looking at dates

library(glmmTMB) # The swiss army knife of modeling packages
library(DHARMa) # inspect model residuals/check assumptions
library(ggeffects) # for extracting predictions and running post hoc tests

library(ggplot2) # Plotting data
library(patchwork) # Arrange multiple plots together
library(visreg) # plot model predictions
library(viridisLite) # colour palettes

library(ggspatial)
library(sf)

# source file where my functions live
source("Code/Functions.R")

# Load abalone data 
abalone <- read_csv("Data/abalone_RLS.csv")

# Where are abalone? -----

# this gives you a dataframe of the coordinates of each site and an average abalone density at each of them in case anyone wants to map this. careful bc some sites had multiple transects, some deep and some shallow and this will just average across those
aba_coords <- abalone %>%
  group_by(site_code) %>%
  mutate(mean_density = mean(density),
         mean_total = mean(transect_total),
         ) %>%
  ungroup() %>%
  select(site_code, site_name, survey_type, latitude, longitude, mean_total, mean_density) %>%
  unique()

write_csv(aba_coords, "Data/abalone_coordinates.csv")


# Mapping -----
# load map
# Load packages -----
library(tidyverse)
library(ggplot2)
library(ggspatial)
library(sf)

# Load functions ----
source("Code/Functions.R") # Length to weight function here!


