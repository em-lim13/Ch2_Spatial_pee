# New R file for making maps of the RLS and Kelp pee data
# Em Lim
# Oct 17, 2023, last updated Dec 2024

# Load packages -----
library(tidyverse)
library(ggplot2)
library(ggspatial)
library(sf)
library(viridis)

# Load functions ----
source("Code/Functions.R") 

# Load data ----
# these are csv files created by the rls and kelp analysis files
kelp_map_data <- read_csv("Output/Output_data/kelp_map_data.csv")

rls_final <- read_csv("Output/Output_data/rls_final.csv")

# load all rls pee data
# this includes the June and July samples, which were excluded from the main RLS analysis
rls_nh4_all <- rbind(read_csv("Output/Output_data/RLS_nh4_2021.csv"),
                     read_csv("Output/Output_data/RLS_nh4_2022.csv"),
                     read_csv("Output/Output_data/RLS_nh4_2023.csv")) %>%
  rename(site_code = site_ID) %>%
  group_by(site_code) %>%
  mutate(nh4_overall_avg_all = mean(nh4_conc)) %>%
  ungroup() %>%
  select(site_code, nh4_overall_avg_all) %>%
  unique()

# Load GREAT shapefile ----
hakai_map <- sf::st_read("Data/Hakai_coast/COAST_TEST2.shp") %>%
  st_sf() %>%
  st_set_crs(4326)

# Make map without pies, just scaling size of point to %
sf_use_s2(FALSE)

# Define colours
blue <- paste("#b9d1df", sep="")


# RLS coords ------
rls_coords <- rls_final %>%
  group_by(site_code) %>%
  # this is the mean of the May samples
  mutate(nh4_avg = mean(nh4_avg)) %>%
  ungroup() %>%
  select(site_code, nh4_avg) %>%
  unique() %>%
  left_join(
    # use real coordinates not the RLS rounded ones
    read_csv("Data/RLS/RLS_data/true_coords.csv"), by = "site_code") %>%
  mutate(Habitat = "Reef") %>%
  left_join(rls_nh4_all, by = "site_code")

# Kelp data coords -----
kelp_coords <- kelp_map_data %>%
  transmute(site_code = site_code,
            site_name = site_name,
            nh4_avg = nh4_avg,
            nh4_overall_avg_all = nh4_avg,
            longitude = longitude,
            latitude = latitude,
            Habitat = "Kelp")

# Both sets of coords ----
all_coords <- rbind(rls_coords, kelp_coords) %>%
  # stupidly manually jitter points that are on top of each other
  mutate(latitude = case_when(site_code == "KCCA12" ~ 48.854448, 
                              site_code == "KCCA22" ~ 48.825960,
                              site_code == "KCCA19" ~ 48.860970,
                              site_code == "KCCA9" ~ 48.855602,
                              site_code == "BMSC8" ~ 48.956866,
                              site_code == "BMSC11" ~ 48.857512,
                              site_code == "BMSC12" ~ 48.858478,
                              site_code == "BMSC4" ~ 48.815043,
                              site_code == "KCCA14" ~ 48.878305,
                              site_code == "KCCA7" ~ 48.837319,
                              TRUE ~ as.numeric(latitude)),
         longitude = case_when(site_code == "KCCA12" ~ -125.168718, 
                               site_code == "KCCA22" ~ -125.19822,
                               site_code == "KCCA19" ~ -125.159304,
                               site_code == "BMSC11" ~ -125.157674,
                               site_code == "KCCA1" ~ -125.156543,
                               site_code == "BMSC8" ~ -125.151189,
                               site_code == "BMSC12" ~ -125.161694,
                               site_code == "BMSC4" ~ -125.177346,
                               site_code == "KCCA7" ~ -125.212146,
                               TRUE ~ as.numeric(longitude))) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  # shrink the one site over 2 uM to 2 uM so the scale is nicer
  mutate(Habitat = factor(as.factor(Habitat), levels = c("Reef", "Kelp")),
         nh4_avg = ifelse(nh4_avg > 2, 2, nh4_avg),
         nh4_overall_avg_all = ifelse(nh4_overall_avg_all > 2, 2, nh4_overall_avg_all),
         dummy = 2) 



# Make maps! ------

# Fig. 1d -----
# All sites
map_daddy(long_min = -125.375,
          long_max = -125.025, 
          lat_min = 48.801, 
          lat_max = 48.965, 
          coord_data = all_coords, 
          nh4_var = nh4_overall_avg_all, 
          kelp_var = Habitat,
          point_size = 6.5, 
          map_file = hakai_map,
          invert = FALSE,
          white_background = TRUE)

# size adjusted
# ggsave("Output/Pub_figs/ppt_shame_folder/Fig.1d.png", device = "png", height = 9, width = 12.606299, dpi = 400)

# old size
# ggsave("Output/Pub_figs/ppt_shame_folder/Fig.1.png", device = "png", height = 9, width = 16, dpi = 400)

#  5 inches in width and 6 inches in height

# Other maps ----
# Barkley Sound map
map_daddy_np(long_min= -127,
             long_max= -123, 
             lat_min  = 48.5, 
             lat_max  = 51, 
             map_file = hakai_map,
             invert = FALSE) +
  # add rectangle for zoomed in part
  geom_rect(aes(xmin = -125.4, xmax = -125.0, ymin = 48.80, ymax = 49),
            color = "red", fill = NA, inherit.aes = FALSE) 

# ggsave("Output/Figures/barkley_sound_map.png", device = "png", height = 9, width = 16, units = "cm", dpi = 400)


# RLS site map
map_daddy(long_min= -125.4,
          long_max= -125.0, 
          lat_min  = 48.75, 
          lat_max  = 49, 
          coord_data = all_coords %>% filter(Habitat == "Reef"), 
          nh4_var = dummy, 
          point_size = 4, 
          kelp_var = Habitat,
          map_file = hakai_map) +
  guides(pch = "none",
         fill = "none")

# ggsave("Output/Figures/rls_site_map.png", device = "png", height = 9, width = 16, dpi = 400)


# RLS nh4 data!
map_daddy(long_min= -125.4,
          long_max= -125.0, 
          lat_min  = 48.80, 
          lat_max  = 49, 
          coord_data = all_coords %>% filter(Habitat == "Reef"), 
          nh4_var = nh4_avg, 
          point_size = 9, 
          kelp_var = Habitat,
          map_file = hakai_map,
          invert = FALSE,
          white_background = FALSE) +
  guides(pch = "none")

# ggsave("Output/Figures/rls_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)


# Kelp site map
map_daddy(long_min= -125.4,
          long_max= -125.0, 
          lat_min  = 48.80, 
          lat_max  = 49, 
          coord_data = all_coords %>% filter(Habitat == "Kelp"), 
          nh4_var = dummy, 
          point_size = 4, 
          kelp_var = Habitat,
          map_file = hakai_map,
          invert = FALSE,
          white_background = FALSE) +
  guides(pch = "none",
         fill = "none")

# ggsave("Output/Figures/kelp_site_map.png", device = "png", height = 9, width = 16, dpi = 400)


# Kelp nh4 data!
map_daddy(long_min= -125.3,
          long_max= -125.0, 
          lat_min  = 48.78, 
          lat_max  = 48.93, 
          coord_data = all_coords %>% filter(Habitat == "Kelp"), 
          nh4_var = nh4_avg, 
          point_size = 9, 
          kelp_var = Habitat,
          map_file = hakai_map,
          invert = FALSE,
          white_background = FALSE) +
  guides(pch = "none")

#ggsave("Output/Figures/kelp_nh4_map2.png", device = "png", height = 9, width = 16, dpi = 400)


# distance between coords  -----
library(geodist)
distance_matrix <- geodist(rls_coords, measure = 'geodesic' )/1000 #converting it to km
colnames(distance_matrix) <- rls_coords$site_name
rownames(distance_matrix) <- rls_coords$site_name

distance_matrix_lon <- as.data.frame.table(distance_matrix, responseName = "value") %>%
  filter(value != 0)

max(distance_matrix_lon$value) # greatest distance is between Hosie and Wouwer = 24.24 km
min(distance_matrix_lon$value) # smallest distance is Baeria Rocks North Island Southside to Baeria Rocks North Island Northside = 65.55 m
