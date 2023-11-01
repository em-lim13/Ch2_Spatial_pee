# New R file for making maps of the RLS and Kelp pee data
# Em Lim
# Oct 17, 2023

# Load packages -----
library(tidyverse)
library(ggplot2)
library(ggspatial)
library(sf)

# Load functions ----
source("Code/Functions.R") # Length to weight function here!

# Load data ----
# these are csv files created by the rls and kelp analysis files

data_s <- read_csv("Output/Output_data/kelp_final.csv")

rls_final <- read_csv("Output/Output_data/rls_final.csv")

kelp_rls <- read_csv("Output/Output_data/kelp_rls.csv")
 

# Load GREAT shapefile ----
load("~/Documents/PhD/Ch2_Spatial_pee/Data/bc_map.Rdata")
bc_map <- slice # rename

# Load not great shapefile that renders quicker
potato_map <- sf::st_read("Data/Potato_shapefiles/eez.shp") %>%
  st_sf() %>%
  st_set_crs(4326)

# Make map without pies, just scaling size of point to %
sf_use_s2(FALSE)

# Define colours
blue <- paste("#b9d1df", sep="")


# RLS coords ------
rls_coords <- rls_final %>%
  select(site_code, survey_id, year, nh4_avg, abundance, tide_cat) %>%
  left_join(
    # use real coordinates not the RLS rounded ones
    read_csv("Data/RLS/RLS_data/true_coords.csv") 
  ) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  group_by(site_code) %>%
  mutate(nh4_overall_avg = mean(nh4_avg),
         nh4_min = min(nh4_avg),
         nh4_max = max(nh4_avg)) %>%
  ungroup() %>%
  mutate(Habitat = "Reef")

# just slack and ebb
coords_slack <- rls_coords %>%
  filter(tide_cat != "Flood") %>%
  group_by(site_code) %>%
  mutate(nh4_avg = mean(nh4_avg)) %>%
  ungroup()


# Kelp data coords -----
kelp_coords <- data_s %>%
  select(site_code, nh4_in_avg, nh4_out_avg, nh4_avg, kelp_sp, tide_cat) %>%
  left_join(
    kelp_rls %>%
      mutate(Habitat = "Kelp")
  ) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  group_by(site_code)


# Both sets of coords ----
all_coords <- rls_coords %>%
  transmute(site_code = site_code,
            nh4_avg = nh4_overall_avg,
            geometry = geometry,
            Habitat = Habitat) %>%
  unique() %>%
  rbind(
    (kelp_coords %>% transmute(nh4_avg = nh4_out_avg,
                               geometry = geometry,
                               Habitat = Habitat)) ) %>%
  # shrink the two sites over 1 uM to 2 uM so the scale is nicer
  mutate(Habitat = factor(as.factor(Habitat), levels = c("Reef", "Kelp")),
         nh4_avg = ifelse(nh4_avg > 2, 2, nh4_avg)) 



# Make maps! ------

# try functions for inset map

# just rls sites
barkley_rls <- site_map(lat_min = -125.6, lat_max = -124.95, long_min = 48.75, long_max = 49.1,
                    coord_data = all_coords %>%
                      filter(Habitat == "Reef"), map_data = bc_map) +
  guides(pch = "none")

van_isle_rls <- inset_map(rect_xmin = -125.6, rect_xmax = -124.95, rect_ymin = 48.75, rect_ymax = 49.1,
                      map_data = bc_map)

  
barkley_rls + 
  inset_element(
    van_isle_rls, 
    left = 0, 
    bottom = 0.5, 
    right = 0.5, 
    top = 1.03,
    align_to = 'panel'
  )

# ggsave("Output/Figures/rls_site_map.png", device = "png", height = 9, width = 16, dpi = 400)

# just kelp sites
barkley_kelp <- site_map(lat_min = -125.6, lat_max = -124.95, long_min = 48.75, long_max = 49.1,
                        coord_data = all_coords %>%
                          filter(Habitat == "Kelp"), map_data = bc_map) +
  guides(pch = "none")

van_isle_kelp <- inset_map(rect_xmin = -125.6, rect_xmax = -124.95, rect_ymin = 48.75, rect_ymax = 49.1,
                          map_data = bc_map)


barkley_kelp + 
  inset_element(
    van_isle_kelp, 
    left = 0, 
    bottom = 0.5, 
    right = 0.5, 
    top = 1.03,
    align_to = 'panel'
  )

 ggsave("Output/Figures/kelp_site_map.png", device = "png", height = 9, width = 16, dpi = 400)


# RLS site averages including all data
all_coords %>%
  filter(Habitat == "Reef") %>%
  map_daddy(nh4_avg, Habitat, potato_map)

#ggsave("Output/Figures/rls_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)


# Kelp only map
all_coords %>%
  filter(Habitat == "Kelp") %>%
  map_daddy(nh4_avg, Habitat, potato_map)

#ggsave("Output/Figures/kelp_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)


# Kelp + RLS map -----
map_daddy(all_coords, nh4_avg, Habitat, potato_map) 

#ggsave("Output/Figures/all_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)


# just avg of the slack and ebb measurements for RLS sites
map_daddy(coords_slack, nh4_avg, Habitat, potato_map) 

#ggsave("Output/Figures/nh4_slack_map.png", device = "png", height = 9, width = 16, dpi = 400)


# Cursed pimple map -----

# need to use full nh4 data, go run the rls analysis script to generate "rls_nh4"
coords_nh4 <- rls_nh4 %>%
  select(site_code, year, nh4_conc) %>%
  left_join(
    # use real coordinates not the RLS rounded ones
    read_csv("Data/RLS/RLS_data/true_coords.csv") 
  ) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  group_by(site_code) %>%
  mutate(nh4_avg = mean(nh4_conc),
         nh4_min = min(nh4_conc),
         nh4_max = max(nh4_conc)) %>%
  ungroup()


# Can I do something cursed? 
ggplot() +
  geom_sf(data = potato_map, fill = blue, colour = "white") +
  geom_sf(data = coords_nh4, 
          colour = "black",
          pch = 21,
          alpha = 0.9,
          size = 13.5,
          aes(fill = nh4_min)) +
  geom_sf(data = coords_nh4, 
          colour = "black",
          pch = 21,
          alpha = 0.9,
          size = 11,
          aes(fill = nh4_avg)) +
  geom_sf(data = coords_nh4, 
          colour = "black",
          pch = 21,
          alpha = 0.9,
          size = 4.5,
          aes(fill = nh4_max)) +
  coord_sf(xlim = c(-125.4, -125.0), ylim = c(48.80, 49), expand = FALSE)  +
  theme_black() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "white")) +
  viridis::scale_fill_viridis(option="magma", direction = -1,
                              guide = guide_colorbar(frame.colour = "white", ticks.colour = "white")) +
  labs(x = "Longitude", y = "Latitude",
       fill = expression(paste("NH"[4]^" +",(mu*M)))) +
  scale_x_continuous(breaks = seq(-125.4, -125.0, by = 0.1))

#ggsave("Output/Figures/cursed_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)
