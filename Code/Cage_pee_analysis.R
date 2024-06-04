# Code to analyze both cage experiments
# Em Lim
# Oct 17, 2023

# Load packages ----
library(tidyverse)
library(ggplot2)
library(TMB)
library(glmmTMB)
library(patchwork)
library(ggeffects)
library(DHARMa)
library(png)
library(ggtext)

renv::restore()
# Load functions
source("Code/Functions.R")

# Load cuke data ----
cuke_pee <- read_csv("Output/Output_data/cuke_cages.csv") %>%
  rename(nh4_avg = nh4_conc) %>%
  mutate(cukes = factor(cukes, levels = c("Control", "Mid", "High"),
                        labels = c("Control", "Medium", "Large")),
         depth_center = c(scale(depth, scale = FALSE))) %>%
  as.data.frame()

# load crab data ------
crab_pee <- read_csv("Data/Cage_experiment/crab_cage_pee.csv") %>%
  mutate(week = c(rep("one", 12), rep("two", 11)),
         treatment = factor(as.factor(treatment), 
                            levels = c("control", "mid", "large"),
                            labels = c("Control", "Medium", "Large"))) %>%
  pivot_longer( cols = c(nh4_conc3, nh4_conc2, nh4_conc1), names_to = "measurement", values_to = "nh4_conc") %>%
  mutate(day = c(rep(3:1, times = 12), rep(6:4, times = 11)))%>%
  as.data.frame()
# this csv is created in the Ch4_Crabs/Code/Crab_cages.R file

# Load and label images for plots
# Load cuke images
cuke <- readPNG("Images/cuke_no_background.png")
cuke_two <- readPNG("Images/two_cukes_no_background.png")

# label them
cuke_labels <- c(Control = "Control",
                 Medium = "<img src='Images/cuke_no_background.png' width='100' />",
                 Large = "<img src='Images/two_cukes_no_background.png' width='100' />")

# crab image
crab <- readPNG("Images/red_crab_no_background.png")

# label it
crab_labels <- c(Control = "Control",
                 Medium = "<img src='Images/red_crab_no_background.png' width='65' />",
                 Large = "<img src='Images/red_crab_no_background.png' width='100' />")


# Cuke stats -----
# start with regular gaussian, gamma is worse!
mod_cu <- glmmTMB(nh4_avg ~ cukes + depth_center,
               cuke_pee)

plot(simulateResiduals(mod_cu))

# add interactions?
mod_cu2 <- glmmTMB(nh4_avg ~ cukes * depth_center * line
                   - cukes:depth_stand:line,
                   cuke_pee)

# random effect of line?
mod_cu3 <- glmmTMB(nh4_avg ~ cukes + depth_center + (1|line),
                   cuke_pee)

AIC(mod_cu, mod_cu2, mod_cu3) # fewer interactions is better

# when we take out the random effect of line, the residuals look better.
# it doesn't change the model output

# look at model output
summary(mod_cu)


# Crab stats ----
# gamma (positive and continuous) is better than gaussian!
mod_cr_gamma <- glmmTMB(nh4_conc ~ treatment + (1|day) +(1|week), 
                 family = Gamma(link = "log"),
                 crab_pee)

mod_cr_gamma2 <- glmmTMB(nh4_conc ~ treatment + (1|week), 
                        family = Gamma(link = "log"),
                        crab_pee)

plot(simulateResiduals(mod_cr_gamma))

AIC(mod_cr_gamma, mod_cr_gamma2)
# OK so let's use the gamma!

summary(mod_cr_gamma)


#Graphing time folks-------
set.seed(1991)
# palettes
pal_cu <- viridis::viridis(6)
pal_c <- c(pal_cu[1], pal_cu[3], pal_cu[5])

pal_cu2 <- viridis::viridis(11)
pal_c2 <- c(pal_cu2[11], pal_cu2[10], pal_cu2[7])

# use ggpredict to get estimates for the cuke model
sum_cukes <- ggpredict(mod_cu, terms = "cukes") %>% 
  dplyr::rename(cukes = x,
                nh4_avg = predicted) %>% 
  as_tibble()

# use ggpredict to get estimates for crab model
sum_crabs <- ggpredict(mod_cr_gamma, terms = "treatment") %>% 
  dplyr::rename(treatment = x,
                nh4_avg = predicted) %>% 
  as_tibble()

# Make the plots
# plot cukes
cuke_plot <- dot_whisker(sum_data = sum_cukes, 
                         all_data = cuke_pee, 
                         x_var = cukes, 
                         y_var = nh4_avg, 
                         labels = cuke_labels,
                         pal = pal_c,
                         theme = "white") +
  ylim(c(0, 3.8)) +
  place_label("(a)")

# plot crabs
crab_plot <- dot_whisker(sum_data = sum_crabs, 
                         all_data = crab_pee,
                         x_var = treatment,
                         y_var = nh4_avg,
                         labels = crab_labels,
                         pal = pal_c,
                         theme = "white") +
  ylim(c(0, 3.8)) +
  place_label("(b)")

# plot together
cuke_plot + crab_plot & 
  plot_annotation(theme = theme(plot.background = 
      element_rect(color = "white", fill = "white")))

# ggsave("Output/Pres_figs/Fig5.png", device = "png", height = 6, width = 12, dpi = 400)
# og is 9 x 16

# Fig 4 white background for pub ----
# ggsave("Output/Pub_figs/Fig6.png", device = "png", height = 9, width = 16, dpi = 400)





# Map -----
library(sf)
library(osmdata)
library(ggimage) # for the pictures!

# Load not great shapefile ----
potato_map <- sf::st_read("Data/Potato_shapefiles/eez.shp") %>%
  st_sf() %>%
  st_set_crs(4326)

# Load GREAT shapefile
load("~/Documents/PhD/Ch2_Spatial_pee/Data/bc_map.Rdata")
bc_map <- slice # rename

# do this so i can trim the map margins
sf_use_s2(FALSE)

# Define colours
blue <- paste("#b9d1df", sep="")

# Choose coordinates 
site <- c("Scott's Bay", "Nova Harvest")
lat <- c(48.835583, 48.829500)
long <- c(-125.146361, -125.136250)
image <- c("Images/cuke_no_background.png", "Images/red_crab_no_background.png")

coords <- data.frame(site, lat, long, image) %>%
  st_as_sf(coords = c("long", "lat")) %>%
  st_set_crs(4326) 


# Make new map
ggplot() +
  geom_sf(data = potato_map, fill = "white", colour = blue) +
  geom_sf(data = coords, 
          alpha = 0.9,
          size = 9,
          aes(colour = site)) +
  coord_sf(xlim = c(-125.4, -125.0), ylim = c(48.80, 49), expand = FALSE)  +
  theme_black() +
  theme(panel.background = element_rect(fill = blue),
        panel.grid.major = element_line(color = blue),
        legend.position = "null") +
  labs(x = "Longitude", y = "Latitude",
       fill = expression(paste("NH"[4]^" +",(mu*M)))) +
  scale_x_continuous(breaks = seq(-125.4, -125.0, by = 0.1)) 



# make the map!
ggplot() +
  geom_sf(data = bc_map, fill = "white", colour = blue) +
  geom_sf(data = coords, 
          alpha = 0.9,
          size = 9,
          aes(colour = site)) +
  coord_sf(xlim = c(-125.25, -125.1), ylim = c(48.80, 48.9), expand = FALSE)  +
  theme_black() +
  theme(panel.background = element_rect(fill = blue),
        panel.grid.major = element_line(color = blue),
        legend.position = "null") +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-125.3, -125.1, by = 0.1))

# ggsave("Output/Figures/cage_site_map.png", device = "png", height = 9, width = 16, dpi = 400)

# Barkley sound for schematic
google_blue <- "#9bbff4"
google_green <- "#bbdaa4"

ggplot() +
  geom_sf(data = bc_map, fill = "white", colour = google_blue) +
  coord_sf(xlim = c(-125.3, -125), ylim = c(48.75, 48.95), expand = FALSE)  +
  theme_black() +
  theme(panel.background = element_rect(fill = google_blue),
        panel.grid.major = element_line(color = google_blue),
        legend.position = "null") +
  scale_x_continuous(breaks = seq(-125.3, -125, by = 0.1))

ggsave("Output/Figures/schematic_map.png", device = "png", height = 9, width = 16, dpi = 400)


# try to replace points with images!
# Help here: https://www.simoncoulombe.com/2020/11/animated-ships/


ggplot() +
  geom_sf(data = coastline_data$osm_lines, fill = blue, colour = "black") +
  # new images
  ggimage::geom_image(data = coords %>%
                        mutate(
                          proj_x= map_dbl( geometry, ~st_coordinates(.x)[1]), # trouver les coordonnées projetées
                          proj_y= map_dbl( geometry, ~st_coordinates(.x)[2])
                        ) %>% st_drop_geometry(), # dropper la géométrie
                      aes(x = proj_x, y = proj_y),
                      size = 0.25,
                      image = image) +
  coord_sf(xlim = c(-125.18, -125.1), ylim = c(48.82, 48.85), expand = FALSE)  +
  theme_black() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "white"),
        legend.position = "null") +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-125.18, -125.1, by = 0.05))


# Graveyard -----
sailboat <- pnw_palette("Sailboat")
csee_pal <- pnw_palette("Starfish")

# Old cuke x depth plot ----
# Effect of depth
ggplot(cuke_pee, aes(depth, nh4_avg)) +
  geom_point(aes(colour = cukes), size = 4) +
  geom_smooth(method = lm, aes(linetype = line), 
              colour = "black", alpha = 0.2) +
  scale_color_brewer(palette = "YlOrRd") +
  scale_y_reverse() +
  labs(y = "Ammonium concentration (umol)", x = "Depth (m)", 
       colour = "Sea cucumbers", linetype = "Line") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15)) 
