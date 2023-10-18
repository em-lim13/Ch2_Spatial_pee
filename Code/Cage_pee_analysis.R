# Code to analyze both cage experiments
# Em Lim
# Oct 17, 2023

# Load packages ----
library(tidyverse)
library(ggplot2)
library(TMB)
library(glmmTMB)
library(patchwork)
library(PNWColors)

# Load functions
source("Code/Functions.R")

# Load cuke data ----
cuke_pee <- read_csv("Output/Output_data/cuke_cages.csv") %>%
  rename(nh4_avg = nh4_conc) %>%
  mutate(cukes = factor(cukes, levels = c("Control", "Mid", "High"),
                        labels = c("Control", "Medium", "Large")))

# load crab data ------
crab_pee <- read_csv("Data/crab_cage_pee.csv") %>%
  mutate(week = c(rep("one", 12), rep("two", 11)),
         treatment = factor(as.factor(treatment), 
                            levels = c("control", "mid", "large"),
                            labels = c("Control", "Medium", "Large")))
# this csv is created in the Ch4_Crabs/Code/Crab_cages.R file

# Cuke stats -----
# start with regular gaussian 
mod_cu <- glmmTMB(nh4_avg ~ cukes + depth_stand * line, 
               cuke_pee)

# plot(simulateResiduals(mod_cu))

# try gamma (positive and continuous)
mod_cu_gamma <- glmmTMB(nh4_avg ~ cukes + depth_stand * line, 
                     family = Gamma(link = "log"),
                     cuke_pee)

# plot(simulateResiduals(mod_cu_gamma))
# OK so let's use the gaussian! 

summary(mod_cu)


# Crab stats ----
# start with regular gaussian 
mod_cr <- glmmTMB(nh4_avg ~ treatment + (1|week), 
                     crab_pee)

# plot(simulateResiduals(mod_cr))

# try gamma (positive and continuous)
mod_cr_gamma <- glmmTMB(nh4_avg ~ treatment + (1|week), 
                 family = Gamma(link = "log"),
                 crab_pee)

# plot(simulateResiduals(mod_cr_gamma))
# OK so let's use the gamma!

summary(mod_cr_gamma)

#Graphing time folks-------

# palettes
sailboat <- pnw_palette("Sailboat")
csee_pal <- pnw_palette("Starfish")


# use ggpredict to get estimates for the cuke model
sum_cukes <- ggpredict(mod_cu, terms = c("cukes")) %>% 
  dplyr::rename(cukes = x,
                nh4_avg = predicted) %>% 
  as_tibble()

# use ggpredict to get estimates for crab model
sum_crabs <- ggpredict(mod_cr_gamma, terms = c("treatment")) %>% 
  dplyr::rename(treatment = x,
                nh4_avg = predicted) %>% 
  as_tibble()


# Make the plots
# plot cukes
cuke_plot <- dot_whisker(sum_data = sum_cukes, all_data = cuke_pee, 
                         x_var = cukes, y_var = nh4_avg) +
  labs(x = "", title = "Sea cucumbers")

# plot crabs
crab_plot <- dot_whisker(sum_data = sum_crabs, all_data = crab_pee, 
            x_var = treatment, y_var = nh4_avg) +
  labs(x = "", title = "Crabs")

# plot together
cuke_plot + crab_plot

ggsave("Output/Figures/both_cages.png", device = "png", height = 9, width = 16, dpi = 400)


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

