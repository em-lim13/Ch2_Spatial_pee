# Code to analyze the May 25 - 27 cage experiment
#May 28, 2021
#By: Emily Leedham and Em Lim

#Load packages and data ------

#Load packages
library(tidyverse)
library(visreg)
library(ggplot2)
library(lme4)
library(lmerTest)
library(RColorBrewer)
library(PNWColors)
source("Code/theme_black.R")

theme_set(theme_bw())

#Load cage data

cages <- read_csv("Data/Cage_experiment/2021_05_28_cage_samples.csv")

#Load fluorometry data

glow <-  read_csv("Data/Cage_experiment/2021_05_28_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

#Load standard curve data
standard <- read_csv("Data/Cage_experiment/2021_05_28_standard curve.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

# # STANDARD CURVE---------------------------------------------
#do the calculations for the standard curve

standard_f <- standard %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading

#linear model between the fluorometer reading and actual concentration of NH4
sc_mod1 <- lm(nh4_conc_final_umol_L ~ mean_FLU, data = standard_f)
#check this summary for the R squared- the relationship should be TIGHT, so if it's lower than like 0.98, there might be some contamination
summary(sc_mod1)
#visualize curve to make sure it looks right
ggplot(standard_f, aes(mean_FLU, nh4_conc_final_umol_L)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#extract intercept and slope to use later to determine sample concentrations
#based on fluorometry readings
int <- coef(sc_mod1)[1]
slope <- coef(sc_mod1)[2]


# Calculate matrix effects --------------------------------------

#joining cages and glow
cages_a<- cages %>%
  left_join(glow, by = c("bottle", "cage_ID")) %>%
  mutate(Fsm_zero = mean_FLU)

#getting the matrix data
cages_matrix <- cages_a %>% slice_tail(n = 4) %>%
  transmute(
  line = line,
  mean_FLU = mean_FLU,
  sample = c("sample", "spike", "sample", "spike")
)

cages_spike <- cages_matrix %>%
  filter(sample == "spike")%>%
transmute(line = line,
          Fsm_spike = mean_FLU)

cages_sample <- cages_matrix %>% 
  filter(sample == "sample") %>% 
  transmute(line = line, 
            Fsm_zero = mean_FLU)

#smooshing back together

cages_matrix_2 <- cages_spike %>%
  left_join(cages_sample, by = "line")

Fst_zero <- standard_f$mean_FLU[standard_f$nh4_vol_uL == "0"]
Fst_spike <- standard_f$mean_FLU[standard_f$nh4_vol_uL == "200"]

cages_matrix_3 <- cages_matrix_2 %>% 
    mutate(
    Fst_zero = Fst_zero,
    Fst_spike = Fst_spike,
    matrix = 100 * ((Fst_spike - Fst_zero - (Fsm_spike - Fsm_zero))/
                      (Fst_spike - Fst_zero))
    ) %>% 
  select(c(line, matrix))

cages_f <- cages_a %>%
  left_join(cages_matrix_3, by = "line") %>%
  slice(1:(n()-4)) %>%
  mutate(Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
         Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100)),
         int = int, 
         slope = slope,
         nh4_conc = int + slope * Fsm_cor_Taylor,
         cukes = factor(cukes, levels = c("0", "1", "2")),
         cuke_vol_cm3 = (pi*cuke_len1*(cuke_width1/(2*pi))^2) +
           (pi*cuke_len2*(cuke_width2/(2*pi))^2)
         )

#Graphing time folks--------

line_lab <- c("Line A", "Line B")
names(line_lab) <- c("A", "B")

# Effect of cuke treatment
ggplot(cages_f, aes(cukes, nh4_conc, fill = cukes)) +
  geom_boxplot() +
  facet_wrap(~line, 
             labeller = labeller(line = line_lab)) +
  scale_fill_brewer(palette = "YlOrRd") +
  theme(legend.position = "none",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 16, colour = "black"),
        strip.text.x = element_text(size = 15) ) +
labs(x = "Sea cucumber treatment", y = "Ammonium concentration (umol)")

# in black
sailboat <- pnw_palette("Sailboat")
csee_pal <- pnw_palette("Starfish")

ggplot(cages_f, aes(cukes, nh4_conc, fill = cukes, colour = cukes)) +
  geom_boxplot(colour = "white", alpha = 0.7) +
  geom_point(alpha = 0.9, size = 4) +
  scale_fill_manual(values = (sailboat)) +
  scale_colour_manual(values = (sailboat)) +
  theme_black() +
  theme(legend.position = "none") +
  labs(x = "Sea cucumber treatment", y = "Ammonium concentration (umol/L)")

#ggsave("Output/Figures/Cage_exp_cuke_nh4.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# Effect of cuke vol
ggplot(cages_f, aes(cuke_vol_cm3, nh4_conc, colour = cukes)) +
  geom_point(size = 4) +
  facet_wrap(~line, 
             labeller = labeller(line = line_lab)) +
  scale_color_brewer(palette = "YlOrRd") +
  labs(x = "Total cuke volume (cm^3)", y = "Ammonium concentration (umol)", 
       colour = "Cukes") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
          strip.text.x = element_text(size = 15)
  )

#ggsave("Output/Figures/Cage_exp_cuke_vol_nh4.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# Effect of depth
ggplot(cages_f, aes(nh4_conc, depth)) +
  geom_point(aes(colour = cukes), size = 4) +
  geom_smooth(method = lm, aes(linetype = line), 
              colour = "black", alpha = 0.2) +
  scale_color_brewer(palette = "YlOrRd") +
  scale_y_reverse() +
  labs(x = "Ammonium concentration (umol)", y = "Depth (m)", 
       colour = "Sea cucumbers", linetype = "Line") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15)) 

# in black
ggplot(cages_f, aes(nh4_conc, depth)) +
  geom_point(aes(colour = cukes), size = 4) +
  geom_smooth(method = lm, aes(linetype = line), 
              colour = "white", alpha = 0.2) +
  scale_color_brewer(palette = "YlOrRd") +
  scale_y_reverse() +
  labs(x = "Ammonium concentration (umol)", y = "Depth (m)", 
       colour = "Cukes", linetype = "Line") +
  theme_black()

#ggsave("Output/Figures/Cage_exp_depth_nh4.png", device = "png",
#       height = 9, width = 16, dpi = 400)


#Lets do some stats :/-----
model <- lm(nh4_conc ~ cukes + depth + line, cages_f)
summary(model)
visreg(model)

