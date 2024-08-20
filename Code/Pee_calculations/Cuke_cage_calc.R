# Code to calculate nh4 from the May 25 - 27 cage experiment
#May 28, 2021
#By: Emily Leedham and Em Lim

#Load packages ------

#Load packages
library(tidyverse)
source("Code/Functions.R")

theme_set(theme_bw())

# Load data -----
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

# final data ------
cages_f <- cages_a %>%
  left_join(cages_matrix_3, by = "line") %>%
  slice(1:(n()-4)) %>%
  mutate(Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
         Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100)),
         int = int, 
         slope = slope,
         nh4_conc = int + slope * Fsm_cor_Taylor,
         cukes = factor(cukes, levels = c("0", "1", "2"),
                        labels = c("Control", "Mid", "High")),
         cuke_vol_cm3 = (pi*cuke_len1*(cuke_width1/(2*pi))^2) +
           (pi*cuke_len2*(cuke_width2/(2*pi))^2),
         depth_stand = c(scale(depth, scale = FALSE))
         )

cages_for_csv <- cages_f %>%
  select(cage_ID, line, depth, cukes, nh4_conc, depth_stand) %>%
  as.data.frame()

# write_csv(cages_for_csv, "Output/Output_data/cuke_cages.csv")