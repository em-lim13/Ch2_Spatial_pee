##Code to analyse Tubs 2 and 1: electric boogaloo
##July 2, 2021
##By: Emily Leedham and Em

##Load packages and data -----

#Load packages
library(tidyverse)
library(visreg)
library(ggplot2)
library(TMB)
library(glmmTMB)
library(PNWColors)
library(viridis)
library(patchwork)
library(ggeffects)
source("Code/Functions.R")

theme_set(theme_bw())

##Code to analyze The Experiment of Many Tubs June 21, 2021----
tubs <- read_csv("Data/Emily_flow/2021_06_21_emily_flow_exp.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(flu_1, flu_2, flu_3)),
         flow_s = as.factor(flow_s), 
         cukes_num = as.factor(cukes_num),
         mean_flow = rowMeans(cbind(flow_s, flow_e)))

#Load standard curve data
standard <- read_csv("Data/Emily_flow/2021_06_23_standard_curve_data.csv") 

##Making the standard curve -----
standard_f <- standard %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading %>%


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

##Calculating the matrix effects -----
Fst_zero <- standard_f$mean_FLU[standard_f$nh4_vol_uL == "0"]
Fst_spike <- standard_f$mean_FLU[standard_f$nh4_vol_uL == "200"]

matrix <- tubs %>%
  filter(sample_matrix == "matrix") %>%
  transmute(Fsm_spike = mean_FLU, tub = tub)%>%
  left_join(tubs, by = "tub") %>%
  filter(sample_matrix == "sample") %>%
  transmute(Fsm_spike = Fsm_spike,
            Fsm_zero = mean_FLU,
            Fst_spike = Fst_spike,
            Fst_zero = Fst_zero, 
            matrix = 100 * ((Fst_spike - Fst_zero - (Fsm_spike - Fsm_zero))/
                              (Fst_spike - Fst_zero)) ) %>%
  select(matrix)

##Calculating NH4 conc
tubs_nh4 <- cbind(tubs, matrix) %>%
  mutate(Fsm_zero = mean_FLU,
         Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100)),
         int = int, 
         slope = slope,
         NH4_conc = int + slope * Fsm_cor_Taylor,
         date = "one") %>%
  filter(sample_matrix != c("matrix")) %>%
  mutate(total_mass = cuke_mass1 +cuke_mass2 +cuke_mass3 + cuke_mass4)

#calculate how much pee was added
tubs_controls <- tubs_nh4 %>%
  filter(cukes_num == "0") %>%
  transmute(NH4_conc_control = mean(NH4_conc), 
            sample_matrix = sample_matrix) %>%
  slice_head()

tubs_nh4_added1 <- tubs_nh4 %>%
  left_join(tubs_controls, by = "sample_matrix") %>%
  mutate(nh4_added = NH4_conc - NH4_conc_control,
         date = "one") %>%
  select(cukes_num, flow_s, mean_flow, total_mass, NH4_conc, nh4_added, date)


# FLOW CODE TAKE 2 -----------
# The fluorometry is split into am and pm because Emily had to go to the nurse part way through fluorometry and by the time Em got back to the lab to finish the OPA had sat for long enough that the standards had to be rerun
# Therefore there are two sets of standards and matrix calculations

#load am and pm data for fluoro
tubs_am <- read.csv("Data/Emily_flow/2021_06_30_emily_flow_exp2_am.csv")%>%
  mutate(mean_FLU = rowMeans(cbind(flu_1, flu_2, flu_3)),
         flow_s = as.factor(flow_s), 
         cukes_num = as.factor(cukes_num),
         flow_e = (((500/flow_e_sec) - 13.594)/3.56),
         mean_flow = rowMeans(cbind(flow_s, flow_e)))

tubs_pm <- read.csv("Data/Emily_flow/2021_06_30_emily_flow_exp2_pm.csv")%>%
  mutate(mean_FLU = rowMeans(cbind(flu_1, flu_2, flu_3)),
         flow_s = as.factor(flow_s), 
         cukes_num = as.factor(cukes_num),
         flow_e = (((500/flow_e_sec) - 13.594)/3.56),
         mean_flow = rowMeans(cbind(flow_s, flow_e)))


#Load standard curve data for am and pm fluoro
standard_am <- read.csv("Data/Emily_flow/2021_06_30_emily_flow_epx2_standard_am.csv")

standard_pm <- read.csv("Data/Emily_flow/2021_06_30_emily_flow_epx2_standard_pm.csv")

##Making the standard curve FOR AM DATA -----
standard_am_f <- standard_am %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading %>%


#linear model between the fluorometer reading and actual concentration of NH4
sc_mod1_am <- lm(nh4_conc_final_umol_L ~ mean_FLU, data = standard_am_f)
#check this summary for the R squared- the relationship should be TIGHT, so if it's lower than like 0.98, there might be some contamination
summary(sc_mod1_am)
#visualize curve to make sure it looks right
ggplot(standard_am_f, aes(mean_FLU, nh4_conc_final_umol_L)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#extract intercept and slope to use later to determine sample concentrations
#based on fluorometry readings
int_am <- coef(sc_mod1_am)[1]
slope_am <- coef(sc_mod1_am)[2]


##Calculating the matrix effects FOR THE AM DATA-----
Fst_zero_am <- standard_am_f$mean_FLU[standard_am_f$nh4_vol_uL == "0"]
Fst_spike_am <- standard_am_f$mean_FLU[standard_am_f$nh4_vol_uL == "200"]

matrix_am <- tubs_am %>%
  filter(sample_matrix == "matrix") %>%
  transmute(Fsm_spike = mean_FLU, tub = tub)%>%
  left_join(tubs_am, by = "tub") %>%
  filter(sample_matrix == "sample") %>%
  transmute(Fsm_spike = Fsm_spike,
            Fsm_zero = mean_FLU,
            Fst_spike_am = Fst_spike_am,
            Fst_zero_am = Fst_zero_am, 
            matrix = 100 * ((Fst_spike_am - Fst_zero_am - (Fsm_spike - Fsm_zero))/
                              (Fst_spike_am - Fst_zero_am)) ) %>%
  select(matrix) #### matrix is negative and it shouldn't be

##Calculating NH4 conc FOR THE AM DATA
tubs_nh4_am <- cbind(tubs_am, matrix_am) %>%
  mutate(Fsm_zero = mean_FLU,
         Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100)),
         int_am = int_am, 
         slope_am = slope_am,
         NH4_conc = int_am + slope_am * Fsm_cor_Taylor) %>%
  filter(sample_matrix != c("matrix")) %>%
  mutate(total_mass = cuke_mass1 +cuke_mass2 +cuke_mass3 + cuke_mass4,
         time = "am") %>%
  select(- c(int_am, slope_am))

##Making the standard curve FOR PM DATA -----

standard_pm_f <- standard_pm %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading %>%


#linear model between the fluorometer reading and actual concentration of NH4
sc_mod1_pm <- lm(nh4_conc_final_umol_L ~ mean_FLU, data = standard_pm_f)
#check this summary for the R squared- the relationship should be TIGHT, so if it's lower than like 0.98, there might be some contamination
summary(sc_mod1_pm)
#visualize curve to make sure it looks right
ggplot(standard_pm_f, aes(mean_FLU, nh4_conc_final_umol_L)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#extract intercept and slope to use later to determine sample concentrations
#based on fluorometry readings
int_pm <- coef(sc_mod1_pm)[1]
slope_pm <- coef(sc_mod1_pm)[2]


##Calculating the matrix effects FOR THE PM DATA-----
Fst_zero_pm <- standard_pm_f$mean_FLU[standard_pm_f$nh4_vol_uL == "0"]
Fst_spike_pm <- standard_pm_f$mean_FLU[standard_pm_f$nh4_vol_uL == "200"]

matrix_pm <- tubs_pm %>%
  filter(sample_matrix == "matrix") %>%
  transmute(Fsm_spike = mean_FLU, tub = tub)%>%
  left_join(tubs_pm, by = "tub") %>%
  filter(sample_matrix == "sample") %>%
  transmute(Fsm_spike = Fsm_spike,
            Fsm_zero = mean_FLU,
            Fst_spike_pm = Fst_spike_pm,
            Fst_zero_pm = Fst_zero_pm, 
            matrix = 100 * ((Fst_spike_pm - Fst_zero_pm - (Fsm_spike - Fsm_zero))/
                              (Fst_spike_pm - Fst_zero_pm)) ) %>%
  select(matrix)

##Calculating NH4 conc FOR THE PM DATA
tubs_nh4_pm <- cbind(tubs_pm, matrix_pm) %>%
  mutate(Fsm_zero = mean_FLU,
         Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100)),
         int_pm = int_pm, 
         slope_pm = slope_pm,
         NH4_conc = int_pm + slope_pm * Fsm_cor_Taylor) %>%
  filter(sample_matrix != c("matrix")) %>%
  mutate(total_mass = cuke_mass1 +cuke_mass2 +cuke_mass3 + cuke_mass4,
         time = "pm") %>%
  select(- c(int_pm, slope_pm)) %>%
  filter(bottle != "D01")


# stick the am and pm dfs together ------
tubs_nh4_2 <- rbind(tubs_nh4_am, tubs_nh4_pm)%>%
  mutate(date = "two")

#calculate how much pee was added
tubs_controls2 <- tubs_nh4_2 %>%
  filter(cukes_num == "0") %>%
  transmute(NH4_conc_control = mean(NH4_conc), 
            sample_matrix = sample_matrix) %>%
  slice_head()

tubs_nh4_added2 <- tubs_nh4_2 %>%
  left_join(tubs_controls2, by = "sample_matrix") %>%
  mutate(nh4_added = NH4_conc - NH4_conc_control,
         date = "two") %>%
  select(cukes_num, flow_s, mean_flow, total_mass, NH4_conc, nh4_added, date)


# smoosh the two dfs together (week 1 + week 2) -----
tubs_nh4_added_final <- rbind(tubs_nh4_added1, tubs_nh4_added2) %>%
  filter(!(date == "one" & flow_s == 15)) %>% # something weird happened in that tub
  mutate(flow2 = ifelse(cukes_num == 0, "Control", as.character(flow_s)),
         flow2 = factor(flow2, levels = c("Control", 10, 15, 20, 25)),
         cukes = as.factor(ifelse(cukes_num == "0", "Zero", "Four")))

# stats ----
model <- lm(nh4_added ~ flow2 + date, data = tubs_nh4_added_final)
summary(model)
visreg(model)

# model with flow as a continuous
model2 <- glmmTMB(nh4_added ~ mean_flow*cukes + (1|date), data = tubs_nh4_added_final)
summary(model2)
visreg(model2, "mean_flow", by = "cukes")


# Plot ------
csee_pal <- pnw_palette("Starfish")

# palette
pal_flow <- viridis::viridis(10)[1:5]

pal <- viridis::viridis(10)
pal_flow2 <- c(pal[5], pal[8])

# continuous plot
predict_flow <- ggpredict(model2, terms = c("mean_flow", "cukes")) %>% 
  dplyr::rename(mean_flow = x,
                cukes = group)

# use plot pred function
plot_pred(raw_data = tubs_nh4_added_final, 
          predict_data = predict_flow, 
          plot_type = "flow",
          x_var = mean_flow, 
          y_var = nh4_added, 
          lty_var = cukes,
          pch_var = NULL,
          x_axis_lab = NULL,
          pal = pal_flow2,
          theme = "white")

# ggsave("Output/Pub_figs/Supp2Fig3.png", device = "png", height = 9, width = 16, dpi = 400)

# use ggpredict to get estimates
sum_stats <- ggpredict(model, terms = c("flow2")) %>% 
  dplyr::rename(flow2 = x,
                nh4_added = predicted) %>% 
  as_tibble()

# now make dot whisker plots
ggplot() +
  geom_point(data = sum_stats,
             aes(x = flow2, y = nh4_added, colour = flow2),
             size = 6) +
  geom_errorbar(data = sum_stats,
                aes(x = flow2,
                    y = nh4_added,
                    ymin = conf.low,
                    ymax = conf.high, 
                    colour = flow2),
                width = 0.4,
                linewidth = 1.5) +
  geom_jitter(data = tubs_nh4_added_final, 
              aes(x = flow2, y = nh4_added, colour = flow2), 
              size = 3, alpha = 0.5, height=0) +
  theme_white() + 
  labs(x = "Flow rate (cm/s)", y = "Change in Ammonium (umol/L)") +
  theme(legend.position = "none") +
  scale_colour_manual(values = (pal_flow))

#ggsave("Output/Figures/Flow_rates.png", device = "png",
#       height = 9, width = 16, dpi = 400)


# ? ------

# Average for each treatment
mean(tubs_f$nh4_conc[tubs_f$cukes_num == 0]) # 1.999316
mean(tubs_f$nh4_conc[tubs_f$cukes_num == 4]) # 2.342993

(2.342993-1.999316)/1.999316 * 100 # percent increase of 17.18973%
2.342993/1.999316 # cuke tub is 1.171897 x higher

##Code to make the flow conversion line
##June 8, 2021
##By: Emily Leedham

flowline <- read.csv("Data/Emily_flow/2021_06_08_flowline.csv") %>%
  mutate(mean_flow = rowMeans(cbind(mL_flow1, mL_flow2, mL_flow3)))

flowline_2 <- flowline %>%
  mutate(period = as.factor(period))

##Making a graph -----

flow_plot <- ggplot(flowline_2, aes(cm_flow, mean_flow)) +
  geom_point(aes(color = period)) +
  geom_smooth(method = lm)
flow_plot
 
flow_line <- lm(mean_flow ~ cm_flow, data = flowline)
summary(flow_line)
coef(flow_line)

#wednesday R^2 was 0.755 (rows 1 to 5)
#June 14 R^2 was 0.7623 (rows 1 to 10)
##From rows 1 to 15, R^2 was 0.6936
