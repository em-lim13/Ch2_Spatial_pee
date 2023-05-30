# Code to create functions to calculate concentration of ammonium via fluorometric methods
# Using the Taylor protocol 1 + 2 methods

# First load packages
library(tidyverse)
library(patchwork)

# function Nikola built for me!!!
pee_calc <- function(data_path) {
  data_bf <- read_csv(data_path, show_col_types = FALSE) %>%
    filter(sample == "bf") %>% # just the BF sample
    transmute(bf = (FLU1 + FLU2 + FLU3)/3,
              site_code = site_code) 
  
  # Load the standards and the sample data
  data <- read_csv(data_path, show_col_types = FALSE) %>%
    as.data.frame() %>%
    filter(sample != "bf") %>% # cut the bf sample
    left_join(data_bf, by = "site_code") %>% # Add the background FLU as a column
    mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
           mean_FLU = mean_FLU_unadj - bf) # correct the samples with this bf value
  
  # Build the standard curve for protocol 1!
  # We need to estimate how much ammonium is in the standard water
  standard <- data %>%
    filter(sample == "standard") %>% 
    mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
           total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
           nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
  #of NH4 in seawater sample
  
  #model for protocol 1
  mod <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard)
  summary(mod)
  print(summary(mod)$adj.r.squared) # Check Adjusted R-squared value, we're looking for at least 0.98
  
  #graph protocol 2 standard curve
  proto1plot<- (ggplot(standard, aes(nh4_conc_final_umol_L, mean_FLU)) +
                  geom_point() +
                  geom_smooth(method = lm, se = FALSE))
  # Inspect this curve to make sure nothing is wonky!
  
  #save coefficients
  int <- coef(mod)[1]
  slope <- coef(mod)[2]
  
  # calculate ammonium in standards using these coefficients
  standard_b <- standard %>%
    mutate(int = int, 
           slope = slope, 
           nh4_conc = abs(int/slope),
           true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)
  
  # Now move into protocol 2! Use the standards to calculate ammonium in samples
  # BF-corrected FLU of standard curve against calculated NH4 concentration
  mod2 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b)
  summary(mod2)
  
  #graph protocol 2 standard curve
  proto2plot<- (ggplot(standard_b, aes(true_nh4_conc, mean_FLU)) +
          geom_point() +
          geom_smooth(method = lm, se = FALSE))
  # Inspect this curve to make sure nothing is wonky!
  
  # inspect standard curves
  print(proto1plot + proto2plot)

  # pull coefficients from model
  int2 <- coef(mod2)[1]
  slope2 <- coef(mod2)[2]
  
  # calculate how much ammonium is in the samples
  samples <- data %>%
    filter(nh4_vol_uL == "0") %>%
    mutate(nh4_conc = (mean_FLU - int2)/slope2)
}








pee_calc2 <- function(data_path) {
  data_bf <- read_csv(data_path, show_col_types = FALSE) %>%
    filter(sample == "bf") %>% # just the BF sample
    transmute(bf = (FLU1 + FLU2 + FLU3)/3,
              site_code = site_code) 
  
  # Load the standards and the sample data
  data <- read_csv(data_path, show_col_types = FALSE) %>%
    as.data.frame() %>%
    filter(sample != "bf") %>% # cut the bf sample
    left_join(data_bf, by = "site_code") %>% # Add the background FLU as a column
    mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
           mean_FLU = mean_FLU_unadj - bf) # correct the samples with this bf value
  
  # Build the standard curve for protocol 1!
  # We need to estimate how much ammonium is in the standard water
  standard <- data %>%
    filter(sample == "standard") %>% 
    mutate(nh4_added_umol = nh4_vol_uL/1e6 * 2014, #amount of NH4
           total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
           nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
  #of NH4 in seawater sample
  
  #model for protocol 1
  mod <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard)
  summary(mod)
  print(summary(mod)$adj.r.squared) # Check Adjusted R-squared value, we're looking for at least 0.98
  
  #graph protocol 2 standard curve
  proto1plot<- (ggplot(standard, aes(nh4_conc_final_umol_L, mean_FLU)) +
                  geom_point() +
                  geom_smooth(method = lm, se = FALSE))
  # Inspect this curve to make sure nothing is wonky!
  
  #save coefficients
  int <- coef(mod)[1]
  slope <- coef(mod)[2]
  
  # calculate ammonium in standards using these coefficients
  standard_b <- standard %>%
    mutate(int = int, 
           slope = slope, 
           nh4_conc = abs(int/slope),
           true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)
  
  # Now move into protocol 2! Use the standards to calculate ammonium in samples
  # BF-corrected FLU of standard curve against calculated NH4 concentration
  mod2 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b)
  summary(mod2)
  
  #graph protocol 2 standard curve
  proto2plot<- (ggplot(standard_b, aes(true_nh4_conc, mean_FLU)) +
                  geom_point() +
                  geom_smooth(method = lm, se = FALSE))
  # Inspect this curve to make sure nothing is wonky!
  
  # inspect standard curves
  print(proto1plot + proto2plot)
  
  # pull coefficients from model
  int2 <- coef(mod2)[1]
  slope2 <- coef(mod2)[2]
  
  # calculate how much ammonium is in the samples
  samples <- data %>%
    filter(sample != "standard") %>%
    mutate(nh4_conc = (mean_FLU - int2)/slope2)
}
