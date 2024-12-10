# Code for functions
# Em Lim 


# Code to create functions to calculate concentration of ammonium via fluorometric methods
# Using the Taylor protocol 1 + 2 methods

# First load packages
library(tidyverse)
library(patchwork)
library(renv)

# First let's get our packages sorted out using renv

# renv::init()
# this adds three new files and directories:
# renv/library = the project library
# renv.lock = lockfile
# .Rprofile 

# Theme black ----
theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "white"),  
      axis.text.x = element_text(size = base_size*2, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*2, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", linewidth  =  0.2),  
      axis.title.x = element_text(size = base_size*2.5, color = "white", margin = ggplot2::margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size*2.5, color = "white", angle = 90, margin = ggplot2::margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "black",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*2, color = "white"),  
      legend.title = element_text(size = base_size*2.5, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "black"),  
      panel.grid.minor = element_line(color = "black"),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*2.5, color = "white"),  
      strip.text.y = element_text(size = base_size*2.5, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines"),
      plot.tag = element_text(size = 30)
      
    )
  
}


theme_white = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "black"),  
      axis.text.x = element_text(size = base_size*2, color = "black", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*2, color = "black", lineheight = 0.9),  
      axis.ticks = element_line(color = "black", linewidth  =  0.2),  
      axis.title.x = element_text(size = base_size*2.5, color = "black", margin = ggplot2::margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size*2.5, color = "black", angle = 90, margin = ggplot2::margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "white"),  
      legend.key = element_rect(color = "white",  fill = "white"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*1.8, color = "black"),  
      legend.title = element_text(size = base_size*2, face = "bold", hjust = 0, color = "black"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "white", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "black"),  
      panel.grid.major = element_line(color = "white"),  
      panel.grid.minor = element_line(color = "white"),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*2.5, color = "black"),  
      strip.text.y = element_text(size = base_size*2.5, color = "black",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "white", fill = "white"),  
      plot.title = element_text(size = base_size*1.2, color = "black"),  
      plot.margin = unit(rep(1, 4), "lines"),
      # changes the size of the patchwork annotations
      plot.tag = element_text(size = 30)
      
    )
  
}

# Publication figure theme ------
pub_theme = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "black", linewidth = 0.1),  
      axis.text.x = element_text(size = base_size*0.9, color = "black", lineheight = 0.2),  
      axis.text.y = element_text(size = base_size*0.9, color = "black", lineheight = 0.2),  
      axis.ticks = element_line(color = "black", linewidth  =  0.2),  
      axis.title.x = element_text(size = base_size*1, color = "black", margin = ggplot2::margin(0, 1, 0, 0)),  
      axis.title.y = element_text(size = base_size*1, color = "black", angle = 90, margin = ggplot2::margin(0, 1, 0, 0)),  
      axis.ticks.length = unit(0.1, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "white"),  
      legend.key = element_rect(color = "white",  fill = "white"),  
      legend.key.size = unit(0.6, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*1, color = "black"),  
      legend.title = element_text(size = base_size*1, face = "bold", hjust = 0, color = "black"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      # Specify panel options
      panel.background = element_rect(fill = "white", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),  
      panel.grid.major = element_line(color = "white"),  
      panel.grid.minor = element_line(color = "white"),  
      panel.spacing = unit(0.1, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*1, color = "black"),  
      strip.text.y = element_text(size = base_size*1, color = "black",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "white", fill = "white"),  
      plot.title = element_text(size = base_size*1, color = "black"),  
      plot.margin = unit(rep(0.5, 4), "lines"),
      # changes the size of the patchwork annotations
      plot.tag = element_text(size = 13)
      
    )
  
}

# Calculate nh4+ -----
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

# Clean up messy taxonomy -----
clean_sp_names <- function(datafile) {
  new_data <- datafile %>%
    mutate(
      # species
      species_name = case_when(
        species_name == "Montereina nobilis" ~ "Peltodoris nobilis",
        species_name == "Parastichopus californicus" ~ "Apostichopus californicus",
        species_name == "Berthella californica" ~ "Berthella chacei",
        species_name == "Henricia leviuscula" ~ "Henricia spp.",
        TRUE ~ as.character(species_name)))
}

# then clean phylo
clean_phylo_names <- function(datafile){
  new_data <- datafile %>%
    clean_sp_names() %>%
    mutate(
      # Phylum
      phylum = case_when(species_name == "Artedius lateralis" ~ "Chordata",
                         species_name == "Ascelichthys rhodorus" ~ "Chordata",
                         species_name == "Clupea pallasii" ~ "Chordata",
                         species_name == "Sebastes paucispinis" ~ "Chordata",
                         species_name == "Orthonopias triacis" ~ "Chordata",
                         species_name == "Blepsias cirrhosus" ~ "Chordata",
                         species_name == "Gobiesox maeandricus" ~ "Chordata",
                         species_name == "Leptocottus armatus" ~ "Chordata",
                         species_name == "Rimicola muscarum" ~ "Chordata",
                         species_name == "Triopha spp." ~ "Mollusca",
                         species_name == "Melibe leonina" ~ "Mollusca",
                         species_name == "Octopus rubescens" ~ "Mollusca",
                         species_name == "Heptacarpus stylus" ~ "Arthropoda",
                         species_name == "Cryptolithodes sitchensis" ~ "Arthropoda",
                         species_name == "Pandulus spp." ~ "Arthropoda",
                         TRUE ~ as.character(phylum)),
      # class
      class = case_when(class == "Actinopteri" | class == "Teleostei" ~ "Actinopterygii",
                        species_name == "Artedius lateralis" ~ "Actinopterygii",
                        species_name == "Ascelichthys rhodorus" ~ "Actinopterygii",
                        species_name == "Clupea pallasii" ~ "Actinopterygii",
                        species_name == "Sebastes paucispinis" ~ "Actinopterygii",
                        species_name == "Orthonopias triacis" ~ "Actinopterygii",
                        species_name == "Blepsias cirrhosus" ~ "Actinopterygii",
                        species_name == "Leptocottus armatus" ~ "Actinopterygii",
                        species_name == "Gobiesox maeandricus" ~ "Actinopterygii",
                        species_name == "Rimicola muscarum" ~ "Actinopterygii",
                        species_name == "Triopha spp." ~ "Gastropoda",
                        species_name == "Melibe leonina" ~ "Gastropoda",
                        species_name == "Octopus rubescens" ~ "Cephalopoda",
                        species_name == "Heptacarpus stylus" ~ "Malacostraca",
                        species_name == "Cryptolithodes sitchensis" ~ "Malacostraca",
                        species_name == "Pandulus spp." ~ "Malacostraca",
                        TRUE ~ as.character(class)),
      # order
      order = case_when(species_name == "Opalia wroblewskyi" ~ "Caenogastropoda", 
                        species_name == "Apostichopus californicus" ~ "Synallactida",
                        family == "Fissurellidae" ~ "Lepetellida",
                        family == "Turbinidae" | 
                          family == "Calliostomatidae" | 
                          family == "Tegulidae" ~ "Trochida",
                        order == "Ovalentaria incertae sedis" ~ "Perciformes",
                        order == "Scorpaeniformes" ~ "Perciformes",
                        order == "Gasterosteiformes" ~ "Perciformes",
                        species_name == "Leptocottus armatus" ~ "Perciformes",
                        species_name == "Gobiesox maeandricus" ~ "Gobiesociformes",
                        species_name == "Rimicola muscarum" ~ "Gobiesociformes",
                        species_name == "Unidentified shrimp" ~ "Decapoda",
                        family == "Pectinidae" ~ "Pectinida",
                        family == "Haliotidae" ~ "Lepetellida",
                        family == "Acmaeidae" ~ "Patellogastropoda", # no order so going to subclass
                        family == "Lottiidae" ~ "Patellogastropoda", # no order so going to subclass
                        species_name == "Artedius lateralis" ~ "Perciformes",
                        species_name == "Ascelichthys rhodorus" ~ "Perciformes",
                        species_name == "Clupea pallasii" ~ "Clupeiformes",
                        species_name == "Sebastes paucispinis" ~ "Perciformes",
                        species_name == "Orthonopias triacis" ~ "Perciformes",
                        species_name == "Blepsias cirrhosus" ~ "Perciformes",
                        species_name == "Triopha spp." ~ "Nudibranchia",
                        species_name == "Melibe leonina" ~ "Nudibranchia",
                        species_name == "Octopus rubescens" ~ "Octopoda",
                        species_name == "Heptacarpus stylus" ~ "Decapoda",
                        species_name == "Cryptolithodes sitchensis" ~ "Decapoda",
                        species_name == "Pandulus spp." ~ "Decapoda",
                        TRUE ~ as.character(order)),
      # family
      family = case_when(species_name == "Paguroidea spp." ~ "Paguridae", 
                         family == "Pycnopodiidae" ~ "Asteriidae",
                         family == "Hapalogastridae" ~ "Lithodidae",
                         species_name == "Hermissenda crassicornis" ~ "Myrrhinidae",
                         species_name == "Artedius lateralis" ~ "Cottidae",
                         species_name == "Ascelichthys rhodorus" ~ "Cottidae",
                         species_name == "Clupea pallasii" ~ "Clupeidae",
                         species_name == "Sebastes paucispinis" ~ "Sebastidae",
                         species_name == "Orthonopias triacis" ~ "Cottidae",
                         species_name == "Blepsias cirrhosus" ~ "Hemitripteridae",
                         species_name == "Leptocottus armatus" ~ "Cottidae",
                         species_name == "Gobiesox maeandricus" ~ "Gobiesocidae",
                         species_name == "Rimicola muscarum" ~ "Gobiesocidae",
                         species_name == "Triopha spp." ~ "Polyceridae",
                         species_name == "Melibe leonina" ~ "Tethydidae",
                         species_name == "Octopus rubescens" ~ "Octopodidae",
                         species_name == "Heptacarpus stylus" ~ "Thoridae",
                         species_name == "Pandulus spp." ~ "Pandalidae",
                         species_name == "Cryptolithodes sitchensis" ~ "Lithodidae",
                         TRUE ~ as.character(family)) ) %>%
    # Cut the non-animal species
    filter(species_name != "Debris - Metal") %>%
    filter(species_name != "Debris - Other") %>%
    filter(species_name != "Debris - Wood") %>%
    filter(species_name != "Debris - Glass") %>%
    filter(species_name != "Debris - Fishing Gear") %>%
    filter(species_name != "Debris - zero")  
}

# Length to weight ----
# all coeffs for fish are from fishbase!
# I confirmed the formula and the units (cm and g)

length_to_weight <- function(datafile) {
  new_data <- datafile %>%
    mutate(
      # then convert length to weight
      weight_per_indiv_g = case_when(
        
        # Blennies
        family == "Clinidae" ~ 0.00513*size_class^3.06,
        
        # Gobies
        species_name == "Rhinogobiops nicholsii" ~ 0.01047*size_class^3.03,
        family == "Gobiidae" ~ 0.01047*size_class^3.03,
        
        # Greenlings
        species_name == "Hexagrammos decagrammus" ~ 0.00813*size_class^3.13,
        species_name == "Hexagrammos stelleri" ~ 0.00692*size_class^3.16,
        species_name == "Oxylebius pictus" ~ 0.01122*size_class^3.04,
        species_name == "Ophiodon elongatus" ~ 0.00389*size_class^3.12,
        species_name == "Hexagrammos spp." ~ 0.00813*size_class^3.13,
        family == "Hexagrammidae" ~ 0.00813*size_class^3.13,
        
        # Rockfish
        family == "Sebastidae" ~ 0.01000*size_class^3.09,
        
        # Sculpins
        species_name == "Jordania zonope" ~ 0.00389*size_class^3.12,
        species_name == "Artedius harringtoni" ~ 0.00631*size_class^3.15,
        species_name == "Artedius lateralis" ~ 0.00631*size_class^3.15,
        species_name == "Artedius fenestralis" ~ 0.00631*size_class^3.15,
        species_name == "Hemilepidotus hemilepidotus" ~ 0.00631*size_class^3.15,
        species_name == "Cottidae spp." ~ 0.00631*size_class^3.15,
        species_name == "Enophrys bison" ~ 0.00794*size_class^3.13,
        species_name == "Rhamphocottus richardsonii" ~ 0.01995*size_class^3.01,
        species_name == "Scorpaenichthys marmoratus" ~ 0.00389*size_class^3.12,
        species_name == "Oligocottus maculosus" ~ 0.00631*size_class^3.15,
        species_name == "Leptocottus armatus" ~ 0.01096*size_class^3.19,
        species_name == "Blepsias cirrhosus" ~ 0.00631*size_class^3.14,
        species_name == "Myoxocephalus polyacanthocephalus" ~ 0.00832*size_class^3.14,
        species_name == "Myoxocephalus aenaeus" ~ 0.00832*size_class^3.14,
        species_name == "Asemichthys taylori" ~ 0.00631*size_class^3.15,
        species_name == "Ascelichthys rhodorus" ~ 0.00676*size_class^3.17,
        species_name == "Orthonopias triacis" ~ 0.01000*size_class^3.04,
        
        #Perch
        species_name == "Brachyistius frenatus" ~ 0.01318*size_class^3.05,
        species_name == "Embiotoca lateralis" ~ 0.01950*size_class^2.97,
        species_name == "Rhacochilus vacca" ~ 0.01950*size_class^2.97,
        species_name == "Cymatogaster aggregata" ~ 0.01950*size_class^2.97,
        species_name == "Embiotocidae spp." ~ 0.01950*size_class^2.97,
        species_name == "Percidae spp." ~ 0.01950*size_class^2.97,
        family == "Embiotocidae" ~ 0.01950*size_class^2.97,
        
        # Gunnels + gunnel-like fish
        species_name == "Anarrhichthys ocellatus" ~ 0.00398*size_class^3.17,
        species_name == "Apodichthys flavidus" ~ 0.00102*size_class^3.06,
        species_name == "Pholis ornata" ~ 0.00162*size_class^3.19,
        species_name == "Pholis laeta" ~ 0.00162*size_class^3.19,
        species_name == "Pholis clemensi" ~ 0.00162*size_class^3.19,
        species_name == "Pholis spp." ~ 0.00162*size_class^3.19,
        species_name == "Lumpenus sagitta" ~ 0.00129*size_class^2.99,
        species_name == "Chirolophis nugator" ~ 0.00372*size_class^3.16,
        
        #Misc
        species_name == "Liparis florae" ~ 0.00525*size_class^3.15, #snailfish
        species_name == "Gobiesox maeandricus" ~ 0.00617*size_class^3.15, #clingfish
        species_name == "Aulorhynchus flavidus" ~ 0.00263*size_class^3.14, #tubesnout
        species_name == "Syngnathus leptorhynchus" ~ 0.00028*size_class^3.18, #pipefish
        species_name == "Clupea pallasii" ~ 0.00603*size_class^3.13, #herring
        species_name == "Gasterosteus aculeatus" ~ 0.00977*size_class^3.09, #stickleback
        species_name == "Porichthys notatus" ~ 0.00562*size_class^3.16, #plainfin
        species_name == "Gibbonsia metzi" ~ 0.00513*size_class^3.06,
        species_name == "Citharichthys stigmaeus" ~ 0.00759*size_class^3.15,
        species_name == "Rimicola muscarum" ~ 0.00617*size_class^3.15, #kelp clingfish
        
        # NOW DO THE INVERTS
        
        # ARTHROPODA
        # Cancridae family
        species_name == "Cancer productus" ~ 200, # mean wet tissue from my data
        family == "Cancridae" ~ 3, # Hines 1982 small crabs avg
        # Epialtidae family
        species_name == "Pugettia producta" ~ 46, # Hines 1982
        species_name == "Scyra acutifrons" ~ 2, # Hines 1982
        family == "Epialtidae" ~ 1.235, # Hines, for Pugettia richii
        # Lithodidae
        species_name == "Cryptolithodes sitchensis" ~ 3, # Hines 1982 small crabs avg
        species_name == "Cryptolithodes typicus" ~ 3, # Hines 1982 small crabs avg
        family == "Lithodidae" ~ 65, # Stewart et al 2015 for Paralithodes rathbuni and Phyllolithodes papillosus
        
        # Oregoniidae family
        family == "Oregoniidae" ~ 3, # Hines 1982 small crabs avg
        # Panopeidae
        family == "Panopeidae" ~ 3, # Hines 1982 small crabs avg
        # Paguridae family
        family == "Paguridae" ~ 0.43, # McKinney et al 2004
        # Pandalidae 
        family == "Pandalidae" ~ 0.11, # McKinney et al 2004 for Palaemonetes pugio 
        species_name == "Heptacarpus stylus" ~ 0.11, # McKinney et al 2004 for Palaemonetes pugio 
        species_name == "Unidentified shrimp" ~ 0.11, # McKinney et al 2004 for Palaemonetes pugio
        # Porcellanidae family
        family == "Porcellanidae" ~ 4.25, # Stillman and Somero 1996 for Petrolisthes
        # Misc decapod
        species_name == "Brachyura spp." ~ 3, # Hines 1982 small crabs avg
        
        # CNIDARIA
        phylum == "Cnidaria" ~ 0.01, # Båmstedt 2015
        
        # CTENOPHORA
        phylum == "Ctenophora" ~ 0.01, # Båmstedt 2015
        
        # ECHINODERMATA
        # Asteroidea class
        # Asteriidae 
        species_name == "Evasterias troschelii" ~ 66.5, # lower end of O'Clair 1985 range
        species_name == "Leptasterias hexactis" ~ 5.5, # Menge 1975
        species_name == "Orthasterias koehleri" ~ 66.5, # lower end of O'Clair 1985 range
        species_name == "Pisaster ochraceus" ~ 128, # Sanford 2002
        species_name == "Pisaster brevispinus" ~ 146.18, # Peters et al 2019 for Pisaster giganteus
        species_name == "Stylasterias forreri" ~ 66.5, # lower end of O'Clair 1985 range
        species_name == "Pycnopodia helianthoides" ~ exp(-3.9989)*size_class^3.133, # Lee 2016
        family == "Asteriidae" ~ 5.5, # Menge 1975
        # Asterinidae
        species_name == "Patiria miniata" ~ 26.97, # 6.83 g dry weight Peters et al 2019
        family == "Asterinidae" ~ 26.97, # 6.83 g dry weight Peters
        # Echinasteridae
        family == "Echinasteridae" ~ 10, # Menge 1975
        # Asteropseidae
        family == "Asteropseidae" ~ 92, # Montgomery 2014
        # Goniasteridae
        family == "Goniasteridae" ~ 10, # Menge 1975 (Henricia)
        # Solasteridae
        family == "Solasteridae" ~ 486, # Montgomery 2014 for Solaster stimpsoni
        # Pterasteridae
        family == "Pterasteridae" ~ 10, # Menge 1975 (Henricia)
        
        # Echinoidea class
        species_name == "Mesocentrotus franciscanus" ~ 29.51, # Schuster and Bates 2023, tissue  # 163.41 whole wet weight Peters et al 2019
        family == "Strongylocentrotidae" ~  20, #Stewart et al 2015 tissue Strongylocentrotus polyacanthus
        # 50 grams from Siikavuopio
        
        # Holothuroidea
        species_name == "Apostichopus californicus" ~ 319.31, # Peters et al 2019 for Apostichopus parvimensis
        # 829 was the cuke avg from my measurements but I'm worried that's just a ton of water weight
        
        # MOLLUSCA
        # Bivalves
        species_name == "Crassadoma gigantea" ~ exp(-3.25924)*size_class^2.39442, #MacDonald 1991
        family == "Pectinidae" ~ 2.5, #MacDonald 1991 rough Chlamys avg
        # Cephalopoda
        species_name == "Enteroctopus dofleini" ~ 137.5, # # Osborn 1995 thesis, upper limit for Octopus rubescens
        species_name == "Octopus rubescens" ~ 80, # Osborn 1995 thesis
        
        # Gastropoda
        species_name == "Pomaulax gibberosus" ~ 31, # Schuster and Bates 2023
        species_name == "Haliotis kamtschatkana" ~ 0.0000578*(size_class*10)^3.2, # Zhang 2007
        order == "Nudibranchia" ~  0.54, # McKinney et al 2004 gastropods
        species_name == "Triopha spp." ~  0.54, # McKinney et al 2004 gastropods
        class == "Gastropoda" ~ 0.91, # Palmer 1982 for Nucella sp
        # Flatworm
        species_name == "Eurylepta leoparda" ~ 0.54, # McKinney et al 2004 gastropods
        
        # Everything else
        TRUE ~ as.numeric(0.5)),
      # set the really big wolf eel weight manually to largest record weight
      # otherwise the calc thinks it's MASSSSSIVE
      weight_per_indiv_g = if_else(size_class == 187.5, 18400, weight_per_indiv_g),
      weight_per_indiv_kg = weight_per_indiv_g/1000,
      weight_size_class_sum = (weight_per_indiv_kg*survey_den)) 
}



# Invert length to weight ------

# Mean cuke wet weight = 829 g
# Mean cuke dry weight = 39

# Calculate scallop weight from length from MacDonald et al 1991
# scallops length to weight
# scallop <- read_csv("Data/Size_data/scallop_l_w.csv") %>%
#  mutate(shell_height_cm = shell_height*0.1,
#         log_weight = log(body_weight),
#         log_shell_height = log(shell_height_cm))
#
# mod <- lm(log_weight ~ log_shell_height, data = scallop)
# summary(mod)
# intercept = log(a) = -3.25924
# slope = b = 2.39442
# mass = exp(-3.25924)*size_class^2.39442

# Calculate mean wet weight from Peters et al 2019
# peters <- read_csv("Data/Size_data/peters_invertebrate_NH4_excretions.csv") %>%
#   unite("species", TAXON_GENUS:TAXON_SPECIES, remove = FALSE)
# 
# peters_avg <- peters %>%
#   group_by(species) %>%
#   summarise(mean_weight = mean(WM_g)) %>%
#   mutate(mean_weight = round(mean_weight, 2))


# Calculate means from Jasmin's data
# schuster <- read_csv("Data/Size_data/Schuster_Grazer_Weights.csv") %>%
#   mutate(shell_g = DryWeight_g - AshedWeight_g,
#     flesh_g = WetWeight_g - shell_g) %>%
#   group_by(Species) %>%
#   summarise(avg_flesh = mean(flesh_g, na.rm=TRUE))

# Calculate diameter to weight from Montgomery 2014
# pycno <- read_csv("Data/Size_data/Montgomery_pycno.csv") %>%
#   mutate(mass_g = 10^log_wet_mass_g,
#          arm_length_mm = 10^log_arm_length_mm,
#          size_class_cm = arm_length_mm*2/10,
#          log_size = log(size_class_cm),
#          log_mass = log(mass_g))
# 
# mod <- lm(log_mass ~ log_size, data = pycno)
# summary(mod)
# intercept = log(a) = -4.368530
# slope = b = 3.301604
# mass = exp(-4.368530)*size_class^3.301604

# Evasterias and pisaster from Sharon Kay
# kay <- read_csv("Data/Size_data/Kay_TableS3_FigS3_data.csv")
# kay_avg <- kay
#   group_by(species) %>%
#   summarise(mean_start = mean(start_weight),
#             mean_end = mean(end_weight))
# pisaster <- kay %>% filter(species == "Evasterias")
# hist(pisaster$start_weight)



# Fish home ranges ----
home_range <- function(datafile){
  new_data <- datafile %>%
    mutate(range = case_when(
      # Gobies
      species_name == "Rhinogobiops nicholsii" ~ 0.05, # guessed from painted greenling
      # Greenlings
      species_name == "Hexagrammos decagrammus" ~ 0.3, # 0.1 - 0.5
      species_name == "Hexagrammos stelleri" ~ 0.3, # from kelp greenling 0.1 - 0.5
      species_name == "Oxylebius pictus" ~ 0.02,
      species_name == "Ophiodon elongatus" ~ 28.3, # 3.3 - 498
      species_name == "Hexagrammos spp." ~ 0.3, # from kelp greenling 0.1 - 0.5
      
      # Rockfish
      species_name == "Sebastes melanops" ~ 0.5,
      species_name == "Sebastes caurinus" ~ 0.15,
      species_name == "Sebastes flavidus" ~ 500, # 140-1400 range
      species_name == "Sebastes maliger" ~ 0.05,
      species_name == "Sebastes nebulosus" ~ 0.01,
      species_name == "Sebastes spp." ~ 0.18, # avg of small range rockfish
      species_name == "Sebastes spp. juv" ~ 0.18, # avg of small range rockfish
      species_name == "Sebastes pinniger" ~ 700, # up to 700
      
      # Sculpins, all guesses based on cabezon + red irish
      species_name == "Jordania zonope" ~ 0.05,
      species_name == "Artedius harringtoni" ~ 0.05,
      species_name == "Artedius lateralis" ~ 0.05,
      species_name == "Artedius fenestralis" ~ 0.05,
      species_name == "Hemilepidotus hemilepidotus" ~ 0.05,# actual estimate
      species_name == "Cottidae spp." ~ 0.05,
      species_name == "Enophrys bison" ~ 0.05,
      species_name == "Rhamphocottus richardsonii" ~ 0.05,
      species_name == "Scorpaenichthys marmoratus" ~ 0.05,# actual estimate
      species_name == "Oligocottus maculosus" ~ 0.05,
      species_name == "Leptocottus armatus" ~ 0.05,
      species_name == "Blepsias cirrhosus" ~0.05,
      species_name == "Myoxocephalus polyacanthocephalus" ~ 0.05,
      species_name == "Myoxocephalus aenaeus" ~ 0.05,
      species_name == "Asemichthys taylori" ~ 0.05,
      
      #Perch
      species_name == "Embiotoca lateralis" ~ 1, # inferred
      species_name == "Rhacochilus vacca" ~ 1, # unknown, inferred from other perch
      species_name == "Brachyistius frenatus" ~ 1, # inferred
      species_name == "Cymatogaster aggregata" ~ 1, # unknown, inferred from other perch
      species_name == "Embiotocidae spp." ~ 1, # unknown, inferred from other perch
      species_name == "Percidae spp." ~ 1, # unknown, inferred from other perch
      
      # Gunnels + gunnel-like fish
      species_name == "Anarrhichthys ocellatus" ~ 0.05,
      species_name == "Apodichthys flavidus" ~ 0.05, # inferred from wolf eel
      species_name == "Pholis ornata" ~ 0.05, # inferred from wolf eel
      species_name == "Pholis laeta" ~ 0.05, # inferred from wolf eel
      species_name == "Pholis clemensi" ~ 0.05, # inferred from wolf eel
      species_name == "Pholis spp." ~ 0.05, # inferred from wolf eel
      species_name == "Lumpenus sagitta" ~ 0.05, # inferred from wolf eel
      species_name == "Chirolophis nugator" ~ 0.05, # inferred from wolf eel
      
      #Misc
      species_name == "Liparis florae" ~ 0.05, # inferred from wolf eel
      species_name == "Aulorhynchus flavidus" ~ 1, # inferred from perch
      species_name == "Syngnathus leptorhynchus" ~ 1, # inferred from perch
      species_name == "Clupea pallasii" ~ 525, # 50 - 1000
      species_name == "Gasterosteus aculeatus" ~ 1, # inferred from perch
      species_name == "Porichthys notatus" ~ 300, # midshipman migrate at night 0 - 366 m deep
      species_name == "Gibbonsia metzi" ~ 1, # inferred from perch
      species_name == "Citharichthys stigmaeus" ~ 100, # flatfish seem to move a lot, avg from them
      TRUE ~ as.numeric(0.01)), # to give all the inverts a range weight of 1
      range_weight = 1/(range*100), # scale smallest range to a weight of 1, everything else is fraction
      weight_weighted = weight_size_class_sum*range_weight # weight weights by range
    )
  
}



# Joining surveys up by depth -----

# this is a function to only keep RLS data that corresponds to the survey that I took a nh4 measurement on!

# tricky ones:
### BMSC1 2021: two survey depths, one nh4 sample 
### BMSC6 2022: two nh4 samples at 2 depths
### BMSC1 2022: two survey depths, one nh4 sample 
### BMSC5 2022: two surveys at same depth, one nh4 sample
# I was on the later shallower survey
### BMSC6 2022: two surveys at same depth, two nh4 samples!
### BMSC11 2022: two survey depths, one nh4 sample 
### BMSC12 2022: two survey depths, one nh4 sample 
### BMSC11 2023: two survey depths, one nh4 sample 
# I tried to take nh4 samples between the shallower and deeper surveys
### BMSC12 2023: two survey depths, one nh4 sample 
# I tried to take nh4 samples between the shallower and deeper surveys
# BMSC24 2023: two survey depths, one nh4 sample
# I tried to take nh4 samples between the shallower and deeper surveys
# BMSC25 2023: two survey depths, one nh4 sample
# these two were back to back, don't average
# BMSC26 2023: two survey depths, one nh4  (only 2 pee reps tho)
# these two were back to back, don't average
# I was with the 10 m team that got in first
# BMSC27 2023: two survey depths, one nh4 sample
# also back to back, don't average
# BMSC1 2023: two survey depths, one nh4 sample
# BMSC5 2023: two survey depths, one nh4 sample
# BMSC6 2023: two survey depths, one nh4 sample
# the shallower team was way shallower
# BMSC8 2023: two survey depths, one nh4 sample

# So I need to decide what to do about repeated surveys
# Is the nh4 sample I took specific to the transect I took it on?
# Which means I should cut the RLS surveys that aren't the transects where the pee samples were taken
# Or would I want to average the two transects and take the mean biomass from the two to relate to the overall pee sample

# For simplicity I think I just want to keep the RLS survey from the transect where the pee is from

depth_function <- function(datafile){
  new_data <- datafile %>%
    mutate(correct = case_when(
      site_code== "BMSC6" &year =="2022" &depth == "8.5" &survey_depth== "6.5" ~ "no",
      site_code== "BMSC6" &year =="2022" &depth == "5.5" &survey_depth== "9" ~ "no",
      site_code== "BMSC6" &year =="2022" &depth == "6" &survey_depth== "9" ~ "no",
      site_code== "BMSC1" &year == "2021" & survey_depth == "6.5" ~ "no",
      site_code== "BMSC1" &year == "2022" & survey_depth == "4.7" ~ "no",
      site_code== "BMSC5" &year == "2022" & survey_depth == "6.2" ~ "no",
      site_code== "BMSC11" &year == "2022" & survey_depth == "8.5" ~ "no",
      site_code== "BMSC12" &year == "2022" & survey_depth == "9" ~ "no",
      site_code== "BMSC11" &year == "2023" & survey_depth == "5.5" ~ "no",
      site_code== "BMSC12" &year == "2023" & survey_depth == "6.5" ~ "no",
      site_code== "BMSC24" &year == "2023" & survey_depth == "7.5" ~ "no",
      site_code== "BMSC25" &year == "2023" & survey_depth == "5.5" ~ "no",
      site_code== "BMSC26" &year == "2023" & survey_depth == "9.5" ~ "no",
      site_code== "BMSC27" &year == "2023" & survey_depth == "7" ~ "no",
      site_code== "BMSC1" &year == "2023" & survey_depth == "8" ~ "no",
      site_code== "BMSC5" &year == "2023" & survey_depth == "8" ~ "no",
      site_code== "BMSC6" &year == "2023" & survey_depth == "5.5" ~ "no",
      site_code== "BMSC8" &year == "2023" & survey_depth == "7.5" ~ "no",
      TRUE ~ as.character("yes")
      # cut out the rls transects I didn't directly measure nh4 on 
    )) %>%
    filter(correct == "yes") %>%
    select(-c(correct, hour))
}


# Mapping -----
map_daddy <- function(lat_min, lat_max, long_min, long_max, 
                      coord_data, nh4_var, kelp_var, point_size, map_file, 
                      white_background = TRUE, invert = FALSE) {
  
  # invert land and sea colours if I use the potato map
  
  if(invert == FALSE){
    sea <- blue
    land <- "white"
  }
  if(invert == TRUE){
    sea <- "white"
    land <- blue
  }
  
  # specify black or white background
  if(white_background == TRUE){
    background <- "white"
    features <- "black"
  }
  if(white_background == FALSE){
    background <- "black"
    features <- "white"
  }
  
  ggplot() +
    geom_sf(data = map_file, fill = land, colour = sea) +
    # add points
    geom_sf(data = coord_data, 
            colour = "black",
            alpha = 0.9,
            size = point_size,
            aes(fill = {{nh4_var}},
                pch = {{kelp_var}})) +
    viridis::scale_fill_viridis(option="magma", direction = -1,
                                limits = c(0, 2),
                                guide = guide_colorbar(frame.colour = features, ticks.colour = features)) +
    coord_sf(xlim = c(lat_min, lat_max), ylim = c(long_min, long_max), expand = FALSE)  +
    labs(fill = expression(paste("NH"[4]^" +",(mu*M)))) +
    scale_shape_manual(values = c(21, 25), drop = F) +
    guides(pch = guide_legend(override.aes = 
                                list(colour = features))) +
    # Themes
    theme_bw() +
    theme(
      # panel stuff
      panel.background = element_rect(fill = sea),
      panel.grid.major = element_line(color = sea),
      panel.border = element_rect(fill = NA, colour = features),
      # remove axis
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      plot.title = NULL,
      plot.margin=grid::unit(c(0,0,0,0), "mm"),
      # Specify legend options
      legend.background = element_rect(color = NA, fill = background),  
      legend.key = element_rect(color = background,  fill = background),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = 22, color = features),  
      legend.title = element_text(size = 22, face = "bold", hjust = 0, color = features),  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      # if you want the legend back on the side unhash these and hash the inside bit
#      legend.box = NULL,
#      legend.position = "right", 
      
      # try to get legend boxes inside plot
      legend.position = "inside",
      legend.position.inside = c(0.99, 0.15),
      legend.justification.inside = c(0.99, 0.15),
      legend.box.background = element_rect(color = features, linewidth = 0),
      
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = 30, color = features),  
      strip.text.y = element_text(size = 30, color = features,angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = background, fill = background),  
    ) +
    annotation_scale(location = "br", width_hint = 0.4, text_cex = 1.75) +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
}


map_daddy_np <- function(lat_min, lat_max, long_min, long_max, 
                         map_file, invert) {
  
  if(invert == FALSE){
    sea <- blue
    land <- "white"
  }
  if(invert == TRUE){
    sea <- "white"
    land <- blue
  }
  
  ggplot() +
    geom_sf(data = map_file, fill = land, colour = sea) +
    coord_sf(xlim = c(lat_min, lat_max), ylim = c(long_min, long_max), expand = FALSE)  +
    # Themes
    theme_bw() +
    theme(
      # panel stuff
      panel.background = element_rect(fill = sea),
      panel.grid.major = element_line(color = sea),
      panel.border = element_rect(fill = NA, colour = "black"),
      # remove axis
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      plot.title = NULL,
      plot.margin=grid::unit(c(0,0,0,0), "mm"),
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "black",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = 20, color = "white"),  
      legend.title = element_text(size = 20, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL,
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = 30, color = "white"),  
      strip.text.y = element_text(size = 30, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
    ) +
    annotation_scale(location = "br", width_hint = 0.4) +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
}



# Plotting -----

# Annotation -----
place_label <- function(label, size = 4.5, ...) {
  annotate("text", label = label, x = -Inf, y = Inf, 
           vjust = 1.4, hjust = -0.15, size = size, ...)
}

# Dot Whisker Plot -----
dot_whisker <- function(sum_data, all_data, x_var, y_var, pch_var = NULL, 
                        labels, pal, theme ="white"){
  
  if(theme == "white"){
    theme <- pub_theme()
    features <- "black"
    dot_var = 3
    jitter_var = 2
    line_var = 0.75
  } else {
    theme <- theme_black()
    features <- "white"
    dot_var = 8
    jitter_var = 5
    line_var = 1.5
  }
  
  ggplot() +
    geom_point(data = {{sum_data}},
               aes(x = {{x_var}}, y = {{y_var}}, colour = {{x_var}}, pch = {{pch_var}}),
               size = dot_var) +
    geom_errorbar(data = {{sum_data}},
                  aes(x = {{x_var}},
                      y = {{y_var}},
                      ymin = conf.low,
                      ymax = conf.high, 
                      colour = {{x_var}}),
                  width = 0.4,
                  linewidth = line_var) +
    geom_jitter(data = {{all_data}}, 
                aes(x = {{x_var}}, y = {{y_var}}, colour = {{x_var}}, pch = {{pch_var}}), 
                size = jitter_var, alpha = 0.4, height=0, width = 0.2) +
    theme + 
    theme(legend.position = "none",
          plot.title = element_text(size = 30)) +
    scale_colour_manual(values = rev(pal)) +
    labs(y = expression(paste("Ammonium"~(mu*M))), x = " ") +
    scale_x_discrete(labels = {{labels}}) +
    theme(axis.text.x = ggtext::element_markdown())
}


# Coeff plots -----
coeff_plot <- function(coeff_df, pal, theme = "white"){
  
  if(theme == "white"){
    theme <- pub_theme()
    features <- "black"
    size_var = 5
    line_var = 1
  } else {
    theme <- theme_black()
    features <- "white"
    size_var = 10
    line_var = 3
  }
  
  ggplot(coeff_df, aes(x = estimate, y = variable, 
                       xmin = lower_CI, xmax = upper_CI, 
                       colour = variable)) +
    geom_point(size = size_var) +
    geom_errorbar(width = 0, linewidth = line_var) +
    geom_vline(xintercept = 0, color = features, linetype = "dashed", linewidth = line_var*0.25) +
    labs(x = "Coefficient", y = " ") +
    scale_y_discrete(limits = rev(levels(coeff_df$variable))) +
    theme +
    theme(legend.position = "none") + 
    scale_colour_manual(values = pal)
}

# Plot model predictions ----

plot_pred <- function(raw_data, predict_data, 
                      plot_type,
                      x_var, y_var, 
                      lty_var = NULL,
                      pch_var = NULL,
                      x_axis_lab = NULL,
                      pal,
                      theme = "white"){
  if(theme == "white"){
    theme <- pub_theme()
    features <- "black"
    size_var = 2
    line_var = 1
  } else {
    theme <- theme_black()
    features <- "white"
    size_var = 4
    line_var = 2
  }
  
  base_pred_plot <-  ggplot() + 
    geom_point(data = raw_data, 
               aes(x = {{x_var}}, y = {{y_var}}, 
                   colour = {{lty_var}}, fill = {{lty_var}},
                   pch = {{pch_var}}), 
               alpha = 0.5, size = size_var) +
    geom_line(data = predict_data,
              aes(x = {{x_var}}, y = predicted,
                  colour = {{lty_var}}),
              linewidth = line_var) +
    geom_ribbon(data = predict_data,
                aes(x = {{x_var}}, y = predicted, fill = {{lty_var}},
                    ymin = conf.low, ymax = conf.high), 
                alpha = 0.15) +
    labs(colour = "Tide", fill = "Tide", pch = "Tide") +
    theme +
    scale_colour_manual(values = (pal)) +
    scale_fill_manual(values = (pal)) +
    guides(lty = guide_legend(override.aes = list(linewidth = line_var/3)),
           size = guide_legend(override.aes = list(colour = features)),
           colour = guide_legend(override.aes = list(size = size_var*0.75, linewidth = line_var/3)))
  
  # then add bells and whistles for new kelp plot
  if(plot_type == "new_kelp"){
    new_plot <- base_pred_plot +
      geom_hline(yintercept= 0, linetype = "dashed", color = features, linewidth = line_var*0.25) +
      labs(y = expression(paste(Delta, " Ammonium ", (mu*M))), 
           x = x_axis_lab) +
      theme(plot.margin = unit(rep(0.1, 4), "lines"))
    #      scale_x_continuous(breaks = c(-1.17905227, -0.1, 1, 2.05),
    #                         labels = c("0", "0.6", "1.2", "1.8"))
  }
  
  # extras for rls plot
  if(plot_type == "rls"){
    new_plot <- base_pred_plot +
      #      scale_x_continuous(breaks = c(-1.85, -0.9, 0.05, 1, 1.95),
      #                         labels = c("300", "600", "900", "1200", "1500")) +
      #      ylim(0, 3.671) +
      labs(y = expression(paste("Ammonium ", (mu*M))), 
           x = expression(paste("Animal abundance/m"^2))) +
      scale_fill_manual(values = pal, drop = FALSE) +
      scale_colour_manual(values = (pal), drop = FALSE) +
      theme(legend.position = c(0.83, 0.865))
  }
    
    if(plot_type == "flow"){
      new_plot <- base_pred_plot +
        labs(x = "Flow rate (m/s)", 
             y = expression(paste(Delta, " Ammonium ", (mu*M)))) +
        labs(colour = "Cucumbers", fill = "Cucumbers", pch = "Cucumbers", lty = "Cucumbers")
    }
    

  
  print(new_plot)
}


# Standardize variables -----

scale_vars <- function(datafile){
  {{datafile}} %>%
    mutate(
      # kelp forest level variables
      forest_biomass = BiomassM*Area_m2,
      den_scale = c(scale(DensityM)), # make sure density is right
      kelp_bio_scale = c(scale(BiomassM)),
      log_kelp = log(BiomassM + 0.001),
      log_kelp_scale = c(scale(log_kelp)),
      forest_bio_scale = c(scale(forest_biomass)),
      area_scale = c(scale(Area_m2)),
      # transect level variables
      bio_tran_scale = c(scale(biomass_trans_mean)),
      log_kelp_tran = log(biomass_trans_mean + 0.001),
      den_tran_scale = c(scale(kelp_den)), # make sure density is right
      # log the pee diff?
      log_pee_diff = log(in_minus_out + 1),
      # the biomass + abundance variables
      weight_sum_scale = c(scale(weight_sum)),
      all_weighted_scale = c(scale(all_weight_weighted)),
      abundance_scale = c(scale(abundance)),
      # biodiversity variables
      rich_scale = c(scale(species_richness)),
      shannon_scale = c(scale(shannon)),
      simpson_scale = c(scale(simpson)),
      # abiotic variables I should control for
      depth_scale = c(scale(depth_avg)),
      tide_scale = c(scale(avg_exchange_rate)),
      tide_cat = factor(as.factor(ifelse(avg_exchange_rate < -0.1897325, "Ebb",
                                         ifelse(avg_exchange_rate < 0.1897325, "Slack", "Flood"))),
                        levels = c("Ebb", "Slack", "Flood")),
      # only slack and flood
      
      # try to center instead of scaling
      kelp_bio_center = c(scale(BiomassM, scale = FALSE)),
      tide_center = c(scale(avg_exchange_rate, scale = FALSE)),
      weight_sum_center = c(scale(weight_sum, scale = FALSE)),
      abundance_center = c(scale(abundance, scale = FALSE)),
      shannon_center = c(scale(shannon, scale = FALSE)),
      simpson_center = c(scale(simpson, scale = FALSE)),
      depth_center = c(scale(depth_avg, scale = FALSE))
    ) 
}


# Family function for RLS Blitz ------
# Family function to get predictions for RLS Blitz
fam_fun_combo <- function(df, family) {
  
  df_fam <- df %>% filter(family == {{family}})
  
  # set levels for ggpredict
  min_fam <- min(df_fam$fam_den)
  max_fam <- max(df_fam$fam_den)
  spread_fam <- (max_fam - min_fam)/100
  v_fam <- seq(min_fam, max_fam, spread_fam)
  
  # model
  mod_fam <- glmmTMB(nh4_avg ~ fam_den*tide_scale + shannon_scale + depth_avg_scale  + (1|year) + (1|site_code), 
                     family = Gamma(link = 'log'),
                     data = df_fam)
  
  # ggpredict
  predict_fam <- ggpredict(mod_fam, terms = c("fam_den[v_fam]")) %>%
    as.data.frame() %>%
    mutate(fam_den = x,
           family = as.character({{family}}),  # Ensure `family` is coerced to character if necessary
           r2 = as.numeric(performance::r2(mod_fam, tolerance = 0.0000000000001)[1])) 
  
}

# Diagnose fam models
diagnose_fun <- function(df, family){
  # Print family name so I know which results are which
  print({{family}})
  
  df_fam <- df %>% filter(family == {{family}})
  
  # model
  mod_fam <- glmmTMB(nh4_avg ~ fam_den_scale*tide_scale + shannon_scale + depth_avg_scale  + (1|year) + (1|site_code), 
                     family = Gamma(link = 'log'),
                     data = df_fam)
  
  print(summary(mod_fam))
  print(plot(DHARMa::simulateResiduals(mod_fam)))
  print(performance::r2(mod_fam, tolerance = 0.0000000000001))
}


# Family function for the kelp pee model
fam_fun_kelp_combo <- function(df, family) {
  
  df_fam <- df %>% filter(family == {{family}})
  
  # set levels for ggpredict
  min_fam <- min(df_fam$weight_fam_sum_g)
  max_fam <- max(df_fam$weight_fam_sum_g)
  spread_fam <- (max_fam - min_fam)/100
  v_fam <- seq(min_fam, max_fam, spread_fam)
  
  # model
  mod_fam <- glmmTMB(in_out_avg ~ kelp_sp + 
                       kelp_bio_scale*tide_scale +
                       kelp_bio_scale*weight_fam_sum_g +
                       weight_fam_sum_g*tide_scale +
                       shannon_scale + depth_scale, 
                     family = 'gaussian',
                     data = df_fam)
  
  # ggpredict
  predict_fam <- ggpredict(mod_fam, terms = c("weight_fam_sum_g[v_fam]")) %>%
    as.data.frame() %>%
    mutate(weight_fam_sum_g = x,
           family = {{family}},
           r2 = as.numeric(performance::r2(mod_fam, tolerance = 0.0000000000001)[1])) 
  
}

# Diagnose kelp pee model
diagnose_kelp_fun <- function(df, family){
  # Print family name so I know which results are which
  print({{family}})
  
  df_fam <- df %>% filter(family == {{family}})
  
  # model
  mod_fam <- glmmTMB(in_out_avg ~ kelp_sp + 
                       kelp_bio_scale*tide_scale +
                       kelp_bio_scale*weight_den_fam_scale +
                       weight_den_fam_scale*tide_scale +
                       shannon_scale + depth_scale, 
                     family = 'gaussian',
                     data = df_fam)
  
  print(summary(mod_fam))
  print(plot(DHARMa::simulateResiduals(mod_fam)))
  print(performance::r2(mod_fam, tolerance = 0.0000000000001))
}

