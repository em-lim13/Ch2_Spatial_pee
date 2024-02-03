# Graveyard ----

# RLS graveyard -----

# plotting graveyard
# plot mean for each year
rls_nh4 %>%
  group_by(site, year) %>%
  mutate(nh4_avg = mean(nh4_conc)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x = reorder(site, -nh4_avg), 
                 nh4_avg, colour = year, fill = year, pch = month),
             size = 3, alpha = 0.75) +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5)) +
  labs(x = "Site", y = "NH4+ Concentration (umol/L)") +
  scale_colour_manual(values = pal) 


# Plot the ranking of each site year to year
# Plot rank
ggplot(data = rank) +
  geom_point(aes(reorder(site, -avg_grade), grade2021, colour = "2021"), size = 3) +
  geom_point(aes(reorder(site, -avg_grade), grade2022, colour = "2022"), size = 3) +
  geom_point(aes(reorder(site, -avg_grade), grade2023, colour = "2023"), size = 3) +
  labs(x= "Site", y = "Rank", colour = "Year") +
  scale_colour_manual(values = pal) + 
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5)) 

# ggsave("Output/Figures/rank_all_years.png", device = "png",
#        height = 5, width = 8, dpi = 400)


# Site vs pee
ggplot(rls_nh4, aes(nh4_conc, site, fill = site)) +
  geom_boxplot(colour = "white") +
  theme_black() +
  labs(x= "Ammonium concentration (umol/L)", y = "Site") +
  theme(legend.position = "none")

#ggsave("Output/Figures/RLS_sites_pee.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# just urchins
urchins <- rls %>%
  filter(species_name == "Mesocentrotus franciscanus") %>%
  group_by(survey_id) %>%
  summarize(urchins = sum(total)) %>%
  mutate(Composition = "None") %>%
  rename(site_code = survey_id) %>%
  rbind(urchins_kelp)

ggplot(urchins, aes(Composition, urchins, colour = Composition)) +
  geom_boxplot() +
  geom_jitter()


#  mutate(xjit = ifelse(site_code == "BMSC6" | 
#                         site_code == "BMSC3" |
#                         site_code == "BMSC21" |
#                         site_code == "BMSC15" |
#                         site_code == "BMSC11" |
#                         site_code == "BMSC19", -0.015, 0.015))
#
## Create legend
#site_names <- as.list(paste(coords$site_num, coords$`Site name`, sep = ": "))


#  + geom_text(data = rls_coords,
#            aes(x = survey_longitude, y = survey_latitude, 
#                label = site_code),
#            size = 3,
#            nudge_x = rls_coords$xjit)


# filter inverts and weirdos for rls
filter(species_name != "Bolinopsis infundibulum") %>%
  filter(species_name != "Pleuronichthys coenosus") %>%
  filter(species_name != "Pleurobrachia bachei") %>%
  filter(species_name != "Polyorchis penicillatus") %>% # Filter inverts
  filter(species_name != "Actinopterygii spp.") %>% # Remove unidentif fish
  
  
  # Put all predictors in a model
  mod_full <- lmer(nh4_avg ~ scale(weight_sum) * scale(shannon) * scale(abundance) * scale(avg_exchange_rate) + (1|site_code), rls_final)
summary(mod_full)

# Put all predictors in a model + depth
# Cut the triple + interactions
mod_fuller <- lmer(nh4_avg ~ scale(weight_sum) + scale(simpson) + scale(abundance) + scale(avg_exchange_rate) + scale(depth_avg) +
                     scale(weight_sum):scale(simpson) + scale(weight_sum):scale(abundance) + 
                     scale(weight_sum):scale(avg_exchange_rate) + scale(weight_sum):scale(depth_avg) +
                     scale(simpson):scale(abundance) + scale(simpson):scale(avg_exchange_rate) + scale(simpson):scale(depth_avg) +
                     scale(abundance):scale(avg_exchange_rate) + scale(abundance):scale(depth_avg) +
                     scale(avg_exchange_rate):scale(depth_avg) +
                     (1|site_code), rls_final)
summary(mod_fuller)


# dredge to compare all models
options(na.action = "na.fail")
dred <- dredge(mod_fuller)
# best model is just the intercept, LOL
# 2nd best (<delta AIC 2) is just richness
# 3rd is just weight, but that's > delta AIC 5
# top 3 models and delta AIC is the same when full mod has all interactions vs when it has none

# When I add depth and only the 2-way interactions top mod = intercept, 2nd = richness (delta 0.73), 3rd is just depth (delta 5.6)

# look at this "best mod"
mod_rich <- lmer(nh4_avg ~ species_richness + (1|site_code), rls_final)
summary(mod_rich)
visreg(mod_rich)
# Negative relationship between species richness and NH4+...... cool cool cool cool cool



# What if I put allll the data together?????
mod_all <- lm(nh4_avg ~ scale(weight_sum) + scale(avg_exchange_rate)  + scale(BiomassM) + 
                scale(weight_sum):avg_exchange_rate + scale(weight_sum):scale(BiomassM), big_rls)

summary(kelp_mod_best)
visreg(kelp_mod_best)



# reduce the df so i can join it with kelp rls data
rls_final_reduced <- rls_final %>%
  mutate(BiomassM = 0) %>%
  select(site, site_code, survey_id, nh4_avg, depth_avg, avg_exchange_rate, BiomassM, weight_sum, species_richness, abundance, avg_exchange_rate,  depth_avg) %>%
  rename(site_code = site_code)

big_rls <- rbind(rls_final_reduced, data_s_reduced) %>%
  mutate(weight_stand = scale(weight_sum),
         rich_stand = scale(species_richness),
         abundance_stand = scale(abundance),
         tide_stand = scale(avg_exchange_rate),
         depth_stand = scale(depth_avg),
         biomass_mean_stand = scale(BiomassM))

# Just look at temp 
temp_df <- rls_final %>%
  drop_na(temp_avg)

mod_temp <- lm(nh4_avg ~ scale(weight_sum) + scale(species_richness) + scale(abundance) + scale(avg_exchange_rate) + scale(temp_avg) + scale(depth_avg), temp_df)
summary(mod_temp)
# Interesting. In 2021, there was a slightly negative relationship between temperature and nh4+ concentration


# Dot and whisker?
sum_stats_pee <- ggpredict(simple_model, terms = c("site_code", "year")) %>% 
  #and then we'll just rename one of the columns so it's easier to plot
  rename(site_code = x,
         nh4_conc = predicted,
         year = group)
# plot
ggplot() +
  geom_point(data = sum_stats_pee, 
             aes(y = site_code, x = nh4_conc, colour = year),
             size = 4) +
  geom_errorbar(data = sum_stats_pee, 
                aes(y = site_code,
                    x = nh4_conc,
                    colour = year,
                    # and you can decide which type of error to show here
                    # we're using 95% CI
                    xmin = conf.low,
                    xmax = conf.high),
                width = 0.2,
                size = 1.2)  +
  geom_point(data = rls_nh4, aes (y = site_code, x = nh4_conc, colour = year), alpha = 0.5, height = 0, size = 2) +
  labs(x= "Ammonium concentration (umol/L)", y = "Site") +
  theme_black() +
  theme(legend.position="none") 

#ggsave("Output/Figures/RLS_pee_black.png", device = "png",
#       height = 9, width = 16, dpi = 400)


# Kelp pee graveyard -----

# kelp pee prediction plot

# get original axis
#ggplot(data, aes(BiomassM, kelp_bio_scale)) +
#  geom_smooth() +
#  geom_vline(xintercept= 1.8,  color = "red", linewidth = 0.5) +
#  scale_y_continuous(n.breaks = 50)

# # Ammonium at the site level?
data_map %>%
  mutate(site = fct_reorder(site, in_out_avg, .fun='median')) %>%
  ggplot() +
  geom_point(aes(x = reorder(site, in_out_avg), nh4_out_avg, 
                 colour = "blue")) +
  geom_point(aes(x = reorder(site, in_out_avg), nh4_in_avg, 
                 colour = "green")) +
  labs(y = "NH4 concentration", x = "Site") + 
  theme(axis.text.x=element_text(angle = 65, hjust = 1))+
  scale_color_identity(name = "Sample",
                       breaks = c("blue", "green"),
                       labels = c("Outside kelp", "Inside kelp"),
                       guide = "legend")
# Ok this is sort of neat, you can see the base ammonium levels and then how different they are, by site. might be neat to arrange these with kelp density instead of site on the x

ggplot(data) +
  geom_point(aes(x = kelp_bio_scale, nh4_out_avg, 
                 colour = "blue")) +
  geom_point(aes(x = kelp_bio_scale, nh4_in_avg, 
                 colour = "green")) +
  labs(y = "NH4 concentration", x = "Kelp biomass") + 
  theme(axis.text.x=element_text(angle = 65, hjust = 1))+
  scale_color_identity(name = "Sample",
                       breaks = c("blue", "green"),
                       labels = c("Outside kelp", "Inside kelp"),
                       guide = "legend")

# are the continuous predictors the same for all kelp species :|
# what if I make nereo first instead of macro?
n_data <- data %>%
  mutate(kelp_sp = factor(kelp_sp, levels = c("nereo", "macro", "mixed", "none")))

mod_in_out_nereo <- glmmTMB(in_minus_out ~ -1 + kelp_sp + kelp_bio_scale*tide_scale*weight_sum_scale + shannon_scale + depth_scale - kelp_bio_scale:tide_scale:weight_sum_scale +
                              (1|site_code), 
                            family = 'gaussian',
                            data = n_data)

summary(mod_in_out_nereo)

# plot model
tt <- (list(no_mix = mod_in_out2, mix = mod_in_out_nereo)
       %>% purrr::map_dfr(tidy, effects = "fixed", conf.int = TRUE,
                          .id = "model")
       %>% select(model, component, term, estimate, conf.low, conf.high)
       ## create new 'term' that combines component and term
       %>% mutate(term_orig = term,
                  term = forcats::fct_inorder(paste(term, component, sep = "_")))
)

dwplot(tt)

# it doesn't seem to matter how I order the kelp_sp factor, all of the estimates are EXACTLY the same



# old kelp pee coeff plot using mod with intercept
int <- confint(mod_in_out, estimate = TRUE)[1,3]

df2 <- confint(mod_in_out, level = 0.95, method = c("wald"), component = c("all", "cond", "zi", "other"), estimate = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(variable = rowname,
         lower_CI = `2.5 %`,
         upper_CI = `97.5 %`,
         estimate = Estimate) %>%
  head(- 1)  %>%
  mutate(variable = factor(as.factor(variable), 
                           levels = c("(Intercept)", "kelp_spnereo", "kelp_spnone", "shannon_scale", "depth_scale", "kelp_bio_scale", "weight_sum_scale", "tide_scale", "kelp_bio_scale:tide_scale", "kelp_bio_scale:weight_sum_scale", "tide_scale:weight_sum_scale"),
                           labels = c("Macro", "Nereo", "No kelp", "Biodiversity", "Depth",  "Kelp biomass", "Animal biomass", "Tide", "Kelp:tide", "Kelp:animals", "Tide:animals")),
         adj_estimate = case_when(variable == "Nereo" ~ estimate + int,
                                  variable == "No kelp" ~ estimate + int,
                                  TRUE ~ estimate),
         se = (upper_CI - estimate)/1.96,
         ci.lb_adjust = adj_estimate - (1.96*se),
         ci.up_adjust = adj_estimate + (1.96*se)
  )


# just urchins
urchins_kelp <- kelp_rls %>%
  filter(species_name == "Mesocentrotus franciscanus") %>%
  group_by(site_code) %>%
  summarize(urchins = sum(total)) %>%
  left_join(kelp %>% select(site_code, Composition) %>% unique())


# try making a more simple model???
mod_gam_kelp3 <- glmmTMB(nh4_out_avg ~ weight_sum_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp3)) # two zig zags but no red

mod_gam_kelp4 <- glmmTMB(nh4_out_avg ~ tide_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp4)) # super fucked

mod_gam_kelp5 <- glmmTMB(nh4_out_avg ~ kelp_bio_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp5)) # three zigzags, 2 red

mod_gam_kelp6 <- glmmTMB(nh4_out_avg ~ shannon_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp6)) # super fucked

mod_gam_kelp7 <- glmmTMB(nh4_out_avg ~ kelp_sp, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp7)) # no signif problems, useable

mod_gam_kelp8 <- glmmTMB(nh4_out_avg ~ depth_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp8))  # a little wavy but ok

# mod with just the variables that didn't throw errors for nh4_outside?
mod_gam_kelp9 <- glmmTMB(nh4_out_avg ~ weight_sum_scale + kelp_sp + depth_scale, 
                         family = Gamma(link = 'log'),
                         data = data_s)
plot(simulateResiduals(mod_gam_kelp9)) # looks ok

aic <- AIC(mod_gam_kelp, mod_gam_kelp2, mod_gam_kelp3, mod_gam_kelp4, mod_gam_kelp5, mod_gam_kelp6, mod_gam_kelp7, mod_gam_kelp8, mod_gam_kelp9)
# ok so the two model with all the terms and lots of interactions are preferred


# so i accidentally ran the models with nh4_in_avg, and kelp_bio_scale*tide_scale + kelp_sp was the "best model" with low residual issues and AIC..... 

# mod with just the variables that didn't throw errors for nh4_inside?
mod_gam_kelp_in <- glmmTMB(nh4_in_avg ~ kelp_bio_scale*tide_scale + kelp_sp, 
                           family = Gamma(link = 'log'),
                           data = data_s)
summary(mod_gam_kelp_in) 
plot(simulateResiduals(mod_gam_kelp_in)) # ehhh? U-shape but less fucked
# dropping the interaction makes it worse
# forest_bio_scale makes it worse
# dropping bio mean scale makes it worse
# dropping tide makes it worse
# dropping kelp_sp gets rid of red line but all results non signif
# can't do tide*kelp_sp interaction or kelp_bio_scale*kelp_sp


# try to plot visreg output nicer
p = visreg(mod_best, "kelp_bio_scale", by="weight_sum_scale",
           overlay = TRUE, partial = FALSE, rug = FALSE,
           plot=FALSE)

ggplot(p$fit, aes(kelp_bio_scale, visregFit, 
                  fill = factor(weight_sum_scale))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha = 0.2) +
  geom_line(aes(colour = factor(weight_sum_scale)), linewidth = 1) +
  labs(linetype="weight_sum_scale", fill="weight_sum_scale", colour = "weight_sum_scale")


# Old stats -----
# I don't think this is linear I'm going to have to fit some kind of curve to it.....

# linear model though for fun
kelp_pee_mod <- lmer(in_minus_out ~ kelp_den + (1|site), data = data)
summary(kelp_pee_mod)
visreg(kelp_pee_mod)

# full model for dredging
model_all <- lmer(in_minus_out ~ den_scale + bio_tran_scale + kelp_bio_scale + area_scale + forest_bio_scale + kelp_sp + (1|site), data = data, na.action = na.fail)
summary(model_all)
visreg(model_all)

dredge <- as.data.frame(dredge(model_all)) %>%
  filter(delta < 3)

# the best model has the mean biomass of the whole forest
# best has mean bio + total forest biomass
# second best is just mean bio
# third is mean bio + total forest biomass + total area

bio_mod <- lmer(in_minus_out ~ kelp_bio_scale + (1|site), data = data)
summary(bio_mod)
visreg(bio_mod)
# Forests with more mean biomass/m2 retain more pee!

# redo dredge with the rest of the vars
model_all <- lmer(in_minus_out ~ den_scale + bio_tran_scale + kelp_bio_scale + area_scale + forest_bio_scale + kelp_sp +
                    weight_stand + rich_scale + abundance_scale + tide_scale + depth_stand + (1|site), data = data, na.action = na.fail)
summary(model_all)

dredge <- as.data.frame(dredge(model_all)) %>%
  filter(delta < 3)

# plot?
ggplot(data, aes(kelp_bio_scale, in_minus_out)) +
  geom_point() +
  geom_smooth(method = lm) 

# Let's be intelligent and build the model I think would be best
# I'd guess the biomass/m2 (density x biomass) = bio_tran_scale would matter bc that's the info about the kelp on the transect the samples were taken on
# And I'd guess the mean biomass/m2 (den x bio) of the kelp forest (how much kelp is around) will matter
mod1 <- lmer(in_minus_out ~ bio_tran_scale* kelp_bio_scale  + (1|site), data = data, na.action = na.fail)
summary(mod1)
visreg(mod1, "bio_tran_scale", by = "kelp_bio_scale")
# It looks like at high mean biomass, the relationship levels out
# So maybe there's an asymptote!

# let's try taking the log of the  in_minus_out and see if that improves fit

mod2 <- lmer(log_pee_diff ~ bio_tran_scale* kelp_bio_scale  + (1|site), data = data, na.action = na.fail)
summary(mod1)
visreg(mod1, "bio_tran_scale", by = "kelp_bio_scale")


# go back to the preferred model with just mean biomass
bio_mod2 <- lmer(log_pee_diff ~ kelp_bio_scale  + (1|site), data = data)
summary(bio_mod2)
visreg(bio_mod2)

# check resid
mod_cont_resids <- simulateResiduals(bio_mod)
plot(mod_cont_resids)
# weight fucking shit



# Use AIC to see if this improves things
AIC(bio_mod, bio_mod2, mod1)
# Yes bio_mod2 is preferred, taking the log of the response variable helps a lot



# Site level model!

kelp_mod_full <- lm(nh4_out_avg ~ weight_stand + rich_scale + abundance_scale + tide_scale + depth_stand + kelp_bio_scale + 
                      weight_stand:tide_scale + weight_stand:kelp_bio_scale +
                      rich_scale:tide_scale + rich_scale:kelp_bio_scale +
                      abundance_scale:tide_scale +
                      abundance_scale:kelp_bio_scale +
                      tide_scale:kelp_bio_scale, data_s)
summary(kelp_mod_full)

options(na.action = "na.fail")
dredge <- as.data.frame(dredge(kelp_mod_full)) %>%
  filter(delta < 3)
# nh4 avg second best mod was just tide
# nh4 inside avg second best mod was just bio mean, then just tide
# nh4 outside avg second best mod just tide

kelp_mod_mean <- lm(nh4_avg ~ weight_stand + tide_scale  + kelp_bio_scale + 
                      weight_stand:tide_scale + weight_stand:kelp_bio_scale, data_s)


summary(kelp_mod_best)

visreg(kelp_mod_best)


# Plots ----

ggplot(data, aes(kelp_bio_scale, log_pee_diff)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_smooth(method = lm, colour = "blue")+
  geom_smooth(method = loess, colour = "red")+
  
  labs(y = "Inside - outside kelp forest ammonium (log(uM))", x = "Forest biomass") 


# each point is a single transect with a pee difference and a kelp density 
# linear model
ggplot(data, aes(BiomassM, in_minus_out)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_hline(yintercept= 0, linetype = "dashed", color = "red", size = 1.5) +
  labs(y = "Inside - outside kelp forest ammonium (uM)", x = "Kelp Density") +
  geom_smooth(method = lm,
              alpha = 0.25)

# try plotting the same data but with a better asymptote?
ggplot(data, aes(kelp_den, in_minus_out)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_hline(yintercept= 0, linetype = "dashed", color = "red", size = 1.5) +
  labs(y = "Inside - outside kelp forest ammonium (uM)", x = "Kelp Density") +
  geom_smooth(method = "loess",
              span = 1,
              alpha = 0.25)

# Percent difference 
# With a truly heinious curve, thanks loess
ggplot(data, aes(kelp_den, percent_diff)) +
  geom_point(aes(pch = kelp_sp, colour = site_code), size =3 )+ 
  geom_hline(yintercept= 0, linetype = "dashed", color = "red", size = 1.5) +
  labs(y = "Inside - outside kelp forest ammonium (uM)", x = "Kelp Density") +
  geom_smooth(method = "loess",
              span = 1,
              alpha = 0.25)

#ggsave("Output/Figures/in_out_kelp_pee.png", device = "png",
#       height = 9, width = 16, dpi = 400)



# Claire has area, density, biomass
# Got average individual biomass for each kelp, which you can scale up to the average biomass of each transect


# Old Kelp density data -----
# Also Claire's biomass estimates
# Basically it would be nice to have a single metric of how much kelp is in each forest
kcca_summary <- kcca_final %>%
  group_by(site_code) %>%
  summarise(kelp_den = mean(kelp_den),
            in_minus_out = mean(in_minus_out),
            kelp_sp = kelp_sp)

# Old map code -----

# OLD map functions! -----

# Site map
old_site_map <- function(lat_min, lat_max, long_min, long_max, 
                         coord_data, map_data, 
                         add_points, add_annotate){
  
  # make basic map plot
  gg <- ggplot() +
    geom_sf(data = {{map_data}}, fill = "white", colour = blue) +
    coord_sf(xlim = c({{lat_min}}, {{lat_max}}), ylim = c({{long_min}}, {{long_max}}), expand = FALSE) +
    # add aesthetic elements
    theme_black() +
    theme(panel.background = element_rect(fill = blue),
          panel.grid.major = element_line(color = blue)) +
    # axis things
    labs(x = "Longitude", y = "Latitude",
         fill = "Habitat") +
    scale_x_continuous(breaks = seq({{lat_min}}, {{lat_max}}, by = 0.1))
  
  # add bells and whistles
  if(add_points == TRUE){
    gg + geom_sf(data = {{coord_data}}, 
                 colour = "black",
                 fill = "black",
                 alpha = 0.5,
                 size = 8,
                 aes(pch = Habitat)) +
      scale_shape_manual(values = c(21, 25), drop = F) +
      guides(pch = guide_legend(override.aes = 
                                  list(colour = "black")))} 
  
  # add arrow + scale bar
  if(add_annotate == TRUE){
    gg +
      annotation_scale(location = "br", width_hint = 0.4) +
      annotation_north_arrow(location = "br", which_north = "true", 
                             pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                             style = north_arrow_fancy_orienteering)
  }
}

# Site map from Nikola
site_map <- function(lat_min,
                     lat_max,
                     long_min,
                     long_max,
                     coord_data,
                     map_data,
                     add_points,
                     add_annotate) {
  # make basic map plot
  plot <- ggplot() +
    geom_sf(data = map_data,
            fill = "white",
            colour = blue) +
    coord_sf(
      xlim = c(lat_min, lat_max),
      ylim = c(long_min, long_max),
      expand = FALSE
    ) +
    # add aesthetic elements
    theme_black() +
    theme(
      panel.background = element_rect(fill = blue),
      panel.grid.major = element_line(color = blue)
    ) +
    # axis things
    labs(x = "Longitude", y = "Latitude",
         fill = "Habitat") +
    scale_x_continuous(breaks = seq(lat_min, lat_max, by = 0.1))
  
  if (add_points) {
    plot <- plot + geom_sf(
      data = coord_data,
      colour = "black",
      fill = "black",
      alpha = 0.5,
      size = 8,
      aes(pch = Habitat)
    ) +
      scale_shape_manual(values = c(21, 25), drop = F) +
      guides(pch = guide_legend(override.aes =
                                  list(colour = "black")))
  }
  if (add_annotate) {
    plot <- plot + annotation_scale(location = "br", width_hint = 0.4) +
      annotation_north_arrow(
        location = "br",
        which_north = "true",
        pad_x = unit(0.0, "in"),
        pad_y = unit(0.2, "in"),
        style = north_arrow_fancy_orienteering
      )
  }
  print(plot)
}


# big inset map
inset_map <- function(rect_xmin, rect_xmax, rect_ymin, rect_ymax, 
                      map_data){
  
  ggplot() +
    geom_sf(data = {{map_data}}, fill = "white", colour = blue) +
    coord_sf(xlim = c(-128.5, -123), ylim = c(48.25, 51), expand = FALSE) +
    # add rectangle for zoomed in part
    geom_rect(aes(xmin = {{rect_xmin}}, xmax = {{rect_xmax}}, ymin = {{rect_ymin}}, ymax = {{rect_ymax}}), color = "red", fill = NA, inherit.aes = FALSE) +
    # add aesthetic elements
    theme_bw() +
    theme(panel.background = element_rect(fill = blue),
          panel.grid.major = element_line(color = blue),
          panel.border = element_rect(fill = NA, colour = "black"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(0, "pt"),
          plot.title = NULL,
          plot.margin=grid::unit(c(0,0,0,0), "mm"))
}

# add label for Vancouver
# geom_text(aes(x = -123.120694, y = 49.282694, 
#               label = "Vancouver",
#               fontface = "bold")) 


# load the osm data map from nikola - ----
load("Data/bc_map.Rdata")
bc_map <- slice # rename


# just slack and ebb
coords_slack <- rls_coords %>%
  filter(tide_cat != "Flood") %>%
  group_by(site_code) %>%
  mutate(nh4_avg = mean(nh4_avg)) %>%
  ungroup()

# use site map! -----
site_map(
  lat_min = -125.3,
  lat_max = -125,
  long_min = 48.8,
  long_max = 48.9,
  coord_data = all_coords,
  map_data = potato_map,
  add_points = TRUE,
  add_annotate = TRUE
)


# just rls sites
barkley_rls <- site_map(lat_min = -125.6, lat_max = -124.95, long_min = 48.75, long_max = 49.1,
                        coord_data = all_coords %>%
                          filter(Habitat == "Reef"), map_data = potato_map,
                        add_points = TRUE, add_annotate = TRUE) +
  guides(pch = "none")

van_isle_rls <- inset_map(rect_xmin = -125.6, rect_xmax = -124.95, rect_ymin = 48.75, rect_ymax = 49.1,
                          map_data = potato_map)


barkley_rls + 
  inset_element(
    van_isle_rls, 
    left = 0, 
    bottom = 0.5, 
    right = 0.5, 
    top = 1.03,
    align_to = 'panel'
  )



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
