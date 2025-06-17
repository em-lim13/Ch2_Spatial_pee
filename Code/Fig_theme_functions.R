# Figure functions
# Em Lim
# Updated June 2025

# This one script should be the same between all of my projects so I can make really consistant figures

# Load packages ----------------------------------------------------------------
library(tidyverse)
library(patchwork)

# Themes -----------------------------------------------------------------------

## Theme black ----------------------------------------------------------------
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

## Theme white ----------------------------------------------------------------
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

## Publication figure theme -----------------------------------------------------
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

# Plot functions ---------------------------------------------------------------

## add annotation label (a, b, c) ----------------------------------------------
place_label <- function(label, size = 4.5, ...) {
  annotate("text", label = label, x = -Inf, y = Inf, 
           vjust = 1.4, hjust = -0.15, size = size, ...)
}