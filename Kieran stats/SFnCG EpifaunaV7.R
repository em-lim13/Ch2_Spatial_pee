
###############################################################################################################
# Cox et al. Intertidal Resource Cultivation Over Millennia Structures Marine Biodiversity  
# Correspondence to: kierancx@gmail.com
# 1. Data
# 2. Manuscript Figures 
# 3. Supplemental Analyses and Figures
### Note: The raw data needed to conduct the analysis is included in the folder. 
###############################################################################################################


###############################################################################################################

## 1. Data Structure
## 2 years, 24 sites, 3 regions, 3 site types, high-mid-low intertidal (quadrat density varies considerably) 
## Epifauna communities 
## Main comparisons Gardens-Reference, Aquaculture-References
### Focus mid intertidal CG-C, CG-Q, SF-Q, SF-B- only

## 2. Manuscript Figures 
## Figrue 1: Map, Coastal Map, 4 regions sub figures, legend
## Figure 2. Diversity Profiles curves
## Figure 3: nMDS (maybe supplemental) 
## Figure 4:  SIMPER plots, tables in supplemenal
## Supplemental: Permanova tables
## Figure 5: Multivariate Regression Tree's Benthic Diversity, Multivariate Random Forest
## Supplemental: Site Complexity models
## Figure 6: 3D model analysis. Complexity influnece on site diversity



## TBD
## 4 themes: benthic composition, 3D complexity
## 4 tools: HMSC, Regression Trees, Mixed effect models, Random Forest

## Current thinking
## Figure XXX: Benthic composition, regression tree's, maybe train on calvert longterm data plus others?
## Figure XXX: 3D Complexity- HMSC species by phylogeny of species, coloured habitat importance values


### 1. Data ####
rm(list = ls(all = TRUE))

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")
SFandCGData=read.csv("SFandCGBiodiversityDataNov2020 V5 .csv",header=T)

SFandCGData
#View(SFandCGData)
head(SFandCGData)
names(SFandCGData)
levels(SFandCGData$Area)
levels(SFandCGData$Site)
levels(SFandCGData$Quad)
SFandCGData$Year=as.factor(SFandCGData$Year)
levels(SFandCGData$Year)
levels(SFandCGData$Zone)
levels(SFandCGData$Comparison)

library(data.table)
library(tidyverse)
library(vegan)
library(pryr)
library(ggplot2)
library(colordistance)
library(cowplot)
library(dplyr)
library(magick)

## 2. Manuscript Figures 
#### Figrue 1: Map, Coastal Map, 4 regions sub figures, legend ####
# packages
library(rgdal)
library(raster)
library(maps) 
library(mapdata)
library(maptools)
library(rgeos)
library(rgdal)
library(ggplot2)
library(ggsn)
library(tidyverse)
library(here)

# Map file
BC.shp <- readOGR(dsn="~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Map", layer="COAST_TEST2")

## Map cropping (optional)
### chose the lat/long extent you want to show
Ncalvert <- extent( c(-128.18, -127.94, 51.61, 51.78))

### if using cropped shapefile polygons to the extent defined
# takes a moment to run (patience grasshopper)
BC.shp2 <- crop(BC.shp,Ncalvert)

### project and fortify (i.e. turn into a dataframe)
BC.df <- fortify(BC.shp2)
BC.dfWhole <- fortify(BC.shp)

## Site Data
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Map")
SiteLocations=read.csv("Biodiversity Site Locations.csv", header=T)
SiteLocations

## Order
SiteLocations$AreaType
SiteLocations$AreaType = factor(SiteLocations$AreaType, levels = c("Clam Garden","Shellfish Farm","Garden Reference","Farm Reference"))

## Color
GardenColor=("#7CAE00")   #### Green 
FarmColor=("#0d7dd9")     #### Blue 
ReferenceGardenColor=("#fa6920") #### Orange
ReferenceFarmColor=("#ff1100") #### Red

col=c(GardenColor,FarmColor,ReferenceGardenColor,ReferenceFarmColor)

## Shape
GardenShape=21   
FarmShape=22     
ReferenceGardenShape=23 
ReferenceFarmShape=24

shape=c(GardenShape,FarmShape,ReferenceGardenShape,ReferenceFarmShape)

## Region: Calvert V1
CalvertMap =ggplot()+ theme_bw()+
  geom_polygon(data= BC.dfWhole, aes(x=long,y=lat,group= group),
               colour= "black", size=0.5, fill='grey84')+
  coord_cartesian(xlim = c(-128.175, -127.95), ylim=c(51.636, 51.73)) +
  geom_point(data=SiteLocations, aes(x=Long, y=Lat, shape=AreaType,color=AreaType,fill=AreaType), size=10, stroke=1.5)+
  scale_shape_manual(values=shape)+scale_color_manual(values=col)+scale_fill_manual(values=col)+
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=10), axis.text.x = element_text(size=10))+guides(color=FALSE,shape=FALSE,fill=FALSE)


## Region: North Quadra V1
NorthQuadraMap=ggplot()+ theme_bw()+
  geom_polygon(data= BC.dfWhole, aes(x=long,y=lat,group= group),
               colour= "black", size=0.5, fill='grey84')+
  coord_cartesian(xlim = c(-125.388, -125.27), ylim=c(50.23, 50.283)) +
  geom_point(data=SiteLocations, aes(x=Long, y=Lat, shape=AreaType,color=AreaType,fill=AreaType), size=10, stroke=1.5)+
  scale_shape_manual(values=shape)+scale_color_manual(values=col)+scale_fill_manual(values=col)+
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=10), axis.text.x = element_text(size=10))+guides(color=FALSE,shape=FALSE,fill=FALSE)

## Region South Quadra
SouthQuadraMap=ggplot()+ theme_bw()+
  geom_polygon(data= BC.dfWhole, aes(x=long,y=lat,group= group),
               colour= "black", size=0.5, fill='grey84')+
  coord_cartesian(xlim = c(-125.24, -125.03), ylim=c(50.08, 50.215)) +
  geom_point(data=SiteLocations, aes(x=Long, y=Lat, shape=AreaType,color=AreaType,fill=AreaType), size=10, stroke=1.5)+
  scale_shape_manual(values=shape)+scale_color_manual(values=col)+scale_fill_manual(values=col)+
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=10), axis.text.x = element_text(size=10))+guides(color=FALSE,shape=FALSE,fill=FALSE)

## Region Baynes Sound
BaynesMap=ggplot()+ theme_bw()+
  geom_polygon(data= BC.dfWhole, aes(x=long,y=lat,group= group),
               colour= "black", size=0.5, fill='grey84')+
  coord_cartesian(xlim = c(-124.94, -124.68), ylim=c(49.445, 49.625)) +
  geom_point(data=SiteLocations, aes(x=Long, y=Lat, shape=AreaType,color=AreaType,fill=AreaType), size=10, stroke=1.5)+
  scale_shape_manual(values=shape)+scale_color_manual(values=col)+scale_fill_manual(values=col)+
  scalebar(BC.dfWhole, dist = 2, st.size=3, height=0.0003, dist_unit = "km",transform = TRUE, model = "WGS84",anchor= c(x = -124.68, y = 49.5382)) +
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=10), axis.text.x = element_text(size=10))+guides(color=FALSE,shape=FALSE,fill=FALSE)

Legend=ggplot()+ theme_bw()+
  geom_polygon(data= BC.dfWhole, aes(x=long,y=lat,group= group),
               colour= "black", size=0.5, fill='grey84')+
  coord_cartesian(xlim = c(-124.94, -124.68), ylim=c(49.445, 49.625)) +
  geom_point(data=SiteLocations, aes(x=Long, y=Lat, shape=AreaType,color=AreaType,fill=AreaType), size=4, stroke=1.5)+
  scale_shape_manual(values=shape)+scale_color_manual(values=col)+scale_fill_manual(values=col)+
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=10), axis.text.x = element_text(size=10))#guides(color=FALSE,shape=FALSE,fill=FALSE)

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Map")
jpeg(filename = "Legend.jpeg", width = 20, height = 20, units = "cm", pointsize = 15, quality = 100, res = 300)  
Legend
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


BaynesMap+theme(axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                rect = element_blank())+   scalebar(BC.dfWhole, dist = 2, st.size=3, height=0.0003, dist_unit = "km",transform = TRUE, model = "WGS84",anchor= c(x = -124.68, y = 49.5382)) 


SouthQuadraMap+theme(axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks = element_blank(),
                     rect = element_blank())+ scalebar(BC.dfWhole, dist = 2, st.size=3, height=0.00027, dist_unit = "km",transform = TRUE, model = "WGS84",anchor= c(x = -125.03, y = 50.075))


NorthQuadraMap+theme(axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks = element_blank(),
                     rect = element_blank())+scalebar(BC.dfWhole, dist = 2, st.size=3, height=0.00008, dist_unit = "km",transform = TRUE, model = "WGS84",anchor= c(x = -125.275, y = 50.128))


CalvertMap+theme(axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 rect = element_blank())+scalebar(BC.dfWhole, dist = 2, st.size=3, height=0.00015, dist_unit = "km",transform = TRUE, model = "WGS84",anchor= c(x = -127.96, y = 51.633))

## Coastal 
CoastMap=ggplot()+ theme_bw()+
  geom_polygon(data= BC.dfWhole, aes(x=long,y=lat,group= group),
               colour= "black", size=0.5, fill='grey84')+
  coord_cartesian(xlim = c(-129.31943, -123.45), ylim=c(48.4, 52.050804)) +
  north(data = BC.dfWhole, scale = 1, symbol = 3, anchor= c(x = -126.15, y = 52)) +
  north(data = BC.dfWhole, scale = 1, symbol = 3, anchor= c(x = -125.15, y = 51)) +
  north(data = BC.dfWhole, scale = 1, symbol = 3, anchor= c(x = -124.15, y = 51.5)) +
  scale_shape_manual(values=shape)+scale_color_manual(values=col)+scale_fill_manual(values=col)+
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=10), axis.text.x = element_text(size=10))+guides(color=FALSE,shape=FALSE,fill=FALSE)
CoastMap

CoastMap=ggplot()+ theme_bw()+
  geom_polygon(data= BC.dfWhole, aes(x=long,y=lat,group= group),
               colour= "black", size=0.5, fill='grey84')+
  coord_cartesian(xlim = c(-170.31943, 0), ylim=c(0, 70.050804)) +
  scale_shape_manual(values=shape)+scale_color_manual(values=col)+scale_fill_manual(values=col)+
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=10), axis.text.x = element_text(size=10))+guides(color=FALSE,shape=FALSE,fill=FALSE)
CoastMap

world <- map_data("world",c("usa", "Canada","Mexico","Guatemala","Belize","Costa Rica","El Salvador","Honduras","Nicaragua","Panama"))

NorthAmerica=ggplot() +
  geom_map(data = world, map = world, aes(x = long, y = lat, map_id = region), color="black", fill = "grey69") +
   xlab("") + ylab("") +
  scale_x_continuous(limits = c(-180, -40), breaks = seq()) +
  scale_y_continuous(limits = c(6, 85), breaks = seq()) +
  scale_color_manual(values=col)+
  theme_bw() 
NorthAmerica

CalvertMapFinal=CalvertMap+theme(axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 rect = element_blank())+scalebar(BC.dfWhole, dist = 2, st.size=3, height=0.00015, dist_unit = "km",transform = TRUE, model = "WGS84",anchor= c(x = -127.96, y = 51.633))

NorthQuadraMapFinal=NorthQuadraMap+theme(axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks = element_blank(),
                     rect = element_blank())+scalebar(BC.dfWhole, dist = 2, st.size=3, height=0.00008, dist_unit = "km",transform = TRUE, model = "WGS84",anchor= c(x = -125.275, y = 50.228))

SouthQuadraMapFinal=SouthQuadraMap+theme(axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks = element_blank(),
                     rect = element_blank())+ scalebar(BC.dfWhole, dist = 2, st.size=3, height=0.00027, dist_unit = "km",transform = TRUE, model = "WGS84",anchor= c(x = -125.03, y = 50.075))

BaynesMapFinal=BaynesMap+theme(axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                rect = element_blank())+   scalebar(BC.dfWhole, dist = 2, st.size=3, height=0.0003, dist_unit = "km",transform = TRUE, model = "WGS84",anchor= c(x = -124.68, y = 49.4382)) 

CoastMapFinal=CoastMap+theme(axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks = element_blank(),
               rect = element_blank())

NorthAmerica

## Exporting images
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Map")
jpeg(filename = "CalvertMapFinal.jpeg", width = 30, height = 20, units = "cm", pointsize = 15, quality = 100, res = 300)  
CalvertMapFinal
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Map")
jpeg(filename = "NorthQuadraMapFinal.jpeg", width = 30, height = 20, units = "cm", pointsize = 15, quality = 100, res = 300)  
NorthQuadraMapFinal
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Map")
jpeg(filename = "SouthQuadraMapFinal.jpeg", width = 30, height = 20, units = "cm", pointsize = 15, quality = 100, res = 300)  
SouthQuadraMapFinal
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Map")
jpeg(filename = "BaynesMapFinal.jpeg", width = 30, height = 20, units = "cm", pointsize = 15, quality = 100, res = 300)  
BaynesMapFinal
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Map")
jpeg(filename = "CoastMapFinal6.jpeg", width = 60, height = 40, units = "cm", pointsize = 15, quality = 100, res = 300)  
CoastMap
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Map")
jpeg(filename = "CoastMapFinal5.jpeg", width = 60, height = 40, units = "cm", pointsize = 15, quality = 100, res = 300)  
CoastMapFinal
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Map")
jpeg(filename = "NorthAmerica.jpeg", width = 30, height = 20, units = "cm", pointsize = 15, quality = 100, res = 300)  
NorthAmerica
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


#### Figure 2. Diversity Profiles Curves ####
### mid intertidal focused figures 
SFandCGData   ## data
SFandCGData %>% count(Site, Year, Zone)
SFandCGDataMid <- subset(SFandCGData, SFandCGData$Zone == "Mid")
SFandCGDataMid %>% count(Site, Year,Zone)  # 3-5 Quadrats per mid zone
options(digits=5)

### Epifauna Data
names(SFandCGDataMid)
MidEpifaunaData <- SFandCGDataMid[,33:92]
MidDataStructure=SFandCGDataMid[,1:7]
EpifaunaSpeciesData2=cbind(MidDataStructure,MidEpifaunaData)

#### Diversity Profiles,Epifaunal, CG-C-Epi - Clam Gardens Calvert Island ####
CalvertMidEpifaunaData <- subset(EpifaunaSpeciesData2, EpifaunaSpeciesData2$Area == "Calvert")
# view(CalvertMidEpifaunaData)
str(CalvertMidEpifaunaData)
names(CalvertMidEpifaunaData)

###Remove zero columns
CalvertMidEpifaunaData2=CalvertMidEpifaunaData[, colSums(CalvertMidEpifaunaData != 0) > 0]
str(CalvertMidEpifaunaData2)
names(CalvertMidEpifaunaData2)

### Removing NAs   not needed as no NAs, keep for consistency across code
CalvertMidEpifaunaData3=(CalvertMidEpifaunaData2[complete.cases(CalvertMidEpifaunaData2), ])
str(CalvertMidEpifaunaData3)
names(CalvertMidEpifaunaData3)

CalvertMidEpifaunaData3$Area
CalvertMidEpifaunaData3$Site
CalvertMidEpifaunaData3$Quad
CalvertMidEpifaunaData3$SiteType

CalvertMidEpifaunaData3 %>% count(SiteType,Site, Year)  # 3-5 Quadrats per mid zone

## Diversity Profile
names(CalvertMidEpifaunaData3)
CalvertMidEpifaunaData4=aggregate(.~SiteType+Year+Site, CalvertMidEpifaunaData3, mean)
str(CalvertMidEpifaunaData4)
CalvertMidEpifaunaData5=aggregate(.~SiteType+Site, CalvertMidEpifaunaData4, mean)

CalvertMidEpifaunaData6=aggregate(.~SiteType, CalvertMidEpifaunaData5, mean)
# view(CalvertMidEpifaunaData6)
names(CalvertMidEpifaunaData6)

### Splitting by SiteType
MidEpifaunaCG= CalvertMidEpifaunaData6[1,]
MidEpifaunaReference= CalvertMidEpifaunaData6[2,]

### Removing zero's
MidEpifaunaCG2  =  MidEpifaunaCG[, colSums(MidEpifaunaCG != 0) > 0]
MidEpifaunaRef2  =  MidEpifaunaReference[, colSums(MidEpifaunaReference != 0) > 0]

names(MidEpifaunaCG2)
names(MidEpifaunaRef2)

MidEpifaunaCG3 = MidEpifaunaCG2[,8:35]
MidEpifaunaRef3 = MidEpifaunaRef2[,8:25]

### Equation
divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

MidEpifaunaRef4=divprof(MidEpifaunaRef3)

### Checking Adult Control Data
exp(diversity(MidEpifaunaRef3, index="shannon"))
###   2.3386
1 / (1-(diversity(MidEpifaunaRef3, index="simpson")))
###  1.4812

MidEpifaunaCG4=divprof(MidEpifaunaCG3)

### Checking Adult Treatment Data
exp(diversity(MidEpifaunaCG3, index="shannon"))
###   3.7876
1 / (1-(diversity(MidEpifaunaCG3, index="simpson")))
###  2.7498

## Building the plot

### Adding colors
GardenColor=alpha("#7CAE00",0.8)   #### Green 
FarmColor=alpha("#0d7dd9",0.8)     #### Blue 
ReferenceFarmColor=alpha("#ff1100",0.8) #### Red
ReferenceGardenColor=alpha("#fa6920",0.8) #### Orange

plot(MidEpifaunaCG4[,1],MidEpifaunaCG4[,2], col=GardenColor,cex=2, pch=16,xlab="",ylab="",ylim = c(0,30))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(MidEpifaunaRef4[,1],MidEpifaunaRef4[,2],pch=16,cex=2,col= ReferenceGardenColor)
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

points(MidEpifaunaCG4[1,1],MidEpifaunaCG4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaRef4[1,1],MidEpifaunaRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.5,29.08, expression( "Richness"["CG"] ))
text(0.5,19.08, expression( "Richness"["GR"] ))

### Now plotting Control shannon and Simpson
points(MidEpifaunaRef4[10,1],MidEpifaunaRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaRef4[19,1],MidEpifaunaRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(MidEpifaunaCG4[10,1],MidEpifaunaCG4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaCG4[19,1],MidEpifaunaCG4[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.45,2.8, expression( "Shannon"["GR"] ))
text(2.3,1.5, expression( "Simpson"["GR"] ))

text(1.45,4, expression( "Shannon"["CG"] ))
text(2.4,3.5, expression( "Simpson"["CG"] ))

### Plot it
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "CalvertEpifaunaOrderPlot6.jpeg", width = 22, height = 18, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(MidEpifaunaCG4[,1],MidEpifaunaCG4[,2], col=GardenColor,cex=2, pch=16,xlab="",ylab="",ylim = c(-1,30))
title(line = 2.5,xlab = list("Order q", cex = 1.3))
title(line = 2.5,ylab = list("Diversity", cex = 1.3))

points(MidEpifaunaRef4[,1],MidEpifaunaRef4[,2],pch=16,cex=2,col= ReferenceGardenColor)

points(MidEpifaunaCG4[1,1],MidEpifaunaCG4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaRef4[1,1],MidEpifaunaRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.49,28, expression( "Richness"["CG"] ))
text(0.49,18, expression( "Richness"["GR"] ))

points(MidEpifaunaRef4[10,1],MidEpifaunaRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaRef4[19,1],MidEpifaunaRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(MidEpifaunaCG4[10,1],MidEpifaunaCG4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaCG4[19,1],MidEpifaunaCG4[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.37,5.05, expression( "Shannon"["CG"] ))
text(1.37,0.5, expression( "Shannon"["GR"] ))

text(2.36,4, expression( "Simpson"["CG"] ))
text(2.36,0, expression( "Simpson"["GR"] ))

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")

### Setting consistent sizes
cex.axis1=1.5
cex1=2.2
cex2=1.7
cex3=2
cexsummary=1.85
cexText=1.2


CalvertCGEpi %<a-% {
  plot(MidEpifaunaCG4[,1],MidEpifaunaCG4[,2], cex.axis=cex.axis1,cex=cex1,col=GardenColor, pch=16,xlab="",ylab="",ylim = c(-1,30))
  title(line = 2.5,xlab = list("Order q", cex = cex2))
  title(line = 2.5,ylab = list("Diversity", cex = cex2))
  mtext(side=3, line=0.5,padj=-0.2,adj=0, "B)", col="Black", font=1, cex=cex3)
  
  points(MidEpifaunaRef4[,1],MidEpifaunaRef4[,2],pch=16,cex=2,col= ReferenceGardenColor)
  
  points(MidEpifaunaCG4[1,1],MidEpifaunaCG4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  points(MidEpifaunaRef4[1,1],MidEpifaunaRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  
  text(0.61,28, cex=cexText,expression( "Richness"["CG"] ))
  text(0.61,18, cex=cexText,expression( "Richness"["GR"] ))
  
  points(MidEpifaunaRef4[10,1],MidEpifaunaRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(MidEpifaunaRef4[19,1],MidEpifaunaRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  points(MidEpifaunaCG4[10,1],MidEpifaunaCG4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(MidEpifaunaCG4[19,1],MidEpifaunaCG4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  text(1.42,5.08,cex=cexText, expression( "Shannon"["CG"] ))
  text(1.42,0.5, cex=cexText,expression( "Shannon"["GR"] ))
  
  text(2.5,4, cex=cexText,expression( "Simpson"["CG"] ))
  text(2.5,0, cex=cexText,expression( "Simpson"["GR"] ))}

##### CG-C Done

#### Diversity Profiles,Epifaunal, CG-Q-Epi - Clam Gardens Quadra Island ####
EpifaunaSpeciesData2$Area
QuadraMidEpifaunaData <- subset(EpifaunaSpeciesData2, EpifaunaSpeciesData2$Area == "Quadra")
#  view(QuadraMidEpifaunaData)
QuadraCGMidEpifaunaData2 <- subset(QuadraMidEpifaunaData, QuadraMidEpifaunaData$Comparison == "Garden")
#  view(QuadraMidEpifaunaData2)

### Remove zero columns
QuadraCGMidEpifaunaData3=QuadraCGMidEpifaunaData2[, colSums(QuadraCGMidEpifaunaData2 != 0) > 0]
str(QuadraCGMidEpifaunaData3)

### Removing NAs
QuadraCGMidEpifaunaData4=(QuadraCGMidEpifaunaData3[complete.cases(QuadraCGMidEpifaunaData3), ])
str(QuadraCGMidEpifaunaData4)

QuadraCGMidEpifaunaData4$Area
QuadraCGMidEpifaunaData4$Site
QuadraCGMidEpifaunaData4$Quad
QuadraCGMidEpifaunaData4 %>% count(Site, Year)  # 4-5 Quadrats per mid zone

## Diversity Profile
names(QuadraCGMidEpifaunaData4)
QuadraCGMidEpifaunaData5=aggregate(.~SiteType+Year+Site, QuadraCGMidEpifaunaData4, mean)
str(QuadraCGMidEpifaunaData5)
QuadraCGMidEpifaunaData6=aggregate(.~SiteType+Site, QuadraCGMidEpifaunaData5, mean)
QuadraCGMidEpifaunaData7=aggregate(.~SiteType, QuadraCGMidEpifaunaData6, mean)
names(QuadraCGMidEpifaunaData7)

### Splitting by SiteType
QuadraMidEpifaunaCG= QuadraCGMidEpifaunaData7[1,]
QuadraMidEpifaunaCGRef= QuadraCGMidEpifaunaData7[2,]

### Removing zero's
QuadraMidEpifaunaCG2  =  QuadraMidEpifaunaCG[, colSums(QuadraMidEpifaunaCG != 0) > 0]
QuadraMidEpifaunaCGRef2  =  QuadraMidEpifaunaCGRef[, colSums(QuadraMidEpifaunaCGRef != 0) > 0]

names(QuadraMidEpifaunaCG2)
names(QuadraMidEpifaunaCGRef2)

QuadraMidEpifaunaCG3 = QuadraMidEpifaunaCG2[,8:26]
QuadraMidEpifaunaCGRef3 = QuadraMidEpifaunaCGRef2[,8:27]

### Equation
divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

QuadraMidEpifaunaCGRef4=divprof(QuadraMidEpifaunaCGRef3)

### Checking Adult Control Data
exp(diversity(QuadraMidEpifaunaCGRef3, index="shannon"))
###   2.4258
1 / (1-(diversity(QuadraMidEpifaunaCGRef3, index="simpson")))
###  1.5482

QuadraMidEpifaunaCG4=divprof(QuadraMidEpifaunaCG3)

### Checking Adult Treatment Data
exp(diversity(QuadraMidEpifaunaCG3, index="shannon"))
###   4.3068
1 / (1-(diversity(QuadraMidEpifaunaCG3, index="simpson")))
###  2.8063

## Building the plot
plot(QuadraMidEpifaunaCGRef4[,1],QuadraMidEpifaunaCGRef4[,2], col=ReferenceGardenColor,cex=2, pch=16,xlab="",ylab="",ylim = c(0,21))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(QuadraMidEpifaunaCG4[,1],QuadraMidEpifaunaCG4[,2],pch=16,cex=2,col= GardenColor)
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

## 
points(QuadraMidEpifaunaCG4[1,1],QuadraMidEpifaunaCG4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaCGRef4[1,1],QuadraMidEpifaunaCGRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.5,19.08, cex=cexText,expression( "Richness"["CG"] ))
text(0.5,20.08, cex=cexText,expression( "Richness"["GR"] ))

### Now plotting Control shannon and Simpson
points(QuadraMidEpifaunaCGRef4[10,1],QuadraMidEpifaunaCGRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaCGRef4[19,1],QuadraMidEpifaunaCGRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(QuadraMidEpifaunaCG4[10,1],QuadraMidEpifaunaCG4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaCG4[19,1],QuadraMidEpifaunaCG4[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.45,1.4,cex=cexText, expression( "Shannon"["GR"] ))
text(2.4,0.9, cex=cexText,expression( "Simpson"["GR"] ))

text(1.45,4.7,cex=cexText, expression( "Shannon"["CG"] ))
text(2.4,3.3,cex=cexText, expression( "Simpson"["CG"] ))

### Plot it
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "QuadraGardensEpifaunaOrderPlot5.jpeg", width = 22, height = 18, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(QuadraMidEpifaunaCGRef4[,1],QuadraMidEpifaunaCGRef4[,2], col=ReferenceGardenColor,cex=2, pch=16,xlab="",ylab="",ylim = c(0,21))

title(line = 2.5,xlab = list("Order q", cex = 1.3))
title(line = 2.5,ylab = list("Diversity", cex = 1.3))

points(QuadraMidEpifaunaCG4[,1],QuadraMidEpifaunaCG4[,2],pch=16,cex=2,col= GardenColor)

points(QuadraMidEpifaunaCG4[1,1],QuadraMidEpifaunaCG4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaCGRef4[1,1],QuadraMidEpifaunaCGRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.49,19.08, cex=cexText,expression( "Richness"["CG"] ))
text(0.49,20.08, cex=cexText,expression( "Richness"["GR"] ))

points(QuadraMidEpifaunaCGRef4[10,1],QuadraMidEpifaunaCGRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaCGRef4[19,1],QuadraMidEpifaunaCGRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(QuadraMidEpifaunaCG4[10,1],QuadraMidEpifaunaCG4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaCG4[19,1],QuadraMidEpifaunaCG4[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.36,1, cex=cexText,expression( "Shannon"["GR"] ))
text(2.37,0.5, cex=cexText,expression( "Simpson"["GR"] ))

text(1.36,5.1, cex=cexText,expression( "Shannon"["CG"] ))
text(2.37,3.65, cex=cexText,expression( "Simpson"["CG"] ))

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


QuadraCGEpi %<a-% {
  plot(QuadraMidEpifaunaCGRef4[,1],QuadraMidEpifaunaCGRef4[,2], col=ReferenceGardenColor,cex.axis=cex.axis1,cex=cex1, pch=16,xlab="",ylab="",ylim = c(0,21))
  
  title(line = 2.5,xlab = list("Order q", cex = cex2))
  title(line = 2.5,ylab = list("Diversity", cex = cex2))
  
  mtext(side=3, line=0.5,padj=-0.2,adj=0, "C)", col="Black", font=1, cex=cex3)
  
  points(QuadraMidEpifaunaCG4[,1],QuadraMidEpifaunaCG4[,2],pch=16,cex=2,col= GardenColor)
  
  points(QuadraMidEpifaunaCG4[1,1],QuadraMidEpifaunaCG4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  points(QuadraMidEpifaunaCGRef4[1,1],QuadraMidEpifaunaCGRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  
  text(0.61,19.08, cex=cexText,expression( "Richness"["CG"] ))
  text(0.61,20.08, cex=cexText,expression( "Richness"["GR"] ))
  
  points(QuadraMidEpifaunaCGRef4[10,1],QuadraMidEpifaunaCGRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(QuadraMidEpifaunaCGRef4[19,1],QuadraMidEpifaunaCGRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  points(QuadraMidEpifaunaCG4[10,1],QuadraMidEpifaunaCG4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(QuadraMidEpifaunaCG4[19,1],QuadraMidEpifaunaCG4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  text(1.36,1, cex=cexText,expression( "Shannon"["GR"] ))
  text(2.47,0.48, cex=cexText,expression( "Simpson"["GR"] ))
  
  text(1.47,5.1, cex=cexText,expression( "Shannon"["CG"] ))
  text(2.47,3.65, cex=cexText,expression( "Simpson"["CG"] ))}

##### CG-Q Done

#### Diversity Profiles,Epifaunal, SF-Q-Epi - Shellfish Farms Quadra Island ####
QuadraMidEpifaunaData <- subset(EpifaunaSpeciesData2, EpifaunaSpeciesData2$Area == "Quadra")
#  view(QuadraMidEpifaunaData)
QuadraSFMidEpifaunaData2 <- subset(QuadraMidEpifaunaData, QuadraMidEpifaunaData$Comparison == "Aquaculture")
# view(QuadraSFMidEpifaunaData2)
str(QuadraSFMidEpifaunaData2)
names(QuadraSFMidEpifaunaData2)

### Remove zero columns
QuadraSFMidEpifaunaData3=QuadraSFMidEpifaunaData2[, colSums(QuadraSFMidEpifaunaData2 != 0) > 0]
str(QuadraSFMidEpifaunaData3)
names(QuadraSFMidEpifaunaData3)

### Removing NAs
QuadraSFMidEpifaunaData4=(QuadraSFMidEpifaunaData3[complete.cases(QuadraSFMidEpifaunaData3), ])
str(QuadraSFMidEpifaunaData4)
names(QuadraSFMidEpifaunaData4)

QuadraSFMidEpifaunaData4$Area
QuadraSFMidEpifaunaData4$Site
QuadraSFMidEpifaunaData4$Quad
QuadraSFMidEpifaunaData4 %>% count(SiteType, Site, Year)  # 4-5 Quadrats per mid zone
QuadraSFMidEpifaunaData4 %>% count(SiteType, Site)  

## Diversity Profile
names(QuadraSFMidEpifaunaData4)
QuadraSFMidEpifaunaData5=aggregate(.~SiteType+Year+Site, QuadraSFMidEpifaunaData4, mean)
str(QuadraSFMidEpifaunaData5)
QuadraSFMidEpifaunaData6=aggregate(.~SiteType+Site, QuadraSFMidEpifaunaData5, mean)
QuadraSFMidEpifaunaData7=aggregate(.~SiteType, QuadraSFMidEpifaunaData6, mean)
names(QuadraSFMidEpifaunaData7)

### Splitting by SiteType
QuadraMidEpifaunaSF= QuadraSFMidEpifaunaData7[1,]
QuadraMidEpifaunaSFRef= QuadraSFMidEpifaunaData7[2,]

### Removing zero's
QuadraMidEpifaunaSF2  =  QuadraMidEpifaunaSF[, colSums(QuadraMidEpifaunaSF != 0) > 0]
QuadraMidEpifaunaSFRef2  =  QuadraMidEpifaunaSFRef[, colSums(QuadraMidEpifaunaSFRef != 0) > 0]

names(QuadraMidEpifaunaSF2)
names(QuadraMidEpifaunaSFRef2)

QuadraMidEpifaunaSF3 = QuadraMidEpifaunaSF2[,8:25]
QuadraMidEpifaunaSFRef3 = QuadraMidEpifaunaSFRef2[,8:28]

### Equation
divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

QuadraMidEpifaunaSFRef4=divprof(QuadraMidEpifaunaSFRef3)

### Checking Adult Control Data
exp(diversity(QuadraMidEpifaunaSFRef3, index="shannon"))
###   3.0678
1 / (1-(diversity(QuadraMidEpifaunaSFRef3, index="simpson")))
###  2.1119

QuadraMidEpifaunaSF4=divprof(QuadraMidEpifaunaSF3)

### Checking Adult Treatment Data
exp(diversity(QuadraMidEpifaunaSF3, index="shannon"))
###   3.3991
1 / (1-(diversity(QuadraMidEpifaunaSF3, index="simpson")))
###  2.5667

## Building the plot
plot(QuadraMidEpifaunaSFRef4[,1],QuadraMidEpifaunaSFRef4[,2], col=ReferenceFarmColor,cex=2, pch=16,xlab="",ylab="",ylim = c(1,22))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(QuadraMidEpifaunaSF4[,1],QuadraMidEpifaunaSF4[,2],pch=16,cex=2,col= FarmColor)
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

## 
points(QuadraMidEpifaunaSFRef4[1,1],QuadraMidEpifaunaSFRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaSF4[1,1],QuadraMidEpifaunaSF4[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.5,18.08, cex=cexText,expression( "Richness"["SF"] ))
text(0.5,21.08, cex=cexText,expression( "Richness"["FR"] ))

### Now plotting Control shannon and Simpson
points(QuadraMidEpifaunaSFRef4[10,1],QuadraMidEpifaunaSFRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaSFRef4[19,1],QuadraMidEpifaunaSFRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(QuadraMidEpifaunaSF4[10,1],QuadraMidEpifaunaSF4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaSF4[19,1],QuadraMidEpifaunaSF4[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.45,1.4,cex=cexText, expression( "Shannon"["FR"] ))
text(2.4,0.9, cex=cexText,expression( "Simpson"["FR"] ))

text(1.45,4.7, cex=cexText,expression( "Shannon"["SF"] ))
text(2.4,3.3, cex=cexText,expression( "Simpson"["SF"] ))

### Plot it
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "QuadraFarmsEpifaunaOrderPlot2.jpeg", width = 22, height = 18, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(QuadraMidEpifaunaSFRef4[,1],QuadraMidEpifaunaSFRef4[,2], col=ReferenceGardenColor,cex=2, pch=16,xlab="",ylab="",ylim = c(1,22))

title(line = 2.5,xlab = list("Order q", cex = 1.3))
title(line = 2.5,ylab = list("Diversity", cex = 1.3))

points(QuadraMidEpifaunaSF4[,1],QuadraMidEpifaunaSF4[,2],pch=16,cex=2,col= FarmColor)

points(QuadraMidEpifaunaSFRef4[1,1],QuadraMidEpifaunaSFRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaSF4[1,1],QuadraMidEpifaunaSF4[1,2],col= alpha("black",0.9),pch=16,cex=2)

QuadraMidEpifaunaSF4

text(0.49,18.08, cex=cexText,expression( "Richness"["SF"] ))
text(0.49,21.08, cex=cexText,expression( "Richness"["FR"] ))

points(QuadraMidEpifaunaSFRef4[10,1],QuadraMidEpifaunaSFRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaSFRef4[19,1],QuadraMidEpifaunaSFRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(QuadraMidEpifaunaSF4[10,1],QuadraMidEpifaunaSF4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidEpifaunaSF4[19,1],QuadraMidEpifaunaSF4[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.41,1.6, cex=cexText,expression( "Shannon"["FR"] ))
text(2.4,1, cex=cexText,expression( "Simpson"["FR"] ))

text(1.41,4.2, cex=cexText,expression( "Shannon"["SF"] ))
text(2.4,3.4, cex=cexText,expression( "Simpson"["SF"] ))

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


QuadraSFEpi %<a-% {
  plot(QuadraMidEpifaunaSFRef4[,1],QuadraMidEpifaunaSFRef4[,2], col=ReferenceFarmColor,cex.axis=cex.axis1,cex=cex1, pch=16,xlab="",ylab="",ylim = c(1,22))
  
  title(line = 2.5,xlab = list("Order q", cex = cex2))
  title(line = 2.5,ylab = list("Diversity", cex = cex2))
  
  mtext(side=3, line=0.5,padj=-0.2,adj=0, "D)", col="Black", font=1, cex=cex3)
  
  points(QuadraMidEpifaunaSF4[,1],QuadraMidEpifaunaSF4[,2],pch=16,cex=2,col= FarmColor)
  
  points(QuadraMidEpifaunaSFRef4[1,1],QuadraMidEpifaunaSFRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  points(QuadraMidEpifaunaSF4[1,1],QuadraMidEpifaunaSF4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  
  QuadraMidEpifaunaSF4
  
  text(0.61,18.08, cex=cexText,expression( "Richness"["SF"] ))
  text(0.61,21.08, cex=cexText,expression( "Richness"["FR"] ))
  
  points(QuadraMidEpifaunaSFRef4[10,1],QuadraMidEpifaunaSFRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(QuadraMidEpifaunaSFRef4[19,1],QuadraMidEpifaunaSFRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  points(QuadraMidEpifaunaSF4[10,1],QuadraMidEpifaunaSF4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(QuadraMidEpifaunaSF4[19,1],QuadraMidEpifaunaSF4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  text(1.39,1.55, cex=cexText,expression( "Shannon"["FR"] ))
  text(1.39,4.3, cex=cexText,expression( "Shannon"["SF"] ))

  text(2.5,1, cex=cexText,expression( "Simpson"["FR"] ))
  text(2.5,3.4, cex=cexText,expression( "Simpson"["SF"] ))}


#SF-Q Done
####

#### Diversity Profiles,Epifaunal, SF-B-Epi - Shellfish Farms Baynes Sound ####
BaynesMidEpifaunaData <- subset(EpifaunaSpeciesData2, EpifaunaSpeciesData2$Area == "Baynes")
# view(BaynesMidEpifaunaData)
BaynesSFMidEpifaunaData2 <- subset(BaynesMidEpifaunaData, BaynesMidEpifaunaData$Comparison == "Aquaculture")
# view(BaynesSFMidEpifaunaData2)
str(BaynesSFMidEpifaunaData2)
names(BaynesSFMidEpifaunaData2)

### Remove zero columns
BaynesSFMidEpifaunaData3=BaynesSFMidEpifaunaData2[, colSums(BaynesSFMidEpifaunaData2 != 0) > 0]
str(BaynesSFMidEpifaunaData3)
names(BaynesSFMidEpifaunaData3)

### Removing NAs
BaynesSFMidEpifaunaData4=(BaynesSFMidEpifaunaData3[complete.cases(BaynesSFMidEpifaunaData3), ])
str(BaynesSFMidEpifaunaData4)
names(BaynesSFMidEpifaunaData4)

BaynesSFMidEpifaunaData4$Area
BaynesSFMidEpifaunaData4$Site
BaynesSFMidEpifaunaData4$Quad
BaynesSFMidEpifaunaData4 %>% count(SiteType, Site, Year)  # 5 Quadrats per mid zone
BaynesSFMidEpifaunaData4 %>% count(SiteType, Site)  

## Diversity Profile
names(BaynesSFMidEpifaunaData4)
BaynesSFMidEpifaunaData5=aggregate(.~SiteType+Year+Site, BaynesSFMidEpifaunaData4, mean)
str(BaynesSFMidEpifaunaData5)
BaynesSFMidEpifaunaData6=aggregate(.~SiteType+Site, BaynesSFMidEpifaunaData5, mean)
BaynesSFMidEpifaunaData7=aggregate(.~SiteType, BaynesSFMidEpifaunaData6, mean)
names(BaynesSFMidEpifaunaData7)

### Splitting by SiteType
BaynesMidEpifaunaSF= BaynesSFMidEpifaunaData7[1,]
BaynesMidEpifaunaSFRef= BaynesSFMidEpifaunaData7[2,]

### Removing zero's
BaynesMidEpifaunaSF2  =  BaynesMidEpifaunaSF[, colSums(BaynesMidEpifaunaSF != 0) > 0]
BaynesMidEpifaunaSFRef2  =  BaynesMidEpifaunaSFRef[, colSums(BaynesMidEpifaunaSFRef != 0) > 0]

names(BaynesMidEpifaunaSF2)
names(BaynesMidEpifaunaSFRef2)

BaynesMidEpifaunaSF3 = BaynesMidEpifaunaSF2[,8:32]
BaynesMidEpifaunaSFRef3 = BaynesMidEpifaunaSFRef2[,8:24]

### Equation
divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

BaynesMidEpifaunaSFRef4=divprof(BaynesMidEpifaunaSFRef3)

### Checking Adult Control Data
exp(diversity(BaynesMidEpifaunaSFRef3, index="shannon"))
###   3.1326
1 / (1-(diversity(BaynesMidEpifaunaSFRef3, index="simpson")))
###  1.8479

BaynesMidEpifaunaSF4=divprof(BaynesMidEpifaunaSF3)

### Checking Adult Treatment Data
exp(diversity(BaynesMidEpifaunaSF3, index="shannon"))
###   4.3831
1 / (1-(diversity(BaynesMidEpifaunaSF3, index="simpson")))
###  2.3645

### Plot it
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "BaynesFarmsEpifaunaOrderPlot2.jpeg", width = 22, height = 18, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(BaynesMidEpifaunaSF4[,1],BaynesMidEpifaunaSF4[,2], col=FarmColor,cex=2, pch=16,xlab="",ylab="",ylim = c(0,26))

title(line = 2.5,xlab = list("Order q", cex = 1.3))
title(line = 2.5,ylab = list("Diversity", cex = 1.3))

points(BaynesMidEpifaunaSFRef4[,1],BaynesMidEpifaunaSFRef4[,2],pch=16,cex=2,col= ReferenceFarmColor)

points(BaynesMidEpifaunaSFRef4[1,1],BaynesMidEpifaunaSFRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(BaynesMidEpifaunaSF4[1,1],BaynesMidEpifaunaSF4[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.48,25.08, cex=cexText,expression( "Richness"["SF"] ))
text(0.48,17.08, cex=cexText,expression( "Richness"["FR"] ))

points(BaynesMidEpifaunaSFRef4[10,1],BaynesMidEpifaunaSFRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(BaynesMidEpifaunaSFRef4[19,1],BaynesMidEpifaunaSFRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(BaynesMidEpifaunaSF4[10,1],BaynesMidEpifaunaSF4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(BaynesMidEpifaunaSF4[19,1],BaynesMidEpifaunaSF4[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.41,1.2,cex=cexText, expression( "Shannon"["FR"] ))
text(1.41,5.1, cex=cexText,expression( "Shannon"["SF"] ))

text(2.4,3.4, cex=cexText,expression( "Simpson"["SF"] ))
text(2.4,0.7, cex=cexText,expression( "Simpson"["FR"] ))


dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")

##
BaynesSFEpi %<a-% {
  plot(BaynesMidEpifaunaSF4[,1],BaynesMidEpifaunaSF4[,2], col=FarmColor,cex.axis=cex.axis1,cex=cex1, pch=16,xlab="",ylab="",ylim = c(0,26))
  
  title(line = 2.5,xlab = list("Order q", cex = cex2))
  title(line = 2.5,ylab = list("Diversity", cex = cex2))
  
  mtext(side=3, line=0.5,padj=-0.2,adj=0, "E)", col="Black", font=1, cex=cex3)
  
  points(BaynesMidEpifaunaSFRef4[,1],BaynesMidEpifaunaSFRef4[,2],pch=16,cex=2,col= ReferenceFarmColor)
  
  points(BaynesMidEpifaunaSFRef4[1,1],BaynesMidEpifaunaSFRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  points(BaynesMidEpifaunaSF4[1,1],BaynesMidEpifaunaSF4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  
  text(0.61,25.08, cex=cexText,expression( "Richness"["SF"] ))
  text(0.61,17.08,cex=cexText, expression( "Richness"["FR"] ))
  
  points(BaynesMidEpifaunaSFRef4[10,1],BaynesMidEpifaunaSFRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(BaynesMidEpifaunaSFRef4[19,1],BaynesMidEpifaunaSFRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  points(BaynesMidEpifaunaSF4[10,1],BaynesMidEpifaunaSF4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(BaynesMidEpifaunaSF4[19,1],BaynesMidEpifaunaSF4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  text(1.41,1.15, cex=cexText,expression( "Shannon"["FR"] ))
  text(1.49,5.1, cex=cexText,expression( "Shannon"["SF"] ))
  
  text(2.5,3.4, cex=cexText,expression( "Simpson"["SF"] ))
  text(2.5,0.62, cex=cexText,expression( "Simpson"["FR"] ))}

### SF-B Done

#### Figure 2: Combining Diversity Profiles ####
par(mfrow = c(1, 4))
par(cex = 0.6)
CalvertCGEpi
QuadraCGEpi
QuadraSFEpi
BaynesSFEpi

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "CombinedPlot2.jpeg", width = 35, height = 10, units = "cm", pointsize = 15, quality = 100, res = 800)  

par(mfrow = c(1, 4))
par(cex = 0.45)

CalvertCGEpi
QuadraCGEpi
QuadraSFEpi
BaynesSFEpi

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")

par(mfrow = c(1, 1))

## Epifaunal diversity summary  (across sites/regions)
EpifaunaSpeciesData2
names(EpifaunaSpeciesData2)
EpifaunaSpeciesData2$Area
EpifaunaSpeciesData2$Comparison
EpifaunaSpeciesData2$Site
EpifaunaSpeciesData2$Quad
EpifaunaSpeciesData2 %>% count(Comparison,Area,SiteType, Site, Year)  # 3-5 Quadrats per mid zone
EpifaunaSpeciesData2 %>% count(Comparison,Area,SiteType, Site)  # 3-5 Quadrats per mid zone
EpifaunaSpeciesData2 %>% count(Comparison,SiteType)  # 3-5 Quadrats per mid zone

EpifaunaSpeciesData2 %>% count(Comparison,SiteType, Site)  

## Diversity Profile

EpifaunaSpeciesData3=aggregate(.~Comparison+Area+SiteType+Site+Year, EpifaunaSpeciesData2, mean)
EpifaunaSpeciesData4=aggregate(.~Comparison+Area+SiteType+Site, EpifaunaSpeciesData3, mean)
EpifaunaSpeciesData5=aggregate(.~Comparison+SiteType, EpifaunaSpeciesData4, mean)

### Splitting by SiteType
EpifaunaSpeciesData5
names(EpifaunaSpeciesData5)
GardenEpifauna= EpifaunaSpeciesData5[1,]
AquacultureEpifauna= EpifaunaSpeciesData5[2,]
AquaReferenceEpifauna= EpifaunaSpeciesData5[3,]
GardenReferenceEpifauna= EpifaunaSpeciesData5[4,]

### Removing zero's
GardenReferenceEpifauna2  =  GardenReferenceEpifauna[, colSums(GardenReferenceEpifauna != 0) > 0]
AquaReferenceEpifauna2  =  AquaReferenceEpifauna[, colSums(AquaReferenceEpifauna != 0) > 0]
AquacultureEpifauna2  =  AquacultureEpifauna[, colSums(AquacultureEpifauna != 0) > 0]
GardenEpifauna2  =  GardenEpifauna[, colSums(GardenEpifauna != 0) > 0]

names(GardenReferenceEpifauna2)
names(AquaReferenceEpifauna2)
names(AquacultureEpifauna2)
names(GardenEpifauna2)

GardenEpifauna3 = GardenEpifauna2[,8:42]
AquacultureEpifauna3= AquacultureEpifauna2[,8:35]
AquaReferenceEpifauna3=AquaReferenceEpifauna2 [,8:33]
GardenReferenceEpifauna3= GardenReferenceEpifauna2[,8:34]

### Equation
divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

GardenEpifauna4=divprof(GardenEpifauna3)
AquacultureEpifauna4=divprof(AquacultureEpifauna3)
AquaReferenceEpifauna4=divprof(AquaReferenceEpifauna3)
GardenReferenceEpifauna4=divprof(GardenReferenceEpifauna3)

ReferenceFarmColor
ReferenceGardenColor
FarmColor
GardenColor

### Plot it
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "EpiDiversitySummaryPlot8.jpeg", width = 22, height = 18, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(GardenEpifauna4[,1],GardenEpifauna4[,2], col=GardenColor,cex.axis=cex.axis1,cex=cex1, pch=16,xlab="",ylab="",ylim = c(-0.75,35))
title(line = 2.5,xlab = list("Order q", cex = cex2))
title(line = 2.5,ylab = list("Diversity", cex = cex2))

mtext(side=3, line=0.5,padj=-0.2,adj=0, "A)", col="Black", font=1, cex=cex3)

points(AquacultureEpifauna4[,1],AquacultureEpifauna4[,2],pch=16,cex=2,col= FarmColor)
points(AquaReferenceEpifauna4[,1],AquaReferenceEpifauna4[,2],pch=16,cex=2,col= ReferenceFarmColor)
points(GardenReferenceEpifauna4[,1],GardenReferenceEpifauna4[,2],pch=16,cex=2,col= ReferenceGardenColor)

points(GardenEpifauna4[1,1],GardenEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(AquacultureEpifauna4[1,1],AquacultureEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(AquaReferenceEpifauna4[1,1],AquaReferenceEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(GardenReferenceEpifauna4[1,1],GardenReferenceEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.48,35, cex=cexText,expression( "Richness"["CG"] ))
text(0.48,27, cex=cexText,expression( "Richness"["GR"] ))
text(0.48,25.7, cex=cexText,expression( "Richness"["FR"]))
text(0.48,28.3, cex=cexText,expression( "Richness"["SF"]))

points(GardenEpifauna4[10,1],GardenEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(GardenEpifauna4[19,1],GardenEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(AquacultureEpifauna4[10,1],AquacultureEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(AquacultureEpifauna4[19,1],AquacultureEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(AquaReferenceEpifauna4[10,1],AquaReferenceEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(AquaReferenceEpifauna4[19,1],AquaReferenceEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(GardenReferenceEpifauna4[10,1],GardenReferenceEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(GardenReferenceEpifauna4[19,1],GardenReferenceEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.45,7.1, cex=cexText,expression(  "Shannon"["SF"]))
text(1.45,5.75,cex=cexText, expression(  "Shannon"["CG"]))
text(0.57,1.93,cex=cexText, expression(  "Shannon"["FR"]))
text(1.35,0.6, cex=cexText,expression(  "Shannon"["GR"]))

text(2.4,6.2, cex=cexText,expression(  "Simpson"["FR"]))
text(2.4,4.8, cex=cexText,expression(  "Simpson"["CG"]))
text(2.4,0.2, cex=cexText,expression( "Simpson"["SF"]))
text(2.4,-1.2,cex=cexText, expression( "Simpson"["GR"]))

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")
###

#### Summary Epifauna Done, Venn Diagram ####
library(VennDiagram)
library(jpeg)

GardenEpifauna3 = GardenEpifauna2[,8:42]
AquacultureEpifauna3= AquacultureEpifauna2[,8:35]
AquaReferenceEpifauna3=AquaReferenceEpifauna2 [,8:33]
GardenReferenceEpifauna3= GardenReferenceEpifauna2[,8:34]


GardenEpifaunaList=(colnames(GardenEpifauna3))
AquacultureEpifaunaList=(colnames(AquacultureEpifauna3))
AquaReferenceEpifaunaList=(colnames(AquaReferenceEpifauna3))
GardenReferenceEpifaunaList=(colnames(GardenReferenceEpifauna3))


### Venn Diagram Species Lists, export and combine in a table
EpifaunaSpeciesData5List  =  EpifaunaSpeciesData5[, colSums(EpifaunaSpeciesData5 != 0) > 0]
EpifaunaSpeciesList=(colnames(EpifaunaSpeciesData5List))
write.csv(GardenEpifaunaList,"GardenEpifaunaList.csv")
write.csv(AquacultureEpifaunaList,"AquacultureEpifaunaList.csv")
write.csv(AquaReferenceEpifaunaList,"AquaReferenceEpifaunaList.csv")
write.csv(GardenReferenceEpifaunaList,"GardenReferenceEpifaunaList.csv")
write.csv(EpifaunaSpeciesList,"EpifaunaSpeciesList.csv")


GardenColor=("#7CAE00")   #### Green 
FarmColor=("#0d7dd9")     #### Blue 
ReferenceFarmColor=("#ff0000") #### Red
ReferenceGardenColor=("#fa7f2d") #### Orange

colors= c(GardenColor, ReferenceGardenColor, FarmColor,ReferenceFarmColor)

venn.diagram(x = list(GardenEpifaunaList,GardenReferenceEpifaunaList,AquacultureEpifaunaList,AquaReferenceEpifaunaList),
             category.names = c("Clam Gardens ","Garden Reference","Shellfish Farm","Farm Reference"),
             filename = "VennDiagramFinal.jpeg",
             col = "black",
             lwd = 2,
             cex = 3,
             height = 8000, width = 8000,
             fill = colors,
             label.col= "black",
             fontfamily = "serif",
             cat.cex = 3,
             cat.fontfamily = "serif",
             margin=0.12,
             resolution=500
)

VennDiagram <-readJPEG("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/VennDiagramFinal.jpg")
VennDiagram
plot(1:10,ty="n")
rasterImage(VennDiagram,6,7,8,8)
plotImage(VennDiagram)

#### Summary Epifauna Done ####

SummaryEpifauna %<a-% {
  plot(GardenEpifauna4[,1],GardenEpifauna4[,2], col=GardenColor,cex.axis=cex.axis1,cex=cex1, pch=16,xlab="",ylab="",ylim = c(-0.75,35))
  title(line = 2.5,xlab = list("Order q", cex = cex2))
  title(line = 2.5,ylab = list("Diversity", cex = cex2))
  
  points(AquacultureEpifauna4[,1],AquacultureEpifauna4[,2],pch=16,cex=2,col= FarmColor)
  points(AquaReferenceEpifauna4[,1],AquaReferenceEpifauna4[,2],pch=16,cex=2,col= ReferenceFarmColor)
  points(GardenReferenceEpifauna4[,1],GardenReferenceEpifauna4[,2],pch=16,cex=2,col= ReferenceGardenColor)
  
  points(GardenEpifauna4[1,1],GardenEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  points(AquacultureEpifauna4[1,1],AquacultureEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  points(AquaReferenceEpifauna4[1,1],AquaReferenceEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  points(GardenReferenceEpifauna4[1,1],GardenReferenceEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  
  mtext(side=3, line=0.5,padj=-0.2,adj=0, "A)", col="Black", font=1, cex=cex3)
  text(0.28,35, cex=cexText,expression( "Richness"["CG"] ))
  text(0.28,27, cex=cexText,expression( "Richness"["GR"] ))
  text(0.28,25.7, cex=cexText,expression( "Richness"["FR"]))
  text(0.28,28.3, cex=cexText,expression( "Richness"["SF"]))
  
  points(GardenEpifauna4[10,1],GardenEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(GardenEpifauna4[19,1],GardenEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  points(AquacultureEpifauna4[10,1],AquacultureEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(AquacultureEpifauna4[19,1],AquacultureEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  points(AquaReferenceEpifauna4[10,1],AquaReferenceEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(AquaReferenceEpifauna4[19,1],AquaReferenceEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  points(GardenReferenceEpifauna4[10,1],GardenReferenceEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(GardenReferenceEpifauna4[19,1],GardenReferenceEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  text(1.23,6.8, cex=cexText,expression(  "Shannon"["SF"]))
  text(1.23,5.63, cex=cexText,expression(  "Shannon"["CG"]))
  text(0.77,1.83, cex=cexText,expression(  "Shannon"["FR"]))
  text(1.2,0.62, cex=cexText,expression(  "Shannon"["GR"]))
  
  text(2.2,5.57,cex=cexText, expression(  "Simpson"["FR"]))
  text(2.2,4.43,cex=cexText, expression(  "Simpson"["CG"]))
  text(2.2,0.42, cex=cexText,expression( "Simpson"["SF"]))
  text(2.2,-0.6, cex=cexText,expression( "Simpson"["GR"]))}


#### Test 
SummaryEpifauna %<a-% {
  plot(GardenEpifauna4[,1],GardenEpifauna4[,2], col=GardenColor,cex.axis=cex.axis1,cex=cex1, pch=16,xlab="",ylab="",ylim = c(-0.75,35))
  title(line = 2.5,xlab = list("Order q", cex = cex2))
  title(line = 2.5,ylab = list("Diversity", cex = cex2))
  
  points(AquacultureEpifauna4[,1],AquacultureEpifauna4[,2],pch=16,cex=2,col= FarmColor)
  points(AquaReferenceEpifauna4[,1],AquaReferenceEpifauna4[,2],pch=16,cex=2,col= ReferenceFarmColor)
  points(GardenReferenceEpifauna4[,1],GardenReferenceEpifauna4[,2],pch=16,cex=2,col= ReferenceGardenColor)
  
  points(GardenEpifauna4[1,1],GardenEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  points(AquacultureEpifauna4[1,1],AquacultureEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  points(AquaReferenceEpifauna4[1,1],AquaReferenceEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  points(GardenReferenceEpifauna4[1,1],GardenReferenceEpifauna4[1,2],col= alpha("black",0.9),pch=16,cex=2)
  
  mtext(side=3, line=0.5,padj=-0.2,adj=0, "A)", col="Black", font=1, cex=cex3)
  text(0.28,35, cex=cexText,expression( "Richness"["CG"] ))
  text(0.28,27, cex=cexText,expression( "Richness"["GR"] ))
  text(0.28,25.7, cex=cexText,expression( "Richness"["FR"]))
  text(0.28,28.3, cex=cexText,expression( "Richness"["SF"]))
  
  points(GardenEpifauna4[10,1],GardenEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(GardenEpifauna4[19,1],GardenEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  points(AquacultureEpifauna4[10,1],AquacultureEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(AquacultureEpifauna4[19,1],AquacultureEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  points(AquaReferenceEpifauna4[10,1],AquaReferenceEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(AquaReferenceEpifauna4[19,1],AquaReferenceEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  points(GardenReferenceEpifauna4[10,1],GardenReferenceEpifauna4[10,2],col= alpha("black",0.9),pch=16,cex=2)
  points(GardenReferenceEpifauna4[19,1],GardenReferenceEpifauna4[19,2],col= alpha("black",0.9),pch=16,cex=2)
  
  text(1.23,6.8, cex=cexText,expression(  "Shannon"["SF"]))
  text(1.23,5.63, cex=cexText,expression(  "Shannon"["CG"]))
  text(0.77,1.83, cex=cexText,expression(  "Shannon"["FR"]))
  text(1.2,0.62, cex=cexText,expression(  "Shannon"["GR"]))
  
  text(2.2,5.57,cex=cexText, expression(  "Simpson"["FR"]))
  text(2.2,4.43,cex=cexText, expression(  "Simpson"["CG"]))
  text(2.2,0.42, cex=cexText,expression( "Simpson"["SF"]))
  text(2.2,-0.6, cex=cexText,expression( "Simpson"["GR"]))
  rasterImage(VennDiagram,2.5,12, 5,35)
  text(2.8,31, cex=cex3,expression( "F)"))}


#### Figure 2: Final Combining Diversity Profiles ####

### Two legend options
# Legend1 <-loadImage(("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/LegendV8.jpg"), lower = c(0.97, 0.97, .97), upper = c(1, 1, 1),hsv = FALSE)
Legend2 <-loadImage(("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/LegendV9.jpg"), lower = c(0.97, 0.97, .97), upper = c(1, 1, 1),hsv = FALSE)
plotImage(Legend2)


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "CombinedPlotRevisedFINAL.jpeg", width = 70, height = 50, units = "cm", pointsize = 15, quality = 100, res = 1200)  
layout(matrix(c(1, 1, 2,
                1, 1, 3,
                4, 5, 6), nrow=3, byrow=TRUE),widths=c(1,1,1), heights=c(1,1,1))
par(cex = 1.03)

par(mai = c(1, 1.1, 0.8, 0.7))
SummaryEpifauna
par(mai = c(1, 1, 0.8, 0.9))
CalvertCGEpi
QuadraCGEpi
par(mai = c(1.1, 1.1, 0.7, 0.7))
QuadraSFEpi
BaynesSFEpi
par(cex = 1.2)
par(mai = c(0.5, 0.5, 0.5, 0.5))
plotImage(Legend2)

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")



#### Figure 3: nMDS  ####

## To Make
# 1 nMDS Epifauna- Types/Sites plotted by Comparison, plus Summary (all together)
### Make this, if summary orinitation is the same, present that, and (maybe) add year together plot ###

# 1 nMDS Epifauna- Types/Region/Sites plotted by Comparison, plus Summary (all together)
#### Data Structure ####
### Epifauna Data
rm(list = ls(all = TRUE))

SFandCGDataMid2=read.csv("SFandCGDataMid2.csv", header=T)
SFandCGDataMid2
SFandCGDataMid2 %>% count(SiteType,SiteType2,Site, Year)

SFandCGDataMid2$SiteType2

names(SFandCGDataMid2)
MidEpifaunaData2 <- SFandCGDataMid2[,34:93]
names(MidEpifaunaData2)
MidDataStructure2=SFandCGDataMid2[,1:8]
names(MidDataStructure2)
EpifaunaSpeciesData3=cbind(MidDataStructure2,MidEpifaunaData2)

### Removing NAs
EpifaunaSpeciesData4=(EpifaunaSpeciesData3[complete.cases(EpifaunaSpeciesData3), ])
str(EpifaunaSpeciesData4)
names(EpifaunaSpeciesData4)
## Removing zero's
EpifaunaSpeciesData5=EpifaunaSpeciesData4[, colSums(EpifaunaSpeciesData4 != 0) > 0]
str(EpifaunaSpeciesData5)
names(EpifaunaSpeciesData5)

#### Epifauna nMDS ####

# 1 nMDS Epifauna- Types/Region/Sites plotted by Comparison, plus Summary (all together)
EpifaunaSpeciesData5
str(EpifaunaSpeciesData5)
names(EpifaunaSpeciesData5)
EpifaunaSpeciesData6=aggregate(.~Comparison+Area+SiteType+SiteType2+Year+Site, EpifaunaSpeciesData5, mean)
EpifaunaSpeciesData7=aggregate(.~Comparison+Area+SiteType+SiteType2+Site, EpifaunaSpeciesData6, mean)

str(EpifaunaSpeciesData7)
names(EpifaunaSpeciesData7)
EpifaunaSpeciesDataStructure = EpifaunaSpeciesData7[,1:8]
EpifaunaSpeciesData8   = EpifaunaSpeciesData7[,9:59]

### Dummy Variable
EpifaunaSpeciesData9=((EpifaunaSpeciesData8+1))

### nMDS
nMDSEpiCommunityMid<- metaMDS(EpifaunaSpeciesData9, distance = "bray", k = 2) 
nMDSEpiCommunityMid$stress  
stressplot(nMDSEpiCommunityMid)

#extract NMDS scores (x and y coordinates)
data.scoresEpi = as.data.frame(scores(nMDSEpiCommunityMid))
# Add data structure
names(EpifaunaSpeciesDataStructure)
data.scoresEpi$Comparison = EpifaunaSpeciesDataStructure$Comparison
data.scoresEpi$Site = EpifaunaSpeciesDataStructure$Site
data.scoresEpi$SiteType2 = EpifaunaSpeciesDataStructure$SiteType2
data.scoresEpi$SiteType = EpifaunaSpeciesDataStructure$SiteType
data.scoresEpi$Area = EpifaunaSpeciesDataStructure$Area
str(data.scoresEpi)
head(data.scoresEpi)

## Plot
GardenColor=("#7CAE00")   #### Green 
FarmColor=("#0d7dd9")     #### Blue 
ReferenceFarmColor=("#ff0000") #### Red
ReferenceGardenColor=("#fa7f2d") #### Orange

hjust=0.85
vjust=-0.1
size=6
stroke = 1.4
nMDSshapes=c(5,12)
AxisSize=18
LegendSize=18    
textsize=6
                           
EpifaunanMDS=ggplot(data.scoresEpi, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Area, colour = SiteType2))+ 
  theme(axis.text.y = element_text(colour = "black", size = AxisSize, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = AxisSize), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = LegendSize), 
        axis.title.x = element_text(face = "bold", size = LegendSize, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "", y = "NMDS2", shape = "")  + 
  scale_colour_manual(values = c("#7CAE00", "#ff0000","#fa7f2d","#0d7dd9"))
EpifaunanMDS

### Adding ordiellipse or connections, 2 options
### Option 1
grp.a <- data.scoresEpi[data.scoresEpi$SiteType2 == "Clam Garden", ][chull(data.scoresEpi[data.scoresEpi$SiteType2 == 
                                                                   "Clam Garden", c("NMDS1", "NMDS2")]), ]

grp.b <- data.scoresEpi[data.scoresEpi$SiteType2 == "Farm Reference", ][chull(data.scoresEpi[data.scoresEpi$SiteType2 == 
                                                                                            "Farm Reference", c("NMDS1", "NMDS2")]), ]
grp.c <- data.scoresEpi[data.scoresEpi$SiteType2 == "Garden Reference", ][chull(data.scoresEpi[data.scoresEpi$SiteType2 == 
                                                                                            "Garden Reference", c("NMDS1", "NMDS2")]), ]
grp.d <- data.scoresEpi[data.scoresEpi$SiteType2 == "Shellfish Farm", ][chull(data.scoresEpi[data.scoresEpi$SiteType2 == 
                                                                                            "Shellfish Farm", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b,grp.c,grp.d)  #combine
hull.data

EpiFaunanMDSPlotFinal=EpifaunanMDS+geom_polygon (data=hull.data,aes (x=NMDS1,y=NMDS2,fill=SiteType2,group=SiteType2),alpha=0.30)+
  scale_fill_manual(values = c("#7CAE00", "#ff0000","#fa7f2d","#0d7dd9"))+guides(fill = FALSE)+
  geom_text(x=0.215, y=0.159, label="Stress = 0.12",size=textsize)
EpiFaunanMDSPlotFinal

### without legend
EpiFaunanMDSPlotFinalnoLegend= EpifaunanMDS+geom_polygon (data=hull.data,aes (x=NMDS1,y=NMDS2,fill=SiteType2,group=SiteType2),alpha=0.30)+
  scale_fill_manual(values = c("#7CAE00", "#ff0000","#fa7f2d","#0d7dd9"))+
  geom_text(x=0.202, y=-0.133, label="Stress = 0.12",size=textsize)+guides(colour = F,shape=F,fill=F)
EpiFaunanMDSPlotFinalnoLegend


####  nMDS Plots ####
# nMDS Epifauna show Type, Site, Region, Year

### Data Structure
EpifaunaSpeciesData5

## Epifauna
CalvertEpifauna <- subset(EpifaunaSpeciesData5, EpifaunaSpeciesData5$Area == "Calvert")
QuadraEpifauna <- subset(EpifaunaSpeciesData5, EpifaunaSpeciesData5$Area == "Quadra")
BaynesEpifauna <- subset(EpifaunaSpeciesData5, EpifaunaSpeciesData5$Area == "Baynes")
QuadraEpifaunaGarden <- subset(QuadraEpifauna, QuadraEpifauna$Comparison == "Garden")
QuadraEpifaunaAquaculture <- subset(QuadraEpifauna, QuadraEpifauna$Comparison == "Aquaculture")

## Averaging 
CalvertEpifaunaAverage=aggregate(.~Comparison+Area+SiteType+SiteType2+Year+Site, CalvertEpifauna, mean)
BaynesEpifaunaAverage=aggregate(.~Comparison+Area+SiteType+SiteType2+Year+Site, BaynesEpifauna, mean)
QuadraEpifaunaGardenAverage=aggregate(.~Comparison+Area+SiteType+SiteType2+Year+Site, QuadraEpifaunaGarden, mean)
QuadraEpifaunaAquacultureAverage=aggregate(.~Comparison+Area+SiteType+SiteType2+Year+Site, QuadraEpifaunaAquaculture, mean)

## Removing zero's
CalvertEpifaunaAverage=CalvertEpifaunaAverage[, colSums(CalvertEpifaunaAverage != 0) > 0]
BaynesEpifaunaAverage=BaynesEpifaunaAverage[, colSums(BaynesEpifaunaAverage != 0) > 0]
QuadraEpifaunaGardenAverage=QuadraEpifaunaGardenAverage[, colSums(QuadraEpifaunaGardenAverage != 0) > 0]
QuadraEpifaunaAquacultureAverage=QuadraEpifaunaAquacultureAverage[, colSums(QuadraEpifaunaAquacultureAverage != 0) > 0]

## Splitting up, adding dummy variable, nMDS, stress plot
names(CalvertEpifaunaAverage)
CalvertEpifaunaAverageStructure = CalvertEpifaunaAverage[,1:8]
CalvertEpifaunaAverageData   = CalvertEpifaunaAverage[,9:41]
CalvertEpifaunaAverageData2=((CalvertEpifaunaAverageData)+1) 
### nMDS
nMDSCalvertEpifaunaMid<- metaMDS(CalvertEpifaunaAverageData2, distance = "bray", k = 2) 
nMDSCalvertEpifaunaMid$stress  
stressplot(nMDSCalvertEpifaunaMid)
data.scoresCalvertEpifauna = as.data.frame(scores(nMDSCalvertEpifaunaMid))
#extract NMDS scores (x and y coordinates)
# Add structure
names(CalvertEpifaunaAverageStructure)
data.scoresCalvertEpifauna$Comparison = CalvertEpifaunaAverageStructure$Comparison
data.scoresCalvertEpifauna$Site =       CalvertEpifaunaAverageStructure$Site
data.scoresCalvertEpifauna$SiteType2 =  CalvertEpifaunaAverageStructure$SiteType2
data.scoresCalvertEpifauna$SiteType =   CalvertEpifaunaAverageStructure$SiteType
data.scoresCalvertEpifauna$Area =       CalvertEpifaunaAverageStructure$Area
data.scoresCalvertEpifauna$Year =       CalvertEpifaunaAverageStructure$Year
data.scoresCalvertEpifauna$Year=as.factor(data.scoresCalvertEpifauna$Year)
str(data.scoresCalvertEpifauna)
head(data.scoresCalvertEpifauna)
## data.scoresCalvertEpifauna   Stress 0.0995

names(BaynesEpifaunaAverage)
BaynesEpifaunaAverageStructure = BaynesEpifaunaAverage[,1:8]
BaynesEpifaunaAverageData   = BaynesEpifaunaAverage[,9:35]
BaynesEpifaunaAverageData2=((BaynesEpifaunaAverageData)+1)
### nMDS
nMDSBaynesEpifaunaMid<- metaMDS(BaynesEpifaunaAverageData2, distance = "bray", k = 2) 
nMDSBaynesEpifaunaMid$stress  
stressplot(nMDSBaynesEpifaunaMid)
data.scoresBaynesEpifauna = as.data.frame(scores(nMDSBaynesEpifaunaMid))
# Add structure
names(BaynesEpifaunaAverageStructure)
data.scoresBaynesEpifauna$Comparison = BaynesEpifaunaAverageStructure$Comparison
data.scoresBaynesEpifauna$Site =       BaynesEpifaunaAverageStructure$Site
data.scoresBaynesEpifauna$SiteType2 =  BaynesEpifaunaAverageStructure$SiteType2
data.scoresBaynesEpifauna$SiteType =   BaynesEpifaunaAverageStructure$SiteType
data.scoresBaynesEpifauna$Area =       BaynesEpifaunaAverageStructure$Area
data.scoresBaynesEpifauna$Year =       BaynesEpifaunaAverageStructure$Year
data.scoresBaynesEpifauna$Year=as.factor(data.scoresBaynesEpifauna$Year)
str(data.scoresBaynesEpifauna)
head(data.scoresBaynesEpifauna)
## data.scoresBaynesEpifauna   Stress 0.080152

names(QuadraEpifaunaGardenAverage)
QuadraEpifaunaGardenAverageStructure = QuadraEpifaunaGardenAverage[,1:8]
QuadraEpifaunaGardenAverageData   = QuadraEpifaunaGardenAverage[,9:32]
QuadraEpifaunaGardenAverageData2=((QuadraEpifaunaGardenAverageData)+1)
### nMDS
nMDSQuadraGardenEpifaunaMid<- metaMDS(QuadraEpifaunaGardenAverageData2, distance = "bray", k = 2) 
nMDSQuadraGardenEpifaunaMid$stress  
stressplot(nMDSQuadraGardenEpifaunaMid)
data.scoresnQuadraGardenEpifauna = as.data.frame(scores(nMDSQuadraGardenEpifaunaMid))
# Add structure
names(QuadraEpifaunaGardenAverageStructure)
data.scoresnQuadraGardenEpifauna$Comparison = QuadraEpifaunaGardenAverageStructure$Comparison
data.scoresnQuadraGardenEpifauna$Site =       QuadraEpifaunaGardenAverageStructure$Site
data.scoresnQuadraGardenEpifauna$SiteType2 =  QuadraEpifaunaGardenAverageStructure$SiteType2
data.scoresnQuadraGardenEpifauna$SiteType =   QuadraEpifaunaGardenAverageStructure$SiteType
data.scoresnQuadraGardenEpifauna$Area =       QuadraEpifaunaGardenAverageStructure$Area
data.scoresnQuadraGardenEpifauna$Year =       QuadraEpifaunaGardenAverageStructure$Year
data.scoresnQuadraGardenEpifauna$Year=as.factor(data.scoresnQuadraGardenEpifauna$Year)
str(data.scoresnQuadraGardenEpifauna)
head(data.scoresnQuadraGardenEpifauna)
## data.scoresnQuadraGardenEpifauna   Stress 0.13988

names(QuadraEpifaunaAquacultureAverage)
QuadraEpifaunaAquacultureAverageStructure = QuadraEpifaunaAquacultureAverage[,1:8]
QuadraEpifaunaAquacultureAverageData   = QuadraEpifaunaAquacultureAverage[,9:31]
QuadraEpifaunaAquacultureAverageData2=((QuadraEpifaunaAquacultureAverageData)+1)
### nMDS
nMDSQuadraAquacultureEpifaunaMid<- metaMDS(QuadraEpifaunaAquacultureAverageData2, distance = "bray", k = 2) 
nMDSQuadraAquacultureEpifaunaMid$stress  
stressplot(nMDSQuadraAquacultureEpifaunaMid)
data.scoresQuadraAquacultureEpifauna = as.data.frame(scores(nMDSQuadraAquacultureEpifaunaMid))
# Add structure
names(QuadraEpifaunaAquacultureAverageStructure)
data.scoresQuadraAquacultureEpifauna$Comparison = QuadraEpifaunaAquacultureAverageStructure$Comparison
data.scoresQuadraAquacultureEpifauna$Site =       QuadraEpifaunaAquacultureAverageStructure$Site
data.scoresQuadraAquacultureEpifauna$SiteType2 =  QuadraEpifaunaAquacultureAverageStructure$SiteType2
data.scoresQuadraAquacultureEpifauna$SiteType =   QuadraEpifaunaAquacultureAverageStructure$SiteType
data.scoresQuadraAquacultureEpifauna$Area =       QuadraEpifaunaAquacultureAverageStructure$Area
data.scoresQuadraAquacultureEpifauna$Year =       QuadraEpifaunaAquacultureAverageStructure$Year
data.scoresQuadraAquacultureEpifauna$Year=as.factor(data.scoresQuadraAquacultureEpifauna$Year)
str(data.scoresQuadraAquacultureEpifauna)
head(data.scoresQuadraAquacultureEpifauna)
## data.scoresQuadraAquacultureEpifauna   Stress 0.096336

#### Plotting nMDS ####
# Epifauna nMDS
## data.scoresCalvertEpifauna   Stress 0.0995
## data.scoresBaynesEpifauna   Stress 0.080716
## data.scoresnQuadraGardenEpifauna   Stress 0.13988
## data.scoresQuadraAquacultureEpifauna   Stress 0.096336

## Pipeline
GardenColors=c("#7CAE00", "#fa7f2d")
FarmColors=c("#ff0000", "#0d7dd9")

### hjust higher is left
### vjust higher is down

# Calvert Epifauna
nMDSscores=data.scoresCalvertEpifauna
CustomColors=GardenColors     #    FarmColors  or  GardenColors

Plot=ggplot(nMDSscores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,stroke = stroke, aes( shape = Year, colour = SiteType2))+ 
  theme(axis.text.y = element_text(colour = "black", size = AxisSize, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = AxisSize), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = LegendSize), 
        axis.title.x = element_text(face = "bold", size = LegendSize, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "", y = "NMDS2", shape = "")  + scale_colour_manual(values = CustomColors)+scale_shape_manual(values=nMDSshapes)

Plot
Plot+geom_text(aes(label=Site),hjust=0, vjust=0) ## Check Sites

### Adding connections, only relavant will be selected
grp.a <- nMDSscores[nMDSscores$SiteType2 == "Clam Garden", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                            "Clam Garden", c("NMDS1", "NMDS2")]), ]
grp.b <- nMDSscores[nMDSscores$SiteType2 == "Farm Reference", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                               "Farm Reference", c("NMDS1", "NMDS2")]), ]
grp.c <- nMDSscores[nMDSscores$SiteType2 == "Garden Reference", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                                 "Garden Reference", c("NMDS1", "NMDS2")]), ]
grp.d <- nMDSscores[nMDSscores$SiteType2 == "Shellfish Farm", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                               "Shellfish Farm", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b,grp.c,grp.d)
hull.data

CalvertEpifaunaMDSPlot=Plot+geom_polygon (data=hull.data,aes (x=NMDS1,y=NMDS2,fill=SiteType2,group=SiteType2),alpha=0.30)+
  scale_fill_manual(values = CustomColors)+guides(fill = F,colour=F,shape=F)+
  annotate("text", y= min(nMDSscores$NMDS2), x =max(nMDSscores$NMDS1),label="Stress = 0.10",hjust=hjust,size = size)
CalvertEpifaunaMDSPlot

## BaynesEpifauna   
nMDSscores=data.scoresBaynesEpifauna
CustomColors=FarmColors   #    FarmColors  or  GardenColors

Plot=ggplot(nMDSscores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Year,stroke = stroke, colour = SiteType2))+ 
  theme(axis.text.y = element_text(colour = "black", size = AxisSize, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = AxisSize), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = LegendSize), 
        axis.title.x = element_text(face = "bold", size = LegendSize, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "", y = "NMDS2", shape = "")  + scale_colour_manual(values = CustomColors)+scale_shape_manual(values=nMDSshapes)
Plot
Plot+geom_text(aes(label=Site),hjust=0, vjust=0) ## Check Sites

### Adding connections, only relavant will be selected
grp.a <- nMDSscores[nMDSscores$SiteType2 == "Clam Garden", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                "Clam Garden", c("NMDS1", "NMDS2")]), ]
grp.b <- nMDSscores[nMDSscores$SiteType2 == "Farm Reference", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                   "Farm Reference", c("NMDS1", "NMDS2")]), ]
grp.c <- nMDSscores[nMDSscores$SiteType2 == "Garden Reference", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                     "Garden Reference", c("NMDS1", "NMDS2")]), ]
grp.d <- nMDSscores[nMDSscores$SiteType2 == "Shellfish Farm", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                   "Shellfish Farm", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b,grp.c,grp.d)
hull.data

BaynesEpifaunaMDSPlot=Plot+geom_polygon (data=hull.data,aes (x=NMDS1,y=NMDS2,fill=SiteType2,group=SiteType2),alpha=0.30)+
  scale_fill_manual(values = CustomColors)+guides(fill = FALSE,colour=F,shape=F) +
  annotate("text", y= min(nMDSscores$NMDS2), x =max(nMDSscores$NMDS1),label="Stress = 0.08",hjust=hjust, size = size)
BaynesEpifaunaMDSPlot

## QuadraGardenEpifauna  
nMDSscores=data.scoresnQuadraGardenEpifauna
CustomColors=GardenColors         #    FarmColors  or  GardenColors

Plot=ggplot(nMDSscores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Year,stroke = stroke, colour = SiteType2))+ 
  theme(axis.text.y = element_text(colour = "black", size = AxisSize, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = AxisSize), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = LegendSize), 
        axis.title.x = element_text(face = "bold", size = LegendSize, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "", y = "NMDS2", shape = "")  + scale_colour_manual(values = CustomColors)+scale_shape_manual(values=nMDSshapes)
Plot
Plot+geom_text(aes(label=Site),hjust=0, vjust=0) ## Check Sites

### Adding connections, only relavant will be selected
grp.a <- nMDSscores[nMDSscores$SiteType2 == "Clam Garden", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                "Clam Garden", c("NMDS1", "NMDS2")]), ]
grp.b <- nMDSscores[nMDSscores$SiteType2 == "Farm Reference", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                   "Farm Reference", c("NMDS1", "NMDS2")]), ]
grp.c <- nMDSscores[nMDSscores$SiteType2 == "Garden Reference", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                     "Garden Reference", c("NMDS1", "NMDS2")]), ]
grp.d <- nMDSscores[nMDSscores$SiteType2 == "Shellfish Farm", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                   "Shellfish Farm", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b,grp.c,grp.d)
hull.data

QuadraGardenEpifaunaMDSPlot=Plot+geom_polygon (data=hull.data,aes (x=NMDS1,y=NMDS2,fill=SiteType2,group=SiteType2),alpha=0.30)+
  scale_fill_manual(values = CustomColors)+guides(fill = F,colour=F,shape=F) +
  annotate("text", y= min(nMDSscores$NMDS2), x =max(nMDSscores$NMDS1),label="Stress = 0.14",hjust=hjust,size = size)
QuadraGardenEpifaunaMDSPlot

### hjust higher is left
### vjust higher is down

# setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Thesis Chapter, Nov 2020/Analysis/Results /nMDS/Regional Figure Construction")
# jpeg(filename = "QuadraGardenEpifaunaMDSPlot.jpeg", width = 25, height = 15, units = "cm", pointsize = 15, quality = 100, res = 600)  
# QuadraGardenEpifaunaMDSPlot
# dev.off()
# setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")

## QuadraAquacultureEpifauna   
nMDSscores=data.scoresQuadraAquacultureEpifauna
CustomColors=FarmColors         #    FarmColors  or  GardenColors

Plot=ggplot(nMDSscores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Year,stroke = stroke, colour = SiteType2))+ 
  theme(axis.text.y = element_text(colour = "black", size = AxisSize, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = AxisSize), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = LegendSize), 
        axis.title.x = element_text(face = "bold", size = LegendSize, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "", y = "NMDS2", shape = "")  + scale_colour_manual(values = CustomColors)+scale_shape_manual(values=nMDSshapes)
Plot
Plot+geom_text(aes(label=Site),hjust=0, vjust=0) ## Check Sites

### Adding connections, only relavant will be selected
grp.a <- nMDSscores[nMDSscores$SiteType2 == "Clam Garden", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                "Clam Garden", c("NMDS1", "NMDS2")]), ]
grp.b <- nMDSscores[nMDSscores$SiteType2 == "Farm Reference", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                   "Farm Reference", c("NMDS1", "NMDS2")]), ]
grp.c <- nMDSscores[nMDSscores$SiteType2 == "Garden Reference", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                     "Garden Reference", c("NMDS1", "NMDS2")]), ]
grp.d <- nMDSscores[nMDSscores$SiteType2 == "Shellfish Farm", ][chull(nMDSscores[nMDSscores$SiteType2 == 
                                                                                   "Shellfish Farm", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b,grp.c,grp.d)
hull.data

QuadraAquacultureEpifaunaMDSPlot=Plot+geom_polygon (data=hull.data,aes (x=NMDS1,y=NMDS2,fill=SiteType2,group=SiteType2),alpha=0.30)+
  scale_fill_manual(values = CustomColors)+guides(fill = F,colour=F,shape=F) +
  annotate("text", y= min(nMDSscores$NMDS2), x =max(nMDSscores$NMDS1),label="Stress = 0.10",hjust=hjust,size = size)
QuadraAquacultureEpifaunaMDSPlot

### hjust higher is left
### vjust higher is down

# setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Thesis Chapter, Nov 2020/Analysis/Results /nMDS/Regional Figure Construction")
# jpeg(filename = "QuadraAquacultureEpifaunaMDSPlotMDSPlot.jpeg", width = 25, height = 15, units = "cm", pointsize = 15, quality = 100, res = 600)  
# QuadraAquacultureEpifaunaMDSPlot
# dev.off()
# setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")

margins=unit(c(0.6, 0.6, 0.6, 0.6), "cm")

CalvertEpifaunaMDSPlot2=CalvertEpifaunaMDSPlot+guides(shape = FALSE,colour=FALSE)+theme(plot.margin = margins)
BaynesEpifaunaMDSPlot2=BaynesEpifaunaMDSPlot+guides(shape = FALSE,colour=FALSE)+theme(plot.margin = margins)
QuadraGardenEpifaunaMDSPlot2=QuadraGardenEpifaunaMDSPlot+guides(shape = FALSE,colour=FALSE)+theme(plot.margin = margins)
QuadraAquacultureEpifaunaMDSPlot2=QuadraAquacultureEpifaunaMDSPlot+guides(shape = FALSE,colour=FALSE)+theme(plot.margin = margins)

plot_grid(CalvertEpifaunaMDSPlot2,QuadraGardenEpifaunaMDSPlot2,QuadraAquacultureEpifaunaMDSPlot2,BaynesEpifaunaMDSPlot2, 
          labels=c('A)','B)','C)','D)','E)'),nrow = 1,ncol=4)

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Thesis Chapter, Nov 2020/Analysis/Results /nMDS/Regional Figure Construction")
jpeg(filename = "nMDSCombinedRegionYear4.jpeg", width = 70, height = 32, units = "cm", pointsize = 15, quality = 100, res = 600)  

plot_grid(CalvertEpifaunaMDSPlot2,QuadraGardenEpifaunaMDSPlot2,QuadraAquacultureEpifaunaMDSPlot2,BaynesEpifaunaMDSPlot2, 
          CalvertBivalveMDSPlot2,QuadraBivalveGardenMDSPlot2,QuadraBivalveAquacultureMDSPlot2,BaynesBivalveMDSPlot2,
          labels=c('A)','B)','C)','D)','E)','F)','G)','H)'),label_size=19,nrow = 2,ncol=4,align="hv")

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


#### Figure 3: nMDS plotting it together #### 

EpiFaunanMDSPlotFinal
EpiFaunanMDSPlotFinalnoLegend

EpiFaunanMDSPlotFinalnoLegend

CalvertEpifaunaMDSPlot
BaynesEpifaunaMDSPlot
QuadraGardenEpifaunaMDSPlot
QuadraAquacultureEpifaunaMDSPlot

Image7 <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Thesis Chapter, Nov 2020/Analysis/Results /nMDS/Regional Figure Construction/Legend7.tiff")
Image7
Legend7 <- image_ggplot(Image7)


layout(matrix(c(1, 1,
                2, 3,
                4, 5,
                6,6), nrow=4, byrow=TRUE),widths=c(1,1,1,1), heights=c(1.50,1,1,0.2))
par(cex = 1.045)
layout.show(n=11)

par(mai = c(0.98, 1, 0.75, 0.5))
EpiFaunanMDSPlotFinalnoLegend

par(mai = c(0.98, 1, 0.9, 0.5))
CalvertEpifaunaMDSPlot
QuadraGardenEpifaunaMDSPlot
QuadraAquacultureEpifaunaMDSPlot
BaynesEpifaunaMDSPlot

par(mai = c(1, 1, 0.9, 0.5))
CalvertBivalveMDSPlot
QuadraBivalveGardenMDSPlot
QuadraBivalveAquacultureMDSPlot
BaynesBivalveMDSPlot

library(gridExtra)
library(egg)
library(grid)
library(gtable)
library(gtable)

margins1 = unit(c(0.55,0.4 ,0.4  ,0.4), "cm")
margins2 = unit(c(0   ,0.4 ,0.15 ,0.4), "cm")
LabelSize=25

###
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Thesis Chapter, Nov 2020/Analysis/Results /nMDS/Regional Figure Construction")
jpeg(filename = "nMDSCombinedLegendFinal3.jpeg", width = 34, height = 55, units = "cm", pointsize = 15, quality = 100, res = 600)  

grid.arrange(EpiFaunanMDSPlotFinalnoLegend+theme(plot.margin =margins1)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black")) +labs(subtitle = "A)"),
             
             CalvertEpifaunaMDSPlot+theme(plot.margin =margins2)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black")) +labs(subtitle = "B)"),
             QuadraGardenEpifaunaMDSPlot+theme(plot.margin =margins2)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black")) +labs(subtitle = "C)"),
             QuadraAquacultureEpifaunaMDSPlot+theme(plot.margin =margins2)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black")) +labs(subtitle = "D)"),
             BaynesEpifaunaMDSPlot+theme(plot.margin =margins2)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black")) +labs( subtitle = "E)"),
             Legend7,
             heights=c(1.75,1,1,0.2),      
             widths = c(1,1),
             layout_matrix = rbind(c(1, 1),
                                   c(2, 3),
                                   c(4, 5),
                                   c(6, 6))
)
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")




##### Supplemental Benthic Composition nMDS ####

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")
EpifaunaDataRegression=read.csv("SFandCGBiodiversityDataNov2020 RegressionTree Family V8.csv", header=T)
EpifaunaDataRegression

### Few concerns- 
# lots of zero's
# averaging at site level
# mid intertidal vs low,mid,high

### nMDS 1: Mid data not averaged
BenthicDataMid <- subset(EpifaunaDataRegression, EpifaunaDataRegression$Zone == "Mid")
names(BenthicDataMid)

BenthicDataMid2 <- BenthicDataMid[,10:23]
names(BenthicDataMid2)
BenthicDataStructure=BenthicDataMid[,1:8]
names(BenthicDataStructure)
BenthicDataMid3=cbind(BenthicDataMid2,BenthicDataStructure)
str(BenthicDataMid3)
BenthicDataMid3

### Removing NAs
BenthicDataMid4=(BenthicDataMid3[complete.cases(BenthicDataMid3), ])
str(BenthicDataMid4)
names(BenthicDataMid4)
## Removing zero's
BenthicDataMid5=BenthicDataMid4[, colSums(BenthicDataMid4 != 0) > 0]
str(BenthicDataMid5)
names(BenthicDataMid5)

### Average 
str(BenthicDataMid5)
names(BenthicDataMid5)
BenthicDataMid5$Year=as.factor(BenthicDataMid5$Year)

BenthicDataMid6=aggregate(.~Area+Comparison+SiteType+SiteType2+Year+Site, BenthicDataMid5, mean)

BenthicDataMid7=aggregate(.~Area+Comparison+SiteType+SiteType2+Site, BenthicDataMid6, mean)

str(BenthicDataMid7)
names(BenthicDataMid7)
BenthicDataMid8   = BenthicDataMid7[,7:20]
BenthicDataMidStructure2 = BenthicDataMid7[,1:5]

###       Optional Dummy Variable
# BenthicData6=((BenthicData5+1))

### nMDS
nMDSEpiBenthicComposition<- metaMDS(BenthicDataMid8, distance = "bray", k = 3) 
nMDSEpiBenthicComposition$stress  
stressplot(nMDSEpiBenthicComposition)

#extract NMDS scores (x and y coordinates)
data.scoresEpi = as.data.frame(scores(nMDSEpiBenthicComposition))
# Add data structure
names(BenthicDataMidStructure2)
data.scoresEpi$Comparison = BenthicDataMidStructure2$Comparison
data.scoresEpi$Site = BenthicDataMidStructure2$Site
data.scoresEpi$SiteType2 = BenthicDataMidStructure2$SiteType2
data.scoresEpi$SiteType = BenthicDataMidStructure2$SiteType
data.scoresEpi$Area = BenthicDataMidStructure2$Area
str(data.scoresEpi)
head(data.scoresEpi)

## Plot
GardenColor=("#7CAE00")   #### Green 
FarmColor=("#0d7dd9")     #### Blue 
ReferenceFarmColor=("#ff0000") #### Red
ReferenceGardenColor=("#fa7f2d") #### Orange

hjust=0.85
vjust=-0.1
size=6
stroke = 1.4
nMDSshapes=c(5,12)
AxisSize=18
LegendSize=18    
textsize=6

BenthicComMDSPlot=ggplot(data.scoresEpi, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Area, colour = SiteType2))+ 
  theme(axis.text.y = element_text(colour = "black", size = AxisSize, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = AxisSize), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = LegendSize), 
        axis.title.x = element_text(face = "bold", size = LegendSize, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "", y = "NMDS2", shape = "")  + 
  scale_colour_manual(values = c("#7CAE00", "#ff0000","#fa7f2d","#0d7dd9"))
BenthicComMDSPlot



### Adding ordiellipse or connections, 2 options
### Option 1
grp.a <- data.scoresEpi[data.scoresEpi$SiteType2 == "Clam Garden", ][chull(data.scoresEpi[data.scoresEpi$SiteType2 == 
                                                                                            "Clam Garden", c("NMDS1", "NMDS2")]), ]

grp.b <- data.scoresEpi[data.scoresEpi$SiteType2 == "Farm Reference", ][chull(data.scoresEpi[data.scoresEpi$SiteType2 == 
                                                                                               "Farm Reference", c("NMDS1", "NMDS2")]), ]
grp.c <- data.scoresEpi[data.scoresEpi$SiteType2 == "Garden Reference", ][chull(data.scoresEpi[data.scoresEpi$SiteType2 == 
                                                                                                 "Garden Reference", c("NMDS1", "NMDS2")]), ]
grp.d <- data.scoresEpi[data.scoresEpi$SiteType2 == "Shellfish Farm", ][chull(data.scoresEpi[data.scoresEpi$SiteType2 == 
                                                                                               "Shellfish Farm", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b,grp.c,grp.d)  #combine
hull.data

EpiFaunanMDSPlotFinal=BenthicComMDSPlot+geom_polygon (data=hull.data,aes (x=NMDS1,y=NMDS2,fill=SiteType2,group=SiteType2),alpha=0.30)+
  scale_fill_manual(values = c("#7CAE00", "#ff0000","#fa7f2d","#0d7dd9"))+guides(fill = FALSE)+
  geom_text(x=0.215, y=0.159, label="Stress = 0.12",size=textsize)
EpiFaunanMDSPlotFinal


## Adding a new column for plotting, region and type
# write.csv(data.scoresEpi,"data.scoresEpi.csv")
data.scoresEpi2=read.csv("data.scoresEpi2.csv",header=T)
names(data.scoresEpi2)
levels(data.scoresEpi2$SiteType3)


BenthicComMDSPlot2=ggplot(data.scoresEpi2, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Area, colour = SiteType2))+ 
  theme(axis.text.y = element_text(colour = "black", size = AxisSize, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = AxisSize), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = LegendSize), 
        axis.title.x = element_text(face = "bold", size = LegendSize, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "", y = "NMDS2", shape = "")  + 
  scale_colour_manual(values = c("#7CAE00", "#ff0000","#fa7f2d","#0d7dd9"))
BenthicComMDSPlot2


### Adding ordiellipse or connections, 2 options
### 
grp.a <- data.scoresEpi[data.scoresEpi$SiteType2 == "Clam Garden", ][chull(data.scoresEpi[data.scoresEpi$SiteType2 == 
                                                                                            "Clam Garden", c("NMDS1", "NMDS2")]), ]



grp.a <- data.scoresEpi2[data.scoresEpi2$SiteType3 == "Clam Garden, Calvert", ][chull(data.scoresEpi2[data.scoresEpi2$SiteType3 == 
                                                                                                        "Clam Garden, Calvert", c("NMDS1", "NMDS2")]), ]
grp.b <- data.scoresEpi2[data.scoresEpi2$SiteType3 == "Clam Garden, Quadra", ][chull(data.scoresEpi2[data.scoresEpi2$SiteType3 == 
                                                                                                       "Clam Garden, Quadra", c("NMDS1", "NMDS2")]), ]
grp.c <- data.scoresEpi2[data.scoresEpi2$SiteType3 == "Farm Reference, Baynes", ][chull(data.scoresEpi2[data.scoresEpi2$SiteType3 == 
                                                                                                          "Farm Reference, Baynes", c("NMDS1", "NMDS2")]), ]
grp.d <- data.scoresEpi2[data.scoresEpi2$SiteType3 == "Farm Reference, Quadra", ][chull(data.scoresEpi2[data.scoresEpi2$SiteType3 == 
                                                                                                          "Farm Reference, Quadra", c("NMDS1", "NMDS2")]), ]

grp.e <- data.scoresEpi2[data.scoresEpi2$SiteType3 == "Garden Reference, Calvert", ][chull(data.scoresEpi2[data.scoresEpi2$SiteType3 == 
                                                                                                             "Garden Reference, Calvert", c("NMDS1", "NMDS2")]), ]

grp.f <- data.scoresEpi2[data.scoresEpi2$SiteType3 == "Garden Reference, Quadra", ][chull(data.scoresEpi2[data.scoresEpi2$SiteType3 == 
                                                                                                            "Garden Reference, Quadra", c("NMDS1", "NMDS2")]), ]

grp.g <- data.scoresEpi2[data.scoresEpi2$SiteType3 == "Shellfish Farm, Baynes", ][chull(data.scoresEpi2[data.scoresEpi2$SiteType3 == 
                                                                                                          "Shellfish Farm, Baynes", c("NMDS1", "NMDS2")]), ]

grp.h <- data.scoresEpi2[data.scoresEpi2$SiteType3 == "Shellfish Farm, Quadra", ][chull(data.scoresEpi2[data.scoresEpi2$SiteType3 == 
                                                                                                          "Shellfish Farm, Quadra", c("NMDS1", "NMDS2")]), ]

hull.data <- rbind(grp.a, grp.b,grp.c,grp.d,grp.e,grp.f,grp.g,grp.h)  #combine
hull.data

levels(hull.data$SiteType3)
Col2= c("#7CAE00", "#7CAE00","#ff1100","#ff1100","#fa6920","#fa6920","#0d7dd9","#0d7dd9")

## Plotting Site type by region ##
BenthicComMDSPlotFinal=BenthicComMDSPlot2+geom_polygon (data=hull.data,aes (x=NMDS1,y=NMDS2,fill=SiteType3,group=SiteType3),alpha=0.30)+
  scale_fill_manual(values = Col2)+guides(fill = FALSE)+
  geom_text(x=0.625, y=0.815, label="Stress = 0.11",size=textsize)
BenthicComMDSPlotFinal



## Export ##
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Final Figures V1/Final Figures V2")

jpeg(filename = "nMDSBenthicComposition.jpeg",width = 25, height = 20, units = "cm", pointsize = 9, quality = 100,res=800)

BenthicComMDSPlotFinal

dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")















## End nMDS ##
#### 

####  Supplemental Permanova Analysis ####
SFandCGDataMid2=read.csv("SFandCGDataMid2.csv", header=T)
SFandCGDataMid2
SFandCGDataMid2 %>% count(SiteType,SiteType2,Site, Year)
names(SFandCGDataMid2)

names(SFandCGDataMid2)
MidEpifaunaData2 <- SFandCGDataMid2[,34:93]
names(MidEpifaunaData2)
MidDataStructure2=SFandCGDataMid2[,1:8]
names(MidDataStructure2)
EpifaunaSpeciesData3=cbind(MidDataStructure2,MidEpifaunaData2)

CalvertEpifauna <- subset(EpifaunaSpeciesData5, EpifaunaSpeciesData5$Area == "Calvert")
QuadraEpifauna <- subset(EpifaunaSpeciesData5, EpifaunaSpeciesData5$Area == "Quadra")
BaynesEpifauna <- subset(EpifaunaSpeciesData5, EpifaunaSpeciesData5$Area == "Baynes")
QuadraEpifaunaGarden <- subset(QuadraEpifauna, QuadraEpifauna$Comparison == "Garden")
QuadraEpifaunaAquaculture <- subset(QuadraEpifauna, QuadraEpifauna$Comparison == "Aquaculture")

## Summary data
EpifaunaSpeciesData3
EpifaunaSpeciesData4=(EpifaunaSpeciesData3[complete.cases(EpifaunaSpeciesData3), ]) ### Removing NAs
EpifaunaSpeciesData5=EpifaunaSpeciesData4[, colSums(EpifaunaSpeciesData4 != 0) > 0] ## Removing zero's
names(EpifaunaSpeciesData5)
EpifaunaSpeciesDataStructure = EpifaunaSpeciesData5[,1:8]
EpifaunaSpeciesData8   = EpifaunaSpeciesData5[,9:59]
EpifaunaSpeciesData9=((EpifaunaSpeciesData8+1))


## Summary Data Aquaculture
EpifaunaSpeciesData3
EpifaunaSpeciesAquaculture <- subset(EpifaunaSpeciesData3, EpifaunaSpeciesData3$Comparison == "Aquaculture")
EpifaunaSpeciesAquaculture2=(EpifaunaSpeciesAquaculture[complete.cases(EpifaunaSpeciesAquaculture), ]) ### Removing NAs
EpifaunaSpeciesAquaculture3=EpifaunaSpeciesAquaculture2[, colSums(EpifaunaSpeciesAquaculture2 != 0) > 0] ## Removing zero's
names(EpifaunaSpeciesAquaculture3)
EpifaunaSpeciesAquacultureDataStructure = EpifaunaSpeciesAquaculture3[,1:8]
EpifaunaSpeciesAquacultureData   = EpifaunaSpeciesAquaculture3[,9:41]
EpifaunaSpeciesAquacultureData2=((EpifaunaSpeciesAquacultureData+1))


## Summary Data Garden
EpifaunaSpeciesData3
EpifaunaSpeciesGarden <- subset(EpifaunaSpeciesData3, EpifaunaSpeciesData3$Comparison == "Garden")
EpifaunaSpeciesGarden2=(EpifaunaSpeciesGarden[complete.cases(EpifaunaSpeciesGarden), ]) ### Removing NAs
EpifaunaSpeciesGarden3=EpifaunaSpeciesGarden2[, colSums(EpifaunaSpeciesGarden2 != 0) > 0] ## Removing zero's
names(EpifaunaSpeciesGarden3)
EpifaunaSpeciesGardenDataStructure = EpifaunaSpeciesGarden3[,1:8]
EpifaunaSpeciesGardenData   = EpifaunaSpeciesGarden3[,9:49]
EpifaunaSpeciesGardenData2=((EpifaunaSpeciesGardenData+1))


## CalvertEpifauna
str(CalvertEpifauna)  
CalvertEpifauna2=(CalvertEpifauna[complete.cases(CalvertEpifauna), ]) ### Removing NAs
CalvertEpifauna3=CalvertEpifauna2[, colSums(CalvertEpifauna2 != 0) > 0] ## Removing zero's
names(CalvertEpifauna3)
str(CalvertEpifauna3)  
CalvertEpifaunaDataStructure = CalvertEpifauna3[,1:8]
CalvertEpifaunaData8   = CalvertEpifauna3[,9:41]
CalvertEpifaunaData9=((CalvertEpifaunaData8+1))


## BaynesEpifauna
str(BaynesEpifauna)
BaynesEpifauna2=(BaynesEpifauna[complete.cases(BaynesEpifauna), ]) ### Removing NAs
BaynesEpifauna3=BaynesEpifauna2[, colSums(BaynesEpifauna2 != 0) > 0] ## Removing zero's
names(BaynesEpifauna3)
BaynesEpifaunaDataStructure = BaynesEpifauna3[,1:8]
BaynesEpifaunaData8   = BaynesEpifauna3[,9:35]
BaynesEpifaunaData9=((BaynesEpifaunaData8+1))


## QuadraEpifaunaGarden
str(QuadraEpifaunaGarden)
QuadraEpifaunaGarden2=(QuadraEpifaunaGarden[complete.cases(QuadraEpifaunaGarden), ]) ### Removing NAs
QuadraEpifaunaGarden3=QuadraEpifaunaGarden2[, colSums(QuadraEpifaunaGarden2 != 0) > 0] ## Removing zero's
names(QuadraEpifaunaGarden3)
QuadraEpifaunaGardenDataStructure = QuadraEpifaunaGarden3[,1:8]
QuadraEpifaunaGardenData8   = QuadraEpifaunaGarden3[,9:32]
QuadraEpifaunaGardenData9=((QuadraEpifaunaGardenData8+1))


## QuadraEpifaunaAquaculture
str(QuadraEpifaunaAquaculture)
QuadraEpifaunaAquaculture2=(QuadraEpifaunaAquaculture[complete.cases(QuadraEpifaunaAquaculture), ]) ### Removing NAs
QuadraEpifaunaAquaculture3=QuadraEpifaunaAquaculture2[, colSums(QuadraEpifaunaAquaculture2 != 0) > 0] ## Removing zero's
names(QuadraEpifaunaAquaculture3)
QuadraEpifaunaAquacultureDataStructure = QuadraEpifaunaAquaculture3[,1:8]
QuadraEpifaunaAquacultureData8   = QuadraEpifaunaAquaculture3[,9:32]
QuadraEpifaunaAquacultureData9=((QuadraEpifaunaAquacultureData8+1))

### Data
#  EpifaunaSpeciesData9                 EpifaunaSpeciesDataStructure
#  CalvertEpifaunaData9                 CalvertEpifaunaDataStructure
#  QuadraEpifaunaGardenData9            QuadraEpifaunaGardenDataStructure
#  QuadraEpifaunaAquacultureData9       QuadraEpifaunaAquacultureDataStructure
#  BaynesEpifaunaData9                  BaynesEpifaunaDataStructure      

##   EpifaunaSpeciesAquacultureData2        EpifaunaSpeciesAquacultureDataStructure
##   EpifaunaSpeciesGardenData2             EpifaunaSpeciesGardenDataStructure


EpifaunaSpeciesDataStructure$Year=as.factor(EpifaunaSpeciesDataStructure$Year)
CalvertEpifaunaDataStructure$Year=as.factor(CalvertEpifaunaDataStructure$Year)
QuadraEpifaunaGardenDataStructure$Year=as.factor(QuadraEpifaunaGardenDataStructure$Year)
QuadraEpifaunaAquacultureDataStructure$Year=as.factor(QuadraEpifaunaAquacultureDataStructure$Year)
BaynesEpifaunaDataStructure$Year=as.factor(BaynesEpifaunaDataStructure$Year)

## PERMANOVA- Adonis
## Summary D4
Summary1=adonis(EpifaunaSpeciesData9 ~ SiteType2/Site+Year,strata=EpifaunaSpeciesDataStructure$Area, data = EpifaunaSpeciesDataStructure,method = "bray")
Summary1
print(Summary1)

## Calvert Gardens D3        Ideal
CalvertGarden1=adonis(CalvertEpifaunaData9 ~ SiteType2/Site/Quad+Year,data = CalvertEpifaunaDataStructure,method = "bray")
CalvertGarden1
print(CalvertGarden1)

## Quadra Gardens D1    Ideal
QuadraGarden1=adonis(QuadraEpifaunaGardenData9 ~ SiteType2/Site/Quad+Year, data = QuadraEpifaunaGardenDataStructure,method = "bray")
QuadraGarden1
print(QuadraGarden1)

## Quadra Aquaculture   Ideal
QuadraAquaculture1=adonis(QuadraEpifaunaAquacultureData9 ~ SiteType2/Site/Quad+Year,  data = QuadraEpifaunaAquacultureDataStructure,method = "bray")
QuadraAquaculture1
print(QuadraAquaculture1)

## Baynes Aquaculture
BaynesAquaculture2=adonis(BaynesEpifaunaData9 ~ SiteType2/Site/Quad+Year, data = BaynesEpifaunaDataStructure,method = "bray")
BaynesAquaculture2
print(BaynesAquaculture2)
### Make into table. 


#### Figure 3:  SIMPER plots, tables in supplemenal ####
#### Epifauna SIMPER #### 
## Data set up ##
EpifaunaSpeciesData5
str(EpifaunaSpeciesData5)
names(EpifaunaSpeciesData5)

## Merging Barnacle Spp and Balanus.glandula to be BarnacleMerged
EpifaunaSpeciesData5$Balanus.glandula
EpifaunaSpeciesData5$Barnacle.spp.
EpifaunaSpeciesData5$BarnacleMerged =(EpifaunaSpeciesData5$Barnacle.spp.) + (EpifaunaSpeciesData5$Balanus.glandula)
EpifaunaSpeciesData5$BarnacleMerged
EpifaunaSpeciesData6 = subset(EpifaunaSpeciesData5, select = -c(Balanus.glandula,Barnacle.spp.) )
names(EpifaunaSpeciesData6)

#Naming species for plot
#Epifauna
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Corophiidae.Amphipod."] = "Corophiidae"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Anenome.SubstrateCounts"] = "Actiniaria"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="batillaria.attramentaria"] = "B. attramentaria"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Neostylidium.eschrichtii"] = "N. eschrichtii"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Magallana.gigas"] = "M. gigas"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Lottia.scutum"] = "L. scutum"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Lottia.spp."] = "Lottia spp."
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Lottia.pelta"] = "L. pelta"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Lottia.digitalis"] = "L. digitalis"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Lottia.Persona"] = "L. Persona"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Idotea.spp."] = "Idotea spp."
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Littorina.scutulata"] = "L. scutulata"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Littorina.spp."] = "Littorina spp."
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Littorina.sitkana"] = "L. sitkana"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Chthamalus.dalli"] = "C. dalli"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Balanus.Nubilus"] = "B. nubilus"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Semibalanus.cariosus"] = "S. cariosus"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Balanus.crenatus"] = "B. crenatus"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Mytilus.spp."] = "Mytilus spp."
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Pagurus.hirsutiusculus"] = "P. hirsutiusculus"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Pagurus.spp."] = "Pagurus spp."
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Hermissenda.spp."] = "Hermissenda spp."
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Dendraster.excentricus"] = "D. excentricus"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Nucella.spp."] = "Nucella spp."
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Strongylocentrotus.spp."] = "Strongylocentrotus spp."
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Upogebia.pugettensis"] = "U. pugettensis"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Gnorimosphaeroma.oregonensis"] = "G. oregonensis"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Lirabuccinum.dirum"] = "L. dirum"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Saxidomus.gigantea"] = "S. gigantea"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Nutricola.tantilla"] = "N. tantilla"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Leukoma.staminea"] = "L. staminea"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Bivalvia"] = "Bivalvia"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Ruditapes.philippinarum"] = "R. philippinarum"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Clinocardium.nuttallii"] = "C. nuttallii"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Leptasterias.hexactis"] = "L. hexactis"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Pisaster.ochraceus"] = "P. ochraceus"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Tonicella.lineata"] = "T. lineata"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Bryozoa"] = "Bryozoa"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Haminoea.vesicula"] = "H. vesicula"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Polychaeta"] = "Polychaeta"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Nemertea"] = "Nemertea"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Balanophyllia.elegans"] = "B. elegans"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Caridean.Shrimps"] = "Caridea"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="BarnacleMerged"] = "Barnacle Spp."
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Hemigrapsus.nudus"] = "H. nudus"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Hemigrapsus.oregonensis"] = "H. oregonensis"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Cancer.productus"] = "C. productus"
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Pugettia.spp."] = "Pugettia spp."
colnames(EpifaunaSpeciesData6)[colnames(EpifaunaSpeciesData6)=="Hemigrapsus.spp."] = "Hemigrapsus spp."

EpifaunaSpeciesData6
names(EpifaunaSpeciesData6)

## Split into comparisions
EpifaunaSpeciesData6$Comparison
GardenEpifauna <- subset(EpifaunaSpeciesData6, EpifaunaSpeciesData6$Comparison == "Garden")
AquacultureEpifauna <- subset(EpifaunaSpeciesData6, EpifaunaSpeciesData6$Comparison == "Aquaculture")

## Removing zeros
GardenEpifauna2=GardenEpifauna[, colSums(GardenEpifauna != 0) > 0]
AquacultureEpifauna2=AquacultureEpifauna[, colSums(AquacultureEpifauna != 0) > 0]

## data sets
GardenEpifauna2
AquacultureEpifauna2

#### Figure 3: Epifauna Summary Gardens ####
### Splitting the data
names(GardenEpifauna2)
GardenEpifaunaStructure = GardenEpifauna2[,1:8]
GardenEpifaunaData   = GardenEpifauna2[,9:48]
GardenEpifaunaData2=(GardenEpifaunaData+1)

### SIMPER
SpeciesSIMPER<-(simper(GardenEpifaunaData2, GardenEpifaunaStructure$SiteType, permutations=999))
summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$`Clam Garden_Reference`[,1]), na.rm = TRUE)
SIMPERscore=sum(as.numeric(SpeciesSIMPERSum$`Clam Garden_Reference`[,1]), na.rm = TRUE)

### communities are 67.7% dissimilar 
summary(SpeciesSIMPER) 
SIMPERsummary=SpeciesSIMPERSum$`Clam Garden_Reference`

SIMPERsummary$averagePercent =(SIMPERsummary$average)*100
SIMPERsummary$Directionality =(SIMPERsummary$ava) - (SIMPERsummary$avb)
  
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="ava"] = "GardenDensity"
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="avb"] = "ReferenceDensity"

SIMPERsummary$Directionality =(SIMPERsummary$GardenDensity) - (SIMPERsummary$ReferenceDensity)
names(SIMPERsummary)
SIMPERsummary$Directionality2[SIMPERsummary$Directionality > 0 ] <- "1"
SIMPERsummary$Directionality2[SIMPERsummary$Directionality < 0 ] <- "-1"  
SIMPERsummary$Directionality2=as.numeric(SIMPERsummary$Directionality2)
SIMPERsummary$Average.Dissim = (SIMPERsummary$averagePercent)*(SIMPERsummary$Directionality2)

GardenEpifaunaSIMPER=SIMPERsummary
#write.csv(GardenEpifaunaSIMPER,"GardenEpifaunaSIMPER.csv")
GardenEpifaunaSIMPER2=read.csv("GardenEpifaunaSIMPER.csv",header=T)
colnames(GardenEpifaunaSIMPER2)[colnames(GardenEpifaunaSIMPER2)=="X"] = "Species"
names(GardenEpifaunaSIMPER2)

## plot
GardenColor=alpha("#7CAE00")   #### Green 
FarmColor=alpha("#0d7dd9")     #### Blue 
ReferenceFarmColor=alpha("#ff1100") #### Red
ReferenceGardenColor=alpha("#fa6920") #### Orange
col=c(ReferenceGardenColor,GardenColor)
italic.text <- element_text(face = "italic")

Gardencol1=c(ReferenceGardenColor,GardenColor)
Gardencol2=c(GardenColor, ReferenceGardenColor)
Farmcol1=c(ReferenceFarmColor, FarmColor)
Farmcol2=c(FarmColor,ReferenceFarmColor)

GardenEpifaunaSIMPER2$Species <- factor(GardenEpifaunaSIMPER2$Species, levels = GardenEpifaunaSIMPER2$Species[order(GardenEpifaunaSIMPER2$Average.Dissim)])
GardenEpifaunaSIMPER2$Directionality2=as.factor(GardenEpifaunaSIMPER2$Directionality2)
size2=5

GardenEpifaunaSIMPERPlot=ggplot(GardenEpifaunaSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-50,50)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Gardencol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
GardenEpifaunaSIMPERPlot

GardenEpifaunaSIMPERPlot2=GardenEpifaunaSIMPERPlot+geom_vline(xintercept = nlevels(GardenEpifaunaSIMPER2$Species)+0.6,size=0.6)
GardenEpifaunaSIMPERPlot3=GardenEpifaunaSIMPERPlot2+annotate("text", y= -40, x =nlevels(GardenEpifaunaSIMPER2$Species),label="Dissimilarity 67.66%",size = size2)
GardenEpifaunaSIMPERPlot4=GardenEpifaunaSIMPERPlot3 + theme(axis.text.y = italic.text)
GardenEpifaunaSIMPERPlot4


#### Figure 3: Epifauna Summary Aquaculture ####

AquacultureEpifauna <- subset(EpifaunaSpeciesData6, EpifaunaSpeciesData6$Comparison == "Aquaculture")
AquacultureEpifauna2=AquacultureEpifauna[, colSums(AquacultureEpifauna != 0) > 0]

### Splitting the data
names(AquacultureEpifauna2)
AquacultureEpifaunaStructure = AquacultureEpifauna2[,1:8]
AquacultureEpifaunaData   = AquacultureEpifauna2[,9:40]
AquacultureEpifaunaData2=(AquacultureEpifaunaData+1)

### SIMPER
SpeciesSIMPER<-(simper(AquacultureEpifaunaData2, AquacultureEpifaunaStructure$SiteType2, permutations=999))
summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$`Farm Reference_Shellfish Farm`[,1]), na.rm = TRUE)

### communities are 74.86% dissimilar 
summary(SpeciesSIMPER)
SIMPERsummary=SpeciesSIMPERSum$`Farm Reference_Shellfish Farm`

SIMPERsummary$averagePercent =(SIMPERsummary$average)*100
SIMPERsummary
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="ava"] = "ReferenceDensity"
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="avb"] = "FarmDensity"
SIMPERsummary$Directionality =(SIMPERsummary$FarmDensity) - (SIMPERsummary$ReferenceDensity)

SIMPERsummary$Directionality2[SIMPERsummary$Directionality > 0 ] <- "1"
SIMPERsummary$Directionality2[SIMPERsummary$Directionality < 0 ] <- "-1"  
SIMPERsummary$Directionality2=as.numeric(SIMPERsummary$Directionality2)
SIMPERsummary$Average.Dissim = (SIMPERsummary$averagePercent)*(SIMPERsummary$Directionality2)

AquacultureEpifaunaSIMPER=SIMPERsummary
# write.csv(AquacultureEpifaunaSIMPER,"AquacultureEpifaunaSIMPER.csv")
AquacultureEpifaunaSIMPER2=read.csv("AquacultureEpifaunaSIMPER.csv",header=T)
colnames(AquacultureEpifaunaSIMPER2)[colnames(AquacultureEpifaunaSIMPER2)=="X"] = "Species"
names(AquacultureEpifaunaSIMPER2)

## plot
GardenColor 
FarmColor 
ReferenceFarmColor
ReferenceGardenColor
col=c(ReferenceFarmColor,FarmColor)

AquacultureEpifaunaSIMPER2$Species <- factor(AquacultureEpifaunaSIMPER2$Species, levels = AquacultureEpifaunaSIMPER2$Species[order(AquacultureEpifaunaSIMPER2$Average.Dissim)])
AquacultureEpifaunaSIMPER2$Directionality2=as.factor(AquacultureEpifaunaSIMPER2$Directionality2)
size2=5

AquacultureEpifaunaSIMPERPlot=ggplot(AquacultureEpifaunaSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-50,50)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Farmcol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
AquacultureEpifaunaSIMPERPlot

AquacultureEpifaunaSIMPERPlot2=AquacultureEpifaunaSIMPERPlot+geom_vline(xintercept = nlevels(AquacultureEpifaunaSIMPER2$Species)+0.6,size=0.6)
AquacultureEpifaunaSIMPERPlot3=AquacultureEpifaunaSIMPERPlot2+annotate("text", y= -40, x =nlevels(AquacultureEpifaunaSIMPER2$Species),label="Dissimilarity 74.86%",size = size2)
AquacultureEpifaunaSIMPERPlot4=AquacultureEpifaunaSIMPERPlot3 + theme(axis.text.y = italic.text)
AquacultureEpifaunaSIMPERPlot4

#### Epifauna Calvert Gardens SIMPER #### 
EpifaunaSpeciesData6
GardenEpifauna <- subset(EpifaunaSpeciesData6, EpifaunaSpeciesData6$Comparison == "Garden")
GardenEpifaunaCalvert <- subset(GardenEpifauna, GardenEpifauna$Area == "Calvert")
GardenEpifaunaCalvert
GardenEpifaunaCalvert2=GardenEpifaunaCalvert[, colSums(GardenEpifaunaCalvert != 0) > 0]
GardenEpifaunaCalvert2

### Splitting the data
names(GardenEpifaunaCalvert2)
GardenEpifaunaCalvertStructure = GardenEpifaunaCalvert2[,1:8]
GardenEpifaunaCalvertData   = GardenEpifaunaCalvert2[,9:40]
GardenEpifaunaCalvertData2=(GardenEpifaunaCalvertData+1)

### SIMPER
SpeciesSIMPER<-(simper(GardenEpifaunaCalvertData2, GardenEpifaunaCalvertStructure$SiteType2, permutations=999))
summary(SpeciesSIMPER)
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$`Clam Garden_Garden Reference`[,1]), na.rm = TRUE)

### communities are 67.53% dissimilar 
summary(SpeciesSIMPER)
SIMPERsummary=SpeciesSIMPERSum$`Clam Garden_Garden Reference`

SIMPERsummary$averagePercent =(SIMPERsummary$average)*100
SIMPERsummary
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="ava"] = "GardenDensity"
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="avb"] = "ReferenceDensity"
SIMPERsummary$Directionality =(SIMPERsummary$GardenDensity) - (SIMPERsummary$ReferenceDensity)

SIMPERsummary$Directionality2[SIMPERsummary$Directionality > 0 ] <- "1"
SIMPERsummary$Directionality2[SIMPERsummary$Directionality < 0 ] <- "-1"  
SIMPERsummary$Directionality2=as.numeric(SIMPERsummary$Directionality2)
SIMPERsummary$Average.Dissim = (SIMPERsummary$averagePercent)*(SIMPERsummary$Directionality2)

GardenEpifaunaCalvertSIMPER=SIMPERsummary
#write.csv(GardenEpifaunaCalvertSIMPER,"GardenEpifaunaCalvertSIMPER.csv")
GardenEpifaunaCalvertSIMPER2=read.csv("GardenEpifaunaCalvertSIMPER.csv",header=T)
colnames(GardenEpifaunaCalvertSIMPER2)[colnames(GardenEpifaunaCalvertSIMPER2)=="X"] = "Species"
names(GardenEpifaunaCalvertSIMPER2)

## plot
#  Gardencol1,Gardencol2,Farmcol1,Farmcol2

GardenEpifaunaCalvertSIMPER2$Species <- factor(GardenEpifaunaCalvertSIMPER2$Species, levels = GardenEpifaunaCalvertSIMPER2$Species[order(GardenEpifaunaCalvertSIMPER2$Average.Dissim)])
GardenEpifaunaCalvertSIMPER2$Directionality2=as.factor(GardenEpifaunaCalvertSIMPER2$Directionality2)

GardenEpifaunaCalvertSIMPERPlot=ggplot(GardenEpifaunaCalvertSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-44,44)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Gardencol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
GardenEpifaunaCalvertSIMPERPlot

GardenEpifaunaCalvertSIMPERPlot2=GardenEpifaunaCalvertSIMPERPlot+geom_vline(xintercept = nlevels(GardenEpifaunaCalvertSIMPER2$Species)+0.6,size=0.6)
GardenEpifaunaCalvertSIMPERPlot3=GardenEpifaunaCalvertSIMPERPlot2+annotate("text", y= -40, x =nlevels(GardenEpifaunaCalvertSIMPER2$Species),label="Dissimilarity 67.53%",size = size2)
GardenEpifaunaCalvertSIMPERPlot4=GardenEpifaunaCalvertSIMPERPlot3 + theme(axis.text.y = italic.text)
GardenEpifaunaCalvertSIMPERPlot4

#### Epifauna Quadra Gardens SIMPER #### 
EpifaunaSpeciesData6

GardenEpifauna <-  subset(EpifaunaSpeciesData6, EpifaunaSpeciesData6$Comparison == "Garden")
GardenEpifaunaQuadra <- subset(GardenEpifauna, GardenEpifauna$Area == "Quadra")
GardenEpifaunaQuadra
GardenEpifaunaQuadra2=GardenEpifaunaQuadra[, colSums(GardenEpifaunaQuadra != 0) > 0]
GardenEpifaunaQuadra2

### Splitting the data
names(GardenEpifaunaQuadra2)
GardenEpifaunaQuadraStructure = GardenEpifaunaQuadra2[,1:8]
GardenEpifaunaQuadraData   = GardenEpifaunaQuadra2[,9:31]
GardenEpifaunaQuadraData2=(GardenEpifaunaQuadraData+1)

### SIMPER
SpeciesSIMPER<-(simper(GardenEpifaunaQuadraData2, GardenEpifaunaQuadraStructure$SiteType2, permutations=999))
summary(SpeciesSIMPER)
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$`Clam Garden_Garden Reference`[,1]), na.rm = TRUE)

### communities are 57.33% dissimilar 
summary(SpeciesSIMPER)
SIMPERsummary=SpeciesSIMPERSum$`Clam Garden_Garden Reference`

SIMPERsummary$averagePercent =(SIMPERsummary$average)*100
SIMPERsummary
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="ava"] = "GardenDensity"
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="avb"] = "ReferenceDensity"
SIMPERsummary$Directionality =(SIMPERsummary$GardenDensity) - (SIMPERsummary$ReferenceDensity)

SIMPERsummary$Directionality2[SIMPERsummary$Directionality > 0 ] <- "1"
SIMPERsummary$Directionality2[SIMPERsummary$Directionality < 0 ] <- "-1"  
SIMPERsummary$Directionality2=as.numeric(SIMPERsummary$Directionality2)
SIMPERsummary$Average.Dissim = (SIMPERsummary$averagePercent)*(SIMPERsummary$Directionality2)

GardenEpifaunaQuadraSIMPER=SIMPERsummary
#write.csv(GardenEpifaunaQuadraSIMPER,"GardenEpifaunaQuadraSIMPER.csv")
GardenEpifaunaQuadraSIMPER2=read.csv("GardenEpifaunaQuadraSIMPER.csv",header=T)
colnames(GardenEpifaunaQuadraSIMPER2)[colnames(GardenEpifaunaQuadraSIMPER2)=="X"] = "Species"
names(GardenEpifaunaQuadraSIMPER2)

## plot
#  Gardencol1,Gardencol2,Farmcol1,Farmcol2

GardenEpifaunaQuadraSIMPER2$Species <- factor(GardenEpifaunaQuadraSIMPER2$Species, levels = GardenEpifaunaQuadraSIMPER2$Species[order(GardenEpifaunaQuadraSIMPER2$Average.Dissim)])
GardenEpifaunaQuadraSIMPER2$Directionality2=as.factor(GardenEpifaunaQuadraSIMPER2$Directionality2)

GardenEpifaunaQuadraSIMPERPlot=ggplot(GardenEpifaunaQuadraSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-40,40)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Gardencol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
GardenEpifaunaQuadraSIMPERPlot

GardenEpifaunaQuadraSIMPERPlot2=GardenEpifaunaQuadraSIMPERPlot+geom_vline(xintercept = nlevels(GardenEpifaunaQuadraSIMPER2$Species)+0.6,size=0.6)
GardenEpifaunaQuadraSIMPERPlot3=GardenEpifaunaQuadraSIMPERPlot2+annotate("text", y= -40, x =nlevels(GardenEpifaunaQuadraSIMPER2$Species),label="Dissimilarity 57.33%",size = size2)
GardenEpifaunaQuadraSIMPERPlot4=GardenEpifaunaQuadraSIMPERPlot3 + theme(axis.text.y = italic.text)
GardenEpifaunaQuadraSIMPERPlot4

#### Epifauna Quadra Aquaculture SIMPER #### 

EpifaunaSpeciesData6
AquacultureEpifauna <- subset(EpifaunaSpeciesData6, EpifaunaSpeciesData6$Comparison == "Aquaculture")
AquacultureEpifaunaQuadra <- subset(AquacultureEpifauna, AquacultureEpifauna$Area == "Quadra")
AquacultureEpifaunaQuadra
AquacultureEpifaunaQuadra2=AquacultureEpifaunaQuadra[, colSums(AquacultureEpifaunaQuadra != 0) > 0]
AquacultureEpifaunaQuadra2

### Splitting the data
names(AquacultureEpifaunaQuadra2)
AquacultureEpifaunaQuadraStructure = AquacultureEpifaunaQuadra2[,1:8]
AquacultureEpifaunaQuadraData   = AquacultureEpifaunaQuadra2[,9:31]
AquacultureEpifaunaQuadraData2=(AquacultureEpifaunaQuadraData+1)

### SIMPER
SpeciesSIMPER<-(simper(AquacultureEpifaunaQuadraData2, AquacultureEpifaunaQuadraStructure$SiteType2, permutations=999))
summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$`Farm Reference_Shellfish Farm`[,1]), na.rm = TRUE)

### communities are 77.30% dissimilar 
summary(SpeciesSIMPER)
SIMPERsummary=SpeciesSIMPERSum$`Farm Reference_Shellfish Farm`

SIMPERsummary$averagePercent =(SIMPERsummary$average)*100
SIMPERsummary
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="ava"] = "ReferenceDensity"
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="avb"] = "FarmDensity"
SIMPERsummary$Directionality =(SIMPERsummary$FarmDensity) - (SIMPERsummary$ReferenceDensity)

SIMPERsummary$Directionality2[SIMPERsummary$Directionality > 0 ] <- "1"
SIMPERsummary$Directionality2[SIMPERsummary$Directionality < 0 ] <- "-1"  
SIMPERsummary$Directionality2=as.numeric(SIMPERsummary$Directionality2)
SIMPERsummary$Average.Dissim = (SIMPERsummary$averagePercent)*(SIMPERsummary$Directionality2)

AquacultureEpifaunaQuadraSIMPER=SIMPERsummary
#write.csv(AquacultureEpifaunaQuadraSIMPER,"AquacultureEpifaunaSIMPER.csv")
AquacultureEpifaunaQuadraSIMPER2=read.csv("AquacultureEpifaunaSIMPER.csv",header=T)
colnames(AquacultureEpifaunaQuadraSIMPER2)[colnames(AquacultureEpifaunaQuadraSIMPER2)=="X"] = "Species"
names(AquacultureEpifaunaQuadraSIMPER2)

## plot
#  Gardencol1,Gardencol2,Farmcol1,Farmcol2

AquacultureEpifaunaQuadraSIMPER2$Species <- factor(AquacultureEpifaunaQuadraSIMPER2$Species, levels = AquacultureEpifaunaQuadraSIMPER2$Species[order(AquacultureEpifaunaQuadraSIMPER2$Average.Dissim)])
AquacultureEpifaunaQuadraSIMPER2$Directionality2=as.factor(AquacultureEpifaunaQuadraSIMPER2$Directionality2)

AquacultureEpifaunaQuadraSIMPERPlot=ggplot(AquacultureEpifaunaQuadraSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-50,50)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Farmcol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
AquacultureEpifaunaQuadraSIMPERPlot

AquacultureEpifaunaQuadraSIMPERPlot2=AquacultureEpifaunaQuadraSIMPERPlot+geom_vline(xintercept = nlevels(AquacultureEpifaunaQuadraSIMPER2$Species)+0.6,size=0.6)
AquacultureEpifaunaQuadraSIMPERPlot3=AquacultureEpifaunaQuadraSIMPERPlot2+annotate("text", y= -40, x =nlevels(AquacultureEpifaunaQuadraSIMPER2$Species),label="Dissimilarity 77.30%",size = size2)
AquacultureEpifaunaQuadraSIMPERPlot4=AquacultureEpifaunaQuadraSIMPERPlot3 + theme(axis.text.y = italic.text)
AquacultureEpifaunaQuadraSIMPERPlot4


#### Epifauna Baynes Aquaculture SIMPER #### 
EpifaunaSpeciesData6
AquacultureEpifauna <- subset(EpifaunaSpeciesData6, EpifaunaSpeciesData6$Comparison == "Aquaculture")
AquacultureEpifaunaBaynes <- subset(AquacultureEpifauna, AquacultureEpifauna$Area == "Baynes")
AquacultureEpifaunaBaynes
AquacultureEpifaunaBaynes2=AquacultureEpifaunaBaynes[, colSums(AquacultureEpifaunaBaynes != 0) > 0]
AquacultureEpifaunaBaynes2

### Splitting the data
names(AquacultureEpifaunaBaynes2)
AquacultureEpifaunaBaynesStructure = AquacultureEpifaunaBaynes2[,1:8]
AquacultureEpifaunaBaynesData   = AquacultureEpifaunaBaynes2[,9:34]
AquacultureEpifaunaBaynesData2=(AquacultureEpifaunaBaynesData+1)

### SIMPER
SpeciesSIMPER<-(simper(AquacultureEpifaunaBaynesData2, AquacultureEpifaunaBaynesStructure$SiteType2, permutations=999))
summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$`Farm Reference_Shellfish Farm`[,1]), na.rm = TRUE)

### communities are 72.14% dissimilar 
summary(SpeciesSIMPER)
SIMPERsummary=SpeciesSIMPERSum$`Farm Reference_Shellfish Farm`

SIMPERsummary$averagePercent =(SIMPERsummary$average)*100
SIMPERsummary
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="ava"] = "ReferenceDensity"
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="avb"] = "FarmDensity"
SIMPERsummary$Directionality =(SIMPERsummary$FarmDensity) - (SIMPERsummary$ReferenceDensity)

SIMPERsummary$Directionality2[SIMPERsummary$Directionality > 0 ] <- "1"
SIMPERsummary$Directionality2[SIMPERsummary$Directionality < 0 ] <- "-1"  
SIMPERsummary$Directionality2=as.numeric(SIMPERsummary$Directionality2)
SIMPERsummary$Average.Dissim = (SIMPERsummary$averagePercent)*(SIMPERsummary$Directionality2)

AquacultureEpifaunaBaynesSIMPER=SIMPERsummary
#write.csv(AquacultureEpifaunaBaynesSIMPER,"AquacultureEpifaunaBaynesSIMPER.csv")
AquacultureEpifaunaBaynesSIMPER2=read.csv("AquacultureEpifaunaBaynesSIMPER.csv",header=T)
colnames(AquacultureEpifaunaBaynesSIMPER2)[colnames(AquacultureEpifaunaBaynesSIMPER2)=="X"] = "Species"
names(AquacultureEpifaunaBaynesSIMPER2)

## plot
#  Gardencol1,Gardencol2,Farmcol1,Farmcol2

AquacultureEpifaunaBaynesSIMPER2$Species <- factor(AquacultureEpifaunaBaynesSIMPER2$Species, levels = AquacultureEpifaunaBaynesSIMPER2$Species[order(AquacultureEpifaunaBaynesSIMPER2$Average.Dissim)])
AquacultureEpifaunaBaynesSIMPER2$Directionality2=as.factor(AquacultureEpifaunaBaynesSIMPER2$Directionality2)

AquacultureEpifaunaBaynesSIMPERPlot=ggplot(AquacultureEpifaunaBaynesSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-42,42)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Farmcol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
AquacultureEpifaunaBaynesSIMPERPlot

AquacultureEpifaunaBaynesSIMPERPlot2=AquacultureEpifaunaBaynesSIMPERPlot+geom_vline(xintercept = nlevels(AquacultureEpifaunaBaynesSIMPER2$Species)+0.6,size=0.6)
AquacultureEpifaunaBaynesSIMPERPlot3=AquacultureEpifaunaBaynesSIMPERPlot2+annotate("text", y= -35, x =nlevels(AquacultureEpifaunaBaynesSIMPER2$Species),label="Dissimilarity 72.14%",size = size2)
AquacultureEpifaunaBaynesSIMPERPlot4=AquacultureEpifaunaBaynesSIMPERPlot3 + theme(axis.text.y = italic.text)
AquacultureEpifaunaBaynesSIMPERPlot4


#### Aquacultures and Gardens ####
EpifaunaSpeciesData6

EpifaunaGardenFarm=(subset(EpifaunaSpeciesData6, SiteType2 %in% c("Clam Garden","Shellfish Farm")))
EpifaunaGardenFarm$SiteType2

EpifaunaGardenFarm2=EpifaunaGardenFarm[, colSums(EpifaunaGardenFarm != 0) > 0]
EpifaunaGardenFarm2

#### Epifauna Aquacultures and Gardens ####
EpifaunaGardenFarm2
names(EpifaunaGardenFarm2)
EpifaunaGardenFarmStructure = EpifaunaGardenFarm2[,1:8]
EpifaunaGardenFarmData   = EpifaunaGardenFarm2[,9:52]
EpifaunaGardenFarmData2=(EpifaunaGardenFarmData+1)

### SIMPER
SpeciesSIMPER<-(simper(EpifaunaGardenFarmData2, EpifaunaGardenFarmStructure$SiteType, permutations=999))
summary(SpeciesSIMPER) 
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$`Farm_Clam Garden`[,1]), na.rm = TRUE)
SIMPERscore=sum(as.numeric(SpeciesSIMPERSum$`Farm_Clam Garden`[,1]), na.rm = TRUE)

### communities are 71.80% dissimilar 
summary(SpeciesSIMPER) 
SIMPERsummary=SpeciesSIMPERSum$`Farm_Clam Garden`

SIMPERsummary$averagePercent =(SIMPERsummary$average)*100
SIMPERsummary$Directionality =(SIMPERsummary$ava) - (SIMPERsummary$avb)

colnames(SIMPERsummary)[colnames(SIMPERsummary)=="ava"] = "FarmDensity"
colnames(SIMPERsummary)[colnames(SIMPERsummary)=="avb"] = "GardenDensity"
SIMPERsummary$Directionality =(SIMPERsummary$FarmDensity) - (SIMPERsummary$GardenDensity)
names(SIMPERsummary)
SIMPERsummary$Directionality2[SIMPERsummary$Directionality > 0 ] <- "1"
SIMPERsummary$Directionality2[SIMPERsummary$Directionality < 0 ] <- "-1"  
SIMPERsummary$Directionality2=as.numeric(SIMPERsummary$Directionality2)
SIMPERsummary$Average.Dissim = (SIMPERsummary$averagePercent)*(SIMPERsummary$Directionality2)

EpifaunaGardenFarmSIMPER=SIMPERsummary
#write.csv(EpifaunaGardenFarmSIMPER,"EpifaunaGardenFarmSIMPER.csv")
EpifaunaGardenFarmSIMPER2=read.csv("EpifaunaGardenFarmSIMPER.csv",header=T)
colnames(EpifaunaGardenFarmSIMPER2)[colnames(EpifaunaGardenFarmSIMPER2)=="X"] = "Species"
names(EpifaunaGardenFarmSIMPER2)

## plot
GardenFarmcol2=c(FarmColor,GardenColor)
GardenFarmcol1=c(FarmColor,GardenColor)

EpifaunaGardenFarmSIMPER2$Species <- factor(EpifaunaGardenFarmSIMPER2$Species, levels = EpifaunaGardenFarmSIMPER2$Species[order(EpifaunaGardenFarmSIMPER2$Average.Dissim)])
EpifaunaGardenFarmSIMPER2$Directionality2=as.factor(EpifaunaGardenFarmSIMPER2$Directionality2)

EpifaunaGardenFarmSIMPERPlot=ggplot(EpifaunaGardenFarmSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-43,43)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=GardenFarmcol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
EpifaunaGardenFarmSIMPERPlot

EpifaunaGardenFarmSIMPERPlot2=EpifaunaGardenFarmSIMPERPlot+geom_vline(xintercept = nlevels(EpifaunaGardenFarmSIMPER2$Species)+0.6,size=0.6)
EpifaunaGardenFarmSIMPERPlot3=EpifaunaGardenFarmSIMPERPlot2+annotate("text", y= -40, x =nlevels(EpifaunaGardenFarmSIMPER2$Species),label="Dissimilarity 71.80%",size = size2)
EpifaunaGardenFarmSIMPERPlot4=EpifaunaGardenFarmSIMPERPlot3 + theme(axis.text.y = italic.text)
EpifaunaGardenFarmSIMPERPlot4

#### SIMPER All Plots Code ####

## Epifauna Garden Summary
GardenEpifaunaSIMPERPlot=ggplot(GardenEpifaunaSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-50,50)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Gardencol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
GardenEpifaunaSIMPERPlot

GardenEpifaunaSIMPERPlot2=GardenEpifaunaSIMPERPlot+geom_vline(xintercept = 40.6,size=0.6)
GardenEpifaunaSIMPERPlot3=GardenEpifaunaSIMPERPlot2+annotate("text", y= -40, x =nlevels(GardenEpifaunaSIMPER2$Species),label="Dissimilarity 67.66%",size = size2)
GardenEpifaunaSIMPERPlot4=GardenEpifaunaSIMPERPlot3 + theme(axis.text.y = italic.text)
GardenEpifaunaSIMPERPlot4

## Epifauna Aquaculture Summary
AquacultureEpifaunaSIMPERPlot=ggplot(AquacultureEpifaunaSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-50,50)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Farmcol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
AquacultureEpifaunaSIMPERPlot

AquacultureEpifaunaSIMPERPlot2=AquacultureEpifaunaSIMPERPlot+geom_vline(xintercept = nlevels(AquacultureEpifaunaSIMPER2$Species)+0.6,size=0.6)
AquacultureEpifaunaSIMPERPlot3=AquacultureEpifaunaSIMPERPlot2+annotate("text", y= -40, x =nlevels(AquacultureEpifaunaSIMPER2$Species),label="Dissimilarity 74.86%",size = size2)
AquacultureEpifaunaSIMPERPlot4=AquacultureEpifaunaSIMPERPlot3 + theme(axis.text.y = italic.text)
AquacultureEpifaunaSIMPERPlot4

## Epifauna Calvert Garden
GardenEpifaunaCalvertSIMPERPlot=ggplot(GardenEpifaunaCalvertSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-44,44)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Gardencol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
GardenEpifaunaCalvertSIMPERPlot

GardenEpifaunaCalvertSIMPERPlot2=GardenEpifaunaCalvertSIMPERPlot+geom_vline(xintercept = nlevels(GardenEpifaunaCalvertSIMPER2$Species)+0.6,size=0.6)
GardenEpifaunaCalvertSIMPERPlot3=GardenEpifaunaCalvertSIMPERPlot2+annotate("text", y= -40, x =nlevels(GardenEpifaunaCalvertSIMPER2$Species),label="Dissimilarity 67.53%",size = size2)
GardenEpifaunaCalvertSIMPERPlot4=GardenEpifaunaCalvertSIMPERPlot3 + theme(axis.text.y = italic.text)
GardenEpifaunaCalvertSIMPERPlot4

## Epifauna Quadra Gardens
GardenEpifaunaQuadraSIMPERPlot=ggplot(GardenEpifaunaQuadraSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-40,40)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Gardencol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
GardenEpifaunaQuadraSIMPERPlot

GardenEpifaunaQuadraSIMPERPlot2=GardenEpifaunaQuadraSIMPERPlot+geom_vline(xintercept = nlevels(GardenEpifaunaQuadraSIMPER2$Species)+0.6,size=0.6)
GardenEpifaunaQuadraSIMPERPlot3=GardenEpifaunaQuadraSIMPERPlot2+annotate("text", y= -40, x =nlevels(GardenEpifaunaQuadraSIMPER2$Species),label="Dissimilarity 57.33%",size = size2)
GardenEpifaunaQuadraSIMPERPlot4=GardenEpifaunaQuadraSIMPERPlot3 + theme(axis.text.y = italic.text)
GardenEpifaunaQuadraSIMPERPlot4

## Epifauna Quadra Aquaculture
AquacultureEpifaunaQuadraSIMPERPlot=ggplot(AquacultureEpifaunaQuadraSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-50,50)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Farmcol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
AquacultureEpifaunaQuadraSIMPERPlot

AquacultureEpifaunaQuadraSIMPERPlot2=AquacultureEpifaunaQuadraSIMPERPlot+geom_vline(xintercept = nlevels(AquacultureEpifaunaQuadraSIMPER2$Species)+0.6,size=0.6)
AquacultureEpifaunaQuadraSIMPERPlot3=AquacultureEpifaunaQuadraSIMPERPlot2+annotate("text", y= -40, x =nlevels(AquacultureEpifaunaQuadraSIMPER2$Species),label="Dissimilarity 77.30%",size = size2)
AquacultureEpifaunaQuadraSIMPERPlot4=AquacultureEpifaunaQuadraSIMPERPlot3 + theme(axis.text.y = italic.text)
AquacultureEpifaunaQuadraSIMPERPlot4

## Epifauna Baynes  Aquaculture
AquacultureEpifaunaBaynesSIMPERPlot=ggplot(AquacultureEpifaunaBaynesSIMPER2, aes(x = Species, y = Average.Dissim,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity")+ylim(-42,42)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Farmcol1)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
AquacultureEpifaunaBaynesSIMPERPlot

AquacultureEpifaunaBaynesSIMPERPlot2=AquacultureEpifaunaBaynesSIMPERPlot+geom_vline(xintercept = nlevels(AquacultureEpifaunaBaynesSIMPER2$Species)+0.6,size=0.6)
AquacultureEpifaunaBaynesSIMPERPlot3=AquacultureEpifaunaBaynesSIMPERPlot2+annotate("text", y= -35, x =nlevels(AquacultureEpifaunaBaynesSIMPER2$Species),label="Dissimilarity 72.14%",size = size2)
AquacultureEpifaunaBaynesSIMPERPlot4=AquacultureEpifaunaBaynesSIMPERPlot3 + theme(axis.text.y = italic.text)
AquacultureEpifaunaBaynesSIMPERPlot4

#### Contructing Figure ####
Image5 <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Thesis Chapter, Nov 2020/Analysis/Results /SIMPER/Figure Construction/LegendV5.tiff")
Image5
Legend5 <- image_ggplot(Image5)

Image6 <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Thesis Chapter, Nov 2020/Analysis/Results /SIMPER/Figure Construction/SIMPERLegendV6.tiff")
Image6
Legend6 <- image_ggplot(Image6)

Image8 <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Thesis Chapter, Nov 2020/Analysis/Results /SIMPER/Figure Construction/SIMPERLegendV8.tiff")
Image8
Legend8 <- image_ggplot(Image8)


Adistance=-6.52
CFdistance=-5.1

size3=14
GJDistance=-4.5

## Epifauna Garden Summary
## Transform Axis ##
GardenEpifaunaSIMPER2$Average.Dissim
sqrt(GardenEpifaunaSIMPER2$averagePercent)

GardenEpifaunaSIMPER2$Directionality3[GardenEpifaunaSIMPER2$Directionality > 0 ] <- "1"
GardenEpifaunaSIMPER2$Directionality3[GardenEpifaunaSIMPER2$Directionality < 0 ] <- "-1"  

GardenEpifaunaSIMPER2$Directionality3=as.numeric(GardenEpifaunaSIMPER2$Directionality3)
GardenEpifaunaSIMPER2$Average.Dissim2 = ((sqrt(GardenEpifaunaSIMPER2$averagePercent)) * (GardenEpifaunaSIMPER2$Directionality3) )
## Plot
GardenEpifaunaSIMPERPlot=ggplot(GardenEpifaunaSIMPER2, aes(x = Species, y = Average.Dissim2,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab(expression(sqrt("Average Dissimilarity")))+ylim(-7,7)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Gardencol1)+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
GardenEpifaunaSIMPERPlot

GardenEpifaunaSIMPERPlot2=GardenEpifaunaSIMPERPlot+geom_vline(xintercept = 40.6,size=0.6)
GardenEpifaunaSIMPERPlot3=GardenEpifaunaSIMPERPlot2+annotate("text", y= -6.51, x =nlevels(GardenEpifaunaSIMPER2$Species),label="Dissimilarity 67.66%",size = size2)
GardenEpifaunaSIMPERPlot4=GardenEpifaunaSIMPERPlot3 + theme(axis.text.y = italic.text)
GardenEpifaunaSIMPERPlot4

## Epifauna Aquaculture Summary
## Transform Axis ##
AquacultureEpifaunaSIMPER2$Average.Dissim
sqrt(AquacultureEpifaunaSIMPER2$averagePercent)

AquacultureEpifaunaSIMPER2$Directionality3[AquacultureEpifaunaSIMPER2$Directionality > 0 ] <- "1"
AquacultureEpifaunaSIMPER2$Directionality3[AquacultureEpifaunaSIMPER2$Directionality < 0 ] <- "-1"  

AquacultureEpifaunaSIMPER2$Directionality3=as.numeric(AquacultureEpifaunaSIMPER2$Directionality3)
AquacultureEpifaunaSIMPER2$Average.Dissim2 = ((sqrt(AquacultureEpifaunaSIMPER2$averagePercent)) * (AquacultureEpifaunaSIMPER2$Directionality3) )

## Plot
AquacultureEpifaunaSIMPERPlot=ggplot(AquacultureEpifaunaSIMPER2, aes(x = Species, y = Average.Dissim2,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab(expression(sqrt("Average Dissimilarity")))+ylim(-7,7)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Farmcol1)+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
AquacultureEpifaunaSIMPERPlot

AquacultureEpifaunaSIMPERPlot2=AquacultureEpifaunaSIMPERPlot+geom_vline(xintercept = nlevels(AquacultureEpifaunaSIMPER2$Species)+0.6,size=0.6)
AquacultureEpifaunaSIMPERPlot3=AquacultureEpifaunaSIMPERPlot2+annotate("text", y= -6.47, x =nlevels(AquacultureEpifaunaSIMPER2$Species),label="Dissimilarity 74.86%",size = size2)
AquacultureEpifaunaSIMPERPlot4=AquacultureEpifaunaSIMPERPlot3 + theme(axis.text.y = italic.text)
AquacultureEpifaunaSIMPERPlot4

## Epifauna Calvert Garden
## Transform Axis ##
GardenEpifaunaCalvertSIMPER2$Directionality3[GardenEpifaunaCalvertSIMPER2$Directionality > 0 ] <- "1"
GardenEpifaunaCalvertSIMPER2$Directionality3[GardenEpifaunaCalvertSIMPER2$Directionality < 0 ] <- "-1"  

GardenEpifaunaCalvertSIMPER2$Directionality3=as.numeric(GardenEpifaunaCalvertSIMPER2$Directionality3)
GardenEpifaunaCalvertSIMPER2$Average.Dissim2 = ((sqrt(GardenEpifaunaCalvertSIMPER2$averagePercent)) * (GardenEpifaunaCalvertSIMPER2$Directionality3) )

## Plot
GardenEpifaunaCalvertSIMPERPlot=ggplot(GardenEpifaunaCalvertSIMPER2, aes(x = Species, y = Average.Dissim2,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab(expression(sqrt("Average Dissimilarity")))+ylim(-7,7)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Gardencol1)+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
GardenEpifaunaCalvertSIMPERPlot

GardenEpifaunaCalvertSIMPERPlot2=GardenEpifaunaCalvertSIMPERPlot+geom_vline(xintercept = nlevels(GardenEpifaunaCalvertSIMPER2$Species)+0.6,size=0.6)
GardenEpifaunaCalvertSIMPERPlot3=GardenEpifaunaCalvertSIMPERPlot2+annotate("text", y= -4.96, x =nlevels(GardenEpifaunaCalvertSIMPER2$Species),label="Dissimilarity 67.53%",size = size2)
GardenEpifaunaCalvertSIMPERPlot4=GardenEpifaunaCalvertSIMPERPlot3 + theme(axis.text.y = italic.text)
GardenEpifaunaCalvertSIMPERPlot4


## Epifauna Quadra Gardens
## Transform Axis ##
GardenEpifaunaQuadraSIMPER2$Directionality3[GardenEpifaunaQuadraSIMPER2$Directionality > 0 ] <- "1"
GardenEpifaunaQuadraSIMPER2$Directionality3[GardenEpifaunaQuadraSIMPER2$Directionality < 0 ] <- "-1"  

GardenEpifaunaQuadraSIMPER2$Directionality3=as.numeric(GardenEpifaunaQuadraSIMPER2$Directionality3)
GardenEpifaunaQuadraSIMPER2$Average.Dissim2 = ((sqrt(GardenEpifaunaQuadraSIMPER2$averagePercent)) * (GardenEpifaunaQuadraSIMPER2$Directionality3) )

## Plot
GardenEpifaunaQuadraSIMPERPlot=ggplot(GardenEpifaunaQuadraSIMPER2, aes(x = Species, y = Average.Dissim2,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab(expression(sqrt("Average Dissimilarity")))+ylim(-7,7)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Gardencol1)+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
GardenEpifaunaQuadraSIMPERPlot

GardenEpifaunaQuadraSIMPERPlot2=GardenEpifaunaQuadraSIMPERPlot+geom_vline(xintercept = nlevels(GardenEpifaunaQuadraSIMPER2$Species)+0.6,size=0.6)
GardenEpifaunaQuadraSIMPERPlot3=GardenEpifaunaQuadraSIMPERPlot2+annotate("text", y= -4.96, x =nlevels(GardenEpifaunaQuadraSIMPER2$Species),label="Dissimilarity 57.33%",size = size2)
GardenEpifaunaQuadraSIMPERPlot4=GardenEpifaunaQuadraSIMPERPlot3 + theme(axis.text.y = italic.text)
GardenEpifaunaQuadraSIMPERPlot4

## Epifauna Quadra Aquaculture
## Transform Axis ##
AquacultureEpifaunaQuadraSIMPER2$Directionality3[AquacultureEpifaunaQuadraSIMPER2$Directionality > 0 ] <- "1"
AquacultureEpifaunaQuadraSIMPER2$Directionality3[AquacultureEpifaunaQuadraSIMPER2$Directionality < 0 ] <- "-1"  

AquacultureEpifaunaQuadraSIMPER2$Directionality3=as.numeric(AquacultureEpifaunaQuadraSIMPER2$Directionality3)
AquacultureEpifaunaQuadraSIMPER2$Average.Dissim2 = ((sqrt(AquacultureEpifaunaQuadraSIMPER2$averagePercent)) * (AquacultureEpifaunaQuadraSIMPER2$Directionality3) )

## Plot
AquacultureEpifaunaQuadraSIMPERPlot=ggplot(AquacultureEpifaunaQuadraSIMPER2, aes(x = Species, y = Average.Dissim2,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab(expression(sqrt("Average Dissimilarity")))+ylim(-7,7)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Farmcol1)+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
AquacultureEpifaunaQuadraSIMPERPlot

AquacultureEpifaunaQuadraSIMPERPlot2=AquacultureEpifaunaQuadraSIMPERPlot+geom_vline(xintercept = nlevels(AquacultureEpifaunaQuadraSIMPER2$Species)+0.6,size=0.6)
AquacultureEpifaunaQuadraSIMPERPlot3=AquacultureEpifaunaQuadraSIMPERPlot2+annotate("text", y= -4.6, x =nlevels(AquacultureEpifaunaQuadraSIMPER2$Species),label="Dissimilarity 77.30%",size = size2)
AquacultureEpifaunaQuadraSIMPERPlot4=AquacultureEpifaunaQuadraSIMPERPlot3 + theme(axis.text.y = italic.text)
AquacultureEpifaunaQuadraSIMPERPlot4


## Epifauna Baynes  Aquaculture
## Transform Axis ##
AquacultureEpifaunaBaynesSIMPER2$Directionality3[AquacultureEpifaunaBaynesSIMPER2$Directionality > 0 ] <- "1"
AquacultureEpifaunaBaynesSIMPER2$Directionality3[AquacultureEpifaunaBaynesSIMPER2$Directionality < 0 ] <- "-1"  

AquacultureEpifaunaBaynesSIMPER2$Directionality3=as.numeric(AquacultureEpifaunaBaynesSIMPER2$Directionality3)
AquacultureEpifaunaBaynesSIMPER2$Average.Dissim2 = ((sqrt(AquacultureEpifaunaBaynesSIMPER2$averagePercent)) * (AquacultureEpifaunaBaynesSIMPER2$Directionality3) )

## Plot
AquacultureEpifaunaBaynesSIMPERPlot=ggplot(AquacultureEpifaunaBaynesSIMPER2, aes(x = Species, y = Average.Dissim2,colour=Directionality2))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab(expression(sqrt("Average Dissimilarity")))+ylim(-7,7)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=Farmcol1)+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)+coord_flip()
AquacultureEpifaunaBaynesSIMPERPlot

AquacultureEpifaunaBaynesSIMPERPlot2=AquacultureEpifaunaBaynesSIMPERPlot+geom_vline(xintercept = nlevels(AquacultureEpifaunaBaynesSIMPER2$Species)+0.6,size=0.6)
AquacultureEpifaunaBaynesSIMPERPlot3=AquacultureEpifaunaBaynesSIMPERPlot2+annotate("text", y= -4.96, x =nlevels(AquacultureEpifaunaBaynesSIMPER2$Species),label="Dissimilarity 72.14%",size = size2)
AquacultureEpifaunaBaynesSIMPERPlot4=AquacultureEpifaunaBaynesSIMPERPlot3 + theme(axis.text.y = italic.text)
AquacultureEpifaunaBaynesSIMPERPlot4


## plots
#A
GardenEpifaunaSIMPERPlot4
AquacultureEpifaunaSIMPERPlot4

#B-
GardenEpifaunaCalvertSIMPERPlot4
GardenEpifaunaQuadraSIMPERPlot4
AquacultureEpifaunaQuadraSIMPERPlot4
AquacultureEpifaunaBaynesSIMPERPlot4

margins1 = unit(c(  0.45  ,0.2, 0.5 ,0.5), "cm")   #A
margins2 = unit(c(  0.45  ,0.2, 0.5 ,0.5), "cm") #B
margins5 = unit(c(  0.15  ,0.2, 0,    0.5), "cm")  #C
margins6 = unit(c(  0.15  ,0.2, 0,    0.5), "cm")  #D
margins7 = unit(c(  0.15  ,0.2, 0,    0.5), "cm")  #E
margins8 = unit(c(  0.15  ,0.2, 0,    0.5), "cm")  #F
#       top, right,bottom,left

#library('addinslist')

## Figure 4 Reduced SIMPER Plot 

GardenEpifaunaSIMPERPlot2=GardenEpifaunaSIMPERPlot+geom_vline(xintercept = 40.6,size=0.6)
GardenEpifaunaSIMPERPlot3=GardenEpifaunaSIMPERPlot2+annotate("text", y= -6.19, x =nlevels(GardenEpifaunaSIMPER2$Species),label="Dissimilarity 67.66%",size = size2)
GardenEpifaunaSIMPERPlot4=GardenEpifaunaSIMPERPlot3 + theme(axis.text.y = italic.text)
GardenEpifaunaSIMPERPlot4


AquacultureEpifaunaSIMPERPlot2=AquacultureEpifaunaSIMPERPlot+geom_vline(xintercept = nlevels(AquacultureEpifaunaSIMPER2$Species)+0.6,size=0.6)
AquacultureEpifaunaSIMPERPlot3=AquacultureEpifaunaSIMPERPlot2+annotate("text", y= -6.1, x =nlevels(AquacultureEpifaunaSIMPER2$Species),label="Dissimilarity 74.86%",size = size2)
AquacultureEpifaunaSIMPERPlot4=AquacultureEpifaunaSIMPERPlot3 + theme(axis.text.y = italic.text)
AquacultureEpifaunaSIMPERPlot4

LabelSize=21
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Thesis Chapter, Nov 2020/Analysis/Results /SIMPER/Figure Construction")
jpeg(filename = "EpifaunalSIMPERFINALFinalReduced4.jpeg", width = 58, height = 30, units = "cm", pointsize = 15, quality = 100, res = 600)  


tiff(filename = "EpifaunalSIMPERFINALFinalReduced5.tiff",
     width = 58, height = 30, units = "cm", pointsize = 15,
     compression = c("none"),
     bg = "white", res = 600)


grid.arrange(GardenEpifaunaSIMPERPlot4+theme(plot.margin =margins1)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black", hjust = -0.1)) +labs(title = NULL, x = NULL, subtitle = "A)"), 
             
             AquacultureEpifaunaSIMPERPlot4+theme(plot.margin =margins2)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black", hjust = -0.1)) +labs(title = NULL, x = NULL, subtitle = "B)"),
             Legend8,
             heights=c(1.3,0.08),      
             widths = c(1,1,1,1),
             layout_matrix = rbind(c(1, 1, 2, 2),
                                   c(3,3,3,3)))
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")


## Full Simper Plot, for supplemental
LabelSize=21
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Thesis Chapter, Nov 2020/Analysis/Results /SIMPER/Figure Construction")
jpeg(filename = "EpifaunalSIMPERFINALFinal.jpeg", width = 72, height = 50, units = "cm", pointsize = 15, quality = 100, res = 300)  

grid.arrange(GardenEpifaunaSIMPERPlot4+theme(plot.margin =margins1)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black", hjust = -0.1)) +labs(title = NULL, x = NULL, subtitle = "A)"), 
             
             AquacultureEpifaunaSIMPERPlot4+theme(plot.margin =margins2)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black", hjust = -0.1)) +labs(title = NULL, x = NULL, subtitle = "B)"),
             
             GardenEpifaunaCalvertSIMPERPlot4+theme(plot.margin =margins5)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black", hjust = -0.28)) +labs(title = NULL, x = NULL, subtitle = "C)"),
             GardenEpifaunaQuadraSIMPERPlot4+theme(plot.margin =margins6)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black", hjust = -0.28)) +labs(title = NULL, x = NULL, subtitle = "D)"),
             AquacultureEpifaunaQuadraSIMPERPlot4+theme(plot.margin =margins7)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black", hjust = -0.28)) +labs(title = NULL, x = NULL, subtitle = "E)"),
             AquacultureEpifaunaBaynesSIMPERPlot4+theme(plot.margin =margins8)+theme(plot.subtitle = element_text(size = LabelSize, colour = "black", hjust = -0.28)) +labs(title = NULL, x = NULL, subtitle = "F)"),
             Legend8,
             heights=c(1.3,1,0.08),      
             widths = c(1,1,1,1),
             layout_matrix = rbind(c(1, 1, 2, 2),
                                   c(3, 4, 5, 6),
                                   c(7,7,7,7)))
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")

#### Done




















# KIERAN SAID TO START HERE ------
#### Figure 4: Regression Tree's Benthic Diversity ####

#### Epifauna Multivaraite Regression Trees with reduced data ####
#### Manuscript Approach
#### Setting up with reduced data ####
rm(list = ls(all = TRUE))
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")
EpifaunaDataRegression=read.csv("SFandCGBiodiversityDataNov2020 RegressionTree Family V8.csv", header=T)
EpifaunaDataRegression

#Packages 
library(tidyverse)
library(devtools)
library(mvpart)
library(MVPARTwrap)
library(ggplot2)
library(Hmisc)
library(ggsci)


#### reduced taxa data format ###
names(EpifaunaDataRegression)   ## reduced ~60 taxa to 24 taxa
str(EpifaunaDataRegression[1:48])
names(EpifaunaDataRegression[25:48]) #24 taxa

## Reordering taxa  
EpifaunaDataRegression2=EpifaunaDataRegression

col_order <- c("Amphilepidida", "Amphipoda", "Anthozoa","Asteriidae","Batillariidae","Bryozoa",
               "Cardiida","Cerithiidae","Chitonida","Decapod","Dendrasteridae","Haminoeidae",
               "Littorinidae","Lottiidae","Myrrhinidae","Mytilidae","Nemertea","Neogastropoda",
               "Ostreidae","Polychaeta","Sessilia","Sphaeromatidae","Strongylocentrotidae","Venerida")

EpifaunaDataRegression2 <- EpifaunaDataRegression2[,col_order]
EpifaunaDataRegression2
names(EpifaunaDataRegression2)
names(EpifaunaDataRegression[1:24])

EpifaunaDataRegression4=cbind(EpifaunaDataRegression[1:24],EpifaunaDataRegression2)
names(EpifaunaDataRegression4)
EpifaunaDataRegression4$Strongylocentrotidae
EpifaunaDataRegression$Strongylocentrotidae

colnames(EpifaunaDataRegression4)[colnames(EpifaunaDataRegression4)=="Boulder..25cm.490cm2."]="Boulder"
colnames(EpifaunaDataRegression4)[colnames(EpifaunaDataRegression4)=="Cobble.6.5.25.cm..33.18cm2."]="Cobble"
colnames(EpifaunaDataRegression4)[colnames(EpifaunaDataRegression4)=="Gravel...2cm...6.5cm."]="Gravel"
colnames(EpifaunaDataRegression4)[colnames(EpifaunaDataRegression4)=="Sand...0125cm....2cm."]="Sand"

GardenData <- subset(EpifaunaDataRegression4, EpifaunaDataRegression4$SiteType == "Clam Garden")
FarmData <- subset(EpifaunaDataRegression4, EpifaunaDataRegression4$SiteType == "Farm")
ReferenceData <- subset(EpifaunaDataRegression4, EpifaunaDataRegression4$SiteType == "Reference")

##Removing zeros
GardenData2=GardenData[, colSums(GardenData != 0) > 0]
FarmData2=FarmData[, colSums(FarmData != 0) > 0]
ReferenceData2=ReferenceData[, colSums(ReferenceData != 0) > 0]
EpifaunaDataRegression2=EpifaunaDataRegression4[, colSums(EpifaunaDataRegression4 != 0) > 0]

## datasets
GardenData2
FarmData2
ReferenceData2
EpifaunaDataRegression2
names(EpifaunaDataRegression2)
EpifaunaDataRegression=EpifaunaDataRegression2

#### All Data: Multivaraite Regression ####
### Workflow development  

names(EpifaunaDataRegression)
### Form  
## Transform, log10 + 1 dummy to reduce the impact of hyper abundant species 
FormAllData = (log10(data.matrix(EpifaunaDataRegression[,24:47])+1))~ (Zone+Cobble+Gravel+Sand+Mud.silt+Bivalve.Shells+Oyster.Shell+Sanddollar.dead+Gracilaria+Brown.Algae+Eelgrass+Ulva+Sargassum+Fuchus+Mastocarpus)

#create MRT and plot tree performance vs. tree size to determine best tree size
#set very small complexity parameter to see how error evolves with tree size
stop.rule = 0.001

mrt.error = mvpart(FormAllData, data = EpifaunaDataRegression, pretty = T, xv = "min", minauto =F,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                   xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = stop.rule, xval = 100), #100-fold cross validation
                   plot.add = F)

plotcp(mrt.error, upper = "size") # size is for # leaves, splits is for # splits
points(mrt.error$cptable[,"rel error"], col = "red")
legend(x = "topright", legend = c("rel error", "xvalidated rel error"), 
       col = c("red", "blue"), pch = c(1, 1))

full.table = as.data.frame(mrt.error$cptable)
print(full.table)

num.leaves = 7
final.cp = 0.01

revised.stop = full.table$CP[full.table$nsplit == (num.leaves-1)]

#create a tree to get coordinates for plotting boxplots later
pre.plot = mvpart(FormAllData, data = EpifaunaDataRegression, pretty = T, xv = "min", minauto = F,
                  which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                  xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                  bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                  plot.add = T)

getxy = as.data.frame(cbind(plot(pre.plot, uniform = F)$x, plot(pre.plot, uniform = F)$y))

#create final tree with boxplot leaves
#comment svg and dev.off to show tree in plot window instead of saving to file

final.mrt = mvpart(FormAllData, data = EpifaunaDataRegression, pretty = T, xv = "min", minauto = F,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                   xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                   plot.add = T)

# Taxa Colours
library(ggsci)
pal_ucscgb("default")(24) 

## assign each taxa a color to make it easy to keep leaf plots consistent across trees
names(EpifaunaDataRegression[24:47]) #24 taxa

Amphilepidida ="#FF0000FF"
Amphipoda ="#FF9900FF"
Anthozoa ="#FFCC00FF"
Asteriidae ="#00FF00FF"
Batillariidae="#6699FFFF"
Bryozoa="#CC33FFFF"
Cardiida="#99991EFF"
Cerithiidae="#999999FF"
Chitonida="#FF00CCFF"
Decapod="#CC0000FF"
Dendrasteridae="#FFCCCCFF"
Haminoeidae="#FFFF00FF"
Littorinidae="#CCFF00FF"
Lottiidae="#358000FF"
Myrrhinidae="#0000CCFF"
Mytilidae="#99CCFFFF"
Nemertea="#00FFFFFF"
Neogastropoda="#CCFFFFFF"
Ostreidae="#9900CCFF"
Polychaeta="#CC99FFFF"
Sessilia="#996600FF"
Sphaeromatidae="#666600FF"
Strongylocentrotidae="#666666FF"
Venerida="#CCCCCCFF"

Taxalist= c(Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,Chitonida,Decapod,Dendrasteridae,Haminoeidae,
            Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae,Venerida)   


mvpart(FormAllData, data = EpifaunaDataRegression, pretty = T, xv = "min", minauto = F,
       which = 4, bord = T, uniform = F, text.add = T, branch = 1,
       xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
       bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
       plot.add = T)

legend(x=0.15, y=1,legend = colnames(EpifaunaDataRegression[,24:47]),
       fill = c(Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,Chitonida,Decapod,Dendrasteridae,Haminoeidae,
                Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae,Venerida),
       bty = "n",cex = 0.75)

#add column to df that records which leave each portfolio belongs in
leaf.assign = cbind(as.matrix(final.mrt$where), EpifaunaDataRegression[24:47])
colnames(leaf.assign) = c("Leaf", colnames(EpifaunaDataRegression[24:47]))
leaf.assign$Leaf
leaf.assign$Leaf=as.factor(leaf.assign$Leaf)
#specify coordinates at which leaves will be plotted
leaves = as.numeric(levels(as.factor(pre.plot$where)))
names(leaf.assign)
#add boxplots to terminal nodes

# test=  (log10(data.matrix(EpifaunaDataRegression[,24:47])+1))

for (i in 1:length(leaves)){
  subplot(boxplot.matrix(as.matrix(  log10(EpifaunaDataRegression[which(leaf.assign$Leaf==leaves[i]),24:47]+1)), 
                         col = c(Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,Chitonida,Decapod,Dendrasteridae,Haminoeidae,
                                 Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae,Venerida),
                         xaxt = "n", yaxt = "n",
                         cex = .09, labels = T, bg = "white"),
          x = 0.98*getxy$V1[leaves[i]], 
          y = 0.9738*getxy$V2[leaves[i]],
          size = c(1.6, .9))}

## Plotting

### Exporting as big and bold as possible, resize later
### Epifauna All data Done
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data /Constructing Figure")

jpeg(filename = "MRTEpifaunaAllDataFinal11.jpeg",width = 40, height = 20, units = "cm", pointsize = 9, quality = 100,res=800)

final.mrt = mvpart(FormAllData, data = EpifaunaDataRegression, pretty = T, xv = "min", minauto = T,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,digits = 3,
                   xadj = 0.7, yadj = 1.2, use.n = T, margin =0.072, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                   plot.add = T)

for (i in 1:length(leaves)){
  subplot(boxplot.matrix(as.matrix(  log10(EpifaunaDataRegression[which(leaf.assign$Leaf==leaves[i]),24:47]+1)), 
                         col = c(Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,Chitonida,Decapod,Dendrasteridae,Haminoeidae,
                                 Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae,Venerida),
                         xaxt = "n", yaxt = "n",
                         cex = .09, labels = T, bg = "white"),
          x = 0.98*getxy$V1[leaves[i]], 
          y = 0.9738*getxy$V2[leaves[i]],
          size = c(1.8, .9))}

legend(x=.44, y=1.025,legend = colnames(EpifaunaDataRegression[,24:47]),
       fill = c(Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,Chitonida,Decapod,Dendrasteridae,Haminoeidae,
                Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae,Venerida),
       bty = "n",cex = 1.2)
for (i in 1:length(leaves)){
  text(0.98*getxy$V1[leaves[i]],0.945*getxy$V2[leaves[i]], paste("Leaf",i))}
dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")


### Plotting Legend
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data /Constructing Figure")

jpeg(filename = "MRTEpifaunaLegend1.jpeg",width = 40, height = 20, units = "cm", pointsize = 9, quality = 100,res=800)

final.mrt = mvpart(FormAllData, data = EpifaunaDataRegression, pretty = T, xv = "min", minauto = T,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,digits = 3,
                   xadj = 0.7, yadj = 1.2, use.n = T, margin =1, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                   plot.add = T)

legend(x=-4, y=1.23,legend = colnames(EpifaunaDataRegression[,24:47]),
       fill = c(Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,Chitonida,Decapod,Dendrasteridae,Haminoeidae,
                Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae,Venerida),
       bty = "n",cex = 1.8)

dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")




#### Find discriminant species with MRT results ####
library(party)
library(partykit)
library(vegan)
library(rdaTest)
library(labdsv)
library(plyr)
library(MASS)
library(mvpart)



final.mrt = mvpart(FormAllData, data = EpifaunaDataRegression, pretty = T, xv = "min", minauto = T,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,digits = 3,
                   xadj = 0.7, yadj = 1.2, use.n = T, margin =0.072, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                   plot.add = T)

TaxaData=EpifaunaDataRegression[,24:47]

doubs.mrt <- mvpart((log10(data.matrix(EpifaunaDataRegression[,24:47])+1))~ (Zone+Cobble+Gravel+Sand+Mud.silt+Bivalve.Shells+Oyster.Shell+Sanddollar.dead+Gracilaria+Brown.Algae+Eelgrass+Ulva+Sargassum+Fuchus+Mastocarpus),
                    data = EpifaunaDataRegression, legend=FALSE, margin=0.01, cp=0.01177389, xv="min", 
                    which=4)

doubs.mrt.wrap<-MRT(doubs.mrt, percent=10, species=colnames(TaxaData))
summary(doubs.mrt.wrap)

### Make above node information into a table

# Extract indval p-values
doubs.mrt.indval<-indval(TaxaData,doubs.mrt$where)
doubs.mrt.indval$pval

# Extract indicator species of each node, with its indval
doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval<=0.05)]
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval<=0.05)]






#### Garden Data Multivaraite Regression ####
### Workflow development      
names(GardenData2)
names(GardenData2[21:38])

FormGardenData =   (log10(data.matrix(GardenData2[,21:38])+1))    ~(Zone+Cobble+Gravel+Sand+Mud.silt+Bivalve.Shells+Oyster.Shell+Gracilaria+Ulva+Sargassum+Fuchus+Mastocarpus)

stop.rule = 0.001

mrt.error = mvpart(FormGardenData, data = GardenData2, pretty = T, xv = "min", minauto =T,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                   xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = stop.rule, xval = 10), #100-fold cross validation
                   plot.add = F)

plotcp(mrt.error, upper = "size") # size is for # leaves, splits is for # splits
points(mrt.error$cptable[,"rel error"], col = "red")
legend(x = "topright", legend = c("rel error", "xvalidated rel error"), 
       col = c("red", "blue"), pch = c(1, 1))

full.table = as.data.frame(mrt.error$cptable)
print(full.table)

num.leaves = 7
final.cp = 0.01

revised.stop = full.table$CP[full.table$nsplit == (num.leaves-1)]

#create a tree to get coordinates for plotting boxplots later
pre.plot = mvpart(FormGardenData, data = GardenData2, pretty = T, xv = "min", minauto = T,
                  which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                  xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                  bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 100),
                  plot.add = T)

getxy = as.data.frame(cbind(plot(pre.plot, uniform = F)$x, plot(pre.plot, uniform = F)$y))
dev.off()

#create final tree with boxplot leaves
#comment svg and dev.off to show tree in plot window instead of saving to file

final.mrt = mvpart(FormGardenData, data = GardenData2, pretty = T, xv = "min", minauto = T,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                   xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 100),
                   plot.add = T)

# Taxa Colours
GardenData2
names(GardenData2)
names(GardenData2[21:38])

taxalist = c(Amphilepidida,Amphipoda,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,Chitonida,Decapod,Littorinidae,
             Lottiidae,Myrrhinidae,Mytilidae,Neogastropoda,Polychaeta,Sessilia,Sphaeromatidae,Venerida)     

mvpart(FormGardenData, data = GardenData2, pretty = T, xv = "min", minauto = T,
       which = 4, bord = T, uniform = F, text.add = T, branch = 1,
       xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
       bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 100),
       plot.add = T)

legend("topleft", legend = colnames(GardenData2[,21:38]),
       fill = taxalist,
       bty = "n",cex = 0.75)

### legend matches previous colour codes ###

#add column to df that records which leave each portfolio belongs in
leaf.assign = cbind(as.matrix(final.mrt$where), GardenData2[,21:38])
colnames(leaf.assign) = c("Leaf", colnames(GardenData2[,21:38]))
leaf.assign$Leaf
leaf.assign$Leaf=as.factor(leaf.assign$Leaf)
#specify coordinates at which leaves will be plotted
leaves = as.numeric(levels(as.factor(pre.plot$where)))
names(leaf.assign)
#add boxplots to terminal nodes

# test=  (log10(data.matrix(EpifaunaDataRegression[,24:47])+1))

for (i in 1:length(leaves)){
  subplot(boxplot.matrix(as.matrix(  log10(GardenData2[which(leaf.assign$Leaf==leaves[i]),21:38]+1)), 
                         col = taxalist,
                         xaxt = "n", yaxt = "n",
                         cex = .09, labels = T, bg = "white"),
          x = 0.98*getxy$V1[leaves[i]], 
          y = 0.9738*getxy$V2[leaves[i]],
          size = c(1.6, .9))}

## Plotting

### Exporting as big and bold as possible, resize later
### Epifauna All data Done
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data /Constructing Figure")

jpeg(filename = "MRTEpifaunaGardenDataFinal3.jpeg",width = 40, height = 20, units = "cm", pointsize = 9, quality = 100,res=800)

final.mrt = mvpart(FormGardenData, data = GardenData2, pretty = T, xv = "min", minauto = T,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,digits = 3,
                   xadj = 0.7, yadj = 1.2, use.n = T, margin =0.05, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 100),
                   plot.add = T)

for (i in 1:length(leaves)){
  subplot(boxplot.matrix(as.matrix(  log10(GardenData2[which(leaf.assign$Leaf==leaves[i]),21:38]+1)), 
                         col = taxalist,
                         xaxt = "n", yaxt = "n",
                         cex = .09, labels = T, bg = "white"),
          x = 0.98*getxy$V1[leaves[i]], 
          y = 0.943*getxy$V2[leaves[i]],
          size = c(1.8, .9))}

for (i in 1:length(leaves)){
  text(0.98*getxy$V1[leaves[i]],0.85*getxy$V2[leaves[i]], paste("Leaf",i))}
dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")



#### Farm Data Multivaraite Regression ####
names(FarmData2)
names(FarmData2[23:38])

FormFarmData =   (log10  (data.matrix(FarmData2[,23:38])+1)  )    ~(Zone+Cobble+Gravel+Sand+Mud.silt+Bivalve.Shells+Oyster.Shell+Gracilaria+Brown.Algae+Eelgrass+Ulva+Sargassum+Fuchus+Mastocarpus)

stop.rule = 0.001

mrt.error = mvpart(FormFarmData, data = FarmData2, pretty = T, xv = "min", minauto =T,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                   xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = stop.rule, xval = 10), #100-fold cross validation
                   plot.add = F)

plotcp(mrt.error, upper = "size") # size is for # leaves, splits is for # splits
points(mrt.error$cptable[,"rel error"], col = "red")
legend(x = "topright", legend = c("rel error", "xvalidated rel error"), 
       col = c("red", "blue"), pch = c(1, 1))

full.table = as.data.frame(mrt.error$cptable)
print(full.table)

num.leaves = 4
final.cp = 0.01

revised.stop = full.table$CP[full.table$nsplit == (num.leaves-1)]

#create a tree to get coordinates for plotting boxplots later
pre.plot = mvpart(FormFarmData, data = FarmData2, pretty = T, xv = "min", minauto = T,
                  which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                  xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                  bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                  plot.add = T)

getxy = as.data.frame(cbind(plot(pre.plot, uniform = F)$x, plot(pre.plot, uniform = F)$y))
dev.off()

#create final tree with boxplot leaves
#comment svg and dev.off to show tree in plot window instead of saving to file

final.mrt = mvpart(FormFarmData, data = FarmData2, pretty = T, xv = "min", minauto = F,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                   xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                   plot.add = T)

# Taxa Colours
FarmData2
names(FarmData2)
names(FarmData2[23:38])

taxalist = c(Amphipoda,Anthozoa,Asteriidae,Batillariidae,Chitonida,Decapod,Dendrasteridae,Littorinidae,Lottiidae,Mytilidae,
             Neogastropoda,Ostreidae,Polychaeta,Sessilia,Strongylocentrotidae,Venerida)

mvpart(FormFarmData, data = FarmData2, pretty = T, xv = "min", minauto = T,
       which = 4, bord = T, uniform = F, text.add = T, branch = 1,
       xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
       bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
       plot.add = T)

legend("topleft", legend = colnames(FarmData2[,23:38]),
       fill = taxalist,
       bty = "n",cex = 0.75)

### legend matches previous colour codes ###

#add column to df that records which leave each portfolio belongs in
leaf.assign = cbind(as.matrix(final.mrt$where), FarmData2[,23:38])
colnames(leaf.assign) = c("Leaf", colnames(FarmData2[,23:38]))
leaf.assign$Leaf
leaf.assign$Leaf=as.factor(leaf.assign$Leaf)
#specify coordinates at which leaves will be plotted
leaves = as.numeric(levels(as.factor(pre.plot$where)))
names(leaf.assign)
#add boxplots to terminal nodes

for (i in 1:length(leaves)){
  subplot(boxplot.matrix(as.matrix(  log10(FarmData2[which(leaf.assign$Leaf==leaves[i]),23:38]+1)), 
                         col = taxalist,
                         xaxt = "n", yaxt = "n",
                         cex = .09, labels = T, bg = "white"),
          x = 0.98*getxy$V1[leaves[i]], 
          y = 0.9738*getxy$V2[leaves[i]],
          size = c(1.6, .9))}

## Plotting

### Exporting as big and bold as possible, resize later
### Epifauna All data Done
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data /Constructing Figure")

jpeg(filename = "MRTEpifaunaFarmDataFinal8.jpeg",width = 40, height = 20, units = "cm", pointsize = 9, quality = 100,res=800)

final.mrt = mvpart(FormFarmData, data = FarmData2, pretty = T, xv = "min", minauto = T,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,digits = 3,
                   xadj = 0.7, yadj = 1.2, use.n = T, margin =0.09, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                   plot.add = T)

# plot(final.mrt); text(final.mrt)

for (i in 1:length(leaves)){
  subplot(boxplot.matrix(as.matrix(  log10(FarmData2[which(leaf.assign$Leaf==leaves[i]),23:38]+1)), 
                         col = taxalist,
                         xaxt = "n", yaxt = "n",
                         cex = .09, labels = T, bg = "white"),
          x = 0.98*getxy$V1[leaves[i]], 
          y = 0.974*getxy$V2[leaves[i]],
          size = c(1.8, .9))}

for (i in 1:length(leaves)){
  text(0.98*getxy$V1[leaves[i]],0.94*getxy$V2[leaves[i]], paste("Leaf",i))}
dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")



#### Reference Data Multivaraite Regression ####
names(ReferenceData2)
names(ReferenceData2[23:44]) ## 22 taxa

FormReferenceData =   (log10  (data.matrix(ReferenceData2[,23:44])+1)  )    ~(Zone+Cobble+Gravel+Sand+Mud.silt+Bivalve.Shells+Oyster.Shell+Sanddollar.dead+Gracilaria+Brown.Algae+Eelgrass+Ulva+Fuchus+Mastocarpus)

stop.rule = 0.001

mrt.error = mvpart(FormReferenceData, data = ReferenceData2, pretty = T, xv = "min", minauto =T,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                   xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = stop.rule, xval = 10), #100-fold cross validation
                   plot.add = F)

plotcp(mrt.error, upper = "size") # size is for # leaves, splits is for # splits
points(mrt.error$cptable[,"rel error"], col = "red")
legend(x = "topright", legend = c("rel error", "xvalidated rel error"), 
       col = c("red", "blue"), pch = c(1, 1))

full.table = as.data.frame(mrt.error$cptable)
print(full.table)

num.leaves = 5
final.cp = 0.01

revised.stop = full.table$CP[full.table$nsplit == (num.leaves-1)]

#create a tree to get coordinates for plotting boxplots later
pre.plot = mvpart(FormReferenceData, data = ReferenceData2, pretty = T, xv = "min", minauto = T,
                  which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                  xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                  bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                  plot.add = T)

getxy = as.data.frame(cbind(plot(pre.plot, uniform = F)$x, plot(pre.plot, uniform = F)$y))

#create final tree with boxplot leaves
#comment svg and dev.off to show tree in plot window instead of saving to file

final.mrt = mvpart(FormReferenceData, data = ReferenceData2, pretty = T, xv = "min", minauto = F,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,
                   xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                   plot.add = T)

# Taxa Colours
ReferenceData2
names(ReferenceData2)
names(ReferenceData2[23:44])

taxalist = c(Amphipoda,Anthozoa,Asteriidae,Batillariidae,Cardiida,Cerithiidae,Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,           
             Myrrhinidae,Mytilidae,Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae,Venerida)

mvpart(FormReferenceData, data = ReferenceData2, pretty = T, xv = "min", minauto = T,
       which = 4, bord = T, uniform = F, text.add = T, branch = 1,
       xadj = .7, yadj = 1.2, use.n = T, margin = 0.05, keep.y = F,
       bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
       plot.add = T)

legend("topleft", legend = colnames(ReferenceData2[,23:44]),
       fill = taxalist,
       bty = "n",cex = 0.75)

### legend matches previous colour codes ###

#add column to df that records which leave each portfolio belongs in
leaf.assign = cbind(as.matrix(final.mrt$where), ReferenceData2[,23:44])
colnames(leaf.assign) = c("Leaf", colnames(ReferenceData2[,23:44]))
leaf.assign$Leaf
leaf.assign$Leaf=as.factor(leaf.assign$Leaf)
#specify coordinates at which leaves will be plotted
leaves = as.numeric(levels(as.factor(pre.plot$where)))
names(leaf.assign)
#add boxplots to terminal nodes

for (i in 1:length(leaves)){
  subplot(boxplot.matrix(as.matrix(  log10(ReferenceData2[which(leaf.assign$Leaf==leaves[i]),23:44]+1)), 
                         col = taxalist,
                         xaxt = "n", yaxt = "n",
                         cex = .09, labels = T, bg = "white"),
          x = 0.98*getxy$V1[leaves[i]], 
          y = 0.9738*getxy$V2[leaves[i]],
          size = c(1.6, .9))}

## Plotting

### Exporting as big and bold as possible, resize later
### Epifauna All data Done
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data /Constructing Figure")

jpeg(filename = "MRTEpifaunaReferenceDataFinal1.jpeg",width = 40, height = 20, units = "cm", pointsize = 9, quality = 100,res=800)

final.mrt = mvpart(FormReferenceData, data = ReferenceData2, pretty = T, xv = "min", minauto = T,
                   which = 4, bord = T, uniform = F, text.add = T, branch = 1,digits = 3,
                   xadj = 0.7, yadj = 1.2, use.n = T, margin =0.09, keep.y = F,
                   bars = F, all.leaves = F, control = rpart.control(cp = revised.stop, xval = 10),
                   plot.add = T)

# plot(final.mrt); text(final.mrt)

for (i in 1:length(leaves)){
  subplot(boxplot.matrix(as.matrix(  log10(ReferenceData2[which(leaf.assign$Leaf==leaves[i]),23:44]+1)), 
                         col = taxalist,
                         xaxt = "n", yaxt = "n",
                         cex = .09, labels = T, bg = "white"),
          x = 0.98*getxy$V1[leaves[i]], 
          y = 0.96*getxy$V2[leaves[i]],
          size = c(1.8, .9))}

for (i in 1:length(leaves)){
  text(0.98*getxy$V1[leaves[i]],0.92*getxy$V2[leaves[i]], paste("Leaf",i))}
dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")





#### Figure 4 B: Multivaraite Random Forest ####
library(randomForestSRC)
EpifaunaDataRegression
EpifaunaDataRegression$Year=as.factor(EpifaunaDataRegression$Year)
EpifaunaDataRegression2


names(EpifaunaDataRegression)
TaxaData= (log10((EpifaunaDataRegression[,24:47])+1))  
names(EpifaunaDataRegression[,9:22])
names(EpifaunaDataRegression[,1:3])
names(EpifaunaDataRegression[5])
names(EpifaunaDataRegression[8])
TaxaData3=cbind(TaxaData,(EpifaunaDataRegression[,9:22]), (EpifaunaDataRegression[,1:3]),(EpifaunaDataRegression[5]),(EpifaunaDataRegression[8]))

MRT.FormAllData2 <- rfsrc(Multivar(Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,
                                   Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,
                                   Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae,
                                   Venerida) ~SiteType2+Area+Site+Zone+Cobble+Gravel+Bivalve.Shells+Mastocarpus+Ulva+Oyster.Shell+Mud.silt+Sanddollar.dead+Gracilaria, 
                          data = TaxaData3,ntree=6000,block.size=1,nodedepth=6, importance = TRUE) 

# Forest terminal node size     =  nodesize
# ntree=3000,nodesize=20,nodedepth=3,block.size=1,nsplit=3,
## switch to 10,000 trees

err2 <- get.mv.error(MRT.FormAllData2)
vmp2 <- get.mv.vimp(MRT.FormAllData2)
pred2 <- get.mv.predicted(MRT.FormAllData2)
err.std2 <- get.mv.error(MRT.FormAllData2, standardize = TRUE)
vmp.std2 <- get.mv.vimp(MRT.FormAllData2, standardize = TRUE)

plot(MRT.FormAllData2,m.target="Amphilepidida")
plot(MRT.FormAllData2,m.target="Amphipoda")
plot(MRT.FormAllData2,m.target="Anthozoa")
plot(MRT.FormAllData2,m.target="Asteriidae")
plot(MRT.FormAllData2,m.target="Batillariidae")
plot(MRT.FormAllData2,m.target="Bryozoa")
plot(MRT.FormAllData2,m.target="Cardiida")
plot(MRT.FormAllData2,m.target="Cerithiidae")
plot(MRT.FormAllData2,m.target="Chitonida")
plot(MRT.FormAllData2,m.target="Decapod")
plot(MRT.FormAllData2,m.target="Dendrasteridae")
plot(MRT.FormAllData2,m.target="Haminoeidae")
plot(MRT.FormAllData2,m.target="Littorinidae")
plot(MRT.FormAllData2,m.target="Lottiidae")
plot(MRT.FormAllData2,m.target="Myrrhinidae")
plot(MRT.FormAllData2,m.target="Mytilidae")
plot(MRT.FormAllData2,m.target="Nemertea")
plot(MRT.FormAllData2,m.target="Neogastropoda")
plot(MRT.FormAllData2,m.target="Ostreidae")
plot(MRT.FormAllData2,m.target="Polychaeta")
plot(MRT.FormAllData2,m.target="Sessilia")
plot(MRT.FormAllData2,m.target="Sphaeromatidae")
plot(MRT.FormAllData2,m.target="Strongylocentrotidae")
plot(MRT.FormAllData2,m.target="Venerida")


### Export and Plot
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data /Constructing Figure")

jpeg(filename = "ErrorRateAmphipoda.jpeg",width = 30, height = 15, units = "cm", pointsize = 9, quality = 100,res=600)
plot(MRT.FormAllData2,m.target="Amphipoda")
dev.off()

jpeg(filename = "ErrorRateAsteriidae.jpeg",width = 30, height = 15, units = "cm", pointsize = 9, quality = 100,res=600)
plot(MRT.FormAllData2,m.target="Asteriidae")
dev.off()

jpeg(filename = "ErrorRateBatillariidae.jpeg",width = 30, height = 15, units = "cm", pointsize = 9, quality = 100,res=600)
plot(MRT.FormAllData2,m.target="Batillariidae")
dev.off()

jpeg(filename = "ErrorRateDecapod.jpeg",width = 30, height = 15, units = "cm", pointsize = 9, quality = 100,res=600)
plot(MRT.FormAllData2,m.target="Decapod")
dev.off()

jpeg(filename = "ErrorRateLottiidae.jpeg",width = 30, height = 15, units = "cm", pointsize = 9, quality = 100,res=600)
plot(MRT.FormAllData2,m.target="Lottiidae")
dev.off()

jpeg(filename = "ErrorRateSessilia.jpeg",width = 30, height = 15, units = "cm", pointsize = 9, quality = 100,res=600)
plot(MRT.FormAllData2,m.target="Sessilia")
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")


#### Extract Each Explained Variance and Test Set Error Rate ####

#### Plotting VariableImportance ####

vmp2 <- get.mv.vimp(MRT.FormAllData2)
# write.csv(vmp2,"VariableImportanceAverageFinal.csv")
# write.csv(vmp.std2,"VariableImportanceAverageStandardFinal.csv")

VariableImportance=read.csv("VariableImportanceAverageStandardFinalV2.csv",header=T)
VariableImportance
names(VariableImportance)
head(VariableImportance)
VariableImportance$Variable
str(VariableImportance)
names(VariableImportance)


### species plots ###
size3=12

# Amphilepidida
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Amphilepidida)])
Amphilepidida.Plot=ggplot(VariableImportance, aes(x = Variable, y = Amphilepidida,fill=AmphilepididaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Amphilepidida.Plot

# Amphipoda
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Amphipoda)])
Amphipoda.Plot=ggplot(VariableImportance, aes(x = Variable, y = Amphipoda,fill=AmphipodaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Amphipoda.Plot

# Anthozoa
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Anthozoa)])
Anthozoa.Plot=ggplot(VariableImportance, aes(x = Variable, y = Anthozoa,fill=AnthozoaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Anthozoa.Plot

# Asteriidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Asteriidae)])
Asteriidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Asteriidae,fill=AsteriidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Asteriidae.Plot

# Batillariidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Batillariidae)])
Batillariidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Batillariidae,fill=BatillariidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3","red2"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Batillariidae.Plot

# Bryozoa
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Bryozoa)])
Bryozoa.Plot=ggplot(VariableImportance, aes(x = Variable, y = Bryozoa,fill=BryozoaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Bryozoa.Plot

# Cardiida
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Cardiida)])
Cardiida.Plot=ggplot(VariableImportance, aes(x = Variable, y = Cardiida,fill=CardiidaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Cardiida.Plot

# Cerithiidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Cerithiidae)])
Cerithiidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Cerithiidae,fill=CerithiidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Cerithiidae.Plot

# Chitonida
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Chitonida)])
Chitonida.Plot=ggplot(VariableImportance, aes(x = Variable, y = Chitonida,fill=ChitonidaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Chitonida.Plot

# Decapod
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Decapod)])
Decapod.Plot=ggplot(VariableImportance, aes(x = Variable, y = Decapod,fill=DecapodDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3","red2"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Decapod.Plot

# Dendrasteridae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Dendrasteridae)])
Dendrasteridae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Dendrasteridae,fill=DendrasteridaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Dendrasteridae.Plot

# Haminoeidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Haminoeidae)])
Haminoeidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Haminoeidae,fill=HaminoeidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Haminoeidae.Plot

# Littorinidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Littorinidae)])
Littorinidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Littorinidae,fill=LittorinidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3","red2"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Littorinidae.Plot

# Lottiidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Lottiidae)])
Lottiidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Lottiidae,fill=LottiidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3","red2"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Lottiidae.Plot

# Myrrhinidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Myrrhinidae)])
Myrrhinidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Myrrhinidae,fill=MyrrhinidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Myrrhinidae.Plot

# Mytilidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Mytilidae)])
Mytilidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Mytilidae,fill=MytilidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Mytilidae.Plot

# Nemertea
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Nemertea)])
Nemertea.Plot=ggplot(VariableImportance, aes(x = Variable, y = Nemertea,fill=NemerteaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Nemertea.Plot

# Neogastropoda
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Neogastropoda)])
Neogastropoda.Plot=ggplot(VariableImportance, aes(x = Variable, y = Neogastropoda,fill=NeogastropodaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Neogastropoda.Plot

# Ostreidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Ostreidae)])
Ostreidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Ostreidae,fill=OstreidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Ostreidae.Plot

# Polychaeta
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Polychaeta)])
Polychaeta.Plot=ggplot(VariableImportance, aes(x = Variable, y = Polychaeta,fill=PolychaetaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3","red2"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Polychaeta.Plot

# Sessilia
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Sessilia)])
Sessilia.Plot=ggplot(VariableImportance, aes(x = Variable, y = Sessilia,fill=SessiliaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3","red2"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Sessilia.Plot

# Sphaeromatidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Sphaeromatidae)])
Sphaeromatidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Sphaeromatidae,fill=SphaeromatidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Sphaeromatidae.Plot

# Strongylocentrotidae
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Strongylocentrotidae)])
Strongylocentrotidae.Plot=ggplot(VariableImportance, aes(x = Variable, y = Strongylocentrotidae,fill=StrongylocentrotidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Strongylocentrotidae.Plot

# Venerida
VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$Venerida)])
Venerida.Plot=ggplot(VariableImportance, aes(x = Variable, y = Venerida,fill=VeneridaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Venerida.Plot



## Averaged
VariableImportance=read.csv("VariableImportanceAverageStandardFinalV2.csv",header=T)

VariableImportance
names(VariableImportance)
## Standardized ##

VariableImportance$Variable <- factor(VariableImportance$Variable, levels = VariableImportance$Variable[order(VariableImportance$AveragePercent)])
AveragePercent.Plot=ggplot(VariableImportance, aes(x = Variable, y = AveragePercent,fill=Directionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=3)+ylab("Variable Importance (%)")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
AveragePercent.Plot

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data /Constructing Figure")
jpeg(filename = "VariableImportanceAveragedStandardized5.jpeg",width = 28, height = 25, units = "cm", pointsize = 9, quality = 100,res=600)
AveragePercent.Plot
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")

###
Amphilepidida.Plot
Amphipoda.Plot
Anthozoa.Plot
Asteriidae.Plot
Batillariidae.Plot
Bryozoa.Plot
Cardiida.Plot
Cerithiidae.Plot
Chitonida.Plot
Decapod.Plot
Dendrasteridae.Plot
Haminoeidae.Plot
Littorinidae.Plot
Lottiidae.Plot
Myrrhinidae.Plot
Mytilidae.Plot
Nemertea.Plot
Neogastropoda.Plot
Ostreidae.Plot
Polychaeta.Plot
Sessilia.Plot
Sphaeromatidae.Plot
Strongylocentrotidae.Plot
Venerida.Plot



library(cowplot)

## Supplemental plot ##
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data /Constructing Figure")
jpeg(filename = "VariableImportance6.jpeg",width = 50, height = 80, units = "cm", pointsize = 9, quality = 100,res=800)
plot_grid(Amphipoda.Plot,        Amphilepidida.Plot,        Chitonida.Plot, 
          Lottiidae.Plot,        Batillariidae.Plot,         Strongylocentrotidae.Plot,   
          Cardiida.Plot,         Mytilidae.Plot,             Neogastropoda.Plot,
          Decapod.Plot,          Littorinidae.Plot,         Cerithiidae.Plot,
          Polychaeta.Plot,      Venerida.Plot,            Asteriidae.Plot,
          Ostreidae.Plot,        Bryozoa.Plot,           Dendrasteridae.Plot,
          Sphaeromatidae.Plot,   Myrrhinidae.Plot,              Haminoeidae.Plot,
          Sessilia.Plot,         Anthozoa.Plot,             Nemertea.Plot,
          labels = c('A) Amphipoda',        'B) Amphilepidida',       'C) Chitonida',
                     'D) Lottiidae',       'E) Batillariidae',        'F) Strongylocentrotidae',
                     'G) Cardiida',         'H) Mytilidae',            'I) Neogastropoda',
                     'J) Decapod',          'K) Littorinidae',        'L) Cerithiidae',
                     'M) Polychaeta',      'N) Venerida',           'O) Asteriidae',
                     'P) Ostreidae',        'Q) Bryozoa',          'R) Dendrasteridae',
                     'S) Sphaeromatidae',   'T) Myrrhinidae',        'U) Haminoeidae',
                     'V) Sessilia',         'W) Anthozoa',            'X) Nemertea'),nrow = 8, align="hv")
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")





MRT.FormAllData2
summary(MRT.FormAllData2)



#### Varability Explained and Random Forest Predictions ####

## Error
MRT.FormAllData2 <- rfsrc(Multivar(Bryozoa,Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Cardiida,Cerithiidae,
                                   Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,
                                   Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae,
                                   Venerida) ~SiteType2+Area+Site+Zone+Cobble+Gravel+Bivalve.Shells+Mastocarpus+Ulva+Oyster.Shell+Mud.silt+Sanddollar.dead+Gracilaria, 
                          data = TaxaData3,ntree=6000,block.size=1,nodedepth=6, importance = TRUE) 
MRT.FormAllData2

err2 <- get.mv.error(MRT.FormAllData2)
vmp2 <- get.mv.vimp(MRT.FormAllData2)
pred2 <- get.mv.predicted(MRT.FormAllData2)
err.std2 <- get.mv.error(MRT.FormAllData2, standardize = TRUE)
vmp.std2 <- get.mv.vimp(MRT.FormAllData2, standardize = TRUE)

### Explained Variance
MRT.FormAllData2 <- rfsrc(Multivar(Venerida,Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,
                                   Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,
                                   Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae) ~SiteType2+Area+Site+Zone+Cobble+Gravel+Bivalve.Shells+Mastocarpus+Ulva+Oyster.Shell+Mud.silt+Sanddollar.dead+Gracilaria, 
                          data = TaxaData3,ntree=6000,block.size=1,nodedepth=6, importance = TRUE) 
MRT.FormAllData2

### Prediction
# train <- sample(1:nrow(TaxaData3), round(nrow(TaxaData3) * 0.80))
MRT.FormAllDataTrain <- rfsrc(Multivar(Venerida,Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,
                                       Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,
                                       Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae) ~SiteType2+Area+Site+Zone+Cobble+Gravel+Bivalve.Shells+Mastocarpus+Ulva+Oyster.Shell+Mud.silt+Sanddollar.dead+Gracilaria, 
                              data = TaxaData3[train, ],ntree=6000,block.size=1,nodedepth=6, importance = TRUE) 

MRT.FormAllDataPredict<- predict(MRT.FormAllDataTrain, TaxaData3[-train , ])
print(MRT.FormAllDataTrain)
print(MRT.FormAllDataPredict)


### Supplemental figure: Conceptual Diagram 
# made in external program


#### Supplemental: Site Complexity models ####
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos")
library(magick)
library(ggplot2)
library(colordistance)
library(cowplot)
library(ggpubr)

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final")
list.files()

GardenColor=("#7CAE00")   #### Green 
FarmColor=("#0d7dd9")     #### Blue 
ReferenceGardenColor=("#fa6920") #### Orange
ReferenceFarmColor=("#ff1100") #### Red

BoardSize=2


### Baynes Sound SF
# Reid
ReidSF <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Baynes-SF-Reid 2.jpeg")
ReidSF
ReidSF <- image_ggplot(ReidSF)
ReidSF

ReidSF=ReidSF+theme(panel.border = element_rect(color = "#0d7dd9",
                                          fill = NA,
                                          size = BoardSize))
ReidSF

# Taylor
TaylorSF <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Baynes-SF-Taylor.jpeg")
TaylorSF
TaylorSF <- image_ggplot(TaylorSF)
TaylorSF

TaylorSF=TaylorSF+theme(panel.border = element_rect(color = "#0d7dd9",
                                                fill = NA,
                                                size = BoardSize))
TaylorSF

# TRAN
TRANSF <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Baynes-SF-TRAN 2.jpeg")
TRANSF
TRANSF <- image_ggplot(TRANSF)
TRANSF

TRANSF=TRANSF+theme(panel.border = element_rect(color = "#0d7dd9",
                                                    fill = NA,
                                                    size = BoardSize))
TRANSF


### Baynes Sound SF Ref

# DBR
DBRSFref <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Baynes-SFref-DBR.jpeg")
DBRSFref
DBRSFref <- image_ggplot(DBRSFref)
DBRSFref

DBRSFref=DBRSFref+theme(panel.border = element_rect(color = "#ff1100",
                                                fill = NA,
                                                size = BoardSize))
DBRSFref

# NOT
NOTSFref <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Baynes-SFref-NOT.jpeg")
NOTSFref
NOTSFref <- image_ggplot(NOTSFref)
NOTSFref

NOTSFref=NOTSFref+theme(panel.border = element_rect(color = "#ff1100",
                                              fill = NA,
                                              size = BoardSize))
NOTSFref


# RES
RESSFref <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Baynes-SFref-RES 2.jpeg")
RESSFref
RESSFref <- image_ggplot(RESSFref)
RESSFref

RESSFref=RESSFref+theme(panel.border = element_rect(color = "#ff1100",
                                                    fill = NA,
                                                    size = BoardSize))
RESSFref


### Calvert CG

# CGN
CGNcg <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Calvert-CG-CGN.jpeg")
CGNcg
CGNcg <- image_ggplot(CGNcg)
CGNcg

CGNcg=CGNcg+theme(panel.border = element_rect(color = "#7CAE00",
                                                    fill = NA,
                                                    size = BoardSize))
CGNcg

# NUXI
NUXIcg <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Calvert-CG-NUXI.jpeg")
NUXIcg
NUXIcg <- image_ggplot(NUXIcg)
NUXIcg

NUXIcg=NUXIcg+theme(panel.border = element_rect(color = "#7CAE00",
                                              fill = NA,
                                              size = BoardSize))
NUXIcg


# Twin
TWINcg <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Calvert-CG-Twin.jpeg")
TWINcg
TWINcg <- image_ggplot(TWINcg)
TWINcg

TWINcg=TWINcg+theme(panel.border = element_rect(color = "#7CAE00",
                                                fill = NA,
                                                size = BoardSize))
TWINcg


### Calvert CG Ref
# HCGN
HCGNcgREF <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Calvert-CGref-HCGN.jpeg")
HCGNcgREF
HCGNcgREF <- image_ggplot(HCGNcgREF)
HCGNcgREF

HCGNcgREF=HCGNcgREF+theme(panel.border = element_rect(color = "#fa6920",
                                                fill = NA,
                                                size = BoardSize))
HCGNcgREF


# Martin
MartincgREF <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Calvert-CGref-Martin 2.jpeg")
MartincgREF
MartincgREF <- image_ggplot(MartincgREF)
MartincgREF

MartincgREF=MartincgREF+theme(panel.border = element_rect(color = "#fa6920",
                                                      fill = NA,
                                                      size = BoardSize))
MartincgREF


# Piles
PilesCGREF <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Calvert-CGref-PBR.jpeg")
PilesCGREF
PilesCGREF <- image_ggplot(PilesCGREF)
PilesCGREF

PilesCGREF=PilesCGREF+theme(panel.border = element_rect(color = "#fa6920",
                                                          fill = NA,
                                                          size = BoardSize))
PilesCGREF




### Quadra CG
# DoubleWallCG
DoubleWallCG <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Qaudra-CG-DoubleWall.jpeg")
DoubleWallCG
DoubleWallCG <- image_ggplot(DoubleWallCG)
DoubleWallCG

DoubleWallCG=DoubleWallCG+theme(panel.border = element_rect(color = "#7CAE00",
                                                        fill = NA,
                                                        size = BoardSize))
DoubleWallCG


# Garden
GardenCG <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Quadra-CG-Garden.jpeg")
GardenCG
GardenCG <- image_ggplot(GardenCG)
GardenCG

GardenCG=GardenCG+theme(panel.border = element_rect(color = "#7CAE00",
                                                            fill = NA,
                                                            size = BoardSize))
GardenCG


# KB11
KB11CG <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Quadra-CG-KB11.jpeg")
KB11CG
KB11CG <- image_ggplot(KB11CG)
KB11CG

KB11CG=KB11CG+theme(panel.border = element_rect(color = "#7CAE00",
                                                    fill = NA,
                                                    size = BoardSize))
KB11CG


### Quadra CG Ref
# KB08
KB08cgRef <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Quadra-CGref-KB08.jpeg")
KB08cgRef
KB08cgRef <- image_ggplot(KB08cgRef)
KB08cgRef

KB08cgRef=KB08cgRef+theme(panel.border = element_rect(color = "#fa6920",
                                                fill = NA,
                                                size = BoardSize))
KB08cgRef



# KB11ref
KB11cgRef <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Quadra-CGref-KB11.jpeg")
KB11cgRef
KB11cgRef <- image_ggplot(KB11cgRef)
KB11cgRef

KB11cgRef=KB11cgRef+theme(panel.border = element_rect(color = "#fa6920",
                                                      fill = NA,
                                                      size = BoardSize))
KB11cgRef


# SWC
SWCcgRef <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Quadra-CGref-SWC 2.jpeg")
SWCcgRef
SWCcgRef <- image_ggplot(SWCcgRef)
SWCcgRef

SWCcgRef=SWCcgRef+theme(panel.border = element_rect(color = "#fa6920",
                                                      fill = NA,
                                                      size = BoardSize))
SWCcgRef


### Quadra SF
# SMB
SMBsf <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Qaudra-SF-SMB 2.jpeg")
SMBsf
SMBsf <- image_ggplot(SMBsf)
SMBsf

SMBsf=SMBsf+theme(panel.border = element_rect(color = "#0d7dd9",
                                                    fill = NA,
                                                    size = BoardSize))
SMBsf


# FAB
FABsf <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Qaudra-SF-FAB.jpeg")
FABsf
FABsf <- image_ggplot(FABsf)
FABsf

FABsf=FABsf+theme(panel.border = element_rect(color = "#0d7dd9",
                                              fill = NA,
                                              size = BoardSize))
FABsf


# LHB
LHBsf <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Qaudra-SF-LHB 2.jpeg")
LHBsf
LHBsf <- image_ggplot(LHBsf)
LHBsf

LHBsf=LHBsf+theme(panel.border = element_rect(color = "#0d7dd9",
                                              fill = NA,
                                              size = BoardSize))
LHBsf


### Quadra SF Ref

# HBR
HBRsfREF <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Quadra-SFref-HBR.jpeg")
HBRsfREF
HBRsfREF <- image_ggplot(HBRsfREF)
HBRsfREF

HBRsfREF=HBRsfREF+theme(panel.border = element_rect(color = "#ff1100",
                                              fill = NA,
                                              size = BoardSize))
HBRsfREF


# LEN
LENsfREF <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Quadra-SFref-LEN.jpeg")
LENsfREF
LENsfREF <- image_ggplot(LENsfREF)
LENsfREF

LENsfREF=LENsfREF+theme(panel.border = element_rect(color = "#ff1100",
                                                    fill = NA,
                                                    size = BoardSize))
LENsfREF


# WWK
WWKsfREF <-image_read("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Complexity Data/Project Photos/Model Photos/Approved Images/Final/Quadra-SFref-WWK 2.jpeg")
WWKsfREF
WWKsfREF <- image_ggplot(WWKsfREF)
WWKsfREF

WWKsfREF=WWKsfREF+theme(panel.border = element_rect(color = "#ff1100",
                                                    fill = NA,
                                                    size = BoardSize))
WWKsfREF


## Baynes
ReidSF       #SF
TaylorSF     #SF
TRANSF       #SF
DBRSFref     #SFref
NOTSFref     #SFref
RESSFref     #SFref
## Calvert
CGNcg        #CG
NUXIcg       #CG
TWINcg       #CG
MartincgREF  #CGref
PilesCGREF   #CGref
HCGNcgREF    #CGref
## Quadra 
DoubleWallCG #CG
GardenCG     #CG
KB11CG       #CG
KB08cgRef    #CGref
KB11cgRef    #CGref
SWCcgRef     #CGref
SMBsf        #SF
FABsf        #SF
LHBsf        #SF
HBRsfREF     #SFref
LENsfREF     #SFref
WWKsfREF     #SFref


## Arrangement is Garden, Garden Reference, Farm, Farm Reference
##       Calvert above Quadra, Quadra above Baynes?

plot_grid(
  CGNcg,        DoubleWallCG,   SMBsf,      TRANSF,
  NUXIcg,       GardenCG,       FABsf,      TaylorSF,
  TWINcg,       KB11CG,         LHBsf,      ReidSF,
  MartincgREF,  KB08cgRef,      WWKsfREF,   DBRSFref,
  PilesCGREF,   KB11cgRef,      HBRsfREF,   NOTSFref, 
  HCGNcgREF,    SWCcgRef,       LENsfREF,   RESSFref, 
  ncol = 4,
  label_size = 12,
  align = "hv")

### Export
## Option 1
tiff(filename = "SupFigSiteComplexModelsV4.tiff",
     width = 30, height = 30, units = "cm", pointsize = 20,
     compression = c("none"),
     bg = "white", res = 100,
     type = c("quartz"))

plot_grid(
  CGNcg,        DoubleWallCG,   SMBsf,      TRANSF,
  NUXIcg,       GardenCG,       FABsf,      TaylorSF,
  TWINcg,       KB11CG,         LHBsf,      ReidSF,
  MartincgREF,  KB08cgRef,      WWKsfREF,   DBRSFref,
  PilesCGREF,   KB11cgRef,      HBRsfREF,   NOTSFref, 
  HCGNcgREF,    SWCcgRef,       LENsfREF,   RESSFref, 
  ncol = 4,
  labels="AUTO",
  label_size = 12,
  align = "hv")

dev.off()

## Option 2
tiff(filename = "SupFigSiteComplexModelsV8.tiff",
     width = 23, height = 30, units = "cm", pointsize = 20,
     compression = c("none"),
     bg = "white", res = 300,
     type = c("quartz"))
ggarrange(
  CGNcg,        DoubleWallCG,   SMBsf,      TRANSF,
  NUXIcg,       GardenCG,       FABsf,      TaylorSF,
  TWINcg,       KB11CG,         LHBsf,      ReidSF,
  MartincgREF,  KB08cgRef,      WWKsfREF,   DBRSFref,
  PilesCGREF,   KB11cgRef,      HBRsfREF,   NOTSFref, 
  HCGNcgREF,    SWCcgRef,       LENsfREF,   RESSFref,  
  ncol = 4, widths = 1.1,heights = 1.1,align = c("hv"),labels="AUTO",
  nrow = 6)
dev.off()
     




#### Supplemental Complexity Data  ####

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")

ComplexityData=read.csv("SfandCGSiteComplexity.csv", header=T)
ComplexityData2=read.csv("SfandCGSiteComplexity2.csv", header=T)
ComplexityData3=read.csv("SfandCGSiteComplexity3.csv", header=T)

ComplexityData
names(ComplexityData)
head(ComplexityData)
str(ComplexityData)

ComplexityData2
names(ComplexityData2)
head(ComplexityData2)
str(ComplexityData2)

###Objectives
## All plots are raw data and logged data
## Plot all data
## Summarize by the four regions/types
## Summarize by types

## Color
ggplot(ComplexityData2, aes(x= RugoistyType, y= Complexity, fill=SiteType) )+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, position = position_dodge(width = 1.0))+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) 

levels(ComplexityData2$RugoistyType)

AverageRugoisty=(subset(ComplexityData2, RugoistyType %in% c("AverageRugoisty")))

ggplot(AverageRugoisty, aes(x= SiteType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, position = position_dodge(width = 1.0))+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")

ggplot(AverageRugoisty, aes(x= SiteType, y= LogComplexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, position = position_dodge(width = 1.0))+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")


Rugoisty=(subset(ComplexityData2, RugoistyType %in% c("Rugoisty1" ,"Rugoisty2","Rugoisty3","Rugoisty4","Rugoisty5")))

ggplot(Rugoisty, aes(x= RugoistyType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, position = position_dodge(width = 1.0))+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")

names(Rugoisty)
ggplot(Rugoisty, aes(x= RugoistyType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, position = position_dodge(width = 1.0))+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")+facet_grid(Region~.)


ggplot(ComplexityData3, aes(x= RugoistyType, y= log10(Complexity), fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, position = position_dodge(width = 1.0))+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")+facet_grid(Region~.)


## Gotta make rug and FD seperately

## FD
names(ComplexityData3)
levels(ComplexityData3$RugoistyType)

FDaverages=(subset(ComplexityData3, RugoistyType %in% c("FractalD50Average" ,"FractalD25Average","FractalD5Average","FractalD2.5Average")))
AverageRugoisty=(subset(ComplexityData3, RugoistyType %in% c("AverageRugoisty")))


## Average 
GardenColor=("#7CAE00")   #### Green 
FarmColor=("#0d7dd9")     #### Blue 
ReferenceGardenColor=("#fa6920") #### Orange
ReferenceFarmColor=("#ff1100") #### Red

col2=c(GardenColor,ReferenceGardenColor,FarmColor,ReferenceFarmColor)



##        Averaged by Site Type
levels(AverageRugoisty$SiteType)[levels(AverageRugoisty$SiteType)=="CG"]<- "Clam Garden"
levels(AverageRugoisty$SiteType)[levels(AverageRugoisty$SiteType)=="CGRef"]<- "Garden Reference"
levels(AverageRugoisty$SiteType)[levels(AverageRugoisty$SiteType)=="SF"]<- "Shellfish Farm"
levels(AverageRugoisty$SiteType)[levels(AverageRugoisty$SiteType)=="SFRef"]<- "Farm Reference"

##   Rugoisty
ggplot(AverageRugoisty, aes(x= SiteType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = 0.3,width = 0.5, position = position_dodge(width = 0.5))+scale_fill_manual(values=col2)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("Linear Rugosity")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")

##  FDaverages
FDaverages$RugoistyType <- factor(FDaverages$RugoistyType, levels=c("FractalD2.5Average","FractalD5Average","FractalD25Average","FractalD50Average"))

levels(FDaverages$RugoistyType)[levels(FDaverages$RugoistyType)=="FractalD2.5Average"]<- "2.5-5"
levels(FDaverages$RugoistyType)[levels(FDaverages$RugoistyType)=="FractalD5Average"]  <- "5-25"
levels(FDaverages$RugoistyType)[levels(FDaverages$RugoistyType)=="FractalD25Average"] <- "25-50"
levels(FDaverages$RugoistyType)[levels(FDaverages$RugoistyType)=="FractalD50Average"] <- "50-100"

ggplot(FDaverages, aes(x= RugoistyType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) +geom_boxplot(size = 0.3,width = 0.5, position = position_dodge(width = 0.5))+scale_fill_manual(values=col2)+ylim(1.99,2.09)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Spatial Scale (cm)")+ ylab("Fractal Dimension (D)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) 

#### Supplemental Complexity by site type and region   ####

## Average by region and site type
FDaverages$Region
BaynesFD=(subset(FDaverages, Region %in% c("Baynes")))
CalvertFD=(subset(FDaverages, Region %in% c("Calvert")))
QuadraFD=(subset(FDaverages, Region %in% c("Quadra")))
QuadraGardensFD=(subset(QuadraFD, SiteType %in% c("CG","CGRef")))
QuadraFarmsFD=(subset(QuadraFD, SiteType %in% c("SF","SFRef")))

AverageRugoisty
BaynesRug=(subset(AverageRugoisty, Region %in% c("Baynes")))
CalvertRug=(subset(AverageRugoisty, Region %in% c("Calvert")))
QuadraRug=(subset(AverageRugoisty, Region %in% c("Quadra")))
QuadraGardensRug=(subset(QuadraRug, SiteType %in% c("CG","CGRef")))
QuadraFarmsRug=(subset(QuadraRug, SiteType %in% c("SF","SFRef")))

##
FarmColors=c(FarmColor,ReferenceFarmColor)
GardenColors=c(GardenColor,ReferenceGardenColor)


## Fractal Dimension
BaynesFDPlot=ggplot(BaynesFD, aes(x= RugoistyType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = 0.3,width = 0.5, position = position_dodge(width = 0.5))+scale_fill_manual(values=FarmColors)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Spatial Scale Range (cm)")+ ylab("Fractal Dimension (D)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) 

QuadraFarmsFDPlot=ggplot(QuadraFarmsFD, aes(x= RugoistyType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, width = 0.5, position = position_dodge(width = 0.5))+scale_fill_manual(values=FarmColors)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Spatial Scale Range (cm)")+ ylab("Fractal Dimension (D)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) 

CalvertFDPlot=ggplot(CalvertFD, aes(x= RugoistyType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, width = 0.5, position = position_dodge(width = 0.5))+scale_fill_manual(values=GardenColors)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Spatial Scale Range (cm)")+ ylab("Fractal Dimension (D)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) 

QuadraGardensFDPlot=ggplot(QuadraGardensFD, aes(x= RugoistyType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, width = 0.5, position = position_dodge(width = 0.5))+scale_fill_manual(values=GardenColors)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Spatial Scale Range (cm)")+ ylab("Fractal Dimension (D)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) 

## Rugoisty
BaynesRug=(subset(AverageRugoisty, Region %in% c("Baynes")))
CalvertRug=(subset(AverageRugoisty, Region %in% c("Calvert")))
QuadraRug=(subset(AverageRugoisty, Region %in% c("Quadra")))
QuadraGardensRug=(subset(QuadraRug, SiteType %in% c("Clam Garden","Garden Reference")))
QuadraFarmsRug=(subset(QuadraRug, SiteType %in% c("Shellfish Farm","Farm Reference")))

BaynesRugPlot=ggplot(BaynesRug, aes(x= SiteType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, width = 0.2, position = position_dodge(width = 0.5))+scale_fill_manual(values=FarmColors)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("Linear Rugosity")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")

CalvertRugPlot=ggplot(CalvertRug, aes(x= SiteType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, width = 0.2, position = position_dodge(width = 0.5))+scale_fill_manual(values=GardenColors)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("Linear Rugosity")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")

QuadraGardensRugPlot=ggplot(QuadraGardensRug, aes(x= SiteType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, width = 0.2, position = position_dodge(width = 0.5))+scale_fill_manual(values=GardenColors)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("Linear Rugosity")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")

QuadraFarmsRugPlot=ggplot(QuadraFarmsRug, aes(x= SiteType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) + geom_boxplot(size = .3, width = 0.2, position = position_dodge(width = 0.5))+scale_fill_manual(values=FarmColors)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab("Linear Rugosity")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")

#### 

BaynesFDPlot
QuadraFarmsFDPlot
CalvertFDPlot
QuadraGardensFDPlot

BaynesRugPlot
CalvertRugPlot
QuadraGardensRugPlot
QuadraFarmsRugPlot
##
library(cowplot)
CalvertPlots=plot_grid(CalvertFDPlot, CalvertRugPlot)
QuadraGardenPlots=plot_grid(QuadraGardensFDPlot, QuadraGardensRugPlot)
QuadraFarmsPlots=plot_grid(QuadraFarmsFDPlot, QuadraFarmsRugPlot)
BaynesPlots=plot_grid(BaynesFDPlot, BaynesRugPlot)


tiff(filename = "ComplexitySiteTypeRegion.tiff",
     width = 60, height = 35, units = "cm",
     compression = c("none"),
     bg = "white", res = 100,
     type = c("quartz"))
plot_grid(CalvertPlots,QuadraGardenPlots,QuadraFarmsPlots,BaynesPlots,
          align = c("hv"),ncol=2, label_size = 25, labels="AUTO",scale = 0.85)
dev.off()
### 

ggarrange(CalvertPlots,QuadraGardenPlots,QuadraFarmsPlots,QuadraFarmsPlots, 
          ncol = 2,nrow = 2, labels = c("A", "B", "C", "D")) 





#### Figure 5: 3D model analysis. Complexity influnece on site diversity ####

#### Complexity of regions and site types ####

##  FDaverages

ggplot(FDaverages, aes(x= RugoistyType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 15) +geom_boxplot(size = 0.3,width = 0.5, position = position_dodge(width = 0.5))+scale_fill_manual(values=col2)+ylim(1.99,2.09)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Spatial Scale (cm)")+ ylab("Fractal Dimension (D)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) +
  guides(colour = FALSE, fill = FALSE) 

library(plotrix)
FDaveragesMeanSE=(aggregate(Complexity ~ RugoistyType+SiteType, data = FDaverages, 
          FUN = function(x) c(mean = mean(x), se = std.error(x))))

FDaveragesMeanSE$mean=FDaveragesMeanSE$Complexity[,1]
FDaveragesMeanSE$es=FDaveragesMeanSE$Complexity[,2]
names(FDaveragesMeanSE)


FDaveragesMeanSEPlot=ggplot(FDaveragesMeanSE, aes(x= RugoistyType, y=mean, color=SiteType,group=SiteType))+
  theme_classic(base_size = 20)+
  geom_errorbar(aes(ymin=mean-es, ymax=mean+es),size=1.5, width=0.6,position = position_dodge(width = 0.5))+
  scale_color_manual(values=col2)+ylim(1.995,2.05)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18)) +
  xlab("Spatial Scale (cm)")+ ylab("Fractal Dimension (D)")+ 
  theme(axis.text.y = element_text(colour = "black", size = 18), axis.title.y = element_text(size = 22)) +
  guides(colour = FALSE, fill = FALSE)+
  geom_line(linetype="dashed", size=1.2,position = position_dodge(width = 0.5))+
  geom_point(size = 4 ,position = position_dodge(width = 0.5))+
  theme(plot.margin = unit(c(1.05,1.05,1.05,1.05), "lines"))

##   Rugoisty
AverageRugoistyPlot=ggplot(AverageRugoisty, aes(x= SiteType, y= Complexity, fill=SiteType))+
  theme_classic(base_size = 20) + geom_boxplot(size = 0.2,width = 0.3, position = position_dodge(width = 0.5))+scale_fill_manual(values=col2)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18)) +
  xlab("")+ ylab("Linear Rugosity")+ theme(axis.text.y = element_text(colour = "black", size = 18), axis.title.y = element_text(size = 22)) +
  guides(colour = FALSE, fill = FALSE) +scale_y_continuous(trans = "reverse")+
  theme(plot.margin = unit(c(1.05,1.05,1.05,1.05), "lines"))

plot_grid(FDaveragesMeanSEPlot, AverageRugoistyPlot)

tiff(filename = "TypeRegionFDRegPlot.tiff",
     width = 60, height = 20, units = "cm",
     compression = c("none"),
     bg = "white", res = 200,
     type = c("quartz"))
plot_grid(FDaveragesMeanSEPlot, AverageRugoistyPlot,
          align = c("hv"),ncol=2, label_size = 25,scale = 0.85)
dev.off()





#### Figure 5 MRF: Complexity and biodiversity #####

## Multivariate Random Forest on higher order taxa (24 taxa), same as previous analyses

### Objectives 

# Supplemental Figure: Variable Importance for all taxa

# Figure 5: Averaged Importance

# Supplemental tables: Explained variablity, error and predicted variability table

### Data
options(digits=5)
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Regression Tree Data ")
## Biodiversity data. 24 taxa
EpifaunaDataRegression=read.csv("SFandCGBiodiversityDataNov2020 RegressionTree Family V8.csv", header=T)
EpifaunaDataRegression

# Mid intertidal. Average within and then among years
EpifaunaDataRegression   ## data
EpifaunaDataRegression %>% count(Site, Year, Zone)
EpifaunaDataMid <- subset(EpifaunaDataRegression, EpifaunaDataRegression$Zone == "Mid")
names(EpifaunaDataMid)
EpifaunaDataMid %>% count(Site, Year,Zone, SiteType,SiteType2)  # 2-5 Quadrats per mid zone

names(EpifaunaDataRegression)
EpifaunaDataRegression[,25:48] #24 taxa
EpifaunaDataRegression[,1:3] 
EpifaunaDataRegression[,7:8] 

EpifaunaDataComplexity=cbind(EpifaunaDataRegression[,1:3],EpifaunaDataRegression[,7:8],EpifaunaDataRegression[,25:48])
EpifaunaDataComplexity
names(EpifaunaDataComplexity)

## Average within year
EpifaunaDataComplexityAverage1=aggregate(. ~ Year+Area+Site+SiteType+SiteType2, data = EpifaunaDataComplexity, mean)

## Average among years
EpifaunaDataComplexityAverage2=aggregate(. ~ Area+Site+SiteType+SiteType2, data = EpifaunaDataComplexityAverage1, mean)

##Data
EpifaunaDataComplexityAverage3 = EpifaunaDataComplexityAverage2[,!(names(EpifaunaDataComplexityAverage2) %in% c("Year"))]

## Reordering taxa
names(EpifaunaDataComplexityAverage3)

col_order <- c("Area","Site","SiteType","SiteType2","Amphilepidida", "Amphipoda", "Anthozoa","Asteriidae","Batillariidae","Bryozoa",
               "Cardiida","Cerithiidae","Chitonida","Decapod","Dendrasteridae","Haminoeidae",
               "Littorinidae","Lottiidae","Myrrhinidae","Mytilidae","Nemertea","Neogastropoda",
               "Ostreidae","Polychaeta","Sessilia","Sphaeromatidae","Strongylocentrotidae","Venerida")

EpifaunaDataComplexityAverage3 <- EpifaunaDataComplexityAverage3[,col_order]
EpifaunaDataComplexityAverage3

### Adding Complexity Data 
# Expoort biodiversity
# write.csv(EpifaunaDataComplexityAverage3,"EpifaunaDataComplexityAverage.csv")
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")
# Read in complexity and biodiversity combined
EpifaunaDataComplexityCombined=read.csv("EpifaunaDataComplexityFinal.csv",header=T)
str(EpifaunaDataComplexityCombined)
EpifaunaDataComplexityCombined$Sphaeromatidae=as.numeric(EpifaunaDataComplexityFinal$Sphaeromatidae)
?tune.nodesize
library(randomForestSRC)
EpifaunaDataComplexityCombined
str(EpifaunaDataComplexityCombined)
names(EpifaunaDataComplexityCombined)

EpifaunaDataComplexityCombined

## Previus taxa data was log +1 transformed, keeping the same
names(EpifaunaDataComplexityCombined)

names(EpifaunaDataComplexityCombined[,5:28])

EpifaunaDataComplexityCombined2= (log10((EpifaunaDataComplexityCombined[,5:28])+1))  
names(EpifaunaDataComplexityCombined[,1:4])
names(EpifaunaDataComplexityCombined[,35:39])
names(EpifaunaDataRegression[5])
names(EpifaunaDataRegression[8])

EpifaunaDataComplexityFinal=cbind((EpifaunaDataComplexityCombined[,1:4]), (log10( (EpifaunaDataComplexityCombined[,5:28]) +1)), (EpifaunaDataComplexityCombined[,35:39])  )

names(EpifaunaDataComplexityFinal)
colnames(EpifaunaDataComplexityFinal)[colnames(EpifaunaDataComplexityFinal)=="SiteType"] = "SiteTypePrevious"
colnames(EpifaunaDataComplexityFinal)[colnames(EpifaunaDataComplexityFinal)=="SiteType2"] = "SiteType"
names(EpifaunaDataComplexityFinal)
EpifaunaDataComplexityFinal$SiteType

## Node size tuned
tune.nodesize(Multivar(Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,
                       Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,
                       Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia, Sphaeromatidae,Strongylocentrotidae,
                       Venerida) ~SiteType+Area+AverageRugosity+FractalD50Average+FractalD25Average+FractalD5Average+FractalD2.5Average, 
              data = EpifaunaDataComplexityFinal,ntree=6000,block.size=1,nodedepth=6, importance = TRUE) 

MRT.Complexity <- rfsrc(Multivar(Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,
                                  Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,
                                  Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia, Sphaeromatidae,Strongylocentrotidae,
                                  Venerida) ~SiteType+Area+AverageRugosity+FractalD50Average+FractalD25Average+FractalD5Average+FractalD2.5Average, 
                         data = EpifaunaDataComplexityCombined,ntree=6000,nodesize=3,block.size=1,nodedepth=6, importance = TRUE) 

plot(get.tree(MRT.Complexity, 5))
plot(get.tree(MRT.Complexity, 25))
plot(get.tree(MRT.Complexity, 50))
plot(get.tree(MRT.Complexity, 500))
plot(get.tree(MRT.Complexity, 5000))
### Note 
# ntree   Number of trees   -- Determined by looking at if error rates stablize

# nodesize    Minumum size of terminal node. The defaults are regression (5), mixed outcomes (3), unsupervised (3)
               ## Danger here is overfitting

# nodedepth   Maximum depth to which a tree should be grown. 
##            The default behaviour is that this parameter is ignored.

# block.size  Cumulative error rate be calculated on how many trees. 
##            VIMP is calculated in "blocks" of size equal to block.size,
#             thus resulting in a useful compromise between ensemble and permutation VIMP. 
#             The default action is to use 10 trees.
            
# bootstrap default is from root without replacement, which is fine



# Supplemental Figure: Variable Importance for all taxa

# Figure 6: Averaged Importance

# Supplemental tables: Explained variablity, error and predicted variability table

err <- get.mv.error(MRT.Complexity)
vmp <- get.mv.vimp(MRT.Complexity)
pred <- get.mv.predicted(MRT.Complexity)
err.std <- get.mv.error(MRT.Complexity, standardize = TRUE)
vmp.std <- get.mv.vimp(MRT.Complexity, standardize = TRUE)


plot(MRT.Complexity,m.target="Amphilepidida")
plot(MRT.Complexity,m.target="Amphipoda")
plot(MRT.Complexity,m.target="Anthozoa")
plot(MRT.Complexity,m.target="Asteriidae")
plot(MRT.Complexity,m.target="Batillariidae")
plot(MRT.Complexity,m.target="Bryozoa")
plot(MRT.Complexity,m.target="Cardiida")
plot(MRT.Complexity,m.target="Cerithiidae")
plot(MRT.Complexity,m.target="Chitonida")
plot(MRT.Complexity,m.target="Decapod")
plot(MRT.Complexity,m.target="Dendrasteridae")
plot(MRT.Complexity,m.target="Haminoeidae")
plot(MRT.Complexity,m.target="Littorinidae")
plot(MRT.Complexity,m.target="Lottiidae")
plot(MRT.Complexity,m.target="Myrrhinidae")
plot(MRT.Complexity,m.target="Mytilidae")
plot(MRT.Complexity,m.target="Nemertea")
plot(MRT.Complexity,m.target="Neogastropoda")
plot(MRT.Complexity,m.target="Ostreidae")
plot(MRT.Complexity,m.target="Polychaeta")
plot(MRT.Complexity,m.target="Sessilia")
plot(MRT.Complexity,m.target="Sphaeromatidae")
plot(MRT.Complexity,m.target="Strongylocentrotidae")
plot(MRT.Complexity,m.target="Venerida")


print(MRT.Complexity,outcome.target="Amphilepidida")
print(MRT.Complexity,outcome.target="Amphipoda")
print(MRT.Complexity,outcome.target="Anthozoa")
print(MRT.Complexity,outcome.target="Asteriidae")
print(MRT.Complexity,outcome.target="Batillariidae")
print(MRT.Complexity,outcome.target="Bryozoa")
print(MRT.Complexity,outcome.target="Cardiida")
print(MRT.Complexity,outcome.target="Cerithiidae")
print(MRT.Complexity,outcome.target="Chitonida")
print(MRT.Complexity,outcome.target="Decapod")
print(MRT.Complexity,outcome.target="Dendrasteridae")
print(MRT.Complexity,outcome.target="Haminoeidae")
print(MRT.Complexity,outcome.target="Littorinidae")
print(MRT.Complexity,outcome.target="Lottiidae")
print(MRT.Complexity,outcome.target="Myrrhinidae")
print(MRT.Complexity,outcome.target="Mytilidae")
print(MRT.Complexity,outcome.target="Nemertea")
print(MRT.Complexity,outcome.target="Neogastropoda")
print(MRT.Complexity,outcome.target="Ostreidae")
print(MRT.Complexity,outcome.target="Polychaeta")
print(MRT.Complexity,outcome.target="Sessilia")
print(MRT.Complexity,outcome.target="Sphaeromatidae")
print(MRT.Complexity,outcome.target="Strongylocentrotidae")
print(MRT.Complexity,outcome.target="Venerida")


#### Extract Each Explained Variance and Test Set Error Rate ####

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/MRF Complexity Data")

## Explained Variance
vmp.std <- get.mv.vimp(MRT.Complexity, standardize = TRUE)

ComplexityError <- get.mv.error(MRT.Complexity)
# write.csv(ComplexityError,"ComplexityError.csv")

ComplexityError.std <- get.mv.error(MRT.Complexity, standardize = TRUE)
# write.csv(ComplexityError.std,"ComplexityError.std.csv")

#### Plotting VariableImportance ####
MRT.Complexity.VIMP.Std <- get.mv.vimp(MRT.Complexity, standardize = TRUE)

# write.csv(MRT.Complexity.VIMP.Std,"MRT.Complexity.VIMP.Std.csv")
ComplexityVariableImportance=read.csv("MRT.Complexity.VIMP.Std.Final.csv",header=T)

ComplexityVariableImportance
names(ComplexityVariableImportance)
head(ComplexityVariableImportance)
ComplexityVariableImportance$Variable
str(ComplexityVariableImportance)
names(ComplexityVariableImportance)

### species plots ###
size3=12

# Amphilepidida
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Amphilepidida)])
Amphilepidida.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Amphilepidida,fill=AmphilepididaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Amphilepidida.Plot

# Amphipoda
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Amphipoda)])
Amphipoda.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Amphipoda,fill=AmphipodaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Amphipoda.Plot

# Anthozoa
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Anthozoa)])
Anthozoa.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Anthozoa,fill=AnthozoaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Anthozoa.Plot

# Asteriidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Asteriidae)])
Asteriidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Asteriidae,fill=AsteriidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Asteriidae.Plot

# Batillariidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Batillariidae)])
Batillariidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Batillariidae,fill=BatillariidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3","red2"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Batillariidae.Plot

# Bryozoa
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Bryozoa)])
Bryozoa.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Bryozoa,fill=BryozoaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Bryozoa.Plot

# Cardiida
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Cardiida)])
Cardiida.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Cardiida,fill=CardiidaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Cardiida.Plot

# Cerithiidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Cerithiidae)])
Cerithiidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Cerithiidae,fill=CerithiidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Cerithiidae.Plot

# Chitonida
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Chitonida)])
Chitonida.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Chitonida,fill=ChitonidaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Chitonida.Plot

# Decapod
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Decapod)])
Decapod.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Decapod,fill=DecapodDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Decapod.Plot

# Dendrasteridae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Dendrasteridae)])
Dendrasteridae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Dendrasteridae,fill=DendrasteridaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Dendrasteridae.Plot

# Haminoeidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Haminoeidae)])
Haminoeidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Haminoeidae,fill=HaminoeidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Haminoeidae.Plot

# Littorinidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Littorinidae)])
Littorinidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Littorinidae,fill=LittorinidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Littorinidae.Plot

# Lottiidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Lottiidae)])
Lottiidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Lottiidae,fill=LottiidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Lottiidae.Plot

# Myrrhinidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Myrrhinidae)])
Myrrhinidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Myrrhinidae,fill=MyrrhinidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Myrrhinidae.Plot

# Mytilidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Mytilidae)])
Mytilidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Mytilidae,fill=MytilidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Mytilidae.Plot

# Nemertea
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Nemertea)])
Nemertea.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Nemertea,fill=NemerteaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Nemertea.Plot

# Neogastropoda
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Neogastropoda)])
Neogastropoda.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Neogastropoda,fill=NeogastropodaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Neogastropoda.Plot

# Ostreidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Ostreidae)])
Ostreidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Ostreidae,fill=OstreidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Ostreidae.Plot

# Polychaeta
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Polychaeta)])
Polychaeta.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Polychaeta,fill=PolychaetaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Polychaeta.Plot

# Sessilia
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Sessilia)])
Sessilia.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Sessilia,fill=SessiliaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3","red2"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Sessilia.Plot

# Sphaeromatidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Sphaeromatidae)])
Sphaeromatidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Sphaeromatidae,fill=SphaeromatidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Sphaeromatidae.Plot

# Strongylocentrotidae
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Strongylocentrotidae)])
Strongylocentrotidae.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Strongylocentrotidae,fill=StrongylocentrotidaeDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Strongylocentrotidae.Plot

# Venerida
ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$Venerida)])
Venerida.Plot=ggplot(ComplexityVariableImportance, aes(x = Variable, y = Venerida,fill=VeneridaDirectionality))+ 
  theme_classic() + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=2)+ylab("Variable Importance")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = size3),axis.title.x  = element_text(size = 15)) +scale_fill_manual(values = c("red2","blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = size3), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
Venerida.Plot



## Standardized Average Variable Importance##
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Final Figures V1/Final Figures V2/Complexity Data Figures")


FDaveragesMeanSEPlot=ggplot(FDaveragesMeanSE, aes(x= RugoistyType, y=mean, color=SiteType,group=SiteType))+
  theme_classic(base_size = 20)+
  geom_errorbar(aes(ymin=mean-es, ymax=mean+es),size=1.5, width=0.6,position = position_dodge(width = 0.5))+
  scale_color_manual(values=col2)+ylim(1.995,2.05)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18)) +
  xlab("Spatial Scale (cm)")+ ylab("Fractal Dimension (D)")+ 
  theme(axis.text.y = element_text(colour = "black", size = 18), axis.title.y = element_text(size = 22)) +
  guides(colour = FALSE, fill = FALSE)+
  geom_line(linetype="dashed", size=1.2,position = position_dodge(width = 0.5))+
  geom_point(size = 4 ,position = position_dodge(width = 0.5))+
  theme(plot.margin = unit(c(1.05,1.05,1.05,1.05), "lines"))


ComplexityVariableImportance$Variable <- factor(ComplexityVariableImportance$Variable, levels = ComplexityVariableImportance$Variable[order(ComplexityVariableImportance$AveragePercent)])
AveragePercent.Plot2=ggplot(ComplexityVariableImportance, aes(x = Variable, y = AveragePercent,fill=Directionality))+ 
  theme_classic(base_size = 20) + ggtitle("") +geom_bar(stat = "identity",width = 0.1)+geom_point(size=3)+ylab("Variable Importance (%)")+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 18),axis.title.x  = element_text()) +scale_fill_manual(values = c("blue3"))+
  theme(axis.text.y = element_text(colour = "black", size = 18), axis.title.y = element_text(size = 22))+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "gray5",size = 0.75,linetype=1)+coord_flip()
AveragePercent.Plot2


tiff(filename = "ComplexityVariableImportance3.tiff",
     width = 30, height = 20, units = "cm",
     compression = c("none"),
     bg = "white", res = 200)
AveragePercent.Plot2
dev.off()

citation()
###
Amphilepidida.Plot
Amphipoda.Plot
Anthozoa.Plot
Asteriidae.Plot
Batillariidae.Plot
Bryozoa.Plot
Cardiida.Plot
Cerithiidae.Plot
Chitonida.Plot
Decapod.Plot
Dendrasteridae.Plot
Haminoeidae.Plot
Littorinidae.Plot
Lottiidae.Plot
Myrrhinidae.Plot
Mytilidae.Plot
Nemertea.Plot
Neogastropoda.Plot
Ostreidae.Plot
Polychaeta.Plot
Sessilia.Plot
Sphaeromatidae.Plot
Strongylocentrotidae.Plot
Venerida.Plot


library(cowplot)

## Supplemental plot ##
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results/Final Figures V1/Final Figures V2/Complexity Data Figures")
jpeg(filename = "ComplexityVariableImportance.jpeg",width = 50, height = 80, units = "cm", pointsize = 9, quality = 100,res=800)
plot_grid(
  
  
  Amphilepidida.Plot,    Amphipoda.Plot,             Anthozoa.Plot,
  Asteriidae.Plot,       Batillariidae.Plot,         Bryozoa.Plot,
  Cardiida.Plot,         Cerithiidae.Plot,           Chitonida.Plot,
  Decapod.Plot,          Dendrasteridae.Plot,        Haminoeidae.Plot,
  Littorinidae.Plot,     Lottiidae.Plot,             Myrrhinidae.Plot,
  Mytilidae.Plot,        Nemertea.Plot,              Neogastropoda.Plot,
  Ostreidae.Plot,        Polychaeta.Plot,            Sessilia.Plot,
  Sphaeromatidae.Plot,   Strongylocentrotidae.Plot,  Venerida.Plot,
 
          labels = c('A) Amphilepidida' ,'B) Amphipoda',           'C) Anthozoa',
                     'D) Asteriidae',    'E) Batillariidae',       'F) Bryozoa ',
                     'G) Cardiida',      'H) Cerithiidae',         'I) Chitonida.Plot',
                     'J) Decapod',       'K) Dendrasteridae',      'L) Haminoeidae',
                     'M) Littorinidae',  'N) Lottiidae',           'O) Myrrhinidae',
                     'P) Mytilidae',     'Q) Nemertea',            'R) Neogastropoda',
                     'S) Ostreidae',     'T) Polychaeta',          'U) Sessilia',
                     'V) Sphaeromatidae','W) Strongylocentrotidae','X) Venerida'),nrow = 8, align="hv")
dev.off()




#### Varability Explained and Random Forest Predictions ####

## Error
MRT.FormAllData2 <- rfsrc(Multivar(Bryozoa,Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Cardiida,Cerithiidae,
                                   Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,
                                   Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae,
                                   Venerida) ~SiteType2+Area+Site+Zone+Cobble+Gravel+Bivalve.Shells+Mastocarpus+Ulva+Oyster.Shell+Mud.silt+Sanddollar.dead+Gracilaria, 
                          data = TaxaData3,ntree=6000,block.size=1,nodedepth=6, importance = TRUE) 
MRT.FormAllData2

err2 <- get.mv.error(MRT.FormAllData2)
vmp2 <- get.mv.vimp(MRT.FormAllData2)
pred2 <- get.mv.predicted(MRT.FormAllData2)
err.std2 <- get.mv.error(MRT.FormAllData2, standardize = TRUE)
vmp.std2 <- get.mv.vimp(MRT.FormAllData2, standardize = TRUE)

### Explained Variance
MRT.FormAllData2 <- rfsrc(Multivar(Venerida,Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,
                                   Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,
                                   Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae) ~SiteType2+Area+Site+Zone+Cobble+Gravel+Bivalve.Shells+Mastocarpus+Ulva+Oyster.Shell+Mud.silt+Sanddollar.dead+Gracilaria, 
                          data = TaxaData3,ntree=6000,block.size=1,nodedepth=6, importance = TRUE) 
MRT.FormAllData2

### Prediction
# train <- sample(1:nrow(TaxaData3), round(nrow(TaxaData3) * 0.80))
MRT.FormAllDataTrain <- rfsrc(Multivar(Venerida,Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,
                                       Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,
                                       Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia,Sphaeromatidae,Strongylocentrotidae) ~SiteType2+Area+Site+Zone+Cobble+Gravel+Bivalve.Shells+Mastocarpus+Ulva+Oyster.Shell+Mud.silt+Sanddollar.dead+Gracilaria, 
                              data = TaxaData3[train, ],ntree=6000,block.size=1,nodedepth=6, importance = TRUE) 

MRT.FormAllDataPredict<- predict(MRT.FormAllDataTrain, TaxaData3[-train , ])
print(MRT.FormAllDataTrain)
print(MRT.FormAllDataPredict)

### Prediction
train <- sample(1:nrow(EpifaunaDataComplexityCombined), round(nrow(EpifaunaDataComplexityCombined) * 0.80))


MRT.ComplexityTrain <- rfsrc(Multivar(Amphilepidida,Amphipoda,Anthozoa,Asteriidae,Batillariidae,Bryozoa,Cardiida,Cerithiidae,
                                 Chitonida,Decapod,Dendrasteridae,Haminoeidae,Littorinidae,Lottiidae,Myrrhinidae,Mytilidae,
                                 Nemertea,Neogastropoda,Ostreidae,Polychaeta,Sessilia, Sphaeromatidae,Strongylocentrotidae,
                                 Venerida) ~SiteType+Area+AverageRugosity+FractalD50Average+FractalD25Average+FractalD5Average+FractalD2.5Average, 
                        data = EpifaunaDataComplexityCombined[train, ] ,ntree=6000,nodesize=3,block.size=1,nodedepth=6, importance = TRUE) 


MRT.ComplexityPredict<- predict(MRT.ComplexityTrain, EpifaunaDataComplexityCombined[-train , ])
print(MRT.ComplexityTrain)
print(MRT.ComplexityPredict)


print(MRT.ComplexityPredict,outcome.target="Amphilepidida")
print(MRT.ComplexityPredict,outcome.target="Amphipoda")
print(MRT.ComplexityPredict,outcome.target="Anthozoa")
print(MRT.ComplexityPredict,outcome.target="Asteriidae")
print(MRT.ComplexityPredict,outcome.target="Batillariidae")
print(MRT.ComplexityPredict,outcome.target="Bryozoa")
print(MRT.ComplexityPredict,outcome.target="Cardiida")
print(MRT.ComplexityPredict,outcome.target="Cerithiidae")
print(MRT.ComplexityPredict,outcome.target="Chitonida")
print(MRT.ComplexityPredict,outcome.target="Decapod")
print(MRT.ComplexityPredict,outcome.target="Dendrasteridae")
print(MRT.ComplexityPredict,outcome.target="Haminoeidae")
print(MRT.ComplexityPredict,outcome.target="Littorinidae")
print(MRT.ComplexityPredict,outcome.target="Lottiidae")
print(MRT.ComplexityPredict,outcome.target="Myrrhinidae")
print(MRT.ComplexityPredict,outcome.target="Mytilidae")
print(MRT.ComplexityPredict,outcome.target="Nemertea")
print(MRT.ComplexityPredict,outcome.target="Neogastropoda")
print(MRT.ComplexityPredict,outcome.target="Ostreidae")
print(MRT.ComplexityPredict,outcome.target="Polychaeta")
print(MRT.ComplexityPredict,outcome.target="Sessilia")
print(MRT.ComplexityPredict,outcome.target="Sphaeromatidae")
print(MRT.ComplexityPredict,outcome.target="Strongylocentrotidae")
print(MRT.ComplexityPredict,outcome.target="Venerida")


## 




















#### Old code ####
## Important questions, impacts of averaging data?

CalvertCGMidEpifaunaData3$Quad

CalvertCGMidEpifaunaData4=aggregate(.~SiteType+Year+Site+Quad, CalvertCGMidEpifaunaData3, mean)
str(CalvertCGMidEpifaunaData4)
names(CalvertCGMidEpifaunaData4)
CalvertCGMidEpifaunaData5  =  CalvertCGMidEpifaunaData4[, colSums(CalvertCGMidEpifaunaData4 != 0) > 0]
str(CalvertCGMidEpifaunaData5)
names(CalvertCGMidEpifaunaData5)

names(CalvertCGMidEpifaunaData5)
CalvertCGMidEpifaunaStructure=CalvertCGMidEpifaunaData5[,1:4]
CalvertCGMidEpifaunaSpeciesData <- CalvertCGMidEpifaunaData5[,6:38]
CalvertCGMidEpifaunaSpeciesData2=((CalvertCGMidEpifaunaSpeciesData+1))

###### Simper
SpeciesSIMPER<-(simper(CalvertCGMidEpifaunaSpeciesData2, CalvertCGMidEpifaunaStructure$SiteType, permutations=999))
summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$`Clam Garden_Reference`[,1]), na.rm = TRUE)

### so the communities are 68% dissimilar 
summary=summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response

summary2=summary$`Clam Garden_Reference`

write.csv(summary2,"CalvertEpifaunaSIMPERS2.csv")
EpifaunaSIMPERCalvert=read.csv("CalvertEpifaunaSIMPERSFinal3.csv",header=T)
EpifaunaSIMPERCalvert
names(EpifaunaSIMPERCalvert)
str(EpifaunaSIMPERCalvert)
EpifaunaSIMPERCalvert$Species <- factor(EpifaunaSIMPERCalvert$Species, levels = EpifaunaSIMPERCalvert$Species[order(EpifaunaSIMPERCalvert$Average.Dissim.....)])

col=c( "#F8766D", "#7CAE00" )
# base_size = 15
p=ggplot(EpifaunaSIMPERCalvert, aes(x = Species, y = Average.Dissim..... ,colour=Direction))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity ")+ ylim(-35,35)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=col)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)

AdultSimperPlot=p+coord_flip()
AdultSimperPlot

library(grid)

# Create a text
grob <- grobTree(textGrob("Dissimilarity 68.70%", x=0.005,  y=1.05, hjust=0,
                          gp=gpar(col="Black", fontsize=14)))
AdultSimperPlot + annotation_custom(grob)+geom_line(linetype = 2)+geom_vline(xintercept = 33.6,size=0.6)

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "EpifaunaSimperPlotV1.jpeg", width = 25, height = 18, units = "cm", pointsize = 15, res = 300)  

AdultSimperPlot + annotation_custom(grob)+geom_line(linetype = 2)+geom_vline(xintercept = 33.6,size=0.7)

dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")










#### Supplemental: Permanova tables ####





















#### old code

EpifaunaSpeciesData7

### EpifaunaSpeciesData7=aggregate(.~SiteType+SiteType2+Site, EpifaunaSpeciesData6, mean)
str(EpifaunaSpeciesData7)
# view(EpifaunaSpeciesData7)

names(EpifaunaSpeciesData7)
MidEpifaunaDataStructure2 =  EpifaunaSpeciesData7[,1:6]
MidEpifaunaSpeciesData2   =  EpifaunaSpeciesData7[,9:59]
names(MidEpifaunaDataStructure2)
names(MidEpifaunaSpeciesData2)

### Transform and dummy    ^(1/4) doesn't have to occur at site average
MidEpifaunaSpeciesData3=((MidEpifaunaSpeciesData2+1))

nMDSCommunityMid<- metaMDS(MidEpifaunaSpeciesData3, distance = "bray", k = 2) 
nMDSCommunityMid$stress  
stressplot(nMDSCommunityMid)

names(MidEpifaunaDataStructure)
co=c(GardenColor,ReferenceGardenColor, ReferenceFarmColor,FarmColor)
shape=c(1:24) 
levels(MidEpifaunaDataStructure2$SiteType2)

plot(nMDSCommunityMid$points, col=co[MidEpifaunaDataStructure2$SiteType2],  pch = shape[MidEpifaunaDataStructure2$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")

ordihull(nMDSCommunityMid,groups=MidEpifaunaDataStructure2$SiteType2,draw="polygon",col=co,label=T,alpha = 40)
txt <- c("Clam Garden", "Farm Reference","Garden Reference","Shellfish Farm")
legend('topright', txt , pch=c(20),col=co,cex=0.8, bty = "y")

### Export it ##

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")

jpeg(filename = "EpifaunaMidnMDS2TypeSite.jpeg", width = 22, height = 15, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(nMDSCommunityMid$points, col=co[MidEpifaunaDataStructure2$SiteType2],  pch = shape[MidEpifaunaDataStructure2$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")

ordihull(nMDSCommunityMid,groups=MidEpifaunaDataStructure2$SiteType2,draw="polygon",col=co,label=F,alpha = 40)
txt <- c("Clam Garden", "Farm Reference","Garden Reference","Shellfish Farm")
legend('topright', txt , pch=c(20),col=co,cex=0.7, bty = "y")

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")





# 4 nMDS Epifauna- Types/Sites plotted by region






















###### Bivalve Diversity Profile Calvert  ####

CalvertCGMidBivalveData3 <- subset(CGMidBivalveSpeciesData3, CGMidBivalveSpeciesData3$Area == "Calvert")

CalvertCGMidBivalveData3$Area
CalvertCGMidBivalveData3$Site
CalvertCGMidBivalveData3$Quad

CalvertCGMidBivalveData4=aggregate(.~SiteType+Site, CalvertCGMidBivalveData3, mean)
CalvertCGMidBivalveData5=aggregate(.~SiteType, CalvertCGMidBivalveData4, mean)
names(CalvertCGMidBivalveData5)

MidBivalveCG= CalvertCGMidBivalveData5[1,]
MidBivalveReference= CalvertCGMidBivalveData5[2,]

### Removing zero's
MidBivalveCG2=cbind(MidBivalveCG[,7:9],MidBivalveCG[,13:16],MidBivalveCG[,18])
MidBivalveReference2=cbind(MidBivalveReference[,7:10], MidBivalveReference[,13:16],MidBivalveReference[,18])
MidBivalveReference2=MidBivalveReference2[,2:10]

### So then using this equation
divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

MidBivalveReference3=divprof(MidBivalveReference2)

### Checking Adult Control Data
exp(diversity(MidBivalveReference2, index="shannon"))
###   3.490625
1 / (1-(diversity(MidBivalveReference2, index="simpson")))
###  2.223225

MidBivalveCG3=divprof(MidBivalveCG2)


### Checking Adult Treatment Data
exp(diversity(MidBivalveCG2, index="shannon"))
###   3.771266
1 / (1-(diversity(MidBivalveCG2, index="simpson")))
###  3.042337


### Building the plot
plot(MidBivalveReference3[,1],MidBivalveReference3[,2], col=alpha("#F8766D",0.6),cex=2, pch=16,xlab="",ylab="",ylim = c(1,10.3))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(MidBivalveCG3[,1],MidBivalveCG3[,2],pch=16,cex=2,col= alpha("#7CAE00",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

###
points(MidBivalveCG3[1,1],MidBivalveCG3[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[1,1],MidBivalveReference3[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.5,8.08, expression( "Richness"["T"] ))
text(0.5,9.08, expression( "Richness"["C"] ))

### Now plotting Control shannon and Simpson

points(MidBivalveReference3[10,1],MidBivalveReference3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[19,1],MidBivalveReference3[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(MidBivalveCG3[10,1],MidBivalveCG3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveCG3[19,1],MidBivalveCG3[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.45,2.8, expression( "Shannon"["C"] ))
text(2.3,1.5, expression( "Simpson"["C"] ))

text(1.45,4, expression( "Shannon"["T"] ))
text(2.4,3.5, expression( "Simpson"["T"] ))

####

### Plot it
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "CalvertBivalveOrderPlot4.jpeg", width = 22, height = 18, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(MidBivalveReference3[,1],MidBivalveReference3[,2], col=alpha("#F8766D",0.6),cex=2, pch=16,xlab="",ylab="",ylim = c(1,10.3))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(MidBivalveCG3[,1],MidBivalveCG3[,2],pch=16,cex=2,col= alpha("#7CAE00",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

### 
points(MidBivalveCG3[1,1],MidBivalveCG3[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[1,1],MidBivalveReference3[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.39,8.33, expression( "Richness"["T"] ))
text(0.39,9.33, expression( "Richness"["C"] ))

### Now plotting Control shannon and Simpson

points(MidBivalveReference3[10,1],MidBivalveReference3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[19,1],MidBivalveReference3[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(MidBivalveCG3[10,1],MidBivalveCG3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveCG3[19,1],MidBivalveCG3[19,2],col= alpha("black",0.9),pch=16,cex=2)

#text(1.4,2.8, expression( "Shannon"["C"] ))

text(.68,3.1, expression( "Shannon"["C"] ))
text(1.4,4.1, expression( "Shannon"["T"] ))

text(2.3,1.6, expression( "Simpson"["C"] ))
text(2.3,3.5, expression( "Simpson"["T"] ))

dev.off()
#### Plot End
























###### Old Code


##### Analysis Objectives
# 1. Clam Garden Analysis
### 1a. Split by Year, Region, Community (Bivavle and Epifauna), mid intertidal only
# Diversity Profiles
# nMDS
# PERMANOVA
# SIMPER
##### Combine years if appropriate
##### Combine communities if appropriate
##### Combine regions if appropriate (keep seperate data for Quadra Aquaculture comparison)

# 2. Shellfish Aqauculture Analysis
### 2a. Split by Year, Region, Community (Bivavle and Epifauna), mid intertidal only
# Diversity Profiles
# nMDS
# PERMANOVA
# SIMPER
##### Combine years if appropriate
##### Combine communities if appropriate
##### Combine regions if appropriate (keep seperate data for Quadra Aquaculture comparison)

# 3. Forces Structing communites within sites- Mid, and High, Mid, and Low
# Models
# Regression Trees



# 1. Clam Garden Analysis
SFandCGData %>% count(Site, Year, Zone)
SFandCGDataMid <- subset(SFandCGData, SFandCGData$Zone == "Mid")
SFandCGDataMid
SFandCGDataMid %>% count(Site, Year)

CGDataMidCG <- subset(SFandCGDataMid, SFandCGDataMid$Comparison == "Garden")
CGDataMidCG
CGDataMidCG %>% count(Site, Year)


################ Bivalves
names(CGDataMidCG)
CGMidBivalveSpeciesData <- CGDataMidCG[,241:254]
CGMidDataStructure1=CGDataMidCG[,1:2]
CGMidDataStructure2=CGDataMidCG[,5:7]
CGMidDataStructure3=cbind(CGMidDataStructure1,CGMidDataStructure2)
CGMidBivalveSpeciesData2=cbind(CGMidDataStructure3,CGMidBivalveSpeciesData)
str(CGMidBivalveSpeciesData2)

## Removing NAs
names(CGMidBivalveSpeciesData2)
str(CGMidBivalveSpeciesData2)
CGMidBivalveSpeciesData3=(CGMidBivalveSpeciesData2[complete.cases(CGMidBivalveSpeciesData2), ])
str(CGMidBivalveSpeciesData3)
CGMidBivalveSpeciesData3
names(CGMidBivalveSpeciesData3)


##### Diversity Profiles, Calvert and Quadra, combined


################ Diversity Profile Curves
CGMidBivalveSpeciesData3$Area

CalvertCGMidBivalveData3 <- subset(CGMidBivalveSpeciesData3, CGMidBivalveSpeciesData3$Area == "Calvert")
QuadraCGMidBivalveData3 <- subset(CGMidBivalveSpeciesData3, CGMidBivalveSpeciesData3$Area == "Quadra")

options(digits=5)
###### Bivalve Diversity Profile Calvert  ####

CalvertCGMidBivalveData3 <- subset(CGMidBivalveSpeciesData3, CGMidBivalveSpeciesData3$Area == "Calvert")

str(CalvertCGMidBivalveData3)
CalvertCGMidBivalveData3$Area
CalvertCGMidBivalveData3$Site
CalvertCGMidBivalveData3$Quad

CalvertCGMidBivalveData4=aggregate(.~SiteType+Site, CalvertCGMidBivalveData3, mean)
CalvertCGMidBivalveData5=aggregate(.~SiteType, CalvertCGMidBivalveData4, mean)
names(CalvertCGMidBivalveData5)

MidBivalveCG= CalvertCGMidBivalveData5[1,]
MidBivalveReference= CalvertCGMidBivalveData5[2,]

### Removing zero's
MidBivalveCG2=cbind(MidBivalveCG[,7:9],MidBivalveCG[,13:16],MidBivalveCG[,18])
MidBivalveReference2=cbind(MidBivalveReference[,7:10], MidBivalveReference[,13:16],MidBivalveReference[,18])
MidBivalveReference2=MidBivalveReference2[,2:10]

### So then using this equation
divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

MidBivalveReference3=divprof(MidBivalveReference2)

### Checking Adult Control Data
exp(diversity(MidBivalveReference2, index="shannon"))
###   3.490625
1 / (1-(diversity(MidBivalveReference2, index="simpson")))
###  2.223225

MidBivalveCG3=divprof(MidBivalveCG2)


### Checking Adult Treatment Data
exp(diversity(MidBivalveCG2, index="shannon"))
###   3.771266
1 / (1-(diversity(MidBivalveCG2, index="simpson")))
###  3.042337


##### Building the plot
plot(MidBivalveReference3[,1],MidBivalveReference3[,2], col=alpha("#F8766D",0.6),cex=2, pch=16,xlab="",ylab="",ylim = c(1,10.3))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(MidBivalveCG3[,1],MidBivalveCG3[,2],pch=16,cex=2,col= alpha("#7CAE00",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

##### 
points(MidBivalveCG3[1,1],MidBivalveCG3[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[1,1],MidBivalveReference3[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.5,8.08, expression( "Richness"["T"] ))
text(0.5,9.08, expression( "Richness"["C"] ))

### Now plotting Control shannon and Simpson

points(MidBivalveReference3[10,1],MidBivalveReference3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[19,1],MidBivalveReference3[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(MidBivalveCG3[10,1],MidBivalveCG3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveCG3[19,1],MidBivalveCG3[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.45,2.8, expression( "Shannon"["C"] ))
text(2.3,1.5, expression( "Simpson"["C"] ))

text(1.45,4, expression( "Shannon"["T"] ))
text(2.4,3.5, expression( "Simpson"["T"] ))

####

#### Plot it
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "CalvertBivalveOrderPlot4.jpeg", width = 22, height = 18, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(MidBivalveReference3[,1],MidBivalveReference3[,2], col=alpha("#F8766D",0.6),cex=2, pch=16,xlab="",ylab="",ylim = c(1,10.3))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(MidBivalveCG3[,1],MidBivalveCG3[,2],pch=16,cex=2,col= alpha("#7CAE00",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

##### 
points(MidBivalveCG3[1,1],MidBivalveCG3[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[1,1],MidBivalveReference3[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.39,8.33, expression( "Richness"["T"] ))
text(0.39,9.33, expression( "Richness"["C"] ))

### Now plotting Control shannon and Simpson

points(MidBivalveReference3[10,1],MidBivalveReference3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[19,1],MidBivalveReference3[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(MidBivalveCG3[10,1],MidBivalveCG3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveCG3[19,1],MidBivalveCG3[19,2],col= alpha("black",0.9),pch=16,cex=2)

#text(1.4,2.8, expression( "Shannon"["C"] ))

text(.68,3.1, expression( "Shannon"["C"] ))
text(1.4,4.1, expression( "Shannon"["T"] ))

text(2.3,1.6, expression( "Simpson"["C"] ))
text(2.3,3.5, expression( "Simpson"["T"] ))

dev.off()


##### Bivalve Diversity Profile Quadra  
QuadraCGMidBivalveData3

str(QuadraCGMidBivalveData3)
QuadraCGMidBivalveData3$Area
QuadraCGMidBivalveData3$Site
QuadraCGMidBivalveData3$Quad

QuadraCGMidBivalveData4=aggregate(.~SiteType+Site, QuadraCGMidBivalveData3, mean)
QuadraCGMidBivalveData5=aggregate(.~SiteType, QuadraCGMidBivalveData4, mean)
names(QuadraCGMidBivalveData5)

QuadraMidBivalveCG= QuadraCGMidBivalveData5[1,]
QuadraMidBivalveReference= QuadraCGMidBivalveData5[2,]

### Removing zero's
QuadraMidBivalveCG2=cbind(QuadraMidBivalveCG[,7:11],QuadraMidBivalveCG[,13:16])
QuadraMidBivalveReference2=cbind(QuadraMidBivalveReference[,7:11], QuadraMidBivalveReference[,13:15])
QuadraMidBivalveReference2
write.csv(QuadraMidBivalveReference2,"QuadraMidBivalveReference2.csv")
QuadraMidBivalveReference2=read.csv("QuadraMidBivalveReference2.csv",header=T)

### So then using this equation
divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

QuadraMidBivalveReference3=divprof(QuadraMidBivalveReference2)

### Checking Adult Control Data
exp(diversity(QuadraMidBivalveReference2, index="shannon"))
###   3.512352
1 / (1-(diversity(QuadraMidBivalveReference2, index="simpson")))
###  2.277252

QuadraMidBivalveCG3=divprof(QuadraMidBivalveCG2)

### Checking Adult Treatment Data
exp(diversity(QuadraMidBivalveCG2, index="shannon"))
###   4.166968
1 / (1-(diversity(QuadraMidBivalveCG2, index="simpson")))
###  2.965608

##### Building the plot
plot(QuadraMidBivalveReference3[,1],QuadraMidBivalveReference3[,2], col=alpha("#F8766D",0.6),cex=2, pch=16,xlab="",ylab="",ylim = c(1,10.3))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(QuadraMidBivalveCG3[,1],QuadraMidBivalveCG3[,2],pch=16,cex=2,col= alpha("#7CAE00",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

##### 
points(QuadraMidBivalveCG3[1,1],QuadraMidBivalveCG3[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidBivalveReference3[1,1],QuadraMidBivalveReference3[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.5,8.08, expression( "Richness"["T"] ))
text(0.5,10.08, expression( "Richness"["C"] ))

### Now plotting Control shannon and Simpson

points(QuadraMidBivalveReference3[10,1],QuadraMidBivalveReference3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidBivalveReference3[19,1],QuadraMidBivalveReference3[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(QuadraMidBivalveCG3[10,1],QuadraMidBivalveCG3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidBivalveCG3[19,1],QuadraMidBivalveCG3[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.45,2.8, expression( "Shannon"["C"] ))
text(2.3,1.5, expression( "Simpson"["C"] ))

text(1.45,4, expression( "Shannon"["T"] ))
text(2.4,3.5, expression( "Simpson"["T"] ))

####

#### Plot it

jpeg(filename = "QuadraBivalveOrderPlot1.jpeg", width = 22, height = 18, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(QuadraMidBivalveReference3[,1],QuadraMidBivalveReference3[,2], col=alpha("#F8766D",0.6),cex=2, pch=16,xlab="",ylab="",ylim = c(1,10.3))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(QuadraMidBivalveCG3[,1],QuadraMidBivalveCG3[,2],pch=16,cex=2,col= alpha("#7CAE00",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

##### 
points(QuadraMidBivalveCG3[1,1],QuadraMidBivalveCG3[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidBivalveReference3[1,1],QuadraMidBivalveReference3[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.45,8.08, expression( "Richness"["T"] ))
text(0.45,10.08, expression( "Richness"["C"] ))

### Now plotting Control shannon and Simpson

points(QuadraMidBivalveReference3[10,1],QuadraMidBivalveReference3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidBivalveReference3[19,1],QuadraMidBivalveReference3[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(QuadraMidBivalveCG3[10,1],QuadraMidBivalveCG3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(QuadraMidBivalveCG3[19,1],QuadraMidBivalveCG3[19,2],col= alpha("black",0.9),pch=16,cex=2)

#text(1.4,2.8, expression( "Shannon"["C"] ))

text(.68,3.1, expression( "Shannon"["C"] ))
text(1.4,4.1, expression( "Shannon"["T"] ))

text(2.3,1.6, expression( "Simpson"["C"] ))
text(2.3,3.5, expression( "Simpson"["T"] ))

dev.off()


#### Clam Garden Combined Profiles ####
CGMidBivalveData4=aggregate(.~SiteType+Site, CGMidBivalveSpeciesData3, mean)
CGMidBivalveData5=aggregate(.~SiteType, CGMidBivalveData4, mean)
names(CalvertCGMidBivalveData5)

MidBivalveCG= CGMidBivalveData5[1,]
MidBivalveReference= CGMidBivalveData5[2,]

### Removing zero's
MidBivalveCG2=cbind(MidBivalveCG[,7:11],MidBivalveCG[,13:16],MidBivalveCG[,18])
MidBivalveReference2=cbind(MidBivalveReference[,7:16],MidBivalveReference[,18])
MidBivalveReference2
write.csv(MidBivalveReference2,"MidBivalveReference2.csv")
MidBivalveReference2=read.csv("MidBivalveReference2.csv",header=T)

### So then using this equation
divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

MidBivalveReference3=divprof(MidBivalveReference2)

### Checking Adult Control Data
exp(diversity(MidBivalveReference2, index="shannon"))
###   4.9811
1 / (1-(diversity(MidBivalveReference2, index="simpson")))
###  3.3766

MidBivalveCG3=divprof(MidBivalveCG2)

### Checking Adult Treatment Data
exp(diversity(MidBivalveCG2, index="shannon"))
###   4.4781
1 / (1-(diversity(MidBivalveCG2, index="simpson")))
###  3.2256


##### Building the plot
plot(MidBivalveReference3[,1],MidBivalveReference3[,2], col=alpha("#F8766D",0.6),cex=2, pch=16,xlab="",ylab="",ylim = c(1,12.4))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(MidBivalveCG3[,1],MidBivalveCG3[,2],pch=16,cex=2,col= alpha("#7CAE00",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

##### 
points(MidBivalveCG3[1,1],MidBivalveCG3[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[1,1],MidBivalveReference3[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.5,8.08, expression( "Richness"["T"] ))
text(0.5,10.08, expression( "Richness"["C"] ))

### Now plotting Control shannon and Simpson

points(MidBivalveReference3[10,1],MidBivalveReference3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[19,1],MidBivalveReference3[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(MidBivalveCG3[10,1],MidBivalveCG3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveCG3[19,1],MidBivalveCG3[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.45,2.8, expression( "Shannon"["C"] ))
text(2.3,1.5, expression( "Simpson"["C"] ))

text(1.45,4, expression( "Shannon"["T"] ))
text(2.4,3.5, expression( "Simpson"["T"] ))

####

#### Plot it

jpeg(filename = "ClamGardensBivalveOrderPlot3.jpeg", width = 22, height = 18, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(MidBivalveReference3[,1],MidBivalveReference3[,2], col=alpha("#F8766D",0.6),cex=2, pch=16,xlab="",ylab="",ylim = c(1,10.3))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(MidBivalveCG3[,1],MidBivalveCG3[,2],pch=16,cex=2,col= alpha("#7CAE00",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

##### 
points(MidBivalveCG3[1,1],MidBivalveCG3[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[1,1],MidBivalveReference3[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.45,8.08, expression( "Richness"["T"] ))
text(0.45,10.08, expression( "Richness"["C"] ))

### Now plotting Control shannon and Simpson

points(MidBivalveReference3[10,1],MidBivalveReference3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveReference3[19,1],MidBivalveReference3[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(MidBivalveCG3[10,1],MidBivalveCG3[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidBivalveCG3[19,1],MidBivalveCG3[19,2],col= alpha("black",0.9),pch=16,cex=2)

#text(1.4,2.8, expression( "Shannon"["C"] ))

text(.68,3.1, expression( "Shannon"["C"] ))
text(1.4,4.1, expression( "Shannon"["T"] ))

text(2.3,1.6, expression( "Simpson"["C"] ))
text(2.3,3.5, expression( "Simpson"["T"] ))

dev.off()


### #### So Combined Bivalve Diversity Plots blurr trends
### Region species analysis seems likely


#### Epifauna
CGDataMidCG
CGDataMidCG %>% count(Site, Year)

################ Bivalves
names(CGDataMidCG)
CGMidEpifaunaData <- CGDataMidCG[,33:92]
CGMidDataStructure1=CGDataMidCG[,1:2]
CGMidDataStructure2=CGDataMidCG[,5:7]
CGMidDataStructure3=cbind(CGMidDataStructure1,CGMidDataStructure2)
CGMidEpifaunaSpeciesData2=cbind(CGMidDataStructure3,CGMidEpifaunaData)
str(CGMidEpifaunaSpeciesData2)

CGMidEpifaunaSpeciesData2
CGMidEpifaunaSpeciesData3=CGMidEpifaunaSpeciesData2[, colSums(CGMidEpifaunaSpeciesData2 != 0) > 0]
names(CGMidEpifaunaSpeciesData2)
names(CGMidEpifaunaSpeciesData3)
str(CGMidEpifaunaSpeciesData3)
str(CGMidEpifaunaSpeciesData2)

CGMidEpifaunaSpeciesData3

##### Diversity Profiles, Calvert and Quadra
################ Diversity Profile Curves
CGMidEpifaunaSpeciesData3$Area

CalvertCGMidEpifaunaData3 <- subset(CGMidEpifaunaSpeciesData3, CGMidEpifaunaSpeciesData3$Area == "Calvert")
QuadraCGMidEpifaunaData3 <- subset(CGMidEpifaunaSpeciesData3, CGMidEpifaunaSpeciesData3$Area == "Quadra")

options(digits=5)
###### Epifauna Diversity Profile Calvert  ####

str(CalvertCGMidEpifaunaData3)
CalvertCGMidEpifaunaData3$Area
CalvertCGMidEpifaunaData3$Site
CalvertCGMidEpifaunaData3$Quad
CalvertCGMidEpifaunaData3$Year

CalvertCGMidEpifaunaData3a=aggregate(.~SiteType+Year+Site, CalvertCGMidEpifaunaData3, mean)

CalvertCGMidEpifaunaData4=aggregate(.~SiteType+Site, CalvertCGMidEpifaunaData3a, mean)
CalvertCGMidEpifaunaData5=aggregate(.~SiteType, CalvertCGMidEpifaunaData4, mean)
names(CalvertCGMidEpifaunaData5)

MidEpifaunaCG= CalvertCGMidEpifaunaData5[1,]
MidEpifaunaReference= CalvertCGMidEpifaunaData5[2,]

### Removing zero's

MidEpifaunaCG2  =  MidEpifaunaCG[, colSums(MidEpifaunaCG != 0) > 0]
MidEpifaunaRef2  =  MidEpifaunaReference[, colSums(MidEpifaunaReference != 0) > 0]

names(MidEpifaunaCG2)
names(MidEpifaunaRef2)

MidEpifaunaCG3 = MidEpifaunaCG2[,6:33]
MidEpifaunaRef3 = MidEpifaunaRef2[,6:23]

### So then using this equation
divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

MidEpifaunaRef4=divprof(MidEpifaunaRef3)

### Checking Adult Control Data
exp(diversity(MidEpifaunaRef3, index="shannon"))
###   2.3386
1 / (1-(diversity(MidEpifaunaRef3, index="simpson")))
###  1.4812

MidEpifaunaCG4=divprof(MidEpifaunaCG3)

### Checking Adult Treatment Data
exp(diversity(MidEpifaunaCG3, index="shannon"))
###   3.8328
1 / (1-(diversity(MidEpifaunaCG3, index="simpson")))
###  2.7925


##### Building the plot
plot(MidEpifaunaCG4[,1],MidEpifaunaCG4[,2], col=alpha("#7CAE00",0.6),cex=2, pch=16,xlab="",ylab="",ylim = c(0,30))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(MidEpifaunaRef4[,1],MidEpifaunaRef4[,2],pch=16,cex=2,col= alpha("#F8766D",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

##### 
points(MidEpifaunaCG4[1,1],MidEpifaunaCG4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaRef4[1,1],MidEpifaunaRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.5,29.08, expression( "Richness"["T"] ))
text(0.5,19.08, expression( "Richness"["C"] ))

### Now plotting Control shannon and Simpson

points(MidEpifaunaRef4[10,1],MidEpifaunaRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaRef4[19,1],MidEpifaunaRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(MidEpifaunaCG4[10,1],MidEpifaunaCG4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaCG4[19,1],MidEpifaunaCG4[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.45,2.8, expression( "Shannon"["C"] ))
text(2.3,1.5, expression( "Simpson"["C"] ))

text(1.45,4, expression( "Shannon"["T"] ))
text(2.4,3.5, expression( "Simpson"["T"] ))

#### Plot it
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "CalvertEpifaunaOrderPlot2.jpeg", width = 22, height = 18, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(MidEpifaunaCG4[,1],MidEpifaunaCG4[,2], col=alpha("#7CAE00",0.6),cex=2, pch=16,xlab="",ylab="",ylim = c(0,30))
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

points(MidEpifaunaRef4[,1],MidEpifaunaRef4[,2],pch=16,cex=2,col= alpha("#F8766D",0.4))

points(MidEpifaunaCG4[1,1],MidEpifaunaCG4[1,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaRef4[1,1],MidEpifaunaRef4[1,2],col= alpha("black",0.9),pch=16,cex=2)

text(0.45,28, expression( "Richness"["T"] ))
text(0.45,18, expression( "Richness"["C"] ))

points(MidEpifaunaRef4[10,1],MidEpifaunaRef4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaRef4[19,1],MidEpifaunaRef4[19,2],col= alpha("black",0.9),pch=16,cex=2)

points(MidEpifaunaCG4[10,1],MidEpifaunaCG4[10,2],col= alpha("black",0.9),pch=16,cex=2)
points(MidEpifaunaCG4[19,1],MidEpifaunaCG4[19,2],col= alpha("black",0.9),pch=16,cex=2)

text(1.4,4.7, expression( "Shannon"["T"] ))
text(1.4,0.6, expression( "Shannon"["C"] ))

text(2.35,4, expression( "Simpson"["T"] ))
text(2.35,0.2, expression( "Simpson"["C"] ))

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")

#### Repeat with Quadra
####


#### For Manuscript: Diversity Profile Figure- Panel
## Diversity Profile Bivalves CG- Calvert, Quadra, SF Quadra, Baynes





### Calvert Epifauna nMDS and SIMPER

##  nMDS
CalvertCGMidEpifaunaData3$Quad

CalvertCGMidEpifaunaData4=aggregate(.~SiteType+Year+Site+Quad, CalvertCGMidEpifaunaData3, mean)

# CalvertCGMidEpifaunaData5=aggregate(.~SiteType+Site, CalvertCGMidEpifaunaData4, mean)

CalvertCGMidEpifaunaData4  =  CalvertCGMidEpifaunaData4[, colSums(CalvertCGMidEpifaunaData4 != 0) > 0]

names(CalvertCGMidEpifaunaData4)
CalvertCGMidEpifaunaStructure=CalvertCGMidEpifaunaData4[,1:4]
CalvertCGMidEpifaunaSpeciesData <- CalvertCGMidEpifaunaData4[,6:38]

CalvertCGMidEpifaunaSpeciesData2=((CalvertCGMidEpifaunaSpeciesData+1))
nMDSCommunityMidEpifauna<- metaMDS(CalvertCGMidEpifaunaSpeciesData2, distance = "bray", k = 2) 
nMDSCommunityMidEpifauna$stress  
stressplot(nMDSCommunityMidEpifauna)

names(CalvertCGMidEpifaunaStructure)
co=c("#F8766D", "#7CAE00")
shape=c(1:24) 

plot(nMDSCommunityMidEpifauna$points, col=co[CalvertCGMidEpifaunaStructure$SiteType],  pch = shape[CalvertCGMidEpifaunaStructure$SiteType], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunityMidEpifauna,groups=CalvertCGMidEpifaunaStructure$SiteType,draw="polygon",col=c("#F8766D", "#7CAE00"),label=T,alpha = 40)

#txt <- c("Clam Garden","Farm","Reference")
legend('topleft', txt , pch=c(10),col=co,cex=0.7, bty = "y")


####   ^ Come back to this ^



#### Simper Analysis Epifauna Clam Gardens
#### 

CalvertCGMidEpifaunaData3$Quad

CalvertCGMidEpifaunaData4=aggregate(.~SiteType+Year+Site+Quad, CalvertCGMidEpifaunaData3, mean)
str(CalvertCGMidEpifaunaData4)
names(CalvertCGMidEpifaunaData4)
CalvertCGMidEpifaunaData5  =  CalvertCGMidEpifaunaData4[, colSums(CalvertCGMidEpifaunaData4 != 0) > 0]
str(CalvertCGMidEpifaunaData5)
names(CalvertCGMidEpifaunaData5)

names(CalvertCGMidEpifaunaData5)
CalvertCGMidEpifaunaStructure=CalvertCGMidEpifaunaData5[,1:4]
CalvertCGMidEpifaunaSpeciesData <- CalvertCGMidEpifaunaData5[,6:38]
CalvertCGMidEpifaunaSpeciesData2=((CalvertCGMidEpifaunaSpeciesData+1))

###### Simper
SpeciesSIMPER<-(simper(CalvertCGMidEpifaunaSpeciesData2, CalvertCGMidEpifaunaStructure$SiteType, permutations=999))
summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$`Clam Garden_Reference`[,1]), na.rm = TRUE)

### so the communities are 68% dissimilar 
summary=summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response

summary2=summary$`Clam Garden_Reference`

write.csv(summary2,"CalvertEpifaunaSIMPERS2.csv")
EpifaunaSIMPERCalvert=read.csv("CalvertEpifaunaSIMPERSFinal3.csv",header=T)
EpifaunaSIMPERCalvert
names(EpifaunaSIMPERCalvert)
str(EpifaunaSIMPERCalvert)
EpifaunaSIMPERCalvert$Species <- factor(EpifaunaSIMPERCalvert$Species, levels = EpifaunaSIMPERCalvert$Species[order(EpifaunaSIMPERCalvert$Average.Dissim.....)])

col=c( "#F8766D", "#7CAE00" )
# base_size = 15
p=ggplot(EpifaunaSIMPERCalvert, aes(x = Species, y = Average.Dissim..... ,colour=Direction))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity ")+ ylim(-35,35)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=col)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)

AdultSimperPlot=p+coord_flip()
AdultSimperPlot

library(grid)

# Create a text
grob <- grobTree(textGrob("Dissimilarity 68.70%", x=0.005,  y=1.05, hjust=0,
                          gp=gpar(col="Black", fontsize=14)))
AdultSimperPlot + annotation_custom(grob)+geom_line(linetype = 2)+geom_vline(xintercept = 33.6,size=0.6)

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "EpifaunaSimperPlotV1.jpeg", width = 25, height = 18, units = "cm", pointsize = 15, res = 300)  

AdultSimperPlot + annotation_custom(grob)+geom_line(linetype = 2)+geom_vline(xintercept = 33.6,size=0.7)

dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")



##### Replot without three main species
EpifaunaSIMPERCalvert
EpifaunaSIMPERCalvert$Species
EpifaunaSIMPERCalvert2= EpifaunaSIMPERCalvert[4:33,]


p=ggplot(EpifaunaSIMPERCalvert2, aes(x = Species, y = Average.Dissim..... ,colour=Direction))+ 
  theme_classic() + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity ")+ ylim(-3.5,3.5)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12),axis.title.x  = element_text(size = 15)) +scale_color_manual(values=col)+
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text())+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)

AdultSimperPlot=p+coord_flip()
AdultSimperPlot

library(grid)

# Create a text
grob <- grobTree(textGrob("Dissimilarity 68.70%", x=0.005,  y=1.05, hjust=0,
                          gp=gpar(col="Black", fontsize=14)))
AdultSimperPlot + annotation_custom(grob)+geom_line(linetype = 2)+geom_vline(xintercept = 30.6,size=0.6)

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis/Results")
jpeg(filename = "EpifaunaSimperPlotV2Cut.jpeg", width = 25, height = 18, units = "cm", pointsize = 15, res = 300)  

AdultSimperPlot + annotation_custom(grob)+geom_line(linetype = 2)+geom_vline(xintercept = 30.6,size=0.7)

dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/CG SF /Epifauna Data/Analysis")





##

























#############   Old Code
#nMDS- Site Type, Site Type and Region, Site Type and Region and Year

SFandCGData %>% count(Site, Year, Zone)
SFandCGDataMid <- subset(SFandCGData, SFandCGData$Zone == "Mid")
SFandCGDataMid
SFandCGData2 %>% count(Site, Year)

### nMDS Plots 
SFandCGDataMid$Zone
names(SFandCGDataMid)
MidEpiSpeciesData <- SFandCGDataMid[,32:91]
str(MidEpiSpeciesData)
MidDataStructure1<-SFandCGDataMid[,1:6]
MidDataStructure1<-SFandCGDataMid[,1:2]
MidDataStructure2<-SFandCGDataMid[,4:6]
MidDataStructure=cbind(MidDataStructure1,MidDataStructure2)
str(MidDataStructure)

################ Bivalves
MidBivalveSpeciesData <- SFandCGDataMid[,240:253]
MidDataStructure1=SFandCGDataMid[,1:2]
MidDataStructure2=SFandCGDataMid[,4:6]
MidDataStructure3=cbind(MidDataStructure1,MidDataStructure2)
MidBivalveSpeciesData2=cbind(MidDataStructure3,MidBivalveSpeciesData)
str(MidBivalveSpeciesData2)

## Removing NAs
names(MidBivalveSpeciesData2)
str(MidBivalveSpeciesData2)
MidBivalveSpeciesData3=(MidBivalveSpeciesData2[complete.cases(MidBivalveSpeciesData2), ])
str(MidBivalveSpeciesData3)
names(MidBivalveSpeciesData3)
MidBivalveDataStructure=MidBivalveSpeciesData3[,1:5]
MidBivalveSpeciesData <- MidBivalveSpeciesData3[,6:19]
names(MidBivalveDataStructure)

### Transform and dummy
MidBivalveSpeciesData2=((MidBivalveSpeciesData+1)^(1/4))

nMDSCommunityMid<- metaMDS(MidBivalveSpeciesData2, distance = "bray", k = 2) 
nMDSCommunityMid$stress  ### Not great, maybe average at site
stressplot(nMDSCommunityMid)

plot(nMDSCommunityMid$points, col=co[MidBivalveDataStructure$SiteType],  pch = shape[MidBivalveDataStructure$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunityMid,groups=MidBivalveDataStructure$SiteType,draw="polygon",col=c("#F8766D", "#7CAE00","#0099ff"),label=F,alpha = 40)
txt <- c("Clam Garden","Farm","Reference")
legend('topleft', txt , pch=c(20),col=co,cex=1, bty = "y")
#### To much quadrat mess


######### average at Zone
MidBivalveSpeciesData4=aggregate(.~Area+SiteType+Site+Year,MidBivalveSpeciesData3,mean)

names(MidBivalveSpeciesData4)
MidBivalveDataStructureTest1=MidBivalveSpeciesData4[,1:4]
MidBivalveSpeciesDataTest1 <- MidBivalveSpeciesData4[,6:19]

MidBivalveSpeciesDataTest2=((MidBivalveSpeciesDataTest1+1)^(1/4))
nMDSCommunityMidAverage<- metaMDS(MidBivalveSpeciesDataTest2, distance = "bray", k = 3) 
nMDSCommunityMidAverage$stress  
stressplot(nMDSCommunityMidAverage)

names(MidBivalveDataStructure)
co=c("#F8766D", "#7CAE00","#0099ff")
shape=c(1:24) 

plot(nMDSCommunityMidAverage$points, col=co[MidBivalveDataStructureTest1$SiteType],  pch = shape[MidBivalveDataStructureTest1$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunityMidAverage,groups=MidBivalveDataStructureTest1$SiteType,draw="polygon",col=c("#F8766D", "#7CAE00","#0099ff"),label=T,alpha = 40)

#txt <- c("Clam Garden","Farm","Reference")
legend('topleft', txt , pch=c(10),col=co,cex=0.7, bty = "y")





######### Here, finish this through ############
#legend('topright', "Stress 0.07" ,cex=.8, bty = "n",x.intersp=-.5,y.intersp=-0.2)

plot(nMDSCommunityAdult$points, col=co[AdultDataStructure$Treatment],  pch = shape[AdultDataStructure$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunityAdult,groups=AdultDataStructure$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 40)


plot(nMDSCommunitySeedling$points, col=co[SeedlingDataStructure$Treatment],  pch = shape[SeedlingDataStructure$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunitySeedling,groups=SeedlingDataStructure$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 40)


### plotting
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis/Results")

jpeg(filename = "AdultnMDS.jpeg", width = 23, height = 13, units = "cm", pointsize = 15, quality = 300, res = 300)  
plot(nMDSCommunityAdult$points, col=co[AdultDataStructure$Treatment],  pch = shape[AdultDataStructure$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunityAdult,groups=AdultDataStructure$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 40)
dev.off()



##### That'll do

### Now make Simper Plots


### Old Code, messy

##### So let's look at Similarity Percentages to get specific values for the differences between treatments
### Simper time

#### Getting data sorted
MidBivalveSpeciesData4$Area

MidBivalveSpeciesData5 <- subset(MidBivalveSpeciesData4, MidBivalveSpeciesData4$Area == "Calvert")
MidBivalveSpeciesData6=aggregate(.~SiteType+Site, MidBivalveSpeciesData5, mean)

MidBivalveSpeciesData7=aggregate(.~SiteType, MidBivalveSpeciesData6, mean)
names(MidBivalveSpeciesData7)


MidBivalveSpeciesData5

MidBivalveSpeciesData8= MidBivalveSpeciesData7[,1]
MidBivalveSpeciesData9= MidBivalveSpeciesData7[,6:18]
MidBivalveSpeciesData10=cbind(MidBivalveSpeciesData8,MidBivalveSpeciesData9)

MidBivalveCG= MidBivalveSpeciesData10[1,]
MidBivalveReference= MidBivalveSpeciesData10[2,]


#### Option 1   Site average
MidBivalveSpeciesData6
str(MidBivalveSpeciesData6)
EnviromentData= data.frame((MidBivalveSpeciesData6[,1]))
colnames(EnviromentData)[1] <- "SiteType"

SpeciesData=cbind(MidBivalveSpeciesData6[,7:10], MidBivalveSpeciesData6[,12:16] , MidBivalveSpeciesData6[,18])
colnames(SpeciesData)[10] <- "Nutricola.tantilla.LargeQuad"


#### Option 2 Quad Values

MidBivalveSpeciesData3
MidCalvertBivalveSpeciesData3 <- subset(MidBivalveSpeciesData3, MidBivalveSpeciesData3$Area == "Calvert")
MidCalvertBivalveSpeciesData3
str(MidCalvertBivalveSpeciesData3)

EnviromentData= (MidCalvertBivalveSpeciesData3[,1:5])

SpeciesData=cbind(MidCalvertBivalveSpeciesData3[,7:10], MidCalvertBivalveSpeciesData3[,12:16] , MidCalvertBivalveSpeciesData3[,18])
colnames(SpeciesData)[10] <- "Nutricola.tantilla.LargeQuad"

##### Option 3 Sites and Years
str(MidBivalveSpeciesData5)
EnviromentData= data.frame((MidBivalveSpeciesData5[,2]))
colnames(EnviromentData)[1] <- "SiteType"
SpeciesData=cbind(MidBivalveSpeciesData5[,7:10], MidBivalveSpeciesData5[,12:16] , MidBivalveSpeciesData5[,18])
colnames(SpeciesData)[10] <- "Nutricola.tantilla.LargeQuad"


###### Simper
SpeciesSIMPER<-(simper(SpeciesData, EnviromentData$SiteType, permutations=999))
summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$`Clam Garden_Reference`[,1]), na.rm = TRUE)
SpeciesData

### so the communities are 49% dissimilar 
summary=summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response

summary2=summary$`Clam Garden_Reference`

write.csv(summary2,"CalvertBivalveSIMPERS2.csv")
BivalveSIMPERS=read.csv("CalvertBivalveSIMPERSFinal.csv",header=T)

BivalveSIMPERS

BivalveSIMPERS$Species <- factor(BivalveSIMPERS$Species, levels = BivalveSIMPERS$Species[order(BivalveSIMPERS$Average.Dissim.....)])

col=c( "#F8766D", "#7CAE00" )

p=ggplot(BivalveSIMPERS, aes(x = Species, y = Average.Dissim..... ,colour=Direction))+ 
  theme_classic(base_size = 15) + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity ")+ ylim(-25,25)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 15)) +scale_color_manual(values=col)+
  theme(axis.text.y = element_text(colour = "black", size = 20), axis.title.y = element_text(size = 15))+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)


AdultSimperPlot=p+coord_flip()


AdultSimperPlot+ geom_text(x=3, y=-30, label="Scatter plot")

library(grid)

# Create a text
grob <- grobTree(textGrob("Dissimilarity 49.87%", x=0.005,  y=0.989, hjust=0,
                          gp=gpar(col="Black", fontsize=14)))
AdultSimperPlot + annotation_custom(grob)+geom_line(linetype = 2)+geom_vline(xintercept = 10.33,size=0.6)
                                               

jpeg(filename = "BivalvesSimperPlotV1.jpeg", width = 25, height = 18, units = "cm", pointsize = 15, res = 300)  

AdultSimperPlot + annotation_custom(grob)+geom_line(linetype = 2)+geom_vline(xintercept = 10.33,size=0.6)

dev.off()






# Plot
AdultSimperPlot + annotation_custom(grob)


# Solution 1
sp2 + geom_text(x=3, y=20, label="Scatter plot")
# Solution 2
sp2 + annotate(geom="text", x=3, y=30, label="Scatter plot",
               color="red")


###





































##### Old code #####


##### Ok so the aim is to get the best look at the community differences across the different treatments and sites
## concern being low species diversity, however, this may be an advantage as it will allow for more meaningful plotting to support community composition work
## i.e. 200 species communties can't really be plotted, thus community statistics, but a 6 species community might get the best of both worlds
## lets find out


data=read.csv("ADULTDATA.csv",header=T)
data
data2=read.csv("ADULTDATA2.csv",header=T) ## transposed verson of your data set
data2

data # let's check the data 
names(data)
### Note, I added a Unique ID column
str(data)
levels(data$Plot) ## switch to factor so we can consider it as levels
data$Plot=as.factor(data$Plot)
str(data)
data2$Plot=as.factor(data2$Plot)
str(data2)

library(ggplot2)
library(data.table)
library(tidyverse)
library(dplyr)
library(vegan)
library(cowplot)

### so first things first, let's look at what we are trying to analyze
## we'll use the jitter function to see all the raw data, which will let us know if we have outliers that are an issue
### Across site trends
ggplot(data2, aes(x = Species , y = Density, fill = Treatment)) + theme_classic(base_size = 15) + ggtitle("Flora Density") + geom_jitter(position = position_dodge(width = 1.0))+ ylab("Density")+ xlab("Species") + geom_boxplot(size = 1, position = position_dodge(width = 1.0)) + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, colour = "black", size = 15)) + theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) #+ guides(colour = FALSE, fill = FALSE)

### strong indication that you'll have reduced have reduced densities across your sites/treatments relative to your treatment
### Dots are sites, but lets take a quick look in more detail at site differences
## we can also explore species diversity, but you mentioned you explored diversity + metrics. Diversity might be tough given low diversity across sites 
# i.e. 3sp vs 4sp might might not seem like a lot but the composition could be totally different

RC <- subset(data2, data2$Species == "RC")
SP <- subset(data2, data2$Species == "SP")
SS <- subset(data2, data2$Species == "SS")
WH <- subset(data2, data2$Species == "WH")
YC <- subset(data2, data2$Species == "YC")
YW <- subset(data2, data2$Species == "YW")

ggplot(RC, aes(x = Treatment , y = Density, fill = Site, color=Treatment)) + theme_classic(base_size = 15) + ggtitle("RC Density") + geom_jitter(position = position_dodge(width = 1.0))+ ylab("Density")+ xlab("Species") +scale_color_manual(values=c("#000000", "#8B0000"))+ geom_boxplot(size = 1, position = position_dodge(width = 1.0)) + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) + theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + guides(colour = FALSE, fill = FALSE)
ggplot(SP, aes(x = Treatment , y = Density, fill = Site, color=Treatment)) + theme_classic(base_size = 15) + ggtitle("SP Density") + geom_jitter(position = position_dodge(width = 1.0))+ ylab("Density")+ xlab("Species") +scale_color_manual(values=c("#000000", "#8B0000"))+ geom_boxplot(size = 1, position = position_dodge(width = 1.0)) + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) + theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + guides(colour = FALSE, fill = FALSE)
ggplot(SS, aes(x = Treatment , y = Density, fill = Site, color=Treatment)) + theme_classic(base_size = 15) + ggtitle("SS Density") + geom_jitter(position = position_dodge(width = 1.0))+ ylab("Density")+ xlab("Species") +scale_color_manual(values=c("#000000", "#8B0000"))+ geom_boxplot(size = 1, position = position_dodge(width = 1.0)) + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) + theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + guides(colour = FALSE, fill = FALSE)
ggplot(WH, aes(x = Treatment , y = Density, fill = Site, color=Treatment)) + theme_classic(base_size = 15) + ggtitle("WH Density") + geom_jitter(position = position_dodge(width = 1.0))+ ylab("Density")+ xlab("Species") +scale_color_manual(values=c("#000000", "#8B0000"))+ geom_boxplot(size = 1, position = position_dodge(width = 1.0)) + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) + theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + guides(colour = FALSE, fill = FALSE)
ggplot(YC, aes(x = Treatment , y = Density, fill = Site, color=Treatment)) + theme_classic(base_size = 15) + ggtitle("YC Density") + geom_jitter(position = position_dodge(width = 1.0))+ ylab("Density")+ xlab("Species") +scale_color_manual(values=c("#000000", "#8B0000"))+ geom_boxplot(size = 1, position = position_dodge(width = 1.0)) + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) + theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + guides(colour = FALSE, fill = FALSE)
ggplot(YW, aes(x = Treatment , y = Density, fill = Site, color=Treatment)) + theme_classic(base_size = 15) + ggtitle("YW Density") + geom_jitter(position = position_dodge(width = 1.0))+ ylab("Density")+ xlab("Species") +scale_color_manual(values=c("#000000", "#8B0000"))+ geom_boxplot(size = 1, position = position_dodge(width = 1.0)) + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) + theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + guides(colour = FALSE, fill = FALSE)


#### k so lets consider a few different levels of the data first, we'll use PERMANOVA to look at species composition via distance matrices  
data
names(data)
SpeciesData= data[,5:10]
SpeciesData

EnviromentData= data[,1:4]
EnviromentData
names(EnviromentData)

adonis2(SpeciesData ~ Treatment+Site+Plot, data = EnviromentData,permutations = 999,method = "bray",contr.unordered = "contr.sum")
##plot doesn't seem to be doing much. We'll consider it solo to make sure
adonis2(SpeciesData ~ Treatment+Site, data = EnviromentData,permutations = 999,method = "bray",contr.unordered = "contr.sum")
## could keep or remove site
adonis2(SpeciesData ~ Treatment, data = EnviromentData,permutations = 999,method = "bray",contr.unordered = "contr.sum")
### So treatment level community differences

### We'll look at the community trends 

#Plot ordination so that points are coloured and shaped according to treatment
nMDSCommunity<- metaMDS(SpeciesData, distance = "bray", k = 3) # Stress of ~ 0.07 is good, it means your control/treatment fit well
stressplot(nMDSCommunity)

## euclidean

co=c("#F8766D", "#7CAE00")
shape=c(10:19)  # This will allow us to show the plots within each site as the same shapes

plot(nMDSCommunity$points, col=co[EnviromentData$Treatment],  pch = shape[EnviromentData$Site], 
     cex=1.2, main="Plant Communities",  xlab = "axis 1", ylab = "axis 2")
### few different ways we can look at the assocaited of these plots
ordispider(nMDSCommunity, groups = EnviromentData$Treatment,  label = F,spiders = c("centroid")) ## Default is ,spiders = c("centroid")
ordiellipse(nMDSCommunity, groups = EnviromentData$Treatment,  label = F,col=c("#F8766D", "#7CAE00"),conf=0.5)
ordihull(nMDSCommunity,groups=EnviromentData$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 30)

### As a quick illistrate of how your much the communties within your sites differ, we can look at this 
plot(nMDSCommunity$points, col=co[EnviromentData$Treatment],  pch = shape[EnviromentData$Site], 
     cex=1.2, main="Plant Communities",  xlab = "axis 1", ylab = "axis 2")
ordihull(nMDSCommunity,groups=EnviromentData$Site,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 30)
### the lack of overland across sites indicates distinct communites at the site level

### For presenting I'd go with ordihull so something like
plot(nMDSCommunity$points, col=co[EnviromentData$Treatment],  pch = shape[EnviromentData$Site], 
     cex=1.2, main="Plant Communities",  xlab = "axis 1", ylab = "axis 2")
ordihull(nMDSCommunity,groups=EnviromentData$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 40)
#Add legend
txt <- c("Control","Treatment")
legend('topleft', txt , pch=c(20),col=co,cex=1, bty = "y")
## if you want to label/legend your sites, you can speak to sites that overlap in community structure, and those that don't

##### So let's look at Similarity Percentages to get specific values for the differences between treatments
### Simper time
SpeciesData

SpeciesSIMPER<-(simper(SpeciesData, EnviromentData$Treatment, permutations=999))
summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response
SpeciesSIMPERSum<-summary(SpeciesSIMPER)
sum(as.numeric(SpeciesSIMPERSum$C_T[,1]), na.rm = TRUE)
SpeciesData
### so the communities are 34% dissimilar 

summary=summary(SpeciesSIMPER) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response
summary

SpeciesSIMPERSum$C_T[,1]

list(as.numeric(SpeciesSIMPERSum$C_T[,1]), na.rm = TRUE)


######## Rerunning simper data but averaging at quadrat level
data

dataSIMPER2= data[,2:10]
dataSIMPER2
dataSIMPER2$Plot=as.factor(dataSIMPER2$Plot)
str(dataSIMPER2)
names(dataSIMPER2)
dataSIMPER3=aggregate(.~Treatment+Site,dataSIMPER2,mean)
dataSIMPER3


EnviromentDatadataSIMPER3=dataSIMPER3[,1:3]
EnviromentDatadataSIMPER3
SpeciesdataSIMPER3=dataSIMPER3[,4:9]
SpeciesdataSIMPER3

SpeciesSIMPER2<-(simper(SpeciesdataSIMPER3, EnviromentDatadataSIMPER3$Treatment, permutations=999))
summary(SpeciesSIMPER2) ## These are the dis sim across the species, this is where the 'interesting' information is. Specific species response
SpeciesSIMPERSum<-summary(SpeciesSIMPER2)
sum(as.numeric(SpeciesSIMPERSum$C_T[,1]), na.rm = TRUE)









##### Final Analysis March 3rd, 2020

#Objectives
## Replot nMDS plots in High res, as A&B Plots
## Plot diversity metrics along a continious measurement of q (increasing abundance effects)
## Plot SIMPER Data

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis")

AdultData=read.csv("ADULTDATA.csv",header=T)
AdultData
str(AdultData)
AdultData$Plot=as.factor(AdultData$Plot)


RegenataDiversity=read.csv("REGENDATADIVERSITY.csv",header=T)
RegenataDiversity
str(RegenataDiversity)
RegenataDiversity$Plot=as.factor(RegenataDiversity$Plot)
RegenataDiversity

####### Adult nMDS  
AdultData
AdultDataSpecies <- AdultData[,4:9]
AdultDataSpecies

AdultDataStructure<- AdultData[,1:3]
AdultDataStructure

nMDSCommunityAdult<- metaMDS(AdultDataSpecies, distance = "bray", k = 3) 
nMDSCommunityAdult$stress
stressplot(nMDSCommunityAdult)

AdultDataStructure$Treatment
AdultDataStructure$Site

plot(nMDSCommunityAdult$points, col=co[AdultDataStructure$Treatment],  pch = shape[AdultDataStructure$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunityAdult,groups=AdultDataStructure$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 40)
# txt <- c("Control","Treatment")
# legend('topleft', txt , pch=c(20),col=co,cex=1, bty = "y")

# nMDSCommunityAdult$stress
#legend('topright', "Stress 0.07" ,cex=.8, bty = "n",x.intersp=-.5,y.intersp=-0.2)
#legend('topright', "Stress 0.07" ,cex=.8, bty = "n",x.intersp=-.5,y.intersp=-0.2)

## Seedling nMSD
RegenataDiversity
names(RegenataDiversity)

SeedlingDataSpecies <- RegenataDiversity[,5:9]
SeedlingDataSpecies

SeedlingDataStructure<- RegenataDiversity[,1:4]
SeedlingDataStructure

nMDSCommunitySeedling<- metaMDS(SeedlingDataSpecies, distance = "bray", k = 3) 
stressplot(nMDSCommunitySeedling)

SeedlingDataStructure$Treatment
SeedlingDataStructure$Site
### For presenting I'd go with ordihull so something like
plot(nMDSCommunitySeedling$points, col=co[SeedlingDataStructure$Treatment],  pch = shape[SeedlingDataStructure$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunitySeedling,groups=SeedlingDataStructure$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 40)
#Add legend

# nMDSCommunitySeedling$stress
legend('topleft', "Stress 0.03" ,cex=1, bty = "n",x.intersp=-.8,y.intersp=-0.2)

# txt <- c("Control","Treatment")
# legend('topleft', txt , pch=c(20),col=co,cex=1, bty = "y")
## if you want to label/legend your sites, you can speak to sites that overlap in community structure, and those that don't

plot(nMDSCommunityAdult$points, col=co[AdultDataStructure$Treatment],  pch = shape[AdultDataStructure$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunityAdult,groups=AdultDataStructure$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 40)


plot(nMDSCommunitySeedling$points, col=co[SeedlingDataStructure$Treatment],  pch = shape[SeedlingDataStructure$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunitySeedling,groups=SeedlingDataStructure$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 40)


### plotting
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis/Results")

jpeg(filename = "AdultnMDS.jpeg", width = 23, height = 13, units = "cm", pointsize = 15, quality = 300, res = 300)  
plot(nMDSCommunityAdult$points, col=co[AdultDataStructure$Treatment],  pch = shape[AdultDataStructure$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunityAdult,groups=AdultDataStructure$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 40)
dev.off()

### 
jpeg(filename = "SeedlingnMDS.jpeg", width = 23, height = 13, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(nMDSCommunitySeedling$points, col=co[SeedlingDataStructure$Treatment],  pch = shape[SeedlingDataStructure$Site], 
     cex=1.2, main="",  xlab = "Axis 1", ylab = "Axis 2")
ordihull(nMDSCommunitySeedling,groups=SeedlingDataStructure$Treatment,draw="polygon",col=c("#F8766D", "#7CAE00"),label=F,alpha = 40)
dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis")

### i'm going to export these to add silhouette to them. Well combined, label etc at the same time




################ Diversity Profile Curves
#### k so now to the hills diversity plots

#### Kyle's data, Adults
AdultData
names(AdultData)
AdultData2= AdultData[,1:9]
names(AdultData2)

AdultData3=aggregate(.~Treatment+Site, AdultData2, mean)
AdultData3

AdultData4= AdultData3[,1:2]
AdultData5= AdultData3[,4:9]
AdultData6=cbind(AdultData4,AdultData5)
AdultData6

AdultData7=aggregate(.~Treatment, AdultData6, mean)
AdultData7

AdultData8= AdultData7[,1]
AdultData9= AdultData7[,3:8]
AdultData10=cbind(AdultData8,AdultData9)
AdultData10


AdultControlData2= AdultData10[1,]
AdultTreatmentData2= AdultData10[2,]

AdultControlData2
AdultTreatmentData2

# write.csv(AdultTreatmentData2,"AdultTreatmentData2.csv")
AdultTreatmentData3=read.csv("AdultTreatmentData2.csv",header=T)
AdultTreatmentData3

AdultTreatmentData4= AdultTreatmentData3[,3:8]
AdultTreatmentData4
## remove zero's

AdultTreatmentData5 <- AdultTreatmentData4[ -c(2, 5) ]

AdultControlData3=AdultControlData2[,2:7]
AdultControlData3



### So then using this equation

divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

AdultControlData4=divprof(AdultControlData3)
AdultControlData4



### Checking Adult Control Data
exp(diversity(AdultControlData3, index="shannon"))
###   3.188747
1 / (1-(diversity(AdultControlData3, index="simpson")))
###  2.645768

AdultTreatmentData6=divprof(AdultTreatmentData5)
AdultTreatmentData6

### Checking Adult Treatment Data
exp(diversity(AdultTreatmentData5, index="shannon"))
###   2.626311
1 / (1-(diversity(AdultTreatmentData5, index="simpson")))
###  2.331009



##### Building the plot
plot(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("#7CAE00",0.6),cex=2, pch=16,xlab="",ylab="")

title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(AdultTreatmentData6[,1],AdultTreatmentData6[,2], col=alpha("black",0.3),cex=2, pch=21)
points(AdultControlData4[,1],AdultControlData4[,2],pch=16,cex=2,col= alpha("#F8766D",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

##### 
points(AdultTreatmentData6[1,1],AdultTreatmentData6[1,2],col= alpha("black",1),pch=21,cex=2)

# text(1,AdultControlData4[1,2],labels=c("Richness"[2]))
# text(1,AdultControlData4[1,2],labels=c("Richness"))

text(0.61,6.01, expression( "Richness"["T&C"] ))

### Now plotting Control shannon and Simpson

# Where q=1, diversity is Shannon diversity
AdultControlData4
points(AdultControlData4[10,1],AdultControlData4[10,2],col= alpha("black",1),pch=21,cex=2)
points(AdultControlData4[10,1],AdultControlData4[10,2],col= alpha("#F8766D",0.4),pch=16,cex=2)

text(1.5,3.26, expression( "Shannon"["C"] ))

# Where q=2, diversity is Simpson diversity
AdultControlData4[19,1]
AdultControlData4[19,2]
points(AdultControlData4[19,1],AdultControlData4[19,2],col= alpha("black",1),pch=21,cex=2)
points(AdultControlData4[19,1],AdultControlData4[19,2],col= alpha("#F8766D",0.4),pch=16,cex=2)

text(2.45,2.78, expression( "Simpson"["C"] ))

### Now plotting Treatment shannon and Simpson

# Where q=1, diversity is Shannon diversity
AdultTreatmentData6

points(AdultTreatmentData6[10,1],AdultTreatmentData6[10,2],col= alpha("black",1),pch=21,cex=2)
points(AdultTreatmentData6[10,1],AdultTreatmentData6[10,2],col= alpha("#7CAE00",0.4),pch=16,cex=2)
text(1.39,2.411708, expression( "Shannon"["T"] ))

# Where q=2, diversity is Simpson diversity
AdultTreatmentData6[19,1]
AdultTreatmentData6[19,2]
points(AdultTreatmentData6[19,1],AdultTreatmentData6[19,2],col= alpha("black",1),pch=21,cex=2)
points(AdultTreatmentData6[19,1],AdultTreatmentData6[19,2],col= alpha("#7CAE00",0.4),pch=16,cex=2)

text(2.41,2.145, expression( "Simpson"["C"] ))


##### Seedling Data
RegenataDiversity

names(RegenataDiversity)
RegenataDiversity2= RegenataDiversity[,1:9]
names(RegenataDiversity2)
RegenataDiversity$Plot=as.factor(RegenataDiversity$Plot)

RegenataDiversity3=aggregate(.~Treatment+Site, RegenataDiversity2, mean)
RegenataDiversity3

RegenataDiversity4= RegenataDiversity3[,1:2]
RegenataDiversity5= RegenataDiversity3[,5:9]
RegenataDiversity6=cbind(RegenataDiversity4,RegenataDiversity5)
RegenataDiversity6

RegenataDiversity7=aggregate(.~Treatment, RegenataDiversity6, mean)
RegenataDiversity7

RegenataDiversity8= RegenataDiversity7[,1]
RegenataDiversity9= RegenataDiversity7[,3:7]
RegenataDiversity10=cbind(RegenataDiversity8,RegenataDiversity9)
RegenataDiversity10


SeedlingControlData2= RegenataDiversity10[1,]
SeedlingTreatmentData2= RegenataDiversity10[2,]


#write.csv(SeedlingTreatmentData2,"SeedlingTreatmentData2.csv")
SeedlingTreatmentData3=read.csv("SeedlingTreatmentData2.csv",header=T)
SeedlingTreatmentData3

SeedlingTreatmentData4= SeedlingTreatmentData3[,3:7]
SeedlingTreatmentData4
SeedlingTreatmentData4=SeedlingTreatmentData4[ -c(4)]

SeedlingControlData3=SeedlingControlData2[,2:6]
SeedlingControlData3
SeedlingControlData3=SeedlingControlData3[ -c(5)]


### So then using this equation

divprof=function(community) {
  cbind(
    seq(0,5,by=0.11),
    unlist(lapply(seq(0,5,by=0.11),function(q) sum(apply(community,1,function(x) 
      (x/sum(x))^q))^(1/(1-q))))) }

SeedlingControlData4=divprof(SeedlingControlData3)
SeedlingControlData4

### Checking Adult Control Data
exp(diversity(SeedlingControlData3, index="shannon"))
###   2.118837
1 / (1-(diversity(SeedlingControlData3, index="simpson")))
###  1.874762

SeedlingTreatmentData5=divprof(SeedlingTreatmentData4)
SeedlingTreatmentData5

### Checking Adult Treatment Data
exp(diversity(SeedlingTreatmentData4, index="shannon"))
###   2.614551
1 / (1-(diversity(SeedlingTreatmentData4, index="simpson")))
###  2.39157



##### Building the plot      #### !! Revised, see below under exporting ###
#### Treatment
plot(SeedlingControlData4[,1],SeedlingControlData4[,2],pch=16,cex=2,col= alpha("#F8766D",0.6),xlab="",ylab="")

points(SeedlingTreatmentData5[,1],SeedlingTreatmentData5[,2], col=alpha("#7CAE00",0.6),cex=2, pch=16)

title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

# points(SeedlingTreatmentData5[,1],SeedlingTreatmentData5[,2], col=alpha("black",0.3),cex=2, pch=21)
#### points(SeedlingControlData4[,1],SeedlingControlData4[,2],pch=16,cex=2,col= alpha("#F8766D",0.4))
# points(AdultControlData4[,1],AdultControlData4[,2],pch=21,cex=2,col= alpha("black",0.3))

##### 
points(SeedlingTreatmentData5[1,1],SeedlingTreatmentData5[1,2],col= alpha("black",1),pch=21,cex=2)

# text(1,SeedlingControlData4[1,2],labels=c("Richness"[2]))
# text(1,SeedlingControlData4[1,2],labels=c("Richness"))
SeedlingTreatmentData5[1,1]
text(0.56,4.99, expression( "Richness"["T&C"] ))

### Now plotting Control shannon and Simpson
# Where q=1, diversity is Shannon diversity

points(SeedlingControlData4[10,1],SeedlingControlData4[10,2],col= alpha("black",1),pch=21,cex=2)
points(SeedlingControlData4[10,1],SeedlingControlData4[10,2],col= alpha("#F8766D",0.4),pch=16,cex=2)

text(1.455,2.18, expression( "Shannon"["C"] ))

# Where q=2, diversity is Simpson diversity
SeedlingControlData4[19,1]
SeedlingControlData4[19,2]
points(SeedlingControlData4[19,1],SeedlingControlData4[19,2],col= alpha("black",1),pch=21,cex=2)
points(SeedlingControlData4[19,1],SeedlingControlData4[19,2],col= alpha("#F8766D",0.4),pch=16,cex=2)

text(2.43,1.96, expression( "Simpson"["C"] ))

### Now plotting Treatment shannon and Simpson

# Where q=1, diversity is Shannon diversity
SeedlingTreatmentData5

points(SeedlingTreatmentData5[10,1],SeedlingTreatmentData5[10,2],col= alpha("black",1),pch=21,cex=2)
points(SeedlingTreatmentData5[10,1],SeedlingTreatmentData5[10,2],col= alpha("#7CAE00",0.4),pch=16,cex=2)
text(1.45,2.7, expression( "Shannon"["T"] ))

# Where q=2, diversity is Simpson diversity
SeedlingTreatmentData5[19,1]
SeedlingTreatmentData5[19,2]
points(SeedlingTreatmentData5[19,1],SeedlingTreatmentData5[19,2],col= alpha("black",1),pch=21,cex=2)
points(SeedlingTreatmentData5[19,1],SeedlingTreatmentData5[19,2],col= alpha("#7CAE00",0.4),pch=16,cex=2)

text(2.425,2.5, expression( "Simpson"["C"] ))


### End



######## Here
### Revise plots via refitting points
### Export
### Revise Figure
### Revise Text


### Exporting
# Adults
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis/Results/Diversity Profile")

jpeg(filename = "AdultOrderPlot.jpeg", width = 20, height = 20, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(AdultControlData4[,1],AdultControlData4[,2], col=alpha("#F8766D",0.4),cex=2, pch=16,xlab="",ylab="",ylim=c(2, 6))

title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

points(AdultTreatmentData6[,1],AdultTreatmentData6[,2],pch=16,cex=2,col= alpha("#7CAE00",0.6))

points(AdultTreatmentData6[1,1],AdultTreatmentData6[1,2],col= alpha("black",1),pch=16,cex=2)
points(AdultControlData4[1,1],AdultControlData4[1,2],col= alpha("black",1),pch=16,cex=2)

text(0.73,4.01, expression( "Richness"["T"] ))
text(0.73,6.01, expression( "Richness"["C"] ))

points(AdultControlData4[10,1],AdultControlData4[10,2],col= alpha("black",1),pch=16,cex=2)

text(1.65,3.35, expression( "Shannon"["C"] ))

points(AdultControlData4[19,1],AdultControlData4[19,2],col= alpha("black",1),pch=16,cex=2)

text(2.6,2.85, expression( "Simpson"["C"] ))

points(AdultTreatmentData6[10,1],AdultTreatmentData6[10,2],col= alpha("black",1),pch=16,cex=2)
text(1.2,2.4, expression( "Shannon"["T"] ))

points(AdultTreatmentData6[19,1],AdultTreatmentData6[19,2],col= alpha("black",1),pch=16,cex=2)

text(2.25,2.1, expression( "Simpson"["T"] ))

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis")


# Seedlings
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis/Results/Diversity Profile")

jpeg(filename = "SeedlingOrderPlot.jpeg", width = 20, height = 20, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot(SeedlingControlData4[,1],SeedlingControlData4[,2],pch=16,cex=2,col= alpha("#F8766D",0.6),xlab="",ylab="")
points(SeedlingTreatmentData5[,1],SeedlingTreatmentData5[,2], col=alpha("#7CAE00",0.6),cex=2, pch=16)
title(xlab = list("Order q", cex = 1.3))
title(ylab = list("Diversity", cex = 1.3))

points(SeedlingTreatmentData5[1,1],SeedlingTreatmentData5[1,2],col= alpha("black",1),pch=16,cex=2)

text(0.56,4.99, expression( "Richness"["T&C"] ))

points(SeedlingControlData4[10,1],SeedlingControlData4[10,2],col= alpha("black",1),pch=16,cex=2)

text(1.455,2.18, expression( "Shannon"["C"] ))

points(SeedlingControlData4[19,1],SeedlingControlData4[19,2],col= alpha("black",1),pch=16,cex=2)

text(2.43,1.96, expression( "Simpson"["C"] ))

points(SeedlingTreatmentData5[10,1],SeedlingTreatmentData5[10,2],col= alpha("black",1),pch=16,cex=2)
text(1.45,2.7, expression( "Shannon"["T"] ))

points(SeedlingTreatmentData5[19,1],SeedlingTreatmentData5[19,2],col= alpha("black",1),pch=16,cex=2)
text(2.425,2.5, expression( "Simpson"["C"] ))

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis")






######### Simper plots

###   colors used    #F8766D     #7CAE00   
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis")
list.files()
SIMPERAdultData=read.csv("SIMPERAdultData.csv",header=T)
SIMPERAdultData
names(SIMPERAdultData)

SIMPERAdultData$Species <- factor(SIMPERAdultData$Species, levels = SIMPERAdultData$Species[order(SIMPERAdultData$Average.Dissim.....)])

col=c( "#F8766D", "#7CAE00" )

p=ggplot(SIMPERAdultData, aes(x = Species, y = Average.Dissim..... ,colour=Direction))+ 
  theme_classic(base_size = 15) + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity ")+ ylim(-25,25)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 15)) +scale_color_manual(values=col)+
  theme(axis.text.y = element_text(colour = "black", size = 20), axis.title.y = element_text(size = 15))+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)


AdultSimperPlot=p+coord_flip()

AdultSimperPlot



##### seedling SIMPER

list.files()
SIMPERSeedlingData=read.csv("SIMPERSeedlingData.csv",header=T)
SIMPERSeedlingData
names(SIMPERSeedlingData)

SIMPERSeedlingData$Species <- factor(SIMPERSeedlingData$Species, levels = SIMPERSeedlingData$Species[order(SIMPERSeedlingData$Average.Dissim.....)])

p2=ggplot(SIMPERSeedlingData, aes(x = Species, y = Average.Dissim..... ,colour=Direction))+ 
  theme_classic(base_size = 15) + ggtitle("") +geom_point(size = 6)+ylab("Average Dissimilarity ")+ ylim(-25,25)+
  xlab("") + theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5, colour = "black", size = 15)) +scale_color_manual(values=col)+
  theme(axis.text.y = element_text(colour = "black", size = 20), axis.title.y = element_text(size = 15))+ guides(colour = FALSE, fill = FALSE)+
  geom_hline(yintercept = 0,colour = "grey",size = 1.25,linetype=2)

SIMPERSeedlingPlot=p2+coord_flip()

SIMPERSeedlingPlot



plot_grid(AdultSimperPlot,SIMPERSeedlingPlot, labels=c('A)', 'B)'),nrow = 1,ncol=2)


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis/Results/SIMPER Plot")
jpeg(filename = "SimperPlotV1.jpeg", width = 35, height = 15, units = "cm", pointsize = 15, quality = 300, res = 300)  

plot_grid(AdultSimperPlot,SIMPERSeedlingPlot, labels=c('A)', 'B)'),nrow = 1,ncol=2)

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis")




####### Tree structure plot

list.files()
AdultStructure=read.csv("ADULTSTRUCTUREDATA.csv",header=T)
AdultStructure
names(AdultStructure)
str(AdultStructure)
AdultStructure$Plot=as.factor(AdultStructure$Plot)
levels(AdultStructure$Species)
AdultStructure$Height

# Aim is to how the target species height increase within quirtile groups
### redcedar, Sitka spruce and western hemlock
# RC      SS       WH 

### so few different ways to plot this
# two routes come to mind
# plot all on the same axis but this might get a bit messy
# Plot each species as it's own section, then the control and treatment within the 4Q blocks
# lets go with each on their own

## selecting target species
SpeciesHeights=(subset(AdultStructure, Species %in% c("RC","SS","WH")))
SpeciesHeights

RCHeights=(subset(AdultStructure, Species %in% c("RC")))
SSHeights=(subset(AdultStructure, Species %in% c("SS")))
WHHeights=(subset(AdultStructure, Species %in% c("WH")))

### few clear groupings but not consistent
hist(RCHeights$Height,breaks=40)
hist(SSHeights$Height,breaks=40)
hist(WHHeights$Height,breaks=40)

quantile(RCHeights$Height)
quantile(SSHeights$Height)
quantile(WHHeights$Height)

## break into 4 
RCHeightsQ4= mutate(RCHeights, quantile_rank = ntile(RCHeights$Height,4))
RCHeightsQ4$quantile_rank= as.factor(RCHeightsQ4$quantile_rank)
levels(RCHeightsQ4$quantile_rank)

SSHeightsQ4= mutate(SSHeights, quantile_rank = ntile(SSHeights$Height,4))
SSHeightsQ4$quantile_rank= as.factor(SSHeightsQ4$quantile_rank)
levels(SSHeightsQ4$quantile_rank)

WHHeightsQ4= mutate(WHHeights, quantile_rank = ntile(WHHeights$Height,4))
WHHeightsQ4$quantile_rank= as.factor(WHHeightsQ4$quantile_rank)
levels(WHHeightsQ4$quantile_rank)


## agg data
## mean
RCHeightsQ4mean=aggregate(Height~Treatment+Site+quantile_rank,RCHeightsQ4,mean)
names(RCHeightsQ4mean)[names(RCHeightsQ4mean) == 'Height'] <- 'mean'

SSHeightsQ4Mean=aggregate(Height~Treatment+Site+quantile_rank,SSHeightsQ4,mean)
names(SSHeightsQ4Mean)[names(SSHeightsQ4Mean) == 'Height'] <- 'mean'

WHHeightsQ4Mean=aggregate(Height~Treatment+Site+quantile_rank,WHHeightsQ4,mean)
names(WHHeightsQ4Mean)[names(WHHeightsQ4Mean) == 'Height'] <- 'mean'

###### So was going to go with SD but the 2 plots means SD isn't possible without pooling plots which isn't ideal
## also can't determine it at the treatment level as there aren't always 3 in each height cat.
## sd  
# RCHeightsQ4sd=aggregate(Height~Treatment+Site+quantile_rank,RCHeightsQ4,sd)
# names(RCHeightsQ4sd)[names(RCHeightsQ4sd) == 'Height'] <- 'sd'

# SSHeightsQ4sd=aggregate(Height~Treatment+Site+quantile_rank,SSHeightsQ4,sd)
# names(SSHeightsQ4sd)[names(SSHeightsQ4sd) == 'Height'] <- 'sd'

# WHHeightsQ4sd=aggregate(Height~Treatment+Site+quantile_rank,WHHeightsQ4,sd)
# names(WHHeightsQ4sd)[names(WHHeightsQ4sd) == 'Height'] <- 'sd'

## lets plot it
library(ggplot2)
names(RCHeightsQ4Mean)
col=c( "#F8766D", "#7CAE00" )

RCHeightsPlot=ggplot(RCHeightsQ4mean, aes(x= quantile_rank, y= mean, fill=Treatment) )+
  theme_classic(base_size = 15) +scale_fill_manual(values=c( "#F8766D", "#7CAE00" ))+ geom_boxplot(size = .3, position = position_dodge(width = 1.0)) +ylim(0,30)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Quantile")+ ylab("Height (m)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + 
  guides(colour = FALSE, fill = FALSE) 
RCHeightsPlot

SSHeightsPlot=ggplot(SSHeightsQ4Mean, aes(x= quantile_rank, y= mean, fill=Treatment) )+
  theme_classic(base_size = 15) +scale_fill_manual(values=c( "#F8766D", "#7CAE00" ))+ geom_boxplot(size = .3, position = position_dodge(width = 1.0)) +ylim(0,30)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Quantile")+ ylab("Height (m)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + 
  guides(colour = FALSE, fill = FALSE) 
SSHeightsPlot

WHHeightsPlot=ggplot(WHHeightsQ4Mean, aes(x= quantile_rank, y= mean, fill=Treatment) )+
  theme_classic(base_size = 15) +scale_fill_manual(values=c( "#F8766D", "#7CAE00" ))+ geom_boxplot(size = .3, position = position_dodge(width = 1.0)) +ylim(0,20)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Quantile")+ ylab("Height (m)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + 
  guides(colour = FALSE, fill = FALSE) 
WHHeightsPlot

citation('vegan')

library(cowplot)
plot_grid(RCHeightsPlot,WHHeightsPlot,SSHeightsPlot, labels=c('A)', 'B)','C)'),nrow = 3,ncol=1)


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis/Results")

jpeg(filename = "HeightsPlot.jpeg", width = 15, height = 30, units = "cm", pointsize = 15, quality = 600, res = 600)  

plot_grid(RCHeightsPlot,WHHeightsPlot,SSHeightsPlot, labels=c('A)', 'B)','C)'),nrow = 3,ncol=1)


dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis")





#### Replotting 

#### Option one 3+4

library(dplyr)
RCHeightsQ4= mutate(RCHeights, quantile_rank = ntile(RCHeights$Height,4))
RCHeightsQ4$quantile_rank[RCHeightsQ4$quantile_rank == "4"] <- "3"
RCHeightsQ4$quantile_rank= as.factor(RCHeightsQ4$quantile_rank)
levels(RCHeightsQ4$quantile_rank)

SSHeightsQ4= mutate(SSHeights, quantile_rank = ntile(SSHeights$Height,4))
SSHeightsQ4$quantile_rank[SSHeightsQ4$quantile_rank == "4"] <- "3"
SSHeightsQ4$quantile_rank= as.factor(SSHeightsQ4$quantile_rank)
levels(SSHeightsQ4$quantile_rank)

WHHeightsQ4= mutate(WHHeights, quantile_rank = ntile(WHHeights$Height,4))
WHHeightsQ4$quantile_rank[WHHeightsQ4$quantile_rank == "4"] <- "3"
WHHeightsQ4$quantile_rank= as.factor(WHHeightsQ4$quantile_rank)
levels(WHHeightsQ4$quantile_rank)


## agg data
## mean
RCHeightsQ4mean=aggregate(Height~Treatment+Site+quantile_rank,RCHeightsQ4,mean)
names(RCHeightsQ4mean)[names(RCHeightsQ4mean) == 'Height'] <- 'mean'

SSHeightsQ4Mean=aggregate(Height~Treatment+Site+quantile_rank,SSHeightsQ4,mean)
names(SSHeightsQ4Mean)[names(SSHeightsQ4Mean) == 'Height'] <- 'mean'

WHHeightsQ4Mean=aggregate(Height~Treatment+Site+quantile_rank,WHHeightsQ4,mean)
names(WHHeightsQ4Mean)[names(WHHeightsQ4Mean) == 'Height'] <- 'mean'


#### Sorting out average increase for manuscript statement

aggregate(mean~quantile_rank+Treatment,SSHeightsQ4Mean,mean)
aggregate(mean~quantile_rank+Treatment,RCHeightsQ4mean,mean)
aggregate(mean~quantile_rank+Treatment,WHHeightsQ4Mean,mean)

24.77037-19.96667
11.347927-10.837634
22.382733-17.502911

24.77-19.97
22.38-17.5

###### So was going to go with SD but the 2 plots means SD isn't possible without pooling plots which isn't ideal
## also can't determine it at the treatment level as there aren't always 3 in each height cat.
## sd  
# RCHeightsQ4sd=aggregate(Height~Treatment+Site+quantile_rank,RCHeightsQ4,sd)
# names(RCHeightsQ4sd)[names(RCHeightsQ4sd) == 'Height'] <- 'sd'

# SSHeightsQ4sd=aggregate(Height~Treatment+Site+quantile_rank,SSHeightsQ4,sd)
# names(SSHeightsQ4sd)[names(SSHeightsQ4sd) == 'Height'] <- 'sd'

# WHHeightsQ4sd=aggregate(Height~Treatment+Site+quantile_rank,WHHeightsQ4,sd)
# names(WHHeightsQ4sd)[names(WHHeightsQ4sd) == 'Height'] <- 'sd'

## lets plot it
library(ggplot2)
names(RCHeightsQ4Mean)
col=c( "#F8766D", "#7CAE00" )

RCHeightsPlot3plu4=ggplot(RCHeightsQ4mean, aes(x= quantile_rank, y= mean, fill=Treatment) )+
  theme_classic(base_size = 15) +scale_fill_manual(values=c( "#F8766D", "#7CAE00" ))+ geom_boxplot(size = .3, position = position_dodge(width = 1.0)) +ylim(0,30)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Quantile")+ ylab("Height (m)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + 
  guides(colour = FALSE, fill = FALSE) 
RCHeightsPlot3plu4

SSHeightsPlot3plu4=ggplot(SSHeightsQ4Mean, aes(x= quantile_rank, y= mean, fill=Treatment) )+
  theme_classic(base_size = 15) +scale_fill_manual(values=c( "#F8766D", "#7CAE00" ))+ geom_boxplot(size = .3, position = position_dodge(width = 1.0)) +ylim(0,30)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Quantile")+ ylab("Height (m)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + 
  guides(colour = FALSE, fill = FALSE) 
SSHeightsPlot3plu4

WHHeightsPlot3plu4=ggplot(WHHeightsQ4Mean, aes(x= quantile_rank, y= mean, fill=Treatment) )+
  theme_classic(base_size = 15) +scale_fill_manual(values=c( "#F8766D", "#7CAE00" ))+ geom_boxplot(size = .3, position = position_dodge(width = 1.0)) +ylim(0,20)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Quantile")+ ylab("Height (m)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + 
  guides(colour = FALSE, fill = FALSE) 
WHHeightsPlot3plu4




library(cowplot)
plot_grid(RCHeightsPlot,WHHeightsPlot,SSHeightsPlot, labels=c('A)', 'B)','C)'),nrow = 3,ncol=1)


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis/Results")

jpeg(filename = "HeightsPlotRevised3plus4.jpeg", width = 15, height = 30, units = "cm", pointsize = 15, quality = 600, res = 600)  

plot_grid(RCHeightsPlot3plu4,WHHeightsPlot3plu4,SSHeightsPlot3plu4, labels=c('A)', 'B)','C)'),nrow = 3,ncol=1)

dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis")



### Replotting 

#### Option one 3 even

library(dplyr)
RCHeightsQ3= mutate(RCHeights, quantile_rank = ntile(RCHeights$Height,3))
RCHeightsQ3$quantile_rank= as.factor(RCHeightsQ3$quantile_rank)
levels(RCHeightsQ3$quantile_rank)

SSHeightsQ3= mutate(SSHeights, quantile_rank = ntile(SSHeights$Height,3))
SSHeightsQ3$quantile_rank= as.factor(SSHeightsQ3$quantile_rank)
levels(SSHeightsQ3$quantile_rank)

WHHeightsQ3= mutate(WHHeights, quantile_rank = ntile(WHHeights$Height,3))
WHHeightsQ3$quantile_rank= as.factor(WHHeightsQ3$quantile_rank)
levels(WHHeightsQ3$quantile_rank)

write.csv(WHHeightsQ3,"WHHeightsQ3.csv")
write.csv(RCHeightsQ3,"RCHeightsQ3.csv")
write.csv(SSHeightsQ3,"SSHeightsQ3.csv")


## agg data
## mean
RCHeightsQ3mean=aggregate(Height~Treatment+Site+quantile_rank,RCHeightsQ3,mean)
names(RCHeightsQ3mean)[names(RCHeightsQ3mean) == 'Height'] <- 'mean'

SSHeightsQ3Mean=aggregate(Height~Treatment+Site+quantile_rank,SSHeightsQ3,mean)
names(SSHeightsQ3Mean)[names(SSHeightsQ3Mean) == 'Height'] <- 'mean'

WHHeightsQ3Mean=aggregate(Height~Treatment+Site+quantile_rank,WHHeightsQ3,mean)
names(WHHeightsQ3Mean)[names(WHHeightsQ3Mean) == 'Height'] <- 'mean'

###### So was going to go with SD but the 2 plots means SD isn't possible without pooling plots which isn't ideal
## also can't determine it at the treatment level as there aren't always 3 in each height cat.
## sd  
# RCHeightsQ4sd=aggregate(Height~Treatment+Site+quantile_rank,RCHeightsQ4,sd)
# names(RCHeightsQ4sd)[names(RCHeightsQ4sd) == 'Height'] <- 'sd'

# SSHeightsQ4sd=aggregate(Height~Treatment+Site+quantile_rank,SSHeightsQ4,sd)
# names(SSHeightsQ4sd)[names(SSHeightsQ4sd) == 'Height'] <- 'sd'

# WHHeightsQ4sd=aggregate(Height~Treatment+Site+quantile_rank,WHHeightsQ4,sd)
# names(WHHeightsQ4sd)[names(WHHeightsQ4sd) == 'Height'] <- 'sd'

## lets plot it
library(ggplot2)
names(RCHeightsQ4Mean)
col=c( "#F8766D", "#7CAE00" )

RCHeightsPlot3=ggplot(RCHeightsQ3mean, aes(x= quantile_rank, y= mean, fill=Treatment) )+
  theme_classic(base_size = 15) +scale_fill_manual(values=c( "#F8766D", "#7CAE00" ))+ geom_boxplot(size = .3, position = position_dodge(width = 1.0)) +ylim(0,30)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Quantile")+ ylab("Height (m)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + 
  guides(colour = FALSE, fill = FALSE) 
RCHeightsPlot3

SSHeightsPlot3=ggplot(SSHeightsQ3Mean, aes(x= quantile_rank, y= mean, fill=Treatment) )+
  theme_classic(base_size = 15) +scale_fill_manual(values=c( "#F8766D", "#7CAE00" ))+ geom_boxplot(size = .3, position = position_dodge(width = 1.0)) +ylim(0,30)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Quantile")+ ylab("Height (m)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + 
  guides(colour = FALSE, fill = FALSE) 
SSHeightsPlot3

WHHeightsPlot3=ggplot(WHHeightsQ3Mean, aes(x= quantile_rank, y= mean, fill=Treatment) )+
  theme_classic(base_size = 15) +scale_fill_manual(values=c( "#F8766D", "#7CAE00" ))+ geom_boxplot(size = .3, position = position_dodge(width = 1.0)) +ylim(0,20)+
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("Quantile")+ ylab("Height (m)")+ theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + 
  guides(colour = FALSE, fill = FALSE) 
WHHeightsPlot3




library(cowplot)
plot_grid(RCHeightsPlot3,WHHeightsPlot3,SSHeightsPlot3, labels=c('A)', 'B)','C)'),nrow = 3,ncol=1)


setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis/Results")

jpeg(filename = "HeightsPlotRevised3.jpeg", width = 15, height = 30, units = "cm", pointsize = 15, quality = 600, res = 600)  

plot_grid(RCHeightsPlot3,WHHeightsPlot3,SSHeightsPlot3, labels=c('A)', 'B)','C)'),nrow = 3,ncol=1)

dev.off()

setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis")


####### Nurse Plot

list.files()
nurse=read.csv("NURSELOGGRAPH.csv", header=T)

## original
nurseplot<-ggplot(nurse, aes(x=Location, y=Count, fill=factor(Treatment))) + geom_boxplot() + labs(fill = "Location", x="") + theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ylab(bquote('Average seedling count within (  '* '25' ~ m^2*')'))
nurseplot

## Revised to match theme
WHHeightsPlot
nurseplot2=ggplot(nurse, aes(x= Location, y= Count, fill=Treatment) )+
  theme_classic(base_size = 15) +scale_fill_manual(values=c( "#F8766D", "#7CAE00" ))+ geom_boxplot(size = .3, position = position_dodge(width = 1.0)) +
  theme(plot.title = element_text(size = rel(1))) + theme(legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 15)) +
  xlab("")+ ylab(bquote('Seedling Density ('* 'Ind./25' ~ m^2*')'))+theme(axis.text.y = element_text(colour = "black", size = 17), axis.title.y = element_text(size = 15)) + 
  guides(colour = FALSE, fill = FALSE) 

nurseplot2



setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis/Results")
jpeg(filename = "NursePlotV2.jpeg", width = 20, height = 10, units = "cm", pointsize = 15, quality = 600, res = 600)  

nurseplot2

dev.off()
setwd("~/Documents/Documents/school/UVic, VIU, Hakai Dudas, Juanes/Mauscripts/Kyle Schang, Trant Lab/Manuscript/Analysis")










### could supplement this with a distribution plot as well to see what it would look like
## will be redundant so maybe not worth it but worth knowing about

SSHeightsControl=(subset(SSHeights, Treatment %in% c("Control")))
SSHeightsTreatment=(subset(SSHeights, Treatment %in% c("Habitation")))

hist(SSHeightsControl$Height,breaks=10)
hist(SSHeightsTreatment$Height,breaks=10)

##
RCHeightsControl=(subset(RCHeights, Treatment %in% c("Control")))
RCHeightsTreatment=(subset(RCHeights, Treatment %in% c("Habitation")))

hist(RCHeightsControl$Height,breaks=10)
hist(RCHeightsTreatment$Height,breaks=10)
#
RCHeightsControl=(subset(RCHeights, Treatment %in% c("Control")))
RCHeightsTreatment=(subset(RCHeights, Treatment %in% c("Habitation")))

hist(RCHeightsControl$Height,breaks=10)
hist(RCHeightsTreatment$Height,breaks=10)

#

WHHeightsControl=(subset(WHHeights, Treatment %in% c("Control")))
WHHeightsTreatment=(subset(WHHeights, Treatment %in% c("Habitation")))

hist(WHHeightsControl$Height,breaks=10)
hist(WHHeightsTreatment$Height,breaks=10)

## ya lets add this
names(RCHeights)

ggplot(WHHeights, aes(x=Height,color=Treatment)) + geom_density(size = 1.2)+ theme_classic()+
  guides(colour = FALSE, fill = FALSE) +scale_color_manual(values=c( "#F8766D", "#7CAE00" ))

ggplot(RCHeights, aes(x=Height,color=Treatment)) + geom_density(size = 1.2)+ theme_classic()+
  guides(colour = FALSE, fill = FALSE) +scale_color_manual(values=c( "#F8766D", "#7CAE00" ))

ggplot(SSHeights, aes(x=Height,color=Treatment)) + geom_density(size = 1.2)+ theme_classic()+
  guides(colour = FALSE, fill = FALSE) +scale_color_manual(values=c( "#F8766D", "#7CAE00" ))

### I like these less but good to know

