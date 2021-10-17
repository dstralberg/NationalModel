library(raster)
library(sf)
library(maptools)
library(dplyr)
library(data.table)
library(reshape2)
library(raster)
library(geosphere)
library(sp)
library(rgdal)
require(rgeos)
rasterOptions(chunksize = 1e+04, maxmemory = 1e+06)

w <-"R:/Boreal/RshProjs/PopnStatus/Density_Popn/NationalModelsV3/"

###### Projections and basemaps ########

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LAEA <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
DD <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


####### Buffered BCR shapefiles ########

bcrs <- c("bcr4_100km.shp","bcr5_100km.shp","bcr9_100km.shp","bcr10_100km.shp","bcr11_100km.shp","bcr12_100km.shp","bcr13_100km.shp","bcr14_100km.shp")
bcrs2 <- c("bcr60_100km.shp","bcr61_100km.shp","bcr70_100km.shp","bcr71_100km.shp","bcr80_100km.shp","bcr81_100km.shp", "bcr82_100km.shp","bcr83_100km.shp")
bcrsus <- c("bcrus2_100km.shp", "bcrus4_100km.shp", "bcrus5_100km.shp", "bcrus14_100km.shp", "bcrus13_100km.shp", "bcrus12_100km.shp", "bcrus11_100km.shp", "bcrus10_100km.shp")

####### Fixed layers ########

#Lambert Azimuthal Equal Area
#Human footprint# (100-m)
usroaddist <- raster(paste0(w,"transportation/USroaddistLAZEA.tif"))
canroaddist <- raster(paste0(w,"transportation/canroaddist1.tif"))

#Adaptwest baseline climate variables
cur <- paste0(w,"climate/baseline19812010/")
setwd(cur)
clim <- list.files(cur, pattern =".asc$")
curclim<-stack(clim)

#Topography (100-m)
TPI <- raster(paste0(w,"topoedaphic/TPI.tif"))
TRI <- raster(paste0(w,"topoedaphic/TRI.tif"))
slope <- raster(paste0(w,"topoedaphic/slope.tif"))
roughness <- raster(paste0(w,"topoedaphic/roughness.tif"))
topo <- stack(TPI,TRI,slope,roughness)

#Lambert Conformal Conic
#kNN biomass estimates (250-m)
biomass2001 <- brick(paste0(w,"biomass/bs2001_1km.grd"))
biomass2011 <- brick(paste0(w,"biomass/bs2001_1km.grd"))

landscape2001 <- brick(paste0(w,"biomass/bs2001_750_1km.grd"))
names(landscape2001)[1:4] <- gsub("v1","v1.1",names(landscape2001)[1:4])
landscape2011 <- brick(paste0(w,"biomass/bs2011_750_1km.grd"))
names(landscape2011)[1:4] <- gsub("v1","v1.1",names(landscape2011)[1:4])

#Topography (100-m)
lf <- raster(paste0(w,"topoedaphic/lf_lcc1.tif"))

####### Functions to extract annual layers ########

#Annual climate variables (1-km, DD, North America)
annclim <- function(surveys, yr) {
  cmi <- raster(paste0(w,"Climate/NRCAN/Annual/",yr,"/cmi60_sum.asc"))
  clim <- stack(lapply(grep("bio",list.files(paste0(w,"Climate/NRCAN/Annual/",yr),pattern=".tif$",full.names=TRUE),value=TRUE),FUN=raster))
  cmi <- resample(cmi,clim[[1]])
  clim <- addLayer(clim,cmi)
  data <- cbind(surveys,raster::extract(clim,surveys))
  return(data)
}

#Landcover Priority 1. ABoVE, LANDSAT, 30-m, western boreal/arctic region)
b1 <- brick(paste0(w,"landcover/NASA_ABoVE/1_to_17_mosaic.tif"))
ABOVE <- function(surveys, yr) {
  data <- cbind(surveys,raster::extract(b1[[yr-1983]],surveys))
  names(data)[ncol(data)] <- "ABOVELC"
  return(data)
}

#Landcover Priority 2. IGBP, MODIS, 250-m, LCC, Canada)
CanadaMODISLC <- function (surveys, yr) {
  nalcx <- raster(paste0(w,"landcover/LCTS_2000-2011/LCTS19c",yr,".tif"))
  data <- cbind(surveys,raster::extract(nalcx,surveys))
  names(data)[ncol(data)]<-"LCTS250m"
  return(data)
}

#Landcover Priority 3. NALCMS, MODIS, 500-m, LCC, North America)
MODIS500 <- function(surveys,yr){
  modx <- raster(paste0(w,"landcover/LCCMODIS500/LC500_",yr,".tif"))  
  data <- cbind(surveys,raster::extract(modx,surveys))
  names(data)[ncol(data)]<-"LC500m"
  return(data)
}

#Crosswalks for landcover layers
IGBP <- read.csv(paste0(w,"IGBP_XWalk.csv"))
NALCMS <- read.csv(paste0(w,"NALCMS_XWalk.csv"))
ABOVE <- read.csv(paste0(w,"ABOVE_XWalk.csv"))

