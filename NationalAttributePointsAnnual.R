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
library(MODISTools)
library(tidyr)

w <-"R:/Boreal/RshProjs/PopnStatus/Density_Popn/NationalModelsV3/"

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LAEA <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
DD <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
TM10 <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +datum=NAD83 +units=m +no_defs")

PKEYBBS <- read.csv("D:/BAM/BAMData/BBS/bbs_output/bbs_visit.csv")
PKEYBBS <- separate(PKEYBBS,col=12,into=c("PROJ","SITE","STATION","YEAR","VISIT"),sep=":",remove=FALSE)
PKEYBBS$SS <- paste(PKEYBBS$PROJ,PKEYBBS$SITE,PKEYBBS$STATION,sep=":")
PKEYBBS <- PKEYBBS[,c(4,12,18)]
names(PKEYBBS) <- c("YEAR","PKEY","SS")

SSBBS <- read.csv("D:/BAM/BAMData/BBS/bbs_output/bbs_XY.csv")[,3:5]
names(SSBBS) <- c("SS", "longitude", "latitude")

PKEYWT <- read.csv("D:/BAM/BAMData/WT-PKEY-2021-06-25.csv")
PKEYWT <- PKEYWT[,c(2,5:6)]
SSWT <- read.csv("D:/BAM/BAMData/WT-XY-2021-06-23.csv")

PKEYcombo <- read.csv("D:/BAM/BAMData/BAMVisitv6Feb2021.csv")
PKEYcombo <- separate(PKEYcombo,col=5,into=c("PROJ","SITE","STATION","YEAR","VISIT"),sep=":",remove=FALSE)
PKEYcombo$SS <- paste(PKEYcombo$PROJ,PKEYcombo$SITE,PKEYcombo$STATION,sep=":")
#PKEYcombo <- as.data.frame(cbind(PKEYcombo$SS,PKEYcombo$PKEY_V6,PKEYcombo$survey_year))
PKEYcombo <- as.data.frame(cbind(PKEYcombo$location_name_4,PKEYcombo$PKEY_V6,PKEYcombo$survey_year))
names(PKEYcombo) <- c("SS","PKEY","YEAR")
PKEYcombo <- rbind(PKEYcombo,PKEYBBS,PKEYWT)
PKEYcombo$YEAR <- as.numeric(PKEYcombo$YEAR)

SScombo <- read.csv("D:/BAM/BAMData/BAMXYv6Feb2021.csv")
# SScombo <- SScombo[,c(4:6)]
SScombo <- SScombo[,c(3,5:6)]
names(SScombo) <- c("SS", "latitude", "longitude")
SScombo <- rbind(SScombo,SSBBS)
SScombo <- rbind(SScombo,SSWT)

SSBBS_DD <- SSBBS
coordinates(SSBBS_DD) <- c("longitude", "latitude")
proj4string(SSBBS_DD) <- DD

SScombo_DD <- SScombo
SScombo_DD <- na.omit(SScombo_DD)
coordinates(SScombo_DD) <- c("longitude", "latitude")
proj4string(SScombo_DD) <- DD

SScombo_LCC <- spTransform(SScombo_DD,LCC)
SScomboLCC <- as.data.frame(SScombo_LCC)
SScombo_LAEA <- spTransform(SScombo_DD,LAEA)
SScomboLAEA <- as.data.frame(SScombo_LAEA)
SScombo_TM10 <- spTransform(SScombo_DD,TM10)
SScomboTM10 <- as.data.frame(SScombo_TM10)

#Prepare visit data
PKEYcombo_DD <- left_join(PKEYcombo, SScombo, by="SS")
PKEYcombo_DD <- na.omit(PKEYcombo_DD)
PKEYcomboDD <- as.data.frame(PKEYcombo_DD)
coordinates(PKEYcombo_DD) <- c("longitude", "latitude")
proj4string(PKEYcombo_DD) <- DD

PKEYcombo_TM10 <- left_join(PKEYcombo, SScomboTM10, by="SS")
PKEYcombo_TM10 <- na.omit(PKEYcombo_TM10)
PKEYcomboTM10 <- as.data.frame(PKEYcombo_TM10)
coordinates(PKEYcombo_TM10) <- c("longitude", "latitude")
proj4string(PKEYcombo_TM10) <- TM10

PKEYcombo_LCC <- left_join(PKEYcombo, SScomboLCC, by="SS")
PKEYcombo_LCC <- na.omit(PKEYcombo_LCC)
PKEYcomboLCC <- as.data.frame(PKEYcombo_LCC)
coordinates(PKEYcombo_LCC) <- c("longitude", "latitude")
proj4string(PKEYcombo_LCC) <- LCC

PKEYcombo_LAEA <- left_join(PKEYcombo, SScomboLAEA, by="SS")
PKEYcombo_LAEA <- na.omit(PKEYcombo_LAEA)
PKEYcomboLAEA <- as.data.frame(PKEYcombo_LAEA)
coordinates(PKEYcombo_LAEA) <- c("longitude", "latitude")
proj4string(PKEYcombo_LAEA) <- LAEA


####### Fixed layers ########

# Lambert Azimuthal Equal Area
#Human footprint#
usroaddist <- raster("E:/GIS/transportation/USroaddistLAZEA.tif")
canroaddist <- raster("E:/GIS/transportation/CanadaRoads2019_lite/canroaddist1.tif")

#Landform (100-m)
TPI <- raster("E:/GIS/topoedaphic/TPI.tif")
TRI <- raster("E:/GIS/topoedaphic/TRI.tif")
slope <- raster("E:/GIS/topoedaphic/slope.tif")
roughness <- raster("E:/GIS/topoedaphic/roughness.tif")
topo <- stack(TPI,TRI,slope,roughness)

SScombo_LAEA <- cbind(SScombo_LAEA,raster::extract(topo,SScombo_LAEA))
SScombo_LAEA <- cbind(SScombo_LAEA,raster::extract(canroaddist,SScombo_LAEA))
names(SScombo_LAEA)[ncol(SScombo_LAEA)] <- "CanRoadDist"
SScombo_LAEA <- cbind(SScombo_LAEA,raster::extract(usroaddist,SScombo_LAEA))
names(SScombo_LAEA)[ncol(SScombo_LAEA)] <- "USRoadDist"
countries <- shapefile("E:/GIS/basemaps/NorthAmericaLazea.shp")
SScombo_LAEA <- cbind(SScombo_LAEA,raster::extract(countries,SScombo_LAEA))
SScombo_LAEA$RoadDist <- pmin(SScombo_LAEA$CanRoadDist,SScombo_LAEA$USRoadDist, na.rm=TRUE)
SScombo_LAEA <- SScombo_LAEA[,c(1:8,13,22)]
names(SScombo_LAEA)[9] <- "Country"

write.csv(SScombo_LAEA,file=paste0(w,"SS_TopoRoad.csv"),row.names=FALSE)

lf <- raster("E:/GIS/topoedaphic/NA_Topo/lf_lcc1.tif")
SScombo_LCC <- cbind(SScombo_LCC,raster::extract(lf,SScombo_LCC))
names(SScombo_LCC)[ncol(SScombo_LCC)] <- "lf"

#Adaptwest baseline climate variables
cur <- "E:/Climate/CMIP5/baseline19812010/"
setwd(cur)
clim <- list.files(cur, pattern =".asc$")
curclim<-stack(clim)
SScombo_LCC <- cbind(SScombo_LCC,raster::extract(curclim,SScombo_LCC))

write.csv(SScombo_LCC,file=paste0(w,"SS_LandformClimate.csv"),row.names=FALSE)

#YearTier1 <- c(2001:2011)
#YearTier2 <- c(1991:2000,2012:2020)

####### Biomass layers ######

#kNN biomass layers, 2001 (250-m)
b2001 <- list.files("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2001/",pattern="tif$")
setwd("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2001/")
bs2001 <- stack(raster(b2001[1]))
for (i in 2:length(b2001)) {bs2001 <- addLayer(bs2001, raster(b2001[i]))}
names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
SScombo_LCC2001 <- cbind(SScombo_LCC,raster::extract(bs2001,SScombo_LCC))

bs2001_Gauss750 <- brick("G:/Boreal/NationalModelsV2/bs2001_750.grd")
names(bs2001_Gauss750) <- gsub("SpeciesGroups","Landsc750",names(bs2001_Gauss750))
names(bs2001_Gauss750) <- gsub("Species","Landsc750",names(bs2001_Gauss750))
names(bs2001_Gauss750) <- gsub("Structure","Landsc750",names(bs2001_Gauss750))
names(bs2001_Gauss750) <- gsub("LandCover","Landsc750",names(bs2001_Gauss750))
SScombo_LCC2001 <- cbind(SScombo_LCC2001,raster::extract(bs2001_Gauss750,SScombo_LCC2001))

#kNN biomass layers, 2011 (250-m)
b2011 <- list.files("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/",pattern="tif$")
setwd("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) {bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
SScombo_LCC2011 <- cbind(SScombo_LCC,raster::extract(bs2011,SScombo_LCC))

bs2011_Gauss750 <- brick("G:/Boreal/NationalModelsV2/bs2011_750.grd")
names(bs2011_Gauss750) <- gsub("SpeciesGroups","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Species","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Structure","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("LandCover","Landsc750",names(bs2011_Gauss750))
SScombo_LCC2011 <- cbind(SScombo_LCC2011,raster::extract(bs2011_Gauss750,SScombo_LCC2011))

PKEYcombo_2001 <- left_join(PKEYcombo[PKEYcombo$YEAR<2006,], as.data.frame(SScombo_LCC2001), by="SS")
PKEYcombo_2011 <- left_join(PKEYcombo[PKEYcombo$YEAR>2005,], as.data.frame(SScombo_LCC2011), by="SS")
PKEYcombo_dat <- rbind(PKEYcombo_2001,PKEYcombo_2011)
write.csv(PKEYcombo_dat,file=paste0(w,"PKEYcombo_VegBiomass.csv"),row.names=FALSE)


####### Annual layers ######

#McKenney annual climate variables
annclim <- function(surveys, yr) {
  cmi <- raster(paste0("E:/Climate/NRCAN/Annual/",yr,"/cmi60_sum.asc"))
  clim <- stack(lapply(grep("bio",list.files(paste0("E:/Climate/NRCAN/Annual/",yr),pattern=".tif$",full.names=TRUE),value=TRUE),FUN=raster))
  cmi <- resample(cmi,clim[[1]])
  clim <- addLayer(clim,cmi)
  data <- cbind(surveys,raster::extract(clim,surveys))
  return(data)
}

CanadaMODISLC <- function (surveys, yr) {
  nalcx <- raster(paste0("E:/GIS/landcover/LCTS_2000-2011/LCTS19c",yr,".tif"))
  data <- cbind(surveys,raster::extract(nalcx,surveys))
  names(data)[ncol(data)]<-"LCTS250m"
  return(data)
}

MODIS500 <- function(surveys,yr){
  modx <- raster(paste0("E:/GIS/landcover/LCCMODIS500/merge/LC500_",yr,".tif"))  
  data <- cbind(surveys,raster::extract(modx,surveys))
  names(data)[ncol(data)]<-"LC500m"
  return(data)
}

# MODIS500 <- function (surveys, yr) {
#   df <- surveys
#   df$site_name = df$SS
#   df <- df[,c(6,5,4)]
#   names(df) <- c("site_name","lon","lat")
#   lcdata <- mt_batch_subset(df = df,
#                             product = "MCD12Q1",
#                             band = "LC_Type1",
#                             start = paste0(yr,"-01-01"),
#                             end = paste0(yr,"-12-31"),
#                             out_dir = "E:/GIS/landcover/LCCMODIS500",
#                             internal = TRUE)
#   names(lcdata)[ncol(lcdata)]<-"LCTS500m"
#   return(lcdata)
# }

b1 <- brick("E:/GIS/landcover/NASA_ABoVE/1_to_17_mosaic.tif")
ABOVE <- function(surveys, yr) {
  data <- cbind(surveys,raster::extract(b1[[yr-1983]],surveys))
  names(data)[ncol(data)] <- "ABOVELC"
  return(data)
}
###############################

yearlist <- sort(unique(PKEYcombo$YEAR))

yearlist1 <- yearlist[yearlist < 2015]
yearlist1 <- yearlist1[yearlist1 > 1983]
yrdat <- setNames(data.frame(matrix(ncol = 22, nrow = 0)), c(names(PKEYcombo_TM10),"ABOVELC")) 
for (i in 1:length(yearlist1)) {
  yr1 <- PKEYcombo_TM10[PKEYcombo_TM10$YEAR==yearlist1[i],]
  yr1dat <- ABOVE(yr1, as.integer(yearlist1[i]))
  yrdat <- rbind(yrdat,yr1dat@data)
}
write.csv(yrdat,file=paste0(w,"PKEYannualLCABOVE.csv"),row.names=FALSE)

clim <- stack(lapply(grep("bio",list.files(paste0("E:/Climate/NRCAN/Annual/",2005),pattern=".tif$",full.names=TRUE),value=TRUE),FUN=raster))
cmi <- raster("E:/Climate/NRCAN/Annual/2005/cmi60_sum.asc")
cmi <- resample(cmi,clim[[1]])
clim <- addLayer(clim,cmi)
yearlist2 <- yearlist[yearlist < 2019]
yrdatclim <- setNames(data.frame(matrix(ncol = 41, nrow = 0)), c(names(PKEYcombo_DD),names(clim))) 
for (i in 1:length(yearlist2)) {
  yr1 <- PKEYcombo_DD[PKEYcombo_DD$YEAR==yearlist2[i],]
  yr1dat <- annclim(yr1, as.integer(yearlist2[i]))
  yrdatclim <- rbind(yrdatclim,as.data.frame(yr1dat))
}
write.csv(yrdatclim,file=paste0(w,"PKEYannualclim.csv"),row.names=FALSE)

yearlist3 <- yearlist[yearlist < 2012]
yearlist3 <- yearlist3[yearlist > 1999]
yrdatlc <- setNames(data.frame(matrix(ncol = 22, nrow = 0)), c(names(PKEYcombo_LCC),"LCTS250m")) 
for (i in 1:length(yearlist3)) {
  yr1 <- PKEYcombo_LCC[PKEYcombo_LCC$YEAR==yearlist3[i],]
  yr1dat <- CanadaMODISLC(yr1, as.integer(yearlist3[i]))
  yrdatlc <- rbind(yrdatlc,yr1dat@data)
}
write.csv(yrdatlc,file=paste0(w,"PKEYannualLC250.csv"),row.names=FALSE)

yearlist4 <- yearlist[yearlist > 2000]
yrdatlc <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c(names(PKEYcombo_DD),"LC500m")) 
for (i in 1:length(yearlist4)) {
  yr1 <- PKEYcombo_DD[PKEYcombo_DD$YEAR==yearlist4[i],]
  yr1dat <- MODIS500(yr1, as.integer(yearlist4[i]))
  yrdatlc <- rbind(yrdatlc,yr1dat@data)
}
write.csv(yrdatlc,file=paste0(w,"PKEYannualLC500.csv"),row.names=FALSE)


### Combine Data ####
toproad <- read.csv(paste0(w,"SS_TopoRoad.csv"))
lfclimate <- read.csv(paste0(w,"SS_LandformClimate.csv"))
vegbiomass <- read.csv(paste0(w,"PKEYcombo_VegBiomass.csv"))
annualclim <- read.csv(paste0(w,"PKEYannualclim.csv"))
annualABOVE <- read.csv(paste0(w,"PKEYannualLCABOVE.csv"))
annualMODIS250 <- read.csv(paste0(w,"PKEYannualLC250.csv"))
annualMODIS500 <- read.csv(paste0(w,"PKEYannualLC500.csv"))

datcombo <- left_join(toproad[,c(1,20,21,2:5,19)],lfclimate[,c(1:28)],by="SS")
datcombo1 <- left_join(annualclim[,c(1:23)],vegbiomass[,c(2,4:189)],by="PKEY")
datcombo1 <- left_join(datcombo1,annualABOVE[,c(2,4)],by="PKEY")
datcombo1 <- left_join(datcombo1,annualMODIS250[,c(2,4)],by="PKEY")
datcombo1 <- left_join(datcombo1,annualMODIS500[,c(2,4)],by="PKEY")
datcombo2 <- right_join(datcombo,datcombo1,by="SS")
#datcombo2001 <- datcombo2[datcombo2$YEAR>2000,]

IGBP <- read.csv("R:/Boreal/RshProjs/PopnStatus/Density_Popn/NationalModelsV3/IGBP_XWalk.csv")
NALCMS <- read.csv("R:/Boreal/RshProjs/PopnStatus/Density_Popn/NationalModelsV3/NALCMS_XWalk.csv")
ABOVE <- read.csv("R:/Boreal/RshProjs/PopnStatus/Density_Popn/NationalModelsV3/ABOVE_XWalk.csv")
datcombo2 <- left_join(datcombo2,ABOVE,by = c("ABOVELC"="ABOVEClass"))
datcombo2 <- left_join(datcombo2,NALCMS,by = c("LCTS250m"="NALCMSClass"))
datcombo2 <- left_join(datcombo2,IGBP,by = c("LC500m"="IGBPClass"))
datcombo2$LCX <- ifelse(datcombo2$ABOVEX>0, datcombo2$ABOVEX, ifelse(datcombo2$NALCMSX>0, datcombo2$NALCMSX, datcombo2$IGBPX))
datcombo2$LCX <- ifelse(is.na(datcombo2$LCX), datcombo2$ABOVEX, datcombo2$LCX)
datcombo2$LCX <- ifelse(is.na(datcombo2$LCX), datcombo2$NALCMSX, datcombo2$LCX)
datcombo2$LCX <- ifelse(is.na(datcombo2$LCX), datcombo2$IGBPX, datcombo2$LCX)
write.csv(datcombo2, file=paste0(w,"PKEYannualcombo.csv"),row.names=FALSE)

table(as.factor(datcombo2$YEAR), as.factor(datcombo2$LCX))
table(as.factor(datcombo2$YEAR), as.factor(datcombo2$LC500m))
table(as.factor(datcombo2$YEAR), as.factor(datcombo2$IGBPX))
table(as.factor(datcombo2$YEAR), as.factor(datcombo2$NALCMSX))
