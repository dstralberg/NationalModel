library(raster)
library(sf)
library(maptools)
library(dplyr)
library(data.table)
library(reshape2)
#source("CTI.R")

LCC <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
lazea <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
w <-"G:/Boreal/NationalModelsV2/"

load("F:/BAM/BAMData/BAMdb-patched-xy.RData")	
sslaz <- st_transform(ss,lazea)

#MODIS-based landcover (250-m)
nalc2005 <- raster("F:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/NA_LandCover_2005_LCC.img")
ss <- cbind(ss,"nalc"=extract(nalc2005,st_coordinates(ss)))
urbag <- raster("G:/Boreal/NationalModelsV2/urbag2011_lcc1.tif")
fw750<-focalWeight(x=urbag,d=750,type="Gauss") #Gaussian filter with sigma=750 (tapers off around 2km)
dev750 <- focal(urbag,w=fw750,na.rm=TRUE)
ss <- cbind(ss,"dev750"=extract(dev750,st_coordinates(ss)))
wat <- raster("G:/Boreal/NationalModelsV2/wat2011_lcc1.tif")
led750 <- focal(wat,w=fw750,na.rm=TRUE)
ss <- cbind(ss,"led750"=extract(led750,st_coordinates(ss)))

#kNN biomass layers, 2011 (250-m)
b2011 <- list.files("F:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/",pattern="tif$")
setwd("F:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) {bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
ss <- cbind(ss,extract(bs2011,st_coordinates(ss)))

#Gaussian filter for 2011 biomass layers
# bs2011_Gauss750<-brick(focal(bs2011[[1]],w=fw750,na.rm=TRUE))
# names(bs2011_Gauss750)<-names(bs2011)[[1]]
# for(i in 2:nlayers(bs2011)){
#   bs2011_Gauss750<-addLayer(bs2011_Gauss750,focal(bs2011[[i]],w=fw750,na.rm=TRUE))
#   names(bs2011_Gauss750)[i]<-names(bs2011)[[i]]
# }
# writeRaster(bs2011_Gauss750,file=paste(w,"bs2011_750",sep=""))
bs2011_Gauss750 <- brick(paste(w,"bs2011_750.grd",sep=""))
names(bs2011_Gauss750) <- gsub("SpeciesGroups","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Species","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Structure","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Landcover","Landsc750",names(bs2011_Gauss750))
ss <- cbind(ss,extract(bs2011_Gauss750,st_coordinates(ss)))

#Landform (100-m)
lf <- raster("D:/NorthAmerica/topo/lf_lcc1.tif")
ss <- cbind(ss,"lf"=extract(lf,st_coordinates(ss)))
#dem <- raster("F:/GIS/topoedaphic/nadem100laz.tif")
#terrain <- terrain(dem, opt=c('slope','TPI','TRI','roughness'))
#setwd("F:/GIS/topoedaphic/")
#writeRaster(terrain, filename=c("slope","TPI","TRI","roughness"), bylayer=TRUE, format="GTiff", overwrite=TRUE)
TPI <- raster("F:/GIS/topoedaphic/TPI.tif")
ss <- cbind(ss,"TPI"=extract(TPI,st_coordinates(sslaz)))
TRI <- raster("F:/GIS/topoedaphic/TRI.tif")
ss <- cbind(ss,"TRI"=extract(TRI,st_coordinates(sslaz)))
slope <- raster("F:/GIS/topoedaphic/slope.tif")
ss <- cbind(ss,"slope"=extract(slope,st_coordinates(sslaz)))
roughness <- raster("F:/GIS/topoedaphic/roughness.tif")
ss <- cbind(ss,"roughness"=extract(roughness,st_coordinates(sslaz)))

#ss1 <- CTI(ss[,c(1:4,199)])

#Adaptwest baseline climate variables
cur <- "E:/CMIP5/baseline19812010/"
setwd(cur)
clim <- list.files(cur, pattern =".asc$")
curclim<-stack(clim)
ss <- cbind(ss,extract(curclim,st_coordinates(ss)))

#Roads buffered by 100-m
road <- raster(paste(w,"roadonoff1.tif",sep=""))
mr <- c(1, 2500000, 1,  NA, NA, 0)
rcroad <- matrix(mr, ncol=3, byrow=TRUE)
rrc <- reclassify(road,rcroad)
ss <- cbind(ss,"ROAD"=extract(rrc,st_coordinates(ss)))

write.csv(ss,file=paste(w,"ss_2011attributes.csv",sep=""),row.names==FALSE)
save(ss,file=paste(w,"ss_2011attributes.RData",sep=""))

temp <- load(paste0(w,"ss_2011attributes.RData"))
