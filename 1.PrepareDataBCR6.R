library(raster)
library(maptools)
library(dplyr)
library(data.table)
library(reshape2)

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
w <-"G:/Boreal/NationalModelsV2/BCR6/"
bcr6 <- raster("G:/Boreal/NationalModelsV2/BCR6/bcr6.tif")

load("F:/BAM/BAMData/data_package_2016-04-18.Rdata")	
load("F:/BAM/BAMData/offsets-v3_2016-04-18.Rdata")
coordinates(SS) <- c("X", "Y") 
proj4string(SS) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSBAM <- as.data.frame(spTransform(SS, LCC))
PCBAM <- PCTBL
PKEYBAM <- PKEY

load("F:/BAM/BAMData/atlas_data_processed-20181018.RData")
SS <- na.omit(SS)
coordinates(SS) <- c("X", "Y") 
proj4string(SS) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSAtlas <- as.data.frame(spTransform(SS, LCC))
PCAtlas <- PCTBL
PKEYAtlas <- PKEY
names(PKEYAtlas)[4] <- "YEAR"

load("F:/BAM/BamData/ARU/nwt-wildtrax-offsets-2019-01-16.RData")
SSWT <- unique(dd[dd$Y>0,c(33,39:40)])
SSWT1 <- na.omit(SSWT)
coordinates(SSWT1) <- c("X", "Y") 
proj4string(SSWT1) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSWTLC <- as.data.frame(spTransform(SSWT1, LCC))
PKEYWT <- unique(dd[,c(33,34,36)])
#PCWT <- dd[,c(33,34,36,38,47)]
PCWT <- melt(y)
names(PCWT) <- c("PKEY","SPECIES","ABUND")
offWT <- data.table(melt(off))
names(offWT) <- c("PKEY","SPECIES","logoffset")
offWT$SPECIES <- as.character(offWT$SPECIES)
offWT$PKEY <- as.character(offWT$PKEY)
write.csv(offWT,file="G:/Boreal/NationalModelsV2/BCR6/offwt.csv", row.names=FALSE)

load("F:/BAM/BamData/ARU/nwt-BU-offsets-2019-01-14.RData")
SSBU <- unique(dd[,c(14,20:21)])
coordinates(SSBU) <- c("X", "Y") 
proj4string(SSBU) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSBULC <- as.data.frame(spTransform(SSBU, LCC))
PKEYBU <- unique(dd[,c(14,15,17)])
PCBU <- melt(y)
names(PCBU) <- c("PKEY","SPECIES","ABUND")
offBU <- data.table(melt(off))
names(offBU) <- c("PKEY","SPECIES","logoffset")
offBU$SPECIES <- as.character(offBU$SPECIES)
offBU$PKEY <- as.character(offBU$PKEY)
#write.csv(offBU,file="G:/Boreal/NationalModelsV2/BCR6/offbu.csv", row.names=FALSE)

load("F:/BAM/BamData/ARU/nwt-offsets-2019-02-01.RData")
SSBU2 <- unique(dd[,c(14,20:21)])
coordinates(SSBU2) <- c("X", "Y") 
proj4string(SSBU2) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSBU2LC <- as.data.frame(spTransform(SSBU2, LCC))
SSBULC <- rbind(SSBULC,SSBU2LC)
PKEYBU2 <- unique(dd[,c(14,15,17)])
PKEYBU <- rbind(PKEYBU,PKEYBU2)
PCBU2 <- melt(y)
names(PCBU2) <- c("PKEY","SPECIES","ABUND")
PCBU <- rbind(PCBU,PCBU2)
offBU2 <- data.table(melt(off))
names(offBU2) <- c("PKEY","SPECIES","logoffset")
offBU2$SPECIES <- as.character(offBU2$SPECIES)
offBU2$PKEY <- as.character(offBU2$PKEY)
offBU <- rbind(offBU,offBU2)
write.csv(offBU,file="G:/Boreal/NationalModelsV2/BCR6/offbu.csv", row.names=FALSE)

SScombo <- rbind(SSBAM[,c(2,48,49)],SSAtlas[,c(1,6,7)],SSWTLC,SSBULC)
SSBCR6 <- cbind(SScombo, extract(bcr6, as.matrix(cbind(SScombo$X,SScombo$Y)))) #n=312906
SSBCR6 <- na.omit(SSBCR6) #n=52362
SSBCR6 <- SSBCR6[,1:3]
PKEYcombo <- rbind(PKEYBAM[,c(1,2,8)],PKEYAtlas[,c(1,2,4)],PKEYWT,PKEYBU)
PCcombo <- rbind(PCBAM[,c(3:5)], PCAtlas[,c(2,4:5)], PCWT, PCBU) 

eco <- raster("F:/GIS/ecoregions/CEC/quebececo1.tif")
nalc <- raster("F:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/NA_LandCover_2005_LCC.img")
bcr6 <- raster("G:/Boreal/NationalModelsV2/BCR6/bcr6.tif")

mwat <- c(0, 17.1, 0,  17.9, 18.1, 1,  18.9, 20, 0)
rclwat <- matrix(mwat, ncol=3, byrow=TRUE)
water <- reclassify(nalc,rclwat)
murb <- c(0, 16.1, 0,  16.9, 17.1, 1,  17.9, 20, 0)
rclurb <- matrix(murb, ncol=3, byrow=TRUE)
urban <- reclassify(nalc,rclurb)
mag <- c(0, 14.1, 0,  14.9, 15.1, 1,  15.9, 20, 0)
rclag <- matrix(mag, ncol=3, byrow=TRUE)
ag <- reclassify(nalc,rclag)
x <- stack(urban,ag)
urbag <- max(x)

fw750<-focalWeight(x=urbag,d=750,type="Gauss") #Gaussian filter with sigma=750 (tapers off around 2km)
dev750 <- focal(urbag,w=fw750,na.rm=TRUE)
led750 <- focal(water,w=fw750,na.rm=TRUE)

# ua6 <- crop(urbag,bcr6)
# ua6 <- mask(ua6,bcr6)
# dev25 <- focal(ua6, fun=mean, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)
# wat6 <- crop(wat,bcr6)
# wat6 <- mask(wat6,bcr6)
# led25 <- focal(wat6, fun=mean, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)
vrug <- raster("G:/Boreal/NationalModelsV2/BCR6/vrug_bcr6.tif")
wet <- raster("G:/Boreal/NationalModelsV2/BCR6/HWL_BCR6.tif")
wet250 <- raster("G:/Boreal/NationalModelsV2/BCR6/wet250_BCR6.tif")

lf <- raster("D:/NorthAmerica/topo/lf_lcc1.tif")
lf6 <- crop(lf,bcr6)
lf250 <- resample(lf6, bs2011bcr6, method='ngb')

road <- raster("G:/Boreal/NationalModelsV2/BCR6/lineDensityMap_BCR6_Roads_v0.2.3_FFT_radius10000_t0.tif")
bs<-stack("G:/Boreal/NationalModelsV2/bcr6clim_1km.grd")

# b2011 <- list.files("F:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/",pattern="tif$")
# setwd("F:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/")
# bs2011 <- stack(raster(b2011[1]))
# for (i in 2:length(b2011)) {bs2011 <- addLayer(bs2011, raster(b2011[i]))}
# names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
# bs2011bcr6 <- crop(bs2011,bcr6)
# bs2011bcr6 <- mask(bs2011bcr6,bcr6)
# bs2011bcr6 <- dropLayer(bs2011bcr6, c(1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28,29,30,31,32,33,34,35,37,38,39,40,41,42,43,46,47,48,49,52,53,54,55,56,57,59,60,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,90,91,92,93))
# 
# bs2011bcr6 <- addLayer(bs2011bcr6,wat6)
# names(bs2011bcr6)[nlayers(bs2011bcr6)] <- "wat"
# bs2011bcr6 <- addLayer(bs2011bcr6,led25)
# names(bs2011bcr6)[nlayers(bs2011bcr6)] <- "led25"
# bs2011bcr6 <- addLayer(bs2011bcr6,ua6)
# names(bs2011bcr6)[nlayers(bs2011bcr6)] <- "urbag"
# bs2011bcr6 <- addLayer(bs2011bcr6,dev25)
# names(bs2011bcr6)[nlayers(bs2011bcr6)] <- "dev25"
# bs2011bcr6 <- addLayer(bs2011bcr6,vrug)
# names(bs2011bcr6)[nlayers(bs2011bcr6)] <- "vrug"
# bs2011bcr6 <- addLayer(bs2011bcr6, lf250)
# names(bs2011bcr6)[nlayers(bs2011bcr6)] <- "landform"
# bs2011bcr6 <- addLayer(bs2011bcr6, wet250)
# names(bs2011bcr6)[nlayers(bs2011bcr6)] <- "wet"
# writeRaster(bs2011bcr6,file=paste(w,"bcr6_2011rasters250",sep=""),overwrite=TRUE)
bs2011bcr6 <- stack(paste(w,"bcr6_2011rasters250",sep=""))

# bs2011bcr6_1km <- resample(bs2011bcr6, bs, method="ngb")
# bs2 <- stack(bs,bs2011bcr6_1km[[14:20]])
# writeRaster(bs2,file=paste(w,"bcr6_1km",sep=""),overwrite=TRUE)
bs2 <- stack(paste(w,"bcr6_1km",sep=""))

# b2001 <- list.files("F:/GIS/landcover/Beaudoin/Processed_sppBiomass/2001/",pattern="tif$")
# setwd("F:/GIS/landcover/Beaudoin/Processed_sppBiomass/2001/")
# bs2001 <- stack(raster(b2001[1]))
# for (i in 2:length(b2001)) {bs2001 <- addLayer(bs2001, raster(b2001[i]))}
# names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
# bs2001bcr6 <- crop(bs2001,bcr6)
# bs2001bcr6 <- mask(bs2001bcr6,bcr6)
# bs2001bcr6 <- dropLayer(bs2001bcr6, c(1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28,29,30,31,32,33,34,35,37,38,39,40,41,42,43,46,47,48,49,52,53,54,55,56,57,59,60,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,90,91,92,93))
# 
# bs2001bcr6 <- addLayer(bs2001bcr6,wat6)
# names(bs2001bcr6)[nlayers(bs2001bcr6)] <- "wat"
# bs2001bcr6 <- addLayer(bs2001bcr6,led25)
# names(bs2001bcr6)[nlayers(bs2001bcr6)] <- "led25"
# bs2001bcr6 <- addLayer(bs2001bcr6,ua6)
# names(bs2001bcr6)[nlayers(bs2001bcr6)] <- "urbag"
# bs2001bcr6 <- addLayer(bs2001bcr6,dev25)
# names(bs2001bcr6)[nlayers(bs2001bcr6)] <- "dev25"
# bs2001bcr6 <- addLayer(bs2001bcr6,vrug)
# names(bs2001bcr6)[nlayers(bs2001bcr6)] <- "vrug"
# bs2001bcr6 <- addLayer(bs2001bcr6, lf250)
# names(bs2001bcr6)[nlayers(bs2001bcr6)] <- "landform"
# bs2001bcr6 <- addLayer(bs2001bcr6, wet250)
# names(bs2001bcr6)[nlayers(bs2001bcr6)] <- "wet"
# writeRaster(bs2001bcr6,file=paste(w,"bcr6_2001rasters250",sep=""),overwrite=TRUE)
bs2001bcr6 <- stack(paste(w,"bcr6_2001rasters250",sep=""))

#bs2011bcr6 <- dropLayer(bs2011bcr6,19)
dat2011 <- cbind(SSBCR6, extract(bs2011bcr6,as.matrix(cbind(SSBCR6$X,SSBCR6$Y))))
dat2011 <-cbind(dat2011,extract(nalc,as.matrix(cbind(dat2011$X,dat2011$Y)))) 
names(dat2011)[ncol(dat2011)] <- "LCC"
dat2011 <-cbind(dat2011,extract(lf6,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "landform"
dat2011$SS <- as.character(dat2011$SS) #n=52362
dat2011 <- na.omit(dat2011) #n=51966
dat_2011 <- inner_join(dat2011, PKEYcombo[,2:3], by=c("SS")) #n=166535
dat_2011 <- distinct(dat_2011[dat_2011$YEAR > 2005,1:23]) #n=34022
write.csv(dat_2011,"G:/Boreal/NationalModelsV2/BCR6/bcr6_dat2011_v2.csv",row.names=FALSE)

bs2001bcr6 <- dropLayer(bs2001bcr6,19)
dat2001 <- cbind(SSBCR6, extract(bs2001bcr6,as.matrix(cbind(SSBCR6$X,SSBCR6$Y))))
dat2001 <-cbind(dat2001,extract(nalc,as.matrix(cbind(dat2001$X,dat2001$Y)))) 
names(dat2001)[ncol(dat2001)] <- "LCC"
dat2001 <-cbind(dat2001,extract(lf6,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "landform"
dat2001$SS <- as.character(dat2001$SS)
dat2001 <- na.omit(dat2001)
dat_2001 <- inner_join(dat2001, PKEYcombo[,2:3], by=c("SS")) 
dat_2001 <- distinct(dat_2001[dat_2001$YEAR < 2006,1:24]) #n=18837
write.csv(dat_2001,"G:/Boreal/NationalModelsV2/BCR6/bcr6_dat2001_v2.csv",row.names=FALSE)

#climate + terrain layers
cdat <- cbind(SSBCR6, extract(bs,as.matrix(cbind(SSBCR6$X,SSBCR6$Y))))
cdat1 <- cbind(SSBCR6, extract(bs2011bcr6[[14:20]],as.matrix(cbind(SSBCR6$X,SSBCR6$Y))))
cdat2 <- cbind(cdat,cdat1[,4:ncol(cdat1)])
cdat3 <- inner_join(cdat2, PKEYcombo[,2:3], by=c("SS")) #n=167850
write.csv(cdat3,"G:/Boreal/NationalModelsV2/BCR6/bcr6_cdat_v2.csv",row.names=FALSE)

PC <- inner_join(PCcombo,PKEYcombo,by=c("PKEY")) 
PCBCR6 <- inner_join(PC, SSBCR6[,1:3], by=c("SS")) 
PCBCR6$SS <- as.character(PCBCR6$SS)
PCBCR6$PKEY <- as.character(PCBCR6$PKEY)
PCBCR6$SPECIES <- as.character(PCBCR6$SPECIES)
PCBCR62001 <- PCBCR6[PCBCR6$YEAR < 2006,] 
PCBCR62011 <- PCBCR6[PCBCR6$YEAR > 2005,] 
write.csv(PCBCR62011,"G:/Boreal/NationalModelsV2/BCR6/BCR6PC2011_v2.csv",row.names=FALSE)
write.csv(PCBCR62001,"G:/Boreal/NationalModelsV2/BCR6/BCR6PC2001_v2.csv",row.names=FALSE)
write.csv(PCBCR6,"G:/Boreal/NationalModelsV2/BCR6/BCR6PC_v2.csv",row.names=FALSE)
