library(raster)
library(dismo)
library(rpart)
library(maptools)
library(dplyr)
library(data.table)

w <-"G:/Boreal/NationalModelsV2/Quebec/"
bith <- read.csv("F:/BAM/BAMData/BITH/Aubry_FichierGlobalGRBI_Oct2015.csv")

PC <- bith[,c(2,4:8)]
names(PC) <- c("ABUND","SS","Date","YEAR","Month","Day")
PC$PKEY <- paste(PC$SS,PC$Year,PC$Month,PC$Day,sep=":")
PC$SPECIES <- "BITH"
PC2001 <- PC[PC$YEAR < 2006,] #n=1093
PC2011 <- PC[PC$YEAR > 2005,] #n=1361

write.csv(PC,paste(w,"BITHPC.csv",sep=""),row.names=FALSE)
write.csv(PC2001,paste(w,"BITHPC2001.csv",sep=""),row.names=FALSE)
write.csv(PC2011,paste(w,"BITHPC2011.csv",sep=""),row.names=FALSE)

PKEY <- unique(PC[,c(2,4,7)])
write.csv(PKEY,paste(w,"BITHPKEY.csv",sep=""),row.names=FALSE)


SS <- bith[,c(4,10:11)]
names(SS) <- c("SS","Y","X")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
coordinates(SS) <- c("X", "Y") 
proj4string(SS) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSQC <- as.data.frame(spTransform(SS, LCC))

coordinates(SSQC) <- c("X", "Y") 
proj4string(SSQC) <- LCC
ss <- as.data.frame(spTransform(SSQC, lazea))
ss$SS <- as.character(ss$SS)

#Landform (100-m) Lazea projection
TPI <- raster("E:/GIS/topoedaphic/TPI.tif")
TRI <- raster("E:/GIS/topoedaphic/TRI.tif")
slope <- raster("E:/GIS/topoedaphic/slope.tif")
roughness <- raster("E:/GIS/topoedaphic/roughness.tif")

ss <- cbind(ss,"TPI"=extract(TPI,as.matrix(cbind(ss$X,ss$Y))))
ss <- cbind(ss,"TRI"=extract(TRI,as.matrix(cbind(ss$X,ss$Y))))
ss <- cbind(ss,"slope"=extract(slope,as.matrix(cbind(ss$X,ss$Y))))
ss <- cbind(ss,"roughness"=extract(roughness,as.matrix(cbind(ss$X,ss$Y))))

#Beaudoin biomass and landscape proportion, LCC projection
b2011 <- list.files("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/",pattern="tif$")
setwd("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) {bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
bs2011quebec <- crop(bs2011,quebec)
bs2011quebec <- mask(bs2011quebec,quebec)

bs2011_Gauss750 <- brick("G:/Boreal/NationalModelsV2/bs2011_750.grd")
names(bs2011_Gauss750) <- gsub("SpeciesGroups","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Species","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Structure","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Landcover","Landsc750",names(bs2011_Gauss750))
bs2011quebec_Gauss750 <- crop(bs2011,quebec)
bs2011quebec_Gauss750 <- mask(bs2011quebec_Gauss750,quebec)

b2001 <- list.files("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2001/",pattern="tif$")
setwd("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2001/")
bs2001 <- stack(raster(b2001[1]))
for (i in 2:length(b2001)) {bs2001 <- addLayer(bs2001, raster(b2001[i]))}
names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
bs2001quebec <- crop(bs2001,quebec)
bs2001quebec <- mask(bs2001quebec,quebec)

bs2001_Gauss750 <- brick("G:/Boreal/NationalModelsV2/bs2001_750.grd")
names(bs2001_Gauss750) <- gsub("SpeciesGroups","Landsc750",names(bs2001_Gauss750))
names(bs2001_Gauss750) <- gsub("Species","Landsc750",names(bs2001_Gauss750))
names(bs2001_Gauss750) <- gsub("Structure","Landsc750",names(bs2001_Gauss750))
names(bs2001_Gauss750) <- gsub("Landcover","Landsc750",names(bs2001_Gauss750))
bs2001quebec_Gauss750 <- crop(bs2001,quebec)
bs2001quebec_Gauss750 <- mask(bs2001quebec_Gauss750,quebec)

#landcover and derived variables, LCC projection
nalc <- raster("F:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/NA_LandCover_2005_LCC.img")
lc <- brick("G:/Boreal/NationalModelsV2/quebec/lc.grd")
lf <- raster("E:/GIS/topoedaphic/lf1k.tif")

#Road (on/off)
road <- raster("G:/Boreal/NationalModelsV2/roadonoff1.tif")
mr <- c(1, 2500000, 1,  NA, NA, 0)
rcroad <- matrix(mr, ncol=3, byrow=TRUE)
rrc <- reclassify(road,rcroad)

dat2011 <- cbind(SSQC, extract(bs2011quebec,as.matrix(cbind(SSQC$X,SSQC$Y))))
dat2011 <- cbind(dat2011, extract(bs2011quebec_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y))))
dat2011 <-cbind(dat2011, extract(nalc,as.matrix(cbind(dat2011$X,dat2011$Y)))) 
names(dat2011)[ncol(dat2011)] <- "nalc"
dat2011 <- cbind(dat2011, extract(rrc,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "rrc"
dat2011 <- cbind(dat2011, extract(lf,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "lf"
dat2011 <- cbind(dat2011, extract(lc,as.matrix(cbind(dat2011$X,dat2011$Y))))
dat_2011 <- merge(as.data.frame(dat2011)[,1:ncol(dat2011)], ss[,c(1,4:7)], by=c("SS"))
dat_2011 <- na.omit(dat_2011)
dat_2011 <- inner_join(dat_2011, PKEY[,1:2], by=c("SS")) 
dat_2011 <- distinct(dat_2011[dat_2011$YEAR > 2005,1:(ncol(dat_2011)-1)]) #n=1340
write.csv(dat_2011,"G:/Boreal/NationalModelsV2/quebec/QCBITHdat2011.csv",row.names=FALSE)

dat2001 <- cbind(SSQC, extract(bs2001quebec,as.matrix(cbind(SSQC$X,SSQC$Y))))
dat2001 <- cbind(dat2001, extract(bs2001quebec_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y))))
dat2001 <-cbind(dat2001, extract(nalc,as.matrix(cbind(dat2001$X,dat2001$Y)))) 
names(dat2001)[ncol(dat2001)] <- "nalc"
dat2001 <- cbind(dat2001, extract(rrc,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "rrc"
dat2001 <- cbind(dat2001, extract(lf,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "lf"
dat2001 <- cbind(dat2001, extract(lc,as.matrix(cbind(dat2001$X,dat2001$Y))))
dat_2001 <- merge(as.data.frame(dat2001)[,1:ncol(dat2001)], ss[,c(1,4:7)], by=c("SS")) 
dat_2001 <- na.omit(dat_2001)
dat_2001 <- inner_join(dat_2001, PKEY[,1:2], by=c("SS")) 
dat_2001 <- distinct(dat_2001[dat_2001$YEAR < 2006,1:(ncol(dat_2001)-1)]) #n=934
write.csv(dat_2001,"G:/Boreal/NationalModelsV2/quebec/QCBITHdat2001.csv",row.names=FALSE)


