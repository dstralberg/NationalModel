library(raster)
library(sf)
library(maptools)
library(dplyr)
library(data.table)
library(reshape2)
rasterOptions(chunksize = 1e+04, maxmemory = 1e+06)
#source("CTI.R")

LCC <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
lazea <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
library(raster)
w <-"G:/Boreal/NationalModelsV2/"
load(paste(w,"BAMdb-GNMsubset-2019-03-01.RData",sep=""))

bcrs <- c("bcr4_100km.shp","bcr5_100km.shp","bcr6_100km.shp","bcr7_100km.shp","bcr8_100km.shp","bcr9_100km.shp","bcr10_100km.shp","bcr11_100km.shp","bcr12_100km.shp","bcr13_100km.shp","bcr14_100km.shp")

#Adaptwest baseline climate variables
cur <- "E:/CMIP5/baseline19812010/"
setwd(cur)
clim <- list.files(cur, pattern =".asc$")
curclim<-stack(clim)

# canroad <- shapefile(paste(w,"canroadlcc.shp",sep=""))
# roadlen <- raster(paste(w,"roadlen.asc",sep=""))
road <- raster(paste(w,"roadonoff1.tif",sep=""))
mr <- c(1, 2500000, 1,  NA, NA, 0)
rcroad <- matrix(mr, ncol=3, byrow=TRUE)
rrc <- reclassify(road,rcroad)

#MODIS-based landcover (250-m)
# nalc2005 <- raster("F:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/NA_LandCover_2005_LCC.img")
# nalc <- resample(nalc2005,curclim,method="ngb")
# writeRaster(nalc,paste(w,"nalc-1km.tif",sep=""))
# mwat <- c(0, 17.1, 0,  17.9, 18.1, 1,  18.9, 20, 0)
# rclwat <- matrix(mwat, ncol=3, byrow=TRUE)
# water <- reclassify(nalc2005,rclwat)
# murb <- c(0, 16.1, 0,  16.9, 17.1, 1,  17.9, 20, 0)
# rclurb <- matrix(murb, ncol=3, byrow=TRUE)
# urban <- reclassify(nalc2005,rclurb)
# mag <- c(0, 14.1, 0,  14.9, 15.1, 1,  15.9, 20, 0)
# rclag <- matrix(mag, ncol=3, byrow=TRUE)
# ag <- reclassify(nalc2005,rclag)
# x <- stack(urban,ag)
# urbag <- max(x)
# fw750<-focalWeight(x=urbag,d=750,type="Gauss") #Gaussian filter with sigma=750 (tapers off around 2km)
# dev750 <- focal(urbag,w=fw750,na.rm=TRUE)
# led750 <- focal(water,w=fw750,na.rm=TRUE)
# landcov <- stack(dev750,led750)
# landcover <- resample(landcov,curclim)
# landcover <- addLayer(landcover,nalc)
# names(landcover) <- c("dev750","led750","nalc")
# writeRaster(landcover,file=paste(w,"landcov_1km",sep=""),overwrite=TRUE)
landcover <- brick(paste(w,"landcov_1km.grd",sep=""))

#kNN biomass layers, 2011 (250-m)
# b2011 <- list.files("F:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/",pattern="tif$")
# setwd("F:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/")
# bs2011 <- stack(raster(b2011[1]))
# for (i in 2:length(b2011)) {bs2011 <- addLayer(bs2011, raster(b2011[i]))}
# names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
# 
# climcan <- crop(curclim,bs2011[[1]])
# biomass <- resample(bs2011,climcan)
# writeRaster(biomass, file=paste(w,"bs2011_1km",sep=""))
biomass <- brick(paste(w,"bs2011_1km.grd",sep=""))

#Gaussian filter for 2011 biomass layers
# bs2011_Gauss750 <- brick(paste(w,"bs2011_750.grd",sep=""))
# names(bs2011_Gauss750) <- gsub("SpeciesGroups","Landsc750",names(bs2011_Gauss750))
# names(bs2011_Gauss750) <- gsub("Species","Landsc750",names(bs2011_Gauss750))
# names(bs2011_Gauss750) <- gsub("Structure","Landsc750",names(bs2011_Gauss750))
# names(bs2011_Gauss750) <- gsub("Landcover","Landsc750",names(bs2011_Gauss750))
# landscape <- resample(bs2011_Gauss750,biomass)
# writeRaster(landscape, file=paste(w,"bs2011_750_1km",sep=""))
landscape <- brick(paste(w,"bs2011_750_1km.grd",sep=""))
names(landscape)[1:4] <- gsub("v1","v1.1",names(landscape)[1:4])

# Landform (100-m)
# lf <- raster("F:/GIS/topoedaphic/lf1k.tif")
# TPI <- raster("F:/GIS/topoedaphic/tpi1k.tif")
# TRI <- raster("F:/GIS/topoedaphic/tri1k.tif")
# slope <- raster("F:/GIS/topoedaphic/slope1k.tif")
# roughness <- raster("F:/GIS/topoedaphic/roughness1k.tif")
# topo <- stack(TPI,TRI,slope,roughness)
# topog <- crop(topo,lf)
# topog <- addLayer(topog,lf)
# names(topog) <- c("TPI","TRI","slope","roughness","lf")
# writeRaster(topog,file=paste(w,"topography_1km",sep=""),overwrite=TRUE)
topog <-brick(paste(w,"topography_1km.grd",sep=""))
names(topog) <- c("TPI","TRI","slope","roughness","lf")

setwd(w)
for (i in 1:length(bcrs)){
  vars <- CN[[i]]
  bcr <- shapefile(bcrs[i])
  biomass1 <- mask(crop(biomass,bcr),bcr)
  landscape1 <- mask(crop(landscape,bcr),bcr)
  clim1 <- mask(crop(curclim,biomass1[[1]]),biomass1[[1]])
  landcov1 <- mask(crop(landcover,biomass1[[1]]),biomass1[[1]])
  topo1 <- mask(crop(topog,biomass1[[1]]),biomass1[[1]])
  bcrr <- rasterize(bcr,curclim[[1]])
  bcrc <- mask(crop(bcrr,biomass1[[1]]),biomass1[[1]])
  bs <- stack(bcrc)
  names(bs) <- "bcr"
  for (j in 1:nlayers(clim1)) {
    if (names(clim1)[j] %in% vars) {bs <- addLayer(bs, clim1[[j]])}
  }
  for (j in 1:nlayers(biomass1)) {
    if (names(biomass1)[j] %in% vars) {bs <- addLayer(bs, biomass1[[j]])}
  }
  for (j in 1:nlayers(landscape1)) {
    if (names(landscape1)[j] %in% vars) {bs <- addLayer(bs, landscape1[[j]])}
  }
  for (j in 1:nlayers(landcov1)) {
    if (names(landcov1)[j] %in% vars) {bs <- addLayer(bs, landcov1[[j]])}
  }  
  for (j in 1:nlayers(topo1)) {
    if (names(topo1)[j] %in% vars) {bs <- addLayer(bs, topo1[[j]])}
  }   
  writeRaster(bs,file=paste(w,gsub("_100km.shp","",bcrs[i]),"_1km",sep=""),overwrite=TRUE)
}

setwd(w)
#write categorical rasters
for (i in 1:length(bcrs)){
  bcr <- shapefile(bcrs[i])
  biomass1 <- mask(crop(biomass,bcr),bcr)
  nalc <- mask(crop(landcover[[1]],biomass1[[1]]),biomass1[[1]])
  lf <- mask(crop(topog[[5]],biomass1[[1]]),biomass1[[1]])
  rrc1 <- crop(extend(rrc,biomass1[[1]]),biomass1[[1]])
  ROAD <- mask(rrc1,biomass1[[1]])
  bcrr <- rasterize(bcr,curclim[[1]])
  bcrc <- mask(crop(bcrr,lf),lf)
  bs <- stack(bcrc,nalc,lf,ROAD)
  names(bs) <- c("bcr","nalc","lf","ROAD")
  writeRaster(bs,file=paste(w,gsub("_100km.shp","",bcrs[i]),"cat_1km",sep=""),overwrite=TRUE)
}

setwd(w)
#write climate rasters
for (i in 1:length(bcrs)){
  bcr <- shapefile(bcrs[i])
  biomass1 <- mask(crop(biomass,bcr),bcr)
  clim1 <- mask(crop(curclim,biomass1[[1]]),biomass1[[1]])
  bcrr <- rasterize(bcr,curclim[[1]])
  bcrc <- mask(crop(bcrr,biomass1[[1]]),biomass1[[1]])
  bs <- stack(bcrc)
  names(bs) <- "bcr"
  for (j in 1:nlayers(clim1)) {
    bs <- addLayer(bs, clim1[[j]])
  }
  writeRaster(bs,file=paste(w,gsub("_100km.shp","",bcrs[i]),"clim_1km",sep=""),overwrite=TRUE)
}
