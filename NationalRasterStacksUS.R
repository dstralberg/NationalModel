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
# load(paste(w,"BAMdb-GNMsubset-2019-03-01.RData",sep=""))

# bcrs <- c("bcr4_100km.shp","bcr5_100km.shp","bcr9_100km.shp","bcr10_100km.shp","bcr11_100km.shp","bcr12_100km.shp","bcr13_100km.shp","bcr14_100km.shp")
# bcrs2 <- c("bcr60_100km.shp","bcr61_100km.shp","bcr70_100km.shp","bcr71_100km.shp","bcr80_100km.shp","bcr81_100km.shp", "bcr82_100km.shp","bcr83_100km.shp")
# canada <- shapefile("F:/GIS/basemaps/canadaLCC.shp")
bcrsus <- c("bcrus2_100km.shp", "bcrus4_100km.shp", "bcrus5_100km.shp", "bcrus14_100km.shp", "bcrus13_100km.shp", "bcrus12_100km.shp", "bcrus11_100km.shp", "bcrus10_100km.shp")

#Adaptwest baseline climate variables
cur <- "E:/CMIP5/baseline19812010/"
setwd(cur)
clim <- list.files(cur, pattern =".asc$")
curclim<-stack(clim)

# canroad <- shapefile(paste(w,"canroadlcc.shp",sep=""))
# roadlen <- raster(paste(w,"roadlen.asc",sep=""))
# road <- raster(paste(w,"roadonoff1.tif",sep=""))
# mr <- c(1, 2500000, 1,  NA, NA, 0)
# rcroad <- matrix(mr, ncol=3, byrow=TRUE)
# rrc <- reclassify(road,rcroad)
# writeRaster(rrc,file=paste(w,"road10",sep=""),overwrite=TRUE)
#rrc <- raster(paste(w,"road10.grd",sep=""))

road1k <- raster("E:/GIS/disturbance/VenterEtAlFootprint/RoadsLCC.tif")
roadr <- resample(road1k, topog, method="ngb")

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
for (i in 1:length(bcrsus)){
  bcr <- shapefile(bcrsus[i])
  clim1 <- mask(crop(curclim,bcr),bcr)
  landcov1 <- mask(crop(landcover,bcr),bcr)
  topo1 <- mask(crop(topog,bcr),bcr)
  bcrr <- rasterize(bcr,clim1[[1]])
  bs <- stack(bcrr)
  names(bs) <- "bcr"
  for (j in 1:nlayers(clim1)) {
    bs <- addLayer(bs, clim1[[j]])}
  for (j in 1:nlayers(landcov1)) {
    bs <- addLayer(bs, landcov1[[j]])}
  for (j in 1:nlayers(topo1)) {
    bs <- addLayer(bs, topo1[[j]])}
  rrc1 <- crop(extend(road1k,bcr),bcr)
  ROAD <- resample(rrc1,clim1,method="ngb")
  ROAD <- mask(ROAD,bcr)
  bs <- stack(bs,ROAD)
  names(bs)[nlayers(bs)] <- "ROAD"
  writeRaster(bs,file=paste(w,gsub("_100km.shp","",bcrsus[i]),"all_1km",sep=""),overwrite=TRUE)
}


