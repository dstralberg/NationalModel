library(raster)
library(sf)
library(maptools)
library(dplyr)
library(data.table)
library(reshape2)
rasterOptions(chunksize = 1e+04, maxmemory = 1e+06)

LCC <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
lazea <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
library(raster)
w <-"F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/"

bcrs <- c("bcr4_100km.shp","bcr5_100km.shp","bcr9_100km.shp","bcr10_100km.shp","bcr11_100km.shp","bcr12_100km.shp","bcr13_100km.shp","bcr14_100km.shp")
bcrs2 <- c("bcr60_100km.shp","bcr61_100km.shp","bcr70_100km.shp","bcr71_100km.shp","bcr80_100km.shp","bcr81_100km.shp", "bcr82_100km.shp","bcr83_100km.shp")
canada <- shapefile("E:/GIS/basemaps/canadaLCC.shp")

#kNN biomass layers, 2011 (250-m)
b2011 <- list.files("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/",pattern="tif$")
x <- setwd("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) {bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
biomass<-bs2011
# writeRaster(bs2011, file=paste0(w,"bs2011_250m"))
# biomass <- brick(paste0(w,"bs2011_250m.grd"))

#Gaussian filter for 2011 biomass layers
bs2011_Gauss750 <- brick(paste(w,"bs2011_750.grd",sep=""))
names(bs2011_Gauss750) <- gsub("SpeciesGroups","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Species","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Structure","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Landcover","Landsc750",names(bs2011_Gauss750))
landscape <- bs2011_Gauss750
names(landscape)[1:4] <- gsub("v1","v1.1",names(landscape)[1:4])

#Adaptwest baseline climate variables
cur <- "E:/CMIP5/baseline19812010/"
# setwd(cur)
# clim <- list.files(cur, pattern =".asc$")
# curclim<-stack(clim)
# clim250 <- disaggregate(curclim,fact=c(4,4))
# writeRaster(clim250,file=paste0(cur,"curclim250"))
clim250 <- brick(paste0(cur,"curclim250.grd"))

# canroad <- shapefile(paste(w,"canroadlcc.shp",sep=""))
# roadlen <- raster(paste(w,"roadlen.asc",sep=""))
#road <- raster(paste(w,"roadonoff1.tif",sep=""))
# mr <- c(1, 2500000, 1,  NA, NA, 0)
# rcroad <- matrix(mr, ncol=3, byrow=TRUE)
# rrc <- reclassify(road,rcroad)
# writeRaster(rrc,file=paste(w,"road10",sep=""),overwrite=TRUE)
# rrc <- raster(paste(w,"road10.grd",sep=""))

#MODIS-based landcover (250-m)
# nalc2005 <- raster("E:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/NA_LandCover_2005_LCC.img")
# landcover <- brick(paste(w,"landcov_1km.grd",sep=""))
# dev750_1km <- landcover[[1]]
# dev750 <- resample(dev750_1km, nalc2005)
# led750_1km <- landcover[[2]]
# led750 <- resample(led750_1km, nalc2005)
# # mwat <- c(0, 17.1, 0,  17.9, 18.1, 1,  18.9, 20, 0)
# # rclwat <- matrix(mwat, ncol=3, byrow=TRUE)
# # water <- reclassify(nalc2005,rclwat)
# # #writeRaster(water,file="E:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/water.tif")
# # murb <- c(0, 16.1, 0,  16.9, 17.1, 1,  17.9, 20, 0)
# # rclurb <- matrix(murb, ncol=3, byrow=TRUE)
# # urban <- reclassify(nalc2005,rclurb)
# # writeRaster(urban,file="E:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/urban.tif")
# # mag <- c(0, 14.1, 0,  14.9, 15.1, 1,  15.9, 20, 0)
# # rclag <- matrix(mag, ncol=3, byrow=TRUE)
# # ag <- reclassify(nalc2005,rclag)
# #writeRaster(ag,file="E:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/ag.tif")
# # water <- raster("E:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/water.tif")
# # urban <- raster("E:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/urban.tif")
# # ag <- raster("E:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/ag.tif")
# # x <- stack(urban,ag)
# # urbag <- max(x)
# # fw750<-focalWeight(x=water,d=750,type="Gauss") #Gaussian filter with sigma=750 (tapers off around 2km)
# # dev750 <- focal(urbag,w=fw750,na.rm=TRUE)
# # led750 <- focal(water,w=fw750,na.rm=TRUE)
# landcov <- stack(dev750,led750,nalc2005)
# names(landcov) <- c("dev750","led750","nalc")
# writeRaster(landcov,file=paste(w,"landcov_250m",sep=""),overwrite=TRUE)
landcov <- brick(paste(w,"landcov_250m.grd",sep=""))

# Landform (100-m)
#dem <- raster("E:/GIS/topoedaphic/nadem100laz.tif")
# tri <- raster("E:/GIS/topoedaphic/tri250.tif")
# tpi <- raster("E:/GIS/topoedaphic/tpi250.tif")
# roughness <- raster("E:/GIS/topoedaphic/rough250.tif")
# slope <- raster("E:/GIS/topoedaphic/slope250.tif")
# lf <- raster("E:/GIS/topoedaphic/NA_Topo/lf_lcc1.tif")
# lf <- crop(lf,tri)
# lf1 <- resample(lf,tri,method="ngb")
# topog <- stack(tri,tpi,roughness,slope)
# topog <- addLayer(topog,lf1)
# names(topog) <- c("TRI","TPI","roughness","slope","lf")
# writeRaster(topog,file=paste0(w,"topography_250m"),overwrite=TRUE)
topog <-brick(paste(w,"topography_250m.grd",sep=""))


# road1k <- raster("E:/GIS/disturbance/VenterEtAlFootprint/RoadsLCC.tif")
# road <- resample(road1k,biomass[[1]],method="bilinear")
# writeRaster(road,file="E:/GIS/disturbance/VenterEtAlFootprint/Roads250m.tif",overwrite=TRUE)
road <- raster("E:/GIS/disturbance/VenterEtAlFootprint/Roads250m.tif")

setwd(w)
for (i in 1:length(bcrs2)){
  bcr <- shapefile(bcrs2[i])
  bcr <- crop(bcr,canada)
  biomass1 <- mask(crop(biomass,bcr),bcr)
  landscape1 <- mask(crop(landscape,bcr),bcr)
  clim1 <- mask(crop(clim250,biomass1[[1]]),biomass1[[1]])
  topo1 <- crop(topog,biomass1[[1]])
  topo2 <- resample(topo1,biomass1[[1]],method="ngb")
  topo3 <- mask(topo2,biomass1[[1]])
  landcov1 <- crop(landcov,biomass1[[1]])
  landcov2 <- resample(landcov1,biomass1[[1]], method="ngb")
  landcov3 <- mask(landcov2,biomass1[[1]])
  ROAD <- crop(extend(road,biomass1[[1]]),biomass1[[1]])
  ROAD <- mask(ROAD,biomass1[[1]])
  bcrr <- rasterize(bcr,biomass1[[1]])
  #bcrc <- mask(crop(bcrr,biomass1[[1]]),biomass1[[1]])
  bs <- stack(bcrr)
  names(bs) <- "bcr"
  for (j in 1:nlayers(clim1)) {
    bs <- addLayer(bs, clim1[[j]])}
  for (j in 1:nlayers(biomass1)) {
    bs <- addLayer(bs, biomass1[[j]])}
  for (j in 1:nlayers(landscape1)) {
    bs <- addLayer(bs, landscape1[[j]])}
  for (j in 1:nlayers(topo3)) {
    bs <- addLayer(bs, topo3[[j]])}
  for (j in 1:nlayers(landcov3)) {
    bs <- addLayer(bs, landcov3[[j]])}
  bs <- stack(bs,ROAD)
  names(bs)[nlayers(bs)] <- "ROAD"
  writeRaster(bs,file=paste(w,gsub("_100km.shp","",bcrs2[i]),"all_250m",sep=""),overwrite=TRUE)
}


setwd(w)
for (i in 1:length(bcrs)){
  bcr <- shapefile(bcrs[i])
  bcr <- crop(bcr,canada)
  biomass1 <- mask(crop(biomass,bcr),bcr)
  landscape1 <- mask(crop(landscape,bcr),bcr)
  clim1 <- mask(crop(clim250,biomass1[[1]]),biomass1[[1]])
  topo1 <- crop(topog,biomass1[[1]])
  topo2 <- resample(topo1,biomass1[[1]],method="ngb")
  topo3 <- mask(topo2,biomass1[[1]])
  landcov1 <- crop(landcov,biomass1[[1]])
  landcov2 <- resample(landcov1,biomass1[[1]], method="ngb")
  landcov3 <- mask(landcov2,biomass1[[1]])
  ROAD <- crop(extend(road,biomass1[[1]]),biomass1[[1]])
  ROAD <- mask(ROAD,biomass1[[1]])
  bcrr <- rasterize(bcr,biomass1[[1]])
  #bcrc <- mask(crop(bcrr,biomass1[[1]]),biomass1[[1]])
  bs <- stack(bcrr)
  names(bs) <- "bcr"
  for (j in 1:nlayers(clim1)) {
    bs <- addLayer(bs, clim1[[j]])}
  for (j in 1:nlayers(biomass1)) {
    bs <- addLayer(bs, biomass1[[j]])}
  for (j in 1:nlayers(landscape1)) {
    bs <- addLayer(bs, landscape1[[j]])}
  for (j in 1:nlayers(topo3)) {
    bs <- addLayer(bs, topo3[[j]])}
  for (j in 1:nlayers(landcov3)) {
    bs <- addLayer(bs, landcov3[[j]])}
  bs <- stack(bs,ROAD)
  names(bs)[nlayers(bs)] <- "ROAD"
  writeRaster(bs,file=paste(w,gsub("_100km.shp","",bcrs[i]),"all_250m",sep=""),overwrite=TRUE)
}

