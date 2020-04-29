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
w <-"G:/Boreal/NationalModelsV2/"

bcrs <- c("bcr4_100km.shp","bcr5_100km.shp","bcr9_100km.shp","bcr10_100km.shp","bcr11_100km.shp","bcr12_100km.shp","bcr13_100km.shp","bcr14_100km.shp")
bcrs2 <- c("bcr60_100km.shp","bcr61_100km.shp","bcr70_100km.shp","bcr71_100km.shp","bcr80_100km.shp","bcr81_100km.shp", "bcr82_100km.shp","bcr83_100km.shp")
canada <- shapefile("E:/GIS/basemaps/canadaLCC.shp")

#Adaptwest baseline climate variables
cur <- "E:/CMIP5/baseline19812010/"
setwd(cur)
clim <- list.files(cur, pattern =".asc$")
curclim<-stack(clim)

rrc <- raster(paste(w,"road10.grd",sep=""))

#MODIS-based landcover (250-m)
# for (x in 2001:2011){ 
  nalcx <- raster(paste0("E:/GIS/landcover/LCTS_2000-2011/LCTS19c",x,".tif"))
  nalc <- resample(nalcx,curclim,method="ngb")
  writeRaster(nalc,paste0(w,"nalc-1km_2001",x,".tif",sep=""))
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
  landcov <- stack(dev750,led750)
  landcover <- resample(landcov,curclim)
  landcover <- addLayer(landcover,nalc)
  names(landcover) <- c("dev750","led750","nalc")
  writeRaster(landcover,file=paste0(w,"landcov_1km_2001",sep=""),overwrite=TRUE)
# }
# 
landcover <- brick(paste(w,"landcov_1km_2001.grd",sep=""))

#kNN biomass layers, 2001 (250-m)
# b2001 <- list.files("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2001/",pattern="tif$")
# setwd("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2001/")
# bs2001 <- stack(raster(b2001[1]))
# for (i in 2:length(b2001)) {bs2001 <- addLayer(bs2001, raster(b2001[i]))}
# names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
# 
# climcan <- crop(curclim,bs2001[[1]])
# biomass <- resample(bs2001,climcan)
# writeRaster(biomass, file=paste0(w,"bs2001_1km"))
# writeRaster(bs2001, file=paste0(w,"bs2001_250m"))

biomass <- brick(paste(w,"bs2001_1km.grd",sep=""))

#Gaussian filter for 2001 biomass layers
# bs2001_Gauss750 <- brick(paste(w,"bs2001_750.grd",sep=""))
# names(bs2001_Gauss750) <- gsub("SpeciesGroups","Landsc750",names(bs2001_Gauss750))
# names(bs2001_Gauss750) <- gsub("Species","Landsc750",names(bs2001_Gauss750))
# names(bs2001_Gauss750) <- gsub("Structure","Landsc750",names(bs2001_Gauss750))
# names(bs2001_Gauss750) <- gsub("LandCover","Landsc750",names(bs2001_Gauss750))
# landscape <- resample(bs2001_Gauss750,biomass)
# writeRaster(landscape, file=paste(w,"bs2001_750_1km",sep=""))
landscape <- brick(paste(w,"bs2001_750_1km.grd",sep=""))
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

road1k <- raster("E:/GIS/disturbance/VenterEtAlFootprint/RoadsLCC.tif")

setwd(w)
for (i in 1:length(bcrs)){
  bcr <- shapefile(bcrs[i])
  bcr <- crop(bcr,canada)
  biomass1 <- mask(crop(biomass,bcr),bcr)
  landscape1 <- mask(crop(landscape,bcr),bcr)
  clim1 <- mask(crop(curclim,biomass1[[1]]),biomass1[[1]])
  landcov1 <- mask(crop(landcover,biomass1[[1]]),biomass1[[1]])
  topo1 <- mask(crop(topog,biomass1[[1]]),biomass1[[1]])
  bcrr <- rasterize(bcr,clim1[[1]])
  #bcrc <- mask(crop(bcrr,biomass1[[1]]),biomass1[[1]])
  bs <- stack(bcrr)
  names(bs) <- "bcr"
  for (j in 1:nlayers(clim1)) {
    bs <- addLayer(bs, clim1[[j]])}
  for (j in 1:nlayers(biomass1)) {
    bs <- addLayer(bs, biomass1[[j]])}
  for (j in 1:nlayers(landscape1)) {
    bs <- addLayer(bs, landscape1[[j]])}
  for (j in 1:nlayers(landcov1)) {
    bs <- addLayer(bs, landcov1[[j]])}
  for (j in 1:nlayers(topo1)) {
    bs <- addLayer(bs, topo1[[j]])}
  rrc1 <- crop(extend(road1k,biomass1[[1]]),biomass1[[1]])
  ROAD <- resample(rrc1,biomass1[[1]],method="ngb")
  ROAD <- mask(ROAD,biomass1[[1]])
  bs <- stack(bs,ROAD)
  names(bs)[nlayers(bs)] <- "ROAD"
  writeRaster(bs,file=paste(w,gsub("_100km.shp","",bcrs[i]),"all_2001_1km",sep=""),overwrite=TRUE)
}

setwd(w)
for (i in 1:length(bcrs2)){
  bcr <- shapefile(bcrs2[i])
  bcr <- crop(bcr,canada)
  biomass1 <- mask(crop(biomass,bcr),bcr)
  landscape1 <- mask(crop(landscape,bcr),bcr)
  clim1 <- mask(crop(curclim,biomass1[[1]]),biomass1[[1]])
  landcov1 <- mask(crop(landcover,biomass1[[1]]),biomass1[[1]])
  topo1 <- mask(crop(topog,biomass1[[1]]),biomass1[[1]])
  bcrr <- rasterize(bcr,clim1[[1]])
  #bcrc <- mask(crop(bcrr,biomass1[[1]]),biomass1[[1]])
  bs <- stack(bcrr)
  names(bs) <- "bcr"
  for (j in 1:nlayers(clim1)) {
    bs <- addLayer(bs, clim1[[j]])}
  for (j in 1:nlayers(biomass1)) {
    bs <- addLayer(bs, biomass1[[j]])}
  for (j in 1:nlayers(landscape1)) {
    bs <- addLayer(bs, landscape1[[j]])}
  for (j in 1:nlayers(landcov1)) {
    bs <- addLayer(bs, landcov1[[j]])}
  for (j in 1:nlayers(topo1)) {
    bs <- addLayer(bs, topo1[[j]])}
  rrc1 <- crop(extend(road1k,biomass1[[1]]),biomass1[[1]])
  ROAD <- resample(rrc1,biomass1[[1]],method="ngb")
  ROAD <- mask(ROAD,biomass1[[1]])
  bs <- stack(bs,ROAD)
  names(bs)[nlayers(bs)] <- "ROAD"
  writeRaster(bs,file=paste(w,gsub("_100km.shp","",bcrs2[i]),"all_2001_1km",sep=""),overwrite=TRUE)
}


