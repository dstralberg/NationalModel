library(raster)
library(geosphere)
library(sp)
library(rgdal)
require(rgeos)
#load("E:/BAM/BAMData/BAM_data_package_November2019.RData")
SScombo <- read.csv("I:/My Drive/BAM.SharedDrive/DataStuff/Avian.Data/BAM-V6-Use/BAMXYFeb2021.csv")
SScombo <- na.omit(SScombo[,c(1:6,12)])
coordinates(SScombo) <- c("longitude", "latitude")

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LAEA <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
DD <- CRS("+proj=longlat +datum=NAD83 +ellps=GRS80 +towgs84=0,0,0")
proj4string(SScombo)<-DD
SS_LAEA <- spTransform(SScombo, LAEA)
SS_LCC <- spTransform(SScombo, LCC)
SSDD <- as.data.frame(SScombo)

Canada <- shapefile("E:/GIS/basemaps/CanadaLAEA.shp")
SS_LAEA_Canada <- crop(SS_LAEA,Canada)

Alaska <- shapefile("E:/GIS/basemaps/AlaskaLAZEA.shp")
SS_LAEA_AK <- crop(SS_LAEA,Alaska)

Lower48 <- shapefile("E:/GIS/basemaps/Lower48LAZEA.shp")
SS_LAEA_Lower48 <- crop(SS_LAEA,Lower48)

f5a <- "E:/GIS/landcover/CECLandsat/FocalStats_ASK_5x5/"
setwd(f5a)
lc5a <- list.files(path=f5a, pattern=".tif$")
stack5a <- raster(lc5a[1])
for (i in 2:length(lc5a)) {
  rast <- raster(lc5a[i])
  stack5a <- addLayer(stack5a,rast)
}

f3a <- "E:/GIS/landcover/CECLandsat/FocalStats_ASK_3x3/"
setwd(f3a)
lc3a <- list.files(path=f3a, pattern=".tif$")
stack3a <- raster(lc3a[1])
for (i in 2:length(lc3a)) {
  rast <- raster(lc3a[i])
  stack3a <- addLayer(stack3a,rast)
}

f5 <- "E:/GIS/landcover/CECLandsat/FocalStats_5x5/"
setwd(f5)
lc5 <- list.files(path=f5, pattern=".tif$")
stack5 <- raster(lc5[1])
for (i in 2:length(lc5)) {
  rast <- raster(lc5[i])
  stack5 <- addLayer(stack5,rast)
}

f3 <- "E:/GIS/landcover/CECLandsat/FocalStats_3x3/"
setwd(f3)
lc3 <- list.files(path=f3, pattern=".tif$")
stack3 <- raster(lc3[1])
for (i in 2:length(lc3)) {
  rast <- raster(lc3[i])
  stack3 <- addLayer(stack3,rast)
}

f5u <- "E:/GIS/landcover/CECLandsat/FocalStats_USA_5x5/"
setwd(f5u)
lc5u <- list.files(path=f5u, pattern=".tif$")
stack5u <- raster(lc5u[1])
for (i in 2:length(lc5u)) {
  rast <- raster(lc5u[i])
  stack5u <- addLayer(stack5u,rast)
}

f3u <- "E:/GIS/landcover/CECLandsat/FocalStats_USA_3x3/"
setwd(f3u)
lc3u <- list.files(path=f3u, pattern=".tif$")
stack3u <- raster(lc3u[1])
for (i in 2:length(lc3u)) {
  rast <- raster(lc3u[i])
  stack3u <- addLayer(stack3u,rast)
}

lcvals5 <- extract(stack5,SS_LAEA_Canada)
lcvals3 <- extract(stack3,SS_LAEA_Canada)

lcvals5a <- extract(stack5a,SS_LAEA_AK)
lcvals3a <- extract(stack3a,SS_LAEA_AK)

lcvals5u <- extract(stack5a,SS_LAEA_Lower48)
lcvals3u <- extract(stack3a,SS_LAEA_Lower48)

AKroaddist <- raster("E:/GIS/transportation/USA Roads_oneFC/AKroaddist.tif")
canroaddist <- raster("E:/GIS/transportation/CanadaRoads2019_lite/canroaddist1.tif")

roadvals <- as.data.frame(extract(canroaddist,SS_LAEA_Canada))
roadvalsa <- as.data.frame(extract(AKroaddist,SS_LAEA_AK))
roadvalsu <- extract(usroaddist,SS_LAEA_Lower48)

valsa <- cbind(SS_LAEA_AK@data,"roaddist"=roadvalsa[,1],lcvals3a,lcvals5a)
write.csv(valsa,file="E:/BAM/BAMData/AKvals.csv", row.names=FALSE)

vals <- cbind(SS_LAEA_Canada@data,"roaddist"=roadvals[,1],lcvals3,lcvals5)
write.csv(vals,file="E:/BAM/BAMData/Canadavals.csv", row.names=FALSE)

test <- shapefile("E:/BAM/BAMData/AKN/pointsTestPhoney.shp")
test5 <- extract(stack5,test)
test3 <- extract(stack3,test)
testroad <- extract(roaddist,test)
testvals <- cbind(test@data,"roaddist"=testroad[,1],test5,test3)
write.csv(testvals,file="E:/BAM/BAMData/AKN/testvals.csv", row.names=FALSE)

#usroad <- shapefile("E:/GIS/transportation/USA Roads_oneFC/USA_Roads_oneFCt.shp")



#write.csv(SScombo, file="E:/BAM/BAMData/BAM_SS_November2019.csv")
#write.csv(SS_LAEA, file="E:/BAM/BAMData/BAM_SS_November2019_LAEA.csv")
#canroad <- shapefile("E:/GIS/transportation/CanadaRoads2019_lite/CanadaRoads2019_lite_LAEA.shp")
#write.csv(lcvals3,file="E:/GIS/landcover/CECLandsat/FocalStats_3x3/lcvals3x3.csv", row.names=FALSE)
#write.csv(lcvals5,file="E:/GIS/landcover/CECLandsat/FocalStats_5x5/lcvals5x5.csv", row.names=FALSE)
#dd = gDistance(canroad, as(SS_LAEA,"SpatialPoints"), byid=TRUE)
#roadclass <- raster("E:/GIS/transportation/CanadaRoads2019_lite/canroadclass2")
#roadvals2 <- dist2Line(SS_DD, canroad)

