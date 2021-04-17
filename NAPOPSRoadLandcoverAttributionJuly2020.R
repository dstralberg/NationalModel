library(raster)
library(geosphere)
library(sp)
library(rgdal)
require(rgeos)
load("E:/BAM/BAMData/BAM_data_package_November2019.RData")
coordinates(SScombo) <- c("X", "Y")

proj4string(SScombo)<-LCC
LAEA <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
DD <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SS_LAEA <- spTransform(SScombo, LAEA)
SS_DD <- spTransform(SScombo, DD)
SSDD <- as.data.frame(SS_DD)

Canada <- shapefile("E:/GIS/basemaps/CanadaLAEA.shp")
SS_LAEA <- crop(SS_LAEA,Canada)

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

lcvals5 <- extract(stack5,SS_LAEA)
lcvals3 <- extract(stack3,SS_LAEA)

canroad <- shapefile("E:/GIS/transportation/CanadaRoads2019_lite/CanadaRoads2019_lite_DD.shp")
roaddist <- raster("E:/GIS/transportation/CanadaRoads2019_lite/canroaddist")
roadvals <- extract(roaddist,SS_LAEA)

vals <- cbind(SS_LAEA@data,"roaddist"=roadvals[,1],lcvals3,lcvals5)
write.csv(vals,file="E:/BAM/BAMData/vals.csv", row.names=FALSE)

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

