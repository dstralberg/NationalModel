library(raster)
w <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/artifacts/"
x <- "I:/My Drive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/website/"
z<- "G:/Boreal/NationalModelsV2/KBAs/"
setwd(x)
speclist <- read.csv("speclist.csv")

bcr <- rgdal::readOGR("E:/GIS/basemaps/BCRs/bcrfinallcc.shp")
models <- list.files(paste0(w,speclist[1,1],"/"),pattern="Mean.tif$")
rast <- raster(paste0(w,speclist[1,1],"/",models[1]))
bcrc <- crop(bcr,rast)

clip1 <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.001)
  zmin <- min(zmin,0.01)
  m <- c(0, zmin, 0)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  rc <- reclassify(rast, rclmat)
  writeRaster(rc,file=paste0(z,spec,"_clip1"),format="GTiff",overwrite=TRUE)
}

clip2 <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.005)
  zmin <- min(zmin,0.05)
  m <- c(0, zmin, 0)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  rc <- reclassify(rast, rclmat)
  writeRaster(rc,file=paste0(z,spec,"_clip2"),format="GTiff",overwrite=TRUE)
}


clip3 <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  zmin <- min(prev,0.001)
  m <- c(0, zmin, 0)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  rc <- reclassify(rast, rclmat)
  writeRaster(rc,file=paste0(z,spec,"_clip3"),format="GTiff",overwrite=TRUE)
}

options(warn = -1) 

for (i in 1:nrow(speclist)) {
  models <- list.files(paste0(w,speclist[i,1],"/"),pattern="Mean.tif$")
  rast <- raster(paste0(w,speclist[i,1],"/",models[1]))
  spec <- substr(models[1],6,9)
  map <- speclist[i,2]
  if (map==1) { 
    clip1(rast,spec)
  } else if (map==2) {
    clip2(rast,spec)
  } else {
    clip3(rast,spec)
  }
}