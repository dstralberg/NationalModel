library(raster)
library(gbm)
library(maptools)
library(dplyr)
library(sf)

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
w <-"G:/Boreal/NationalModelsV2/"
x <-"G:/Boreal/NationalModelsV2/Nov2019/"


#generate predictions and plots from models
brtplot <- function (model,bs,bcrname) {
  load(model)
  spec <- substr(model,1,4)
  rast <- raster::predict(bs, out, type="response", n.trees=out$n.trees)
  writeRaster(rast, filename=paste(x,bcrname,"/",spec,"_", bcrname, "_pred1km1",sep=""), format="GTiff",overwrite=TRUE)
}

#mosaic predictions
predmos <- function(spec,bcrs) {
  for (i in 1:length(bcrs)){
    bcrname <- gsub("_100km.shp","",bcrs[i])
    xx <- try(pred <- raster(paste(x,bcrname,"/",spec,"_", bcrname, "_pred1km1.tif",sep="")))
    if(class(xx)!="try-error"){
      if(i==1) {mospred <- pred}
      if(i>1) {mospred <- mosaic(mospred,pred,fun=mean)}
    }
  }
  return(mospred)
}

bcrs <- c("bcr4_100km.shp","bcr5_100km.shp","bcr9_100km.shp","bcr10_100km.shp","bcr11_100km.shp","bcr12_100km.shp","bcr13_100km.shp","bcr14_100km.shp","bcr60_100km.shp","bcr61_100km.shp","bcr70_100km.shp","bcr71_100km.shp","bcr80_100km.shp","bcr81_100km.shp", "bcr82_100km.shp","bcr83_100km.shp")

# for (i in 1:length(bcrs)){
#   bcr <- shapefile(paste(w,bcrs[i],sep=""))
#   bcrname <- gsub("_100km.shp","",bcrs[i])
#   dir.create(paste(x,bcrname,"/",sep=""))
# }

for (i in 1:length(bcrs)){
  bcr <- shapefile(paste(w,bcrs[i],sep=""))
  bcrname <- gsub("_100km.shp","",bcrs[i])
  bs <- brick(paste(w,bcrname,"all_1km.grd",sep=""))
  ARU <- bs[[1]]*0
  bs <- addLayer(bs,ARU)
  names(bs)[nlayers(bs)]<-"ARU"
  setwd(paste(x,bcrname,"/",sep=""))
  models <- list.files(pattern=".RData")
  if(length(models)>0){
    for (j in 1:length(models)) {
      brtplot(models[j],bs,bcrname)
    }
    gc()
  }
}

CAWApred <- predmos("CAWA",bcrs)
writeRaster(CAWApred,file=paste(x,"mosaicpred","CAWApred.tif",sep=""),overwrite=TRUE)
OSFLpred <- predmos("OSFL",bcrs)
writeRaster(OSFLpred,file=paste(x,"mosaicpred","OSFLpred.tif",sep=""),overwrite=TRUE)
