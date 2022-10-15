library(raster)
library(ggplot2)
library(maptools)
library(dplyr)

w <-"H:/My Drive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModelsV4.0/Feb2020/artifacts/"
x <-"H:/My Drive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModelsV4.0/Feb2020/website/"
setwd(w)
specpred <- list.dirs(w, full.names=FALSE)
bcr <- shapefile("G:/Boreal/NationalModelsV2/BCRSubunits.shp")
bcrlu <- unique(as.data.frame(bcr[,c(1:2)]))
names(bcrlu) <- c("BCR","BCRchar")
#canada <- shapefile("E:/GIS/basemaps/canadaLCC.shp")
#bcrcan <- crop(bcr,canada)
setwd(w)
models <- list.files(paste0(w,"ALFL/"),pattern=".tif")
a <- raster(paste0(w,"ALFL/",models[1]))
# lcr <- resample(lc,a,method='ngb')
# writeRaster(lcr,file="G:/Boreal/NationalModelsV2/lcr.tif",overwrite=TRUE)
# bcrr <- rasterize(bcrcan,lcr,field="BCR")
# writeRaster(bcrr,file="G:/Boreal/NationalModelsV2/bcrr.tif",overwrite=TRUE)
landcov <- raster("G:/Boreal/NationalModelsV2/lcr.tif")
lc <- crop(landcov,a)
units <- rasterize(bcr,a,field="BCR")
#bc <- crop(units,a)
lu <- read.csv("G:/Boreal/NationalModelsV2/landcov.csv")
names(lu) <- c("nalc","landcover")

bbcr <- shapefile("I:/My Drive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/BCR_BAMSubunits.shp")
subspp <- shapefile("H:/My Drive/BAM.SharedDrive/RshProjs/PopnStatus/GCTH/C_m_minimus_range_v6_8Apr2022_NoSmallNSIslands_LCC.shp")
subsppr <- trim(rasterize(subspp,a))
m <- c(0, 0, 0,  1, 10000, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
range <- as.integer(reclassify(subsppr, rclmat))
 

sumdens <- function(species,landcov,units){
  setwd(paste0(w,species,"/"))
  models <- list.files(paste0(w,species,"/"),pattern=".tif")
  models <- grep("-SD.tif",models,invert=TRUE,value=TRUE) 
  bcrv <- getValues(units)
  lcv <- getValues(lc)
  bcrlc <- as.data.frame(cbind(bcrv,lcv))
  bcrlc$count <- 1
  names(bcrlc) <- c("BCR","nalc","count")
  lc1 <- aggregate(bcrlc$count,by=list(bcrlc$BCR,bcrlc$nalc),FUN='sum',na.rm=TRUE)
  names(lc1) <- c("BCR","nalc","area")
  bcr1 <- aggregate(bcrlc$count,by=list(bcrlc$BCR),FUN='sum',na.rm=TRUE)
  names(bcr1) <- c("BCR","barea")
  lc1 <- merge(lc1,bcr1,by="BCR")
  lc1$areaprop <- lc1$area/lc1$barea
  lc1<-lc1[lc1$BCR>0,]
  lc1 <- lc1[lc1$nalc %in% list(1,2,5,6,8,10,11,12,14,15),]
  lc1 <- merge(lc1,bcrlu,by="BCR")
  denstable <- data.frame("spec"="","BCR"=0,"BCRname"="","nalc"=0,"landcover"="","area"=0,"barea"=0,"areaprop"=0,"meandens"=0,"Q5"=0,"Q50"=0,"Q95"=0)[0,]
  
  for (i in 1:32){
    rast <- raster(models[i])
    ss <- unlist(strsplit(gsub(".tif","",models[i]),"-"))
    spec <- paste0(ss[2],ss[4],ss[5])    
    pred <- getValues(rast)
    pred1 <- as.data.frame(cbind(bcrv,lcv,pred))
    #pred1 <- pred1[pred1$bcrv>0,]
    densmean <- aggregate(pred1$pred,by=list(pred1$bcrv,pred1$lcv),FUN='mean',na.rm=TRUE)
    names(densmean) <- c("BCR","nalc","mean")
    densquant <- aggregate(pred1$pred,by=list(pred1$bcrv,pred1$lcv),FUN=function(x){quantile(x,c(0.05,0.5,0.95),na.rm=TRUE)})
    densmean <- cbind(densmean,unlist(densquant$x))
    names(densmean)[4:6] <- c("Q5","Q50","Q95")
    densmean$spec <- spec
    densmean <- merge(densmean,lc1,by=c("BCR","nalc"))
    densmean <- merge(densmean,lu,by="nalc")
    denstable <- rbind(denstable,densmean)
  }
  return(denstable)
}

pop <- function(species,range){
  setwd(paste0(w,species,"/"))
  models <- list.files(paste0(w,species,"/"),pattern=".tif")
  models <- grep("-SD.tif",models,invert=TRUE,value=TRUE) 
  poptab <- data.frame("spec"="","pop"=0)
  

  for (i in 1:32){
    rast <- raster(models[i])
    ss <- unlist(strsplit(gsub(".tif","",models[i]),"-"))
    spec <- paste0(ss[2],ss[4],ss[5])  
    rr <- crop(rast,range)
    rrm <- rr*range
    poptab[i,1]= spec
    poptab[i,2]= cellStats(rrm,stat='sum')*100
  }
    
  return(poptab)
}

species<-"GCTH"
poptab <- pop(species,range)
q <- quantile(poptab$pop, probs=c(0.05,0.5,0.95))


  
write.csv(poptab,file=paste0("H:/My Drive/BAM.SharedDrive/RshProjs/PopnStatus/GCTH/",species,"_population_subspp.csv"),row.names=FALSE)

    
