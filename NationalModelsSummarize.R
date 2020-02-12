library(raster)
library(ggplot2)
library(maptools)
library(dplyr)

w <-"F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/feb2020/artifacts/"
#x <-"F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus//NationalModels/Nov2019/mosaics/"
setwd(w)
specpred <- list.dirs(w, full.names=FALSE)
bcr <- shapefile("E:/GIS/basemaps/BCRs/bcrfinallcc.shp")
bcrlu <- unique(as.data.frame(bcr[,c(5,7)]))
names(bcrlu) <- c("BCR","BCRname")
canada <- shapefile("E:/GIS/basemaps/canadaLCC.shp")
bcrcan <- crop(bcr,canada)
landcov <- raster("E:/GIS/landcover/NALC/LandCover_IMG/NA_LandCover_2005/data/NA_LandCover_2005/NA_LandCover_2005_LCC.img")
setwd(w)
a<- raster(specpred[1])
# lc <- crop(landcov,a)
# lcr <- resample(lc,a,method='ngb')
# writeRaster(lcr,file="G:/Boreal/NationalModelsV2/lcr.tif",overwrite=TRUE)
# bcrr <- rasterize(bcrcan,lcr,field="BCR")
# writeRaster(bcrr,file="G:/Boreal/NationalModelsV2/bcrr.tif",overwrite=TRUE)
landcov <- raster("G:/Boreal/NationalModelsV2/lcr.tif")
units <- raster("G:/Boreal/NationalModelsV2/bcrr.tif")
lu <- read.csv("G:/Boreal/NationalModelsV2/landcov.csv")
names(lu) <- c("nalc","landcover")


sumdens <- function(species,landcov,units){
  setwd(paste0(w,species,"/"))
  models <- list.files(paste0(w,species,"/"),pattern=".tif")
  models <- grep("-SD.tif",models,invert=TRUE,value=TRUE) 
  bcrv <- getValues(units)
  lcv <- getValues(landcov)
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
  #write.csv(lc1,file=paste(w,"bcrlandcov.csv",sep=""))
  denstable <- data.frame("spec"="","BCR"=0,"BCRname"="","nalc"=0,"landcover"="","area"=0,"barea"=0,"areaprop"=0,"meandens"=0,"Q5"=0,"Q50"=0,"Q95"=0)[0,]
  
  for (i in 1:length(models)){
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
    #write.csv(densmean,file=paste0(w,specpred,"/",spec,"_densities.csv"),row.names=FALSE)
    denstable <- rbind(denstable,densmean)
  }
  write.csv(denstable,file=paste0(w,species,"/",species,"_densities.csv"),row.names=FALSE)
}


plotdens <- function (models,landcov,units){
for (i in 1:length(models)){
  spec <- paste0(substr(models[i],6,9),substr(models[i],15,20))  
  densmean <- read.csv(paste(w,spec,"_densities.csv",sep=""))
  png(filename=paste0(w,specpred,"/",spec,"_densityplot.png"),width=1800,height=1800,res=216)
  p<-ggplot(densmean,aes(x=landcover,y=mean))+
    geom_bar(aes(fill=landcover),width=densmean$areaprop*2, stat="identity")+
    xlab("Land cover type")+
    ylab("Mean density (males/ha)")+
    theme(legend.position = "none")+
    theme(axis.text.x=element_text(angle=90,vjust=0.4))+
    ylim(0,0.5*max(densmean$mean,na.rm=TRUE))+
    facet_wrap(~BCR,ncol=4)
  print(p)
  dev.off()
}
}
  
boxplotdens <- function(models,landcov,units){
bcrv <- getValues(units)
lcv <- getValues(landcov)
for (i in 1:length(models)){
  rast <- raster(models[i])
  spec <- paste0(substr(models[i],6,9),substr(models[i],15,20))  
  pred <- getValues(rast)
  pred1 <- as.data.frame(cbind(bcrv,lcv,pred))
  pred1 <- pred1[pred1$bcrv>0,]
  pred1 <- na.omit(pred1)
  names(pred1) <- c("BCR","nalc","pred")
  #pred1$lcv <- as.factor(pred1$lcv)
  pred1 <- left_join(pred1,lu)
  png(filename=paste0(w,specpred,"/",spec,"_boxplot.png",sep=""),width=2200,height=1800,res=216)
  p<-ggplot(pred1,aes(x=landcover,y=pred))+
    geom_boxplot()+
    xlab("Land cover type")+
    ylab("Mean density (males/ha)")+
    theme(legend.position = "none")+
    theme(axis.text.x=element_text(angle=90,vjust=0.4))+
    facet_wrap(~BCR,ncol=4)
  print(p)
  dev.off()
}
}

for (j in 2:length(specpred)) {
  species<-specpred[j]
  d <- sumdens(species,landcov,units)
  #boxplotdens(species,landcov,units)
  #plotdens(species,landcov,units)
}


