library(raster)
library(dismo)
library(gbm)
library(maptools)
library(dplyr)
library(sf)

bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.5)
provstate <- rgdal::readOGR("F:/GIS/basemaps/province_state_line.shp")
lakes <- rgdal::readOGR("F:/GIS/hydrology/lakes_lcc.shp")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
w <-"G:/Boreal/NationalModelsV2/"


#generate predictions and plots from models
brtplot <- function (j,bs,bcrname,landcov) {
  load(models[j])
  spec <- substr(models[j],1,4)
  #varimp <- as.data.frame(out$contributions)
  #write.csv(varimp,file=paste(w,speclist[j],"varimp3.csv",sep=""))
  #cvstats <- t(as.data.frame(out$cv.statistics))
  #write.csv(cvstats,file=paste(w,speclist[j],"cvstats3.csv",sep=""))
  #pdf(paste(spec,"_plot1.pdf",sep=""))
  #gbm.plot(out,n.plots=12,smooth=TRUE)
  #dev.off()
  rast <- raster::predict(bs, out, type="response", n.trees=out$n.trees)
  writeRaster(rast, filename=paste(spec,"_", bcrname, "_pred1km1",sep=""), format="GTiff",overwrite=TRUE)
  densities <- zonal(rast,landcov,fun='mean')
  write.csv(densities,file=paste(w,spec,"_",bcrname,"_density_by_landcover.csv",sep=""),row.names=FALSE)
  
  #q99 <- quantile(rast, probs=c(0.99))	
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  p <- crop(provstate,bcr)
  l <- crop(lakes,bcr)
  png(file=paste(spec,"_", bcrname, "_pred1km1.png",sep=""), height=600, width=500)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(spec),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(spec), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.55,0.8,0.82,0.87), axis.args=list(cex.axis=1.3))
  plot(p, col="gray", add=TRUE)
  plot(l, col="gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="dark gray", add=TRUE)
  text(-300000,8900000,"Potential density (males/ha)", cex=1.3)
  dev.off()
  
}


bcrs <- c("bcr4_100km.shp","bcr5_100km.shp","bcr6_100km.shp","bcr7_100km.shp","bcr8_100km.shp","bcr9_100km.shp","bcr10_100km.shp","bcr11_100km.shp","bcr12_100km.shp","bcr13_100km.shp","bcr14_100km.shp")
for (i in 1:length(bcrs)){
  bcr <- shapefile(paste(w,bcrs[i],sep=""))
  bcrname <- gsub("_100km.shp","",bcrs[i])
  bs <- stack(paste(w,bcrname,"_1km.grd",sep=""))
  bs2 <- stack(paste(w,bcrname,"cat_1km.grd",sep=""))
  bs <- stack(bs,bs2)
  landcov <- bs[[148]]
  setwd(paste(w,"March2019","/",bcrname,"/",sep=""))
  models <- list.files(pattern=".RData")
  for (j in 1:length(models)) {
    brtplot(j,bs,bcrname,landcov)
    }
  gc()
}



