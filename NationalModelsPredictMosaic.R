library(raster)
library(dismo)
library(rgdal)
library(gbm)
library(maptools)
library(dplyr)
library(sf)
library(Matrix)
library(reshape2)
library(colorspace)

bluegreen.colors <- colorRampPalette(c("#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.8)
bgtrunc <- colorRampPalette(c("#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=10)

bgy <- sequential_hcl(10, "ag_GrnYl",rev=TRUE)

p<- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("E:/GIS/hydrology/lakes_lcc.shp")
bcr <- rgdal::readOGR("E:/GIS/basemaps/BCRs/bcrfinallcc.shp")
canada <- rgdal::readOGR("E:/GIS/basemaps/canadaLCC.shp")
natureserve <- "E:/GIS/NatureServe/Abbreviated/_lcc/"
LCC <- CRS(projection(canada))
w <-"F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/artifacts/"
specpred <- list.dirs(w, full.names=FALSE)

load("E:/BAM/BAMData/BAM_data_package_November2019.RData")

#summarize mean densities
prev <- data.frame("species"=specpred[2:length(specpred)],"meandens"=0)
for (i in 2:length(specpred)){
  models <- list.files(paste0(w,specpred[i],"/"),pattern="Mean.tif$")
  rast <- raster(paste0(w,specpred[i],"/",models[1]))
  spec <- substr(models[1],6,9)
  meandens <- cellStats(rast, 'mean')
  prev$meandens[i-1] = meandens
}
write.csv(prev,file=paste0(w,"meandensities.csv"),row.names=FALSE)

#generate maps
# brtplot <- function (rast,spec) {
#   prev <- cellStats(rast, 'mean')	
#   q99 <- quantile(rast, probs=c(0.99))
#   max <- max(3*prev,q99)
#   png(file=paste0(w,spec,"/",spec,"_pred1km1.png"), width=2600, height=1600, res=216)
#   par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0), xpd=TRUE)
#   par(mar=c(0,0,0,0))
#   plot(rast, col="blue", axes=FALSE, legend=FALSE)
#   plot(rast, col=bluegreen.colors(15), maxpixels=5000000, zlim=c(0,max), axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
#   plot(p, col="gray", add=TRUE)
#   plot(l, col="gray", border=NA,add=TRUE)
#   plot(bcr, col=NA, border="dark gray", add=TRUE)
#   #text(2200000,9400000,"Potential density (males/ha)", cex=1.3)
#   dev.off()
# }

brtplot <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  min <- max(prev,0.01)
  min <- min(min,0.05)
  #q99 <- quantile(rast, probs=c(0.99))
  #max <- max(3*prev,q99)
  max <- cellStats(rast, 'max')
  png(file=paste0(w,spec,"/",spec,"_pred1km2.png"), width=2600, height=1600, res=216)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0), xpd=TRUE)
  par(mar=c(0,0,0,0))
  plot(rast, col="light yellow", maxpixels=5000000, zlim=c(0,min), axes=FALSE, legend=FALSE)
  #plot(rast, col=bgy[10], zlim=c(max,2), axes=FALSE, add=TRUE, legend=FALSE)
  plot(rast, col=bgy, zlim=c(min,max), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(p, col="gray", add=TRUE)
  plot(l, col="gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="dark gray", add=TRUE)
  #text(2200000,9400000,"Potential density (males/ha)", cex=1.3)
  dev.off()
}

#maps with range boundaries and occurrence points
brtplot2 <- function (rast,spec,range) {
  prev <- cellStats(rast, 'mean')	
  q99 <- quantile(rast, probs=c(0.99))
  max <- max(3*prev,q99)
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  write.csv(as.data.frame(occur), file=paste(x,spec,"_occur.csv",sep=""), row.names=FALSE)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(w,spec,"_pred1km2.png",sep=""), width=2250, height=1600, res=216)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,0,0))
  plot(rast, col="blue", maxpixels=5000000, axes=FALSE, legend=FALSE)
  plot(rast, col=bluegreen.colors(15), maxpixels=5000000, zlim=c(0,max), axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95), axis.args=list(cex.axis=1.3))
  plot(p, col="gray", add=TRUE)
  plot(l, col="gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark red", add=TRUE)
  points(occurcan[,7:8], col = 'black', pch=1, cex=0.4)
  text(2200000,9000000,"Potential density (males/ha)", cex=1.3)
  dev.off()
}

setwd(w)
for (i in 2:length(specpred)){
    models <- list.files(paste0(w,specpred[i],"/"),pattern="Mean.tif$")
    rast <- raster(paste0(w,specpred[i],"/",models[1]))
    #rast <- extend(rast,r)
    spec <- substr(models[1],6,9)
    #x1<-try(range <- shapefile(paste(natureserve,spec,".shp",sep="")))
    brtplot(rast,spec)
    #if (class(x1)!="try-error"){
    #range <- range[range$ORIGIN %in% list(2,1),]
    #try(brtplot2(rast,spec,range))
    #}
}
gc()



