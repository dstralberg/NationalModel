library(raster)
#library(rgeos)
library(dismo)
#library(rgdal)
library(gbm)
#library(maptools)
library(dplyr)
library(sf)
library(Matrix)
library(reshape2)
library(colorspace)
library(terra)

#w <-"F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/artifacts/"
#x <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/feb2020/website/map-images/"
x <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/website/NewMosaicApproach/"

#bluegreen.colors <- colorRampPalette(c("#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.8)
#bgtrunc <- colorRampPalette(c("#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=10)

bgy <- sequential_hcl(10, "ag_GrnYl",rev=TRUE)
#bgy2 <- colorRamp(bgy, bias=0.8)
#blueyellow <- sequential_hcl(10, "BluYl",rev=TRUE)

adjust <- read.csv("G:/Shared drives/BAM_NationalModels4/NationalModels4.0/website/offset-adjustments-2025-04-04.csv")
#specpred <- unique(adjust$spp)
speclist <- read.csv("G:/Shared drives/BAM_NationalModels4/NationalModels4.0/website/speclist.csv")
speclist <- speclist[speclist$drop == 0,]
specpred <- as.vector(speclist[,1])

p<- shapefile("G:/Shared drives/GIS/basemaps/province_state_line.shp")
l <- shapefile("G:/Shared drives/GIS/hydrology/lakes_lcc.shp")
bcr <- shapefile("G:/Shared drives/GIS/basemaps/BCRs/bcrfinallcc.shp")
canada <- shapefile("G:/Shared drives/GIS/basemaps/canadaLCC.shp")
natureserve <- "G:/Shared drives/GIS/NatureServe/Abbreviated/"
LCC <- CRS(projection(canada))
#specpred <- list.dirs(w, full.names=FALSE)
#specpred <- list.files(x,pattern="_TSSRcorrected.tif$")
#specpred <- substr(specpred,1,4)

#models <- list.files(paste0(x,specpred[2],"/"),pattern="Mean.tif$")
#rast <- raster(paste0(w,specpred[2],"/",models[1]))

subr<- raster("G:/Shared drives/BAM_NationalModels4/NationalModels4.0/BCRUnits/BCRSubunits.tif")
models <- list.files(x,pattern=paste0("WeightedMosaic_",specpred[1]))
rast <- raster(paste0(x,models[1]))
rast <- mask(rast,subr)
bcrc <- raster::crop(bcr,rast)

subunits <- shapefile("G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/BCRSubunits/BCRSubunits.shp")
#subr <- rasterize(subunits,rast)

#load("D:/BAM/BAMData/BAMdb-GNMsubset-2020-01-08.RData")

load("G:/Shared drives/BAM_NationalModels4/NationalModels4.0/data/BAMData/BAM_data_package_November2019.RData")
LCC <- CRS(projection(canada))
occur <- left_join(PCmatch,SScombo,by="SS")
occur <- occur[occur$ABUND > 0,]
occur <- na.omit(occur)
occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
occurcan <- occursp[canada,]
detect <- unique(occurcan[,c(1,2)])
detect$count <- 1
countd <- stats::aggregate(x=detect$count,by=list(detect$SPECIES), FUN="sum")
#write.csv(countd, file=paste0(w,"detections.csv"),row.names=FALSE)

#write.csv(occurcan[,2:6],file="G:/Boreal/NationalModelsV2/abund_noxy.csv",row.names=FALSE)
#write.csv(occurcan[,1:6],file="G:/Boreal/NationalModelsV2/GNMabund_noxy.csv",row.names=FALSE)

pkeymodel <- data.table::unique(as.data.frame(occurcan[,c(1,4,5,7,8)]))
#write.csv(pkeymodel, file="G:/Boreal/NationalModelsV2/pkeymodel.csv",row.names=FALSE)

#writeRaster(subr,file = "G:/Boreal/NationalModelsV2/BCRSubunits.tif",overwrite=TRUE)

#plot number of sampling points at 10-km resolution
load("H:/Shared drives/BAM_NationalModels4/NationalModels4.0/data/BAMData/BAMdb-patched-xy.RData")	
subr20k <- aggregate(subr, fact=20)
subr10k <- aggregate(subr, fact=10)
subr5k <- aggregate(subr, fact=5)
samprast <- rasterize(cbind(SScombo$X,SScombo$Y), subr10k, fun='count')
samprast2 <- rasterize(cbind(SScombo$X,SScombo$Y), subr20k, fun='count')
can <- rasterize(canada,samprast2)
sampcan <- mask(samprast2,can)


brtplotdens <- function (rast) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.001)
  zmin <- min(zmin,0.01)
  zmax <- cellStats(rast,'max')
  q99 <- quantile(rast, probs=c(0.99))
  m <- c(q99,zmax,0)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  rast2 <- reclassify(rast,rclmat)
  png(file=paste0(x,"_pointdens20k.png"), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast2, col=bgy, range=c(0,q99), maxcell=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="black", lwd=1.5, add=TRUE)
  plot(canada,col=NA,border="black", lwd=1.5, add=TRUE)
  dev.off()
}
#brtplotdens(sampcan)

#plot total density
setwd(x)
models <- list.files(x,pattern=paste0("WeightedMosaic_",specpred[1]))
rast <- raster(paste0(x,models[1]))
totdens <- mask(rast,subr)
for (i in 2:length(specpred)){
  models <- list.files(x,pattern=paste0("WeightedMosaic_",specpred[i]))
  rast <- raster(paste0(x,models[1]))
  rast <- mask(rast,subr)
  totdens <- totdens + rast
}
#writeRaster(totdens,filename=paste0(x,"_TotDens.tif"))

#rast <- raster(paste0(x,"_TotDens.tif"))
brtplotabund <- function (rast) {
  prev <- cellStats(rast, 'mean')	
  zmin <- cellStats(rast,'min')
  zmax <- cellStats(rast,'max')
  q99 <- quantile(rast, probs=c(0.99))
  rast2 <- clamp(rast,4,q99)
  png(file=paste0(x,"_TotDens.png"), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast2, col=bgy, range=c(6,q99), maxcell=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  #plot(rast, col="#255668", range=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="black", lwd=1.5, add=TRUE)
  plot(canada,col=NA,border="black", lwd=1.5, add=TRUE)
  dev.off()
}
brtplotabund(totdens)


brtplot1 <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.001)
  zmin <- min(zmin,0.01)
  q99 <- quantile(rast, probs=c(0.999))
  #max <- max(3*prev,q99)
  zmax <- cellStats(rast, 'max')
  png(file=paste0(x,spec,"_pred1km1.png"), width=2600, height=1600, res=216)
  par(cex.main=1.8, mfcol=c(1,1), mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="dark gray", add=TRUE)
  plot(p, col="black", add=TRUE)
  dev.off()
}

brtplot2 <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.005)
  zmin <- min(zmin,0.05)
  q99 <- quantile(rast, probs=c(0.999))
  #max <- max(3*prev,q99)
  zmax <- cellStats(rast, 'max')
  png(file=paste0(x,spec,"_pred1km2.png"), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, add=TRUE, legend=FALSE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="dark gray", add=TRUE)
  plot(p, col="black", add=TRUE)
  dev.off()
}

#maps with range boundaries and occurrence points
brtplot3 <- function (rast,spec,range) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.005)
  zmin <- min(zmin,0.05)
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  #write.csv(as.data.frame(occur), file=paste(w,spec,"_occur.csv",sep=""), row.names=FALSE)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km3.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

#maps with range boundaries only (for manuscript)
brtplot3a <- function (rast,spec,range) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.005)
  zmin <- min(zmin,0.05)
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  png(file=paste(x,spec,"_pred1km3a.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  dev.off()
}

#maps with points only
brtplot3b <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.005)
  zmin <- min(zmin,0.05)
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km3b.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(p, col="black", add=TRUE)
  points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

brtplot4 <- function (rast,spec,range) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.001)
  zmin <- min(zmin,0.01)
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  #write.csv(as.data.frame(occur), file=paste(w,spec,"_occur.csv",sep=""), row.names=FALSE)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km4.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", maxpixels=5000000, zlim=c(0,zmin), axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

#maps with range boundaries only (for manuscript)
brtplot4a <- function (rast,spec,range) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.001)
  zmin <- min(zmin,0.01)
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  png(file=paste(x,spec,"_pred1km4a.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", maxpixels=5000000, zlim=c(0,zmin), axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  dev.off()
}

#maps with points only
brtplot4b <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  zmin <- max(prev,0.001)
  zmin <- min(zmin,0.01)
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km4b.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", maxpixels=5000000, zlim=c(0,zmin), axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(p, col="black", add=TRUE)
  points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

brtplot5 <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  zmin <- min(prev,0.001)
  q99 <- quantile(rast, probs=c(0.999))
  #max <- max(3*prev,q99)
  zmax <- cellStats(rast, 'max')
  png(file=paste0(x,spec,"_pred1km5.png"), width=2600, height=1600, res=216)
  par(cex.main=1.8, mfcol=c(1,1), mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", maxpixels=5000000, zlim=c(0,zmin), axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="dark gray", add=TRUE)
  plot(p, col="black", add=TRUE)
  dev.off()
}

brtplot6 <- function (rast,spec,range) {
  prev <- cellStats(rast, 'mean')	
  zmin <- min(prev,0.001)
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  #write.csv(as.data.frame(occur), file=paste(w,spec,"_occur.csv",sep=""), row.names=FALSE)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km6.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

#maps with range boundaries only (for manuscript)
brtplot6a <- function (rast,spec,range) {
  prev <- cellStats(rast, 'mean')	
  zmin <- min(prev,0.001)
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  png(file=paste(x,spec,"_pred1km6a.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  dev.off()
}

#maps with points only
brtplot6b <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  zmin <- min(prev,0.001)
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km6b.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(p, col="black", add=TRUE)
  points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

setwd(x)
for (i in 1:length(specpred)){
  models <- list.files(x,pattern=paste0("WeightedMosaic_",specpred[i]))
  rast <- raster(paste0(x,models[1]))
  rast <- mask(rast,subr)
  x1<-try(range <- shapefile(paste(natureserve,specpred[i],".shp",sep="")))
  brtplot1(rast,specpred[i])
  brtplot2(rast,specpred[i])
  brtplot5(rast,specpred[i])
  if (class(x1)!="try-error"){
    range1 <- spTransform(range,LCC)
    range1 <- range1[range1$ORIGIN %in% list(2,1),]
    try(brtplot3(rast,specpred[i],range1))
    try(brtplot3a(rast,specpred[i],range1))
    try(brtplot4(rast,specpred[i],range1))
    try(brtplot6(rast,specpred[i],range1))
    try(brtplot4a(rast,specpred[i],range1))
    try(brtplot6a(rast,specpred[i],range1))
  }
  gc()
}
gc()

setwd(x)
for (i in 1:length(specpred)){
  models <- list.files(x,pattern=paste0("WeightedMosaic_",specpred[i]))
  rast <- raster(paste0(x,models[1]))
  rast <- mask(rast,subr)
  try(brtplot3b(rast,specpred[i]))
  try(brtplot4b(rast,specpred[i]))
  try(brtplot6b(rast,specpred[i]))
}
gc()


for (i in 1:length(specpred)){
  models <- list.files(x,pattern=paste0(specpred[i],"_TSSRcorrected.tif$"))
  rast <- raster(paste0(x,models[1]))
  rast <- mask(rast,subr)
  writeRaster(rast,file = paste0(x,specpred[i],"_MeanCrop.tif"),overwrite=TRUE)
}


#CAWA map for manuscript
bcr60 <- subunits[subunits$BCRChar=="6-0",]
bcr12 <- subunits[subunits$BCRChar=="12",]
rast <- raster(paste0(x,"WeightedMosaic_CAWA.tiff"))
rast <- mask(rast,subr)
x1<-try(range <- shapefile(paste0(natureserve,"CAWA.shp",sep="")))
range <- range[range$ORIGIN %in% list(2,1),]
range1 <- spTransform(range,LCC)
prev <- cellStats(rast, 'mean')	
zmin <- max(prev,0.005)
zmin <- min(zmin,0.05)
zmax <- cellStats(rast, 'max')
q99 <- quantile(rast, probs=c(0.999))
png(file=paste(x,"CAWA_pred1km_ms.png",sep=""), width=2600, height=1600, res=216)
par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
plot(bcrc, col=NA, border=NA, axes=FALSE)
plot(bcr, col="white", border=NA, add=TRUE)
plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
plot(l, col="light gray", border=NA,add=TRUE)
plot(range1, col=NA, border="#3366cc", lwd=1.5, add=TRUE)
plot(p, col="black", add=TRUE)
plot(bcr60, col=NA, border="red", lwd=2, add=TRUE)
plot(bcr12, col=NA, border="red", lwd=2, add=TRUE)
dev.off()

#CONW map for manuscript
bcr60 <- subunits[subunits$BCRChar=="6-0",]
bcr81 <- subunits[subunits$BCRChar=="8-1",]
rast <- raster(paste0(x,"WeightedMosaic_CONW.tiff"))
rast <- mask(rast,subr)
x1<-try(range <- shapefile(paste0(natureserve,"CONW.shp",sep="")))
range <- range[range$ORIGIN %in% list(2,1),]
range1 <- spTransform(range,LCC)
prev <- cellStats(rast, 'mean')	
zmin <- max(prev,0.005)
zmin <- min(zmin,0.05)
zmax <- cellStats(rast, 'max')
q99 <- quantile(rast, probs=c(0.999))
png(file=paste(x,"CONW_pred1km_ms.png",sep=""), width=2600, height=1600, res=216)
par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
plot(bcrc, col=NA, border=NA, axes=FALSE)
plot(bcr, col="white", border=NA, add=TRUE)
plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
plot(l, col="light gray", border=NA,add=TRUE)
plot(range1, col=NA, border="#3366cc", lwd=1.5, add=TRUE)
plot(p, col="black", add=TRUE)
plot(bcr60, col=NA, border="red", lwd=2, add=TRUE)
plot(bcr81, col=NA, border="red", lwd=2, add=TRUE)
dev.off()

