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

#bluegreen.colors <- colorRampPalette(c("#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.8)
#bgtrunc <- colorRampPalette(c("#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=10)

bgy <- sequential_hcl(10, "ag_GrnYl",rev=TRUE)
#bgy2 <- colorRamp(bgy, bias=0.8)
#blueyellow <- sequential_hcl(10, "BluYl",rev=TRUE)

load("E:/BAM/BAMData/BAM_data_package_November2019.RData")

p<- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("E:/GIS/hydrology/lakes_lcc.shp")
bcr <- rgdal::readOGR("E:/GIS/basemaps/BCRs/bcrfinallcc.shp")
canada <- rgdal::readOGR("E:/GIS/basemaps/canadaLCC.shp")
natureserve <- "E:/GIS/NatureServe/Abbreviated/_lcc/"
LCC <- CRS(projection(canada))
w <-"F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/artifacts/"
x <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/feb2020/website/map-images/"
specpred <- list.dirs(w, full.names=FALSE)

models <- list.files(paste0(w,specpred[2],"/"),pattern="Mean.tif$")
rast <- raster(paste0(w,specpred[2],"/",models[1]))
bcrc <- crop(bcr,rast)

subunits <- shapefile("G:/Boreal/NationalModelsV2/BCRSubunits.shp")
subr <- rasterize(subunits,rast)


brtplotmean <- function (rast,spec,range) {
  prev <- cellStats(rast, 'mean')	
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  # PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  # occur <- left_join(PC1,SScombo,by="SS")
  # occur <- occur[occur$ABUND > 0,]
  # occur <- na.omit(occur)
  # occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  # occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km3.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,zmin), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(zmin,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  #plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  # points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

brtplotmedian <- function (rast,spec) {
  prev <- quantile(rast, probs=c(0.5))	
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  # PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  # occur <- left_join(PC1,SScombo,by="SS")
  # occur <- occur[occur$ABUND > 0,]
  # occur <- na.omit(occur)
  #write.csv(as.data.frame(occur), file=paste(w,spec,"_occur.csv",sep=""), row.names=FALSE)
  # occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  # occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km_median.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,prev), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(prev,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  #plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  # points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

brtplotmedpts <- function (rast,spec,range) {
  prev <- quantile(rast, probs=c(0.5))	
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km_median_pts.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,prev), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(prev,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

brtplot25pts <- function (rast,spec,range) {
  prev <- quantile(rast, probs=c(0.25))	
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km_25%ile_pts.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,prev), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(prev,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

brtplot75pts <- function (rast,spec,range) {
  prev <- quantile(rast, probs=c(0.75))	
  zmax <- cellStats(rast, 'max')
  q99 <- quantile(rast, probs=c(0.999))
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km_75%ile_pts.png",sep=""), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rast, col="#F9FFAF", zlim=c(0,prev), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(rast, col=bgy, zlim=c(prev,q99), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(rast, col="#255668", zlim=c(q99,zmax), maxpixels=5000000, axes=FALSE, legend=FALSE, add=TRUE)
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark blue", add=TRUE)
  plot(p, col="black", add=TRUE)
  points(occurcan[,7:8], col = "#0000007D", pch=20, cex=0.4)
  dev.off()
}

setwd(w)
for (i in 111:151){
    models <- list.files(paste0(w,specpred[i],"/"),pattern="Mean.tif$")
    spec <- substr(models[1],6,9)
    x1<-try(range <- shapefile(paste(natureserve,spec,".shp",sep="")))
    rast <- raster(paste0(w,specpred[i],"/",models[1]))
    rast <- mask(rast,subr)
    brtplotmedian(rast,spec)
    brtplotmedpts(rast,spec,range)
    brtplot25pts(rast,spec,range)
    brtplot75pts(rast,spec,range)
}
gc()

    
