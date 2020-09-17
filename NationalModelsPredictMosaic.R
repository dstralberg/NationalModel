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

#summarize mean densities
# prev <- data.frame("species"=specpred[2:length(specpred)],"meandens"=0)
# for (i in 2:length(specpred)){
#   models <- list.files(paste0(w,specpred[i],"/"),pattern="Mean.tif$")
#   rast <- raster(paste0(w,specpred[i],"/",models[1]))
#   spec <- substr(models[1],6,9)
#   meandens <- cellStats(rast, 'mean')
#   prev$meandens[i-1] = meandens
# }
# write.csv(prev,file=paste0(w,"meandensities.csv"),row.names=FALSE)

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

setwd(w)
for (i in 2:151){
    models <- list.files(paste0(w,specpred[i],"/"),pattern="Mean.tif$")
    rast <- raster(paste0(w,specpred[i],"/",models[1]))
    spec <- substr(models[1],6,9)
    #x1<-try(range <- shapefile(paste(natureserve,spec,".shp",sep="")))
    brtplot1(rast,spec)
    brtplot2(rast,spec)
    brtplot5(rast,spec)
    # if (class(x1)!="try-error"){
    #   range <- range[range$ORIGIN %in% list(2,1),]
    #   try(brtplot3(rast,spec,range))
    # }
    # if (class(x1)!="try-error"){
    #   range <- range[range$ORIGIN %in% list(2,1),]
    #   try(brtplot4(rast,spec,range))
    # }
    # if (class(x1)!="try-error"){
    #   range <- range[range$ORIGIN %in% list(2,1),]
    #   try(brtplot6(rast,spec,range))
    # }
}
gc()

# for (j in 2:length(specpred)) {
#   species<-specpred[j]
#   r <- raster(paste0(w,species,"/pred-",species,"-CAN-Mean.tif"))
#   writeRaster(r,filename=paste0("F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/feb2020/website/pred-",species,"-CAN-Mean.tif"), format="GTiff", overwrite=TRUE)
# }  
    
