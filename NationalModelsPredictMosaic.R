library(raster)
library(dismo)
library(rgdal)
library(gbm)
library(maptools)
library(dplyr)
library(sf)
library(Matrix)
library(reshape2)

bluegreen.colors <- colorRampPalette(c("#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.8)
p<- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("E:/GIS/hydrology/lakes_lcc.shp")
bcr <- rgdal::readOGR("E:/GIS/basemaps/BCRs/bcrfinallcc.shp")
canada <- rgdal::readOGR("E:/GIS/basemaps/canadaLCC.shp")
natureserve <- "E:/GIS/NatureServe/Abbreviated/_lcc/"
LCC <- CRS(projection(canada))
w <-"F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Nov2019/mosaics/"
x <-"F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Nov2019/"

load("F:/BAM/BAMData/BAM_data_package_November2019.RData")

#generate maps
brtplot <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(x,spec,"_pred1km1.png",sep=""), width=2500, height=1800, res=216)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(spec),"current prediction"))
  plot(rast, col=bluegreen.colors(15), maxpixels=5000000, zlim=c(0,max), axes=FALSE, main=as.character(spec), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.70,0.90,0.80,0.85), axis.args=list(cex.axis=1.3))
  plot(p, col="gray", add=TRUE)
  plot(l, col="gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="dark gray", add=TRUE)
  text(2200000,9400000,"Potential density (males/ha)", cex=1.3)
  dev.off()
}

#maps with range boundaries and occurrence points
brtplot2 <- function (rast,spec,range) {
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  PC1 <- PCmatch[PCmatch$SPECIES==spec,]
  occur <- left_join(PC1,SScombo,by="SS")
  occur <- occur[occur$ABUND > 0,]
  occur <- na.omit(occur)
  write.csv(as.data.frame(occur), file=paste(x,spec,"_occur.csv",sep=""), row.names=FALSE)
  occursp <- SpatialPointsDataFrame(coords = occur[,7:8], data = occur, proj4string = LCC)
  occurcan <- occursp[canada,]
  png(file=paste(x,spec,"_pred1km2.png",sep=""), width=2250, height=1800, res=216)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(spec),"current prediction"))
  plot(rast, col=bluegreen.colors(15), maxpixels=5000000, zlim=c(0,max), axes=FALSE, main=as.character(spec), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.70,0.90,0.80,0.85), axis.args=list(cex.axis=1.3))
  plot(p, col="gray", add=TRUE)
  plot(l, col="gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark red", add=TRUE)
  points(occurcan[,7:8], col = 'black', pch=1, cex=0.4)
  text(2200000,9000000,"Potential density (males/ha)", cex=1.3)
  dev.off()
}

setwd(w)
models <- list.files(pattern=".tif$")
r <- raster(models[2])
for (i in 1:length(models)){
    rast <- raster(models[i])
    rast <- extend(rast,r)
    spec <- substr(models[i],8,11)
    x1<-try(range <- shapefile(paste(natureserve,spec,".shp",sep="")))
    #brtplot(rast,spec)
    if (class(x1)!="try-error"){
    range <- range[range$ORIGIN %in% list(2,1),]
    brtplot2(rast,spec,range)
    }
}
gc()



