library(raster)
library(dismo)
library(gbm)
library(maptools)
library(dplyr)
library(sf)

bluegreen.colors <- colorRampPalette(c("#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.8)
bluegreen.colors2 <- colorRampPalette(c("#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.5)

provstate <- rgdal::readOGR("F:/GIS/basemaps/province_state_line.shp")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
w <-"G:/Boreal/NationalModelsV2/BCR6/June2019Results/"
bcr6 <- shapefile("G:/Boreal/NationalModelsV2/BCR6/bcr6.shp")
p<- rgdal::readOGR("F:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("F:/GIS/hydrology/lakes_lcc.shp")
lc <- crop(l,bcr6)

speclist <- read.csv("F:/BAM/BAMDAta/SpeciesClassesModv5.csv")
speclist <- speclist[speclist$NWT==1|speclist$Alberta==1,]
speclist <- speclist[speclist$spp!="CORE",]
speclist <- speclist[,1]

PC2011 <- read.csv("G:/Boreal/NationalModelsV2/BCR6/BCR6PC2011_v2.csv")
PC2001 <- read.csv("G:/Boreal/NationalModelsV2/BCR6/BCR6PC2001_v2.csv")
PC <- read.csv("G:/Boreal/NationalModelsV2/BCR6/BCR6PC_v2.csv")

mapplot <- function (j) {
  rast <- raster(paste(w,"CSpredicted",speclist[j],"year2011.tif",sep=""))
  futrast <- raster(paste(w,"CSpredicted",speclist[j],"year2100.tif",sep=""))
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(w,speclist[j],"_pred2011.png",sep=""), height=800, width=650)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()

  png(file=paste(w,speclist[j],"_pred2100.png",sep=""), width=650, height=800)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(futrast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(futrast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()

}

mapplot2 <- function (j) {
  rast <- raster(paste(w,"CSpredicted",speclist[j],"year2011.tif",sep=""))
  futrast <- raster(paste(w,"CSpredicted",speclist[j],"year2100.tif",sep=""))
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(w,speclist[j],"_pred2011_2.png",sep=""), height=800, width=650)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors2(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()
  
  png(file=paste(w,speclist[j],"_pred2100_2.png",sep=""), width=650, height=800)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(futrast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"2100 projection"))
  plot(futrast, col=bluegreen.colors2(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()
  
}

for (j in 1:length(speclist)) {
    mapplot(j)
}

for (j in 1:length(speclist)) {
  x1 <- try(raster(paste(w,"CSpredicted",speclist[j],"year2011.tif",sep="")))
  if (class(x1) != "try-error") {
  mapplot2(j)
  }
}
