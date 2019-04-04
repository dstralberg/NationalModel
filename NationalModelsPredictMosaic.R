library(raster)
library(dismo)
library(gbm)
library(maptools)
library(dplyr)
library(sf)

bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.5)
p<- rgdal::readOGR("F:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("F:/GIS/hydrology/lakes_lcc.shp")
bcr <- rgdal::readOGR("F:/GIS/basemaps/BCRs/bcrfinallcc.shp")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
w <-"G:/Boreal/NationalModelsV2/March2019/rasters/"
x <-"G:/Boreal/NationalModelsV2/March2019/mosaicpred/"


#generate maps
brtplot <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(x,spec,"_pred1km1.png",sep=""), height=600, width=750)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(spec),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(spec), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.71,0.89,0.82,0.87), axis.args=list(cex.axis=1.3))
  plot(p, col="gray", add=TRUE)
  plot(l, col="gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="dark gray", add=TRUE)
  text(2200000,9400000,"Potential density (males/ha)", cex=1.3)
  dev.off()
}

setwd(w)
models <- list.files(pattern=".tif")
for (i in 1:length(models)){
    rast <- raster(models[i])
    spec <- substr(models[i],8,11)
    brtplot(rast,spec)
    gc()
}



