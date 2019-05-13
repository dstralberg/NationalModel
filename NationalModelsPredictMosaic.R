library(raster)
library(dismo)
library(gbm)
library(maptools)
library(dplyr)
library(sf)
library(Matrix)
library(reshape2)

bluegreen.colors <- colorRampPalette(c("#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.8)
p<- rgdal::readOGR("F:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("F:/GIS/hydrology/lakes_lcc.shp")
bcr <- rgdal::readOGR("F:/GIS/basemaps/BCRs/bcrfinallcc.shp")
natureserve <- "F:/GIS/NatureServe/Abbreviated/_lcc/"
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
w <-"G:/Boreal/NationalModelsV2/March2019/rasters/"
x <-"G:/Boreal/NationalModelsV2/March2019/mosaicpred/"
load("F:/BAM/BAMData/data_package_2016-04-18.Rdata")
PKEYBAM <- PKEY
PCBAM <- PCTBL
load("F:/BAM/BAMData/atlas_data_processed-20181018.RData")
PKEYAtlas <- PKEY
names(PKEYAtlas)[4] <- "YEAR"
PCAtlas <- PCTBL
load("F:/BAM/BamData/ARU/nwt-wildtrax-offsets-2019-01-16.RData")
PKEYWT <- unique(dd[,c(33,34,36)])
PCWT <- melt(y)
PCWT <-PCWT[PCWT$ABUND>0,]
names(PCWT) <- c("PKEY","SPECIES","ABUND")
load("F:/BAM/BamData/ARU/nwt-BU-offsets-2019-01-14.RData")
PKEYBU <- unique(dd[,c(14,15,17)])
PCBU <- melt(y)
PCBU <- PCBU[PCBU$ABUND>0,]
names(PCBU) <- c("PKEY","SPECIES","ABUND")
PCcombo <- rbind(PCBAM[,c(3:5)], PCAtlas[,c(2,4:5)], PCWT, PCBU)
PKEYcombo <- rbind(PKEYBAM[,c(1,2,8)],PKEYAtlas[,c(1,2,4)],PKEYWT,PKEYBU)
PC <- left_join(PCcombo,PKEYcombo)
load("F:/BAM/BAMData/BAMdb-patched-xy.RData")	
#PCSS <- left_join(PC,ss[,c(1,5)])

#generate maps
brtplot <- function (rast,spec) {
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(x,spec,"_pred1km1.png",sep=""), width=2250, height=1800, res=216)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(spec),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(spec), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.69,0.89,0.82,0.87), axis.args=list(cex.axis=1.3))
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
  PC1 <- PC[PC$SPECIES==spec,]
  occur <- inner_join(ss,PC1,by="SS")
  #write.csv(as.data.frame(occur), file=paste(x,spec,"_occur.csv",sep=""), row.names=FALSE)
  png(file=paste(x,spec,"_pred1km2.png",sep=""), width=2250, height=1800, res=216)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(spec),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(spec), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.69,0.89,0.82,0.87), axis.args=list(cex.axis=1.3))
  plot(p, col="gray", add=TRUE)
  plot(l, col="gray", border=NA,add=TRUE)
  plot(range, col=NA, border="dark red", add=TRUE)
  plot(st_geometry(occur), col = 'black', pch=1, cex=0.4, add = TRUE)
  text(2200000,9400000,"Potential density (males/ha)", cex=1.3)
  dev.off()
}

setwd(w)
models <- list.files(pattern=".tif")
for (i in 1:length(models)){
    rast <- raster(models[i])
    spec <- substr(models[i],8,11)
    x1<-try(range <- shapefile(paste(natureserve,spec,".shp",sep="")))
    brtplot(rast,spec)
    if (class(x1)!="try-error"){
    range <- range[range$ORIGIN %in% list(2,1),]
    brtplot2(rast,spec,range)
    }  
    gc()
}



