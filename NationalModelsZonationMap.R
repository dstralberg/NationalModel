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

bgy <- sequential_hcl(5, "ag_GrnYl",rev=FALSE)
#bgy2 <- colorRamp(bgy, bias=0.8)
blueyellow <- sequential_hcl(5, "BluYl",rev=TRUE)


p<- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("E:/GIS/hydrology/lakes_lcc.shp")
bcr <- rgdal::readOGR("E:/GIS/basemaps/BCRs/bcrfinallcc.shp")
canada <- rgdal::readOGR("E:/GIS/basemaps/canadaLCC.shp")
natureserve <- "E:/GIS/NatureServe/Abbreviated/_lcc/"
LCC <- CRS(projection(canada))
w <-"F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/artifacts/"
x <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/feb2020/website/map-images/"

rast <- raster("I:/My Drive/BAM.SharedDrive/RshProjs/ConsPlanning/MultiSpeciesPlanning/Experimentation Work Project (Birds&Caribou)/Zonation/outputs_scenario1/scenario1.CAZ_E.rank.compressed.tif")
bcrc <- crop(bcr,rast)

subunits <- shapefile("G:/Boreal/NationalModelsV2/BCRSubunits.shp")
subr <- rasterize(subunits,rast)
rast <- mask(rast,subr)

png(file=paste0(x,"Zonation_63spp_Combo.png"), width=2600, height=1600, res=216)
par(cex.main=1.8, mfcol=c(1,1), mar=c(0,0,0,0), bg="light gray", bty="n")
plot(bcrc, col=NA, border=NA, axes=FALSE)
plot(bcr, col="white", border=NA, add=TRUE)
plot(rast, col=blueyellow,  maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
plot(l, col="light gray", border=NA,add=TRUE)
plot(bcr, col=NA, border="dark gray", add=TRUE)
plot(p, col="black", add=TRUE)
dev.off()




