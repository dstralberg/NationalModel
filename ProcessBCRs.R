library(sf)
w <-"G:/Boreal/NationalModelsV2/"
setwd(w)
bcr <- st_read("F:/GIS/basemaps/BCRs/bcrfinallcc.shp")

for (i in 4:14) {
  b <- bcr[bcr$BCR==i,]
  #br <- rasterize(b,r1k)
  bb <- st_buffer(b, dist=100000)
  bbb <- st_union(bb)
  st_write(bbb,paste("bcr",i,"_100km.shp",sep=""), delete_layer=TRUE)
}

library(rgeos)
library(raster)
bcrs <- shapefile("F:/GIS/basemaps/BCRs/bcrfinallcc.shp")
bcrc <- crop(bcrs,canada)
n60 <- shapefile("F:/GIS/basemaps/60N_LCC.shp")
prov <- shapefile("F:/GIS/basemaps/province_state_lcc.shp")
canada <- shapefile("F:/GIS/basemaps/canadaLCC.shp")
bcr60 <- gIntersection(bcrc,n60)
bcr60b <- gBuffer(bcr60, width = 1)
bcr1 <- gDifference(bcrc,bcr60b)
bcr2 <- intersect(bcrc,bcr1)
bcr3 <- intersect(bcr2,prov)
shapefile(bcr3,filename="G:/Boreal/NationalModelsV2/BCRunits.shp",overwrite=TRUE)

bcr <- shapefile("G:/Boreal/NationalModelsV2/BCRunits.shp")
bcrlist <- c(61,60,71,70,81,82,80,83)
for (i in bcrlist) {
  b <- bcr[bcr$BCR==i,]
  #br <- rasterize(b,r1k)
  bb <- gBuffer(b, width=100000)
  bbb <- union(bb)
  shapefile(bbb,filename=paste("bcr",i,"_100km.shp",sep=""), overwrite=TRUE)
}