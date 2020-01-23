library(sf)
library(rgeos)
library(raster)

w <-"G:/Boreal/NationalModelsV2/"
setwd(w)
bcr <- st_read("E:/GIS/basemaps/BCRs/bcrfinallcc.shp")
bcrus <- shapefile("E:/GIS/basemaps/BCRs/HemiborealUS_LCC.shp")

for (i in 4:14) {
  b <- bcr[bcr$BCR==i,]
  #br <- rasterize(b,r1k)
  bb <- st_buffer(b, dist=100000)
  bbb <- st_union(bb)
  st_write(bbb,paste("bcr",i,"_100km.shp",sep=""), delete_layer=TRUE)
}

bcrs <- shapefile("E:/GIS/basemaps/BCRs/bcrfinallcc.shp")
n60 <- shapefile("E:/GIS/basemaps/60N_LCC.shp")
prov <- shapefile("E:/GIS/basemaps/province_state_lcc.shp")
canada <- shapefile("E:/GIS/basemaps/canadaLCC.shp")
alaska <- shapefile("E:/GIS/basemaps/AlaskaLCC.shp")
bcra <- crop(bcrs,alaska)
bcrc <- crop(bcrs,canada)
# bcr60 <- gIntersection(bcrc,n60)
# bcr60b <- gBuffer(bcr60, width = 1)
# bcr1 <- gDifference(bcrc,bcr60b)
# bcr2 <- intersect(bcrc,bcr1)
# bcr3 <- intersect(bcr2,prov)
# #shapefile(bcr3,filename="G:/Boreal/NationalModelsV2/BCRunits.shp",overwrite=TRUE)

bcr3 <- shapefile("G:/Boreal/NationalModelsV2/BCRunits.shp")
bcrlist <- c(61,60,71,70,81,82,80,83)
for (i in bcrlist) {
  b <- bcr3[bcr3$BCR==i,]
  #br <- rasterize(b,r1k)
  bb <- gBuffer(b, width=100000)
  bbb <- union(bb)
  shapefile(bbb,filename=paste("bcr",i,"_100km.shp",sep=""), overwrite=TRUE)
}

bcrlistus <- c(10,11,12,13,14)
for (i in bcrlistus) {
  b <- bcrus[bcrus$BCR==i,]
  #br <- rasterize(b,r1k)
  bb <- gBuffer(b, width=100000)
  bbb <- union(bb)
  shapefile(bbb,filename=paste("bcrus",i,"_100km.shp",sep=""), overwrite=TRUE)
}

bcrlistak <- c(2,4,5)
for (i in bcrlistak) {
  b <- bcra[bcra$BCR==i,]
  #br <- rasterize(b,r1k)
  bb <- gBuffer(b, width=100000)
  bbb <- union(bb)
  shapefile(bbb,filename=paste("bcrus",i,"_100km.shp",sep=""), overwrite=TRUE)
}

bcrcombo <- bind(bcr3,bcra,bcrus)
shapefile(bcrcombo,filename="G:/Boreal/NationalModelsV2/BCRunits2.shp",overwrite=TRUE)