library(parallel)
library(gbm)
library(raster)
library(rgdal)
library(colorspace)

PROJ <- "boot"
ROOT1 <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/May2020/"
ROOT2 <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/"
SPP <- list.files(file.path(ROOT2, "out", PROJ))
SPP2 <- c("BBWA","BHVI","BLJA","BLPW","BOCH","CMWA","CONW","GCKI", "GRAJ","LEFL","OVEN","OSFL","PAWA","PHVI","RUBL","SWTH","TEWA","YRWA")

B <- 32
u <- c(60, 61, 70, 71, 80, 81, 82, 83, 9, 10, 11, 12, 13, 14, 4, 5, 
       200, 400, 500, 1000, 1100, 1200, 1300, 1400)

bgy <- sequential_hcl(10, "ag_GrnYl",rev=TRUE)
blueyellow <- sequential_hcl(10, "BluYl",rev=TRUE)
bluered <- diverging_hcl(10, c = c(100, 0), l = c(50, 90), power = 1.3,rev=TRUE)
tropic <- diverging_hcl(7, palette = "Tropic", h2 = 0,rev=TRUE)

p<- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("E:/GIS/hydrology/lakes_lcc.shp")
bcr <- rgdal::readOGR("E:/GIS/basemaps/BCRs/bcrfinallcc.shp")
canada <- rgdal::readOGR("E:/GIS/basemaps/canadaLCC.shp")
natureserve <- "E:/GIS/NatureServe/Abbreviated/_lcc/"
abbcr6 <- rgdal::readOGR("E:/GIS/basemaps/BCRs/AlbertaBCR6LCC.shp")
bcr6 <- bcr[bcr$BCR==6,]
l <- crop(l,bcr6)
LCC <- CRS(projection(canada))

brtplot2 <- function (rast,spp,prev) {
  #if (file.exists(paste0(ROOT1,spp,"_BCR",BCR,".png"))==FALSE) {
  zmin <- max(prev,0.005)
  zmin <- min(zmin,0.05)
  #q99 <- quantile(rast, probs=c(0.99))
  #max <- max(3*prev,q99)
  zmax <- cellStats(rast, 'max')
  png(file=paste0(ROOT1,spp,"_BCR",BCR,".png"), width=1800, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,3,0), bg="white", bty="n")
  #plot(bcrc, col=NA, border=NA, axes=FALSE)
  #plot(bcr, col="white", border=NA, add=TRUE)
  rast <- trim(rast)
  plot(rast, col="#F9FFAF", maxpixels=5000000, zlim=c(0,zmin), axes=FALSE, legend=FALSE, main=spp)
  plot(rast, col=bgy, zlim=c(zmin,zmax), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.58,0.78,0.85,0.9))
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr6, col=NA, border="black", add=TRUE)
  plot(p, col="black", add=TRUE)
  dev.off()
  #}
}


for (b in 1:B) {
  for (BCR in u) {
    for (spp in SPP) {
      r1 <- raster(file.path(ROOT2, "data", "templates", paste0("bcr-template-", BCR, ".grd")))
        fout <- file.path(ROOT1,paste0("pred250-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif"))
        a <- try(raster(fout))
        if(class(a) == "try-error"){
          YS <- stack(file.path(ROOT1, "stacks250m", paste0("bcr", BCR, "all_250m.grd")))
          ARU <- YS[[1]]*0 
          YEAR <- ARU + 2011
          YS <- addLayer(YS,ARU,YEAR)
          names(YS)[(nlayers(YS)-1):(nlayers(YS))] <- c("ARU","YEAR")
          fin <- file.path(ROOT2, "out", PROJ, spp, paste0("BCR_", BCR), paste0("gnmboot-", spp, "-BCR_", BCR, "-", b, ".RData"))
          e <- new.env()
          load(fin, envir=e)
          if(class(e$out) != "try-error") {
            pred <- raster::predict(YS,e$out, n.trees=e$out$n.trees, type="response")
            writeRaster(pred,fout,overwrite=TRUE)
          }
        }
      }
    }
  }
