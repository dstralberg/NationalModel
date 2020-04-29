library(parallel)
library(gbm)
library(raster)
library(rgdal)
library(colorspace)

PROJ <- "boot"
ROOT1 <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/"
ROOT2 <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/"
SPP <- list.files(file.path(ROOT1, "out", PROJ))
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

brtplot2 <- function (rast,spp,YR,prev) {
  #if (file.exists(paste0(ROOT2,"out/","parts/",spp,"/",spp,YR,"_BCR",BCR,".png"))==FALSE) {
  zmin <- max(prev,0.005)
  zmin <- min(zmin,0.05)
  #q99 <- quantile(rast, probs=c(0.99))
  #max <- max(3*prev,q99)
  zmax <- cellStats(rast, 'max')
  png(file=paste0(ROOT2,"out/","parts/",spp,"/",spp,YR,"_BCR",BCR,".png"), width=1800, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,3,0), bg="white", bty="n")
  #plot(bcrc, col=NA, border=NA, axes=FALSE)
  #plot(bcr, col="white", border=NA, add=TRUE)
  rast <- trim(rast)
  plot(rast, col="#F9FFAF", maxpixels=5000000, zlim=c(0,zmin), axes=FALSE, legend=FALSE, main=paste(spp,YR))
  plot(rast, col=bgy, zlim=c(zmin,zmax), maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.58,0.78,0.85,0.9))
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr6, col=NA, border="black", add=TRUE)
  plot(p, col="black", add=TRUE)
  dev.off()
  #}
}

chgplot <- function (rast,spp) {
  #if (file.exists(paste0(ROOT2,"out/","parts/",spp,"/",spp,"_BCR",BCR,"_1996-2016.png"))==FALSE) {
  png(file=paste0(ROOT2,"out/","parts/",spp,"/",spp,"_BCR",BCR,"_1996-2016.png"), width=1800, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,3,0), bg="white", bty="n")
  plot(rast, col="gray", maxpixels=5000000, axes=FALSE, legend=FALSE, main=paste(spp,"1996 - 2016 change"))
  plot(rast, col=bluered, maxpixels=5000000, axes=FALSE, add=TRUE, main=spp, horizontal = TRUE, smallplot = c(0.58,0.78,0.85,0.90))
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr6, col=NA, border="black", add=TRUE)
  plot(p, col="black", add=TRUE)
  dev.off()
  #}
}

year_stack <- function (YR) {
  r1 <- raster(file.path(ROOT1, "data", "templates", paste0("bcr-template-", BCR, ".grd")))
  ND <- stack(file.path(ROOT1, "data", "stacks2011", paste0("bcr", BCR, "all_1km.grd")))
  ARU <- ND[[1]]*0 
  YEAR <- ARU + YR
  ND <- addLayer(ND,ARU,YEAR)
  names(ND)[(nlayers(ND)-1):(nlayers(ND))] <- c("ARU","YEAR")
  return(ND)
  # ND$data$YEAR <- 2011
  # assign("r1", r1, envir=.GlobalEnv)
  # assign("ND", ND, envir=.GlobalEnv)
  # invisible(TRUE)
}

year_stack_2001 <- function (YR) {
  r1 <- raster(file.path(ROOT1, "data", "templates", paste0("bcr-template-", BCR, ".grd")))
  ND <- stack(file.path(ROOT1, "data", "stacks2001", paste0("bcr", BCR, "all_2001_1km.grd")))
  ARU <- ND[[1]]*0 
  YEAR <- ARU + YR
  ND <- addLayer(ND,ARU,YEAR)
  names(ND)[(nlayers(ND)-1):(nlayers(ND))] <- c("ARU","YEAR")
  return(ND)
  # ND$data$YEAR <- 2011
  # assign("r1", r1, envir=.GlobalEnv)
  # assign("ND", ND, envir=.GlobalEnv)
  # invisible(TRUE)
}


for (b in 1:B) {
  for (BCR in u) {
    for (spp in SPP) {
      r1 <- raster(file.path(ROOT1, "data", "templates", paste0("bcr-template-", BCR, ".grd")))
      for (YR in 2006:2019) {
        fout <- file.path(ROOT2, "out", "parts", spp,paste0("pred-", YR, "-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif"))
        a <- try(raster(fout))
        if(class(a) == "try-error"){
          YS <- year_stack(YR)
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
}

for (b in 1:B) {
  for (BCR in u) {
    for (spp in SPP) {
      r1 <- raster(file.path(ROOT1, "data", "templates", paste0("bcr-template-", BCR, ".grd")))
      for (YR in 1992:2005) {
        fout <- file.path(ROOT2, "out", "parts", spp,paste0("pred-", YR, "-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif"))
        a <- try(raster(fout))
        if(class(a) == "try-error"){
          YS <- year_stack_2001(YR)
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
}


#maps
for (b in 1:B) {
  for (BCR in u) {
    for (spp in SPP) {
      z1 <- try(pred96 <- raster(paste0(ROOT2, "out/", "parts/", spp,"/pred-", 1996, "-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif")))
      if (class(z1) != "try-error") {
      prev <- cellStats(pred96, 'mean')
      brtplot2(pred96,spp,1996,prev)
      pred16 <- raster(paste0(ROOT2, "out/", "parts/", spp,"/pred-", 2016, "-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif"))
      brtplot2(pred16,spp,2016,prev)
      predchg <- pred16 - pred96
      chgplot(predchg,spp)
      for (YR in 1997:2015) {
        pred <- raster(paste0(ROOT2, "out/", "parts/", spp,"/pred-", YR, "-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif"))
        brtplot2(pred,spp,YR,prev)
        }
      }
    }
  }  
}

#plots
for (b in 1:B) {
  for (BCR in u) {
    spectrend <- data.frame(SPP,trend=0)
    j<-1
    for (spp in SPP) {
      abund <-data.frame(YEAR=0,SPECIES=spp,abund=0)
      i<-1
      for (YR in 1996:2016) {
        x <- try(pred <- raster(paste0(ROOT2, "out/", "parts/", spp,"/pred-", YR, "-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif")))
        pred <- mask(pred,abbcr6)
        if (class(x) !="try-error") {
        abund[i,3] <- 200*cellStats(pred,stat='sum')
        abund[i,1] <- YR
        abund[i,2] <- spp
        }
        i<-i+1
      }
      if (class(x) != "try-error"){
      tm <- lm(abund$abund ~ abund$YEAR)
      trend <- (tm$coefficients[2]/max(abund$abund))
      png(file=paste0(ROOT2,"out/","parts/",spp,"/",spp,"_BCR",BCR,"_1996-2016_trend.png"), width=2600, height=1600, res=216)
      par(0,0,3,0)
      plot(abund$YEAR, abund$abund, type="l", xlab="YEAR",ylab="Abundance", main=paste0(spp,", trend = ", signif(trend,3)))
      abline(tm$coefficients[1], tm$coefficients[2], col="red")
      #text(1998, 1000, label=paste0("trend = ", signif(trend,3)))
      dev.off()
      spectrend[j,2] <-trend}
      j<-j+1
    }
    write.csv(spectrend, file=paste0(ROOT2,"out/","parts/_BCR",BCR,"_1996-2016_trend.csv"),row.names=FALSE )
  }
}

