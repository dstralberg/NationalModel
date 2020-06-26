library(parallel) #not used
library(gbm)
library(raster)
library(rgdal)

PROJ <- "boot"
ROOT1 <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/May2020/" #input models
ROOT2 <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/" #input BCR templates
ROOT3 <- "I:/My Drive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/May2020/" #output predictions
SPP <- list.files(file.path(ROOT2, "out", PROJ))
SPP2 <- c("BBWA","BHVI","BLJA","BLPW","BOCH","CMWA","CONW","GCKI", "GRAJ","LEFL","OVEN","OSFL","PAWA","PHVI","RUBL","SWTH","TEWA","YRWA")

B <- 32
u <- c(60, 61, 70, 71, 80, 81, 82, 83, 9, 10, 11, 12, 13, 14, 4, 5, 
       200, 400, 500, 1000, 1100, 1200, 1300, 1400)

for (BCR in u) {
  for (b in 1:B) {
    for (spp in SPP2) {
      r1 <- raster(file.path(ROOT2, "data", "templates", paste0("bcr-template-", BCR, ".grd")))
        fout <- file.path(ROOT3,paste0("pred250-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif"))
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
