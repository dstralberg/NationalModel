library(raster)
library(dismo)
library(rgdal)
library(gbm)
library(maptools)
library(dplyr)
library(sf)
library(Matrix)
library(reshape2)
library(mefa)
library(colorspace)

bgy <- sequential_hcl(11, "ag_GrnYl",rev=TRUE)
bgy2 <- sequential_hcl(4, "ag_GrnYl",rev=TRUE)
blueyellow <- sequential_hcl(10, "BluYl",rev=TRUE)

load("E:/BAM/BAMData/BAM_data_package_November2019.RData")
w <-"G:/Boreal/NationalModelsV2/"

p<- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("E:/GIS/hydrology/lakes_lcc.shp")
bcr <- rgdal::readOGR("E:/GIS/basemaps/BCRs/bcrfinallcc.shp")
canada <- rgdal::readOGR("E:/GIS/basemaps/canadaLCC.shp")
natureserve <- "E:/GIS/NatureServe/Abbreviated/_lcc/"
LCC <- CRS(projection(canada))

subunits <- shapefile("G:/Boreal/NationalModelsV2/BCRSubunits.shp")
subr <- rasterize(subunits,rast)

samplesp <- SpatialPointsDataFrame(coords = SScombo[,2:3], data = SScombo, proj4string = LCC)
samplecan <- samplesp[canada,]

#' Evaluate predictor sets based on hist, SD, etc
get_cn <- function(z, rmax=0.9) {
  SD <- apply(z, 2, sd)
  COR <- cor(z[,SD > 0])
  cr <- mefa:::stack.dist(as.dist(COR), dim.names = TRUE)
  cr <- cr[order(abs(cr$dist), decreasing = TRUE),]
  cr[,1] <- as.character(cr[,1])
  cr[,2] <- as.character(cr[,2])
  cr$Landsc1 <- startsWith(cr[,1], "Landsc750_")
  cr$Landsc2 <- startsWith(cr[,2], "Landsc750_")
  cr1 <- cr[cr$Landsc1 | cr$Landsc2,]
  cr2 <- cr[!(cr$Landsc1 | cr$Landsc2),]
  while(any(abs(cr1$dist) > rmax)) {
    i <- if (cr1$Landsc1[1])
      cr1[1,1] else cr1[1,2]
    j <- cr1[,1] == i | cr1[,2] == i
    cr1 <- cr1[!j,]
  }
  cr3 <- rbind(cr1, cr2)
  cr3 <- cr3[order(abs(cr3$dist), decreasing = TRUE),]
  while(any(abs(cr3$dist) > rmax)) {
    i <- if (cr3$Landsc1[1])
      cr3[1,1] else cr3[1,2]
    j <- cr3[,1] == i | cr3[,2] == i
    cr3 <- cr3[!j,]
  }
  vars <- union(as.character(cr3[,1]), as.character(cr3[,2]))
  return(vars)
}

brtplot <- function (rast,type) {
  # varimp <- as.data.frame(model$contributions)
  # write.csv(varimp,file=paste(w,type,"varimp.csv",sep=""))
  # cvstats <- as.data.frame(model$cv.statistics[c(1,3)])
  # cvstats$deviance.null <- model$self.statistics$mean.null
  # cvstats$deviance.exp <- (cvstats$deviance.null-cvstats$deviance.mean)/cvstats$deviance.null
  # write.csv(cvstats,file=paste(w,type,"cvstats.csv",sep=""))
  # pdf(paste(w,type,"_plot.pdf",sep=""))
  # gbm.plot(model,n.plots=12,smooth=TRUE)
  # dev.off()
  
  q <- c(10,20,30,40,50,60,70,80,90,100)
  pp <- quantile(rast, c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1))
  pp0 <- c(0,pp[1:9])
  m <- cbind(pp0,pp,q)
  zmin <- pp[1]
  # ix <- findInterval(getValues(rast), pp)
  # quant <- setValues(rast, ix)
  quant <- reclassify(rast,m)
  quant <- mask(quant,subr)
  bcrc <- crop(bcr,quant)
  
  png(file=paste0(w,type,"_predquant.png"), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(quant, col=bgy, maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(quant, col="#F9FFAF", maxpixels=5000000, zlim=c(0,zmin), axes=FALSE, add=TRUE, legend=FALSE)
  
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="dark gray", add=TRUE)
  plot(p, col="#8B0000", add=TRUE)
  dev.off()
  
}

#plot with sampling locations
brtplot1 <- function (rast,type) { 
  q <- c(10,20,30,40,50,60,70,80,90,100)
  pp <- quantile(rast, c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1))
  pp0 <- c(0,pp[1:9])
  m <- cbind(pp0,pp,q)
  zmin <- pp[1]
  # ix <- findInterval(getValues(rast), pp)
  # quant <- setValues(rast, ix)
  quant <- reclassify(rast,m)
  quant <- mask(quant,subr)
  bcrc <- crop(bcr,quant)
  
  png(file=paste0(w,type,"_predquant_points.png"), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(quant, col=bgy, maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  plot(quant, col="#F9FFAF", maxpixels=5000000, zlim=c(0,zmin), axes=FALSE, add=TRUE, legend=FALSE)
  
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="dark gray", add=TRUE)
  plot(p, col="black", add=TRUE)
  points(samplesp[,2:3], col = "#8B0000", pch=20, cex=0.4)
  
  dev.off()
  
}

#plot mess with sampling locations
brtplot2 <- function (rast) { 
  class <- c(1,2,3,4)
  lower <- c(-3000,-4,-2,0)
  upper <- c(-4,-2,0,2)
  # pp0 <- c(0,pp[1:9])
  m <- cbind(lower,upper,class)
  #zmin <- pp[1]
  # # ix <- findInterval(getValues(rast), pp)
  # # quant <- setValues(rast, ix)
  rastr <- reclassify(rast,m)
  rastr <- mask(rastr,subr)
  bcrc <- crop(bcr,rastr)
  
  png(file=paste0(w,"_mess_points.png"), width=2600, height=1600, res=216)
  par(cex.main=1.8, mar=c(0,0,0,0), bg="light gray", bty="n")
  plot(bcrc, col=NA, border=NA, axes=FALSE)
  plot(bcr, col="white", border=NA, add=TRUE)
  plot(rastr, col=bgy2, maxpixels=5000000, axes=FALSE, add=TRUE, horizontal = TRUE, smallplot = c(0.70,0.90,0.90,0.95))
  #plot(rastr, col="#F9FFAF", maxpixels=5000000, zlim=c(0,zmin), axes=FALSE, add=TRUE, legend=FALSE)
  
  plot(l, col="light gray", border=NA,add=TRUE)
  plot(bcr, col=NA, border="dark gray", add=TRUE)
  plot(p, col="black", add=TRUE)
  points(samplesp[,2:3], col = "#8B0000", pch=20, cex=0.3)
  
  dev.off()
  
}

models <- list.files(paste0("F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/artifacts/ALFL/"),pattern="Mean.tif$")
rast <- raster(paste0("F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/Feb2020/artifacts/ALFL/",models[1]))

# topog <-brick(paste(w,"topography_1km.grd",sep=""))
# names(topog) <- c("TPI","TRI","slope","roughness","lf")
# 
# landcover <- brick(paste(w,"landcov_1km.grd",sep=""))
# 
# road1k <- raster("E:/GIS/disturbance/VenterEtAlFootprint/RoadsLCC.tif")
# ROAD <- crop(road1k,rast)
# ROAD <- resample(ROAD,rast,method="ngb")
# projection(ROAD) <- LCC
# 
# biomass <- brick(paste(w,"bs2011_1km.grd",sep=""))
# landscape <- brick(paste(w,"bs2011_750_1km.grd",sep=""))
# names(landscape) <- gsub("LandCover","Landsc750",names(landscape))
# 
# cur <- "E:/Climate/CMIP5/baseline19812010/"
# setwd(cur)
# clim <- list.files(cur, pattern =".asc$")
# curclim<-stack(clim)
# 
# setwd(w)
# 
# SS <- cbind(SScombo, extract(landcover, SScombo[,2:3]))
# SS <- cbind(SS, extract(road1k, SScombo[,2:3]))
# names(SS)[ncol(SS)] <- "ROAD"
# SS <- cbind(SS, extract(topog, SScombo[,2:3]))
# SS <- cbind(SS, extract(biomass, SScombo[,2:3]))
# SS <- cbind(SS, extract(landscape, SScombo[,2:3]))
# 
# 
# samprast <- rasterize(cbind(SScombo$X,SScombo$Y), rast, field=1, fun='sum', background=0)
# samprast <- mask(samprast,rast)
# bs <- stack(samprast)
# names(bs) <-"samp"
# 
# clim1 <- mask(crop(curclim,rast),rast)
# landcov1 <- mask(crop(landcover,rast),rast)
# topo1 <- mask(crop(topog,rast),rast)
# landscape1 <- mask(crop(landscape,rast),rast)
# biomass1 <- mask(crop(biomass,rast),rast)
# ROAD1 <- mask(ROAD,rast)
# 
# bs <- addLayer(bs, ROAD1)
# names(bs)[nlayers(bs)] <- "ROAD"
# for (j in 1:nlayers(topo1)) {
#   bs <- addLayer(bs, topo1[[j]])}
# for (j in 1:nlayers(landcov1)) {
#   bs <- addLayer(bs, landcov1[[j]])}
# for (j in 1:nlayers(biomass1)) {
#   bs <- addLayer(bs, biomass1[[j]])}
# for (j in 1:nlayers(landscape1)) {
#   bs <- addLayer(bs, landscape1[[j]])}
# for (j in 1:nlayers(clim1)) {
#   bs <- addLayer(bs, clim1[[j]])}
# 
# writeRaster(bs,file=paste0(w,"rasters2011_all.grd"),overwrite=TRUE)
bs <- stack(paste0(w,"rasters2011_all.grd"))

# samp1kk <- sampleRandom(bs, size=100000, sp=TRUE)
# samp1kk$pres <- ifelse(samp1kk$samp > 0, 1, 0)
# samp1kkdf <- as.data.frame(samp1kk)
# write.csv(samp1kkdf,file=paste0(w,"surveysample1kk.csv"),row.names=FALSE)
samp1kkdf <- read.csv(paste0(w,"surveysample1kk.csv"))

potvar <- samp1kkdf[,2:222]
var <- get_cn(potvar)

x1 <- try(brt1 <- gbm.step(samp1kkdf, gbm.y = 1, gbm.x = var, family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5))
x2 <- try(brt2 <- gbm.step(samp1kkdf, gbm.y = 223, gbm.x = var, family = "bernoulli", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5))
save(brt1, file=paste0(w,"survey_poisson.RData"))
save(brt2, file=paste0(w,"survey_binomial.RData"))

predpois <- predict(bs, brt1, n.trees=brt1$n.trees, type="response")
writeRaster(predpois, filename=paste0(w,"survey_poisson"), format="GTiff",overwrite=TRUE)

predprob <- predict(bs, brt2, n.trees=brt2$n.trees, type="response")
writeRaster(predprob, filename=paste0(w,"survey_binomial"), format="GTiff",overwrite=TRUE)

predpois <- raster(paste0(w,"survey_poisson.tif"))
load(paste0(w,"survey_poisson.RData"))
brtplot(predpois, "poisson")

predprob <- raster(paste0(w,"survey_binomial.tif"))
load(paste0(w,"survey_binomial.RData"))
brtplot(predprob, "binomial")

#without roads
potvar1 <- samp1kkdf[,3:222]
var1 <- get_cn(potvar1)
x3 <- try(brt3 <- gbm.step(samp1kkdf, gbm.y = 1, gbm.x = var1, family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5))
x4 <- try(brt4 <- gbm.step(samp1kkdf, gbm.y = 223, gbm.x = var1, family = "bernoulli", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5))
save(brt3, file=paste0(w,"survey_poisson_noroad.RData"))
save(brt4, file=paste0(w,"survey_binomial_noroad.RData"))

# predpois3 <- predict(bs, brt3, n.trees=brt3$n.trees, type="response")
# writeRaster(predpois3, filename=paste0(w,"survey_poisson_noroad"), format="GTiff",overwrite=TRUE)
# 
# predprob4 <- predict(bs, brt4, n.trees=brt4$n.trees, type="response")
# writeRaster(predprob4, filename=paste0(w,"survey_binomial_noroad"), format="GTiff",overwrite=TRUE)

predpois3 <- raster(paste0(w,"survey_poisson_noroad.tif"))
#load(paste0(w,"survey_poisson_noroad.RData"))
brtplot1(predpois3, "poisson_noroad")

predprob4 <- raster(paste0(w,"survey_binomial_noroad.tif"))
#load(paste0(w,"survey_binomial_noroad.RData"))
brtplot1(predprob4, "binomial_noroad")

#MESS analysis
sampled <- samp1kkdf[samp1kkdf$samp>0,]
samprast <- rasterize(sampled[,223:224],bs[[1]],field=sampled[,1])
sampstack <- mask(bs,samprast)
sampsubset <- raster::subset(sampstack,var1,drop=FALSE)

sampsub <- sampled[,var1]

bssubset <- raster::subset(bs,var1,drop=FALSE)

# m1 <- mess(bssubset,sampsub,full=FALSE)
# writeRaster(m1,filename=paste0(w,"survey_mess"),format="GTiff",overwrite=TRUE)
m1 <- raster(paste0(w,"survey_mess.tif"))
brtplot2(m1)

baseline <- mess(sampsubset,sampsub,full=FALSE)
writeRaster(baseline,filename=paste0(w,"baseline_mess"),format="GTiff",overwrite=TRUE)

            
            