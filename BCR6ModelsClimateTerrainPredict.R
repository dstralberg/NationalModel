library(raster)
library(dismo)
library(gbm)
library(maptools)
library(dplyr)
library(sf)

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

bluegreen.colors <- colorRampPalette(c("#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.8)
provstate <- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
w <-"G:/Boreal/NationalModelsV2/"
bcr6 <- shapefile("G:/Boreal/NationalModelsV2/BCR6/bcr6.shp")
p<- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("E:/GIS/hydrology/lakes_lcc.shp")
lc <- crop(l,bcr6)

speclist <- read.csv("F:/BAM/BAMDAta/SpeciesClassesModv5.csv")
speclist <- speclist[speclist$NWT==1|speclist$Alberta==1,]
speclist <- speclist[speclist$spp!="CORE",]
speclist <- speclist[,1]

#speclist <- as.factor(c(as.character(speclist),"CAWA","RUBL"))

# bs2001 <- stack(paste(w,"bcr6_2001rasters250.grd",sep=""))
# bs2011 <- stack(paste(w,"bcr6_2011rasters250.grd",sep=""))
bs<-stack("G:/Boreal/NationalModelsV2/bcr6clim_1km.grd")
bs2 <- stack("G:/Boreal/NationalModelsV2/BCR6/bcr6_2011rasters1km.grd")
bs <- crop(bs,bs2)
bs <- resample(bs,bs2)
bs3 <- stack(bs,bs2)

bsf <-stack("G:/Boreal/NationalModelsV2/BCR6/bcr6_clim2050_CANESM2_RCP45.grd")
bsf <- resample(bsf,bs2)
bs2050 <- stack(bsf,bs2)
bsf <-stack("G:/Boreal/NationalModelsV2/BCR6/bcr6_clim2080_CanESM2_RCP45.grd")
bsf <- resample(bsf,bs2)
bs2080 <- stack(bsf,bs2)

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# offl <- read.csv("G:/Boreal/NationalModelsV2/Quebec/BAMoffsets.csv")
# offla <- read.csv("G:/Boreal/NationalModelsV2/Quebec/Atlasoffsets.csv")
# offlc <- rbind(offl[2:4],offla[2:4])
# offlc$PKEY <- as.character(offlc$PKEY)
# offlc$SPECIES <- as.character(offlc$SPECIES)
# offlb <- read.csv("G:/Boreal/NationalModelsV2/BCR6/offwt.csv")
# offlb$PKEY <- as.character(offlb$PKEY)
# offlb$SPECIES <- as.character(offlb$SPECIES)
# offld <- read.csv("G:/Boreal/NationalModelsV2/BCR6/offbu.csv")
# offld$PKEY <- as.character(offld$PKEY)
# offld$SPECIES <- as.character(offld$SPECIES)
# offcombo <- rbind(offlc,offlb,offld)

cdat <- read.csv("G:/Boreal/NationalModelsV2/BCR6/bcr6_cdat_v3.csv") #n=48500
cdat$SS <- as.character(cdat$SS)
#cdat <- cdat[!duplicated(cdat[, 2:3]), ]
#cdat$count <- 1

PC2011 <- read.csv(paste(w,"BCR6/BCR6PC2011_v3.csv",sep=""))
PC2001 <- read.csv(paste(w,"BCR6/BCR6PC2001_v3.csv",sep=""))
PC <- read.csv("G:/Boreal/NationalModelsV2/BCR6/BCR6PC_v3.csv")

survey2001 <- aggregate(PC2001$ABUND, by=list("PKEY"=PC2001$PKEY,"SS"=PC2001$SS), FUN=sum) 
survey2001 <- survey2001[sample(1:nrow(survey2001)), ]
survey2011 <- aggregate(PC2011$ABUND, by=list("PKEY"=PC2011$PKEY,"SS"=PC2011$SS), FUN=sum) 
survey2011 <- survey2011[sample(1:nrow(survey2011)), ]

#calculating sample weights as inverse of survey effort within 5x5 pixel area
samprast <- rasterize(cbind(cdat$X,cdat$Y), r2, field=1, fun='sum')
sampsum25 <- focal(samprast, w=matrix(1,nrow=5,ncol=5), na.rm=TRUE)
dat_c <- cbind(cdat,extract(sampsum25,cbind(cdat$X,cdat$Y))) #n=48500
names(dat_c)[ncol(dat_c)] <- "sampsum25"
dat_c$wt <- 1/dat_c$sampsum25
dat_c$SS <- as.character(dat_c$SS)

w1 <- "G:/Boreal/NationalModelsV2/BCR6/"
setwd(w1)

#generate current predictions and plots from models
brtplot <- function (j) {
  load(paste(w1,speclist[j],"brt5.R",sep=""))
  varimp <- as.data.frame(brt1$contributions)
  write.csv(varimp,file=paste(w,speclist[j],"varimp5.csv",sep=""))
  cvstats <- t(as.data.frame(brt1$cv.statistics))
  write.csv(cvstats,file=paste(w1,speclist[j],"cvstats5.csv",sep=""))
  pdf(paste(w1,speclist[j],"_plot5.pdf",sep=""))
  gbm.plot(brt1,n.plots=12,smooth=TRUE)
  dev.off()
  rast <- raster::predict(bs3, brt1, type="response", n.trees=brt1$n.trees)
  writeRaster(rast, filename=paste(w1,speclist[j],"_pred1km5",sep=""), format="GTiff",overwrite=TRUE)
  
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(w1,speclist[j],"_pred1km5.png",sep=""), height=800, width=650)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()

  PC1 <- PC[PC$SPECIES==as.character(speclist[j]),]
  PC1 <- PC1[PC1$ABUND>0,]
  xy <- PC1[,c(7,8)]
  spdf <- SpatialPointsDataFrame(coords = xy, data = PC1, proj4string = LCC)
  png(file=paste(w,speclist[j],"_pred1km5_pts.png",sep=""), width=650, height=800)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  plot(spdf, col = 'red', pch=1, cex=PC1$ABUND, add = TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()
}

mapplot <- function (j) {
  rast <- raster(paste(w1,speclist[j],"_pred1km5.tif",sep=""))
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(w1,speclist[j],"_pred1km5.png",sep=""), height=800, width=650)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()
  
  PC1 <- PC[PC$SPECIES==as.character(speclist[j]),]
  PC1 <- PC1[PC1$ABUND>0,]
  xy <- PC1[,c(7,8)]
  spdf <- SpatialPointsDataFrame(coords = xy, data = PC1, proj4string = LCC)
  png(file=paste(w1,speclist[j],"_pred1km5_pts.png",sep=""), width=650, height=800)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  plot(spdf, col = 'red', pch=1, cex=PC1$ABUND, add = TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()
  
  PC1 <- PC[PC$SPECIES==as.character(speclist[j]),]
  PC1 <- PC1[PC1$ABUND>0,]
  xy <- PC1[,c(7,8)]
  spdf <- SpatialPointsDataFrame(coords = xy, data = PC1, proj4string = LCC)
  png(file=paste(w1,speclist[j],"_pred1km5_pts_small.png",sep=""), width=650, height=800)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  plot(spdf, col = 'red', pch=1, cex=0.4, add = TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()
}

#generate future predictions 
futplot <- function (j) {
  load(paste(w1,speclist[j],"brt5.R",sep=""))
  r2080 <- raster::predict(bs2080, brt1, type="response", n.trees=brt1$n.trees)
  writeRaster(r2080, filename=paste(w1,speclist[j],"_pred1km5_2080_CanESM2_RCP45",sep=""), format="GTiff",overwrite=TRUE)
  r2050 <- raster::predict(bs2050, brt1, type="response", n.trees=brt1$n.trees)
  writeRaster(r2050, filename=paste(w1,speclist[j],"_pred1km5_2050_CanESM2_RCP45",sep=""), format="GTiff",overwrite=TRUE)  
  
  rast <- raster(paste(w1,speclist[j],"_pred1km5.tif",sep=""))
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  
  png(file=paste(w1,speclist[j],"_pred1km5_2080_CanESM2_RCP45.png",sep=""), height=850, width=600)
  par(cex.main=1, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(r2080, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"2080 prediction"))
  plot(r2080, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.5))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(2400000,7950000,"Potential density (males/ha)", cex=1)
  dev.off()
  
  png(file=paste(w1,speclist[j],"_pred1km5_2050_CanESM2_RCP45.png",sep=""), height=850, width=600)
  par(cex.main=1, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(r2050, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"2050 prediction"))
  plot(r2050, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.5))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(2400000,7950000,"Potential density (males/ha)", cex=1)
  dev.off()  
  
}

futmapplot <- function (j) {
  r2080 <- raster(paste(w1,speclist[j],"_pred1km5_2080_CanESM2_RCP45.tif",sep=""))
  r2050 <- raster(paste(w1,speclist[j],"_pred1km5_2050_CanESM2_RCP45.tif",sep=""))  
  rast <- raster(paste(w1,speclist[j],"_pred1km5.tif",sep=""))
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  
  png(file=paste(w1,speclist[j],"_pred1km5_2080_CanESM2_RCP45.png",sep=""), height=850, width=600)
  par(cex.main=1, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(r2080, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"2080 prediction"))
  plot(r2080, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.5))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(2400000,7950000,"Potential density (males/ha)", cex=1)
  dev.off()
  
  png(file=paste(w1,speclist[j],"_pred1km5_2050_CanESM2_RCP45.png",sep=""), height=850, width=600)
  par(cex.main=1, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(r2050, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"2050 prediction"))
  plot(r2050, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.84,0.87), axis.args=list(cex.axis=1.5))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(2400000,7950000,"Potential density (males/ha)", cex=1)
  dev.off()  
  
}

cvstatsum <- function (speclist) {
  cvstats <- read.csv(paste(w1,speclist[1],"cvstats5.csv",sep=""))
  cvstatmean <- as.data.frame(cbind(as.character(cvstats[,1]),rowMeans(cvstats[,2:6])))
  names(cvstatmean) <- c("stat",as.character(speclist[1]))
  for (j in 2:length(speclist)) {
    x<-try(cv2 <- read.csv(paste(w1,speclist[j],"cvstats5.csv",sep="")))
    if(class(x) != "try-error") {
      cvstatmean <- as.data.frame(cbind(cvstatmean,rowMeans(cv2[,2:6])))
      names(cvstatmean)[ncol(cvstatmean)] <- as.character(speclist[j])
    }
  }
  return(cvstatmean)
}

for (j in 1:length(speclist)) {
  x<-try(rast <- raster(paste(w1,speclist[j],"_pred1km5.tif",sep="")))
  if(class(x)=="try-error"){
  specoff <- filter(offcombo, SPECIES==as.character(speclist[j]))
  specoff <- distinct(specoff) 
  
  specdat2001 <- filter(PC2001, SPECIES == as.character(speclist[j]))
  specdat2001x <- aggregate(specdat2001$ABUND,by=list("PKEY"=specdat2001$PKEY,"SS"=specdat2001$SS), FUN=sum)
  names(specdat2001x)[3] <- "ABUND"
  dat1 <- right_join(specdat2001x,survey2001,by=c("SS","PKEY")) 
  dat1$SPECIES <- as.character(speclist[j])
  dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
  dat11 <- distinct(dat1,SS,.keep_all=TRUE) #randomly select one survey for analysis 
  s2001 <- left_join(dat11,specoff, by=c("PKEY","SPECIES"))
  d2001 <- left_join(s2001, dat_c, by=c("SS")) 
  
  specdat2011 <- filter(PC2011, SPECIES == as.character(speclist[j])) 
  specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS), FUN=sum)
  names(specdat2011x)[3] <- "ABUND"  
  dat2 <- right_join(specdat2011x,survey2011,by=c("SS","PKEY"))
  dat2$SPECIES <- as.character(speclist[j])
  dat2$ABUND <- as.integer(ifelse(is.na(dat2$ABUND),0,dat2$ABUND)) 
  dat22 <- distinct(dat2,SS,.keep_all=TRUE) #randomly select one survey for analysis 
  s2011 <- left_join(dat22,specoff, by=c("PKEY","SPECIES"))
  d2011 <- left_join(s2011, dat_c, by=c("SS"))

  datcombo <- rbind(d2001,d2011)
  datcombo <- na.omit(datcombo[,c(1:18,20:27,29:30,32:40,42,44)])

  potvar <- datcombo[,c(10:38)]
  var <- get_cn(potvar)
  
  datcombo$wat <- as.factor(datcombo$wat)
  datcombo$urbag <- as.factor(datcombo$urbag)
  #datcombo$landform <- as.factor(datcombo$landform)
  datcombo$wet <- as.factor(datcombo$wet)
  
  x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = var, family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
  if (class(x1) != "NULL") {
    save(brt1,file=paste(w1,speclist[j],"brt5.R",sep=""))
    brtplot(j)
  }
  if(class(x1)=="NULL"){ #retry models that didn't converge with smaller learning rate
    x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = var, family = "poisson", tree.complexity = 3, learning.rate = 0.0001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
    if (class(x1) != "NULL") {
      save(brt1,file=paste(w1,speclist[j],"brt5.R",sep=""))
      brtplot(j)
    }
    if(class(x1)=="NULL"){ #retry models that didn't converge with smaller learning rate
      x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = var, family = "poisson", tree.complexity = 3, learning.rate = 0.00001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
      if (class(x1) != "NULL") {
        save(brt1,file=paste(w1,speclist[j],"brt5.R",sep=""))
        brtplot(j)
      }  
    }
  gc()
   }
  }
}

for (j in 1:length(speclist)) {
  x1 <- try(load(paste(w1,speclist[j],"brt5.R",sep="")))
  if (class(x1) != "try-error") {
  brtplot(j)
  }
}

for (j in 1:length(speclist)) {
  x1 <- try(load(paste(w1,speclist[j],"brt5.R",sep="")))
  if (class(x1) != "try-error") {
    futplot(j)
  }
}

#redo maps
for (j in 1:length(speclist)) {
  x1 <- try(load(paste(w1,speclist[j],"brt5.R",sep="")))
  if (class(x1) != "try-error") {
    mapplot(j)
  }
}

#redo future maps
for (j in 1:length(speclist)) {
  x1 <- try(load(paste(w1,speclist[j],"brt5.R",sep="")))
  if (class(x1) != "try-error") {
    futmapplot(j)
  }
}


cvstats <- cvstatsum(speclist)
write.csv(cvstats,file=paste(w1,"_cvstats5.csv",sep=""))
