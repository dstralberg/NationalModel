library(raster)
library(dismo)
library(gbm)
library(maptools)
library(dplyr)
library(sf)
library(mefa)

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
bcr6 <- shapefile("G:/Boreal/NationalModelsV2/BCR6/bcr6.shp")
p<- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("E:/GIS/hydrology/lakes_lcc.shp")
lc <- crop(l,bcr6)

speclist <- read.csv("E:/BAM/BAMDAta/SpeciesClassesModv5.csv")
speclist <- speclist[speclist$NWT==1,]
speclist <- speclist[,1]
#speclist <- as.factor(c(as.character(speclist),"CAWA","RUBL"))

# bs2001 <- stack(paste(w,"bcr6_2001rasters250.grd",sep=""))
bs2011 <- stack(paste(w,"bcr6_2011rasters250_2.grd",sep=""))
# bs2011_1km <- stack(paste(w,"bcr6_2011rasters1km.grd",sep=""))
bs<-stack("G:/Boreal/NationalModelsV2/bcr6clim_1km.grd")
#bs2 <- stack("G:/Boreal/NationalModelsV2/BCR6/bcr6_2011rasters1km.grd")
bs <- crop(bs,bs2011)
bs <- resample(bs,bs2011)
bs3 <- stack(bs,bs2011)
ARU <- bs3[[1]]*0
names(ARU) <- "ARU"
bs3 <- addLayer(bs3,ARU)

wet250 <- raster("G:/Boreal/NationalModelsV2/BCR6/wet250_BCR6.tif")
wet <- resample(wet250, bs, method="ngb")
writeRaster(wet,file="G:/Boreal/NationalModelsV2/BCR6/wet1km_BCR6.tif")
bs3 <- stack(bs3,wet)
names(bs3[[60]]) <- "wet"

bs2011_1km <- stack(paste(w,"bcr6_2011rasters1km.grd",sep=""))
r2 <- bs2011_1km[[1]]

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
off6 <- read.csv("G:/Boreal/NationalModelsV2/BCR6/BCR6offsets_v3.csv")
off6$PKEY <- as.character(off6$PKEY)
off6$SPECIES <- as.character(off6$SPECIES)

dat2001 <- read.csv("G:/Boreal/NationalModelsV2/BCR6/bcr6_dat2001_v3.csv") #n=40012
dat2001$SS <- as.character(dat2001$SS)
#dat_2001 <- dat2001[!duplicated(dat2001[, 2:3]), ] 
#dat2001$count <- 1

dat2011 <- read.csv("G:/Boreal/NationalModelsV2/BCR6/bcr6_dat2011_V3.csv") #n=33667
dat2011$SS <- as.character(dat2011$SS)
#dat_2011 <- dat2011[!duplicated(dat2011[, 2:3]), ] 
#dat2011$count <- 1

cdat <- read.csv("G:/Boreal/NationalModelsV2/BCR6/bcr6_cdat_v3.csv") 
cdat$SS <- as.character(cdat$SS)

dat2011 <- left_join(dat2011,cdat[,1:30],by=c("SS","X","Y"))
dat2001 <- left_join(dat2001,cdat[,1:30],by=c("SS","X","Y"))

PC2011 <- read.csv("G:/Boreal/NationalModelsV2/BCR6/BCR6PC2011_v3.csv")
PC2011$PKEY <- as.character(PC2011$PKEY)
PC2011$SS <- as.character(PC2011$SS)
PC2001 <- read.csv("G:/Boreal/NationalModelsV2/BCR6/BCR6PC2001_v3.csv")
PC2001$PKEY <- as.character(PC2001$PKEY)
PC2001$SS <- as.character(PC2001$SS)
PC <- read.csv("G:/Boreal/NationalModelsV2/BCR6/BCR6PC_v3.csv")

survey2001 <- aggregate(PC2001$ABUND, by=list("PKEY"=PC2001$PKEY,"SS"=PC2001$SS, "ARU"=PC2001$ARU), FUN=sum) #n=61408
# survey2001 <- survey2001[sample(1:nrow(survey2001)), ]
survey2011 <- aggregate(PC2011$ABUND, by=list("PKEY"=PC2011$PKEY,"SS"=PC2011$SS, "ARU"=PC2011$ARU), FUN=sum) #n=64994
# survey2011 <- survey2011[sample(1:nrow(survey2011)), ]
survey <- rbind(survey2001,survey2011) #n=126402

#calculating sample weights as inverse of survey effort within 5x5 pixel area
samprast2011 <- rasterize(cbind(dat2011$X,dat2011$Y), r2, field=1, fun='sum')
sampsum25 <- focal(samprast2011, w=matrix(1,nrow=5,ncol=5), na.rm=TRUE)
dat_2011 <- cbind(dat2011,extract(sampsum25,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat_2011)[ncol(dat_2011)] <- "sampsum25"
dat_2011$wt <- 1/dat_2011$sampsum25
dat_2011$SS <- as.character(dat_2011$SS) #n=34022

samprast2001 <- rasterize(cbind(dat2001$X,dat2001$Y), r2, field=1, fun='sum')
sampsum25 <- focal(samprast2001, w=matrix(1,nrow=5,ncol=5), na.rm=TRUE)
dat_2001 <- cbind(dat2001,extract(sampsum25,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat_2001)[ncol(dat_2001)] <- "sampsum25"
dat_2001$wt <- 1/dat_2001$sampsum25
dat_2001$SS <- as.character(dat_2001$SS) #n=18837

w <-"G:/Boreal/NationalModelsV2/BCR6/"
setwd(w)

#generate current predictions and plots from models
brtplot <- function (j) {
  load(paste(w,speclist[j],"brt8.R",sep=""))
  varimp <- as.data.frame(brt1$contributions)
  write.csv(varimp,file=paste(w,speclist[j],"varimp8.csv",sep=""))
  cvstats <- as.data.frame(brt1$cv.statistics[c(1,3)])
  cvstats$deviance.null <- brt1$self.statistics$mean.null
  cvstats$deviance.exp <- (cvstats$deviance.null-cvstats$deviance.mean)/cvstats$deviance.null
  write.csv(cvstats,file=paste(w,speclist[j],"cvstats8.csv",sep=""))
  pdf(paste(w,speclist[j],"_plot8.pdf",sep=""))
  gbm.plot(brt1,n.plots=12,smooth=TRUE)
  dev.off()
  rast <- raster::predict(bs3, brt1, type="response", n.trees=brt1$n.trees)
  writeRaster(rast, filename=paste(w,speclist[j],"_pred1km8",sep=""), format="GTiff",overwrite=TRUE)
  
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(w,speclist[j],"_pred1km8.png",sep=""), height=800, width=650)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), maxpixels=5000000, zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(-200000,8800000,"Potential density (males/ha)", cex=1.2)
  dev.off()
  
  PC1 <- PC[PC$SPECIES==as.character(speclist[j]),]
  PC1 <- PC1[PC1$ABUND>0,]
  xy <- PC1[,c(7,8)]
  spdf <- SpatialPointsDataFrame(coords = xy, data = PC1, proj4string = LCC)
  png(file=paste(w,speclist[j],"_pred1km8_pts.png",sep=""), width=650, height=800)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), maxpixels=5000000, zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  plot(spdf, col = 'red', pch=1, cex=0.4, add = TRUE)
  text(-200000,8800000,"Potential density (males/ha)", cex=1.2)
  dev.off()
}

mapplot <- function (j) {
  rast <- raster(paste(w,speclist[j],"_pred1km8.tif",sep=""))
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(w,speclist[j],"_pred1km8.png",sep=""), height=800, width=650)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), maxpixels=5000000, zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(-200000,8800000,"Potential density (males/ha)", cex=1.2)
  dev.off()
  
  PC1 <- PC[PC$SPECIES==as.character(speclist[j]),]
  PC1 <- PC1[PC1$ABUND>0,]
  xy <- PC1[,c(7,8)]
  spdf <- SpatialPointsDataFrame(coords = xy, data = PC1, proj4string = LCC)
  png(file=paste(w,speclist[j],"_pred1km8_pts.png",sep=""), width=650, height=800)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), maxpixels=5000000, zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  plot(spdf, col = 'red', pch=1, cex=0.4, add = TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()
}

cvstatsum <- function (speclist) {
  cvstats <- read.csv(paste(w,speclist[1],"cvstats8.csv",sep=""))
  cvstatmean <- as.data.frame(cbind(as.character(cvstats[,1]),rowMeans(cvstats[,2:6])))
  names(cvstatmean) <- c("stat",as.character(speclist[1]))
  for (j in 2:length(speclist)) {
    x<-try(cv2 <- read.csv(paste(w,speclist[j],"cvstats8.csv",sep="")))
    if(class(x) != "try-error") {
      cvstatmean <- as.data.frame(cbind(cvstatmean,rowMeans(cv2[,2:6])))
      names(cvstatmean)[ncol(cvstatmean)] <- as.character(speclist[j])
    }
  }
  return(cvstatmean)
}

for (j in 1:length(speclist)) {
  x<-try(rast <- raster(paste(w,speclist[j],"_pred1km8.tif",sep="")))
  if(class(x)=="try-error"){
  specoff <- filter(off6, SPECIES==as.character(speclist[j]))
  specoff <- distinct(specoff) 
  
  specdat2001 <- filter(PC2001, SPECIES == as.character(speclist[j]))
  specdat2001x <- aggregate(specdat2001$ABUND,by=list("PKEY"=specdat2001$PKEY,"SS"=specdat2001$SS, "ARU"=specdat2001$ARU), FUN=sum)
  names(specdat2001x)[4] <- "ABUND"
  dat1 <- right_join(specdat2001x,survey2001[,1:3],by=c("SS","PKEY","ARU")) 
  dat1$SPECIES <- as.character(speclist[j])
  dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
  #dat11 <- distinct(dat1,SS,.keep_all=TRUE) #randomly select one survey for analysis
  s2001 <- dplyr::left_join(dat1,specoff, by=c("SPECIES","PKEY"))
  d2001 <- dplyr::left_join(s2001, dat_2001, by=c("SS")) 
  
  specdat2011 <- filter(PC2011, SPECIES == as.character(speclist[j])) 
  specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS, "ARU"=specdat2011$ARU), FUN=sum)
  names(specdat2011x)[4] <- "ABUND"  
  dat2 <- right_join(specdat2011x,survey2011[,1:3],by=c("SS","PKEY","ARU"))
  dat2$SPECIES <- as.character(speclist[j])
  dat2$ABUND <- as.integer(ifelse(is.na(dat2$ABUND),0,dat2$ABUND)) 
  #dat22 <- distinct(dat2,SS,.keep_all=TRUE) 
  s2011 <- left_join(dat2,specoff, by=c("SPECIES","PKEY"))
  d2011 <- left_join(s2011, dat_2011, by=c("SS")) 

  datcombo <- rbind(d2001,d2011)
  datcombo <- na.omit(datcombo[,c(1:6,11:15,18,21:32,35,37,40,42,44:48,50,52,54:59,61:62,64:67,69)])
  
  potvar <- datcombo[,c(3,7:47)]
  var <- get_cn(potvar)
  
  datcombo$wat <- as.factor(datcombo$wat)
  datcombo$urbag <- as.factor(datcombo$urbag)
  #datcombo$landform <- as.factor(datcombo$landform)
  datcombo$wet <- as.factor(datcombo$wet)
  datcombo$ARU <- as.factor(datcombo$ARU)

  x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 4, gbm.x = var, family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
  if (class(x1) != "NULL") {
    save(brt1,file=paste(w,speclist[j],"brt8.R",sep=""))
    brtplot(j)
  }
  if(class(x1)=="NULL"){ #retry models that didn't converge with smaller learning rate
    x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 4, gbm.x = var, family = "poisson", tree.complexity = 3, learning.rate = 0.0001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
    if (class(x1) != "NULL") {
      save(brt1,file=paste(w,speclist[j],"brt8.R",sep=""))
      brtplot(j)
    }
    if(class(x1)=="NULL"){ #retry models that didn't converge with smaller learning rate
      x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 4, gbm.x = var, family = "poisson", tree.complexity = 3, learning.rate = 0.00001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
      if (class(x1) != "NULL") {
        save(brt1,file=paste(w,speclist[j],"brt8.R",sep=""))
        brtplot(j)
      }  
    }
  gc()
    }
  }
}

for (j in 1:length(speclist)) {
  x1 <- try(r<-raster(paste(w,speclist[j],"_pred1km6.tif",sep="")))
  if (class(x1) == "try-error") {
    x1 <- try(load(paste(w,speclist[j],"brt8.R",sep="")))
    if (class(x1) != "try-error") {
      brtplot(j)
    }
  }
}

#redo maps
for (j in 1:length(speclist)) {
  x1 <- try(load(paste(w,speclist[j],"brt8.R",sep="")))
  if (class(x1) != "try-error") {
    mapplot(j)
  }
}

cvstats <- cvstatsum(speclist)
write.csv(cvstats,file=paste(w,"_cvstats8.csv",sep=""))
