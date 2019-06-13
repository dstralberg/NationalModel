library(raster)
library(dismo)
library(gbm)
library(maptools)
library(dplyr)

bluegreen.colors <- colorRampPalette(c("#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.8)
provstate <- rgdal::readOGR("F:/GIS/basemaps/province_state_line.shp")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
w <-"G:/Boreal/NationalModelsV2/BCR6/"
bcr6 <- shapefile("G:/Boreal/NationalModelsV2/BCR6/bcr6.shp")
p<- rgdal::readOGR("F:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("F:/GIS/hydrology/lakes_lcc.shp")
lc <- crop(l,bcr6)

speclist <- read.csv("F:/BAM/BAMDAta/SpeciesClassesModv5.csv")
speclist <- speclist[speclist$NWT==1|speclist$Alberta==1,]
speclist <- speclist[,1]
#speclist <- as.factor(c(as.character(speclist),"CAWA","RUBL"))

bs2001 <- stack(paste(w,"bcr6_2001rasters250.grd",sep=""))
bs2011 <- stack(paste(w,"bcr6_2011rasters250.grd",sep=""))
bs2011_1km <- stack(paste(w,"bcr6_2011rasters1km.grd",sep=""))
r2 <- bs2011_1km[[1]]

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
offl <- read.csv("G:/Boreal/NationalModelsV2/Quebec/BAMoffsets.csv")
offla <- read.csv("G:/Boreal/NationalModelsV2/Quebec/Atlasoffsets.csv")
offlc <- rbind(offl[2:4],offla[2:4])
offlc$PKEY <- as.character(offlc$PKEY)
offlc$SPECIES <- as.character(offlc$SPECIES)
offlb <- read.csv("G:/Boreal/NationalModelsV2/BCR6/offwt.csv")
offlb$PKEY <- as.character(offlb$PKEY)
offlb$SPECIES <- as.character(offlb$SPECIES)
offld <- read.csv("G:/Boreal/NationalModelsV2/BCR6/offbu.csv")
offld$PKEY <- as.character(offld$PKEY)
offld$SPECIES <- as.character(offld$SPECIES)
offcombo <- rbind(offlc,offlb,offld)

dat2001 <- read.csv("G:/Boreal/NationalModelsV2/BCR6/bcr6_dat2001_v2.csv") #n=18946
dat2001$SS <- as.character(dat2001$SS)
dat_2001 <- dat2001[!duplicated(dat2001[, 2:3]), ] #n=18343
dat_2001$count <- 1

dat2011 <- read.csv("G:/Boreal/NationalModelsV2/BCR6/bcr6_dat2011_V2.csv") #n=33859
dat2011$SS <- as.character(dat2011$SS)
dat_2011 <- dat2011[!duplicated(dat2011[, 2:3]), ] #n=30835
dat_2011$count <- 1

PC2011 <- read.csv(paste(w,"BCR6PC2011_v2.csv",sep=""))
PC2001 <- read.csv(paste(w,"BCR6PC2001_v2.csv",sep=""))

survey2001 <- aggregate(PC2001$ABUND, by=list("PKEY"=PC2001$PKEY,"SS"=PC2001$SS), FUN=sum) 
survey2001 <- survey2001[sample(1:nrow(survey2001)), ]
survey2011 <- aggregate(PC2011$ABUND, by=list("PKEY"=PC2011$PKEY,"SS"=PC2011$SS), FUN=sum) 
survey2011 <- survey2011[sample(1:nrow(survey2011)), ]

#calculating sample weights as inverse of survey effort within 5x5 pixel area
samprast2011 <- rasterize(cbind(dat_2011$X,dat_2011$Y), r2, field=1, fun='sum')
sampsum25 <- focal(samprast2011, w=matrix(1,nrow=5,ncol=5), na.rm=TRUE)
dat_2011 <- cbind(dat_2011,extract(sampsum25,as.matrix(cbind(dat_2011$X,dat_2011$Y))))
names(dat_2011)[ncol(dat_2011)] <- "sampsum25"
dat_2011$wt <- 1/dat_2011$sampsum25
dat_2011$SS <- as.character(dat_2011$SS)

samprast2001 <- rasterize(cbind(dat_2001$X,dat_2001$Y), r2, field=1, fun='sum')
sampsum25 <- focal(samprast2001, w=matrix(1,nrow=5,ncol=5), na.rm=TRUE)
dat_2001 <- cbind(dat_2001,extract(sampsum25,as.matrix(cbind(dat_2001$X,dat_2001$Y))))
names(dat_2001)[ncol(dat_2001)] <- "sampsum25"
dat_2001$wt <- 1/dat_2001$sampsum25
dat_2001$SS <- as.character(dat_2001$SS)

setwd(w)

#generate current predictions and plots from models
brtplot <- function (j) {
  load(paste(w1,speclist[j],"brt3.R",sep=""))
  varimp <- as.data.frame(brt1$contributions)
  write.csv(varimp,file=paste(w,speclist[j],"varimp3.csv",sep=""))
  cvstats <- t(as.data.frame(brt1$cv.statistics))
  write.csv(cvstats,file=paste(w,speclist[j],"cvstats3.csv",sep=""))
  pdf(paste(w1,speclist[j],"_plot3.pdf",sep=""))
  gbm.plot(brt1,n.plots=12,smooth=TRUE)
  dev.off()
  rast <- raster::predict(bs2011, brt1, type="response", n.trees=brt1$n.trees)
  writeRaster(rast, filename=paste(w,speclist[j],"_pred1km3",sep=""), format="GTiff",overwrite=TRUE)
  
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(w,speclist[j],"_pred1km3.png",sep=""), height=800, width=650)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()
  
  PC1 <- PC[PC$SPECIES==as.character(speclist[j]),]
  PC1 <- PC1[PC1$ABUND>0,]
  xy <- PC1[,c(6,7)]
  spdf <- SpatialPointsDataFrame(coords = xy, data = PC1, proj4string = LCC)
  png(file=paste(w,speclist[j],"_pred1km3_pts.png",sep=""), width=650, height=800)
  par(cex.main=1.2, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.2))
  plot(bcr6, border="gray", add=TRUE)
  plot(lc, col="gray", border=NA,add=TRUE)
  plot(spdf, col = 'red', pch=1, cex=0.4, add = TRUE)
  text(-200000,8900000,"Potential density (males/ha)", cex=1.2)
  dev.off()
}

cvstatsum <- function (speclist) {
  cvstats <- read.csv(paste(w,speclist[1],"cvstats3.csv",sep=""))
  cvstatmean <- as.data.frame(cbind(as.character(cvstats[,1]),rowMeans(cvstats[,2:6])))
  names(cvstatmean) <- c("stat",as.character(speclist[1]))
  for (j in 2:length(speclist)) {
    x<-try(cv2 <- read.csv(paste(w,speclist[j],"cvstats3.csv",sep="")))
    if(class(x) != "try-error") {
      cvstatmean <- as.data.frame(cbind(cvstatmean,rowMeans(cv2[,2:6])))
      names(cvstatmean)[ncol(cvstatmean)] <- as.character(speclist[j])
    }
  }
  return(cvstatmean)
}

for (j in 1:length(speclist)) {
  x<-try(rast <- raster(paste(w,speclist[j],"_pred1km3.tif",sep="")))
  if(class(x)=="try-error"){
  specoff <- filter(offcombo, SPECIES==as.character(speclist[j]))
  specoff <- distinct(specoff) 
  
  specdat2001 <- filter(PC2001, SPECIES == as.character(speclist[j]))
  specdat2001x <- aggregate(specdat2001$ABUND,by=list("PKEY"=specdat2001$PKEY,"SS"=specdat2001$SS), FUN=sum)
  names(specdat2001x)[3] <- "ABUND"
  dat1 <- right_join(specdat2001x,survey2001[,1:3],by=c("SS","PKEY")) 
  dat1$SPECIES <- as.character(speclist[j])
  dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) #n=61407
  dat11 <- distinct(dat1,SS,.keep_all=TRUE) #n=14075
  s2001 <- left_join(dat11,specoff, by=c("SPECIES","PKEY"))
  d2001 <- left_join(s2001, dat_2001, by=c("SS")) 
  
  specdat2011 <- filter(PC2011, SPECIES == as.character(speclist[j])) 
  specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS), FUN=sum)
  names(specdat2011x)[3] <- "ABUND"  
  dat2 <- right_join(specdat2011x,survey2011[,1:3],by=c("SS","PKEY"))
  dat2$SPECIES <- as.character(speclist[j])
  dat2$ABUND <- as.integer(ifelse(is.na(dat2$ABUND),0,dat2$ABUND)) #n=58040
  dat22 <- distinct(dat2,SS,.keep_all=TRUE) #n=30531
  s2011 <- left_join(dat22,specoff, by=c("SPECIES","PKEY"))
  d2011 <- left_join(s2011, dat_2011, by=c("SS")) 

  datcombo <- rbind(d2001,d2011)
  datcombo <- na.omit(datcombo)
  datcombo$wat <- as.factor(datcombo$wat)
  datcombo$urbag <- as.factor(datcombo$urbag)
  datcombo$landform <- as.factor(datcombo$landform)

  x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = c(9,11:18,21:26,28:29), family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
  if (class(x1) != "NULL") {
    save(brt1,file=paste(w,speclist[j],"brt3.R",sep=""))
    brtplot(j)
  }
  if(class(x1)=="NULL"){ #retry models that didn't converge with smaller learning rate
    x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = c(9,11:18,21:26,28:29), family = "poisson", tree.complexity = 3, learning.rate = 0.0001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
    if (class(x1) != "NULL") {
      save(brt1,file=paste(w,speclist[j],"brt3.R",sep=""))
      brtplot(j)
    }
    if(class(x1)=="NULL"){ #retry models that didn't converge with smaller learning rate
      x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = c(9,11:18,21:26,28:29), family = "poisson", tree.complexity = 3, learning.rate = 0.00001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
      if (class(x1) != "NULL") {
        save(brt1,file=paste(w,speclist[j],"brt3.R",sep=""))
        brtplot(j)
      }  
    }
  gc()
  }
  }
}

for (j in 1:length(speclist)) {
  x1 <- try(load(paste(w,speclist[j],"brt3.R",sep="")))
  if (class(x1) != "try-error") {
  brtplot(j)
  }
}

cvstats <- cvstatsum(speclist)
write.csv(cvstats,file=paste(w,"_cvstats3.csv",sep=""))
