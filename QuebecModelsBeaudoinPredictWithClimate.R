library(raster)
library(dismo)
library(gbm)
library(maptools)
library(dplyr)

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


bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.8)
provstate <- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")

speclist <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCspecies.csv")
speclist <- speclist[,1]

qbs2011_1km <- brick("G:/Boreal/NationalModelsV2/Quebec/QC2011rasters.grd")
r2 <- qbs2011_1km[[1]]

combo2011 <- brick("G:/Boreal/NationalModelsV2/quebec/combo2011.grd")
names(combo2011)[191:194] <- c("TPI","TRI","slope","roughness")

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
quebec <- raster("E:/GIS/basemaps/quebec250m1.tif")
cur <- "E:/CMIP5/baseline19812010/"
setwd(cur)
clim <- list.files(cur, pattern =".asc$")
curclim<-stack(clim)
qclim <- crop(curclim,quebec)
projection(qclim) <- LCC
qclima <- resample(qclim,combo2011[[1]])
pred <- stack(combo2011,qclima)

dat2001 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/quebec_dat2001_v4.csv") #n=5474
dat2001$SS <- as.character(dat2001$SS)

dat2011 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/quebec_dat2011_v4.csv") #n=38481
dat2011$SS <- as.character(dat2011$SS)

#calculating sample weights as inverse of number of survey points within 5x5 pixel radius
samprast2011 <- rasterize(cbind(dat2011$X,dat2011$Y), r2, field=1, fun='sum')
gf <- focalWeight(samprast2011, 25, "Gauss")
sampsum25 <- focal(samprast2011, w=gf, na.rm=TRUE)
dat2011 <- cbind(dat2011,extract(sampsum25,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "sampsum25"
dat2011$wt <- 1/dat2011$sampsum25
dat2011$SS <- as.character(dat2011$SS) #n=42210
rm(samprast2011)
dat_2011 <- cbind(dat2011, extract(qclima,as.matrix(cbind(dat2011$X,dat_2011$Y))))

samprast2001 <- rasterize(cbind(dat2001$X,dat2001$Y), r2, field=1, fun='sum')
gf <- focalWeight(samprast2001, 25, "Gauss")
sampsum25 <- focal(samprast2001, w=gf, na.rm=TRUE)
dat2001 <- cbind(dat2001,extract(sampsum25,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "sampsum25"
dat2001$wt <- 1/dat2001$sampsum25
dat2001$SS <- as.character(dat2001$SS) #n=20765
rm(samprast2001)
dat_2001 <- cbind(dat2001, extract(qclima,as.matrix(cbind(dat2001$X,dat_2001$Y))))

QCPC2011 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCPC2011_v4.csv") #n=555611
QCPC2011$PKEY <- as.character(QCPC2011$PKEY)
QCPC2011$SS <- as.character(QCPC2011$SS)
QCPC2001 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCPC2001_v4.csv") #n=209163
QCPC2001$PKEY <- as.character(QCPC2001$PKEY)
QCPC2001$SS <- as.character(QCPC2001$SS)

survey2001 <- aggregate(QCPC2001$ABUND, by=list("PKEY"=QCPC2001$PKEY,"SS"=QCPC2001$SS), FUN=sum) #n=25646
survey2011 <- aggregate(QCPC2011$ABUND, by=list("PKEY"=QCPC2011$PKEY,"SS"=QCPC2011$SS), FUN=sum) #n=51916

w <- "G:/Boreal/NationalModelsV2/Quebec/"
setwd(w)

for (j in 1:length(speclist)) {
  x<-try(rast <- raster(paste(w,speclist[j],"_pred1km4.tif",sep="")))
  if(class(x)=="try-error"){
  specoff <- filter(offlc, SPECIES==as.character(speclist[j]))
  specoff <- distinct(specoff) 
  
  specdat2001 <- filter(QCPC2001, SPECIES == as.character(speclist[j]))
  specdat2001x <- aggregate(specdat2001$ABUND,by=list("PKEY"=specdat2001$PKEY,"SS"=specdat2001$SS), FUN=sum)
  names(specdat2001x)[3] <- "ABUND"
  dat1 <- right_join(specdat2001x,survey2001[,1:3],by=c("SS","PKEY")) 
  dat1$SPECIES <- as.character(speclist[j])
  dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
  s2001 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
  d2001 <- left_join(s2001, dat_2001, by=c("SS")) 
  
  specdat2011 <- filter(PC2011, SPECIES == as.character(speclist[j])) 
  specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS), FUN=sum)
  names(specdat2011x)[3] <- "ABUND"  
  dat2 <- right_join(specdat2011x,survey2011[,1:3],by=c("SS","PKEY"))
  dat2$SPECIES <- as.character(speclist[j])
  dat2$ABUND <- as.integer(ifelse(is.na(dat2$ABUND),0,dat2$ABUND)) 
  s2011 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
  d2011 <- left_join(s2011, dat_2011, by=c("SS")) 
  d2011 <- na.omit(d2011) #eliminate non-Quebec data 

  datcombo <- rbind(d2001,d2011)
  datcombo <- na.omit(datcombo)

  potvar <- datcombo[,c(103:107,110:118,120:127,129:130,132:135)]
  var <- get_cn(potvar)
  
  datcombo$wat <- as.factor(datcombo$wat)
  datcombo$urbag <- as.factor(datcombo$urbag)
  datcombo$landform <- as.factor(datcombo$landform)

  x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = var, family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
  if (class(x1) != "NULL") {
    save(brt1,file=paste(w,speclist[j],"brtQC4.R",sep=""))
    varimp <- as.data.frame(brt1$contributions)
    write.csv(varimp,file=paste(w,speclist[j],"varimp4.csv",sep=""))
    cvstats <- t(as.data.frame(brt1$cv.statistics))
    write.csv(cvstats,file=paste(w,speclist[j],"cvstats4.csv",sep=""))
    pdf(paste(w,speclist[j],"_plot4.pdf",sep=""))
    gbm.plot(brt1,n.plots=12,smooth=TRUE)
    dev.off()
    rast <- raster::predict(pred, brt1, type="response", n.trees=brt1$n.trees)
    writeRaster(rast, filename=paste(w,speclist[j],"_pred1km4",sep=""), format="GTiff",overwrite=TRUE)
    
    q99 <- quantile(rast, probs=c(0.99))	
    prev <- cellStats(rast, 'mean')	
    max <- 3*prev
    png(file=paste(w,speclist[j],"_pred1km4.png",sep=""), height=600, width=850)
    par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
    par(mar=c(0,0,5,0))
    plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
    plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=paste(as.character(speclist[j]),", 1961-1990"), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.5))
    plot(provstate, col="gray", add=TRUE)
    text(2400000,7950000,"Potential density (males/ha)", cex=1.3)
    dev.off()
  }
  if(class(x1)=="NULL"){ #retry models that didn't converge with smaller learning rate
    x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = var, family = "poisson", tree.complexity = 3, learning.rate = 0.0001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
    if (class(x1) != "NULL") {
      save(brt1,file=paste(w,speclist[j],"brtQC4.R",sep=""))
      varimp <- as.data.frame(brt1$contributions)
      write.csv(varimp,file=paste(w,speclist[j],"varimp4.csv",sep=""))
      cvstats <- t(as.data.frame(brt1$cv.statistics))
      write.csv(cvstats,file=paste(w,speclist[j],"cvstats4.csv",sep=""))
      pdf(paste(w,speclist[j],"_plot4.pdf",sep=""))
      gbm.plot(brt1,n.plots=12,smooth=TRUE)
      dev.off()
      rast <- raster::predict(pred, brt1, type="response", n.trees=brt1$n.trees)
      writeRaster(rast, filename=paste(w,speclist[j],"_pred1km4",sep=""), format="GTiff",overwrite=TRUE)
      
      q99 <- quantile(rast, probs=c(0.99))	
      prev <- cellStats(rast, 'mean')	
      max <- 3*prev
      png(file=paste(w,speclist[j],"_pred1km4.png",sep=""), height=600, width=850)
      par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
      par(mar=c(0,0,5,0))
      plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
      plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=paste(as.character(speclist[j]),", 1961-1990"), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.5))
      plot(provstate, col="gray", add=TRUE)
      text(2400000,7950000,"Potential density (males/ha)", cex=1.3)
      dev.off()
    }
  }
  gc()
  }
}
