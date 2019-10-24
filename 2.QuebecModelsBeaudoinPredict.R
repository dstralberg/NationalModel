library(raster)
library(dismo)
library(gbm)
library(maptools)
library(dplyr)

bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.5)
provstate <- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")

speclist <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCspecies.csv")
speclist <- speclist[,1]

qbs2011_1km <- brick("G:/Boreal/NationalModelsV2/Quebec/QC2011rasters.grd")
r2 <- qbs2011_1km[[1]]

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
offlc <- read.csv("G:/Boreal/NationalModelsV2/quebec/QCoffsets_v4.csv")

dat2001 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/quebec_dat2001_v4.csv") #n=5474
dat2001$SS <- as.character(dat2001$SS)

dat2011 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/quebec_dat2011_v4.csv") #n=38481
dat2011$SS <- as.character(dat2011$SS)

#calculating sample weights as inverse of number of survey points within 5x5 pixel radius
samprast2011 <- rasterize(cbind(dat_2011$X,dat_2011$Y), r2, field=1)
gf <- focalWeight(samprast2011, 25, "Gauss")
sampsum25 <- focal(samprast2011, w=gf, na.rm=TRUE)
dat_2011 <- cbind(dat2011,extract(sampsum25,as.matrix(cbind(dat_2011$X,dat_2011$Y))))
names(dat_2011)[ncol(dat_2011)] <- "sampsum25"
dat_2011$wt <- 1/dat_2011$sampsum25
dat_2011$SS <- as.character(dat_2011$SS) #n=42210
rm(samprast2011)

samprast2001 <- rasterize(cbind(dat_2001$X,dat_2001$Y), r2, field=1)
gf <- focalWeight(samprast2001, 25, "Gauss")
sampsum25 <- focal(samprast2001, w=gf, na.rm=TRUE)
dat_2001 <- cbind(dat2001,extract(sampsum25,as.matrix(cbind(dat_2001$X,dat_2001$Y))))
names(dat_2001)[ncol(dat_2001)] <- "sampsum25"
dat_2001$wt <- 1/dat_2001$sampsum25
dat_2001$SS <- as.character(dat_2001$SS) #n=20765
rm(samprast2001)

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

#generate predictions and plots from models
brtplot <- function (j) {
  load(paste(w,speclist[j],"brtQC5.R",sep=""))
  varimp <- as.data.frame(brt1$contributions)
  write.csv(varimp,file=paste(w,speclist[j],"varimp5.csv",sep=""))
  cvstats <- t(as.data.frame(brt1$cv.statistics))
  write.csv(cvstats,file=paste(w,speclist[j],"cvstats5.csv",sep=""))
  pdf(paste(w,speclist[j],"_plot5.pdf",sep=""))
  gbm.plot(brt1,n.plots=12,smooth=TRUE)
  dev.off()
  rast <- raster::predict(qbs2011_1km, brt1, type="response", n.trees=brt1$n.trees)
  writeRaster(rast, filename=paste(w,speclist[j],"_pred1km5",sep=""), format="GTiff",overwrite=TRUE)
  
  q99 <- quantile(rast, probs=c(0.99))	
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(w,speclist[j],"_pred1km5.png",sep=""), height=600, width=850)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=paste(as.character(speclist[j]),", 1961-1990"), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.5))
  plot(provstate, col="gray", add=TRUE)
  text(2400000,7950000,"Potential density (males/ha)", cex=1.3)
  dev.off()
}

for (j in 1:length(speclist)) {
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
  
  specdat2011 <- filter(QCPC2011, SPECIES == as.character(speclist[j])) 
  specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS), FUN=sum)
  names(specdat2011x)[3] <- "ABUND"  
  dat2 <- right_join(specdat2011x,survey2011[,1:3],by=c("SS","PKEY"))
  dat2$SPECIES <- as.character(speclist[j])
  dat2$ABUND <- as.integer(ifelse(is.na(dat2$ABUND),0,dat2$ABUND)) 
  s2011 <- left_join(dat2,specoff, by=c("SPECIES","PKEY"))
  d2011 <- left_join(s2011, dat_2011, by=c("SS")) 
  d2011 <- na.omit(d2011) #eliminate non-Quebec data 

  datcombo <- rbind(d2001,d2011)
  datcombo <- na.omit(datcombo)
  datcombo$wat <- as.factor(datcombo$wat)
  datcombo$urbag <- as.factor(datcombo$urbag)
  datcombo$lf <- as.factor(datcombo$lf)

  # Beaudoin covariates for model 
  # Species_Abie_Bal_v1                                                                                    
  # Species_Acer_Rub_v1                      
  # Species_Acer_Sah_v1                      
  # Species_Betu_All_v1                  
  # Species_Betu_Pap_v1                      
  # Species_Fagu_Gra_v1                      
  # Species_Frax_Ame_v1                      
  # Species_Lari_Lar_v1                      
  # Species_Pice_Gla_v1                      
  # Species_Pice_Mar_v1                      
  # Species_Pice_Rub_v1                                            
  # Species_Pinu_Ban_v1    
  # Species_Pinu_Res_v1                                           
  # Species_Pinu_Str_v1                      
  # Species_Popu_Bal_v1                      
  # Species_Popu_Tre_v1                                           
  # Species_Quer_Rub_v1                                           
  # Species_Thuj_Occ_v1                                         
  # Species_Tsug_Can_v1                                                                                 
  # Structure_Biomass_TotalLiveAboveGround_v1
  # Structure_Stand_Age_v1  
  # dev750
  # led750
  # water
  # urbag
  # TPI
  # TRI
  # slope
  # roughness

  x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = c(12,14,20,22,28,29,35,36,44,52,53,54,58,62,66,69,76,80,83,96,97,198,199,200,201,202,203,204,205), family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
  if (class(x1) != "NULL") {
    save(brt1,file=paste(w,speclist[j],"brtQC5.R",sep=""))
    brtplot(j)
  }
  if(class(x1)=="NULL"){ #retry models that didn't converge with smaller learning rate
    x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = c(12,14,20,22,28,29,35,36,44,52,53,54,58,62,66,69,76,80,83,96,97,198,199,200,201,202,203,204,205), family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
    if (class(x1) != "NULL") {
      save(brt1,file=paste(w,speclist[j],"brtQC5.R",sep=""))
      brtplot(j)
    }
    if(class(x1)=="NULL"){ #retry models that didn't converge with smaller learning rate
      x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = c(12,14,20,22,28,29,35,36,44,52,53,54,58,62,66,69,76,80,83,96,97,198,199,200,201,202,203,204,205), family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
      if (class(x1) != "NULL") {
        save(brt1,file=paste(w,speclist[j],"brtQC5.R",sep=""))
        brtplot(j)
      }  
    }
  gc()
  }
}
