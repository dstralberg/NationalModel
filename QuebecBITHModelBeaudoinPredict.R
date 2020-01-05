library(raster)
library(dismo)
library(gbm)
library(maptools)
library(dplyr)

bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab")
provstate <- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")

qbs2011_1km <- brick("G:/Boreal/NationalModelsV2/Quebec/QC2011rasters.grd")
r2 <- qbs2011_1km[[1]]

combo2011 <- brick("G:/Boreal/NationalModelsV2/quebec/combo2011.grd")
names(combo2011)[191:194] <- c("TPI","TRI","slope","roughness")
lf <- raster("G:/Boreal/NationalModelsV2/quebec/landform.grd")
combo2011 <- addLayer(combo2011,lf)
names(combo2011)[195]<-"lf"

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

dat2001 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/quebec_dat2001_v4.csv") #n=38256
dat2001$SS <- as.character(dat2001$SS)
names(dat2001)[193:200] <- c("dev750","led750","water","urbag","TPI","TRI","slope","roughness")
#dat_2001 <- dat2001
bdat2001 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCBITHdat2001.csv") #n=934
bdat2001$SS <- as.character(bdat2001$SS)
#bdat_2001 <- bdat2001

d2001 <- rbind(dat2001[,1:200],bdat2001) #n=23219
dat_2001 <- d2001[!duplicated(d2001[, 1]), ] #n=23034

dat2011 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/quebec_dat2011_v4.csv") #n=38256
dat2011$SS <- as.character(dat2011$SS)
names(dat2011)[193:200] <- c("dev750","led750","water","urbag","TPI","TRI","slope","roughness")
bdat2011 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCBITHdat2011.csv") #n=1350
bdat2011$SS <- as.character(bdat2011$SS)
#bdat_2011 <- bdat2011

d2011 <- rbind(dat2011[,1:200],bdat2011) #n=39596
dat_2011 <- d2011[!duplicated(d2011[, 1]), ] #n=39569

#calculating sample weights as inverse of number of survey points within 5x5 pixel radius
samprast2011 <- rasterize(cbind(dat_2011$X,dat_2011$Y), r2, field=1, fun='sum')
gf <- focalWeight(samprast2011, 25, "Gauss")
sampsum25 <- focal(samprast2011, w=gf, na.rm=TRUE)
dat_2011 <- cbind(dat_2011,extract(sampsum25,as.matrix(cbind(dat_2011$X,dat_2011$Y))))
names(dat_2011)[ncol(dat_2011)] <- "sampsum25"
dat_2011$wt <- 1/dat_2011$sampsum25
dat_2011$SS <- as.character(dat_2011$SS) #n=42210
rm(samprast2011)

samprast2001 <- rasterize(cbind(dat_2001$X,dat_2001$Y), r2, field=1, fun='sum')
gf <- focalWeight(samprast2001, 25, "Gauss")
sampsum25 <- focal(samprast2001, w=gf, na.rm=TRUE)
dat_2001 <- cbind(dat_2001,extract(sampsum25,as.matrix(cbind(dat_2001$X,dat_2001$Y))))
names(dat_2001)[ncol(dat_2001)] <- "sampsum25"
dat_2001$wt <- 1/dat_2001$sampsum25
dat_2001$SS <- as.character(dat_2001$SS) #n=20765
rm(samprast2001)

BPC2011 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/BITHPC2011.csv") #n=1361
BPC2011$PKEY <- as.character(BPC2011$PKEY)
BPC2011$SS <- as.character(BPC2011$SS)
QCPC2011 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCPC2011_v4.csv") #n=252792
QCPC2011$PKEY <- as.character(QCPC2011$PKEY)
QCPC2011$SS <- as.character(QCPC2011$SS)

BPC2001 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/BITHPC2001.csv") #n=1361
BPC2001$PKEY <- as.character(BPC2001$PKEY)
BPC2001$SS <- as.character(BPC2001$SS)
QCPC2001 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCPC2001_v4.csv") #n=212901
QCPC2001$PKEY <- as.character(QCPC2001$PKEY)
QCPC2001$SS <- as.character(QCPC2001$SS)

PC2001 <- rbind(BPC2001[,c(1:2,7:8)],QCPC2001[,1:4]) 
PC2001$PKEY <- as.character(PC2001$PKEY)
PC2001$SS <- as.character(PC2001$SS)
PC2011 <- rbind(BPC2011[,c(1:2,7:8)],QCPC2011[,1:4]) #n=882940
PC2011$PKEY <- as.character(PC2011$PKEY)
PC2011$SS <- as.character(PC2011$SS)

survey2001 <- aggregate(PC2001$ABUND, by=list("PKEY"=PC2001$PKEY,"SS"=PC2001$SS), FUN=sum) #n=26161
survey2011 <- aggregate(PC2011$ABUND, by=list("PKEY"=PC2011$PKEY,"SS"=PC2011$SS), FUN=sum) #n=89289

w <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/CC/CCImpacts/QuebecLANDIS/"
setwd(w)

#generate predictions and plots from model
brtplot <- function (brt1) {
   varimp <- as.data.frame(brt1$contributions)
   write.csv(varimp,file=paste(w,"BITHvarimp7.csv",sep=""))
   cvstats <- t(as.data.frame(brt1$cv.statistics))
   write.csv(cvstats,file=paste(w,"BITHcvstats7.csv",sep=""))
   pdf(paste(w,"BITH_plot5.pdf",sep=""))
   gbm.plot(brt1,n.plots=12,smooth=TRUE)
   dev.off()
   rast <- raster::predict(combo2011, brt1, type="response", n.trees=brt1$n.trees)
   writeRaster(rast, filename=paste(w,"BITH_pred1km7",sep=""), format="GTiff",overwrite=TRUE)

   png(file=paste(w,"BITH_pred1km7.png",sep=""), height=800, width=710)
   par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
   par(mar=c(0,0,5,0))
   plot(rast, col="blue", axes=FALSE, legend=FALSE, main="BITH current prediction")
   plot(rast, col=bluegreen.colors(15), zlim=c(0,0.5), axes=FALSE, main=paste(as.character(speclist[j]),", 1961-1990"), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.5))
   plot(provstate, col="gray", add=TRUE)
   text(2400000,7950000,"probability of occurrence", cex=1.3)
   dev.off()
   }

specdat2001 <- filter(PC2001, SPECIES == "BITH")
specdat2001x <- aggregate(specdat2001$ABUND,by=list("PKEY"=specdat2001$PKEY,"SS"=specdat2001$SS), FUN=sum)
names(specdat2001x)[3] <- "ABUND"
dat1 <- right_join(specdat2001x,survey2001[,1:3],by=c("SS","PKEY")) 
dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND))
dat11 <- distinct(dat1,SS,.keep_all=TRUE) #randomly select one survey for analysis
d2001 <- left_join(dat11, dat_2001, by=c("SS"))
  
specdat2011 <- filter(PC2011, SPECIES == "BITH") 
specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS), FUN=sum)
names(specdat2011x)[3] <- "ABUND"  
dat2 <- right_join(specdat2011x,survey2011[,1:3],by=c("SS","PKEY"))
dat2$ABUND <- as.integer(ifelse(is.na(dat2$ABUND),0,dat2$ABUND)) 
dat22 <- distinct(dat2,SS,.keep_all=TRUE) #randomly select one survey for analysis
d2011 <- left_join(dat22, dat_2011, by=c("SS")) 

datcombo <- rbind(d2001,d2011)
datcombo <- na.omit(datcombo) #n=78923
datcombo$water <- as.factor(datcombo$wat)
datcombo$urbag <- as.factor(datcombo$urbag)
datcombo$lf <- as.factor(datcombo$lf)
datcombo$ABUND <- ifelse(datcombo$ABUND>0,1,0)

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
  # Landsc750_Abie_Bal_v1                                                                                    
  # Landsc750_Acer_Rub_v1                      
  # Landsc750_Acer_Sah_v1                      
  # Landsc750_Betu_All_v1                  
  # Landsc750_Betu_Pap_v1                      
  # Landsc750_Fagu_Gra_v1                      
  # Landsc750_Frax_Ame_v1                      
  # Landsc750_Lari_Lar_v1                      
  # Landsc750_Pice_Gla_v1                      
  # Landsc750_Pice_Mar_v1                      
  # Landsc750_Pice_Rub_v1                                            
  # Landsc750_Pinu_Ban_v1    
  # Landsc750_Pinu_Res_v1                                           
  # Landsc750_Pinu_Str_v1                      
  # Landsc750_Popu_Bal_v1                      
  # Landsc750_Popu_Tre_v1                                           
  # Landsc750_Quer_Rub_v1                                           
  # Landsc750_Thuj_Occ_v1                                         
  # Landsc750_Tsug_Can_v1                                                                                 
  # Landsc750_Biomass_TotalLiveAboveGround_v1
  # Landsc750_Stand_Age_v1  
  # dev750
  # led750
  # water
  # urbag
  # landcover
  # roughness

x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = c(12,18,20,26,27,33,34,42,50,51,52,56,60,62,64,67,74,78,81,94,95,105,111,113,119,120,126,127,135,143,144,145,149,153,155,157,160,167,171,174,187,188,195,196,197,198,199,203), family = "bernoulli", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, site.weights=datcombo$wt))
  if (class(x1) != "NULL") {
    save(brt1,file=paste(w,"BITHbrtQC7.R",sep=""))
    brtplot(brt1)
  }


