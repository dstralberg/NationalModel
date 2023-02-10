library(raster)
library(gbm)
library(dismo)
library(maptools)
library(dplyr)
library(data.table)

bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.5)
provstate <- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
#varimpclasses <- read.csv("G:/Boreal/NationalModelsV2/BCR6/varimpclasses.csv")

speclist <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCspecies_longlist.csv")
speclist <- speclist[,1]

qbs2011_1km <- stack("G:/Boreal/NationalModelsV2/Quebec/QC2011rasters.grd")
r2 <- qbs2011_1km[[1]]

combo2011 <- stack("G:/Boreal/NationalModelsV2/quebec/combo2011.grd")
names(combo2011)[191:194] <- c("TPI","TRI","slope","roughness")
lf <- raster("G:/Boreal/NationalModelsV2/quebec/landform.grd")
combo2011 <- addLayer(combo2011,lf)
names(combo2011)[195] <- "lf"
names(combo2011)[194] <- "rough250"
r <- combo2011[[182]]
r[r>100] <- 100
combo2011 <- stack(combo2011[[1:181]],r,combo2011[[183:nlayers(combo2011)]])

provcrop <- crop(provstate, combo2011[[1]])

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
offlc <- read.csv("G:/Boreal/NationalModelsV2/quebec/QCoffsets_v4.csv")

dat2001 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/quebec_dat2001_v4.csv") #n=5495
dat2001$SS <- as.character(dat2001$SS)
dat2001 <- dat2001[dat2001$nalc<15,] #n=4576

dat2011 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/quebec_dat2011_v4.csv") #n=38748
dat2011$SS <- as.character(dat2011$SS)
dat2011 <- dat2011[dat2011$nalc<15,] #n=33075

#calculating sample weights as inverse of number of survey points within 5x5 pixel radius
samprast2011 <- rasterize(cbind(dat2011$X,dat2011$Y), r2, field=1, fun='sum')
gf <- focalWeight(samprast2011, 25, "Gauss")
sampsum25 <- focal(samprast2011, w=gf, na.rm=TRUE)
dat2011 <- cbind(dat2011,extract(sampsum25,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "sampsum25"
dat2011$wt <- 1/dat2011$sampsum25
dat2011$SS <- as.character(dat2011$SS) #n=42210
rm(samprast2011)

samprast2001 <- rasterize(cbind(dat2001$X,dat2001$Y), r2, field=1, fun='sum')
gf <- focalWeight(samprast2001, 25, "Gauss")
sampsum25 <- focal(samprast2001, w=gf, na.rm=TRUE)
dat2001 <- cbind(dat2001,extract(sampsum25,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "sampsum25"
dat2001$wt <- 1/dat2001$sampsum25
dat2001$SS <- as.character(dat2001$SS) #n=20765
rm(samprast2001)

QCPC2011 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCPC2011_v4.csv") #n=555611
QCPC2011$PKEY <- as.character(QCPC2011$PKEY)
QCPC2011$SS <- as.character(QCPC2011$SS) #n=555611
QCPC2011 <- right_join(QCPC2011[,1:6], dat2011[,1:3], by=c("SS")) #n=466831
QCPC2001 <- read.csv("G:/Boreal/NationalModelsV2/Quebec/QCPC2001_v4.csv") #n=209163
QCPC2001$PKEY <- as.character(QCPC2001$PKEY)
QCPC2001$SS <- as.character(QCPC2001$SS) #n=209163
QCPC2001 <- right_join(QCPC2001[,1:6], dat2001[,1:3], by=c("SS")) #n=166123 

QCPC <- rbind(QCPC2001,QCPC2011) #n=632954

survey2001 <- aggregate(QCPC2001$ABUND, by=list("PKEY"=QCPC2001$PKEY,"SS"=QCPC2001$SS), FUN=sum) #n=20712 (25646 with urban/ag) 
survey2011 <- aggregate(QCPC2011$ABUND, by=list("PKEY"=QCPC2011$PKEY,"SS"=QCPC2011$SS), FUN=sum) #n=43395 (51916 with urban/ag)

w <- "H:/My Drive/BAM.SharedDrive/RshProjs/CC/CCImpacts/QuebecLANDIS/"

w2 <- "G:/Boreal/NationalModelsV2/Quebec/"
setwd(w2)

#generate predictions and plots from models
brtplot <- function (j,PC) {
  e <- new.env()
  f <- file.path(w2,paste0(speclist[j],"brtQC8.RData"))
  load(f,env=e)
  varimp <- as.data.frame(e$brt1$contributions)
  write.csv(varimp,file=paste(w,"v8NoUrban/",speclist[j],"varimp8.csv",sep=""))
  cvstats <- as.data.frame(e$brt1$cv.statistics[c(1,3)])
  cvstats$deviance.null <- e$brt1$self.statistics$mean.null
  cvstats$deviance.exp <- (cvstats$deviance.null-cvstats$deviance.mean)/cvstats$deviance.null
  write.csv(cvstats,file=paste(w,"v8NoUrban/",speclist[j],"cvstats8.csv",sep=""))
  pdf(paste(w,"v8NoUrban/",speclist[j],"_plot8.pdf",sep=""))
  gbm.plot(e$brt1,n.plots=12,smooth=TRUE)
  dev.off()
  rast <- raster::predict(object=combo2011, model=e$brt1, type="response", n.trees=e$brt1$n.trees)
  writeRaster(rast, filename=paste(w,"v8NoUrban/",speclist[j],"_pred1km8",sep=""), format="GTiff",overwrite=TRUE)
  
  q99 <- quantile(rast, probs=c(0.99))	
  prev <- cellStats(rast, 'mean')	
  max <- 3*prev
  png(file=paste(w,"v8NoUrban/",speclist[j],"_pred1km8.png",sep=""), height=800, width=710)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=paste(as.character(speclist[j]),", 1961-1990"), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.5))
  plot(provcrop, col="gray", add=TRUE)
  text(2150000,7950000,"Potential density (males/ha)", cex=1.3)
  dev.off()
  
  PC1 <- PC[PC$ABUND>0,]
  xy1 <- PC1[,c(6,7)]
  spdf1 <- SpatialPointsDataFrame(coords = xy1, data = PC1, proj4string = LCC)
  PC2 <- PC[PC$ABUND==0,]
  xy2 <- PC2[,c(6,7)]
  spdf2 <- SpatialPointsDataFrame(coords = xy2, data = PC2, proj4string = LCC)
  png(file=paste(w,"v8NoUrban/",speclist[j],"_pred1km8_pts.png",sep=""), width=710, height=800)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(as.character(speclist[j]),"current prediction"))
  plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE, main=as.character(speclist[j]), add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.2))
  plot(spdf1, col = 'red', pch=1, cex=0.5, add = TRUE)
  plot(spdf2, col = 'black', pch=1, cex=0.2, add = TRUE)
  plot(provcrop, col="gray", add=TRUE)
  text(2150000,7950000,"Potential density (males/ha)", cex=1.3)
  dev.off()
}

for (j in 1:length(speclist)) {
  x1 <- try(load(paste(w2,speclist[j],"brtQC8.RData",sep="")))
  if (class(x1) != "NULL") {
  specoff <- filter(offlc, SPECIES==as.character(speclist[j]))
  specoff <- distinct(specoff) 
  
  specdat2001 <- filter(QCPC2001, SPECIES == as.character(speclist[j]))
  specdat2001x <- aggregate(specdat2001$ABUND,by=list("PKEY"=specdat2001$PKEY,"SS"=specdat2001$SS), FUN=sum)
  names(specdat2001x)[3] <- "ABUND"
  dat1 <- right_join(specdat2001x,survey2001[,1:3],by=c("SS","PKEY")) 
  dat1$SPECIES <- as.character(speclist[j])
  dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
  dat11 <- distinct(dat1,SS,.keep_all=TRUE) #randomly select one survey for analysis
  s2001 <- left_join(dat11,specoff, by=c("SPECIES","PKEY"))
  d2001 <- left_join(s2001, dat2001, by=c("SS")) 
  
  specdat2011 <- filter(QCPC2011, SPECIES == as.character(speclist[j])) 
  specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS), FUN=sum)
  names(specdat2011x)[3] <- "ABUND"  
  dat2 <- right_join(specdat2011x,survey2011[,1:3],by=c("SS","PKEY"))
  dat2$SPECIES <- as.character(speclist[j])
  dat2$ABUND <- as.integer(ifelse(is.na(dat2$ABUND),0,dat2$ABUND)) 
  dat22 <- distinct(dat2,SS,.keep_all=TRUE) #randomly select one survey for analysis
  s2011 <- left_join(dat22,specoff, by=c("SPECIES","PKEY"))
  d2011 <- left_join(s2011, dat2011, by=c("SS")) 
  
  PC <- rbind(dat11,dat22)
  PC <- left_join(PC,QCPC[,c(1,4,7,8)],by=c("PKEY", "SS"))

  datcombo <- rbind(d2001,d2011)
  datcombo$water <- as.factor(datcombo$wat)
  datcombo$urbag <- as.factor(datcombo$urbag)
  datcombo$lf <- as.factor(datcombo$lf)
  datcombo <- na.omit(datcombo)
  datcombo$Structure_Stand_Age_v1 <- ifelse(datcombo$Structure_Stand_Age_v1 < 100, datcombo$Structure_Stand_Age_v1, 100)
  datcombo$Landsc750_Stand_Age_v1 <- ifelse(datcombo$Landsc750_Stand_Age_v1 < 100, datcombo$Structure_Stand_Age_v1, 100)
  
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
  # landform
  # roughness

  x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 3, gbm.x = c(14,20,22,28,29,35,36,44,52,53,54,58,62,66,69,76,80,83,96,97,107,113,115,121,122,128,129,137,145,146,147,151,155,157,159,162,169,173,176,189,190,197,198,199,200,201,205), family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
  if (class(x1) != "NULL") {
    save(brt1,file=paste(w2,speclist[j],"brtQC8.RData",sep=""))
  }
  }
}


for (j in 1:length(speclist)) {
  x1 <- try(load(paste(w2,speclist[j],"brtQC8.RData",sep="")))
  if (class(x1) != "NULL") {
  specoff <- filter(offlc, SPECIES==as.character(speclist[j]))
  specoff <- distinct(specoff) 
  
  specdat2001 <- filter(QCPC2001, SPECIES == as.character(speclist[j]))
  specdat2001x <- aggregate(specdat2001$ABUND,by=list("PKEY"=specdat2001$PKEY,"SS"=specdat2001$SS), FUN=sum)
  names(specdat2001x)[3] <- "ABUND"
  dat1 <- right_join(specdat2001x,survey2001[,1:3],by=c("SS","PKEY")) 
  dat1$SPECIES <- as.character(speclist[j])
  dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
  dat11 <- distinct(dat1,SS,.keep_all=TRUE) #randomly select one survey for analysis
  s2001 <- left_join(dat11,specoff, by=c("SPECIES","PKEY"))
  d2001 <- left_join(s2001, dat2001, by=c("SS")) 
  
  specdat2011 <- filter(QCPC2011, SPECIES == as.character(speclist[j])) 
  specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS), FUN=sum)
  names(specdat2011x)[3] <- "ABUND"  
  dat2 <- right_join(specdat2011x,survey2011[,1:3],by=c("SS","PKEY"))
  dat2$SPECIES <- as.character(speclist[j])
  dat2$ABUND <- as.integer(ifelse(is.na(dat2$ABUND),0,dat2$ABUND)) 
  dat22 <- distinct(dat2,SS,.keep_all=TRUE) #randomly select one survey for analysis
  s2011 <- left_join(dat22,specoff, by=c("SPECIES","PKEY"))
  d2011 <- left_join(s2011, dat2011, by=c("SS")) 
  
  PC <- rbind(dat11,dat22)
  PC <- left_join(PC,QCPC[,c(1,4,7,8)],by=c("PKEY", "SS"))
  brtplot(j,PC)
  }
  gc()
}


varimpsum <- function (speclist) {
  varimp <- read.csv(paste(w,"v8NoUrban/",speclist[1],"varimp8.csv",sep=""))
  varimp$SPEC <- speclist[1]
  for (j in 2:length(speclist)) {
    x<-try(varimp1 <- read.csv(paste(w,"v8NoUrban/",speclist[j],"varimp8.csv",sep="")))
    if(class(x)!="try-error"){
      varimp1$SPEC <- speclist[j]
      varimp <- rbind(varimp,varimp1)
    }
  }
  #varimp <- merge(varimp,varimpclasses,by="var")
  return(varimp)
}

cvstats2 <- function (speclist) {
  x <- try(load(paste(w,"v8NoUrban/",speclist[1],"brtQC8.RData",sep="")))
  varimp <- read.csv(paste(w,"v8NoUrban/",speclist[1],"varimp8.csv",sep=""))
  cvstats <- as.data.frame(brt1$cv.statistics[c(1,3)])
  cvstats$deviance.null <- brt1$self.statistics$mean.null
  cvstats$deviance.exp <- (cvstats$deviance.null-cvstats$deviance.mean)/cvstats$deviance.null
  cvstats$SPEC <- speclist[1]
  for (j in 2:length(speclist)) {
    x <- try(load(paste(w,"v8NoUrban/",speclist[j],"brtQC8.RData",sep="")))
    if(class(x)!="try-error"){
      varimp <- read.csv(paste(w,"v8NoUrban/",speclist[j],"varimp8.csv",sep=""))
      cvstats1 <- as.data.frame(brt1$cv.statistics[c(1,3)])
      cvstats1$deviance.null <- brt1$self.statistics$mean.null
      cvstats1$deviance.exp <- (cvstats1$deviance.null-cvstats1$deviance.mean)/cvstats1$deviance.null
      cvstats1$SPEC <- speclist[j]
      cvstats <- rbind(cvstats,cvstats1)
    }
  }
  return(cvstats)
}

cvstats <- cvstats2(speclist)
varimp <- varimpsum(speclist)

varimpsummary <- aggregate(varimp[,3],by=list(varimp$SPEC,varimp$var),FUN=sum)
names(varimpsummary)<- c("SPEC","var","rel.inf")
varimpwide <- dcast(varimpsummary, SPEC ~ var)
statscombo <- merge(cvstats,varimpwide,by="SPEC")
write.csv(statscombo,file=paste(w,"v8NoUrban/","_statscombo8.csv",sep=""))


