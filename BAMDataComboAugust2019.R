library(raster)
library(maptools)
library(dplyr)
library(plyr)
library(tidyr)
#library(data.table)
library(reshape2)

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

PCBBS <- read.delim("F:/BAM/BAMData/BBS/BBS_ptcount_to2017.txt",sep=",",header=TRUE)
names(PCBBS) <- c("SS", "PKEY", "SPECIES","dur","dis","BEHAVIOUR","ABUND","OLD_SPECIES")

load("F:/BAM/BAMData/data_package_2016-04-18.Rdata")	
load("F:/BAM/BAMData/offsets-v3_2016-04-18.Rdata")
coordinates(SS) <- c("X", "Y") 
proj4string(SS) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSBAM <- as.data.frame(spTransform(SS, LCC)) #n=206066 (198736 unique coordinates)
PCBAM <- PCTBL
PCBAM$dis <- as.integer(PCBAM$dis)
PCBAM$dur <- as.integer(PCBAM$dur)
PKEYBAM <- PKEY #n=921275 (all unique)
PKEYBAM$ARU <- 0

load("F:/BAM/BAMData/atlas_data_processed-20181018.RData")
SS <- na.omit(SS)
coordinates(SS) <- c("X", "Y") 
proj4string(SS) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSAtlas <- as.data.frame(spTransform(SS, LCC)) #n=103148 (95338 unique coordinates)
PCAtlas <- PCTBL
PCAtlas$ABUND <- as.integer(PCAtlas$ABUND)
names(PCAtlas) <- c("PCODE", "PKEY","SS","SPECIES","ABUND","BEH","dis","dur")
PKEYAtlas <- PKEY #n=103470 (all unique)
names(PKEYAtlas)[4] <- "YEAR"
PKEYAtlas$ARU <- 0

load("F:/BAM/BamData/ARU/nwt-wildtrax-offsets-2019-01-16.RData")
SSWT <- unique(dd[dd$Y>0,c(33,39:40)]) #n=619
SSWT1 <- na.omit(SSWT) #n=618
coordinates(SSWT1) <- c("X", "Y") 
proj4string(SSWT1) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSWTLC <- as.data.frame(spTransform(SSWT1, LCC))
PKEYWT <- unique(dd[,c(33,34,36)]) #n=2764
PKEYWT$ARU <- 1
PCWT <- melt(y)
names(PCWT) <- c("PKEY","SPECIES","ABUND")
# offWT <- data.table(melt(off))
# names(offWT) <- c("PKEY","SPECIES","logoffset")
# offWT$SPECIES <- as.character(offWT$SPECIES)
# offWT$PKEY <- as.character(offWT$PKEY)

# load("F:/BAM/BamData/ARU/nwt-BU-offsets-2019-01-14.RData")
# SSBU <- unique(dd[,c(14,20:21)])
# coordinates(SSBU) <- c("X", "Y") 
# proj4string(SSBU) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# SSBULC <- as.data.frame(spTransform(SSBU, LCC)) #n=1311 (all unique)
# PKEYBU <- unique(dd[,c(14,15,17)]) #n=4137
# PKEYBU$ARU <- 1
# PCBU <- melt(y)
# names(PCBU) <- c("PKEY","SPECIES","ABUND")
# offBU <- data.table(melt(off))
# names(offBU) <- c("PKEY","SPECIES","logoffset")
# offBU$SPECIES <- as.character(offBU$SPECIES)
# offBU$PKEY <- as.character(offBU$PKEY)

load("F:/BAM/BamData/ARU/nwt-offsets-2019-02-01.RData")
SSBU2 <- unique(dd[,c(14,20:21)])
coordinates(SSBU2) <- c("X", "Y") 
proj4string(SSBU2) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSBU2LC <- as.data.frame(spTransform(SSBU2, LCC)) #n=1763 (1761 unique locations)
# SSBULC <- unique(rbind(SSBULC,SSBU2LC)) #n=1763
PKEYBU2 <- unique(dd[,c(14,15,17)]) #n=5806
PKEYBU2$ARU <- 1
# PKEYBU <- unique(rbind(PKEYBU,PKEYBU2)) #n=5806
PCBU2 <- melt(y)
names(PCBU2) <- c("PKEY","SPECIES","ABUND")
#PCBU <- rbind(PCBU,PCBU2)
# offBU2 <- data.table(melt(off))
# names(offBU2) <- c("PKEY","SPECIES","logoffset")
# offBU2$SPECIES <- as.character(offBU2$SPECIES)
# offBU2$PKEY <- as.character(offBU2$PKEY)
# offBU <- rbind(offBU,offBU2)

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

offcombo <- read.csv("G:/Boreal/NationalModelsV2/BCR6/offcombo.csv")
offcombo <- unique(offcombo, by="PKEY")#n=134790563
SScombo <- rbind(SSBAM[,c(2,48,49)],SSAtlas[,c(1,6,7)],SSWTLC,SSBU2LC) #n=311595
SScombo <- unique(SScombo, by="SS") #n=288662
#SScombo <- unique(SScombo, by=c("X","Y"))
PKEYcombo <- rbind(PKEYBAM[,c(1,2,8,27)],PKEYAtlas[,c(1,2,4,29)],PKEYWT,PKEYBU2) #n=1033315 (19424 duplicates)
PKEYcombo <- unique(PKEYcombo, by="PKEY") #n=1013891
PCcombo <- rbind(PCBAM[,c(2:7)], PCAtlas[,c(2:5,7:8)], PCBBS[,c(1:5,7)])
PCcombo <- unique(PCcombo)
PCcombo <- rbind(PCcombo[,2:4],PCWT, PCBU2)

PCmatch <- left_join(PCcombo,PKEYcombo,by="PKEY")
PCnomatch <- PCmatch[is.na(PCmatch$SS),]
PCsplitnomatch <- separate(data = PCnomatch, col = PKEY, into = c("PCODE", "Transect1","Transect2","Point","Visit"), sep = "\\:")
unique(PCsplitnomatch$PCODE)
#[1] "CL"    "BBSWT" "BBSWA" "BBSYT" "BBSAB" "BBSAK" "BBSDE" "BBSBC" "BBSCT" "BBSID" "BBSIL" "BBSIN"
#[13] "BBSME" "BBSMB" "BBSNT" "USMT"  "USMI"  "USMN"  "CANB"  "CANL"  "BBSNJ" "BBSNH" "BBSNY" "BBSND"
#[25] "BBSOH" "BBSNS" "BBSNU" "BBSON"

PCmatch <- na.omit(PCmatch)

PKEYmatch <- left_join(PKEYcombo,PCcombo,by="PKEY")
PKEYnomatch <- PKEYmatch[is.na(PKEYmatch$ABUND),]
PKEYnomatch <- unique(PKEYnomatch[,1:4])
PKEYsplitnomatch <- separate(data = PKEYnomatch, col = PKEY, into = c("PCODE", "Transect","Point","Year","Visit"), sep = "\\:")
unique(PKEYsplitnomatch$PCODE)
# [1] "ABCAWAWEST"  "AK_YC"       "AKALMS_DIST" "AKCAFS"      "AKDEN"       "AKOR"       
# [7] "AKOR_CNT"    "ALMS"        "BR"          "CF"          "CHSS"        "CL"         
# [13] "CLSS"        "COCL"        "COMW"        "CW"          "DRAP"        "EHPP"       
# [19] "EM"          "EMB-BS"      "EMCLA"       "EMCLA2014"   "FLPC"        "GPMN"       
# [25] "HEDW"        "HOBBBS"      "HR"          "IMBE"        "JL"          "KH"         
# [31] "KTPC"        "LP"          "LR"          "MB06"        "MC"          "ML"         
# [37] "MM94"        "MM95"        "MNBBA"       "MNNFB"       "NLCTRLAB"    "NLLR"       
# [43] "NLMR"        "NWSS"        "PA"          "PF"          "PN"          "PR"         
# [49] "QCATLAS"     "ROMA"        "RP"          "SC"          "SNPC"        "SRDR"       
# [55] "THIN"        "WAP"         "WF"          "WH"          "WR"          "YKFL"       
# [61] "YKRIPARIAN"  "YKSW"        "YKTeslin"    "BBSAB"       "BBSAK"       "BBSBC"      
# [67] "BBSMB"       "BBSMI"       "BBSMN"       "BBSMT"       "BBSNB"       "BBSNL"      
# [73] "BBSNM"       "BBSNS"       "BBSNT"       "BBSNU"       "BBSON"       "BBSPEI"     
# [79] "BBSQC"       "BBSSD"       "BBSSK"       "BBSVT"       "BBSWI"       "BBSYK"   
PKEYmatch <- na.omit(PKEYmatch)
PKEYmatch <- unique(PKEYmatch[,1:4])

rm(dd)
rm(off)
rm(OFF)
rm(offBU)
rm(offBU2)
rm(offWT)
rm(PCBAM)
rm(PCAtlas)
rm(PCWT)
rm(PCBU)
rm(PCBU2)
rm(PCBBS)
rm(YY)
rm(y)
rm(SSBU)
rm(SSBU2)
rm(SSBULC)
rm(SSBU2LC)
rm(SS)
rm(SSWT)
rm(SSWT1)
rm(SSWTLC)
rm(SSAtlas)
rm(SSBAM)
rm(PKEY)
rm(PKEYBU)
rm(PKEYBU2)
rm(PKEYBAM)
rm(PCTBL)
rm(PKEYAtlas)
rm(PKEYWT)

save.image(file="F:/BAM/BAMData/BAM_data_package_November2019.RData")
write.csv(PKEYmatch, file="F:/BAM/BAMData/BAM_data_package_November2019_PKEYmatch.csv")
write.csv(PKEYnomatch, file="F:/BAM/BAMData/BAM_data_package_November2019_PKEYnomatch.csv")
write.csv(PCnomatch, file="F:/BAM/BAMData/BAM_data_package_November2019_PCnomatch.csv")
write.csv(PCmatch, file="F:/BAM/BAMData/BAM_data_package_November2019_PCmatch.csv")


