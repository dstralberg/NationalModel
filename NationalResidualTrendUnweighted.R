library(raster)
library(broom)
library(mefa)
library(ggplot2)
load("F:/BAM/BAMdb-GNMsubset-2019-06-20.RData")

# install.packages("devtools")
# devtools::install_github("cardiomoon/ggiraphExtra")
# require(ggplot2)
# require(ggiraph)
# require(ggiraphExtra)
# require(plyr)

speclist <- colnames(yy)
coordinates(dd) <- c("X", "Y") 
proj4string(dd) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
SS <- as.data.frame(spTransform(dd, LCC))
r <- "F:/GoogleDrive/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels/rasters/roadfix-mosaic/"
setwd(r)

for (j in 1:length(speclist)) {
  x<-try(p <- raster(paste("mosaic-",speclist[j],"-roadfix.tif",sep="")))
  if (class(x) != "try-error") {
  obs <- yy[,j]
  pred <- cbind(SS,obs)
  pred1 <- cbind(pred,raster::extract(p,as.matrix(cbind(pred$X,pred$Y))))
  names(pred1)[ncol(pred1)] <- "pred"
  offset <- off[,j]
  pred1 <- cbind(pred1,offset)
  pred1$expoff <- exp(pred1$offset)
  pred1$pred1 <- pred1$pred*pred1$expoff
  pred1$resid <- pred1$pred1-pred1$obs
  pred1 <- na.omit(pred1)
  t <- lm(log(resid+0.5) ~ YEAR*bcrsu, data=pred1)
  beta1 <- as.data.frame(summary(t)$coefficients[,1:2])
  beta1$species <- speclist[j]
  betayr <- mefa:::rep.data.frame(beta1[2,1:2], 15)
  beta2 <- beta1[18:32,]
  beta2$effect <- row.names(beta2)
  beta <- as.data.frame(cbind((beta2[,1] + betayr[,1]),sqrt(beta2[,2]^2 + betayr[,2]^2), beta2[,3:4]))
  names(beta) <- c("coeff","se","species","effect")
  beta$effect <- gsub("YEAR:","",beta3$effect)
  bcr4 <- beta1[2,]
  bcr4$effect <- "bcrsu4"
  names(bcr4) <- c("coeff","se","species","effect")
  beta <- rbind(bcr4,beta)
  write.csv(beta,file=paste(speclist[j],"_betas_subregion_adjust_long_unweighted.csv",sep=""), row.names=FALSE)
  }
}  
  
yr <- 0:20

for (j in 1:length(speclist)) {
    x<-try(b <- read.csv(paste(speclist[j],"_betas_subregion_adjust_long_unweighted.csv",sep="")))
    if (class(x) != "try-error") {
    names(b) <- c("trend","se","species","su")
    b$su <- gsub("bcrsu","",b$su)
    b$su <- as.numeric(gsub("bcsu","",b$su))
    b <- b[order(b$su),]
    # b$su <- gsub("YEAR:","",b$su)
    b1 <- data.frame(yr=yr,rep(b,each=length(yr)))
    b1$chg <- 100*(exp(b1$trend*b1$yr)-1)
    b1$min <- 100*(exp((b1$trend-b1$se)*b1$yr)-1)
    b1$max <- 100*(exp((b1$trend+b1$se)*b1$yr)-1)
    b1$su <- as.factor(b1$su)
    p <- ggplot(b1, aes(x=yr, y=chg, group=su)) + ylim(-50,50) + geom_ribbon(aes(ymin=min,ymax=max,fill=su,alpha=0.3)) + geom_line(aes(color=su)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
    pdf(paste(speclist[j],"_trends_unweighted.pdf",sep="")) 
        plot(p + labs(title=speclist[j]))
        dev.off()
    }
}

# slope1=0.1
# slope2=-0.05
# x1 <- 100*(exp(slope1*yr)-1)
# x2 <- 100*(exp(slope2 *yr )-1)
# d <- data.frame(year=yr, pop=c(x1,x2),what=rep(c("s1", "s2"),each=length(yr)))
# p <- ggplot(d, aes(x=year, y=pop, group=what)) + geom_line(aes(color=what))
# ggPredict(t)
# plot(c(0,1),c(0,1), type = "n", xlab = "x", ylab = "y", asp = 1)
# plot(abline(b=betalong$trend[1]))
# 
# mt <- ggplot(betalong, aes(~ trend, colour = factor(su))) + geom_line()
# 
# mt + facet_wrap( ~ species, scales = "free")

