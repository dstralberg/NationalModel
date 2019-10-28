library(raster)
library(dismo)
library(cluster)
library(maptools)
library(dplyr)

bluegreen.colors <- colorRampPalette(c("#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.8)
provstate <- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
w <-"G:/Boreal/NationalModelsV2/BCR6/"
bcr6 <- shapefile("G:/Boreal/NationalModelsV2/BCR6/bcr6.shp")
p<- rgdal::readOGR("E:/GIS/basemaps/province_state_line.shp")
l <- rgdal::readOGR("E:/GIS/hydrology/lakes_lcc.shp")
lc <- crop(l,bcr6)

#bs<-stack("G:/Boreal/NationalModelsV2/bcr6clim_1km.grd")
bs2 <- stack("G:/Boreal/NationalModelsV2/BCR6/bcr6_2011rasters1km.grd")
#bs <- crop(bs,bs2)
#bs <- resample(bs,bs2)
#bs3 <- stack(bs,bs2)
rst <- bs2[[c(1:18)]]

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

rstDF <- values(rst)

# Check NA's in the data
idx <- complete.cases(rstDF)

rstKM <- raster(rst[[1]])
rstCLARA <- raster(rst[[1]])

for(nClust in 5:12){
  
  cat("-> Clustering data for nClust =",nClust,"......")
  
  # Perform K-means clustering
  km <- kmeans(rstDF[idx,], centers = nClust, iter.max = 50)
  
  # Perform CLARA's clustering (using manhattan distance)
  cla <- clara(rstDF[idx, ], k = nClust, metric = "manhattan")
  
  # Create a temporary integer vector for holding cluster numbers
  kmClust <- vector(mode = "integer", length = ncell(rst))
  claClust <- vector(mode = "integer", length = ncell(rst))
  
  # Generate the temporary clustering vector for K-means (keeps track of NA's)
  kmClust[!idx] <- NA
  kmClust[idx] <- km$cluster
  
  # Generate the temporary clustering vector for CLARA (keeps track of NA's too ;-)
  claClust[!idx] <- NA
  claClust[idx] <- cla$clustering
  
  # Create a temporary raster for holding the new clustering solution
  # K-means
  tmpRstKM <- raster(rst[[1]])
  # CLARA
  tmpRstCLARA <- raster(rst[[1]])
  
  # Set raster values with the cluster vector
  # K-means
  values(tmpRstKM) <- kmClust
  # CLARA
  values(tmpRstCLARA) <- claClust
  
  # Stack the temporary rasters onto the final ones
  if(nClust==2){
    rstKM    <- tmpRstKM
    rstCLARA <- tmpRstCLARA
  }else{
    rstKM    <- stack(rstKM, tmpRstKM)
    rstCLARA <- stack(rstCLARA, tmpRstCLARA)
  }
  writeRaster(tmpRstKM,paste("G:/Boreal/NationalModelsV2/BCR6/LT8_PGNP_KMeans_nc",nClust,"-1.tif",sep=""), overwrite=TRUE)
  writeRaster(tmpRstCLARA,paste("G:/Boreal/NationalModelsV2/BCR6/LT8_PGNP_CLARA_nc",nClust,"-1.tif",sep=""), overwrite=TRUE)
  cat(" done!\n\n")
}

# Write the clustering solutions for each algorithm
writeRaster(rstKM,"G:/Boreal/NationalModelsV2/BCR6/LT8_PGNP_KMeans_nc2_12-1.tif", overwrite=TRUE)
writeRaster(rstCLARA,"G:/Boreal/NationalModelsV2/BCR6/LT8_PGNP_CLARA_nc2_12-1.tif", overwrite=TRUE)




