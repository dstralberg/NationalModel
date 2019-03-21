CTI <- function(tempdd) {
w <- "F:/GIS/topoedaphic//GLCFDEM/"
setwd(w)
wet <- list.files(w, pattern="wetness.asc$")
setwd(w)
tempdd$cti90 <- 0
r <- raster(wet[1])
lat <- as.numeric(substr(wet[1],10,11))
lon <- as.numeric(substr(wet[1],13,15))*-1

cti <- tempdd[floor(tempdd$lon) == lon,]
cti <- cti[floor(cti$lat) == lat,]
if (nrow(cti)>0) {cti$cti90 <- extract(r,st_coordinates(x))}
for (i in 2:length(wet)) {
  r <- raster(wet[i])
  lat <- as.numeric(substr(wet[i],10,11))
  lon <- as.numeric(substr(wet[i],13,15))*-1
  x <- tempdd[floor(tempdd$lon) == lon,]
  x <- x[floor(x$lat) == lat,]
  x$cti90 <- extract(r,st_coordinates(x))
  cti <- rbind(cti, x)
}

cti1 <- tempdd[ceiling(tempdd$lon) == lon,]
cti1 <- cti1[floor(cti1$lat) == lat,]
if (nrow(cti1)>0) {cti1$cti90 <- extract(r,st_coordinates(x))}
for (i in 2:length(wet)) {
  r <- raster(wet[i])
  lat <- as.numeric(substr(wet[i],10,11))
  lon <- as.numeric(substr(wet[i],13,15))*-1
  x <- tempdd[ceiling(tempdd$lon) == lon,]
  x <- x[floor(x$lat) == lat,]
  x$cti90 <- extract(r,st_coordinates(x))
  cti1 <- rbind(cti1, x)
}	

cti2 <- tempdd[ceiling(tempdd$lon) == lon,]
cti2 <- cti2[ceiling(cti2$lat) == lat,]
if (nrow(cti2)>0) {cti2$cti90 <- extract(r,st_coordinates(x))}
for (i in 2:length(wet)) {
  r <- raster(wet[i])
  lat <- as.numeric(substr(wet[i],10,11))
  lon <- as.numeric(substr(wet[i],13,15))*-1
  x <- tempdd[ceiling(tempdd$lon) == lon,]
  x <- x[ceiling(x$lat) == lat,]
  x$cti90 <- extract(r,st_coordinates(x))
  cti2 <- rbind(cti2, x)
}		

cti3 <- tempdd[floor(tempdd$lon) == lon,]
cti3 <- cti3[ceiling(cti3$lat) == lat,]
if (nrow(cti3)>0) {cti3$cti90 <- extract(r,st_coordinates(x))}
for (i in 2:length(wet)) {
  r <- raster(wet[i])
  lat <- as.numeric(substr(wet[i],10,11))
  lon <- as.numeric(substr(wet[i],13,15))*-1
  x <- tempdd[floor(tempdd$lon) == lon,]
  x <- x[ceiling(x$lat) == lat,]
  x$cti90 <- extract(r,st_coordinates(x))
  cti3 <- rbind(cti3, x)
}		

cticombo <- rbind(cti, cti1, cti2, cti3)	
cticombo <- cticombo[is.na(cticombo$cti90)==FALSE,]
return(cticombo)
}