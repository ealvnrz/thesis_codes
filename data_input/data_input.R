rm(list=ls())
library(GeoModels)
library(raster)
library(fields)
library(sp)
require(rgdal)
require(mapproj)
library(mgcv)

set.seed(385736)

#setwd("")

##################################################################
##                          Data Input                          ##
##################################################################
raw=raster("maxres.TIF")
#plot(raw)
data1 <- as.data.frame(raw, xy=TRUE)
tt=data1[,3]
sel=!is.na(tt)&tt<99998
temp=data1[sel,]

### Raster Subselection 

sel1=temp[,1]> -97&temp[,1]< -92&temp[,2]<  37&temp[,2]>33 
temp2=temp[sel1,]
coords=temp2[,1:2]
z=0.0001*temp2[,3]

### Sinusoidal proyection 

distance="Eucl"
if(distance=="Eucl"){
  P.sinusoidal <- mapproject(coords[,1],coords[,2],projection="sinusoidal")
  coords<-cbind(P.sinusoidal$x,P.sinusoidal$y)*6371
}

### Sample 

sel<-sample(nrow(coords), 1000)
z=z[sel]
coords=coords[sel,]
maxdist=max(dist(coords))

### Save data

#write.csv(z, file = "data_z.csv")
#write.csv(coords, file= "data_coords.csv")