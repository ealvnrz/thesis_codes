rm(list=ls())
library(GeoModels)
library(mapproj)
library(fields)

################################################################
###
### Example 1. Simulation of a spatial Kumaraswamy on a regular grid
###
###############################################################
corrmodel="Matern"
sparse=TRUE


setwd("~")


# Define the spatial-coordinates of the points:
x <- seq(0,1,0.009)
y <- seq(0,1,0.009)

model="Kumaraswamy"
shape1=0.5  
shape2=0.05 
param=list(smooth=0.5,min=0,max=1,
             mean=0,sill=1,scale=0.1,nugget=0,shape1=shape1,shape2=shape2)
# Simulation of a spatial Gaussian RF with Matern correlation function
set.seed(89)
data1 <- GeoSim(x,y,grid=TRUE, corrmodel=corrmodel, model=model,param=param,sparse=sparse)$data


#par(mfrow=c(1,2))
pdf("55.pdf")
image.plot(x,y,data1,col=terrain.colors(100),main="",xlab="",ylab="")
dev.off()


pdf("66.pdf")
hist(data1,freq=FALSE,main="",xlab="",ylab="")
dev.off()





mean(data1)
mm=shape1*beta(1+1/shape2,shape1)
mm
var(c(data1))
shape1*beta(1+2/shape2,shape1)-mm^2


         




model="Kumaraswamy"
shape1=1.5
shape2=0.5
param=list(smooth=0.5,min=0,max=1,
             mean=0,sill=1,scale=0.1,nugget=0,shape1=shape1,shape2=shape2)
# Simulation of a spatial Gaussian RF with Matern correlation function
set.seed(89)
data1 <- GeoSim(x,y,grid=TRUE, corrmodel=corrmodel, model=model,param=param,sparse=sparse)$data


#par(mfrow=c(1,2))
pdf("11.pdf")
image.plot(x,y,data1,col=terrain.colors(100),main="",xlab="",ylab="")
dev.off()


pdf("22.pdf")
hist(data1,freq=FALSE,main="",xlab="",ylab="")
dev.off()





model="Kumaraswamy"
shape1=0.5  
shape2=0.5
param=list(smooth=0.5,min=0,max=1,
             mean=0,sill=1,scale=0.05,nugget=0,shape1=shape1,shape2=shape2)
# Simulation of a spatial Gaussian RF with Matern correlation function
set.seed(89)
data1 <- GeoSim(x,y,grid=TRUE, corrmodel=corrmodel, model=model,param=param,sparse=sparse)$data


#par(mfrow=c(1,2))
pdf("33.pdf")
image.plot(x,y,data1,col=terrain.colors(100),main="",xlab="",ylab="")
dev.off()


pdf("44.pdf")
hist(data1,freq=FALSE,main="",xlab="",ylab="")
dev.off()






