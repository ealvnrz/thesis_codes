rm(list=ls())
library(GeoModels)
library(mapproj)
library(fields)

# Define the spatial-coordinates of the points:
x <- seq(0,1,0.009)
y <- seq(0,1,0.009)

sparse=FALSE
power2=5
smooth=0
set.seed(261)

data1 <- GeoSimCopula(x,y,grid=TRUE, corrmodel="GenWend", sparse=sparse, 
             model="Gaussian",
             copula="Clayton",
             param=list(power2=power2,smooth=smooth,nu=1,
             mean=0,sill=1,scale=0.15,nugget=0))$data
image.plot(x,y,data1,col=terrain.colors(100),xlab="",ylab="")

set.seed(261)
data2 <- GeoSimCopula(x,y,grid=TRUE, corrmodel="GenWend", 
             model="Gaussian",
             copula="Clayton",
                param=list(power2=power2,smooth=smooth,nu=2,
             mean=0,sill=1,scale=0.15,nugget=0))$data
image.plot(x,y,data2,col=terrain.colors(100),xlab="",ylab="") 

set.seed(261)
data5 <- GeoSimCopula(x,y,grid=TRUE, corrmodel="GenWend", 
             model="Gaussian",
             copula="Clayton",
              param=list(power2=power2,smooth=smooth,nu=5,
             mean=0,sill=1,scale=0.15,nugget=0))$data

setwd("~")

pdf("cl1.pdf")
image.plot(x,y,data1,col=terrain.colors(100),xlab="",ylab="") 
dev.off()
pdf("cl2.pdf") 
image.plot(x,y,data2,col=terrain.colors(100),xlab="",ylab="")
dev.off()  
pdf("cl3.pdf")
image.plot(x,y,data5,col=terrain.colors(100),xlab="",ylab="")   
dev.off()

smooth=1
set.seed(261)
data1 <- GeoSimCopula(x,y,grid=TRUE, corrmodel="GenWend", sparse=sparse, 
             model="Gaussian",
             copula="Clayton",
             param=list(power2=power2,smooth=smooth,nu=1,
             mean=0,sill=1,scale=0.15,nugget=0))$data
image.plot(x,y,data1,col=terrain.colors(100),xlab="",ylab="")


set.seed(261)
data2 <- GeoSimCopula(x,y,grid=TRUE, corrmodel="GenWend", 
             model="Gaussian",
             copula="Clayton",
                param=list(power2=power2,smooth=smooth,nu=2,
             mean=0,sill=1,scale=0.15,nugget=0))$data
image.plot(x,y,data2,col=terrain.colors(100),xlab="",ylab="") 

set.seed(261)
data5 <- GeoSimCopula(x,y,grid=TRUE, corrmodel="GenWend", 
             model="Gaussian",
             copula="Clayton",
              param=list(power2=power2,smooth=smooth,nu=5,
             mean=0,sill=1,scale=0.15,nugget=0))$data
par(mfrow=c(1,3))
pdf("cl4.pdf")
image.plot(x,y,data1,col=terrain.colors(100),xlab="",ylab="")  
dev.off()
pdf("cl5.pdf")
image.plot(x,y,data2,col=terrain.colors(100),xlab="",ylab="")  
dev.off()
pdf("cl6.pdf")
image.plot(x,y,data5,col=terrain.colors(100),xlab="",ylab="") 
dev.off()  

