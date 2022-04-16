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

z=read.csv("data_z.csv")[,2]
coords=read.csv("data_coords.csv")[,2:3]

maxdist=max(dist(coords))
#################################################################
##                      Descriptive plots                      ##
#################################################################
quilt.plot(coords,z)
### Histogram
hist(z)
### Empirical Semivariogram
a<-GeoVariogram(coordx=coords,data=z,maxdist=maxdist/2)
plot(a,pch=20)

##------------------------------------------------------------------------------
##  Is the Beta model a good marginal model?; Analysis assuming independence   -
##------------------------------------------------------------------------------
corrmodel="GenWend"
scale=200;smooth=0;power2=4; nugget=0
copula="Clayton"
nu=6
model="Beta2"
min_1=-1
max_1=1
shape=1
mean=0
optimizer="nlminb"

I=5
fixed<-list(nugget=nugget,sill=1,scale=scale,smooth=smooth,power2=power2,min=-1,max=1,nu=nu)
start<-list(shape=shape,mean=mean)
lower<-list(shape=0,mean=-I)
upper<-list(shape=1000,mean=I)

### Maximum independence likelihood
fit0 <- GeoFit(data=z,coordx=coords,corrmodel=corrmodel, model=model,likelihood="Marginal",type="Independence",
               optimizer=optimizer,lower=lower,upper=upper,copula=copula,
               start=start,fixed=fixed)
GeoQQ(fit0)
GeoQQ(fit0,type="D")

### PIT transformation
aa=GeoPit(fit0,type="Uniform")
hist(fit0$data,freq=FALSE)
GeoScatterplot(aa$data,coords,neighb=c(1,2,3,4))

### Normal score transformation
bb=GeoPit(fit0,type="Gaussian")
hist(bb$data,freq=FALSE)
GeoScatterplot(bb$data,coords,neighb=c(1,2,3,4))

##################################################################
##                  Gaussian Copula Beta Model                  ##
##################################################################

copula="Gaussian"
model="Beta2"
corrmodel="GenWend"
smooth=0

I=7
lower<-list(shape=0,scale=0,mean=-I)
upper<-list(shape=1500,scale=1000000,mean=I)
fixed<-list(sill=1,smooth=smooth,power2=4,min=-1,max=1,nugget=0)
start<-list(mean=0,shape=6,scale=1500)
fit1 <- GeoFit2(data=z,coordx=coords,corrmodel=corrmodel, model=model,
                neighb=2,likelihood="Marginal",type="Pairwise",
                optimizer="nlminb",lower=lower,upper=upper,
                copula=copula,
                start=start,fixed=fixed,sensitivity=TRUE)

##----------------------------------------------------------------
##        Boostrap. Computing std errors and composite AIC       -
##----------------------------------------------------------------

#geovarest1<-GeoVarestbootstrap(fit1,K=250,optimizer="nlminb",lower=lower,upper=upper)
GeoQQ(fit1)


##---------------------------------------------------------------
##                         Correlation                         --
##---------------------------------------------------------------

param=append(fit1$param,fit1$fixed)
x = seq(5,300,5)
corr2= GeoCorrFct_Cop(x=x, corrmodel="GenWend", param=param,copula="Gaussian",model="Beta2", covariance=TRUE)
plot(x,corr2,type="l")
plot(a$centers,a$variograms,ylim=c(0,0.003),pch=20,xlim=c(0,max(a$centers)))

#fix
#vario2= GeoCorrFct_Cop(x=x, corrmodel="GenWend", param=param,copula="Gaussian",model="Beta2", covariance=TRUE,variogram=TRUE)
#lines(x,vario2)


#################################################################
##                  Clayton Copula Beta Model                  ##
#################################################################

copula="Clayton"
optimizer="nlminb"
nu=6   ## nu=2,4,6

I=5
fixed<-list(nu=nu,sill=1,smooth=smooth,power2=4,min=-1,max=1,nugget=0)
start<-list(mean=0,shape=6,scale=200)
lower<-list(shape=0,scale=0,mean=-I)
upper<-list(shape=1500,scale=1000000,mean=I)


fit2 <- GeoFit(data=z,coordx=coords,corrmodel=corrmodel, model=model,
               likelihood="Marginal",type="Pairwise",
               copula=copula,optimizer=optimizer,lower=lower,upper=upper,
               start=start,fixed=fixed,neigh=2,sensitivity = TRUE)

##----------------------------------------------------------------
##        Boostrap. Computing std errors and composite AIC       -
##----------------------------------------------------------------

#geovarest_2<-GeoVarestbootstrap(fit2,K=250,optimizer="nlminb",lower=lower,upper=upper)
GeoQQ(fit2)


##---------------------------------------------------------------
##                         Correlation                         --
##---------------------------------------------------------------

param=append(fit2$param,fit2$fixed)
x = seq(5,300,5)
corr2= GeoCorrFct_Cop(x=x, corrmodel="GenWend", param=param,copula="Clayton",model="Beta2", covariance=TRUE)
plot(x,corr2,type="l")
plot(a$centers,a$variograms,ylim=c(0,0.003),pch=20,xlim=c(0,max(a$centers)))

#fix
#vario2= GeoCorrFct_Cop(x=x, corrmodel="GenWend", param=param,copula="Clayton",model="Beta2", covariance=TRUE,variogram=TRUE)
#lines(x,vario2)
