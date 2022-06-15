rm(list=ls())
library(GeoModels)
library(raster)
library(fields)
library(sp)
require(rgdal)
require(mapproj)
library(mgcv)
library(ggplot2)
library(gcKrig)

set.seed(385736)
seed=385736

setwd("~/")

z=read.csv("data_z.csv")[,2]
coords=read.csv("data_coords.csv")[,2:3]
maxdist=max(dist(coords))
#################################################################
##                      Descriptive plots                      ##
#################################################################

quilt.plot(coords,z)

### Empirical Semivariogram
a<-GeoVariogram(coordx=coords,data=z,maxdist=maxdist/3)
plot(a,pch=20,ylim=c(0,0.0025),xlab="Distance",ylab="Semivariogram")

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

### Histogram
hist(z,xlim=c(-1,1), xlab="NDVI", main="Histogram of NDVI", prob=TRUE, ylim=c(0,16))
GeoQQ(fit0,type="D",ylim=c(0,10))

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

#geovarest1<-GeoVarestbootstrap(fit1,K=500,optimizer="nlminb",lower=lower,upper=upper,seed=seed)
GeoQQ(fit1)

##---------------------------------------------------------------
##                         Correlation                         --
##---------------------------------------------------------------

param=append(fit1$param,fit1$fixed)
x=seq(0.001,300,10)
semi=c(rep(0,length(x)))
exp_rho=(1-x/(112.4101))^4*I(x<=112.4101)
MM=fit0$param$mean
mm=1/(1+exp(-MM))
sh=fit0$param$shape
pmin=fit0$fixed$min;pmax=fit0$fixed$max;
shape1=mm*sh
shape2=(1-mm)*sh
a1=mm*sh;a2=mm*sh
b1=(1-mm)*sh;b2=(1-mm)*sh
e1=a1/(a1+b1); e2=a2/(a2+b2)   
v1=a1*b1/((a1+b1)^2*(a1+b1+1));v2=a2*b2/((a2+b2)^2*(a2+b2+1))
min=pmin;max=pmax
dd=max-min
vs= sqrt(v1*v2)*dd^2
for (i in 1:length(x)){
  corr_g<-corrTG(marg1 =beta.gc(shape1 = shape1, shape2 = shape2) , marg2 = beta.gc(shape1 = shape1, shape2 =shape2) , 
       corrGauss = exp_rho[i], method = "integral")
  semi[i]=vs*(1-corr_g)
}


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

#geovarest_2<-GeoVarestbootstrap(fit2,K=500,optimizer="nlminb",lower=lower,upper=upper,seed=seed)
GeoQQ(fit2)

##---------------------------------------------------------------
##                         Correlation                         --
##---------------------------------------------------------------

param=append(fit2$param,fit2$fixed)
corr2= GeoCorrFct_Cop(x=x, corrmodel="GenWend", param=param,copula="Clayton",model="Beta2", covariance=TRUE)
plot(x,corr2,type="l")
plot(a$centers,a$variograms,ylim=c(0,0.003),pch=20,xlim=c(0,max(a$centers)),xlab="Distance",ylab="Semivariogram")

vario_clayton= GeoCorrFct_Cop(x=x, corrmodel="GenWend", param=param,copula="Clayton",model="Beta2", covariance=TRUE,variogram=TRUE)
lines(x,vario_clayton)
lines(x,semi,col="red")

