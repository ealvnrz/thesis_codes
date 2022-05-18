library(GeoModels)
library(mapproj)
library(fields)
library(hypergeo)

set.seed(261)
NN=5000
x <- runif(NN);y <- runif(NN)
coords=cbind(x,y)

shape1=5
shape2=2
sill=1
# Replace with 0.18 & 0.1
scale=0.18
# Replace with 0.0.5 & 1.5
smooth=0.5
nugget=0
corrmodel="Matern"
min=0;max=1
param=list(smooth=smooth, min=min,max=max,
           mean=0,sill=sill,scale=scale,nugget=0,shape1=shape1,shape2=shape2)
data1 <- GeoSim(coordx=coords, corrmodel=corrmodel, model="Beta",param=param,
                      sparse=TRUE)$data
quilt.plot(coords,data1)
hist(data1^(shape1/2), freq=FALSE, xlab="", main="")


# Correlation function


corr_beta2<-function(nu,alpha,rho){
  c=(nu+alpha)/2
  nu2=nu/2
  alpha2=alpha/2
  rho2=rho^2
  res=0
  sum=0
  k=0
  while(k<=500){
    p1=2*(lgamma(c+k)-lgamma(c)+lgamma(nu2+1+k)-lgamma(nu2+1))
    p2=lgamma(k+1)+(lgamma(nu2+k)-lgamma(nu2))+2*(lgamma(c+1+k)-lgamma(c+1))
    b1=p1-p2
    b2=log(genhypergeo(U=c(c+k,c+k,alpha2), L=c(c+k+1,c+k+1), check_mod=FALSE, z=rho2))
    b3=k*log(rho2)
    sum=exp(b1+b2+b3)
    res=res+sum
    if (sum<1e-10){
      break
    } else{
      A=res
    }
    k=k+1
  }
  corr=((nu*(c+1))/alpha)*((((1-rho2)^c)*A)-1)
  corr[is.nan(corr)]=0
  return(corr)
}

# Parameters
nu=2
alpha=2
h=seq(0.05,1,0.01)
matern1=exp(-h/scale)
matern2=exp(-h/scale)*(1+(h/scale))

# Replace with matern1 & matern2
corr_v=corr_beta2(nu,alpha,matern1)
corr_v=c(1,corr_v)
exp_rho=c(1,matern1)
h=c(0,h)


plot(h,corr_v, type = 'l', ylim = c(0,1),lty=1,lwd=2,xlab="",ylab="")
lines(h,exp_rho, type = 'l', ylim = c(0,1),lty=2,lwd=2)



