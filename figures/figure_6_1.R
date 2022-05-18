rm(list=ls())

library(cubature)
library(hypergeo)
library(mvtnorm)
library(gcKrig)
library(plot3D)


biv2gauss=function(zi,zj,rho){
  x=zi
  y=zj
  aa=1/(sqrt(1-rho^2)*2*pi)
  tst=exp(-0.5*(1-rho^2)^(-1)*(x^2+y^2-2*rho*x*y))
  out=aa*tst
  return(out)
}

biv2gauss_copula=function(zi,zj,rho){
  x=qnorm(zi)
  y=qnorm(zj)
  aa=exp(0.5*(x^2+y^2))/(sqrt(1-rho^2))
  tst=exp(-0.5*(1-rho^2)^(-1)*(x^2+y^2-2*rho*x*y))
  out=aa*tst
  return(out)
}

####################################################

############# bivariate uniform beta
biv2beta=function(zi,zj,rho,nu){
  k=0;a1=0;res0=0.0;RR=0.0;pp1=0.0;
  x=zi^(2/nu)
  y=zj^(2/nu)
  nu2=nu/2
  cc=nu2+1
  c=((1-rho^2)^(-cc))
  aux=(rho^2)*x*y
  aux1=(rho^2)*(1-x)*(1-y)

  while( k<=10000 ){
    #pp1=log(Re(hypergeo(cc+k,cc+k,nu2,aux)))
    pp1=(nu2-2*(cc+k))*log1p(-aux)+log(Re(hypergeo(nu2-(cc+k),nu2-(cc+k),nu2,aux)))
    bb1=pp1+k*log(aux1)+2*(lgamma(cc+k)-lgamma(cc))-lgamma(k+1)-lgamma(1+k)
    a1= a1 + exp(bb1)
    RR=a1/c
    if((abs(RR-res0)<1e-40)  ) {break;}
    else {res0=RR;}
    k=k+1
  }
  return(RR);
}

############# bivariate uniform beta in Gaussian scale
biv2beta_gauss=function(zi,zj,rho,nu){
  k=0;a1=0;res0=0.0;RR=0.0;pp1=0.0;
  x=pnorm(zi)^(2/nu)
  y=pnorm(zj)^(2/nu)
  nu2=nu/2
  cc=nu2+1
  c=((1-(rho^2))^(-cc))
  aux=(rho^2)*x*y
  aux1=(rho^2)*(1-x)*(1-y)
  while( k<=10000 ){
    #pp1=log(Re(hypergeo(cc+k,cc+k,nu2,aux)))
    pp1=(nu2-2*(cc+k))*log1p(-aux)+log(Re(hypergeo(nu2-(cc+k),nu2-(cc+k),nu2,aux)))
    bb1=pp1+k*log(aux1)+2*(lgamma(cc+k)-lgamma(cc))-lgamma(k+1)-lgamma(1+k)
    a1= a1 + exp(bb1)
    RR=(a1/c)
    if((abs(RR-res0)<1e-40)  ) {break;}
    else {res0=RR;}
    k=k+1
  }
  return(RR*(dnorm(zi)*dnorm(zj)));
}
##############################################
##############################################

biv2beta_vec <- Vectorize(biv2beta, c("zi","zj")) ## bivariate uniform
biv2gauss_copula_vec <- Vectorize(biv2gauss_copula, c("zi","zj")) # uniform gaussian copula
###
biv2beta_gauss_vec <- Vectorize(biv2beta_gauss, c("zi","zj")) #bivunif in scale  gaussiana
biv2gauss_vec <- Vectorize(biv2gauss, c("zi","zj")) #biv gaussiana
###########################


xx=0.5;yy=1

### testing simmetry
biv2gauss(xx,yy,0.7);biv2gauss(-xx,-yy,0.7);
biv2gauss(xx,-yy,0.7);biv2gauss(-xx,yy,0.7);

#nu=2
biv2beta_gauss(xx,yy,0.7,2);biv2beta_gauss(-xx,-yy,0.7,2);
biv2beta_gauss(xx,-yy,0.7,2);biv2beta_gauss(-xx,yy,0.7,2);

#nu=1
biv2beta_gauss(xx,yy,0.7,1);biv2beta_gauss(-xx,-yy,0.7,1);
biv2beta_gauss(xx,-yy,0.7,1);biv2beta_gauss(-xx,yy,0.7,1);

#nu=5
biv2beta_gauss(xx,yy,0.7,5);biv2beta_gauss(-xx,-yy,0.7,5);
biv2beta_gauss(xx,-yy,0.7,5);biv2beta_gauss(-xx,yy,0.7,5);




rho=0.9
bbb=4

##### uniform scale
start=0.001
xx=seq(start,1-start,0.025)
BBbeta=outer(xx, xx, biv2beta_vec,rho=rho,nu=bbb)
BBuni=outer(xx, xx, biv2gauss_copula_vec,rho=rho)
image2D(BBbeta, xx, xx,contour=TRUE,breaks=seq(0, max(BBbeta), by = 0.05))
image2D(BBuni, xx, xx,contour=TRUE,breaks=seq(0, max(BBuni), by = 0.05))
##############################

#### normal scale
ff=seq(-3,3,0.05)
BB2gauss=outer(ff, ff, biv2gauss_vec,rho=rho)
BB2betagauss=outer(ff, ff, biv2beta_gauss_vec,rho=rho,nu=bbb)
BB2betagauss2=outer(ff, ff, biv2beta_gauss_vec,rho=rho,nu=4)
image2D(BB2gauss, ff, ff,contour=TRUE)
image2D(BB2betagauss, ff, ff,contour=TRUE)

#pdf("c02.pdf")
#contour(ff,ff,BB2gauss,lwd=2,levels=seq(0.02,.4,0.04))
contour(ff,ff,BB2betagauss, add = TRUE,lwd=2,col="red",levels=seq(0.02,.4,0.04))
#contour(ff,ff,BB2betagauss2, add = TRUE,lwd=2,col="green",levels=seq(0.02,.4,0.04))
dev.off()



#####################################################

hh=seq(0.0001,0.99,0.001)
corr=exp(-3*hh/0.4)
NN=length(corr)
aa=double(NN)
for(i in 1:NN)
aa[i]=corrTG(beta.gc(shape1 = 1, shape2 = 1), 
         beta.gc(shape1 = 1, shape2 = 1), 
         corrGauss = corr[i], method = "mc", nrep = 600000)
plot(c(0,hh),c(1,aa),type="l",ylim=c(0,1))

rho2=corr^2
aa=((2*(rho2*(3*rho2-1)-(rho2-1)^2*log1p(-rho2)))/(rho2)^2)-3
lines(hh,aa,type="l",ylim=c(0,1),lty=2)
lines(hh,corr,type="l",ylim=c(0,1),lty=3)