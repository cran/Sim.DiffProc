## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=4, fig.path='Figs/', fig.show='hold',
                      warning=FALSE, message=FALSE)
library(Sim.DiffProc)

## ------------------------------------------------------------------------
theta = 0.5
f <- expression( (0.5*theta^2*x) )
g <- expression( theta*x )
mod1 <- snssde1d(drift=f,diffusion=g,x0=10,M=500,type="ito") # Using Ito
mod2 <- snssde1d(drift=f,diffusion=g,x0=10,M=500,type="str") # Using Stratonovich 
mod1
mod2

## ------------------------------------------------------------------------
summary(mod1, at = 1)
summary(mod2, at = 1)

## ------------------------------------------------------------------------
x1 <- rsde1d(object = mod1, at = 1)  # X(t=1) | X(0)=x0 (Itô SDE)
x2 <- rsde1d(object = mod2, at = 1)  # X(t=1) | X(0)=x0 (Stratonovich SDE)
summary(data.frame(x1,x2))

## ----01,fig.env='figure*', fig.cap='  '----------------------------------
mu1 = log(10); sigma1= sqrt(theta^2)  # log mean and log variance for mod1 
mu2 = log(10)-0.5*theta^2 ; sigma2 = sqrt(theta^2) # log mean and log variance for mod2
AppdensI <- dsde1d(mod1, at = 1)
AppdensS <- dsde1d(mod2, at = 1)
plot(AppdensI , dens = function(x) dlnorm(x,meanlog=mu1,sdlog = sigma1))
plot(AppdensS , dens = function(x) dlnorm(x,meanlog=mu2,sdlog = sigma2))

## ----02,fig.env='figure*', fig.cap=' ',fig.show='hold'-------------------
plot(mod1,plot.type="single",ylab=expression(X^mod1))
lines(time(mod1),mean(mod1),col=2,lwd=2)
lines(time(mod1),bconfint(mod1,level=0.95)[,1],col=4,lwd=2)
lines(time(mod1),bconfint(mod1,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of",95,"% confidence")),col=c(2,4),lwd=2,cex=0.8)
plot(mod2,plot.type="single",ylab=expression(X^mod2))
lines(time(mod2),mean(mod2),col=2,lwd=2)
lines(time(mod2),bconfint(mod2,level=0.95)[,1],col=4,lwd=2)
lines(time(mod2),bconfint(mod2,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of",95,"% confidence")),col=c(2,4),lwd=2,cex=0.8)

## ------------------------------------------------------------------------
x=5;y=0
mu=3;sigma=0.5
fx <- expression(-(x/mu),x)  
gx <- expression(sqrt(sigma),0)
mod2d <- snssde2d(drift=fx,diffusion=gx,Dt=0.01,M=500,x0=c(x,y),method="smilstein")
mod2d
summary(mod2d)

## ----03,fig.env='figure*', fig.cap=' '-----------------------------------
plot(mod2d)

## ------------------------------------------------------------------------
out <- rsde2d(object = mod2d, at = 10) 
summary(out)

## ----04,fig.env='figure*', fig.cap='  '----------------------------------
denM <- dsde2d(mod2d,pdf="M",at =10)
denM
plot(denM, main="Marginal Density")

## ----05,fig.env='figure*', fig.cap='  '----------------------------------
denJ <- dsde2d(mod2d,pdf="J",at =10)
denJ
plot(denJ,display="contour",main="Bivariate Density")
plot(denJ,display="image",drawpoints=TRUE,col.pt="green",cex=0.25,pch=19,main="Bivariate Density")

## ----06,fig.env='figure*', fig.cap='  '----------------------------------
plot(denJ,main="Bivariate Density")

## ------------------------------------------------------------------------
mu = 4; sigma=0.1
fx <- expression( y ,  (mu*( 1-x^2 )* y - x)) 
gx <- expression( 0 ,2*sigma)
mod2d <- snssde2d(drift=fx,diffusion=gx,N=10000,Dt=0.01,type="str",method="rk1")
mod2d

## ----07,fig.env='figure*', fig.cap='  '----------------------------------
plot2d(mod2d) ## in plane (O,X,Y)
plot(mod2d)   ## back in time

## ------------------------------------------------------------------------
fx <- expression(4*(-1-x)*y , 4*(1-y)*x , 4*(1-z)*y) 
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(x0=c(x=2,y=-2,z=-2),drift=fx,diffusion=gx,N=1000,M=500)
mod3d
summary(mod3d)

## ----08,fig.env='figure*', fig.cap='  '----------------------------------
plot(mod3d,union = TRUE)         ## back in time
plot3D(mod3d,display="persp")    ## in space (O,X,Y,Z)

## ----09,fig.env='figure*', fig.cap='  '----------------------------------
den <- dsde3d(mod3d,pdf="M",at =1)
den
plot(den, main="Marginal Density") 

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  denJ <- dsde3d(mod3d,pdf="J")
#  denJ
#  plot(denJ,display="persp")
#  mtext("Contour surface of kernel density estimate")

## ------------------------------------------------------------------------
K = 4; s = 1; sigma = 0.2
fx <- expression( (-K*x/sqrt(x^2+y^2+z^2)) , (-K*y/sqrt(x^2+y^2+z^2)) , (-K*z/sqrt(x^2+y^2+z^2)) ) 
gx <- rep(expression(sigma),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,N=10000,x0=c(x=1,y=1,z=1))
mod3d

## ----11,fig.env='figure*', fig.cap='  '----------------------------------
plot3D(mod3d,display="persp",col="blue")

## ------------------------------------------------------------------------
fx <- expression(y,0,0) 
gx <- expression(z,1,1)
modtra <- snssde3d(drift=fx,diffusion=gx,M=500)
modtra
summary(modtra)

## ----12,fig.env='figure*', fig.cap='  '----------------------------------
plot(modtra$X,plot.type="single",ylab="X")
lines(time(modtra),mean(modtra)$X,col=2,lwd=2)
lines(time(modtra),bconfint(modtra,level=0.95)$X[,1],col=4,lwd=2)
lines(time(modtra),bconfint(modtra,level=0.95)$X[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of",95,"% confidence")),col=c(2,4),lwd=2,cex=0.8)

## ----13,fig.env='figure*', fig.cap='  '----------------------------------
den <- dsde3d(modtra,at=1)
den$resx
MASS::truehist(den$ech$x,xlab = expression(X[t==1]));box()
lines(den$resx,col="red",lwd=2)
legend("topleft",c("Distribution histogram","Kernel Density"),inset =.01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"),lwd=2,cex=0.8)


