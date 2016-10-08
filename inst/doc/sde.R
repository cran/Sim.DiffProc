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

## ----01,fig.env='figure*', fig.cap='  '----------------------------------
x1 <- rsde1d(mod1,at=1) # X(t) at time t = 1 (Ito SDE)
x2 <- rsde1d(mod2,at=1) # X(t) at time t = 1 (Stratonovich SDE)
require(MASS)
truehist(x1,xlab = expression(X[t==1]^mod1));box()
curve(dlnorm(x,meanlog=log(10),sdlog = sqrt(theta^2)),col="red",lwd=2,add=TRUE)
legend("topright",c("Distribution histogram","Exact Log-normal density"),inset =.01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"),lwd=2,cex=0.8)
truehist(x2,xlab = expression(X[t==1]^mod2));box()
curve(dlnorm(x,meanlog=log(10)-0.5*theta^2,sdlog = sqrt(theta^2)),col="red",lwd=2,add=TRUE)
legend("topright",c("Distribution histogram","Exact Log-normal density"),inset = .01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"),lwd=2,cex=0.8)

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
x0=5;y0=0
mu=3;sigma=0.5
fx <- expression(-(x/mu)) ; gx <- expression(sqrt(sigma))
fy <- expression(x) ; gy <- expression(0)
mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,Dt=0.01,M=500,x0=x0,y0=y0,method="smilstein")
mod2d
summary(mod2d)

## ----03,fig.env='figure*', fig.cap=' '-----------------------------------
plot(mod2d)

## ----04,fig.env='figure*', fig.cap='  '----------------------------------
out <- rsde2d(mod2d,at =10)
summary(out)
cov(out) ## Variance-Covariance matrice
E_x <- function(t) x0*exp(-t/mu)
V_x <- function(t) 0.5*sigma*mu *(1-exp(-2*(t/mu)))
E_y <- function(t) y0+x0*mu*(1-exp(-t/mu))
V_y <- function(t) sigma*mu^3*((t/mu)-2*(1-exp(-t/mu))+0.5*(1-exp(-2*(t/mu))))
truehist(out$x,xlab = expression(X[t==10]));box()
curve(dnorm(x,mean=E_x(10) ,sd =sqrt(V_x(10))),col="red",lwd=2,add=TRUE)
legend("topright",c("Distribution histogram","Exact Normal density"),inset = .01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"),lwd=2,cex=0.8)
truehist(out$y,xlab = expression(Y[t==10]));box()
curve(dnorm(x,mean=E_y(10) ,sd =sqrt(V_y(10))),col="red",lwd=2,add=TRUE)
legend("topright",c("Distribution histogram","Exact Normal density"),inset = .01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"),lwd=2,cex=0.8)

## ----05,fig.env='figure*', fig.cap='  '----------------------------------
library(ggplot2)
m <- ggplot(out, aes(x = x, y = y)) 
m + stat_density_2d(aes(fill = ..level..), geom = "polygon")


## ----06,fig.env='figure*', fig.cap='  '----------------------------------
library(sm)
sm.density(out,display="persp")

## ------------------------------------------------------------------------
mu = 4; sigma=0.1
fx <- expression( y ) ; gx <- expression( 0 )
fy <- expression( (mu*( 1-x^2 )* y - x) ) ; gy <- expression( 2*sigma)
mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,N=10000,
                   Dt=0.01,type="str",method="rk1")
mod2d

## ----07,fig.env='figure*', fig.cap='  '----------------------------------
plot2d(mod2d) ## in plane (O,X,Y)
plot(mod2d)   ## back in time

## ------------------------------------------------------------------------
fx <- expression(4*(-1-x)*y) ; gx <- expression(0.2)
fy <- expression(4*(1-y)*x) ; gy <- expression(0.2)
fz <- expression(4*(1-z)*y) ; gz <- expression(0.2)
mod3d <- snssde3d(x0=2,y0=-2,z0=-2,driftx=fx,diffx=gx,drifty=fy,diffy=gy,
               driftz=fz,diffz=gz,N=1000,M=500)
mod3d
summary(mod3d)

## ----08,fig.env='figure*', fig.cap='  '----------------------------------
plot(mod3d,union = TRUE)         ## back in time
plot3D(mod3d,display="persp")    ## in space (O,X,Y,Z)

## ----09,fig.env='figure*', fig.cap='  '----------------------------------
out <- rsde3d(mod3d,at =1)
summary(out)
cov(out) ## Variance-Covariance matrice
truehist(out$x,xlab = expression(X[t==1]));box()
lines(density(out$x),col="red",lwd=2)
legend("topright",c("Distribution histogram","Kernel Density"),inset =.01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"),lwd=2,cex=0.8)
truehist(out$y,xlab = expression(Y[t==1]));box()
lines(density(out$y),col="red",lwd=2)
legend("topright",c("Distribution histogram","Kernel Density"),inset =.01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"),lwd=2,cex=0.8)
truehist(out$z,xlab = expression(Z[t==1]));box()
lines(density(out$z),col="red",lwd=2)
legend("topright",c("Distribution histogram","Kernel Density"),inset =.01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"),lwd=2,cex=0.8)

## ------------------------------------------------------------------------
K = 4; s = 1; sigma = 0.2
fx <- expression( (-K*x/sqrt(x^2+y^2+z^2)) ) ; gx <- expression(sigma)
fy <- expression( (-K*y/sqrt(x^2+y^2+z^2)) ) ; gy <- expression(sigma)
fz <- expression( (-K*z/sqrt(x^2+y^2+z^2)) ) ; gz <- expression(sigma)
mod3d <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,
                   N=10000,x0=1,y0=1,z0=1)
mod3d

## ----10,fig.env='figure*', fig.cap='  '----------------------------------
plot3D(mod3d,display="persp",col="blue")

## ------------------------------------------------------------------------
fx <- expression(y) ; gx <- expression(z)
fy <- expression(0) ; gy <- expression(1)
fz <- expression(0) ; gz <- expression(1)
modtra <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,M=500)
modtra
summary(modtra)

## ----11,fig.env='figure*', fig.cap='  '----------------------------------
plot(modtra$X,plot.type="single",ylab="X")
lines(time(modtra),mean(modtra)$X,col=2,lwd=2)
lines(time(modtra),bconfint(modtra,level=0.95)$X[,1],col=4,lwd=2)
lines(time(modtra),bconfint(modtra,level=0.95)$X[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of",95,"% confidence")),col=c(2,4),lwd=2,cex=0.8)

## ----12,fig.env='figure*', fig.cap='  '----------------------------------
x <- rsde3d(modtra,at=1)$x
summary(x)
truehist(x,xlab = expression(X[t==1]));box()
lines(density(x),col="red",lwd=2)
legend("topright",c("Distribution histogram","Kernel Density"),inset =.01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"),lwd=2,cex=0.8)


