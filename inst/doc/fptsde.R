## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=4, fig.path='Figs/', fig.show='hold',
                      warning=FALSE, message=FALSE)
library(Sim.DiffProc)

## ------------------------------------------------------------------------
alpha=2
f <- expression( alpha^2 * x )
g <- expression( alpha * x )
mod1d <- snssde1d(drift=f,diffusion=g,x0=0.5,M=1000)

## ------------------------------------------------------------------------
St  <- expression( -5*t+1 )
fpt1d <- fptsde1d(mod1d, boundary = St)
names(fpt1d) ## names of output
summary(fpt1d)

## ----0, echo=FALSE, fig.cap='  ', fig.env='figure*'----------------------
mod <- snssde1d(drift=f,diffusion=g,x0=0.5,M=1)
plot( fptsde1d(mod, boundary = St))

## ----1,fig.env='figure*', fig.cap='  '-----------------------------------
require(MASS)
truehist(fpt1d$fpt,xlab = expression(tau[S(t)])  );box()
lines(density(fpt1d$fpt),col="red",lwd=2)
legend("topright",c("Distribution histogram","Kernel Density"),inset =.01,pch=c(15,NA),lty=c(NA,1),col=c("cyan","red"),lwd=2,cex=0.8)

## ------------------------------------------------------------------------
fx <- expression(5*(-1-y)*x) ; gx <- expression(0.5)
fy <- expression(5*(-1-x)*y) ; gy <- expression(0.5)
mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,x0=2,y0=-2,M=1000)

## ------------------------------------------------------------------------
St <- expression(-3+5*t)
fpt2d <- fptsde2d(mod2d, boundary = St)
names(fpt2d) ## names of output
summary(fpt2d)

## ----00, echo=FALSE, fig.cap='  ', fig.env='figure*'---------------------
mod2 <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,x0=2,y0=-2,M=1)
plot( fptsde2d(mod2, boundary = St))

## ----2,fig.env='figure*', fig.cap='  '-----------------------------------
out <- data.frame(x=fpt2d$fptx,y=fpt2d$fpty)
library(ggplot2)
m <- ggplot(out, aes(x = x, y = y)) 
m + stat_density_2d(aes(fill = ..level..), geom = "polygon")


## ----3, echo=TRUE, fig.cap='  ', fig.env='figure*', message=FALSE, warning=FALSE----
library(sm)
sm.density(out,display="persp")

## ------------------------------------------------------------------------
fx <- expression(4*(-1-x)*y) ; gx <- expression(0.2)
fy <- expression(4*(1-y)*x)  ; gy <- expression(0.2)
fz <- expression(4*(1-z)*y)  ; gz <- expression(0.2) 
mod3d <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,x0=2,y0=-2,z0=0,M=1000)

## ------------------------------------------------------------------------
St <- expression(-3+5*t)
fpt3d <- fptsde3d(mod3d, boundary = St)
names(fpt3d) ## names of output
summary(fpt3d)

## ----000, echo=FALSE, fig.cap='  ', fig.env='figure*'--------------------
mod3 <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,x0=2,y0=-2,z0=0,M=1)
plot( fptsde3d(mod3, boundary = St),pos=3)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  out <- data.frame(x=fpt3d$fptx,y=fpt3d$fpty,z=fpt3d$fptz)
#  library(sm)
#  sm.density(out,display="rgl")

