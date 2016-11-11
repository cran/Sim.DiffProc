## ----setup, include=FALSE,results="asis"---------------------------------
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
fpt1d <- rfptsde1d(mod1d, boundary = St)
summary(fpt1d)

## ----1,fig.env='figure*', fig.cap='  '-----------------------------------
den <- dfptsde1d(mod1d, boundary = St, bw ='ucv')
den 
plot(den)

## ------------------------------------------------------------------------
fx <- expression(5*(-1-y)*x , 5*(-1-x)*y)
gx <- rep(expression(0.5),2)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(x=2,y=-2),M=1000)

## ------------------------------------------------------------------------
St <- expression(-3+5*t)
fpt2d <- rfptsde2d(mod2d, boundary = St)
summary(fpt2d)

## ----2,fig.env='figure*', fig.cap='  '-----------------------------------
denM <- dfptsde2d(mod2d, boundary = St, pdf = 'M')
denM
plot(denM)

## ----3,fig.env='figure*', fig.cap='  '-----------------------------------
denJ <- dfptsde2d(mod2d, boundary = St, pdf = 'J')
denJ
plot(denJ,display="contour",main="Bivariate Density")
plot(denJ,display="image",drawpoints=TRUE,col.pt="green",cex=0.25,pch=19,main="Bivariate Density")

## ----4,webgl=TRUE--------------------------------------------------------
plot(denJ,display="persp",main="Bivariate Density")

## ------------------------------------------------------------------------
fx <- expression(4*(-1-x)*y , 4*(1-y)*x , 4*(1-z)*y) 
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,x0=c(x=2,y=-2,z=0),M=1000)

## ------------------------------------------------------------------------
St <- expression(-3+5*t)
fpt3d <- rfptsde3d(mod3d, boundary = St)
summary(fpt3d)

## ----5,fig.env='figure*', fig.cap='  '-----------------------------------
denM <- dfptsde3d(mod3d, boundary = St)
denM
plot(denM)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  library(sm)
#  sm.density(fpt3d,display="rgl")
#  
#  ##
#  
#  library(ks)
#  fhat <- kde(x=fpt3d)
#  plot(fhat, drawpoints=TRUE)

