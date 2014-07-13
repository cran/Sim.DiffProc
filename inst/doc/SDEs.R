### R code from vignette source 'SDEs.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: SDEs.Rnw:131-132 (eval = FALSE)
###################################################
## install.packages("Sim.DiffProc")


###################################################
### code chunk number 2: packages
###################################################
library(Sim.DiffProc)


###################################################
### code chunk number 3: SDEs.Rnw:140-141 (eval = FALSE)
###################################################
## library(help = "Sim.DiffProc")


###################################################
### code chunk number 4: SDEs.Rnw:188-194
###################################################
f <- expression( (0.5*0.5^2*x) )
g <- expression( 0.5*x )
mod1 <- snssde1d(drift=f,diffusion=g,x0=10,M=1,N=1000)
mod2 <- snssde1d(drift=f,diffusion=g,x0=10,M=1,N=1000,type="str")
mod1
mod2


###################################################
### code chunk number 5: SDEs.Rnw:197-199 (eval = FALSE)
###################################################
## plot(mod1)
## plot(mod2)


###################################################
### code chunk number 6: SDEs.Rnw:204-205
###################################################
plot(mod1)


###################################################
### code chunk number 7: SDEs.Rnw:207-208
###################################################
plot(mod2)


###################################################
### code chunk number 8: SDEs.Rnw:225-229
###################################################
mod1 <- snssde1d(drift=f,diffusion=g,x0=10,M=100,N=1000)
mod2 <- snssde1d(drift=f,diffusion=g,x0=10,M=100,N=1000,type="str")
summary(mod1)
summary(mod2)


###################################################
### code chunk number 9: SDEs.Rnw:234-249 (eval = FALSE)
###################################################
## plot(mod1,plot.type="single")
## lines(time(mod1),mean(mod1),col=2,lwd=2)
## lines(time(mod1),bconfint(mod1,level=0.95)[,1],col=4,lwd=2)
## lines(time(mod1),bconfint(mod1,level=0.95)[,2],col=4,lwd=2)
## legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),
##        inset = .01,col=c(2,4),lwd=2,cex=0.8)
##        
## dev.new()
## 
## plot(mod2,plot.type="single")
## lines(time(mod2),mean(mod2),col=2,lwd=2)
## lines(time(mod2),bconfint(mod2,level=0.95)[,1],col=4,lwd=2)
## lines(time(mod2),bconfint(mod2,level=0.95)[,2],col=4,lwd=2)
## legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),
##        inset = .01,col=c(2,4),lwd=2,cex=0.8)


###################################################
### code chunk number 10: SDEs.Rnw:254-260
###################################################
plot(mod1,ylim=c(0,40),plot.type="single")
lines(time(mod1),mean(mod1),col=2,lwd=2)
lines(time(mod1),bconfint(mod1,level=0.95)[,1],col=4,lwd=2)
lines(time(mod1),bconfint(mod1,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)


###################################################
### code chunk number 11: SDEs.Rnw:262-268
###################################################
plot(mod2,ylim=c(0,40),plot.type="single")
lines(time(mod2),mean(mod2),col=2,lwd=2)
lines(time(mod2),bconfint(mod2,level=0.95)[,1],col=4,lwd=2)
lines(time(mod2),bconfint(mod2,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)


###################################################
### code chunk number 12: SDEs.Rnw:280-286
###################################################
a = 0.5; b=1; sigma=0.1
fx <- expression( -( a*sin(x)+2*b*sin(2*x) ) )
gx <- expression( sigma )
mod <- snssde1d(drift=fx,diffusion=gx, x0=5, M=100, N=1000,Dt=0.002, method="rk3")
mod
summary(mod)


###################################################
### code chunk number 13: SDEs.Rnw:288-289 (eval = FALSE)
###################################################
## plot(mod,plot.type="single")


###################################################
### code chunk number 14: SDEs.Rnw:295-301
###################################################
plot(mod,plot.type="single")
lines(time(mod),mean(mod),col=2,lwd=2)
lines(time(mod),bconfint(mod,level=0.95)[,1],col=4,lwd=2)
lines(time(mod),bconfint(mod,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)


###################################################
### code chunk number 15: SDEs.Rnw:335-341
###################################################
K = 4; s = 1; sigma = 0.2
fx <- expression( ((0.5*sigma^2 *x^(s-1) - K)/ x^s) )
gx <- expression( sigma )
mod <- snssde1d(drift=fx,diffusion=gx, x0=3, M=100, N=1000)
mod
summary(mod)


###################################################
### code chunk number 16: SDEs.Rnw:343-344 (eval = FALSE)
###################################################
## plot(mod,plot.type="single")


###################################################
### code chunk number 17: SDEs.Rnw:349-350
###################################################
plot(mod,plot.type="single")


###################################################
### code chunk number 18: SDEs.Rnw:428-437
###################################################
a1  <- function(t) 2*t
a2  <- function(t) 0.5*t
b1  = b2 <- function(t) 0.1*t
fx    <- expression(a1(t)*x)
gx    <- expression(b1(t))
fy    <- expression(a2(t)*x)
gy    <- expression(b2(t))
mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,N=1000)
mod2d


###################################################
### code chunk number 19: SDEs.Rnw:440-442 (eval = FALSE)
###################################################
## plot(mod2d,plot.type="single")
## plot(mod2d)


###################################################
### code chunk number 20: SDEs.Rnw:447-448
###################################################
plot(mod2d,plot.type="single")


###################################################
### code chunk number 21: SDEs.Rnw:450-451
###################################################
plot(mod2d)


###################################################
### code chunk number 22: SDEs.Rnw:512-520
###################################################
mu = 4; sigma=0.1
fx <- expression( y )
gx <- expression( 0 )
fy <- expression( (mu*( 1-x^2 )* y - x) )
gy <- expression( 2*sigma)
mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,N=10000,Dt=0.01,
                  type="str")
mod2d


###################################################
### code chunk number 23: SDEs.Rnw:523-525 (eval = FALSE)
###################################################
## plot2d(mod2d,type="l")
## plot(mod2d)


###################################################
### code chunk number 24: SDEs.Rnw:530-531
###################################################
plot2d(mod2d,type="l")


###################################################
### code chunk number 25: SDEs.Rnw:533-534
###################################################
plot(mod2d)


###################################################
### code chunk number 26: SDEs.Rnw:542-545
###################################################
mu = .2
mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,N=10000,Dt=0.01,type="str")
plot2d(mod2d,type="l")


###################################################
### code chunk number 27: SDEs.Rnw:547-548
###################################################
plot(mod2d)


###################################################
### code chunk number 28: SDEs.Rnw:636-646
###################################################
K = 4; s = 1; sigma = 0.2
fx <- expression( (-K*x/sqrt(x^2+y^2+z^2)) )
gx <- expression(sigma)
fy <- expression( (-K*y/sqrt(x^2+y^2+z^2)) )
gy <- expression(sigma)
fz <- expression( (-K*z/sqrt(x^2+y^2+z^2)) )
gz <- expression(sigma)
mod3d <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,
                  N=10000,x0=1,y0=1,z0=1)
mod3d


###################################################
### code chunk number 29: SDEs.Rnw:651-653 (eval = FALSE)
###################################################
## plot3D(mod3d,display="persp",col="blue")  ## in space
## plot(mod3d,plot.type="signle")            ## with time


###################################################
### code chunk number 30: SDEs.Rnw:658-659
###################################################
plot3D(mod3d,display="persp",col="blue")


###################################################
### code chunk number 31: SDEs.Rnw:661-662
###################################################
plot(mod3d,plot.type="signle")


###################################################
### code chunk number 32: SDEs.Rnw:701-710
###################################################
fx <- expression((  x - x*y))
gx <- expression(0.03)
fy <- expression(( -y + x*y-y*z ))
gy <- expression(0.03)
fz <- expression(( -z+ y*z ))
gz <- expression(0.03)
mod3d <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,
                  N=10000,T=50,x0=0.5,y0=1,z0=2,type="str")
mod3d


###################################################
### code chunk number 33: SDEs.Rnw:713-715 (eval = FALSE)
###################################################
## plot3D(mod3d,"persp",col="blue")  ## in space 
## plot(mod3d,plot.type="signle")    ## with time


###################################################
### code chunk number 34: SDEs.Rnw:720-721
###################################################
plot3D(mod3d,"persp",col="blue")


###################################################
### code chunk number 35: SDEs.Rnw:723-724
###################################################
plot(mod3d,plot.type="signle")


###################################################
### code chunk number 36: SDEs.Rnw:744-754
###################################################
mu = 2; sigma=0.2
fx <- expression(mu*y)
gx <- expression(sigma*z)
fy <- expression(0)
gy <- expression(1)
fz <- expression(0)
gz <- expression(1)
modtra <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,
                   N=1000)
modtra


###################################################
### code chunk number 37: SDEs.Rnw:757-758 (eval = FALSE)
###################################################
## plot(modtra$XYZ[,1],ylab="X")


###################################################
### code chunk number 38: SDEs.Rnw:763-764
###################################################
plot(modtra$XYZ[,1],ylab="X")


