### R code from vignette source 'SDEs.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: SDEs.Rnw:126-127 (eval = FALSE)
###################################################
## install.packages("Sim.DiffProc")


###################################################
### code chunk number 2: packages
###################################################
library(Sim.DiffProc)


###################################################
### code chunk number 3: SDEs.Rnw:135-136 (eval = FALSE)
###################################################
## library(help = "Sim.DiffProc")


###################################################
### code chunk number 4: SDEs.Rnw:183-189
###################################################
f <- expression( (0.5*0.5^2*x) )
g <- expression( 0.5*x )
mod1 <- snssde1d(drift=f,diffusion=g,x0=10,M=1,N=1000)
mod2 <- snssde1d(drift=f,diffusion=g,x0=10,M=1,N=1000,type="str")
mod1
mod2


###################################################
### code chunk number 5: SDEs.Rnw:192-194 (eval = FALSE)
###################################################
## plot(mod1)
## plot(mod2)


###################################################
### code chunk number 6: SDEs.Rnw:199-200
###################################################
plot(mod1)


###################################################
### code chunk number 7: SDEs.Rnw:202-203
###################################################
plot(mod2)


###################################################
### code chunk number 8: SDEs.Rnw:220-224
###################################################
mod1 <- snssde1d(drift=f,diffusion=g,x0=10,M=50,N=1000)
mod2 <- snssde1d(drift=f,diffusion=g,x0=10,M=50,N=1000,type="str")
summary(mod1,at=1)
summary(mod2,at=1)


###################################################
### code chunk number 9: SDEs.Rnw:229-244 (eval = FALSE)
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
### code chunk number 10: SDEs.Rnw:249-255
###################################################
plot(mod1,ylim=c(0,40),plot.type="single")
lines(time(mod1),mean(mod1),col=2,lwd=2)
lines(time(mod1),bconfint(mod1,level=0.95)[,1],col=4,lwd=2)
lines(time(mod1),bconfint(mod1,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)


###################################################
### code chunk number 11: SDEs.Rnw:257-263
###################################################
plot(mod2,ylim=c(0,40),plot.type="single")
lines(time(mod2),mean(mod2),col=2,lwd=2)
lines(time(mod2),bconfint(mod2,level=0.95)[,1],col=4,lwd=2)
lines(time(mod2),bconfint(mod2,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)


###################################################
### code chunk number 12: SDEs.Rnw:297-303
###################################################
K = 4; s = 1; sigma = 0.2
fx <- expression( ((0.5*sigma^2 *x^(s-1) - K)/ x^s) )
gx <- expression( sigma )
mod <- snssde1d(drift=fx,diffusion=gx, x0=3, M=50, N=1000)
mod
summary(mod,at=0.5)


###################################################
### code chunk number 13: SDEs.Rnw:305-306 (eval = FALSE)
###################################################
## plot(mod,plot.type="single")


###################################################
### code chunk number 14: SDEs.Rnw:311-312
###################################################
plot(mod,plot.type="single")


###################################################
### code chunk number 15: SDEs.Rnw:387-395
###################################################
fx  <- expression(4*(-1-x)*y)
gx  <- expression(0.2)
fy  <- expression(4*(1-y)*x)
gy  <- expression(0.2)
mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,x0=1,y0=-1,M=50,
                  Dt=0.001,method="rk3")
mod2d
summary(mod2d)


###################################################
### code chunk number 16: SDEs.Rnw:398-400 (eval = FALSE)
###################################################
## plot(mod2d,pos=2)
## plot2d(mod2d)


###################################################
### code chunk number 17: SDEs.Rnw:405-406
###################################################
plot(mod2d,pos=2)


###################################################
### code chunk number 18: SDEs.Rnw:408-409
###################################################
plot2d(mod2d)


###################################################
### code chunk number 19: SDEs.Rnw:427-436
###################################################
a1  <- function(t) 2*t
a2  <- function(t) 0.5*t
b1  = b2 <- function(t) 0.1*t
fx    <- expression(a1(t)*x)
gx    <- expression(b1(t))
fy    <- expression(a2(t)*x)
gy    <- expression(b2(t))
mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy)
mod2d


###################################################
### code chunk number 20: SDEs.Rnw:439-440 (eval = FALSE)
###################################################
## plot(mod2d,union=TRUE,pos=3)


###################################################
### code chunk number 21: SDEs.Rnw:445-446
###################################################
plot(mod2d,union=TRUE,pos=3)


###################################################
### code chunk number 22: SDEs.Rnw:477-484
###################################################
mu = 4; sigma=0.1
fx <- expression( y )
gx <- expression( 0 )
fy <- expression( (mu*( 1-x^2 )* y - x) )
gy <- expression( 2*sigma)
mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,type="str",T=100,N=10000)
mod2d


###################################################
### code chunk number 23: SDEs.Rnw:487-489 (eval = FALSE)
###################################################
## plot2d(mod2d)
## plot(mod2d,pos=3)


###################################################
### code chunk number 24: SDEs.Rnw:494-495
###################################################
plot2d(mod2d)


###################################################
### code chunk number 25: SDEs.Rnw:497-498
###################################################
plot(mod2d,pos=3)


###################################################
### code chunk number 26: SDEs.Rnw:506-509
###################################################
mu = .2
mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,type="str",T=100,,N=10000)
plot2d(mod2d)


###################################################
### code chunk number 27: SDEs.Rnw:511-512
###################################################
plot(mod2d,pos=3)


###################################################
### code chunk number 28: SDEs.Rnw:599-609
###################################################
fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)
fz <- expression(4*(1-z)*y)
gz <- expression(0.2)
mod3d <- snssde3d(x0=2,y0=-2,z0=-2,driftx=fx,diffx=gx,drifty=fy,diffy=gy,
                driftz=fz,diffz=gz,N=1000,M=50,method="rk2")
mod3d
summary(mod3d)


###################################################
### code chunk number 29: SDEs.Rnw:613-615 (eval = FALSE)
###################################################
## plot(mod3d,union = TRUE,pos=2)   ## with time
## plot3D(mod3d,display="persp")    ## in space (O,X,Y,Z)


###################################################
### code chunk number 30: SDEs.Rnw:620-621
###################################################
plot(mod3d,union = TRUE,pos=2)


###################################################
### code chunk number 31: SDEs.Rnw:623-624
###################################################
plot3D(mod3d,display="persp")


###################################################
### code chunk number 32: SDEs.Rnw:643-653
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
### code chunk number 33: SDEs.Rnw:658-660 (eval = FALSE)
###################################################
## plot3D(mod3d,display="persp",col="blue")  ## in space
## plot(mod3d,union=TRUE,pos=2)              ## with time


###################################################
### code chunk number 34: SDEs.Rnw:665-666
###################################################
plot3D(mod3d,display="persp",col="blue")


###################################################
### code chunk number 35: SDEs.Rnw:668-669
###################################################
plot(mod3d,union=TRUE,pos=2)


###################################################
### code chunk number 36: SDEs.Rnw:708-717
###################################################
fx <- expression((  x - x*y))
gx <- expression(0.03)
fy <- expression(( -y + x*y-y*z ))
gy <- expression(0.03)
fz <- expression(( -z+ y*z ))
gz <- expression(0.03)
mod3d <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,
                  N=10000,T=20,x0=0.5,y0=1,z0=2,type="str")
mod3d


###################################################
### code chunk number 37: SDEs.Rnw:720-722 (eval = FALSE)
###################################################
## plot3D(mod3d,"persp",col="blue")  ## in space 
## plot(mod3d,union=TRUE)            ## with time


###################################################
### code chunk number 38: SDEs.Rnw:727-728
###################################################
plot3D(mod3d,"persp",col="blue")


###################################################
### code chunk number 39: SDEs.Rnw:730-731
###################################################
plot(mod3d,union=TRUE) 


###################################################
### code chunk number 40: SDEs.Rnw:751-759
###################################################
fx <- expression(2*y)
gx <- expression(0.2*z)
fy <- expression(0)
gy <- expression(1)
fz <- expression(0)
gz <- expression(1)
modtra <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz)
modtra


###################################################
### code chunk number 41: SDEs.Rnw:762-763 (eval = FALSE)
###################################################
## plot(modtra$X,plot.type="single",ylab="X")


###################################################
### code chunk number 42: SDEs.Rnw:768-769
###################################################
plot(modtra$X,plot.type="single",ylab="X")


