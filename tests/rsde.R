library(Sim.DiffProc)

## rsde 1-dim
## Itô sde

f <- expression( 2*(3-x) )
g <- expression( 1 )
res1 <- rsde1d(drift=f,diffusion=g,M=40,N=1000,tau=0.5412)
res1
summary(res1)
bconfint(res1,level=0.95)
moment(res1,order=c(2,3,4,5))
plot(density(res1$x))

## Stratonovich sde

f <- expression(-2*(x<=0)+2*(x>=0))
g <- expression(0.5)
res2 <- rsde1d(drift=f,diffusion=g,M=40,N=1000,tau=0.981232)
res2
summary(res2)
bconfint(res2,level=0.95)
moment(res2,order=c(2,3,4,5))
plot(density(res2$x))

## rsde 2-dim
## Itô sde

fx <- expression(0)
gx <- expression(1)
fy <- expression(0)
gy <- expression(1)
res1 <- rsde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,tau=1,M=40,N=1000)
res1
summary(res1)
bconfint(res1,level=0.95)
moment(res1,order=c(2,3,4,5))


## Stratonovich sde 

fx <- expression(4*(-2-x)*t)
gx <- expression(0.4*y)
fy <- expression(0*(y>0)-2*(y<=0))
gy <- expression(x)
res2 <- rsde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,tau=1,M=40,N=1000,type="str")
res2
summary(res2)
bconfint(res2,level=0.95)
moment(res2,order=c(2,3,4,5))


## rsde 3-dim
## Itô sde

fx <- expression(2*(3-x))
gx <- expression(y+z)
fy <- expression(2*(3-y))
gy <- expression(x+z)
fz <- expression(2*(3-z))
gz <- expression(x+y)

res1 <- rsde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,N=1000,M=40,Dt=0.01,tau=10)
res1
summary(res1)
bconfint(res1,level=0.95)
moment(res1,order=c(2,3,4,5))

