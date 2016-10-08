library(Sim.DiffProc)

## 1-dim SDE
## Example 1: Stratonovich sde
## dX(t) = (-2*(X(t)<=0)+2*(X(t)>=0)) *dt + 0.5 o dW(t)
set.seed(1234)

f <- expression(-2*(x<=0)+2*(x>=0))
g <- expression(0.5)
res1 <- snssde1d(drift=f,diffusion=g,M=25,type="str",T=10)
x <- rsde1d(res1,at=10)
summary(x)
plot(density(x))

## 2-dim SDE's
## Example 2:
## random numbers of two standard Brownian motion W1(t) and W2(t) at time = 0.5
set.seed(1234)

fx <- expression(0)
gx <- expression(1)
fy <- expression(0)
gy <- expression(1)
res <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,M=25)

Data2d <- rsde2d(res,at=0.5)
summary(Data2d)
## library(sm)
## sm.density(Data2d,display="persp")

## 3-dim SDE's

## Example 3: Ito sde 3-dim
## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 * dW1(t) 
## dY(t) = 4*(1-Y(t)) *X(t) dt + 0.2 * dW2(t) 
## dZ(t) = 4*(1-Z(t)) *Y(t) dt + 0.2 * dW3(t)       
## W1(t), W2(t) and W3(t) three independent Brownian motion  
set.seed(1234)

fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)
fz <- expression(4*(1-z)*y)
gz <- expression(0.2)

res <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,
               x0=2,y0=-2,z0=0,M=25)
Data3d <- rsde3d(res,at=0.75)
summary(Data3d)
## library(sm)
## sm.density(Data3d,display="rgl")

