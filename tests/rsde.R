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

fx <- expression(0,0)
gx <- expression(1,1)
res <- snssde2d(drift=fx,diffusion=gx,M=25)

Data2d <- rsde2d(res,at=0.5)
summary(Data2d)
