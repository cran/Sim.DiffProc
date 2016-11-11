library(Sim.DiffProc)

## Itô sde 2-dim
set.seed(1234)

fx <- expression(4*(-1-x)*y , 4*(1-y)*x)
gx <- expression(0.2 , 0.2)
res <- snssde2d(drift=fx,diffusion=gx,x0=c(1,-1),M=50)
res
summary(res,at=0.02548)
plot(res)
dev.new()
plot2d(res) ## in plane (O,X,Y)

## Stratonovich sde 2-dim
set.seed(1234)

fx <- expression( y , (4*( 1-x^2 )* y - x))
gx <- expression( 0 , 0.2)

res1 <- snssde2d(drift=fx,diffusion=gx,type="str",T=100,N=10000)
res1
plot(res1,pos=2)
dev.new()
plot(res1,union = FALSE)
dev.new()
plot2d(res1,type="n") ## in plane (O,X,Y)
points2d(res1,col=rgb(0,100,0,50,maxColorValue=255), pch=16)


