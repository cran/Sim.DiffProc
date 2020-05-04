library(Sim.DiffProc)

## sde 2-dim
set.seed(1234)

fx <- expression( y , (4*( 1-x^2 )* y - x))
gx <- expression( 0 , 0.2)

mod2d2 <- snssde2d(drift=fx,diffusion=gx,type="str",T=100,N=10000)
mod2d2
plot(mod2d2,pos=2)
dev.new()
plot(mod2d2,union = FALSE)
dev.new()
plot2d(mod2d2,type="n") ## in plane (O,X,Y)
points2d(mod2d2,col=rgb(0,100,0,50,maxColorValue=255), pch=16)


## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 dW1(t)
## dY(t) = 4*(1-Y(t))*X(t) dt + 0.2 dW2(t)
set.seed(1234)

fx <- expression(4*(-1-x)*y , 4*(1-y)*x )
gx <- expression(0.25*y,0.2*x)

Sigma <- matrix(c(0.9,0.3,0.3,0.4),nrow=2,ncol=2) # correlation matrix
mod2d1 <- snssde2d(drift=fx,diffusion=gx,corr=Sigma,x0=c(x0=1,y0=-1),M=1000)
mod2d1
summary(mod2d1)
##
dev.new()
plot(mod2d1,type="n")
mx <- apply(mod2d1$X,1,mean)
my <- apply(mod2d1$Y,1,mean)
lines(time(mod2d1),mx,col=1)
lines(time(mod2d1),my,col=2)
legend("topright",c(expression(E(X[t])),expression(E(Y[t]))),lty=1,inset = .01,col=c(1,2),cex=0.95)
##
dev.new()
plot2d(mod2d1) ## in plane (O,X,Y)
lines(my~mx,col=2)

## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 dW1(t)
## dY(t) = 4*(1-Y(t))*X(t) dt + 0.2 dW2(t)
set.seed(1234)

fx <- expression(4*(-1-x)*y , 4*(1-y)*x )
gx <- expression(0.25*y,0.2*x)

mod2d1 <- snssde2d(drift=fx,diffusion=gx,x0=c(x0=1,y0=-1),M=1000)
mod2d1
summary(mod2d1)
##
dev.new()
plot(mod2d1,type="n")
mx <- apply(mod2d1$X,1,mean)
my <- apply(mod2d1$Y,1,mean)
lines(time(mod2d1),mx,col=1)
lines(time(mod2d1),my,col=2)
legend("topright",c(expression(E(X[t])),expression(E(Y[t]))),lty=1,inset = .01,col=c(1,2),cex=0.95)
##
dev.new()
plot2d(mod2d1) ## in plane (O,X,Y)
lines(my~mx,col=2)

## Now W1(t) and W2(t) are correlated.

set.seed(1234)
Sigma <- matrix(c(0.9,0.3,0.3,0.4),nrow=2,ncol=2) # correlation matrix
mod2d1 <- snssde2d(drift=fx,diffusion=gx,corr=Sigma,x0=c(x0=1,y0=-1),M=1000)
mod2d1
summary(mod2d1)
##
dev.new()
plot(mod2d1,type="n")
mx <- apply(mod2d1$X,1,mean)
my <- apply(mod2d1$Y,1,mean)
lines(time(mod2d1),mx,col=1)
lines(time(mod2d1),my,col=2)
legend("topright",c(expression(E(X[t])),expression(E(Y[t]))),lty=1,inset = .01,col=c(1,2),cex=0.95)
##
dev.new()
plot2d(mod2d1) ## in plane (O,X,Y)
lines(my~mx,col=2)


## dX(t) = Y(t)* dt  + 0.2 o dW3(t)         
## dY(t) = (4*( 1-X(t)^2 )* Y(t) - X(t))* dt + 0.2 o dW2(t)
## dZ(t) = (4*( 1-X(t)^2 )* Z(t) - X(t))* dt + 0.2 o dW3(t)
## W1(t), W2(t) and W3(t) are three correlated Brownian motions with Sigma

fx <- expression( y , (4*( 1-x^2 )* y - x), (4*( 1-x^2 )* z - x))
gx <- expression( 0.2 , 0.2, 0.2)
# correlation matrix
Sigma <-matrix(c(1,0.3,0.5,0.3,1,0.2,0.5,0.2,1),nrow=3,ncol=3) 

mod3d2 <- snssde3d(drift=fx,diffusion=gx,N=10000,T=100,type="str",corr=Sigma)
mod3d2
