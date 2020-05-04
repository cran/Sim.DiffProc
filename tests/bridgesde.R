library(Sim.DiffProc)


### 2-dim Ito bridge sde
set.seed(1234)

fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)
res <- bridgesde2d(drift=fx,diffusion=gx,Dt=0.005,M=2000)
res
summary(res) ## Monte-Carlo statistics at time T/2=2.5
summary(res,at=1) ## Monte-Carlo statistics at time 1
summary(res,at=4) ## Monte-Carlo statistics at time 4
##
plot(res,type="n")
lines(time(res),apply(res$X,1,mean),col=3,lwd=2)
lines(time(res),apply(res$Y,1,mean),col=4,lwd=2)
legend("topright",c(expression(E(X[t])),expression(E(Y[t]))),lty=1,inset = .7,col=c(3,4))
##
plot2d(res)


##

fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)
Sigma= matrix(c(1,0.7,0.7,1),nrow=2)
res <- bridgesde2d(drift=fx,diffusion=gx,Dt=0.005,M=500,corr=Sigma)
res

## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 * dW1(t) ; x01 = 0 and y01 = 0
## dY(t) = 4*(1-Y(t)) *X(t) dt + 0.2 * dW2(t) ; x02 = -1 and y02 = -2
## dZ(t) = 4*(1-Z(t)) *Y(t) dt + 0.2 * dW3(t) ; x03 = 0.5 and y03 = 0.5       
## W1(t), W2(t) and W3(t) are three correlated Brownian motions with Sigma
set.seed(1234)

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
# correlation matrix
Sigma <-matrix(c(1,0.3,0.5,0.3,1,0.2,0.5,0.2,1),nrow=3,ncol=3) 

res <- bridgesde3d(x0=c(0,-1,0.5),y=c(0,-2,0.5),drift=fx,diffusion=gx,corr=Sigma,M=200)
res
summary(res) ## Monte-Carlo statistics at time T/2=0.5
summary(res,at=0.25) ## Monte-Carlo statistics at time 0.25
summary(res,at=0.75) ## Monte-Carlo statistics at time 0.75
##
plot(res,type="n")
lines(time(res),apply(res$X,1,mean),col=3,lwd=2)
lines(time(res),apply(res$Y,1,mean),col=4,lwd=2)
lines(time(res),apply(res$Z,1,mean),col=5,lwd=2)
legend("topleft",c(expression(E(X[t])),expression(E(Y[t])),
       expression(E(Z[t]))),lty=1,inset = .01,col=c(3,4,5))
##
plot3D(res,display = "persp",main="3-dim bridge sde")