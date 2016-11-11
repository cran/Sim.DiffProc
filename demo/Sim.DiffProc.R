############################################################################
#                               Demo 1                                     # 
#                              1-dim SDE                                   #
############################################################################        
set.seed(1234)
theta = 0.5
f <- expression( (0.5*theta^2*x) )
g <- expression( theta*x )
mod1 <- snssde1d(drift=f,diffusion=g,x0=10,M=500,type="ito")  
mod1
summary(mod1, at = 1)
x1 <- rsde1d(object = mod1, at = 1) 
summary(x1)
mu1 = log(10); sigma1= sqrt(theta^2) 
AppdensI <- dsde1d(mod1, at = 1)
plot(AppdensI , dens = function(x) dlnorm(x,meanlog=mu1,sdlog = sigma1))
dev.new()
plot(mod1,plot.type="single")
lines(time(mod1),mean(mod1),col=2,lwd=2)
lines(time(mod1),bconfint(mod1,level=0.95)[,1],col=4,lwd=2)
lines(time(mod1),bconfint(mod1,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of",95,"% confidence")),col=c(2,4),lwd=2)

############################################################################
#                               Demo 2                                     # 
#                              2-dim SDE                                   #
############################################################################        
set.seed(1234)
fx <- expression( y ,(4*( 1-x^2 )* y - x))
gx <- expression( 0 , 0.2)

res1 <- snssde2d(drift=fx,diffusion=gx,type="str",T=100,N=10000)
res1
plot(res1,pos=2)
dev.new()
plot(res1,union = FALSE)
dev.new()
plot2d(res1,type="n") ## in plane (O,X,Y)
points2d(res1,col=rgb(0,100,0,50,maxColorValue=255), pch=16)

############################################################################
#                               Demo 3                                     # 
#                              3-dim SDE                                   #
############################################################################        
set.seed(1234)

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
res <- snssde3d(x0=c(2,-2,-2),drift=fx,diffusion=gx,M=100)
res
plot(res,pos=2)
dev.new()
plot3D(res,display="persp")

############################################################################
#                               Demo 4                                     # 
#                          2-dim Bridge SDE                                #
############################################################################
set.seed(1234)
fx <- expression(4*(-1-x)*y, 4*(1-y)*x)
gx <- expression(0.2, 0.2)

res <- bridgesde2d(x0=c(0,-1),y=c(1,0),drift=fx,diffusion=gx,M=50)
res
plot(res)
dev.new()
plot2d(res,type="n")
points2d(res,col=rgb(0,100,0,50,maxColorValue=255), pch=16)

############################################################################
#                               Demo 5                                     # 
#                            Bridge 3-dim SDE                              #
############################################################################  
set.seed(1234)
fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)

res <- bridgesde3d(x0=c(0,-1,0.5),y=c(0,-2,0.5),drift=fx,diffusion=gx,M=20)
res
plot(res,union=TRUE)
dev.new()
plot3D(res,display = "persp",main="3-dim bridge sde")

############################################################################
#                               Demo 6                                     # 
#                              1-dim FPT                                   #
############################################################################        

## X(t) Brownian motion
## S(t) = 0.3+0.2*t (time-dependent boundary)
set.seed(1234)

f <- expression( 0 )
g <- expression( 1 )
St <- expression(0.5-0.5*t) 
mod1 <- snssde1d(drift=f,diffusion=g,M=100)
fptmod1 <- rfptsde1d(mod1,boundary=St)
summary(fptmod1)
plot(dfptsde1d(mod1,boundary=St))



############################################################################
#                               Demo 7                                     # 
#                            Fiting 1-dim SDE                              #
############################################################################  
set.seed(1234)

## Application to real data
## CKLS modele vs CIR modele 
## CKLS (mod1):  dX(t) = (theta1+theta2* X(t))* dt + theta3 * X(t)^theta4 * dW(t)
## CIR  (mod2):  dX(t) = (theta1+theta2* X(t))* dt + theta3 * sqrt(X(t))  * dW(t)
set.seed(1234)

data(Irates)
rates <- Irates[,"r1"]
rates <- window(rates, start=1964.471, end=1989.333)

fx1 <- expression(theta[1]+theta[2]*x)
gx1 <- expression(theta[3]*x^theta[4])
gx2 <- expression(theta[3]*sqrt(x))

fitmod1 <- fitsde(rates,drift=fx1,diffusion=gx1,pmle="euler",start = list(theta1=1,theta2=1,
                  theta3=1,theta4=1),optim.method = "L-BFGS-B")
fitmod2 <- fitsde(rates,drift=fx1,diffusion=gx2,pmle="euler",start = list(theta1=1,theta2=1,
                  theta3=1),optim.method = "L-BFGS-B")	
summary(fitmod1)
summary(fitmod2)
coef(fitmod1)
coef(fitmod2)
confint(fitmod1,parm=c('theta2','theta3'))
confint(fitmod2,parm=c('theta2','theta3'))
AIC(fitmod1)
AIC(fitmod2)	

############################################################################
#                               Demo 8                                     # 
#                 Transition Density of Brownian motion                    #
############################################################################  

f <- expression(0); g <- expression(1)
B <- snssde1d(drift=f,diffusion=g,M=1000)
for (i in seq(B$t0,B$T,by=B$Dt)){
plot(dsde1d(B, at = i),main=paste0('Transition Density \n t = ',i))
}

############################################################################
#                               Demo 9                                     # 
#     Bivariate Transition Density of 2 Brownian motion                    #
############################################################################  

fx <- expression(0,0)
gx <- expression(1,1)
B <- snssde2d(drift=fx,diffusion=gx,M=1000)
B
for (i in seq(B$Dt,B$T,by=B$Dt)){
plot(dsde2d(B, at = i),display="contour",main=paste0('Transition Density \n t = ',i))
}
