############################################################################
#                               Demo 1                                     # 
#                              1-dim SDE                                   #
############################################################################        
set.seed(1234)
f <- expression((2-x)/(1-t))
g <- expression(x)
res1 <- snssde1d(type="str",drift=f,diffusion=g,M=100,x0=1,N=1000)
res1
summary(res1,at=1)
moment(res1,order=c(2,3,4))[which(time(res1)==1),]
plot(res1,plot.type="single")
lines(time(res1),mean(res1),col=2,lwd=2)
lines(time(res1),bconfint(res1,level=0.95)[,1],col=4,lwd=2)
lines(time(res1),bconfint(res1,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),inset = .01,
       col=c(2,4),lwd=2,cex=0.8)
	   
rsde1d(res1,at=1)
plot(density(rsde1d(res1,at=1)))	   

############################################################################
#                               Demo 2                                     # 
#                              2-dim SDE                                   #
############################################################################        
set.seed(1234)
fx <- expression( y )
gx <- expression( 0 )
fy <- expression( (4*( 1-x^2 )* y - x) )
gy <- expression( 0.2)

res1 <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,type="str",T=100,
                 ,N=10000)
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
fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)
fz <- expression(4*(1-z)*y)
gz <- expression(0.2)

res <- snssde3d(x0=2,y0=-2,z0=-2,driftx=fx,diffx=gx,drifty=fy,diffy=gy,
                driftz=fz,diffz=gz,N=1000,M=100)
plot(res,pos=2)
dev.new()
plot3D(res,display="persp")

############################################################################
#                               Demo 4                                     # 
#                          2-dim Bridge SDE                                #
############################################################################
set.seed(1234)
fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)

res <- bridgesde2d(x0=c(0,-1),y=c(1,0),driftx=fx,diffx=gx,drifty=fy,diffy=gy,M=50)
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
fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)
fz <- expression(4*(1-z)*y)
gz <- expression(0.2)

res <- bridgesde3d(x0=c(0,-1,0.5),y=c(0,-2,0.5),driftx=fx,diffx=gx,
                drifty=fy,diffy=gy,driftz=fz,diffz=gz,M=20)
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
mod1 <- snssde1d(drift=f,diffusion=g,M=40)
fptmod1 <- fptsde1d(mod1,boundary=St)
summary(fptmod1)
plot(fptmod1)
plot(density(fptmod1$fpt[!is.na(fptmod1$fpt)]),main="Kernel Density of a First-Passage-Time")


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





