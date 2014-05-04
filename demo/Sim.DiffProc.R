############################################################################
#                               Demo 1                                     # 
#                              1-dim SDE                                   #
############################################################################        

f <- expression((2-x)/(1-t))
g <- expression(x)
res1 <- snssde1d(type="str",drift=f,diffusion=g,M=10,x0=1,N=1000)
res1
summary(res1)
moment(res1,order=c(2,3,4))[which(time(res1)==1),]
plot(res1,plot.type="single")
lines(time(res1),mean(res1),col=2,lwd=2)
lines(time(res1),bconfint(res1,level=0.95)[,1],col=4,lwd=2)
lines(time(res1),bconfint(res1,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),inset = .01,
       col=c(2,4),lwd=2,cex=0.8)

############################################################################
#                               Demo 2                                     # 
#                              2-dim SDE                                   #
############################################################################        

fx <- expression(tan(x-y))
gx <- expression(y-x)
fy <- expression(tan(y-x))
gy <- expression(x-y)

res1 <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,N=10000,type="str",x0=1,y0=-1)
res1
plot(res1)
dev.new()
plot(res1,plot.type="single")
dev.new()
plot2d(res1,type="n",main="2-dim sde") 
points2d(res1,col=rgb(0,100,0,50,maxColorValue=255), pch=16)

############################################################################
#                               Demo 3                                     # 
#                              3-dim SDE                                   #
############################################################################        

fx <- expression(y)
gx <- expression(z)
fy <- expression(0)
gy <- expression(1)
fz <- expression(0)
gz <- expression(1)

res <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,N=10000)
plot3D(res,display="persp") 
plot3D(res,display="rgl")

############################################################################
#                               Demo 4                                     # 
#                              1-dim FPT                                   #
############################################################################        

f <- expression( -3*(1+x) )
g <- expression( 0.5*x )
res <- fptsde1d(drift=f,diffusion=g,x0=1,c=0,M=100,N=1000)
summary(res)
bconfint(res,level=0.95)
moment(res,order=c(2,3,4,5))
dev.new()
plot(density(res$tau[!is.na(res$tau)]))

############################################################################
#                               Demo 5                                     # 
#                              1-dim RN's SDE                              #
############################################################################        

f <- expression( -3*(1+x) )
g <- expression( 0.5*x )
res <- rsde1d(drift=f,diffusion=g,M=100,N=1000,tau=0.5412)
summary(res)
bconfint(res,level=0.95)
moment(res,order=c(2,3,4,5))
dev.new()
plot(density(res$x))


############################################################################
#                               Demo 5                                     # 
#                            Fiting 1-dim SDE                              #
############################################################################  

true <- c(1,-11,2,1,0.5)
pmle <- eval(formals(fitsde.default)$pmle)

fx <- expression(theta[1] + theta[2]*x + theta[3]*x^2)
gx <- expression(theta[4]*x^theta[5])

fres <- lapply(1:4, function(i) fitsde(mydata1,drift=fx,diffusion=gx,
	             pmle=pmle[i],start = list(theta1=1,theta2=1,theta3=1,theta4=1,
				 theta5=1),optim.method = "L-BFGS-B"))
Coef <- data.frame(true,do.call("cbind",lapply(1:4,function(i) coef(fres[[i]]))))
names(Coef) <- c("True",pmle)
Summary <- data.frame(do.call("rbind",lapply(1:4,function(i) logLik(fres[[i]]))),
                      do.call("rbind",lapply(1:4,function(i) AIC(fres[[i]]))),
                      do.call("rbind",lapply(1:4,function(i) BIC(fres[[i]]))),
                      row.names=pmle)
names(Summary) <- c("logLik","AIC","BIC")
Coef	
Summary
