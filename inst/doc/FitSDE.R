### R code from vignette source 'FitSDE.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: packages
###################################################
library(Sim.DiffProc)


###################################################
### code chunk number 2: FitSDE.Rnw:151-155
###################################################
f <- expression( (1+2*x) )   
g <- expression( 0.5*x^0.3 ) 
sim    <- snssde1d(drift=f,diffusion=g,x0=2,M=1,N=1000,Dt=0.001)
mydata <- sim$X


###################################################
### code chunk number 3: FitSDE.Rnw:160-164
###################################################
fx <- expression( theta[1]+theta[2]*x ) ## drift coefficient of model (9)
gx <- expression( theta[3]*x^theta[4] ) ## diffusion coefficient of model (9)
fitmod <- fitsde(data=mydata,drift=fx,diffusion=gx,start = list(theta1=1,
                 theta2=1,theta3=1,theta4=1),pmle="euler")


###################################################
### code chunk number 4: FitSDE.Rnw:167-168
###################################################
coef(fitmod)


###################################################
### code chunk number 5: FitSDE.Rnw:171-172
###################################################
summary(fitmod)


###################################################
### code chunk number 6: FitSDE.Rnw:183-186
###################################################
vcov(fitmod)
AIC(fitmod)
confint(fitmod, level=0.95)


###################################################
### code chunk number 7: FitSDE.Rnw:213-217
###################################################
f <- expression( 3*(2-x) )   
g <- expression( 0.5 ) 
sim <- snssde1d(drift=f,diffusion=g,x0=5,Dt=0.01)
HWV <- sim$X


###################################################
### code chunk number 8: FitSDE.Rnw:221-226
###################################################
fx <- expression( theta[1]*(theta[2]- x) ) ## drift coefficient of model (13)
gx <- expression( theta[3] )           ## diffusion coefficient of model (13)
fitmod <- fitsde(data=HWV,drift=fx,diffusion=gx,start = list(theta1=1,theta2=1,
                 theta3=1),pmle="ozaki")
summary(fitmod)                 


###################################################
### code chunk number 9: FitSDE.Rnw:229-230
###################################################
confint(fitmod,parm=c("theta1","theta2"),level=0.95)


###################################################
### code chunk number 10: FitSDE.Rnw:253-257
###################################################
f <- expression(-2*x*t)
g <- expression(0.2*x)
sim <- snssde1d(drift=f,diffusion=g,N=1000,Dt=0.001,x0=10)
mydata <- sim$X


###################################################
### code chunk number 11: FitSDE.Rnw:260-268
###################################################
fx <- expression( theta[1]*x*t ) ## drift coefficient of model (17)
gx <- expression( theta[2]*x )   ## diffusion coefficient of model (17)
fitmod <- fitsde(data=mydata,drift=fx,diffusion=gx,start = list(theta1=1, 
                 theta2=1),pmle="shoji")
summary(fitmod)
vcov(fitmod)
logLik(fitmod)
confint(fitmod,level=0.9)


###################################################
### code chunk number 12: FitSDE.Rnw:286-290
###################################################
f <- expression(3*t*(sqrt(t)-x))
g <- expression(0.3*t)
sim <- snssde1d(drift=f,diffusion=g,M=1,N=1000,x0=2,Dt=0.001)
mydata <- sim$X


###################################################
### code chunk number 13: FitSDE.Rnw:293-300
###################################################
## drift coefficient of model (20)
fx <- expression( theta[1]*t* ( theta[2]*sqrt(t) - x ) ) 
## diffusion coefficient of model (20)
gx <- expression( theta[3]*t )       
fitmod <- fitsde(data=mydata,drift=fx,diffusion=gx,start = list(theta1=1,
                 theta2=1,theta3=1),pmle="kessler")
summary(fitmod)


###################################################
### code chunk number 14: FitSDE.Rnw:312-317
###################################################
theta1 = 5; theta2 = 1; theta3 = 0.2
f <- expression( ((0.5*theta3^2 *x^(theta2-1) - theta1)/ x^theta2) )
g <- expression( theta3 )
sim <- snssde1d(drift=f,diffusion=g,M=1,N=1000,x0=3,Dt=0.001)
mydata <- sim$X


###################################################
### code chunk number 15: FitSDE.Rnw:320-325
###################################################
fx <- expression( ((0.5*theta[3]^2 *x^(theta[2]-1) - theta[1])/ x^theta[2])  ) 
gx <- expression(theta[3])      
fitmod <- fitsde(mydata,drift=fx,diffusion=gx, start = list(theta1=1,theta2=1,
                theta3=1),pmle="euler")
coef(fitmod)


###################################################
### code chunk number 16: FitSDE.Rnw:328-332
###################################################
true <- c(theta1,theta2,theta3)   ## True parameters
bias <- true-coef(fitmod)
bias
confint(fitmod)


###################################################
### code chunk number 17: FitSDE.Rnw:352-356
###################################################
f <- expression( 2*x )
g <- expression( 0.3*x^0.5 )
sim <- snssde1d(drift=f,diffusion=g,M=1,N=1000,x0=2,Dt=0.001)
mydata <- sim$X


###################################################
### code chunk number 18: FitSDE.Rnw:359-386
###################################################
## True model
fx <- expression( theta[1]*x )
gx <- expression( theta[2]*x^theta[3] )
truemod <- fitsde(data=mydata,drift=fx,diffusion=gx,start = list(theta1=1,
                  theta2=1,theta3=1),pmle="euler")
## competing model 1
fx1 <- expression( theta[1]+theta[2]*x )
gx1 <- expression( theta[3]*x^theta[4] )
mod1 <- fitsde(data=mydata,drift=fx1,diffusion=gx1,start = list(theta1=1,
          theta2=1,theta3=1,theta4=1),pmle="euler")
## competing model 2  
fx2 <- expression( theta[1]+theta[2]*x )
gx2 <- expression( theta[3]*sqrt(x) )
mod2 <- fitsde(data=mydata,drift=fx2,diffusion=gx2,start = list(theta1=1,
          theta2=1,theta3=1),pmle="euler")
## competing model 3      
fx3 <- expression( theta[1] )
gx3 <- expression( theta[2]*x^theta[3] )
mod3 <- fitsde(data=mydata,drift=fx3,diffusion=gx3,start = list(theta1=1,
          theta2=1,theta3=1),pmle="euler")
## Computes AIC
AIC <- c(AIC(truemod),AIC(mod1),AIC(mod2),AIC(mod3))
Test <- data.frame(AIC,row.names = c("True mod","Comp mod1","Comp mod2",
                   "Comp mod3"))
Test
Bestmod <- rownames(Test)[which.min(Test[,1])] 
Bestmod                                              


###################################################
### code chunk number 19: FitSDE.Rnw:389-396
###################################################
Theta1 <- c(coef(truemod)[[1]],coef(mod1)[[1]],coef(mod2)[[1]],coef(mod3)[[1]])
Theta2 <- c(coef(truemod)[[2]],coef(mod1)[[2]],coef(mod2)[[2]],coef(mod3)[[2]])
Theta3 <- c(coef(truemod)[[3]],coef(mod1)[[3]],coef(mod2)[[3]],coef(mod3)[[3]])
Theta4 <- c(NA,coef(mod1)[[4]],NA,NA)
Parms  <- data.frame(Theta1,Theta2,Theta3,Theta4,row.names = c("True mod",
                     "Comp mod1","Comp mod2","Comp mod3"))
Parms                     


###################################################
### code chunk number 20: FitSDE.Rnw:401-404
###################################################
data(Irates)
rates <- Irates[, "r1"]
X <- window(rates, start = 1964.471, end = 1989.333)


###################################################
### code chunk number 21: FitSDE.Rnw:406-407 (eval = FALSE)
###################################################
## plot(X)


###################################################
### code chunk number 22: FitSDE.Rnw:412-413
###################################################
plot(X)


###################################################
### code chunk number 23: FitSDE.Rnw:420-434
###################################################
fx <- expression( theta[1]+theta[2]*x ) ## drift coefficient of model (9)
gx <- expression( theta[3]*x^theta[4] ) ## diffusion coefficient of model (9)
pmle <- eval(formals(fitsde.default)$pmle)
fitres <- lapply(1:4, function(i) fitsde(X,drift=fx,diffusion=gx,pmle=pmle[i],
                 start = list(theta1=1,theta2=1,theta3=1,theta4=1)))
Coef <- data.frame(do.call("cbind",lapply(1:4,function(i) coef(fitres[[i]]))))
Info <- data.frame(do.call("rbind",lapply(1:4,function(i) logLik(fitres[[i]]))),
                   do.call("rbind",lapply(1:4,function(i) AIC(fitres[[i]]))),
                   do.call("rbind",lapply(1:4,function(i) BIC(fitres[[i]]))),
                   row.names=pmle)
names(Coef) <- c(pmle)
names(Info) <- c("logLik","AIC","BIC")
Coef	
Info


###################################################
### code chunk number 24: FitSDE.Rnw:440-445
###################################################
f <- expression( (2.076-0.263*x) ) 
g <- expression( 0.130*x^1.451 )
mod <- snssde1d(drift=f,diffusion=g,x0=X[1],M=100, N=length(X),t0=1964.471,
                 T=1989.333)
mod


###################################################
### code chunk number 25: FitSDE.Rnw:447-454 (eval = FALSE)
###################################################
## plot(mod,plot.type="single",type="n",ylim=c(0,30))
## lines(X,col=4,lwd=2)
## lines(time(mod),mean(mod),col=2,lwd=2)
## lines(time(mod),bconfint(mod,level=0.95)[,1],col=5,lwd=2)
## lines(time(mod),bconfint(mod,level=0.95)[,2],col=5,lwd=2)
## legend("topleft",c("real data","mean path",paste("bound of", 95,"% confidence")),
##        inset = .01,col=c(4,2,5),lwd=2,cex=0.8)


###################################################
### code chunk number 26: FitSDE.Rnw:459-466
###################################################
plot(mod,plot.type="single",type="n",ylim=c(0,30))
lines(X,col=4,lwd=2)
lines(time(mod),mean(mod),col=2,lwd=2)
lines(time(mod),bconfint(mod,level=0.95)[,1],col=5,lwd=2)
lines(time(mod),bconfint(mod,level=0.95)[,2],col=5,lwd=2)
legend("topleft",c("real data","mean path",paste("bound of", 95,"% confidence")),
       inset = .01,col=c(4,2,5),lwd=2,cex=0.8)


