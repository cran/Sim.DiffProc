## Mon Apr 30 23:33:23 2012
## GUIDOUM  Arsalane Chouaib
## Stochastics Oscillators
## stochastic Van der Pol oscillator

N=10000;T=100;x0=1;v0=0;a=3;b=0.3;omega=2.5;sigma=0.1

diffx  <- expression(0)
diffy  <- expression(sigma)
driftx <- expression(y)
drifty <- expression(-(a*y*(b^-2 * x^2 -1)+omega^2 * x))
DSxx   <- D(diffx,"x")
DSyy   <- D(diffy,"y")
Ax     <- function(t,x,y)  eval(driftx)
Ay     <- function(t,x,y)  eval(drifty)
Sx     <- function(t,x,y)  eval(diffx)
DSx    <- function(t,x,y)  eval(DSxx)
Sy     <- function(t,x,y)  eval(diffy)
DSy    <- function(t,x,y)  eval(DSyy)
Dt = T/N
t <- seq(0,T,length=N+1)
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X     <- numeric()
Y     <- numeric()
X[1]  <- x0
Y[1]  <- v0
for (i in 2:(N+1)){
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]+
           0.5 *Sx(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1])*((Dx[i-1])^2 -Dt)
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]+
           0.5 *Sy(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1])*((Dy[i-1])^2 -Dt) 
                  } 
plot(X,(Y/omega),type="l",axes = FALSE,xlab=expression(x[t]),ylab=expression(x[t]*minute/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
points(x0,v0/omega,pch=20,col="red")
mtext(expression("The phase portrait of stochastic Van Der Pol oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second+a*x[t]*minute*(x[t]^2 / b - 1)+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=2.9,cex=1,adj=1,col="blue")
mtext(bquote(b==.(b)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)				  
