## Tue May 22 16:12:16 2012
## GUIDOUM  Arsalane Chouaib
## Stochastic Lotka-Volterra Model

N=5000;t0=0;T=100;x0=1;y0=1;a=1;b=2;c=0.5;d=0.25;sigma=0.0001

diffx  <- expression(sigma)
diffy  <- expression(sigma)
driftx <- expression(a * x - b * x * y)
drifty <- expression(c * x * y - d * y)

DSxx   <- D(diffx,"x")
DSyy   <- D(diffy,"y")
Ax     <- function(t,x,y)  eval(driftx)
Ay     <- function(t,x,y)  eval(drifty)
Sx     <- function(t,x,y)  eval(diffx)
DSx    <- function(t,x,y)  eval(DSxx)
Sy     <- function(t,x,y)  eval(diffy)
DSy    <- function(t,x,y)  eval(DSyy)


Dt = (T-t0)/N
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X    <- numeric()
Y    <- numeric()
X[1] <- x0
Y[1] <- y0
for (i in 2:(N+1)){
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]+
           0.5 *Sx(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1])*((Dx[i-1])^2 -Dt)
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]+
           0.5 *Sy(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1])*((Dy[i-1])^2 -Dt) 
                  } 
plot(X,Y,type="l",xlab=expression(X[t]),ylab=expression(Y[t]),las=1)
points(x0,y0,pch=20,col="red")
mtext(expression("Stochastic Lotka-Volterra Model"),line=3,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*minute-a*x[t]+b*x[t]*y[t]==epsilon[t]),line=2.3,adj=0,cex=1,col="red")
mtext(expression(y[t]*minute-c*x[t]*y[t]+d*y[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.1,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=2,adj=0.78,cex=1,col="blue")
mtext(bquote(y[0]==.(y0)),line=1.2,adj=0.78,cex=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.3,cex=1,adj=0.78,col="blue")
mtext(bquote(a==.(a)),line=2.5,cex=1,adj=1,col="blue")
mtext(bquote(b==.(b)),line=1.8,cex=1,adj=1,col="blue")
mtext(bquote(c==.(c)),line=0.9,cex=1,adj=1,col="blue")
mtext(bquote(d==.(d)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)

