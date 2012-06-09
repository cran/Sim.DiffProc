## Mon Jan 31 11:03:46 2011
## BOUKHETALA Kamal , GUIDOUM Arsalane.
## USTHB, Maths-PS
## Simulation SDE in 3D (O,X,Y,Z) 

N = 2500
t0= 0
x0 = y0 = z0 = 0
Dt = 0.001
driftx <- expression(cos(t*x*y))
drifty <- expression(cos(t))
driftz <- expression(cos(t*x))
diffx <- expression(0.1)
diffy <- expression(0.1)
diffz <- expression(0.1)

Ax    <- function(t,x,y,z)  eval(driftx)
Ay    <- function(t,x,y,z)  eval(drifty)
Az    <- function(t,x,y,z)  eval(driftz)
Sx    <- function(t,x,y,z)  eval(diffx)
Sy    <- function(t,x,y,z)  eval(diffy)
Sz    <- function(t,x,y,z)  eval(diffz)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
wz = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dz    <- diff(wz)
X    <- numeric()
Y    <- numeric()
Z    <- numeric()
X[1] <- x0
Y[1] <- y0
Z[1] <- z0
for (i in 2:(N+1)){
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dx[i-1]
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dy[i-1] 
	Z[i] = Z[i-1] + Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dz[i-1] 
                  } 
G <- data.frame(X,Y,Z)
open3d()
a <- c(0,1,0,0)
b <- c(0,0,1,0)
c <- c(0,0,0,1)
labels <- c("O", "X", "Y", "Z")
i <- c(1,2,1,3,1,4)
segments3d(a[i],b[i],c[i],color = c("black"),lwd= 2.0,box=T)
text3d(a,b,c,labels,adj=0.5,col="red",cex=1.2,family=c("serif"))
text3d(x0,y0,z0,c("(X0,Y0,Z0)"),adj=c(0.5,-0.25),cex=0.8,family=c("serif"),col="blue")
points3d(G[1,],color = c("blue"),size=6)
for (i in 1:N) {lines3d(c(G[i,1],G[i+1,1]),c(G[i,2],G[i+1,2]),c(G[i,3],G[i+1,3]),col="black",from ="lines",lwd=2)}
title3d(family=c("serif"),main="Euler scheme : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='  Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
