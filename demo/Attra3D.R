## Sat Feb 05 16:24:36 2011
## Guidoum Arsalane (PG-PS/USTHB)

N=1000;t0=0;Dt=0.001;T=1;X0=2;Y0=1;Z0=1;v=0.4;K=3.5;Sigma=0.2;

drifx     <- expression( (-K*x) / (x^2 + y^2 + z^2) )
drify     <- expression( (-K*y) / (x^2 + y^2 + z^2) )
drifz     <- expression( (-K*z) / (x^2 + y^2 + z^2) )
diff      <- expression( Sigma ) 

Ax    <- function(t,x,y,z)  eval(drifx)
Ay    <- function(t,x,y,z)  eval(drify)
Az    <- function(t,x,y,z)  eval(drifz)
S     <- function(t,x,y,z)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt= (T-t0)/N 

ux = runif(N,0,1)
ox = rep(1,N)
ox [ which(ux < 0.5) ] = -1
wx = cumsum(c(0,ox))*sqrt((T-t0)/N)
Dx    <- diff(wx)

uy = runif(N,0,1)
oy = rep(1,N)
oy [ which(uy < 0.5) ] = -1
wy = cumsum(c(0,oy))*sqrt((T-t0)/N)
Dy    <- diff(wy)

uz = runif(N,0,1)
oz = rep(1,N)
oz [ which(uz < 0.5) ] = -1
wz = cumsum(c(0,oz))*sqrt((T-t0)/N)
Dz    <- diff(wz)

X    <- numeric()
Y    <- numeric()
Z    <- numeric()
X[1] <- X0
Y[1] <- Y0
Z[1] <- Z0
for (i in 2:(N+1)){     
        X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + S(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dx[i-1]
        Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + S(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dy[i-1]
        Z[i] = Z[i-1] + Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + S(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dz[i-1]
                   } 
R = sqrt((X^2) + (Y^2) + (Z^2))
n = which(R <= v)

if (length(n) > 0){
X    <- X[c(seq(1,min(n),by=1))]
Y    <- Y[c(seq(1,min(n),by=1))]
Z    <- Z[c(seq(1,min(n),by=1))]
time <- t[c(seq(1,min(n),by=1))]}else{X <- X ; Y <- Y ; Z <- Z ; time <- t}

G <- data.frame(X,Y,Z)
V = 1
if ( V  > v ) { V =1 }
if ( V <= v ) { V = v + 1 }

a <- c(0,V,0,0)
b <- c(0,0,V,0)
c <- c(0,0,0,V)
labels <- c("Origin", "X", "Y", "Z")
i <- c(1,2,1,3,1,4)

open3d()
spheres3d(0,0,0,v,color=("white"), shininess = 128,alpha=0.2,front= "line")  
segments3d(c(0,v*cos(-pi/4)*cos(3*pi/4)),c(0,v*sin(-3*pi/4)*cos(pi/4)),c(0,v*sin(pi/4)),color = c("black"),lwd= 2.0)
text3d(0.5*v*cos(-pi/4)*cos(3*pi/4),0.5*v*sin(-3*pi/4)*cos(pi/4),0.5*v*sin(pi/4),c("v"),adj=c(0.5,-0.25),cex=1.2,family=c("serif"))

segments3d(a[i],b[i],c[i],color = c("black"),lwd= 2.0)
text3d(a,b,c,labels,adj=0.5,col="red",cex=1.2,family=c("serif"))
text3d(X0,Y0,Z0,c("(X0,Y0,Z0)"),adj=c(0.5,-0.25),cex=1.2,family=c("serif"))
points3d(G[1,],color = c("blue"),size=6)
title3d(family=c("serif"),main="Simulation Three-Dimensional Attractive Model M(S=1,Sigma)",color = c("black"),cex=1.2)

if (length(n) > 0){
for (i in 1:min(n)) {lines3d(c(G[i,1],G[i+1,1]),c(G[i,2],G[i+1,2]),c(G[i,3],G[i+1,3]),col="red",from ="lines",lwd=2)}}else
{for (i in 1:N) {lines3d(c(G[i,1],G[i+1,1]),c(G[i,2],G[i+1,2]),c(G[i,3],G[i+1,3]),col="red",from ="lines",lwd=2)}}
if (length(n) > 0 ){points3d(X[min(n)],Y[min(n)],Z[min(n)],col="blue",size=8)
                    text3d(X[min(n)],Y[min(n)],Z[min(n)],texts=c("FPT = "),adj=c(0.5,-0.8),color = c("blue"),cex=1.2,family=c("serif"))
                    text3d(X[min(n)],Y[min(n)],Z[min(n)],texts=c(t[min(n)]),adj=c(-0.9,-0.8),color = c("blue"),cex=1.2,family=c("serif"))
                    text3d(V,V,V,c("FPT : First Passage Time"),adj=c(0.5,-0.25),cex=1.2,col="blue",family=c("serif"))
                    }
title3d(family=c("serif"),font=4,sub='  Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)






