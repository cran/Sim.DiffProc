## Sat Feb 05 16:24:36 2011
## Guidoum Arsalane (PG-PS/USTHB)

N=1000;t0=0;Dt=0.001;T=1;X1_0=1;X2_0=0;X3_0=-1;Y1_0=1;
Y2_0=0;Y3_0=1;v=0.05;K=3;m=0.1;Sigma=0.2;

Sigmax <- Sigma
Sigmay <- Sigma
drifx1     <- expression( (-K*(x1-y1)) / (sqrt((x1-y1)^2+(x2-y2)^2 +(x3-y3)^2))^(m+1) )
drifx2     <- expression( (-K*(x2-y2)) / (sqrt((x1-y1)^2+(x2-y2)^2 +(x3-y3)^2))^(m+1) )
drifx3     <- expression( (-K*(x3-y3)) / (sqrt((x1-y1)^2+(x2-y2)^2 +(x3-y3)^2))^(m+1) )
diffx      <- expression( Sigmax ) 
diffy      <- expression( Sigmay )


Ax1    <- function(t,x1,x2,x3,y1,y2,y3)  eval(drifx1)
Ax2    <- function(t,x1,x2,x3,y1,y2,y3)  eval(drifx2)
Ax3    <- function(t,x1,x2,x3,y1,y2,y3)  eval(drifx3)
Sx     <- function(t,x1,x2,x3,y1,y2,y3)  eval(diffx)
Sy     <- function(t,x1,x2,x3,y1,y2,y3)  eval(diffy)


if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt= (T-t0)/N 

ux1 = runif(N,0,1)
ox1 = rep(1,N)
ox1 [ which(ux1 < 0.5) ] = -1
wx1 = cumsum(c(0,ox1))*sqrt((T-t0)/N)
Dx1     <- diff(wx1)

ux2 = runif(N,0,1)
ox2 = rep(1,N)
ox2 [ which(ux2 < 0.5) ] = -1
wx2 = cumsum(c(0,ox2))*sqrt((T-t0)/N)
Dx2     <- diff(wx2)

ux3 = runif(N,0,1)
ox3 = rep(1,N)
ox3 [ which(ux3 < 0.5) ] = -1
wx3 = cumsum(c(0,ox3))*sqrt((T-t0)/N)
Dx3     <- diff(wx3)

uy1 = runif(N,0,1)
oy1 = rep(1,N)
oy1 [ which(uy1 < 0.5) ] = -1
wy1 = cumsum(c(0,oy1))*sqrt((T-t0)/N)
Dy1     <- diff(wy1)

uy2 = runif(N,0,1)
oy2 = rep(1,N)
oy2 [ which(uy2 < 0.5) ] = -1
wy2 = cumsum(c(0,oy2))*sqrt((T-t0)/N)
Dy2     <- diff(wy2)

uy3 = runif(N,0,1)
oy3 = rep(1,N)
oy3 [ which(uy3 < 0.5) ] = -1
wy3 = cumsum(c(0,oy3))*sqrt((T-t0)/N)
Dy3     <- diff(wy3)

X1    <- numeric()
X2    <- numeric()
X3    <- numeric()
Y1    <- numeric()
Y2    <- numeric()
Y3    <- numeric()
D     <- numeric()
X1[1] <- X1_0
X2[1] <- X2_0
X3[1] <- X3_0
Y1[1] <- Y1_0
Y2[1] <- Y2_0
Y3[1] <- Y3_0
for (i in 2:(N+1)){ 
        Y1[i] = Y1[i-1] + Sy(t[i-1],X1[i-1],X2[i-1],X3[i-1],Y1[i-1],Y2[i-1],Y3[i-1]) * Dy1[i-1]
        Y2[i] = Y2[i-1] + Sy(t[i-1],X1[i-1],X2[i-1],X3[i-1],Y1[i-1],Y2[i-1],Y3[i-1]) * Dy2[i-1] 
        Y3[i] = Y3[i-1] + Sy(t[i-1],X1[i-1],X2[i-1],X3[i-1],Y1[i-1],Y2[i-1],Y3[i-1]) * Dy3[i-1]   
   
        X1[i] = X1[i-1] + Ax1(t[i-1],X1[i-1],X2[i-1],X3[i-1],Y1[i-1],Y2[i-1],Y3[i-1])*Dt + 
                2* Sx(t[i-1],X1[i-1],X2[i-1],X3[i-1],Y1[i-1],Y2[i-1],Y3[i-1]) * Dx1[i-1]
        X2[i] = X2[i-1] + Ax2(t[i-1],X1[i-1],X2[i-1],X3[i-1],Y1[i-1],Y2[i-1],Y3[i-1])*Dt + 
                2* Sx(t[i-1],X1[i-1],X2[i-1],X3[i-1],Y1[i-1],Y2[i-1],Y3[i-1]) * Dx2[i-1]
        X3[i] = X3[i-1] + Ax3(t[i-1],X1[i-1],X2[i-1],X3[i-1],Y1[i-1],Y2[i-1],Y3[i-1])*Dt + 
                2* Sx(t[i-1],X1[i-1],X2[i-1],X3[i-1],Y1[i-1],Y2[i-1],Y3[i-1]) * Dx3[i-1]
}
D  = sqrt((X1-Y1)^2 + (X2-Y2)^2 +(X3-Y3)^2)
n = which(D <= v)
if (length(n) > 0){
X1    <- X1[c(seq(1,min(n),by=1))]
X2    <- X2[c(seq(1,min(n),by=1))]
X3    <- X3[c(seq(1,min(n),by=1))]
Y1    <- Y1[c(seq(1,min(n),by=1))]
Y2    <- Y2[c(seq(1,min(n),by=1))]
Y3    <- Y3[c(seq(1,min(n),by=1))]
D     <- D[c(seq(1,min(n),by=1))]
t     <- t[c(seq(1,min(n),by=1))]}else{X1 <- X1 ; X2 <- X2 ; X3 <- X3;Y1 <- Y1;Y2 <- Y2;Y3<-Y3;D<-D ; t <- t}

Gx <- data.frame(X1,X2,X3)
Gy <- data.frame(Y1,Y2,Y3)
V = 1
if ( V  > v ) { V =1 }
if ( V <= v ) { V = v + 1 }

a <- c(0,V,0,0)
b <- c(0,0,V,0)
c <- c(0,0,0,V)
labels <- c("Origin", "X", "Y", "Z")
i <- c(1,2,1,3,1,4)

open3d()
##spheres3d(0,0,0,v,color=("white"), shininess = 128,alpha=0.2,front= "line")  
##segments3d(c(0,v*cos(-pi/4)*cos(3*pi/4)),c(0,v*sin(-3*pi/4)*cos(pi/4)),c(0,v*sin(pi/4)),color = c("black"),lwd= 2.0)
##text3d(0.5*v*cos(-pi/4)*cos(3*pi/4),0.5*v*sin(-3*pi/4)*cos(pi/4),0.5*v*sin(pi/4),c("v"),adj=c(0.5,-0.25),cex=1.2,family=c("serif"))

segments3d(a[i],b[i],c[i],color = c("black"),lwd= 2.0)
text3d(a,b,c,labels,adj=0.5,col="green4",cex=1.2,family=c("serif"))
text3d(X1_0,X2_0,X3_0,c("(X1,X2,X3)"),adj=c(0.5,-0.25),cex=1.2,family=c("serif"))
points3d(Gx[1,],color = c("red"),size=6)
text3d(Y1_0,Y2_0,Y3_0,c("(Y1,Y2,Y3)"),adj=c(0.5,-0.25),cex=1.2,family=c("serif"))
points3d(Gy[1,],color = c("blue"),size=6)
title3d(family=c("serif"),main="3-Dimensional Attractive Model for 2-Diffusion Processes",color = c("black"),cex=1.2)

if (length(n) > 0){
for (i in 1:min(n)) {lines3d(c(Gy[i,1],Gy[i+1,1]),c(Gy[i,2],Gy[i+1,2]),c(Gy[i,3],Gy[i+1,3]),col="blue",from ="lines",lwd=2)
                     lines3d(c(Gx[i,1],Gx[i+1,1]),c(Gx[i,2],Gx[i+1,2]),c(Gx[i,3],Gx[i+1,3]),col="red",from ="lines",lwd=2)}}else
{for (i in 1:N) {lines3d(c(Gy[i,1],Gy[i+1,1]),c(Gy[i,2],Gy[i+1,2]),c(Gy[i,3],Gy[i+1,3]),col="blue",from ="lines",lwd=2)
                 lines3d(c(Gx[i,1],Gx[i+1,1]),c(Gx[i,2],Gx[i+1,2]),c(Gx[i,3],Gx[i+1,3]),col="red",from ="lines",lwd=2)}}

if (length(n) > 0 ){points3d(X1[min(n)],X2[min(n)],X3[min(n)],col="green",size=8)
                    text3d(X1[min(n)],X2[min(n)],X3[min(n)],texts=c("FPT = "),adj=c(0.5,-0.8),color = c("green"),cex=1.2,family=c("serif"))
                    text3d(X1[min(n)],X2[min(n)],X3[min(n)],texts=c(t[min(n)]),adj=c(-0.9,-0.8),color = c("green"),cex=1.2,family=c("serif"))
                    text3d(V,V,V,c("FPT : First Passage Time"),adj=c(0.5,-0.25),cex=1.2,col="green",family=c("serif"))
                    }
title3d(family=c("serif"),font=4,sub='USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria Sat Feb 05 16:34:11 2011',color = c("blue"),cex=0.8)

