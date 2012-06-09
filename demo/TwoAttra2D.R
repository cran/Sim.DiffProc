## Sat Feb 05 16:24:36 2011
## Guidoum Arsalane (PG-PS/USTHB)

N=3000;t0=0;Dt=0.001;T=1;X1_0=2;X2_0=2;Y1_0=-0.5;Y2_0=-1;v=0.05;K=5;m=0.4;Sigma=0.3;

Sigmax <- Sigma
Sigmay <- Sigma
drifx1     <- expression( (-K*(x1-y1)) / (sqrt((x1-y1)^2+(x2-y2)^2))^(m+1) )
drifx2     <- expression( (-K*(x2-y2)) / (sqrt((x1-y1)^2+(x2-y2)^2))^(m+1) )
diffx      <- expression( Sigmax ) 
diffy      <- expression( Sigmay )

Ax1    <- function(t,x1,x2,y1,y2)  eval(drifx1)
Ax2    <- function(t,x1,x2,y1,y2)  eval(drifx2)
Sx     <- function(t,x1,x2,y1,y2)  eval(diffx)
Sy     <- function(t,x1,x2,y1,y2)  eval(diffy)

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

X1    <- numeric()
X2    <- numeric()
Y1    <- numeric()
Y2    <- numeric()
D     <- numeric()
X1[1] <- X1_0
X2[1] <- X2_0
Y1[1] <- Y1_0
Y2[1] <- Y2_0

for (i in 2:(N+1)){ 
        Y1[i] = Y1[i-1] + Sy(t[i-1],X1[i-1],X2[i-1],Y1[i-1],Y2[i-1]) * Dy1[i-1]
        Y2[i] = Y2[i-1] + Sy(t[i-1],X1[i-1],X2[i-1],Y1[i-1],Y2[i-1]) * Dy2[i-1]       
        X1[i] = X1[i-1] + Ax1(t[i-1],X1[i-1],X2[i-1],Y1[i-1],Y2[i-1])*Dt + 
                2 * Sx(t[i-1],X1[i-1],X2[i-1],Y1[i-1],Y2[i-1]) * Dx1[i-1]
        X2[i] = X2[i-1] + Ax2(t[i-1],X1[i-1],X2[i-1],Y1[i-1],Y2[i-1])*Dt + 
                2 * Sx(t[i-1],X1[i-1],X2[i-1],Y1[i-1],Y2[i-1]) * Dx2[i-1]
}
for (i in 1:(N+1)){D[i]  = sqrt((X1[i]-Y1[i])^2 + (X2[i]-Y2[i])^2)}
NN <- which(D <= v)
if(length(NN)>0) {nn <- min(NN)}
if(length(NN)> 0){D <- D[seq(1,nn,by=1)]}
if(length(NN)> 0){X1 <- X1[seq(1,nn,by=1)]}
if(length(NN)> 0){X2 <- X2[seq(1,nn,by=1)]}
if(length(NN)> 0){Y1 <- Y1[seq(1,nn,by=1)]}
if(length(NN)> 0){Y2 <- Y2[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}

plot(X1,X2,type="n",las=1,xlim=c(min(min(X1,Y1)),max(max(X1,Y1))),ylim=c(min(min(X2,Y2)),max(max(X2,Y2))),
     xlab=expression((list(X[t]^(1),Y[t]^(1)))),ylab=expression((list(X[t]^(2),Y[t]^(2)))),cex.lab=0.8)

mtext(expression("2-Dimensional Attractive Model for 2-Diffusion Processes"),line=3.55,adj=0.5,cex=0.8,col="black")
mtext(expression( bolditalic(M[mu(.)]^sigma * (V[t]^(1)) %<-% M[0]^sigma*(V[t]^(2)))),line=2.1,adj=0.5,cex=0.8,col="green4")

mtext(bquote( bolditalic(V[t]^(1)==(list(X[t]^(1),X[t]^(2))) ) ),cex=0.7,adj=0,line=2.4,col="red")
mtext(bquote( bolditalic(V[t]^(2)==(list(Y[t]^(1),Y[t]^(2))) ) ),cex=0.7,adj=0,line=1.4,col="blue")
mtext(bquote( bolditalic(D[t]==V[t]^(1) - V[t]^(2)) ),cex=0.7,adj=0,line=0.4,col="green4")
mtext(bquote(X[t[0]]^(1)==.(X1_0)),line=2,adj=0.85,cex=0.7,col="red")
mtext(bquote(X[t[0]]^(2)==.(X2_0)),line=1,adj=0.85,cex=0.7,col="red")
mtext(bquote(Y[t[0]]^(1)==.(Y1_0)),line=2,adj=1,cex=0.7,col="blue")
mtext(bquote(Y[t[0]]^(2)==.(Y2_0)),line=1,adj=1,cex=0.7,col="blue")
mtext(bquote((list(K,m,sigma,v))==(list(.(K),.(m),.(Sigmax),.(v)))),line=0.2,adj=1,cex=0.7,col="green4")

for (i in 1:length(D)){
lines(c(Y1[i],Y1[i+1]),c(Y2[i],Y2[i+1]),type="l",col="blue",lwd=2)
lines(c(X1[i],X1[i+1]),c(X2[i],X2[i+1]),type="l",col="red",lwd=2)
                 }
n <- which(D <= v)
if (length(n) > 0) {
thoV1V2 <- t[min(n)]
mtext(bquote( bolditalic( tau[group("||",D[t],"||")<=v]^(m)*(list(V[t]^(1),V[t]^(2)))==.(round(thoV1V2,4)) )  ),cex=0.9,adj=0.5,line=0.8,col="black")
points(X1[length(D)],X2[length(D)],pch=19,col="black",cex=1.1)
legend("topleft",bg="gray85",border="gray",expression(group("||",D[t],"||")<=v),col=c("black"),pch=19,cex=0.7)
}
thoV1V2 <- t[min(n)]
mtext(bquote( bolditalic( tau[group("||",D[t],"||")<=v]^(m)*(list(V[t]^(1),V[t]^(2)))==.(round(thoV1V2,4)) )  ),cex=0.9,adj=0.5,line=0.8,col="black")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)


 