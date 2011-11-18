## Mon Jan 31 11:03:46 2011
## BOUKHETALA Kamal , GUIDOUM Arsalane.
## USTHB, Maths-PS
## Simulation SDE in 2D (O,X,Y) 


N=5000;T=1;t0=0;x0=0;y0=0;Dt=0.001;
driftx <- expression(cos(t*x*y))
drifty <- expression(cos(t))
diffx  <- expression(0.1)
diffy  <- expression(0.1)

Ax    <- function(t,x,y)  eval(driftx)
Ay    <- function(t,x,y)  eval(drifty)
Sx    <- function(t,x,y)  eval(diffx)
Sy    <- function(t,x,y)  eval(diffy)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
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
X    <- numeric()
Y    <- numeric()
X[1] <- x0
Y[1] <- y0
for (i in 2:(N+1)){
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1] 
                  } 
plot(X,Y,type="n",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="black",lwd=1)}
mtext(expression("Euler scheme : Simulation SDE Two-Dimensional"),line=3.4,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",bg="gray85",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria Mon Jan 31 11:03:46 2011"),side = 1, line = 4, adj = 0.5, cex = .66)


