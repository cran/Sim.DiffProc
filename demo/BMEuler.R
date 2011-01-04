## 07/04/2010 18:27:55
## Guidoum Arsalane (PG-PS/USTHB)
## Brownian trajectory Approxsime by Euler scheme
## N  Size of process  
## T  Final time
## t0 Initial time
## x0 Initial value


pa <- par (ask=FALSE)

N  = 5000
T  = 10
t0 = 0
x0 = 0
a     <- expression (0)
sigma <- expression(1)
t = seq(t0,T,length=N)
dt = (T-t0)/N
u = runif(N-1,0,1)
o = rep(1,N-1)
o[ which( u < 0.5) ]= -1
w = c(x0,cumsum(o)*sqrt(dt))
plot(t,w,las=1,lwd=3,type="l",xlab="time",ylab=expression(W[t]))
mtext("Brownian Trajectory",line=2,cex=1.2)
mtext("Brownian trajectory Approxsime by Euler scheme",adj=0,col="red",line=0.25,cex=0.8)
mtext(bquote(x[.(t0)]==.(x0)),line=0.1,adj=1,cex=1,col="red")
mtext(bquote(t[.(x0)]==.(t0)),line=0.9,adj=1,cex=1,col="red")
legend("topleft",bg="gray85",border="gray",c("Euler Scheme"),lty=c(1),col=c("red"),lwd=2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 07/04/2010 18:27:55"),
      side = 1, line = 4, adj = 0.5, cex = .66)
D    <- diff(w)
X    <- numeric(N)
X[1] <- x0
A   <- function(t,x)  eval(a)
S   <- function(t,x)  eval(sigma)
op = par(pch = 19)
for (i in 2:N)    {X[i] = X[i-1] + A(t[i-1],X[i-1])*dt + S(t[i-1],X[i-1])*D[i-1]}
for (i in 1:(N-1)){lines(c(t[i],t[i+1]),c(X[i],X[i+1]),type="l",col="red",panel.frist=grid(col="gray"))}

par(pa)