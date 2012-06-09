## Sat Dec 25 16:15:16 2010
## Guidoum Arsalane (PG-PS/USTHB)
## Brownian trajectory in 2D (O,X,Y)
## N  Size of process  
## T  Final time
## t0 Initial time
## x0 Initial value

N=10000;t0=0;x0=0;T=1;
t = seq(t0,T,length=N+1)
dt = (T-t0)/N
u = runif(N,0,1)
o = rep(1,N)
o[ which( u < 0.5) ]= -1
w1 = c(x0,cumsum(o)*sqrt(dt))
u1 = runif(N,0,1)
o1 = rep(1,N)
o1[ which( u1 < 0.5) ]= -1
w2 = c(x0,cumsum(o1)*sqrt(dt))
plot(w1,w2,las=1,lwd=3,type="l",xlab=expression(W[t]),ylab=expression(W[t]))
mtext("Brownian motion in 2D plane",line=2,cex=1.2) 
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)


