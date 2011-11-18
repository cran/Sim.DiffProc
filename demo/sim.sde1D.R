## 23:12 04/04/2010
## Guidoum Arsalane (PG-PS/USTHB)
## Stochastic process
## Simulation Numerical Solution of Stochastic Differential Equation
## N Size of process  
## T Final time
## Dt Discretization  
## t0 Initial time
## x0 Initial value 
## a   = a(t,X(t))   Coefficient of drift 
## sigma=s(t,X(t))   Coefficient of diffusion
## Stochastic differential equation : dX(t) = (0.03*t*x-x^5)*dt + 0.1*dW(t)
## Euler scheme

pa <- par (ask=FALSE)

N  = 1000
t0 = 0
x0 = 0
T  = 100
Dt = (T-t0)/N
a     <- expression( (0.03*t*x-x^5) )
sigma <- expression( (0.1) )
A     <- function(t,x)  eval(a)
S     <- function(t,x)  eval(sigma)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
X    <- numeric(N+1)
X[1] <- x0
for (i in 2:(N+1)){
    X[i] = X[i-1] + A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]
    plot(t,X,type="n",ylab=expression(X[t]),xlab="time",las=1,main=expression(dX[t]== (0.03*t*X[t]-X[t]^5)*dt+0.1*dW[t]))
    points(t,X,type="l",col="blue")
    mtext(paste(date()),adj=0.5,col="blue",line=3.0,cex=1.4)
    mtext(bquote(t[i]==.(round(t[i],2))),adj=0,col="red",line=0.32,cex=1.2)
    mtext(bquote(X[t[i]]==.(round(X[i],2))),adj=0.2,col="red",line=0.2,cex=1.2)
    mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
    mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
    mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 23:12 04/04/2010"),side = 1, line = 4, adj = 0.5, cex = .66)
    }
par(pa)

