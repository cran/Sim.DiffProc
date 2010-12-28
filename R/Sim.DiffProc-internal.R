.ABM <-
function(N,t0,T,x0,theta,sigma,output=FALSE)
         {

if ( N <= 0 ) 
                       stop(tkmessageBox(title="Error",message=paste( "size of process : N >>> 0" ),icon="error"))

if ( t0 >= T || t0 < 0 ) 
                       stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if (sigma < 0 ) 
                       stop(tkmessageBox(title="Error",message=paste( "constant positive : Sigma > 0" ),icon="error"))

temps = seq(t0,T,length=N+1)
dt = (T-t0)/N
u = runif(N,0,1)
                 for (i in 1:length(u))
                 {
                 if ( u[i] >= 0.5)
                 u[i] = +1 
                 else
                 u[i] = -1
                 }
w = cumsum(c(0,u))*sqrt(dt)
X <- vector()
X[1] <- x0
for (i in 1:N){X[i+1] <- X[i]+ theta*dt + sigma*(w[i+1]-w[i])}
plot(temps,X,las=1,type="n",xlab="time",ylab=expression(X[t]))
points(temps,X,type="l",col="black",lwd=1,panel.frist=grid(col="gray"))
mtext("Arithmetic Brownian Motion",line=2,cex=1.2)
mtext(bquote(dX[t]==.(theta)*dt+.(sigma)*dW[t]),line=0.4,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(dt)),line=0.2,cex=1,adj=0,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 15:54 31/03/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "Arithmetic Brownian Motion.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.ABMF <-
function(N,M,t0,T,x0,theta,sigma,output=FALSE)
     {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 2)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 2" ),icon="error"))

if (sigma <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

ABM <- function(N,t0,T,x0,theta,sigma)
         {
temps = seq(t0,T,length=N+1)
dt = (T-t0)/N
u = runif(N,0,1)
                 for (i in 1:length(u))
                 {
                 if ( u[i] >= 0.5)
                 u[i] = +1 
                 else
                 u[i] = -1
                 }
w = cumsum(c(0,u))*sqrt(dt)
X <- vector()
X[1] <- x0
for (i in 1:N){X[i+1] <- X[i]+ theta*dt + sigma*(w[i+1]-w[i])}
X
         }
Q = sapply(rep(N,length=M),ABM,t0=t0,T=T,x0=x0,theta=theta,sigma=sigma)
temps = seq(t0,T,length=N+1)
dt = (T-t0)/N
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(temps,Q[,1],type="n",las=1,xlab="time",ylim=c(r1,r2),ylab=expression(X[t]),cex.lab=1)
mtext("Flow of Arithmetic Brownian Motion",line=2.5,cex=1.2)
mtext(bquote(dX[t]==.(theta)*dt+.(sigma)*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(NbrT==.(M)),line=1.2,cex=0.9,adj=0,col="red")
mtext(bquote(Delta*t==.(dt)),line=0.2,cex=1,adj=0,col="red")
for (i in 1:M){points(temps,Q[,i],type="l",panel.frist=grid(col="gray"))}
if (M >=2) {lines(temps,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 15:54 31/03/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q
time <- temps
X.mean <- Q.mean
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,X,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "ABMF.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.Asys <-
function(lambda,mu,t,T)
         {
if ( lambda <= 0 || mu <= 0 ) 
                       stop(tkmessageBox(title="Error",message=paste( "Lambda > 0 and  mu > 0" ),icon="error"))

if ( t <= 0 || T <= 0 ) 
                       stop(tkmessageBox(title="Error",message=paste( "t > 0 and  T > 0" ),icon="error"))

lambda = 2
mu = 0.5
t = 1
T = 10
Q <- matrix(c(-lambda,mu,lambda,-mu),ncol=2)
V <- eigen(Q)
G <- matrix(c(V$vectors[,1],V$vectors[,2]),ncol=2)
invG <- solve(G)
D1 <- expression ( exp(V$values[1]*t) )
D2 <- expression ( exp(V$values[2]*t) )
D  <- expression ( matrix(c(eval(D1),0,0,eval(D2)),ncol=2) ) 
P_t <- G%*%eval(D)%*%invG
temps =seq(0,T,by=0.1)
f1 <- numeric()
f2 <- numeric()
         for ( i in 1:length(temps))
             {
             f1[i] = (mu/(mu+lambda)) + (lambda/(mu+lambda))*exp(-(lambda+mu)*temps[i])
             f2[i] = (mu/(mu+lambda)) - (mu/(mu+lambda))    *exp(-(lambda+mu)*temps[i])
             }
plot(temps,f1,type="l",ylim=c(0,1),las=1,lwd=2,col="red",xlab="time",ylab=expression(pi[0](t)))
points(temps,f1,type="l",col="red",lwd=2)
points(temps,f2,type="l",col="blue",lwd=2,panel.frist=grid(col="gray"))
mtext("Evolution a telegraphic process in time",line=2,cex=1.2)
mtext(bquote(mu==.(mu)),adj=0,line=0.25,cex=1,col="red")
mtext(bquote(lambda==.(lambda)),adj=0.25,line=0.25,cex=1,col="blue")
abline(h=mu/(mu+lambda),lwd=2,col="gray50",lty=2)
axis(2,at=round(mu/(mu+lambda),2),las=1,col.axis="gray50")
text(T/4,(mu/(mu+lambda))+0.2,c(expression(pi[0]==(list(1,0)))),cex=0.8,col="red")
text(T/4,(mu/(mu+lambda))-0.2,c(expression(pi[0]==(list(0,1)))),cex=0.8,col="blue")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 01:33 31/01/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Result <- data.frame(P_t)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
return(invisible(Result)) 
        }

.BB <-
function(N,t0,T,x0,y,output=FALSE)
    {
if ( t0 >= T || t0 < 0 ) 
                       stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
x = c(0,x)
w = cumsum(x)*sqrt(delta.temps)
X <- x0 + w - (temps-t0)/(T-t0) * (w[N+1]-y +x0)
plot(temps,X,las=1,type="n",xlab="time",ylab=expression(X[t]))
points(temps,X,type="l",col="black",lwd=1,panel.frist=grid(col="gray"))
mtext("Brownian Bridge",line=2,cex=1.2)
mtext(bquote(x[.(0)]==.(x0)),line=0.15,cex=1.2,adj=0,col="red")
mtext(bquote(y==.(y)),line=0.1,cex=1.2,adj=0.2,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.3,cex=1,adj=1,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 00:47 18/03/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "Brownian Bridge.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
    }

.BBF <-
function(N,M,t0,T,x0,y,output=FALSE)
      {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 2)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 2" ),icon="error"))

BB <- function(N,t0,T,x0,y)
    {
temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
x = c(0,x)
w = cumsum(x)*sqrt(delta.temps)
X <- x0 + w - (temps-t0)/(T-t0) * (w[N+1]-y +x0)
}
Q = sapply(rep(N,length=M),BB,t0=t0,T=T,x0=x0,y=y)
temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(temps,Q[,1],type="n",las=1,xlab="time",ylim=c(r1,r2),ylab=expression(X[t]),cex.lab=1)
for (i in 1:M){points(temps,Q[,i],type="l",panel.frist=grid(col="gray"))}
mtext("flow of the Brownian bridge",line=2.5,cex=1.2)
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1.2,adj=0,col="red")
mtext(bquote(y==.(y)),line=0.1,cex=1.2,adj=0.2,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1.2,adj=0.4,col="red")
if ( M >= 2 ) {lines(temps,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2,cex=1.1)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 04:11 18/03/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BBF.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
    }

.Besselp <-
function(N,M,t0,T,x0,alpha,output=FALSE)
       {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( x0 == 0) 
            stop(tkmessageBox(title="Error",message=paste( "x0 =! 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( "N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if (alpha < 2 ) 
            stop(tkmessageBox(title="Error",message=paste( "alpha >= 2" ),icon="error"))

Bes <- function(N,T,t0,x0,alpha)
   {
Dt <- (T-t0)/N
a <- expression((alpha-1)/(2*x))
s <- expression(1)
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
       0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X     
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),Bes,T=T,t0=t0,x0=x0,alpha)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("Bessel Process",adj=0.5,line=3,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==frac(.(alpha-1),2*X[t])*dt+dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 29/08/2010 02:39:15"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "Bessel.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
    }

.BMcov <-
function(N,M,T,C)
       {
if( T <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( " T > 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( " M must be >= 1" ),icon="error"))

if (C <= 0) 
            stop(tkmessageBox(title="Error",message=paste( " C >= 1" ),icon="error"))

temps = seq(0,T,length=N)
delta.temps = T/N
TB = matrix(rnorm((N-1)*M,sd=sqrt(C*delta.temps)),nrow=M)
B = matrix(NA,ncol=N,nrow=M)
for (i in 1:M){B[i,] = c(0,cumsum(TB[i,]))}
B.cov = cov(B) 
filled.contour(temps, temps,B.cov, col = terrain.colors(10),plot.title = 
               title(main = "Empirical Covariance of BM",xlab = "time",
               ylab = "time"),key.title =title(main=bquote(cov(BM[t]))))
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 22:45 20/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
       }

.BMinf <-
function(N,T)
      {
temps = seq(0,T,length=N+1)
delta.temps = T/N
u = runif(N,0,1)
                 for (i in 1:length(u))
                 {
                 if ( u[i] >= 0.5)
                 u[i] = +1
                 else
                 u[i] = -1
                 }
w = (cumsum(c(0,u))*sqrt(delta.temps))
limB1 = w/temps
plot(temps,w,las=1,lwd=1,type="l",xlab="t",ylab=expression(W[t]))
points(temps,w,type="n",panel.frist=grid(col="gray"))
points(temps[-1],limB1[-1],type="l",lwd=2,col="green")
mtext("Standard Brownian Motion has the infinite",line=2,cex=1.2)
legend("topleft",bg="gray",border="gray",c(expression(lim(frac(w[t],t),t%->%+infinity)%~~%0)),lty=c(1),col=c("green"),lwd=3)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 19:46 21/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
}

.BMIrt <-
function(N,T)
      {
temps = seq(0,T,length=N+1)
delta.temps = T/N
u = runif(N,0,1)
                 for (i in 1:length(u))
                 {
                 if ( u[i] >= 0.5)
                 u[i] = +1
                 else
                 u[i] = -1
                 }
w = (cumsum(c(0,u))*sqrt(delta.temps))
i = seq(1,N+1,1)
x =  w[N-i+2]-w[N+1]
r1 = min(min(w),min(x))
r2 = max(max(w),max(x))
plot(temps,w,type="l",ylim=c(r1,r2),col="black",las=1,xlab="time",ylab="B(t)")
points(temps,x,col="red",type="l",panel.frist=grid(col="gray"))
mtext("Brownian Motion invariance by reversal of time",line=2,cex=1.2)
legend("topleft",bg="gray",border="gray",
c("B(t)","B(t)=B(T-t)-B(T)"),lty=c(1,1),col=c("black","red"),lwd=1)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 01:22 21/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
}

.BMIto1 <-
function(N,T,output=FALSE)
       {
temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
Ito <- 0.5*w^2 - 0.5 * temps
Ito2 <- vector()
for (i in 1:(N+1)){ Ito2[i] <- w[i]*(w[i+1]-w[i])}
r1= max(Ito)
r2= min(Ito)
Ito.sum <- cumsum(Ito2)
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Ito.sum,type="l",col="red",panel.frist=grid(col="gray"))
mtext(c((expression("Stochastic Integral":I(w[t])==integral(W[s] * dW[s], 0, t)))),adj=0.5,cex=1.2)
mtext(bquote(Delta*t==.(delta.temps)),line=0.2,cex=1,adj=0,col="red")
legend("topleft",bg="gray",border="gray",c(expression(frac(1,2)*(w[t]^2-t)),expression(sum(w[t[i]]*(w[t[i+1]]-w[t[i]]),i=0,n))),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 18:44 24/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Ito,Ito.sum)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMIto1.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.BMIto2 <-
function(N,T,output=FALSE)
{
temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
Ito <- 0.5*w^2 - 0.5 * temps
Ito.sum <- c(0, sapply (2:(N+1), function (x) {w[x] * (w[x+1]-w[x])} ) )
r1= max(Ito)
r2= min(Ito)
Ito.sum <- cumsum(Ito.sum)
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Ito.sum,type="l",col="red",panel.frist=grid(col="gray"))
mtext(c((expression("Stochastic integral":I(w[t])==integral(W[s] * dW[s], 0, t)))),adj=0.5,cex=1.2)
mtext(bquote(Delta*t==.(delta.temps)),line=0.2,cex=1,adj=0,col="red")
legend("topleft",bg="gray",border="gray",c(expression(frac(1,2)*(w[t]^2-t)),expression(sum(w[t[i]]*(w[t[i+1]]-w[t[i]]),i=0,n))),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 18:44 24/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Ito,Ito.sum)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMIto2.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.BMItoC <-
function(N,T,alpha,output=FALSE)
        {
if ( alpha == 0 )
            stop(tkmessageBox(title="Error",message=paste( "alpha =! 0" ),icon="error"))

temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
Ito <- alpha*w
Ito.sum <- sapply (1:(N+1), function (x) {alpha * (w[x+1]-w[x])} ) 
r1= min(min(Ito,na.rm=T),min(Ito.sum,na.rm=T))
r2= max(max(Ito,na.rm=T),max(Ito.sum,na.rm=T))
plot(temps,Ito,type="l",las=1,col="blue",ylim=c(r1,r2),xlab="time",ylab=expression(I(w[t])),main=bquote("Stochastic Integral":I(w[t])==alpha*integral(dW[s], 0, t)),cex.lab=1.1)
points(temps,cumsum(Ito.sum),type="l",col="red",panel.frist=grid(col="gray"))
mtext(bquote(alpha==.(alpha)),line=0.25,cex=1.2,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=0,col="red")
legend("topleft",bg="gray",border="gray",c(bquote(alpha*w[t]),expression(sum(alpha*(w[t[i+1]]-w[t[i]]),i=0,n))),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 20:27 25/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Ito.sum <- cumsum(Ito.sum)
time <- temps
Result <- data.frame(time,Ito,Ito.sum)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMItoC.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.BMItoP <-
function(N,T,power,output=FALSE)
       {
if ( power  <= 0)  
            stop(tkmessageBox(title="Error",message=paste( "power  > 0 " ),icon="error"))

temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
Ito <- (1/(power+1))*w^(power +1) - (power /2) * cumsum(w^(power -1)*delta.temps)
Ito.sum <- sapply (1:(N+1), function (x) {w[x]^power * (w[x+1]-w[x])} )
r1= min(min(Ito,na.rm=T),min(Ito.sum,na.rm=T))
r2= max(max(Ito,na.rm=T),max(Ito.sum,na.rm=T))
plot(temps,Ito,type="l",las=1,col="blue",ylim=c(r1,r2),ylab=expression(I(w[t])),xlab="time",main=bquote("Stochastic Integral":I(w[t])==integral(W[s]^n * dW[s], 0, t)),cex.lab=1.1)
points(temps,cumsum(Ito.sum),type="l",col="red",panel.frist=grid(col="gray"))
mtext(bquote(n==.(power)),line=0.25,cex=1.2,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=0,col="red")
legend("topleft",bg="gray",border="gray",c(expression(frac(w[t]^(n+1),n+1)-frac(n,2)*integral(W[s]^(n-1) * ds, 0, t)),expression(sum(w[t[i]]^n*(w[t[i+1]]-w[t[i]]),i=0,N))),
      lty=c(1,1),col=c("blue","red"),lwd=2,cex=0.85)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 20:24 26/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Ito.sum <- cumsum(Ito.sum)
time <- temps
Result <- data.frame(time,Ito,Ito.sum)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMItoP.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.BMItoT <-
function(N,T,output=FALSE)
       {
temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
for (i in 1:length(x)){if ( x[i] >= 0.5)
                       x[i] = +1 
                       else
                       x[i] = -1}
w = cumsum(c(0,x))*sqrt(delta.temps)
Ito <- temps*w - cumsum(w*delta.temps)
Ito.sum <- sapply (1:(N+1), function (x) {temps[x]*(w[x+1]-w[x])} )
r1= min(min(Ito,na.rm=T),min(Ito.sum,na.rm=T))
r2= max(max(Ito,na.rm=T),max(Ito.sum,na.rm=T))
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",main=bquote("Stochastic Integral":I(w[t])==integral(s*dW[s], 0, t)),cex.lab=1.1)
points(temps,cumsum(Ito.sum),type="l",col="red",panel.frist=grid(col="gray"))
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=0,col="red")
legend("topleft",bg="gray",border="gray",c(expression(t*w[t]-integral(w[s]*ds, 0, t)),expression(sum(t[i]*(w[t[i+1]]-w[t[i]]),i=0,n))),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 23:03 26/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Ito.sum <- cumsum(Ito.sum)
time <- temps
Result <- data.frame(time,Ito,Ito.sum)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMItoT.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.BMN <-
function(N,t0,T,C,output=FALSE)
                 {

if ( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if (C <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "C >= 1" ),icon="error"))

temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
MB = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(C*delta.temps))))
plot(temps,MB,las=1,lwd=1,type="l",xlab="time",ylab=expression(X[t]))
points(temps,MB,type="n",panel.frist=grid(col="gray"))
mtext("Brownian Motion",line=2,cex=1.2)
mtext("by normal law",line=0.15,cex=1.2,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=1,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 01:36 11/12/2009"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- MB
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMN.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
                 }

.BMNF <-
function(N,M,t0,T,C,output=FALSE)
     {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M <= 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 2" ),icon="error"))

if (C <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "C >= 1" ),icon="error"))

MB <- function(N,t0,T)
     {
temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
BM = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(C*delta.temps))))
      }
Q = sapply(rep(N,length=M),MB,t0=t0,T=T)
temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(temps,Q[,1],type="n",las=1,xlab="time",ylab=expression(X[t]),ylim=c(r1,r2),cex.lab=1)
mtext("flow of Brownian Motion",line=2.5,cex=1.2)
mtext("by normal law",line=0.25,cex=1,adj=0,col="red")
mtext(bquote(C==.(C)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.9,cex=1,adj=1,col="red")
for ( i in 1:M){points(temps,Q[,i],type="l",lwd=1,panel.frist=grid(col="gray"))}
if(M > 1) {lines(temps,Q.mean,lwd=2,col="blue")}
if(M > 1) {legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 02:26 25/01/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q
time <- temps
X.mean <- Q.mean
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,X,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMNF.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
}

.BMP <-
function(N,M,T,C)
     {

if( T <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > 0" ),icon="error"))

if (C <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "C >= 1" ),icon="error"))

if (M <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1 " ),icon="error"))

if (N <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "N >>> 0 " ),icon="error"))

MB <- function(N,T,C)
     {
temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(C*delta.temps)
      }
Q = sapply(rep(N,length=M),MB,T=T,C=C)
temps = seq(0,T,length=N+1)
f1 <- function(t) 2*sqrt(C*t)
plot(temps,f1(temps),ylim=c(min(-f1(temps)),max(f1(temps))),type="l",las=1,lwd=4,col="red",xlab="time",ylab=expression(W[t]),main=bquote("Trajectories Brownian in the curves"%+-%2*sqrt(C*t)))
points(temps,-f1(temps),type="l",col="red",lwd=4,panel.first=grid(col="gray"))
mtext(bquote(C==.(C)),line=0.25,cex=1.2,adj=1,col="red")
mtext(bquote("Numbers of the trajectories"==.(M)),line=0.25,cex=1,adj=0,col="blue")
for ( i in 1:M) { points(temps,Q[,i],type="l",lwd=1)}
legend("topleft",bg="gray",border="gray",c(expression(""%+-%2*sqrt(C*t))),lty=c(1),col=c("red"),lwd=2,cex=1.2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 23:53 20/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
}

.BMRW <-
function(N,t0,T,C,output=FALSE)
      {

if ( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if (C <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "C >= 1" ),icon="error"))

temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
x = c(0,x)
w = cumsum(x)*sqrt(C*delta.temps)
plot(temps,w,las=1,type="n",xlab="time",ylab=expression(X[t]))
points(temps,w,type="l",col="black",lwd=1,panel.frist=grid(col="gray"))
mtext("Brownian Motion",line=2,cex=1.2)
mtext("by a Random Walk",line=0.15,cex=1.2,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=1,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 22:28 31/01/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- w
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMRW.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
    }

.BMRWF <-
function(N,M,t0,T,C,output=FALSE)
     {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( C <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "C >= 1" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M <= 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 2" ),icon="error"))

MB1 <- function(N,t0,T)
     {
temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(C*delta.temps)
      }
Q = sapply(rep(N,length=M),MB1,t0=t0,T=T)
temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(temps,Q[,1],type="n",las=1,xlab="time",ylab=expression(X[t]),ylim=c(r1,r2),cex.lab=1)
mtext("flow of Brownian Motion",line=2.5,cex=1.2)
mtext("by a Random Walk",line=0.25,cex=1,adj=0,col="red")
mtext(bquote(C==.(C)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.9,cex=1,adj=1,col="red")
for ( i in 1:M){points(temps,Q[,i],type="l",lwd=1,panel.frist=grid(col="gray"))}
if(M > 1) {lines(temps,Q.mean,lwd=2,col="blue")}
if(M > 1) {legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 02:26 25/01/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q
time <- temps
X.mean <- Q.mean
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,X,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMRWF.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.BMscal <-
function(N,T,S1,S2,S3,output=FALSE)
      {
temps1 <- (S1^2)*seq(0,T,length=N+1)
temps2 <- (S2^2)*seq(0,T,length=N+1)
temps3 <- (S3^2)*seq(0,T,length=N+1)
r <- max(max(temps1),max(temps2),max(temps3))
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1
                 else
                 x[i] = -1
                 }
w1 = (1/S1)*(cumsum(c(0,x))*sqrt(delta.temps))
w2 = (1/S2)*(cumsum(c(0,x))*sqrt(delta.temps))
w3 = (1/S3)*(cumsum(c(0,x))*sqrt(delta.temps))
r1 <- min(min(w1),min(w2),min(w3))
r2 <- max(max(w1),max(w2),max(w3))
plot(temps1,w1,type="l",ylim=c(r1,r2),xlim=c(0,r),col="black",las=1,xlab="time",ylab=expression(W[t]))
points(temps2,w2,col="red",type="l",panel.frist=grid(col="gray"))
points(temps3,w3,col="blue",type="l")
mtext("Brownian Motion with different scales",line=2,cex=1.2)
legend("topleft",bg="gray",border="gray",
c(paste("S1=",S1),paste("S2=",S2),paste("S3=",S3)),lty=c(1,1,1),col=c("black","red","blue"),lwd=1)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 00:42 21/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Result <- data.frame(w1,w2,w3)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMscal.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.BMStra <-
function(N,T,output=FALSE)
        {
temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
Stra <- 0.5 * w^2
r1= max(Stra)
r2= min(Stra)
plot(temps,Stra,type="n",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red",panel.frist=grid(col="gray"))
mtext(c((expression("Stratonovitch integral":I(w[t])==integral(W[s]*o*dW[s], 0, t)))),adj=0.5,cex=1.2)
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=0,col="red")
legend("topleft",bg="gray",border="gray",c(expression(integral(W[s]*o*dW[s], 0, t)==frac(1,2)*W[t]^2)),
      lty=c(1),col=c("red"),lwd=2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 15:35 05/10/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMStra.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.BMStraC <-
function(N,T,alpha,output=FALSE)
        {

if ( alpha == 0 )
            stop(tkmessageBox(title="Error",message=paste( "alpha =! 0" ),icon="error"))

temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
Stra <- alpha * w
r1= max(Stra)
r2= min(Stra)
plot(temps,Stra,type="n",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red",panel.frist=grid(col="gray"))
mtext(c((expression("Stratonovitch integral":I(w[t])==integral(alpha*o*dW[s], 0, t)))),adj=0.5,cex=1.2)
legend("topleft",bg="gray",border="gray",c(expression(integral(alpha*o*dW[s], 0, t)==alpha*W[t])),
      lty=c(1),col=c("red"),lwd=2)
mtext(bquote(alpha == .(alpha)), line = 0.25, cex = 1.2,adj=0, col = "red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=1,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 15:15 05/10/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMStraC.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.BMStraP <-
function(N,T,power,output=FALSE)
        {

if ( power  <= 0)  
            stop(tkmessageBox(title="Error",message=paste( "power > 0" ),icon="error"))

temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
n <- power
Stra <- numeric(N+1)
for ( i in 1:N){
   Stra[i] <- sum(0.5*(w[i]^(n-1)+w[i+1]^(n-1))*(w[i+1]^2 - w[i]^2))
}
Stra <- cumsum(Stra)
r1= max(Stra)
r2= min(Stra)
plot(temps,Stra,type="n",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red",panel.frist=grid(col="gray"))
mtext(c((expression("Stratonovitch integral":I(w[t])==integral(W[s]^n*o*dW[s], 0, t)))),adj=0.5,cex=1.2)
legend("topleft",bg="gray",border="gray",c(expression(integral(W[s]^n*o*dW[s], 0, t))),
      lty=c(1),col=c("red"),lwd=2)
mtext(bquote(n==.(power)),line=0.25,cex=1.2,col="red",adj=0)
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=1,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 16:03 05/10/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMStraP.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.BMStraT <-
function(N,T,output=FALSE)
        {
temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
Stra <- numeric(N+1)
for ( i in 1:N){
   Stra[i] <- sum(0.5*(temps[i]*(w[i+1]-w[i])+temps[i+1]*(w[i+1]-w[i])))
}
Stra <- cumsum(Stra)
r1= max(Stra)
r2= min(Stra)
plot(temps,Stra,type="n",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red",panel.frist=grid(col="gray"))
mtext(c((expression("Stratonovitch integral":I(w[t])==integral(s*o*dW[s], 0, t)))),adj=0.5,cex=1.2)
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=0,col="red")
legend("topleft",bg="gray85",border="gray",c(expression(integral(s*o*dW[s], 0, t))),
      lty=c(1),col=c("red"),lwd=2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 16:25 05/10/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "BMStraT.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.CEV <-
function(N,M,t0,T,x0,mu,sigma,gamma,output=FALSE)
       {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( x0 == 0) 
            stop(tkmessageBox(title="Error",message=paste( "x0 =! 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if ( gamma <= 0)
            stop(tkmessageBox(title="Error",message=paste( "Gamma > 0" ),icon="error"))

CE <- function(N,T,t0,x0,mu,sigma,gamma)
   {
Dt <- (T-t0)/N
a <- expression(mu*x)
s <- expression(sigma*(x^gamma))
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
       0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X     
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),CE,T=T,t0=t0,x0=x0,mu,sigma,gamma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("Constant Elasticity of Variance process (CEV)",adj=0.5,line=2.5,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==.(mu)*X[t]*dt+.(sigma)*X[t]^.(gamma)*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 29/08/2010 02:39:15"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "CEV.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
    }

.CIR <-
function(N,M,t0,T,x0,theta,r,sigma,output=FALSE)
       {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if (x0 <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "x0 > 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if (sigma <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if ( r <= 0)
            stop(tkmessageBox(title="Error",message=paste( "r > 0" ),icon="error"))

if ( 2*r <= sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*r > (sigma)^2" ),icon="error"))

if (theta < 0 )
            stop(tkmessageBox(title="Error",message=paste( "Theta >= 0" ),icon="error"))


CI <- function(N,T,t0,x0,theta,r,sigma)
   {
Dt <- (T-t0)/N
a <- expression(r-theta*x)
s <- expression(sigma*sqrt(x))
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
       0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X       
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),CI,T=T,t0=t0,x0=x0,theta=theta,r=r,sigma=sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("Cox-Ingersoll-Ross process (CIR)",adj=0.5,line=2.5,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==(.(r)-.(theta)*X[t])*dt+.(sigma)*sqrt(X[t])*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 29/08/2010 01:46:33"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "CIR.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
    }

.CIRhy <-
function(N,M,t0,T,x0,r,sigma,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if (sigma <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if ( r+ ((sigma^2)/2) <= 0)
            stop(tkmessageBox(title="Error",message=paste( "r+(sigma^2)/2 > 0" ),icon="error"))

CHY <- function(N,T,t0,x0,r,sigma)
   {
Dt <- (T-t0)/N
a <- expression(-r*x)
s <- expression(sigma*sqrt(1+x^2))
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
       0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X      
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),CHY,T=T,t0=t0,x0=x0,r=r,sigma=sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("The modified CIR and hyperbolic processes",adj=0.5,line=2.5,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==-.(r)*X[t]*dt+.(sigma)*sqrt(1+X[t]^2)*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 29/08/2010 03:51:12"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "CIRhy.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))   
}

.CKLS <-
function(N,M,t0,T,x0,r,theta,sigma,gamma,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if (sigma <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if ( gamma < 0)
            stop(tkmessageBox(title="Error",message=paste( "gamma > 0" ),icon="error"))

CK <- function(N,T,t0,x0,r,theta,sigma,gamma)
   {
Dt <- (T-t0)/N
a <- expression(r+theta*x)
s <- expression(sigma*(x^gamma))
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
       0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X      
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),CK,T=T,t0=t0,x0=x0,r=r,theta=theta,sigma=sigma,gamma=gamma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("Chan-Karolyi-Longstaff-Sanders (CKLS)",adj=0.5,line=2.5,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==(.(r)+.(theta)*X[t])*dt+.(sigma)*X[t]^.(gamma)*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 29/08/2010 03:11:53"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "CKLS.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result)) 
    }

.diffBridge <-
function(N,t0,T,x,y,drift,diffusion,Output=FALSE) 
             {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if(!is.expression(drift)) 
            stop(tkmessageBox(title="Error",message=paste( " The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(diffusion))
            stop(tkmessageBox(title="Error",message=paste( " The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

Euler <- function(N,T,t0,x0,drift,diffusion)
     {
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diffusion)
t = seq(t0,T,length=N+1)
Dt = (T-t0)/N
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){X[i] = X[i-1] + A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]}      
X <- ts(X, start =t0 ,deltat =Dt)
X
   }
t = seq(t0,T,length=N+1)
Dt = (T-t0)/N
done <- FALSE
while (!done) {
              X1 <-Euler(N=N,T=T,t0=t0,x0=x,drift=drift,diffusion=diffusion)
              X2 <-Euler(N=N,T=T,t0=t0,x0=y,drift=drift,diffusion=diffusion)    
              X3 <-ts(rev(X2),start=start(X2),end=end(X2),deltat=deltat(X2))
G <- Inf
if (X1[1] >= X3[1]) {
                     if (!all(X1 > X3)) 
                     G <- min(which(X1 <= X3)) - 1
                    }
               else {
                    if (!all(X1 < X3)) 
                    G <- min(which(X1 >= X3)) - 1
                    }
                    if (G == 0 || G == length(X1) || G == Inf) {
                    stop(tkmessageBox(title="Information",message=paste( "A crossing has been no realized,trying again (Repeat)..." ),icon="info"))
                    done <- FALSE
                    }
               else {
                    done <- TRUE
                    }
            }
X <- ts(c(X1[1:G],X3[-(1:G)]),start = t0 ,end=end(X1),frequency = frequency(X1) )
plot(time(X),X,type="n",ylab=expression(X[t]),xlab="time",las=1)
points(time(X),X,type="l",col="black",lwd=1,panel.frist=grid(col="gray"))
mtext(expression("Diffusion Bridges Process"),line=3,cex=1.2,adj=0.5)
mtext(expression(dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=1.5,cex=1.2,adj=0.5,col="blue")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=0.2)
mtext(expression(sigma(t,X[t])==.),adj=0.5,col="red",line=0.2)
mtext(drift,adj=0.17,col="red",line=0.17)
mtext(diffusion,adj=0.65,col="red",line=0.17)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.9,adj=1,col="red")
mtext(bquote(Start==.(x)),line=1.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=2.7,cex=0.9,adj=1,col="red")
mtext(bquote(End==.(y)),line=2,cex=0.9,adj=1,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 04/09/2010 00:47:19"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- time(X)
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.xlsx(Result,"diffBridge.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
if ( donc <- TRUE ) {return(invisible(Result))
}
}

.DWP <-
function(N,M,t0,T,x0,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

DW <- function(N,T,t0,x0)
   {
Dt <- (T-t0)/N
a <- expression(x-x^3)
s <- expression(1)
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
       0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X      
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),DW,T=T,t0=t0,x0=x0)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("Double-Well Potential",adj=0.5,line=2.5,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==(X[t]-X[t]^3)*dt+dW[t] ),cex=1,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 30/08/2010 00:58:47"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "DWP.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
    }

.Euler <-
function(N,M,T=1,t0,x0,Dt,a,sigma,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(a))
            stop(tkmessageBox(title="Error",message=paste( "The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(sigma))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

Eul <- function(N,T=1,Dt,t0,x0,a,sigma)
     {

A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(sigma)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){X[i] = X[i-1] + A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]}      
X
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Eul,T=T,t0=t0,x0=x0,Dt=Dt,a=a,sigma=sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Euler scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria",date()),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.xlsx(Result, "Euler.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
    }

.GBM <-
function(N,t0,T,x0,theta,sigma,output=FALSE)
      {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( x0 <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "x0 > 0" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

temps = seq(t0,T,length=N+1)
dt = (T-t0)/N
u = runif(N,0,1)
                 for (i in 1:length(u))
                 {
                 if ( u[i] >= 0.5)
                 u[i] = +1 
                 else
                 u[i] = -1
                 }
w = cumsum(c(0,u))*sqrt(dt)
X <- length(temps)
X <- x0*exp((theta - 0.5*sigma^2)*temps + sigma*w)
plot(temps,X,las=1,type="n",xlab="time",ylab=expression(X[t]))
points(temps,X,type="l",col="black",lwd=1,panel.frist=grid(col="gray"))
mtext(" Brownian Motion Geometrical (Model of Black-Scholes)" ,line=2.5,cex=1.2)
mtext(bquote(dX[t]==.(theta)*X[t]*dt+.(sigma)*X[t]*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(dt)),line=0.4,cex=1,adj=0,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 22:31 31/03/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "GBM.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.GBMF <-
function(N,M,t0,T,x0,theta,sigma,output=FALSE)
        {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( x0 <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "x0 > 0" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M <= 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 2" ),icon="error"))

MBG <- function(N,T,t0,x0,theta,sigma)
      {
temps = seq(t0,T,length=N+1)
dt = (T-t0)/N
u = runif(N,0,1)
                 for (i in 1:length(u))
                 {
                 if ( u[i] >= 0.5)
                 u[i] = +1 
                 else
                 u[i] = -1
                 }
w = cumsum(c(0,u))*sqrt(dt)
X <- length(temps)
X <- x0*exp((theta - 0.5*sigma^2)*temps + sigma*w)
X 
}
Q = sapply(rep(N,length=M),MBG,T=T,t0=t0,x0=x0,theta=theta,sigma=sigma)
temps = seq(t0,T,length=N+1)
dt = (T-t0)/N
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(temps,Q[,1],las=1,type="n",ylim=c(r1,r2),xlab="time",ylab=expression(X[t]),cex.lab=1)
mtext("Flow of Brownian Motion Geometrical",line=2.5,cex=1.2)
mtext(bquote(dX[t]==.(theta)*X[t]*dt+.(sigma)*X[t]*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(dt)),line=0.4,cex=1,adj=0,col="red")
for (i in 1:M){points(temps,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
if (M >=2){lines(temps,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 11/04/2010 20:35:57"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q
time <- temps
X.mean <- Q.mean
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,X,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "GBMF.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.Heun <-
function(N,M,T=1,t0,x0,Dt,a,sigma,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(a))
            stop(tkmessageBox(title="Error",message=paste( "The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(sigma))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

H <- function(N,T=1,Dt,t0,x0,a,sigma)
   {
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(sigma)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
X    <- numeric()
Y    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
         Y[i-1]= X[i-1]+A(t[i-1],X[i-1])*Dt+S(t[i-1],X[i-1])*D[i-1]
         X[i]  = X[i-1]+0.5*Dt*(A(t[i-1],X[i-1])+A(t[i-1],Y[i-1]))+
                 0.5*(S(t[i-1],X[i-1])+S(t[i-1],Y[i-1]))*D[i-1]
               }          
X      
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),H,t0=t0,T=T,x0=x0,Dt=Dt,a=a,sigma=sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Heun scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria",date()),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.xlsx(Result, "Heun.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
    }

.HWV <-
function(N,t0,T,x0,theta,r,sigma,output=FALSE)
     {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if ( r <= 0)
            stop(tkmessageBox(title="Error",message=paste( "r > 0" ),icon="error"))

temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
u = runif(N,0,1)
for (i in 1:length(u)){if ( u[i] >= 0.5)
                       u[i] = +1
                       else
                       u[i] = -1}
w = cumsum(c(0,u))*sqrt(delta.temps)
Ito.sum <- c(0,sapply(1:(N+1),function(x){exp(-r*(temps[x+1]-temps[x]))*(w[x+1]-w[x])}))
X <- sapply(1:(N+1),function(x){theta+(x0-theta)*exp(-r*temps[x])+sigma*sum(Ito.sum[1:x])})
plot(temps,X,type="n",las=1,xlab="time",ylab=expression(X[t]))
points(temps,X,type="l",panel.frist=grid(col="gray"))
mtext("Hull-White/Vasicek (HWV) Gaussian Diffusion Models",cex=1.2,line=2.5)
mtext(bquote(dX[t]==.(r)*(.(theta)-X[t])*dt+.(sigma)*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 22:44 01/04/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "HWV.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.HWVF <-
function(N,M,t0,T,x0,theta,r,sigma,output=FALSE)
        {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if ( r <= 0)
            stop(tkmessageBox(title="Error",message=paste( "r > 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M <= 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 2" ),icon="error"))

OUG <- function(N,t0,T,x0,theta,r,sigma){
temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
u = runif(N,0,1)
for (i in 1:length(u)){if ( u[i] >= 0.5)
                       u[i] = +1
                       else
                       u[i] = -1}
w = cumsum(c(0,u))*sqrt(delta.temps)
Ito.sum <- c(0,sapply(1:(N+1),function(x){exp(-r*(temps[x+1]-temps[x]))*(w[x+1]-w[x])}))
X <- sapply(1:(N+1),function(x){theta+(x0-theta)*exp(-r*temps[x])+sigma*sum(Ito.sum[1:x])})
}
Q = sapply(rep(N,length=M),OUG,t0=t0,T=T,x0=x0,theta=theta,r=r,sigma=sigma)
temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(temps,Q[,1],type="n",las=1,xlab="time",ylim=c(r1,r2),ylab=expression(X[t]))
for (i in 1:M){points(temps,Q[,i],type="l",panel.frist=grid(col="gray"))}
mtext("Gaussian Diffusion Models(Hull-White/Vasicek)",line=2.5,cex=1.2)
mtext(bquote(dX[t]==.(r)(.(theta)-X[t])*dt+.(sigma)*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0,col="red")
if ( M >= 2 ) {lines(temps,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 00:59 02/04/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "HWVF.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
    }

.Hyproc <-
function(N,M,t0,T,x0,theta,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if (theta <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "theta > 0" ),icon="error"))

hy <- function(N,T,t0,x0,theta)
   {
Dt <- (T-t0)/N
a <- expression((-theta*x)/(sqrt(1+x^2)))
s <- expression(1)
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
        0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X      
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),hy,T=T,t0=t0,x0=x0,theta=theta)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("The hyperbolic process",adj=0.5,line=3,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==- frac(.(theta)*X[t],sqrt(1+X[t]^2))*dt+dW[t] ),cex=1,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 29/08/2010 18:46:26"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "Hyproc.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
    }

.Hyprocg <-
function(N,M,t0,T,x0,beta,gamma,theta,mu,sigma,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if ( gamma <= 0)
            stop(tkmessageBox(title="Error",message=paste( "Gamma > 0" ),icon="error"))

if ( theta < 0)
            stop(tkmessageBox(title="Error",message=paste( "theta >= 0" ),icon="error"))

if (abs(beta) < 0 || abs(beta) >= gamma )
            stop(tkmessageBox(title="Error",message=paste( "0 <= abs(beta) < gamma" ),icon="error"))
               
hyg <- function(N,T,t0,x0,beta,gamma,theta,mu,sigma)
   {
Dt <- (T-t0)/N
a <- expression( (0.5*sigma^2)* ( beta- ((gamma*x)/sqrt((theta^2)+(x-mu)^2)) ) )
s <- expression(sigma)
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
          0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X      
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),hyg,T=T,t0=t0,x0=x0,beta=beta,gamma=gamma,theta=theta,mu=mu,sigma=sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("The general hyperbolic diffusion",adj=0.5,line=3.2,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==frac(.(sigma)^2,2)* (.(beta)-frac(.(gamma)*X[t],sqrt(.(theta)^2+(X[t]-.(mu))^2) )) *dt+.(sigma)*dW[t] ),cex=1,adj=0.5,line=0.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 29/08/2010 19:37:25"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "Hyprocg.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
    }

.INFSR <-
function(N,M,t0,T,x0,theta,r,sigma,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( x0 <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "x0 > 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

IFS <- function(N,T,t0,x0,theta,r,sigma)
   {
Dt <- (T-t0)/N
a <- expression( x*(theta-((sigma^3)-theta*r)*x) )
s <- expression( sigma*(x^(3/2)) )
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
        0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X       
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),IFS,T=T,t0=t0,x0=x0,theta=theta,r=r,sigma=sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("Inverse of Feller Square Root Model",adj=0.5,line=2.5,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==X[t]*(.(theta)-(.(theta*r))*X[t])*dt+.(sigma)*X[t]^frac(3,2)*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 30/08/2010 00:31:36"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "INFSR.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
    }

.ItovsStra <-
function(N,T,output=FALSE)
          {
temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
Stra <- 0.5 * w^2
Ito <- 0.5*w^2 - 0.5 * temps
r1= max(max(Ito),max(Stra))
r2= min(min(Ito),min(Stra))
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),ylim=c(r2,r1),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red",panel.frist=grid(col="gray"))
mtext("Stochastic integral",adj=0.5,cex=1.2,line=3)
mtext(c((expression(I(w[t])==integral(W[s] * dW[s], 0, t)))),adj=0,cex=1,col="blue",line=0.4)
mtext(c((expression(I(w[t])==integral(W[s]*o*dW[s], 0, t)))),adj=1,cex=1,col="red",line=0.4)
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0.5,col="black")
legend("topleft",bg="gray85",border="gray",c("Ito Integral","Stratonovitch Integral"),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 21:51 05/10/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Ito,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "ItoStra.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.ItovsStraP <-
function(N,T,power,output=FALSE)
        {
if ( power  <= 0)  
            stop(tkmessageBox(title="Error",message=paste( "power  > 0" ),icon="error"))

temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
n <- power
Stra <- numeric(N+1)
for ( i in 1:N){
   Stra[i] <- sum(0.5*(w[i]^(n-1)+w[i+1]^(n-1))*(w[i+1]^2 - w[i]^2))
}
Stra <- cumsum(Stra)
Ito <- (1/(power+1))*w^(power +1) - (power /2) * cumsum(w^(power -1)*delta.temps)
r1= max(max(Ito),max(Stra))
r2= min(min(Ito),min(Stra))
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),ylim=c(r2,r1),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red",panel.frist=grid(col="gray"))
mtext("Stochastic integral",adj=0.5,cex=1.2,line=3)
mtext(c((expression(I(w[t])==integral(W[s]^n * dW[s], 0, t)))),adj=0,cex=1,col="blue",line=0.4)
mtext(c((expression(I(w[t])==integral(W[s]^n*o*dW[s], 0, t)))),adj=1,cex=1,col="red",line=0.4)
mtext(bquote(n==.(power)),line=1.4,cex=1,col="black",adj=0.5)
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0.5,col="black")
legend("topleft",bg="gray85",border="gray",c("Ito Integral","Stratonovitch Integral"),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 21:30 05/10/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Ito,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "ItoStraP.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.ItovsStraT <-
function(N,T,output=FALSE)
        {
temps = seq(0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:length(x))
                 {
                 if ( x[i] >= 0.5)
                 x[i] = +1 
                 else
                 x[i] = -1
                 }
w = cumsum(c(0,x))*sqrt(delta.temps)
Stra <- numeric(N+1)
for ( i in 1:N){
   Stra[i] <- sum(0.5*(temps[i]*(w[i+1]-w[i])+temps[i+1]*(w[i+1]-w[i])))
}
Stra <- cumsum(Stra)
Ito <- temps*w - cumsum(w*delta.temps)
r1= max(max(Ito),max(Stra))
r2= min(min(Ito),min(Stra))
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),ylim=c(r2,r1),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red",panel.frist=grid(col="gray"))
mtext("Stochastic integral",adj=0.5,cex=1.2,line=3)
mtext(c((expression(I(w[t])==integral(s * dW[s], 0, t)))),adj=0,cex=1,col="blue",line=0.4)
mtext(c((expression(I(w[t])==integral(s*o*dW[s], 0, t)))),adj=1,cex=1,col="red",line=0.4)
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0.5,col="black")
legend("topleft",bg="gray85",border="gray",c("Ito Integral","Stratonovitch Integral"),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 22:10 05/10/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Ito,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "ItoStraT.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.JDP <-
function(N,M,t0,T,x0,theta,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if (theta <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "theta > 0" ),icon="error"))

if ( x0 > 1 || x0 < 0 )
            stop(tkmessageBox(title="Error",message=paste( "0 =< x0 <= 1" ),icon="error"))

JD <- function(N,T,t0,x0,theta)
   {
Dt <- (T-t0)/N
a <- expression(-theta*(x-0.5))
s <- expression(sqrt(theta*x*(1-x)))
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]
              }
X      
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),JD,T=T,t0=t0,x0=x0,theta=theta)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("The Jacobi diffusion process",adj=0.5,line=3,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==-.(theta)*(X[t]-frac(1,2))*dt+sqrt(.(theta)*X[t]*(1-X[t]))*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 29/08/2010 04:09:16"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "JDP.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
    }

.MartExp <-
function(N,t0,T,sigma,output=FALSE)
      {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
u = runif(N,0,1)
                 for (i in 1:length(u))
                 {
                 if ( u[i] >= 0.5)
                 u[i] = +1
                 else
                 u[i] = -1
                 }
w = (cumsum(c(0,u))*sqrt(delta.temps))
X = w^2 - temps
Y = exp(sigma*w-0.5*(sigma^2)*temps)
X11()
plot(temps,X,las=1,lwd=1,type="n",xlab="time",ylab=expression(X[t]))
mtext(c(expression(X[t]==W[t]^2 - t )),cex=1.4,col="red",line=1)
points(temps,X,type="l",panel.frist=grid(col="gray"))
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 23:34 21/03/2010"),
     side = 1, line = 4, adj = 0.5, cex = .66)
X11()
plot(temps,Y,las=1,lwd=1,type="n",xlab="time",ylab=expression(Y[t]))
mtext(c(expression(Y[t]== exp(sigma*W[t] - frac(sigma^2,2)*t) )),cex=1.4,col="red",line=1)
points(temps,Y,type="l",panel.frist=grid(col="gray"))
mtext(bquote(sigma^2==.(sigma)^2),line=0.25,cex=1.2,adj=1,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 23:34 21/03/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,X,Y)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "MartExp.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.Milstein <-
function(N,M,T=1,t0,x0,Dt,a,sigma,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(a))
            stop(tkmessageBox(title="Error",message=paste( "The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(sigma))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

Mils <- function(N,T=1,Dt,t0,x0,a,sigma)
   {
DSx  <- D(sigma,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(sigma)
Sx   <- function(t,x)  eval(DSx)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
       0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X      
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Mils,t0=t0,T=T,x0=x0,Dt=Dt,a=a,sigma=sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Milstein scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria",date()),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.xlsx(Result, "Milstein.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
    }

.MilsteinS <-
function(N,M,T=1,t0,x0,Dt,a,sigma,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(a))
            stop(tkmessageBox(title="Error",message=paste( "The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(sigma))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

MilS <- function(N,T=1,Dt,t0,x0,a,sigma)
   {
DAx  <- D(a,"x")
DAxx <- D(D(a,"x"),"x")
DSx  <- D(sigma,"x")
DSxx <- D(D(sigma,"x"),"x")
A    <- function(t,x)  eval(a)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(sigma)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1] + A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1] +
       0.5 *S(t[i-1],X[i-1]) * Sx(t[i-1],X[i-1])*(D[i-1]^2-Dt)+ 
       Dt^(3 /2)*(0.5 *A(t[i-1],X[i-1])*Sx(t[i-1],X[i-1]) +
       0.5 *Ax(t[i-1],X[i-1])*S(t[i-1],X[i-1])+
       0.25 *(S(t[i -1] ,X[i -1])^2) * Sxx(t[i -1] ,X[i -1]))*D[i -1]+ 
       (Dt^2) * (0.5*A(t[i -1],X[i -1])*Ax(t[i-1],X[i-1])+
       0.25 *Axx(t[i-1],X[i-1])*(S(t[i-1],X[i-1])^2))
               }               
X      
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),MilS,t0=t0,T=T,x0=x0,Dt=Dt,a=a,sigma=sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Second Milstein scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria",date()),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.xlsx(Result, "MilsteinS.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
    }

.OU <-
function(N,t0,T,x0,r,sigma,output=FALSE)
     {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if ( r <= 0)
            stop(tkmessageBox(title="Error",message=paste( "r > 0" ),icon="error"))

temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
u = runif(N,0,1)
for (i in 1:length(u)){
                       if ( u[i] >= 0.5)
                       u[i] = +1
                       else
                       u[i] = -1
                      }
w = cumsum(c(0,u))*sqrt(delta.temps)
Ito.sum <- c(0,sapply(1:(N+1),function(x){exp(-r*(temps[x+1]-temps[x]))*(w[x+1]-w[x])}))
X <- sapply(1:(N+1),function(x){x0*exp(-r*temps[x])+sigma*sum(Ito.sum[1:x])})
plot(temps,X,type="n",las=1,xlab="time",ylab=expression(X[t]))
points(temps,X,type="l",panel.frist=grid(col="gray"))
mtext("Ornstein-Uhlenbeck Process",line=2.5,cex=1.2 )
mtext(bquote(dX[t]==-.(r)*X[t]*dt+.(sigma)*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 19:35 18/03/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "OU.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.OUF <-
function(N,M,t0,T,x0,r,sigma,output=FALSE)
        {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M <= 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 2" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if ( r <= 0)
            stop(tkmessageBox(title="Error",message=paste( "r > 0" ),icon="error"))

OU <- function(N,t0,T,x0,r,sigma)
     {
temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
u = runif(N,0,1)
for (i in 1:length(u)){if ( u[i] >= 0.5)
                       u[i] = +1
                       else
                       u[i] = -1}
w = cumsum(c(0,u))*sqrt(delta.temps)
Ito.sum <- c(0,sapply(1:(N+1),function(x){exp(-r*(temps[x+1]-temps[x]))*(w[x+1]-w[x])}))
X <- sapply(1:(N+1),function(x){x0*exp(-r*temps[x])+sigma*sum(Ito.sum[1:x])})
}
Q = sapply(rep(N,length=M),OU,t0=t0,T=T,x0=x0,r=r,sigma=sigma)
temps = seq(t0,T,length=N+1)
delta.temps = (T-t0)/N
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(temps,Q[,1],type="n",las=1,xlab="time",ylim=c(r1,r2),ylab=expression(X[t]))
mtext("Flow of Ornstein-Uhlenbeck Process",line=2.5,cex=1.2)
mtext(bquote(dX[t]==-.(r)*X[t]*dt+.(sigma)*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0,col="red")
for (i in 1:M){points(temps,Q[,i],type="l",panel.frist=grid(col="gray"))}
if (M >=2) {lines(temps,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 23:12 30/03/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "OUF.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
   }

.PDP <-
function(N,M,t0,T,x0,theta,mu,a,b,c,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if (theta <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "theta > 0" ),icon="error"))
            
pd <- function(N,T,t0,x0,theta,mu,a,b,c)
   {
Dt <- (T-t0)/N
aa <- expression( -theta*(x-mu) )
s <- expression( sqrt(2*theta*(a*x^2 +b*x+c)) )
DSx  <- D(s,"x")
A    <- function(t,x)  eval(aa)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
        0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X      
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),pd,T=T,t0=t0,x0=x0,theta=theta,mu=mu,a=a,b=b,c=c)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("Pearson Diffusions Process",adj=0.5,line=2.5,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]== -.(theta)*(X[t]-.(mu))*dt+sqrt(.(2*theta)*(.(a)*X[t]^2+.(b)*X[t]+.(c)))*dW[t] ),cex=1,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 29/08/2010 21:18:43"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "PDP.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
    }

.PEABM <-
function (X, delta,starts = list(theta = 1 , sigma = 1), leve = 0.95) 
        {
if (length(dim(X)) > 0)
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))
 
CDABM <- function(x,t,x0,theta,sigma,log = FALSE )
         {
ml <- x0 + theta * t
sl <- sigma^2 * t
dnorm (x, mean = ml , sd = sqrt(sl) , log = log  )
      }
ABMlik <-function(theta,sigma) 
       {
n  <- length(X)
dt <- deltat(X)
-sum(CDABM(x=X[2:n],t=dt,x0=X[1:(n-1)],theta,sigma,log = TRUE ))
       }
X <- ts(X,start=0,deltat=delta)
res <- mle(ABMlik, start = starts, method = "L-BFGS-B", lower = c(0, 
          0))
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.PEBS <-
function (X, delta,starts = list(theta = 1, sigma = 1), leve = 0.95) 
        {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

CDBS <- function(x,t,x0,theta,sigma,log = TRUE)
         {
ml <- log(x0) + (theta - ((sigma)^2 /2))*t
sl <- sqrt(t)*sigma
lik <- dlnorm (x, meanlog = ml , sdlog = sl , log = TRUE )
if(!log )
lik <- exp(lik)
lik
      }
BSlik <-function(theta,sigma) 
       {
n  <- length(X)
dt <- deltat(X)
-sum(CDBS(x=X[2:n],t=dt,x0=X[1:(n-1)],theta,sigma,log = TRUE ))
       }
X <- ts(X,start=0,deltat=delta)
res <- mle(BSlik, start = starts, method = "L-BFGS-B", lower = c(0.01, 
          0.01))
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.PEOU <-
function(X,delta,starts=list(r = 1,sigma = 1),leve=0.95)
     {
if (length(dim(X)) > 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

DCOU <-function(x,t,x0,r,sigma,log = FALSE )
      {
Ex <- x0 * exp(-r*t)
Vx <- ((sigma^2)/(2*r))*(1-exp(-2*r*t))
dnorm (x, mean =Ex , sd = sqrt(Vx), log = log )
      }
OUlik <-function ( r , sigma )
       {
n <- length (X)
dt <- deltat (X)
-sum ( DCOU (X [2: n], dt , X [1:(n -1)] , r , sigma , log = TRUE ))
       }
X <- ts(X,start=0,deltat=delta)
res <- mle(OUlik,start = starts,method = "L-BFGS-B",
           lower =c (0 ,0))
{return(print(list(summary=summary(res),coef=coef(res),AIC=AIC(res),vcov=vcov(res),confint=confint(res,level=leve))))}
}

.PEOUexp <-
function(X,delta) 
{
    if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))
    X <- ts(X,start=0,deltat=delta)
    N <- length(X)
    sumtmp <- sum(X[1:(N - 1)] * X[2:N])
    dt <- deltat(X)
    r <- ifelse(sumtmp > 0, -log(sumtmp/sum(X[1:(N - 1)]^2))/dt, 
        NA)
    sigma2 <- 2 * r/((N - 1) * (1 - exp(-2 * dt * r))) * sum((X[2:N] - 
        X[1:(N - 1)] * exp(-dt * r))^2)
    sigma <- sqrt(sigma2)
    {
        return(print(list(r = r, sigma = sigma)))
    }
}

.PEOUG <-
function(X,delta,starts=list(r = 1,theta = 1,sigma = 1),leve=0.95)
     {
if (length(dim(X)) > 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

DCOUG <-function(x,t,x0,r,theta,sigma,log = FALSE )
      {
Ex <- theta+(x0-theta)*exp(-r*t)
Vx <- (sigma^2/(2*r))*(1-exp(-2*r*t))
dnorm (x, mean =Ex , sd = sqrt(Vx), log = log )
      }
OUGlik <-function ( r,theta,sigma )
       {
n <- length (X)
dt <- deltat (X)
-sum ( DCOUG(X[2:n],dt,X[1:(n-1)],r,theta,sigma,log = TRUE ))
       }
X <- ts(X,start=0,deltat=delta)
res <- mle(OUGlik,start = starts,method = "L-BFGS-B",
           lower =c (0,0,0))
{return(print(list(summary=summary(res),coef=coef(res),AIC=AIC(res),vcov=vcov(res),confint=confint(res,level=leve))))}
}

.RK3 <-
function(N,M,T=1,t0,x0,Dt,a,sigma,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(a))
            stop(tkmessageBox(title="Error",message=paste( "The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(sigma))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

rk3 <- function(N,T=1,Dt,t0,x0,a,sigma)
   {
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(sigma)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
              Y    <- numeric()
              Z    <- numeric()
              Y[i-1]=X[i-1]+0.5*Dt*A(t[i-1],X[i-1])+S(t[i-1],X[i-1])*D[i-1]
              Z[i-1]=X[i-1]-A(t[i-1],X[i-1])*Dt+2*Dt*A(t[i-1]+0.5*Dt,Y[i-1])+
                     (2*S(t[i-1]+0.5*Dt,Y[i-1])-S(t[i-1],X[i-1]))*D[i-1]
              X[i] = X[i-1]+(Dt/6)*(A(t[i-1],X[i-1])+4*A(t[i-1]+0.5*Dt,Y[i-1])+A(t[i-1]+Dt,Z[i-1]))+
                     (1/6)*(S(t[i-1],X[i-1])+4*S(t[i-1]+0.5*Dt,Y[i-1])+S(t[i-1]+Dt,Z[i-1]))*D[i-1]
               }          
X      
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),rk3,t0=t0,T=T,x0=x0,Dt=Dt,a=a,sigma=sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Runge-Kutta scheme Order3":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria",date()),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.xlsx(Result, "RK3.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
    }

.ROU <-
function(N,M,t0,T,x0,theta,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( x0 == 0) 
            stop(tkmessageBox(title="Error",message=paste( "x0 =! 0" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

RO <- function(N,T,t0,x0,theta)
   {
Dt <- (T-t0)/N
a <- expression(theta*(x^-1)-x)
s <- expression(1)
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
        0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)
              }
X      
   }
t = seq(t0,T,length=N+1)
Dt <- (T-t0)/N
Q = sapply(rep(N,length=M),RO,T=T,t0=t0,x0=x0,theta=theta)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("Radial Ornstein-Uhlenbeck process",adj=0.5,line=2.5,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(bquote( dX[t]==(.(theta)*X[t]^-1 -X[t] )*dt+dW[t] ),cex=1,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 29/08/2010 20:22:32"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "ROU.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
    }

.snssde <-
function(N,M,T=1,t0,x0,Dt,drift,diffusion,Output=c(FALSE,TRUE),Methods=c("SchEuler","SchMilstein",
                   "SchMilsteinS","SchTaylor","SchHeun",
                   "SchRK3"),...)
        {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(drift))
            stop(tkmessageBox(title="Error",message=paste( "The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(diffusion))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

Methods <- match.arg(Methods)

if ( Methods=="SchEuler" & Output==TRUE)      {R <- .Euler(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=TRUE)}
if ( Methods=="SchEuler" & Output==FALSE)     {R <- .Euler(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=FALSE)}
if ( Methods=="SchMilstein" & Output==TRUE)   {R <- .Milstein(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=TRUE)}
if ( Methods=="SchMilstein" & Output==FALSE)  {R <- .Milstein(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=FALSE)}
if ( Methods=="SchMilsteinS" & Output==TRUE ) {R <- .MilsteinS(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=TRUE)}
if ( Methods=="SchMilsteinS" & Output==FALSE ){R <- .MilsteinS(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=FALSE)}
if ( Methods=="SchTaylor" & Output==TRUE)     {R <- .STS(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=TRUE)}
if ( Methods=="SchTaylor" & Output==FALSE)    {R <- .STS(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=FALSE)}
if ( Methods=="SchHeun" & Output==TRUE)       {R <- .Heun(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=TRUE)}
if ( Methods=="SchHeun" & Output==FALSE)      {R <- .Heun(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=FALSE)}
if ( Methods=="SchRK3" & Output==TRUE)        {R <- .RK3(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=TRUE)}
if ( Methods=="SchRK3" & Output==FALSE)       {R <- .RK3(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output=FALSE)}
return(invisible(R))   
      }

.SRW <-
function(N,t0,T,p=0.5,output=FALSE) 
                 {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( p > 1 | p < 0) 
            stop(tkmessageBox(title="Error",message=paste( "p probability of choosing X in [0,1]" ),icon="error"))

temps = seq(t0,T,length=N+1)
delta.temps = T/N
x = runif(N,0,1)
                 for (i in 1:N)
                 {
                 if ( x[i] >= p)
                 x[i] = +1
                 else
                 x[i] = -1
                 }
x = c(0,x)
M = cumsum(x)
plot(temps,M,las=1,lwd=1,type="n",xlab="time",ylab=expression(R[t]))
points(temps,M,type="s",panel.frist=grid(col="gray"))
mtext("Simulation a Random Walk",line=2.5,cex=1.2)
mtext(bquote(P(X[t]==+1)==.(p)),line=0.25,cex=1.2,adj=0,col="red")
mtext(bquote(P(X[t]==-1)==.(1-p)),line=0.25,cex=1.2,adj=0.5,col="red")
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 22:13 15/01/2010"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,M)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "SRW.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
                 }

.Stgamma <-
function(N,t0,T,alpha,beta,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rgamma(u,alpha,beta)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(0,y)*sqrt(dt)
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :Gamma(.(alpha),.(beta))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1,panel.frist=grid(col="gray"))
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 16:11 30/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "Stgamma.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.STS <-
function(N,M,T=1,t0,x0,Dt,a,sigma,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(a))
            stop(tkmessageBox(title="Error",message=paste( "The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(sigma))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

st <- function(N,T=1,Dt,t0,x0,a,sigma)
   {
DAx  <- D(a,"x")
DAxx <- D(D(a,"x"),"x")
DSx  <- D(sigma,"x")
DSxx <- D(D(sigma,"x"),"x")
A    <- function(t,x)  eval(a)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(sigma)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
DZ= rnorm(N,0,sqrt((1/3)*Dt^3))
X    <- numeric()
X[1] <- x0
for (i in 2:(N+1)){
X[i]=X[i-1]+A(t[i-1],X[i-1])*Dt+S(t[i-1],X[i-1])*D[i-1]+
     0.5*S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1]^2)-Dt)+
     Ax(t[i-1],X[i-1])*S(t[i-1],X[i-1])*DZ[i-1]+0.5*(A(t[i-1],X[i-1])*Ax(t[i-1],X[i-1])+
     0.5*(S(t[i-1],X[i-1])^2)*Axx(t[i-1],X[i-1]))*(Dt^2)+(A(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])+
     0.5*(S(t[i-1],X[i-1])^2)*Sxx(t[i-1],X[i-1]))*(D[i-1]*Dt-DZ[i-1])+
     0.5*S(t[i-1],X[i-1])*(S(t[i-1],X[i-1])*Sxx(t[i-1],X[i-1])+
     (Sx(t[i-1],X[i-1])^2))*((1/3)*(D[i-1]^2)-Dt)*D[i-1]
                }            
X      
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),st,t0=t0,T=T,x0=x0,Dt=Dt,a=a,sigma=sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Strong Taylor Scheme Order 1.5":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria",date()),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.xlsx(Result, "STS.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
    }

.Stst <-
function(N,t0,T,n,output=FALSE)
     {
if ( n <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "degrees of freedom n > 0" ),icon="error"))

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

temps <- seq(t0,T,length=N+1)
dt <- T/N
u <- runif(N)
x <- rt(u,n)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procst <- c(0,y)*sqrt(dt)
plot(temps,procst,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :St[.(n)]),font.main=2)
points(temps,procst,type="l",col="black",lwd=1,panel.frist=grid(col="gray"))
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 16:24 30/03/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
X <- procst
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "Stst.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
                }
return(invisible(Result))
}

.Telegproc <-
function(t0,x0=1,T,lambda,output=FALSE )
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( x0 > 1 | x0 < -1) 
            stop(tkmessageBox(title="Error",message=paste( "x0={-1 Or +1}" ),icon="error"))

temps =seq(t0,T,by=0.1)
u   <- runif(length(temps))
tho <- numeric(length(temps))
x   <- numeric(length(temps)-1)
for ( i in 1:(length(temps)-1))
    {
    tho[i]=(-1/lambda)*log(u[i])   
    if ( temps[i] >= sum(tho[i-1]) && temps[i] < sum(tho[i]) )
    x[i] = +1
    else 
    x[i] = -1
    }
x <- c(x0,x)
THO=cumsum(tho)
plot(temps,x,type="n",xlab="time",ylab=bquote("Space States":E=={-1 | +1}),las=1,lwd=2,col="blue",font.main=1,ylim=c(-2,2),main="Realization a Telegraphic Process")
points(temps,x,type="s",lty=2,col="blue",lwd=2,panel.frist=grid(col="gray"))
mtext(bquote(x[0]==.(x0)),adj=0,line=0.25,cex=1,col="red")
mtext(bquote(lambda==.(lambda)),adj=0.25,line=0.25,cex=1,col="red")
for (i in 1:length(x)){lines(c(temps[i],temps[i+1]),c(x[i],x[i]),type="s",lwd=2,col="red")}
states <- x
time <- temps
X_t <- THO
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 23:45 29/01/2010"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Result <- data.frame(states,time,X_t,tho)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "Telegproc.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
         }

.WNG <-
function(N,t0,T,m,sigma2,output = FALSE)
   {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( sigma2 <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma2 >= 0" ),icon="error"))

temps = seq(t0,T,length=N+1)
bbG <- rnorm(N+1,mean=m,sd=sqrt(sigma2))
par(mfrow=c(2,2))
plot(temps,bbG,type="l",las=1,ylab=expression(epsilon[t]),xlab="time"
     ,cex.lab=1.3,main=paste("white Noise G(m=",m,",var=",sigma2,")" ) )
points(temps,bbG,type="n",panel.frist=grid(col="gray"))
acf(bbG,lag.max=N/5,plot=TRUE,main ="Autocovariance Function",xlab ="Decalage temporel",
    ylab="ACF",type = c("covariance"),panel.frist=grid(col="gray"))
spectrum(bbG,method=c("ar"),las=1,lwd=2,main="Spectral Density Estimation")
mtext("from AR Fit",col="red")
points(bbG,type="n",panel.frist=grid(col="gray"))
spectrum(bbG,method=c("pgram"),las=1,lwd=1,main="Spectral Density Estimation")
mtext("by a Smoothed Periodogram ",col="red")
points(bbG,type="n",panel.frist=grid(col="gray"))
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria 20:51 23/03/2010"),
      side = 4, line = 1, adj = 0, cex = .66)
time <- temps
WNG <- bbG
Result <- data.frame(time,WNG)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "WNG.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
}

.PredCorr <-
function(N,M,T=1,t0,x0,Dt,alpha=0.5,mu=0.5,drift,diffusion,output)
       {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(drift))
            stop(tkmessageBox(title="Error",message=paste( "The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(diffusion))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

if ( alpha > 1 || alpha < 0 )
            stop(tkmessageBox(title="Error",message=paste( "0 <= alpha <= 1" ),icon="error"))

if ( mu > 1 || mu < 0 )
            stop(tkmessageBox(title="Error",message=paste( "0 <= mu <= 1" ),icon="error"))


PC <- function(N,T=1,Dt,t0,x0,alpha=0.5,mu=0.5,drift,diffusion)
     {
DSx  <- D(diffusion,"x")
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diffusion)
Sx   <- function(t,x)  eval(DSx)
SS   <- function(t,x)  eval(drift) - mu * eval(diffusion) * eval(DSx)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt(Dt)
D    <- diff(w)
X    <- numeric()
Y    <- numeric()
X[1] <- x0
Y[1] <- x0
for (i in 2:(N+1)){Y[i] = Y[i-1] + A(t[i-1],Y[i-1])*Dt + S(t[i-1],Y[i-1])*D[i-1]} 
for (i in 2:(N+1)){X[i] = X[i-1] +(alpha*SS(t[i],Y[i])+(1-alpha)*SS(t[i-1],X[i-1]))*Dt+
                          (mu*S(t[i],Y[i])+(1-mu)*S(t[i-1],X[i-1]))*D[i-1]}
X
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),PC,t0=t0,T=T,x0=x0,Dt=Dt,alpha=alpha,mu=mu,drift=drift,diffusion=diffusion)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Predictor-Corrector Method":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1,panel.frist=grid(col="gray"))}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diffusion,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(alpha==.(alpha)),line=0.2,cex=0.9,adj=0.50,col="blue")
mtext(bquote(mu==.(mu)),line=1,cex=0.9,adj=0.50,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="blue")
legend("topleft",bg="gray85",border="gray",c("Average trajectory"),lty=c(1),col=c("blue"),lwd=2)}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria",date()),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.xlsx(Result, "PredCorr.xlsx", sheetName="Sheet 1", formatTemplate=NULL,
           col.names=TRUE, row.names=FALSE, append=FALSE)
}
return(invisible(Result))
}

