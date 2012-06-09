# Mon May 28 19:06:18 2012

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

.ABM <-
function(N,t0,T,x0,theta,sigma,output=FALSE)
         {

if ( N <= 0 ) 
                       stop(tkmessageBox(title="Error",message=paste( "size of process : N >>> 0" ),icon="error"))

if ( t0 >= T || t0 < 0 ) 
                       stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if (sigma < 0 ) 
                       stop(tkmessageBox(title="Error",message=paste( "constant positive : sigma > 0" ),icon="error"))

temps = seq(t0,T,length=N+1)
dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(dt))))
X <- vector()
X[1] <- x0
for (i in 1:N){X[i+1] <- X[i]+ theta*dt + sigma*(w[i+1]-w[i])}
plot(temps,X,las=1,type="n",xlab="time",ylab=expression(X[t]))
points(temps,X,type="l",col="black",lwd=1)
mtext("Arithmetic Brownian Motion",line=2,cex=1.2)
mtext(bquote(dX[t]==.(theta)*dt+.(sigma)*dW[t]),line=0.4,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(dt)),line=0.2,cex=1,adj=0,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
ABM  <- X
Result <- data.frame(time,ABM)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "ABM.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(dt))))
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
for (i in 1:M){points(temps,Q[,i],type="l")}
if (M >=2) {lines(temps,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
ABMF <- Q
time <- temps
X.mean <- Q.mean
Result <- data.frame(time,ABMF)
if (M >=2) {Result <- data.frame(time,ABMF,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "ABMF.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
points(temps,f2,type="l",col="blue",lwd=2)
mtext("Evolution a telegraphic process in time",line=2,cex=1.2)
mtext(bquote(mu==.(mu)),adj=0,line=0.25,cex=1,col="red")
mtext(bquote(lambda==.(lambda)),adj=0.25,line=0.25,cex=1,col="blue")
abline(h=mu/(mu+lambda),lwd=2,col="gray50",lty=2)
axis(2,at=round(mu/(mu+lambda),2),las=1,col.axis="gray50")
text(T/4,(mu/(mu+lambda))+0.2,c(expression(pi[0]==(list(1,0)))),cex=0.8,col="red")
text(T/4,(mu/(mu+lambda))-0.2,c(expression(pi[0]==(list(0,1)))),cex=0.8,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Result <- data.frame(P_t)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
attach(Result)
        }

.BB <-
function(N,t0,T,x0,y,output=FALSE)
    {
if ( t0 >= T || t0 < 0 ) 
                       stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

temps = seq(t0,T,length=N+1)
dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(dt))))
X <- x0 + w - (temps-t0)/(T-t0) * (w[N+1]-y +x0)
plot(temps,X,las=1,type="n",xlab="time",ylab=expression(X[t]))
points(temps,X,type="l",col="black",lwd=1)
mtext("Brownian Bridge",line=2,cex=1.2)
mtext(bquote(x[.(0)]==.(x0)),line=0.15,cex=1.2,adj=0,col="red")
mtext(bquote(y==.(y)),line=0.1,cex=1.2,adj=0.2,col="red")
mtext(bquote(Delta*t==.(dt)),line=0.3,cex=1,adj=1,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
BB   <- X
Result <- data.frame(time,BB)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BB.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(dt))))
X <- x0 + w - (temps-t0)/(T-t0) * (w[N+1]-y +x0)
}
Q = sapply(rep(N,length=M),BB,t0=t0,T=T,x0=x0,y=y)
temps = seq(t0,T,length=N+1)
dt = (T-t0)/N
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(temps,Q[,1],type="n",las=1,xlab="time",ylim=c(r1,r2),ylab=expression(X[t]),cex.lab=1)
for (i in 1:M){points(temps,Q[,i],type="l")}
mtext("flow of the Brownian bridge",line=2.5,cex=1.2)
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1.2,adj=0,col="red")
mtext(bquote(y==.(y)),line=0.1,cex=1.2,adj=0.2,col="red")
mtext(bquote(Delta*t==.(dt)),line=0.4,cex=1.2,adj=0.4,col="red")
if ( M >= 2 ) {lines(temps,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2,cex=1.1)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
X.mean <- Q.mean
BBF <- Q
Result <- data.frame(time,BBF)
if (M >=2) {Result <- data.frame(time,BBF,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BBF.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==frac(.(alpha-1),2*X[t])*dt+dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "Bessel.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
       }

.BMinf <-
function(N,T)
      {
temps = seq(0,T,length=N+1)
Dt = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
limB1 = w/temps
plot(temps,w,las=1,lwd=1,type="l",xlab="t",ylab=expression(W[t]))
points(temps,w,type="n")
points(temps[-1],limB1[-1],type="l",lwd=2,col="green")
mtext("Standard Brownian Motion has the infinite",line=2,cex=1.2)
legend("topleft",border="gray",c(expression(lim(frac(w[t],t),t%->%+infinity)%~~%0)),lty=c(1),col=c("green"),lwd=3)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
}

.BMIrt <-
function(N,T)
      {
temps = seq(0,T,length=N+1)
Dt = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
i = seq(1,N+1,1)
x =  w[N-i+2]-w[N+1]
r1 = min(min(w),min(x))
r2 = max(max(w),max(x))
plot(temps,w,type="l",ylim=c(r1,r2),col="black",las=1,xlab="time",ylab="B(t)")
points(temps,x,col="red",type="l")
mtext("Brownian Motion invariance by reversal of time",line=2,cex=1.2)
legend("topleft",border="gray",c("B(t)","B(t)=B(T-t)-B(T)"),lty=c(1,1),col=c("black","red"),lwd=1)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
}

.BMIto1 <-
function(N,T,output=FALSE)
       {
temps = seq(0,T,length=N+1)
Dt = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Ito <- 0.5*w^2 - 0.5 * temps
Ito2 <- vector()
for (i in 1:(N+1)){ Ito2[i] <- w[i]*(w[i+1]-w[i])}
r1= max(Ito)
r2= min(Ito)
Ito.sum <- cumsum(Ito2)
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Ito.sum,type="l",col="red")
mtext(c((expression("Stochastic Integral":I(w[t])==integral(W[s] * dW[s], 0, t)))),adj=0.5,cex=1.2)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=0,col="red")
legend("topleft",border="gray",c(expression(frac(1,2)*(w[t]^2-t)),expression(sum(w[t[i]]*(w[t[i+1]]-w[t[i]]),i=0,n))),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Ito,Ito.sum)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMIto1.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.BMIto2 <-
function(N,T,output=FALSE)
{
temps = seq(0,T,length=N+1)
Dt = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Ito <- 0.5*w^2 - 0.5 * temps
Ito.sum <- c(0, sapply (2:(N+1), function (x) {w[x] * (w[x+1]-w[x])} ) )
r1= max(Ito)
r2= min(Ito)
Ito.sum <- cumsum(Ito.sum)
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Ito.sum,type="l",col="red")
mtext(c((expression("Stochastic integral":I(w[t])==integral(W[s] * dW[s], 0, t)))),adj=0.5,cex=1.2)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=0,col="red")
legend("topleft",border="gray",c(expression(frac(1,2)*(w[t]^2-t)),expression(sum(w[t[i]]*(w[t[i+1]]-w[t[i]]),i=0,n))),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Ito,Ito.sum)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMIto2.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.BMItoC <-
function(N,T,alpha,output=FALSE)
        {
if ( alpha == 0 )
            stop(tkmessageBox(title="Error",message=paste( "alpha =! 0" ),icon="error"))

temps = seq(0,T,length=N+1)
Dt = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Ito <- alpha*w
Ito.sum <- sapply (1:(N+1), function (x) {alpha * (w[x+1]-w[x])} ) 
r1= min(min(Ito,na.rm=T),min(Ito.sum,na.rm=T))
r2= max(max(Ito,na.rm=T),max(Ito.sum,na.rm=T))
plot(temps,Ito,type="l",las=1,col="blue",ylim=c(r1,r2),xlab="time",ylab=expression(I(w[t])),main=bquote("Stochastic Integral":I(w[t])==alpha*integral(dW[s], 0, t)),cex.lab=1.1)
points(temps,cumsum(Ito.sum),type="l",col="red")
mtext(bquote(alpha==.(alpha)),line=0.25,cex=1.2,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.25,cex=1,adj=0,col="red")
legend("topleft",border="gray",c(bquote(alpha*w[t]),expression(sum(alpha*(w[t[i+1]]-w[t[i]]),i=0,n))),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Ito.sum <- cumsum(Ito.sum)
time <- temps
Result <- data.frame(time,Ito,Ito.sum)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMItoC.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.BMItoP <-
function(N,T,power,output=FALSE)
       {
if ( power  <= 0)  
            stop(tkmessageBox(title="Error",message=paste( "power  > 0 " ),icon="error"))

temps = seq(0,T,length=N+1)
delta.temps = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
Ito <- (1/(power+1))*w^(power +1) - (power /2) * cumsum(w^(power -1)*delta.temps)
Ito.sum <- sapply (1:(N+1), function (x) {w[x]^power * (w[x+1]-w[x])} )
r1= min(min(Ito,na.rm=T),min(Ito.sum,na.rm=T))
r2= max(max(Ito,na.rm=T),max(Ito.sum,na.rm=T))
plot(temps,Ito,type="l",las=1,col="blue",ylim=c(r1,r2),ylab=expression(I(w[t])),xlab="time",main=bquote("Stochastic Integral":I(w[t])==integral(W[s]^n * dW[s], 0, t)),cex.lab=1.1)
points(temps,cumsum(Ito.sum),type="l",col="red")
mtext(bquote(n==.(power)),line=0.25,cex=1.2,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=0,col="red")
legend("topleft",border="gray",c(expression(frac(w[t]^(n+1),n+1)-frac(n,2)*integral(W[s]^(n-1) * ds, 0, t)),expression(sum(w[t[i]]^n*(w[t[i+1]]-w[t[i]]),i=0,N))),
      lty=c(1,1),col=c("blue","red"),lwd=2,cex=0.85)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Ito.sum <- cumsum(Ito.sum)
time <- temps
Result <- data.frame(time,Ito,Ito.sum)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMItoP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.BMItoT <-
function(N,T,output=FALSE)
       {
temps = seq(0,T,length=N+1)
Dt = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Ito <- temps*w - cumsum(w*Dt)
Ito.sum <- sapply (1:(N+1), function (x) {temps[x]*(w[x+1]-w[x])} )
r1= min(min(Ito,na.rm=T),min(Ito.sum,na.rm=T))
r2= max(max(Ito,na.rm=T),max(Ito.sum,na.rm=T))
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",main=bquote("Stochastic Integral":I(w[t])==integral(s*dW[s], 0, t)),cex.lab=1.1)
points(temps,cumsum(Ito.sum),type="l",col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.25,cex=1,adj=0,col="red")
legend("topleft",border="gray",c(expression(t*w[t]-integral(w[s]*ds, 0, t)),expression(sum(t[i]*(w[t[i+1]]-w[t[i]]),i=0,n))),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Ito.sum <- cumsum(Ito.sum)
time <- temps
Result <- data.frame(time,Ito,Ito.sum)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMItoT.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
points(temps,MB,type="n")
mtext("Brownian Motion",line=2,cex=1.2)
mtext("by normal law",line=0.15,cex=1.2,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=1,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- MB
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMN.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
for ( i in 1:M){points(temps,Q[,i],type="l",lwd=1)}
if(M > 1) {lines(temps,Q.mean,lwd=2,col="red")}
if(M > 1) {legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q
time <- temps
X.mean <- Q.mean
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,X,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMNF.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(C*delta.temps))))
      }
Q = sapply(rep(N,length=M),MB,T=T,C=C)
temps = seq(0,T,length=N+1)
f1 <- function(t) 2*sqrt(C*t)
plot(temps,f1(temps),ylim=c(min(-f1(temps)),max(f1(temps))),type="l",las=1,lwd=4,col="red",xlab="time",ylab=expression(W[t]),main=bquote("Trajectories Brownian in the curves"%+-%2*sqrt(C*t)))
points(temps,-f1(temps),type="l",col="red",lwd=4)
mtext(bquote(C==.(C)),line=0.25,cex=1.2,adj=1,col="red")
mtext(bquote("Numbers of the trajectories"==.(M)),line=0.25,cex=1,adj=0,col="blue")
for ( i in 1:M) { points(temps,Q[,i],type="l",lwd=1)}
legend("topleft",border="gray",c(expression(""%+-%2*sqrt(C*t))),lty=c(1),col=c("red"),lwd=2,cex=1.2)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
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
points(temps,w,type="l",col="black",lwd=1)
mtext("Brownian Motion",line=2,cex=1.2)
mtext("by a Random Walk",line=0.15,cex=1.2,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=1,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- w
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMRW.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
for ( i in 1:M){points(temps,Q[,i],type="l",lwd=1)}
if(M > 1) {lines(temps,Q.mean,lwd=2,col="red")}
if(M > 1) {legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q
time <- temps
X.mean <- Q.mean
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,X,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMRWF.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.BMscal <-
function(N,T,S1,S2,S3,output=FALSE)
      {
temps1 <- (S1^2)*seq(0,T,length=N+1)
temps2 <- (S2^2)*seq(0,T,length=N+1)
temps3 <- (S3^2)*seq(0,T,length=N+1)
r <- max(max(temps1),max(temps2),max(temps3))
delta.temps = T/N
w = (1/S1)*c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
w1 = (1/S1)* w
w2 = (1/S2)* w
w3 = (1/S3)* w
r1 <- min(min(w1),min(w2),min(w3))
r2 <- max(max(w1),max(w2),max(w3))
plot(temps1,w1,type="l",ylim=c(r1,r2),xlim=c(0,r),col="black",las=1,xlab="time",ylab=expression(W[t]))
points(temps2,w2,col="red",type="l")
points(temps3,w3,col="blue",type="l")
mtext("Brownian Motion with different scales",line=2,cex=1.2)
legend("topleft",border="gray",
c(paste("S1=",S1),paste("S2=",S2),paste("S3=",S3)),lty=c(1,1,1),col=c("black","red","blue"),lwd=1)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Result <- data.frame(w1,w2,w3)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMscal.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.BMStra <-
function(N,T,output=FALSE)
        {
temps = seq(0,T,length=N+1)
delta.temps = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
Stra <- 0.5 * w^2
r1= max(Stra)
r2= min(Stra)
plot(temps,Stra,type="n",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red")
mtext(c((expression("Stratonovitch integral":I(w[t])==integral(W[s]*o*dW[s], 0, t)))),adj=0.5,cex=1.2)
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=0,col="red")
legend("topleft",border="gray",c(expression(integral(W[s]*o*dW[s], 0, t)==frac(1,2)*W[t]^2)),
      lty=c(1),col=c("red"),lwd=2)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMStra.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.BMStraC <-
function(N,T,alpha,output=FALSE)
        {

if ( alpha == 0 )
            stop(tkmessageBox(title="Error",message=paste( "alpha =! 0" ),icon="error"))

temps = seq(0,T,length=N+1)
delta.temps = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
Stra <- alpha * w
r1= max(Stra)
r2= min(Stra)
plot(temps,Stra,type="n",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red")
mtext(c((expression("Stratonovitch integral":I(w[t])==integral(alpha*o*dW[s], 0, t)))),adj=0.5,cex=1.2)
legend("topleft",border="gray",c(expression(integral(alpha*o*dW[s], 0, t)==alpha*W[t])),
      lty=c(1),col=c("red"),lwd=2)
mtext(bquote(alpha == .(alpha)), line = 0.25, cex = 1.2,adj=0, col = "red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=1,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMStraC.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.BMStraP <-
function(N,T,power,output=FALSE)
        {

if ( power  <= 0)  
            stop(tkmessageBox(title="Error",message=paste( "power > 0" ),icon="error"))

temps = seq(0,T,length=N+1)
delta.temps = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
n <- power
Stra <- numeric(N+1)
for ( i in 1:N){
   Stra[i] <- sum(0.5*(w[i]^(n-1)+w[i+1]^(n-1))*(w[i+1]^2 - w[i]^2))
}
Stra <- cumsum(Stra)
r1= max(Stra)
r2= min(Stra)
plot(temps,Stra,type="n",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red")
mtext(c((expression("Stratonovitch integral":I(w[t])==integral(W[s]^n*o*dW[s], 0, t)))),adj=0.5,cex=1.2)
legend("topleft",border="gray",c(expression(integral(W[s]^n*o*dW[s], 0, t))),
      lty=c(1),col=c("red"),lwd=2)
mtext(bquote(n==.(power)),line=0.25,cex=1.2,col="red",adj=0)
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=1,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMStraP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.BMStraT <-
function(N,T,output=FALSE)
        {
temps = seq(0,T,length=N+1)
delta.temps = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
Stra <- numeric(N+1)
for ( i in 1:N){
   Stra[i] <- sum(0.5*(temps[i]*(w[i+1]-w[i])+temps[i+1]*(w[i+1]-w[i])))
}
Stra <- cumsum(Stra)
r1= max(Stra)
r2= min(Stra)
plot(temps,Stra,type="n",las=1,col="blue",ylab=expression(I(w[t])),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red")
mtext(c((expression("Stratonovitch integral":I(w[t])==integral(s*o*dW[s], 0, t)))),adj=0.5,cex=1.2)
mtext(bquote(Delta*t==.(delta.temps)),line=0.25,cex=1,adj=0,col="red")
legend("topleft",border="gray",c(expression(integral(s*o*dW[s], 0, t))),
      lty=c(1),col=c("red"),lwd=2)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "BMStraT.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==.(mu)*X[t]*dt+.(sigma)*X[t]^.(gamma)*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "CEV.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==(.(r)-.(theta)*X[t])*dt+.(sigma)*sqrt(X[t])*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "CIR.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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

if ( r + ((sigma^2)/2) <= 0)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==-.(r)*X[t]*dt+.(sigma)*sqrt(1+X[t]^2)*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "CIRhy.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)  
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==(.(r)+.(theta)*X[t])*dt+.(sigma)*X[t]^.(gamma)*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "CKLS.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
points(time(X),X,type="l",col="black",lwd=1)
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
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- time(X)
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "diffBridge.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
if ( donc <- TRUE ) {attach(Result)}
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==(X[t]-X[t]^3)*dt+dW[t] ),cex=1,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "DWP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "Euler.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(dt))))
X <- length(temps)
X <- x0*exp((theta - 0.5*sigma^2)*temps + sigma*w)
plot(temps,X,las=1,type="n",xlab="time",ylab=expression(X[t]))
points(temps,X,type="l",col="black",lwd=1)
mtext("Brownian Motion Geometrical (Model of Black-Scholes)" ,line=2.5,cex=1.2)
mtext(bquote(dX[t]==.(theta)*X[t]*dt+.(sigma)*X[t]*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(dt)),line=0.4,cex=1,adj=0,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "GMB.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(dt))))
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
for (i in 1:M){points(temps,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(temps,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q
time <- temps
X.mean <- Q.mean
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,X,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "GMBF.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "Heun.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
Ito.sum <- c(0,sapply(1:(N+1),function(x){exp(-r*(temps[x+1]-temps[x]))*(w[x+1]-w[x])}))
X <- sapply(1:(N+1),function(x){theta+(x0-theta)*exp(-r*temps[x])+sigma*sum(Ito.sum[1:x])})
plot(temps,X,type="n",las=1,xlab="time",ylab=expression(X[t]))
points(temps,X,type="l")
mtext("Hull-White/Vasicek (HWV) Gaussian Diffusion Models",cex=1.2,line=2.5)
mtext(bquote(dX[t]==.(r)*(.(theta)-X[t])*dt+.(sigma)*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "HWV.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
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
for (i in 1:M){points(temps,Q[,i],type="l")}
mtext("Gaussian Diffusion Models(Hull-White/Vasicek)",line=2.5,cex=1.2)
mtext(bquote(dX[t]==.(r)(.(theta)-X[t])*dt+.(sigma)*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0,col="red")
if ( M >= 2 ) {lines(temps,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "HWVF.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==- frac(.(theta)*X[t],sqrt(1+X[t]^2))*dt+dW[t] ),cex=1,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "Hyproc.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==frac(.(sigma)^2,2)* (.(beta)-frac(.(gamma)*X[t],sqrt(.(theta)^2+(X[t]-.(mu))^2) )) *dt+.(sigma)*dW[t] ),cex=1,adj=0.5,line=0.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "Hyprocg.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==X[t]*(.(theta)-(.(theta*r))*X[t])*dt+.(sigma)*X[t]^frac(3,2)*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "INFSR.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
    }

.ItovsStra <-
function(N,T,output=FALSE)
          {
temps = seq(0,T,length=N+1)
delta.temps = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
Stra <- 0.5 * w^2
Ito <- 0.5*w^2 - 0.5 * temps
r1= max(max(Ito),max(Stra))
r2= min(min(Ito),min(Stra))
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),ylim=c(r2,r1),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red")
mtext("Stochastic integral",adj=0.5,cex=1.2,line=3)
mtext(c((expression(I(w[t])==integral(W[s] * dW[s], 0, t)))),adj=0,cex=1,col="blue",line=0.4)
mtext(c((expression(I(w[t])==integral(W[s]*o*dW[s], 0, t)))),adj=1,cex=1,col="red",line=0.4)
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0.5,col="black")
legend("topleft",border="gray",c("Ito Integral","Stratonovitch Integral"),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Ito,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "ItoStra.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.ItovsStraP <-
function(N,T,power,output=FALSE)
        {
if ( power  <= 0)  
            stop(tkmessageBox(title="Error",message=paste( "power  > 0" ),icon="error"))

temps = seq(0,T,length=N+1)
delta.temps = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
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
points(temps,Stra,type="l",col="red")
mtext("Stochastic integral",adj=0.5,cex=1.2,line=3)
mtext(c((expression(I(w[t])==integral(W[s]^n * dW[s], 0, t)))),adj=0,cex=1,col="blue",line=0.4)
mtext(c((expression(I(w[t])==integral(W[s]^n*o*dW[s], 0, t)))),adj=1,cex=1,col="red",line=0.4)
mtext(bquote(n==.(power)),line=1.4,cex=1,col="black",adj=0.5)
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0.5,col="black")
legend("topleft",border="gray",c("Ito Integral","Stratonovitch Integral"),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Ito,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "ItoStraP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.ItovsStraT <-
function(N,T,output=FALSE)
        {
temps = seq(0,T,length=N+1)
delta.temps = T/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
Stra <- numeric(N+1)
for ( i in 1:N){
   Stra[i] <- sum(0.5*(temps[i]*(w[i+1]-w[i])+temps[i+1]*(w[i+1]-w[i])))
}
Stra <- cumsum(Stra)
Ito <- temps*w - cumsum(w*delta.temps)
r1= max(max(Ito),max(Stra))
r2= min(min(Ito),min(Stra))
plot(temps,Ito,type="l",las=1,col="blue",ylab=expression(I(w[t])),ylim=c(r2,r1),xlab="time",cex.lab=1.1)
points(temps,Stra,type="l",col="red")
mtext("Stochastic integral",adj=0.5,cex=1.2,line=3)
mtext(c((expression(I(w[t])==integral(s * dW[s], 0, t)))),adj=0,cex=1,col="blue",line=0.4)
mtext(c((expression(I(w[t])==integral(s*o*dW[s], 0, t)))),adj=1,cex=1,col="red",line=0.4)
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0.5,col="black")
legend("topleft",border="gray",c("Ito Integral","Stratonovitch Integral"),
      lty=c(1,1),col=c("blue","red"),lwd=2)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
time <- temps
Result <- data.frame(time,Ito,Stra)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "ItoStraT.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==-.(theta)*(X[t]-frac(1,2))*dt+sqrt(.(theta)*X[t]*(1-X[t]))*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "JDP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
X = w^2 - temps
Y = exp(sigma*w-0.5*(sigma^2)*temps)
X11()
plot(temps,X,las=1,lwd=1,type="n",xlab="time",ylab=expression(X[t]))
mtext(c(expression(X[t]==W[t]^2 - t )),cex=1.4,col="red",line=1)
points(temps,X,type="l")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
     side = 1, line = 4, adj = 0.5, cex = .66)
X11()
plot(temps,Y,las=1,lwd=1,type="n",xlab="time",ylab=expression(Y[t]))
mtext(c(expression(Y[t]== exp(sigma*W[t] - frac(sigma^2,2)*t) )),cex=1.4,col="red",line=1)
points(temps,Y,type="l")
mtext(bquote(sigma^2==.(sigma)^2),line=0.25,cex=1.2,adj=1,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,X,Y)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "MartExp.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "Milstein.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "MilsteinS.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
Ito.sum <- c(0,sapply(1:(N+1),function(x){exp(-r*(temps[x+1]-temps[x]))*(w[x+1]-w[x])}))
X <- sapply(1:(N+1),function(x){x0*exp(-r*temps[x])+sigma*sum(Ito.sum[1:x])})
plot(temps,X,type="n",las=1,xlab="time",ylab=expression(X[t]))
points(temps,X,type="l")
mtext("Ornstein-Uhlenbeck Process",line=2.5,cex=1.2 )
mtext(bquote(dX[t]==-.(r)*X[t]*dt+.(sigma)*dW[t]),line=0.25,cex=1.2,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(delta.temps)),line=0.4,cex=1,adj=0,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "OU.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(delta.temps))))
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
for (i in 1:M){points(temps,Q[,i],type="l")}
if (M >=2) {lines(temps,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "OUF.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]== -.(theta)*(X[t]-.(mu))*dt+sqrt(.(2*theta)*(.(a)*X[t]^2+.(b)*X[t]+.(c)))*dW[t] ),cex=1,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "PDP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "RK3.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==(.(theta)*X[t]^-1 -X[t] )*dt+dW[t] ),cex=1,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "ROU.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
    }

.snssde <-
function(N,M,T=1,t0,x0,Dt,drift,diffusion,Output=FALSE,Methods=c("SchEuler","SchMilstein",
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

if ( Methods=="SchEuler")     {R <- .Euler(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output)}
if ( Methods=="SchMilstein")  {R <- .Milstein(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output)}
if ( Methods=="SchMilsteinS") {R <- .MilsteinS(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output)}
if ( Methods=="SchTaylor")    {R <- .STS(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output)}
if ( Methods=="SchHeun")      {R <- .Heun(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output)}
if ( Methods=="SchRK3")       {R <- .RK3(N,M,T=1,t0,x0,Dt,a=drift,sigma=diffusion,Output)}
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
points(temps,M,type="s")
mtext("Simulation a Random Walk",line=2.5,cex=1.2)
mtext(bquote(P(X[t]==+1)==.(p)),line=0.25,cex=1.2,adj=0,col="red")
mtext(bquote(P(X[t]==-1)==.(1-p)),line=0.25,cex=1.2,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
time <- temps
Result <- data.frame(time,M)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SRW.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "Stgamma.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
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
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "PredCorr.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(a,adj=0.17,col="red",line=1.8)
mtext(sigma,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "STS.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
points(temps,procst,type="l",col="black",lwd=1)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
X <- procst
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "Stst.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
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
points(temps,x,type="s",lty=2,col="blue",lwd=2)
mtext(bquote(x[0]==.(x0)),adj=0,line=0.25,cex=1,col="red")
mtext(bquote(lambda==.(lambda)),adj=0.25,line=0.25,cex=1,col="red")
for (i in 1:length(x)){lines(c(temps[i],temps[i+1]),c(x[i],x[i]),type="s",lwd=2,col="red")}
states <- x
time <- temps
X_t <- THO
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.6, cex = .66)
Result <- data.frame(states,time,X_t,tho)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "Telegproc.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
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
points(temps,bbG,type="n")
acf(bbG,lag.max=N/5,plot=TRUE,main ="Autocovariance Function",xlab ="Decalage temporel",
    ylab="ACF",type = c("covariance"))
spectrum(bbG,method=c("ar"),las=1,lwd=2,main="Spectral Density Estimation")
mtext("from AR Fit",col="red")
points(bbG,type="n")
spectrum(bbG,method=c("pgram"),las=1,lwd=1,main="Spectral Density Estimation")
mtext("by a Smoothed Periodogram ",col="red")
points(bbG,type="n")
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 4, line = 1, adj = 0, cex = .66)
time <- temps
WNG <- bbG
Result <- data.frame(time,WNG)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "WNG.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}


.Euler2D <-
function(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y)" ),icon="error"))

Ax    <- function(t,x,y)  eval(driftx)
Ay    <- function(t,x,y)  eval(drifty)
Sx    <- function(t,x,y)  eval(diffx)
Sy    <- function(t,x,y)  eval(diffy)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
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
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1] 
                  } 

if(Step==FALSE){
plot(X,Y,type="l",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
mtext(expression("Euler scheme : Simulation SDE Two-Dimensional"),line=3.4,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}

if(Step==TRUE){
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
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2 ,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
time <- t
X1   <- X
X2   <- Y
Result <- data.frame(time,X1,X2)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "Euler2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.Milstein2D <-
function(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y)" ),icon="error"))


DSxx   <- D(diffx,"x")
DSyy   <- D(diffy,"y")
Ax    <- function(t,x,y)  eval(driftx)
Ay    <- function(t,x,y)  eval(drifty)
Sx    <- function(t,x,y)  eval(diffx)
DSx   <- function(t,x,y)  eval(DSxx)
Sy    <- function(t,x,y)  eval(diffy)
DSy   <- function(t,x,y)  eval(DSyy)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
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
				  
if(Step==FALSE){
plot(X,Y,type="l",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
mtext(expression("Milstein scheme : Simulation SDE Two-Dimensional"),line=3.4,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}

if(Step==TRUE){
plot(X,Y,type="n",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="black",lwd=1)}
mtext(expression("Milstein scheme : Simulation SDE Two-Dimensional"),line=3.4,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2 ,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}				  
time <- t
X1   <- X
X2   <- Y
Result <- data.frame(time,X1,X2)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "Milstein2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.MilsteinS2D <-
function(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y)" ),icon="error"))



DDAx   <- D(driftx,"x")
DDDAx  <- D(D(driftx,"x"),"x")

DDSx   <- D(diffx,"x")
DDDSx  <- D(D(diffx,"x"),"x")

Ax     <- function(t,x,y)  eval(driftx)
DAxx   <- function(t,x,y)  eval(DDAx)
DAxxx  <- function(t,x,y)  eval(DDDAx)
Sx     <- function(t,x,y)  eval(diffx)
DSx    <- function(t,x,y)  eval(DDSx)
DSxx   <- function(t,x,y)  eval(DDDSx)

DDAy   <- D(drifty,"y")
DDDAy  <- D(D(drifty,"y"),"y")
DDSy   <- D(diffy,"y")
DDDSy  <- D(D(diffy,"y"),"y")
Ay     <- function(t,x,y)  eval(drifty)
DAyy   <- function(t,x,y)  eval(DDAy)
DAyyy  <- function(t,x,y)  eval(DDDAy)
Sy     <- function(t,x,y)  eval(diffy)
DSy    <- function(t,x,y)  eval(DDSy)
DSyy   <- function(t,x,y)  eval(DDDSy)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
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
X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1] +
       0.5 *Sx(t[i-1],X[i-1],Y[i-1]) * DSx(t[i-1],X[i-1],Y[i-1])*(Dx[i-1]^2-Dt)+ 
       Dt^(3/2)*(0.5 *Ax(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1]) +
       0.5 *DAxx(t[i-1],X[i-1],Y[i-1])*Sx(t[i-1],X[i-1],Y[i-1])+
       0.25 *(Sx(t[i-1],X[i-1],Y[i-1])^2) * DSxx(t[i-1],X[i-1],Y[i-1]))*Dx[i -1]+ 
       (Dt^2) * (0.5*Ax(t[i-1],X[i-1],Y[i-1])*DAxx(t[i-1],X[i-1],Y[i-1])+
       0.25 *DAxxx(t[i-1],X[i-1],Y[i-1])*(Sx(t[i-1],X[i-1],Y[i-1])^2))

Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1] +
       0.5 *Sy(t[i-1],X[i-1],Y[i-1]) * DSy(t[i-1],X[i-1],Y[i-1])*(Dy[i-1]^2-Dt)+ 
       Dt^(3/2)*(0.5 *Ay(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1]) +
       0.5 *DAyy(t[i-1],X[i-1],Y[i-1])*Sy(t[i-1],X[i-1],Y[i-1])+
       0.25 *(Sy(t[i-1],X[i-1],Y[i-1])^2) * DSyy(t[i-1],X[i-1],Y[i-1]))*Dy[i -1]+ 
       (Dt^2) * (0.5*Ay(t[i-1],X[i-1],Y[i-1])*DAyy(t[i-1],X[i-1],Y[i-1])+
       0.25 *DAyyy(t[i-1],X[i-1],Y[i-1])*(Sy(t[i-1],X[i-1],Y[i-1])^2))
                  } 

if(Step==FALSE){
plot(X,Y,type="l",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
mtext(expression("Second Milstein scheme : Simulation SDE Two-Dimensional"),line=3.4,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
plot(X,Y,type="n",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="black",lwd=1)}
mtext(expression("Second Milstein scheme : Simulation SDE Two-Dimensional"),line=3.4,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2 ,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}					  
time <- t
X1   <- X
X2   <- Y
Result <- data.frame(time,X1,X2)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "MilsteinS2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.Heun2D <-
function(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y)" ),icon="error"))

Ax    <- function(t,x,y)  eval(driftx)
Ay    <- function(t,x,y)  eval(drifty)
Sx    <- function(t,x,y)  eval(diffx)
Sy    <- function(t,x,y)  eval(diffy)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X    <- numeric()
Y    <- numeric()
XX   <- numeric()
YY   <- numeric()
X[1] <- x0
Y[1] <- y0
for (i in 2:(N+1)){
         XX[i-1]= X[i-1]+Ax(t[i-1],X[i-1],Y[i-1])*Dt+Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]
         YY[i-1]= Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1])*Dt+Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]
         X[i]   = X[i-1]+0.5*Dt*(Ax(t[i-1],X[i-1],Y[i-1])+Ax(t[i-1],XX[i-1],Y[i-1]))+
                 0.5*(Sx(t[i-1],X[i-1],Y[i-1])+Sx(t[i-1],XX[i-1],Y[i-1]))*Dx[i-1]
         Y[i]   = Y[i-1]+0.5*Dt*(Ay(t[i-1],X[i-1],Y[i-1])+Ay(t[i-1],X[i-1],YY[i-1]))+
                 0.5*(Sy(t[i-1],X[i-1],Y[i-1])+Sy(t[i-1],X[i-1],YY[i-1]))*Dy[i-1]
                  } 
				  
if(Step==FALSE){
plot(X,Y,type="l",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
mtext(expression("Heun scheme : Simulation SDE Two-Dimensional"),line=3.2,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
plot(X,Y,type="n",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="black",lwd=1)}
mtext(expression("Heun scheme : Simulation SDE Two-Dimensional"),line=3.2,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2 ,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}					  
time <- t
X1   <- X
X2   <- Y
Result <- data.frame(time,X1,X2)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "Heun2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.RK32D <-
function(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y)" ),icon="error"))

Ax    <- function(t,x,y)  eval(driftx)
Ay    <- function(t,x,y)  eval(drifty)

Sx    <- function(t,x,y)  eval(diffx)
Sy    <- function(t,x,y)  eval(diffy)


if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X    <- numeric()
Y    <- numeric()
XX   <- numeric()
YY   <- numeric()
XXX  <- numeric()
YYY  <- numeric()
X[1] <- x0
Y[1] <- y0
for (i in 2:(N+1)){
              XX[i-1]=X[i-1]+0.5*Dt*Ax(t[i-1],X[i-1],Y[i-1])+Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]
              YY[i-1]=Y[i-1]+0.5*Dt*Ay(t[i-1],X[i-1],Y[i-1])+Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]
              XXX[i-1]=X[i-1]-Ax(t[i-1],X[i-1],Y[i-1])*Dt+2*Dt*Ax(t[i-1]+0.5*Dt,XX[i-1],Y[i-1])+
                     (2*Sx(t[i-1]+0.5*Dt,XX[i-1],Y[i-1])-Sx(t[i-1],X[i-1],Y[i-1]))*Dx[i-1]
              YYY[i-1]=Y[i-1]-Ay(t[i-1],X[i-1],Y[i-1])*Dt+2*Dt*Ay(t[i-1]+0.5*Dt,X[i-1],YY[i-1])+
                     (2*Sy(t[i-1]+0.5*Dt,X[i-1],YY[i-1])-Sy(t[i-1],X[i-1],Y[i-1]))*Dy[i-1]
              X[i] = X[i-1]+(Dt/6)*(Ax(t[i-1],X[i-1],Y[i-1])+4*Ax(t[i-1]+0.5*Dt,XX[i-1],Y[i-1])+Ax(t[i-1]+Dt,XXX[i-1],Y[i-1]))+
                     (1/6)*(Sx(t[i-1],X[i-1],Y[i-1])+4*Sx(t[i-1]+0.5*Dt,XX[i-1],Y[i-1])+Sx(t[i-1]+Dt,XXX[i-1],Y[i-1]))*Dx[i-1]
              Y[i] = Y[i-1]+(Dt/6)*(Ay(t[i-1],X[i-1],Y[i-1])+4*Ay(t[i-1]+0.5*Dt,X[i-1],YY[i-1])+Ay(t[i-1]+Dt,X[i-1],YYY[i-1]))+
                     (1/6)*(Sy(t[i-1],X[i-1],Y[i-1])+4*Sy(t[i-1]+0.5*Dt,X[i-1],YY[i-1])+Sy(t[i-1]+Dt,X[i-1],YYY[i-1]))*Dy[i-1]
                  } 

if(Step==FALSE){
plot(X,Y,type="l",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
mtext(expression("Runge-Kutta scheme Order3: Simulation SDE Two-Dimensional"),line=3.2,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
plot(X,Y,type="n",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="black",lwd=1)}
mtext(expression("Runge-Kutta scheme Order3: Simulation SDE Two-Dimensional"),line=3.2,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2 ,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}					  
time <- t
X1   <- X
X2   <- Y
Result <- data.frame(time,X1,X2)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "RK32D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.STS2D <-
function(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y)" ),icon="error"))



DDAx   <- D(driftx,"x")
DDDAx  <- D(D(driftx,"x"),"x")

DDSx   <- D(diffx,"x")
DDDSx  <- D(D(diffx,"x"),"x")
Ax     <- function(t,x,y)  eval(driftx)
DAxx   <- function(t,x,y)  eval(DDAx)
DAxxx  <- function(t,x,y)  eval(DDDAx)
Sx     <- function(t,x,y)  eval(diffx)
DSx    <- function(t,x,y)  eval(DDSx)
DSxx   <- function(t,x,y)  eval(DDDSx)

DDAy   <- D(drifty,"y")
DDDAy  <- D(D(drifty,"y"),"y")
DDSy   <- D(diffy,"y")
DDDSy  <- D(D(diffy,"y"),"y")
Ay     <- function(t,x,y)  eval(drifty)
DAyy   <- function(t,x,y)  eval(DDAy)
DAyyy  <- function(t,x,y)  eval(DDDAy)
Sy     <- function(t,x,y)  eval(diffy)
DSy    <- function(t,x,y)  eval(DDSy)
DSyy   <- function(t,x,y)  eval(DDDSy)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
DZx= rnorm(N,0,sqrt((1/3)*Dt^3))
DZy= rnorm(N,0,sqrt((1/3)*Dt^3))
X    <- numeric()
Y    <- numeric()
X[1] <- x0
Y[1] <- y0
for (i in 2:(N+1)){
X[i]=X[i-1]+Ax(t[i-1],X[i-1],Y[i-1])*Dt+Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]+
     0.5*Sx(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1])*((Dx[i-1]^2)-Dt)+
     DAxx(t[i-1],X[i-1],Y[i-1])*Sx(t[i-1],X[i-1],Y[i-1])*DZx[i-1]+0.5*(Ax(t[i-1],X[i-1],Y[i-1])*DAxx(t[i-1],X[i-1],Y[i-1])+
     0.5*(Sx(t[i-1],X[i-1],Y[i-1])^2)*DAxxx(t[i-1],X[i-1],Y[i-1]))*(Dt^2)+(Ax(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1])+
     0.5*(Sx(t[i-1],X[i-1],Y[i-1])^2)*DSxx(t[i-1],X[i-1],Y[i-1]))*(Dx[i-1]*Dt-DZx[i-1])+
     0.5*Sx(t[i-1],X[i-1],Y[i-1])*(Sx(t[i-1],X[i-1],Y[i-1])*DSxx(t[i-1],X[i-1],Y[i-1])+
     (DSx(t[i-1],X[i-1],Y[i-1])^2))*((1/3)*(Dx[i-1]^2)-Dt)*Dx[i-1]
Y[i]=Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1])*Dt+Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]+
     0.5*Sy(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1])*((Dy[i-1]^2)-Dt)+
     DAyy(t[i-1],X[i-1],Y[i-1])*Sy(t[i-1],X[i-1],Y[i-1])*DZy[i-1]+0.5*(Ay(t[i-1],X[i-1],Y[i-1])*DAyy(t[i-1],X[i-1],Y[i-1])+
     0.5*(Sy(t[i-1],X[i-1],Y[i-1])^2)*DAyyy(t[i-1],X[i-1],Y[i-1]))*(Dt^2)+(Ay(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1])+
     0.5*(Sy(t[i-1],X[i-1],Y[i-1])^2)*DSyy(t[i-1],X[i-1],Y[i-1]))*(Dy[i-1]*Dt-DZy[i-1])+
     0.5*Sy(t[i-1],X[i-1],Y[i-1])*(Sy(t[i-1],X[i-1],Y[i-1])*DSyy(t[i-1],X[i-1],Y[i-1])+
     (DSy(t[i-1],X[i-1],Y[i-1])^2))*((1/3)*(Dy[i-1]^2)-Dt)*Dy[i-1]
                  } 
				  
if(Step==FALSE){
plot(X,Y,type="l",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
mtext(expression("Strong Taylor Scheme Order 1.5: Simulation SDE Two-Dimensional"),line=3.2,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
plot(X,Y,type="n",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="black",lwd=1)}
mtext(expression("Strong Taylor Scheme Order 1.5: Simulation SDE Two-Dimensional"),line=3.2,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2 ,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}					  
time <- t
X1   <- X
X2   <- Y
Result <- data.frame(time,X1,X2)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "STS2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.snssde2D <-
function(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step=FALSE,Output=FALSE,Methods=c("SchEuler","SchMilstein",
                   "SchMilsteinS","SchTaylor","SchHeun",
                   "SchRK3"),...)
        {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y)" ),icon="error"))

Methods <- match.arg(Methods)

if ( Methods=="SchEuler" )    {R <- .Euler2D(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step,Output)}
if ( Methods=="SchMilstein")  {R <- .Milstein2D(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step,Output)}
if ( Methods=="SchMilsteinS") {R <- .MilsteinS2D(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step,Output)}
if ( Methods=="SchTaylor")    {R <- .STS2D(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step,Output)}
if ( Methods=="SchHeun")      {R <- .Heun2D(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step,Output)}
if ( Methods=="SchRK3")       {R <- .RK32D(N,T=1,t0,x0,y0,Dt,driftx,drifty,diffx,diffy,Step,Output)}
      }

.PredCorr2D <-
function(N,T=1,t0,x0,y0,Dt,alpha=0.5,mu=0.5,driftx,drifty,diffx,diffy,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if ( alpha > 1 || alpha < 0 )
            stop(tkmessageBox(title="Error",message=paste( "0 <= alpha <= 1" ),icon="error"))

if ( mu > 1 || mu < 0 )
            stop(tkmessageBox(title="Error",message=paste( "0 <= mu <= 1" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y)" ),icon="error"))


DSx  <- D(diffx,"x")
Ax    <- function(t,x,y)  eval(driftx)
Sx    <- function(t,x,y)  eval(diffx)
Sxx   <- function(t,x,y)  eval(DSx)
SSx   <- function(t,x,y)  eval(driftx) - mu * eval(diffx) * eval(DSx)

DSy  <- D(diffy,"y")
Ay    <- function(t,x,y)  eval(drifty)
Sy    <- function(t,x,y)  eval(diffy)
Syy   <- function(t,x,y)  eval(DSy)
SSy   <- function(t,x,y)  eval(drifty) - mu * eval(diffy) * eval(DSy)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X    <- numeric()
Y    <- numeric()
XX    <- numeric()
YY    <- numeric()
X[1] <- XX[1] <- x0
Y[1] <- YY[1] <- y0
for (i in 2:(N+1)){
XX[i] = XX[i-1] + Ax(t[i-1],XX[i-1],YY[i-1])*Dt + Sx(t[i-1],XX[i-1],YY[i-1])*Dx[i-1]
YY[i] = YY[i-1] + Ay(t[i-1],XX[i-1],YY[i-1])*Dt + Sy(t[i-1],XX[i-1],YY[i-1])*Dy[i-1]
}

for (i in 2:(N+1)){
X[i] = X[i-1] +(alpha*SSx(t[i],XX[i],YY[i])+(1-alpha)*SSx(t[i-1],X[i-1],Y[i-1]))*Dt+
       (mu*Sx(t[i],XX[i],YY[i])+(1-mu)*Sx(t[i-1],X[i-1],Y[i-1]))*Dx[i-1]
Y[i] = Y[i-1] +(alpha*SSy(t[i],XX[i],YY[i])+(1-alpha)*SSy(t[i-1],X[i-1],Y[i-1]))*Dt+
       (mu*Sy(t[i],XX[i],YY[i])+(1-mu)*Sy(t[i-1],X[i-1],Y[i-1]))*Dy[i-1]
} 
if(Step==FALSE){
plot(X,Y,type="l",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
mtext(expression("Predictor-Corrector Method: Simulation SDE Two-Dimensional"),line=3.2,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
plot(X,Y,type="n",xlab=expression(X[t]^1),ylab=expression(X[t]^2),las=1)
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="black",lwd=1)}
mtext(expression("Predictor-Corrector Method: Simulation SDE Two-Dimensional"),line=3.2,adj=0.5,cex=1,col="black")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=1.6,col="red3")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="red3")
mtext(bquote(X[t[0]]^1==.(x0)),line=1.6,adj=0.78,cex=1,col="blue")
mtext(bquote(X[t[0]]^2==.(y0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(T==.(T)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=1,col="blue")
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X11(6,6)
par( mar =c(3 ,3 ,2 ,1))
par(mfrow=c(2,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2)*dt + sigma[1](t,X[t]^1,X[t]^2) *d*W[t]^1),cex=1,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2)*dt + sigma[2](t,X[t]^1,X[t]^2) *d*W[t]^2),cex=1,adj=0,line=0.1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria") ,side = 1, line = 4, adj = 0.5, cex = .66)
}					  
time <- t
X1   <- X
X2   <- Y
Result <- data.frame(time,X1,X2)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "PredCorr2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.AnaSimX <-
function(N,M,t0,Dt,T=1,X0,v,drift,diff,Output=FALSE,
                    Methods = c("Euler", "Milstein", "MilsteinS", 
                                "Ito-Taylor", "Heun", "RK3"), ...)
{

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M length of first time passege must be >>> 30" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if ( (N*Dt) < v || t0 >= v )
            stop(tkmessageBox(title="Error",message=paste( " v = k * Dt with k integer, 1 <= k <= N " ),icon="error"))


if(!is.expression(drift))
            stop(tkmessageBox(title="Error",message=paste( "The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(diff))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))


Methods <- match.arg(Methods)

if(Methods=="Euler")
                  {
Eul <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {

A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
X    <- numeric()
X[1] <- X0
for (i in 2:(N+1)){X[i] = X[i-1] + A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]}      
X
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Eul,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Euler scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(v=v,lwd=2,col="red")
axis(1,at=v,las=1,labels=expression(v[t]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(v[t]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q[which(t==v),]
Anay_Euler <- data.frame(X)
showData(Anay_Euler, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_Euler, file = "AnaSimX.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_Euler)
}

if(Methods=="Milstein")
                  {
Mil <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
DSx  <- D(diff,"x")
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
X    <- numeric()
X[1] <- X0
for (i in 2:(N+1)){
            X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
                   0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)}     
X
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Mil,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Milstein scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(v=v,lwd=2,col="red")
axis(1,at=v,las=1,labels=expression(v[t]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(v[t]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q[which(t==v),]
Anay_Milstein <- data.frame(X)
showData(Anay_Milstein, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_Milstein, file = "AnaSimX.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_Milstein)
}

if(Methods=="MilsteinS")
                  {
MilS <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
DAx  <- D(drift,"x")
DAxx <- D(D(drift,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drift)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
X    <- numeric()
X[1] <- X0
for (i in 2:(N+1)){
             X[i] = X[i-1] + A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1] +
                    0.5 *S(t[i-1],X[i-1]) * Sx(t[i-1],X[i-1])*(D[i-1]^2-Dt)+ 
                    Dt^(3 /2)*(0.5 *A(t[i-1],X[i-1])*Sx(t[i-1],X[i-1]) +
                    0.5 *Ax(t[i-1],X[i-1])*S(t[i-1],X[i-1])+
                    0.25 *(S(t[i -1] ,X[i -1])^2) * Sxx(t[i -1] ,X[i -1]))*D[i -1]+ 
                   (Dt^2) * (0.5*A(t[i -1],X[i -1])*Ax(t[i-1],X[i-1])+
                    0.25 *Axx(t[i-1],X[i-1])*(S(t[i-1],X[i-1])^2))}     
X
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),MilS,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Second Milstein scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(v=v,lwd=2,col="red")
axis(1,at=v,las=1,labels=expression(v[t]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(v[t]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q[which(t==v),]
Anay_MilsteinS <- data.frame(X)
showData(Anay_MilsteinS, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_MilsteinS, file = "AnaSimX.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_MilsteinS)
}

if(Methods=="Ito-Taylor")
                  {
ITO <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
DAx  <- D(drift,"x")
DAxx <- D(D(drift,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drift)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
DZ= rnorm(N,0,sqrt((1/3)*Dt^3))
X    <- numeric()
X[1] <- X0
for (i in 2:(N+1)){
                 X[i]=X[i-1]+A(t[i-1],X[i-1])*Dt+S(t[i-1],X[i-1])*D[i-1]+
                      0.5*S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1]^2)-Dt)+
                      Ax(t[i-1],X[i-1])*S(t[i-1],X[i-1])*DZ[i-1]+0.5*(A(t[i-1],X[i-1])*Ax(t[i-1],X[i-1])+
                      0.5*(S(t[i-1],X[i-1])^2)*Axx(t[i-1],X[i-1]))*(Dt^2)+(A(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])+
                      0.5*(S(t[i-1],X[i-1])^2)*Sxx(t[i-1],X[i-1]))*(D[i-1]*Dt-DZ[i-1])+
                      0.5*S(t[i-1],X[i-1])*(S(t[i-1],X[i-1])*Sxx(t[i-1],X[i-1])+
                     (Sx(t[i-1],X[i-1])^2))*((1/3)*(D[i-1]^2)-Dt)*D[i-1]}     
X
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),ITO,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Strong Taylor Scheme Order 1.5":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(v=v,lwd=2,col="red")
axis(1,at=v,las=1,labels=expression(v[t]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(v[t]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q[which(t==v),]
Anay_STS <- data.frame(X)
showData(Anay_STS, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_STS, file = "AnaSimX.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_STS)
}

if(Methods=="Heun")
                  {
Heu <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
X    <- numeric()
X[1] <- X0
for (i in 2:(N+1)){
         Y[i-1]= X[i-1]+A(t[i-1],X[i-1])*Dt+S(t[i-1],X[i-1])*D[i-1]
         X[i]  = X[i-1]+0.5*Dt*(A(t[i-1],X[i-1])+A(t[i-1],Y[i-1]))+
                 0.5*(S(t[i-1],X[i-1])+S(t[i-1],Y[i-1]))*D[i-1]}     
X
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Heu,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Heun scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(v=v,lwd=2,col="red")
axis(1,at=v,las=1,labels=expression(v[t]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(v[t]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q[which(t==v),]
Anay_Heun <- data.frame(X)
showData(Anay_Heun, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_Heun, file = "AnaSimX.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_Heun)
}

if(Methods=="RK3")
                  {
RK <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
X    <- numeric()
X[1] <- X0
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
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),RK,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Runge-Kutta scheme Order3":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(v=v,lwd=2,col="red")
axis(1,at=v,las=1,labels=expression(v[t]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(v[t]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
X <- Q[which(t==v),]
Anay_RK3 <- data.frame(X)
showData(Anay_RK3, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_RK3, file = "AnaSimX.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_RK3)
}
}

.AnaSimFPT <- function(N,M,t0,Dt,T=1,X0,v,drift,diff,ELRENA=c("No","Yes","Mean","Median"),Output=FALSE,
                    Methods = c("Euler", "Milstein", "MilsteinS", 
                                "Ito-Taylor", "Heun", "RK3"), ...)
{

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M length of first time passege must be >>> 30" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))


if(!is.expression(drift))
            stop(tkmessageBox(title="Error",message=paste( "The coefficient of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(diff))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))


Methods <- match.arg(Methods)
ELRENA <- match.arg(ELRENA)

if(Methods=="Euler")
                  {
Eul <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
X    <- numeric()
X[1] <- X0
for (i in 2:(N+1)){X[i] = X[i-1] + A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]}      
X
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Eul,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Euler scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(h=v,lwd=2,col="red")
axis(2,at=v,las=1,labels=expression(X[v]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(X[v]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
Y <- rep(NA,M)
if (X0 > v){
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}else{
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}
Anay_Euler <- data.frame(tau)
showData(Anay_Euler, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_Euler, file = "AnaSimFTP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_Euler)
names(Anay_Euler)
}

if(Methods=="Milstein")
                  {
Mil <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
DSx  <- D(diff,"x")
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
X    <- numeric()
X[1] <- X0
for (i in 2:(N+1)){
            X[i] = X[i-1]+ A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1]+ 
                   0.5 *S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1])^2 -Dt)}     
X
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Mil,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Milstein scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(h=v,lwd=2,col="red")
axis(2,at=v,las=1,labels=expression(X[v]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(X[v]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
Y <- rep(NA,M)
if (X0 > v){
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}else{
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}
Anay_Milstein <- data.frame(tau)
showData(Anay_Milstein, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_Milstein, file = "AnaSimFTP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_Milstein)
names(Anay_Milstein)
}

if(Methods=="MilsteinS")
                  {
MilS <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
DAx  <- D(drift,"x")
DAxx <- D(D(drift,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drift)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
X    <- numeric()
X[1] <- X0
for (i in 2:(N+1)){
             X[i] = X[i-1] + A(t[i-1],X[i-1])*Dt + S(t[i-1],X[i-1])*D[i-1] +
                    0.5 *S(t[i-1],X[i-1]) * Sx(t[i-1],X[i-1])*(D[i-1]^2-Dt)+ 
                    Dt^(3 /2)*(0.5 *A(t[i-1],X[i-1])*Sx(t[i-1],X[i-1]) +
                    0.5 *Ax(t[i-1],X[i-1])*S(t[i-1],X[i-1])+
                    0.25 *(S(t[i -1] ,X[i -1])^2) * Sxx(t[i -1] ,X[i -1]))*D[i -1]+ 
                   (Dt^2) * (0.5*A(t[i -1],X[i -1])*Ax(t[i-1],X[i-1])+
                    0.25 *Axx(t[i-1],X[i-1])*(S(t[i-1],X[i-1])^2))}     
X
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),MilS,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Second Milstein scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(h=v,lwd=2,col="red")
axis(2,at=v,las=1,labels=expression(X[v]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(X[v]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
Y <- rep(NA,M)
if (X0 > v){
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}else{
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}
Anay_MilsteinS <- data.frame(tau)
showData(Anay_MilsteinS, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_MilsteinS, file = "AnaSimFTP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_MilsteinS)
names(Anay_MilsteinS)
}

if(Methods=="Ito-Taylor")
                  {
ITO <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
DAx  <- D(drift,"x")
DAxx <- D(D(drift,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drift)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
DZ= rnorm(N,0,sqrt((1/3)*Dt^3))
X    <- numeric()
X[1] <- X0
for (i in 2:(N+1)){
                 X[i]=X[i-1]+A(t[i-1],X[i-1])*Dt+S(t[i-1],X[i-1])*D[i-1]+
                      0.5*S(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])*((D[i-1]^2)-Dt)+
                      Ax(t[i-1],X[i-1])*S(t[i-1],X[i-1])*DZ[i-1]+0.5*(A(t[i-1],X[i-1])*Ax(t[i-1],X[i-1])+
                      0.5*(S(t[i-1],X[i-1])^2)*Axx(t[i-1],X[i-1]))*(Dt^2)+(A(t[i-1],X[i-1])*Sx(t[i-1],X[i-1])+
                      0.5*(S(t[i-1],X[i-1])^2)*Sxx(t[i-1],X[i-1]))*(D[i-1]*Dt-DZ[i-1])+
                      0.5*S(t[i-1],X[i-1])*(S(t[i-1],X[i-1])*Sxx(t[i-1],X[i-1])+
                     (Sx(t[i-1],X[i-1])^2))*((1/3)*(D[i-1]^2)-Dt)*D[i-1]}     
X
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),ITO,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Strong Taylor Scheme Order 1.5":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(h=v,lwd=2,col="red")
axis(2,at=v,las=1,labels=expression(X[v]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(X[v]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
Y <- rep(NA,M)
if (X0 > v){
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}else{
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}
Anay_STS <- data.frame(tau)
showData(Anay_STS, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_STS, file = "AnaSimFTP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_STS)
names(Anay_STS)
}

if(Methods=="Heun")
                  {
Heu <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
X    <- numeric()
X[1] <- X0
for (i in 2:(N+1)){
         Y[i-1]= X[i-1]+A(t[i-1],X[i-1])*Dt+S(t[i-1],X[i-1])*D[i-1]
         X[i]  = X[i-1]+0.5*Dt*(A(t[i-1],X[i-1])+A(t[i-1],Y[i-1]))+
                 0.5*(S(t[i-1],X[i-1])+S(t[i-1],Y[i-1]))*D[i-1]}     
X
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Heu,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Heun scheme":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(h=v,lwd=2,col="red")
axis(2,at=v,las=1,labels=expression(X[v]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(X[v]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
Y <- rep(NA,M)
if (X0 > v){
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}else{
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}
Anay_Heun <- data.frame(tau)
showData(Anay_Heun, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_Heun, file = "AnaSimFTP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_Heun)
names(Anay_Heun)
}

if(Methods=="RK3")
                  {
RK <- function(N,T=1,Dt,t0,X0,v,drift,diff)
     {
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
X    <- numeric()
X[1] <- X0
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
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),RK,T=T,t0=t0,X0=X0,Dt=Dt,v=v,drift=drift,diff=diff)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(expression("Runge-Kutta scheme Order3":dX[t]== a(t,X[t])*dt+sigma(t,X[t])*dW[t]),line=3,cex=1.2,adj=0.5)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
abline(h=v,lwd=2,col="red")
axis(2,at=v,las=1,labels=expression(X[v]),tick = FALSE,col.axis="red")
mtext(expression(a(t,X[t])==.),adj=0,col="red",line=1.8)
mtext(expression(sigma(t,X[t])==.),adj=0,col="red",line=0.4)
mtext(drift,adj=0.17,col="red",line=1.8)
mtext(diff,adj=0.17,col="red",line=0.4)
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=1,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1.1,cex=1,adj=1,col="blue")
mtext(bquote(x[.(0)]==.(X0)),line=0.1,cex=1,adj=0.75,col="blue")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=1,adj=0.75,col="blue")
mtext(bquote(X[v]==.(v)),line=0.1,cex=1,adj=0.5,col="red")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
Y <- rep(NA,M)
if (X0 > v){
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] <= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}else{
if(ELRENA=="No"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]}
if(ELRENA == "Yes"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau <- tau[c(which(!is.na(tau)))]
                 tau}
if(ELRENA=="Mean"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=mean(tau,na.rm=TRUE)
                 tau}
if(ELRENA=="Median"){for (i in 1:M){Y[i] = min(which(Q[,i] >= v))}
                 tau <- t[Y]
                 tau[c(which(is.na(tau)))]=median(tau,na.rm=TRUE)
                 tau}
}
Anay_RK3 <- data.frame(tau)
showData(Anay_RK3, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Anay_RK3, file = "AnaSimFTP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Anay_RK3)
names(Anay_RK3)
}
}


.hist_meth <-
function(X,Breaks,Prob = c("TRUE","FALSE"))
           {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

if (Breaks=='Sturges')   {Breaks <- nclass.Sturges(X)}
if (Breaks=='scott')     {Breaks <- nclass.scott(X)}
if (Breaks=='FD')        {Breaks <- nclass.FD(X)}

hist(X,breaks = Breaks,probability = Prob,col="light blue",border="dark blue",
     xlab = expression(X),las=1,main="")
box()
mtext(expression("Histogram for the Random Variable X"),line=2.5,adj=0.5,cex=1,col="black")
mtext(bquote(nclass==.(Breaks)),adj=0,line=0.2,cex=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}

.fctrep_Meth <-
function(X)
       {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

n <- length(X)
y <- sort(X)
f <- numeric()
for ( i in 1:n) {f[i] <- (i/n)}
plot(y,f,type="p",pch="*",las=1,xlab = expression(X),ylab="Frequence")
points(y,f,type="n")
mtext(expression("Empirical Distribution for the Random Variable X "),line=2.5,adj=0.5,cex=1,col="black")
legend("topleft",bg="gray85",border="gray",c("Empirical Distr"),pch=c("*"),col=c("black"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}

.Kern_meth <-
function(X,bw,k)
           {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

if (bw=='Irt')       {bw <- bw.nrd0(X)}
if (bw=='scott')     {bw <- bw.nrd(X)}
if (bw=='Ucv')       {bw <- bw.ucv(X)}
if (bw=='Bcv')       {bw <- bw.bcv(X)}
if (bw=='SJ')        {bw <- bw.SJ(X)}

plot(density(X, bw, kernel=k),xlab=expression(X),main="",las=1)
mtext(bquote("Estimated Density of the Random Variable X\nUsing The Kernel Method"),line=2.1,adj=0.5,cex=0.9,col="black")
mtext(bquote(Kernel== .(k)),line=1.1,adj=0,cex=0.9,col="red")
mtext(bquote(Bandwidth== .(round(bw,4))),line=0.25,adj=0,cex=0.9,col="red")
mtext(bquote(N == .(length(X))),line=0.25,adj=0.5,cex=0.9,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
   }

.Kern_general <-
function(Data,bw,k,
             Law=c("exp","GAmma","chisq","Beta","fisher",
                   "student","weibull","Normlog","Norm"))
           {
if (length(dim(Data)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

X <- Data
x <- seq(min(X),max(X),by=0.001)

if (bw=='Irt')       {bw <- bw.nrd0(X)}
if (bw=='scott')     {bw <- bw.nrd(X)}
if (bw=='Ucv')       {bw <- bw.ucv(X)}
if (bw=='Bcv')       {bw <- bw.bcv(X)}
if (bw=='SJ')        {bw <- bw.SJ(X)}

Law <- match.arg(Law)

if (Law=="exp")
      {
res  <- Ajdexp(X,starts = list(lambda = 1))
res1 <-ks.test(X,"dexp",res$coef[1])
plot(density(X, bw, kernel=k),xlab=expression(X),main="",las=1)
curve(dexp(x,res$coef[1]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Estimated Density of the Random Variable X\nUsing The Kernel Method "),line=2.5,adj=0.5,cex=0.9,col="black")
mtext(bquote(Kernel== .(k)),line=1.1,adj=0,cex=0.9,col="red")
mtext(bquote(Bandwidth== .(round(bw,4))),line=0.25,adj=0,cex=0.9,col="red")
mtext(bquote(N == .(length(X))),line=0.25,adj=0.5,cex=0.9,col="blue")
mtext(bquote("Law" :exp(lambda==.(round(res$coef[1],3)))),adj=0.5,line=1.1,cex=0.8,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("estimated density","exp Law"),col=c("black","red"),lwd=c(2,2),lty=c(1,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
   }

if (Law=="GAmma")
      {
res  <- Ajdgamma(X,starts = list(shape =1, rate=1))
res1 <-ks.test(X,"dgamma",res$coef[1],res$coef[2])
plot(density(X, bw, kernel=k),xlab=expression(X),main="",las=1)
curve(dgamma(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Estimated Density of the Random Variable X\nUsing The Kernel Method "),line=2.5,adj=0.5,cex=0.9,col="black")
mtext(bquote(Kernel== .(k)),line=1.1,adj=0,cex=0.9,col="red")
mtext(bquote(Bandwidth== .(round(bw,4))),line=0.25,adj=0,cex=0.9,col="red")
mtext(bquote(N == .(length(X))),line=0.25,adj=0.5,cex=0.9,col="blue")
mtext(bquote("Law" :Gamma( list(shape==.(round(res$coef[1],3)),rate==.(round(res$coef[2],3)))  )  ),adj=0.5,line=1.1,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("estimated density","gamma Law"),col=c("black","red"),lwd=c(2,2),lty=c(1,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="chisq")
      {
res  <- Ajdchisq(X,starts = list(df = 1))
res1 <-ks.test(X,"dchisq",res$coef[1])
plot(density(X, bw, kernel=k),xlab=expression(X),main="",las=1)
curve(dchisq(x,res$coef[1]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Estimated Density of the Random Variable X\nUsing The Kernel Method "),line=2.5,adj=0.5,cex=0.9,col="black")
mtext(bquote(Kernel== .(k)),line=1.1,adj=0,cex=0.9,col="red")
mtext(bquote(Bandwidth== .(round(bw,4))),line=0.25,adj=0,cex=0.9,col="red")
mtext(bquote(N == .(length(X))),line=0.25,adj=0.5,cex=0.9,col="blue")
mtext(bquote("Law" :chi[2](df%~~%.(round(res$coef[1],0)))),adj=0.5,line=1.1,cex=0.8,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("estimated density","chisq Law"),col=c("black","red"),lwd=c(2,2),lty=c(1,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="Beta")
      {
res  <- Ajdbeta(X,starts = list(shape1 =1, shape2=1))
res1 <-ks.test(X,"dbeta",res$coef[1],res$coef[2])
plot(density(X, bw, kernel=k),xlab=expression(X),main="",las=1)
curve(dbeta(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Estimated Density of the Random Variable X\nUsing The Kernel Method "),line=2.5,adj=0.5,cex=0.9,col="black")
mtext(bquote(Kernel== .(k)),line=1.1,adj=0,cex=0.9,col="red")
mtext(bquote(Bandwidth== .(round(bw,4))),line=0.25,adj=0,cex=0.9,col="red")
mtext(bquote(N == .(length(X))),line=0.25,adj=0.5,cex=0.9,col="blue")
mtext(bquote("Law" :beta( list(shape1==.(round(res$coef[1],3)),shape2==.(round(res$coef[2],3)))  )  ),adj=0.5,line=1.1,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("estimated density","beta Law"),col=c("black","red"),lwd=c(2,2),lty=c(1,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="fisher")
      {
res  <- Ajdf(X,starts = list(df1=1, df2=1))
res1 <-ks.test(X,"df",res$coef[1],res$coef[2])
plot(density(X, bw, kernel=k),xlab=expression(X),main="",las=1)
curve(df(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Estimated Density of the Random Variable X\nUsing The Kernel Method "),line=2.5,adj=0.5,cex=0.9,col="black")
mtext(bquote(Kernel== .(k)),line=1.1,adj=0,cex=0.9,col="red")
mtext(bquote(Bandwidth== .(round(bw,4))),line=0.25,adj=0,cex=0.9,col="red")
mtext(bquote(N == .(length(X))),line=0.25,adj=0.5,cex=0.9,col="blue")
mtext(bquote("Law" :F( list(df1%~~%.(round(res$coef[1],0)),df2%~~%.(round(res$coef[2],0)))  )  ),adj=0.5,line=1.1,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("estimated density","fisher Law"),col=c("black","red"),lwd=c(2,2),lty=c(1,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="student")
      {
res  <- Ajdt(X,starts = list(df=1))
res1 <-ks.test(X,"dt",res$coef[1])
plot(density(X, bw, kernel=k),xlab=expression(X),main="",las=1)
curve(dt(x,res$coef[1]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Estimated Density of the Random Variable X\nUsing The Kernel Method "),line=2.5,adj=0.5,cex=0.9,col="black")
mtext(bquote(Kernel== .(k)),line=1.1,adj=0,cex=0.9,col="red")
mtext(bquote(Bandwidth== .(round(bw,4))),line=0.25,adj=0,cex=0.9,col="red")
mtext(bquote(N == .(length(X))),line=0.25,adj=0.5,cex=0.9,col="blue")
mtext(bquote("Law" :St( list(df%~~%.(round(res$coef[1],0))))  ),adj=0.5,line=1.1,cex=0.8,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("estimated density","student Law"),col=c("black","red"),lwd=c(2,2),lty=c(1,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="weibull")
      {
res  <- Ajdweibull(X,starts = list(shape =1, scale=1))
res1 <-ks.test(X,"dweibull",res$coef[1],res$coef[2])
plot(density(X, bw, kernel=k),xlab=expression(X),main="",las=1)
curve(dweibull(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Estimated Density of the Random Variable X\nUsing The Kernel Method "),line=2.5,adj=0.5,cex=0.9,col="black")
mtext(bquote(Kernel== .(k)),line=1.1,adj=0,cex=0.9,col="red")
mtext(bquote(Bandwidth== .(round(bw,4))),line=0.25,adj=0,cex=0.9,col="red")
mtext(bquote(N == .(length(X))),line=0.25,adj=0.5,cex=0.9,col="blue")
mtext(bquote("Law" :W( list(shape==.(round(res$coef[1],3)),scale==.(round(res$coef[2],3)))  )  ),adj=0.5,line=1.1,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("estimated density","weibull Law"),col=c("black","red"),lwd=c(2,2),lty=c(1,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="Normlog")
      {
res  <- Ajdlognorm(X,starts = list(meanlog =1, sdlog=1))
res1 <-ks.test(X,"dlnorm",res$coef[1],res$coef[2])
plot(density(X, bw, kernel=k),xlab=expression(X),main="",las=1)
curve(dlnorm(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Estimated Density of the Random Variable X\nUsing The Kernel Method "),line=2.5,adj=0.5,cex=0.9,col="black")
mtext(bquote(Kernel== .(k)),line=1.1,adj=0,cex=0.9,col="red")
mtext(bquote(Bandwidth== .(round(bw,4))),line=0.25,adj=0,cex=0.9,col="red")
mtext(bquote(N == .(length(X))),line=0.25,adj=0.5,cex=0.9,col="blue")
mtext(bquote("Law" :logN( list(mean==.(round(res$coef[1],3)),sd==.(round(res$coef[2],3)))  )  ),adj=0.5,line=1.1,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("estimated density","log normal Law"),col=c("black","red"),lwd=c(2,2),lty=c(1,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="Norm")
      {
res  <- Ajdnorm(X,starts = list(mean =1, sd=1))
res1 <-ks.test(X,"dnorm",res$coef[1],res$coef[2])
plot(density(X, bw, kernel=k),xlab=expression(X),main="",las=1)
curve(dnorm(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Estimated Density of the Random Variable X\nUsing The Kernel Method "),line=2.5,adj=0.5,cex=0.9,col="black")
mtext(bquote(Kernel== .(k)),line=1.1,adj=0,cex=0.9,col="red")
mtext(bquote(Bandwidth== .(round(bw,4))),line=0.25,adj=0,cex=0.9,col="red")
mtext(bquote(N == .(length(X))),line=0.25,adj=0.5,cex=0.9,col="blue")
mtext(bquote("Law" :N( list(mean==.(round(res$coef[1],3)),sd==.(round(res$coef[2],3)))  )  ),adj=0.5,line=1.1,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("estimated density","normal Law"),col=c("black","red"),lwd=c(2,2),lty=c(1,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}
}

.fctgeneral <-
function(Data,Law=c("exp","GAmma","chisq","Beta","fisher","student","weibull","Normlog","Norm"))
       {
if (length(dim(Data)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

Law <- match.arg(Law)
X <- Data
x <- seq(min(X),max(X),by=0.001)
if (Law=="exp")
      {
res  <- Ajdexp(X,starts = list(lambda = 1))
res1 <-ks.test(X,"pexp",res$coef[1])
n <- length(X)
y <- sort(X)
f <- numeric()
for ( i in 1:n)
     {
     f[i] <- (i/n)
     }
plot(y,f,type="p",pch="*",las=1,xlab = expression(X),ylab="Frequence")
curve(pexp(x,res$coef[1]), col = 2, lwd = 2, add = TRUE)
mtext(expression("Empirical Distribution for the Random Variable X "),line=2.8,adj=0.5,cex=1,col="black")
mtext(bquote("Law" :exp(lambda==.(round(res$coef[1],3)))),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topleft",border="gray",c("empirical Distr","exp Law"),pch=c("*",""),col=c("black","red"),lwd=c("",2),lty=c(0,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="GAmma")
      {
res  <- Ajdgamma(X,starts = list(shape =1, rate=1))
res1 <-ks.test(X,"pgamma",res$coef[1],res$coef[2])
n <- length(X)
y <- sort(X)
f <- numeric()
for ( i in 1:n)
     {
     f[i] <- (i/n)
     }
plot(y,f,type="p",pch="*",las=1,xlab = expression(X),ylab="Frequence")
curve(pgamma(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(expression("Empirical Distribution for the Random Variable X "),line=2.8,adj=0.5,cex=1,col="black")
mtext(bquote("Law" :Gamma( list(shape==.(round(res$coef[1],3)),rate==.(round(res$coef[2],3)))  )  ),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topleft",border="gray",c("empirical Distr","gamma Law"),pch=c("*",""),col=c("black","red"),lwd=c("",2),lty=c(0,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="chisq")
      {
res  <- Ajdchisq(X,starts = list(df = 1))
res1 <-ks.test(X,"pchisq",res$coef[1])
n <- length(X)
y <- sort(X)
f <- numeric()
for ( i in 1:n)
     {
     f[i] <- (i/n)
     }
plot(y,f,type="p",pch="*",las=1,xlab = expression(X),ylab="Frequence")
curve(pchisq(x,res$coef[1]), col = 2, lwd = 2, add = TRUE)
mtext(expression("Empirical Distribution for the Random Variable X "),line=2.8,adj=0.5,cex=1,col="black")
mtext(bquote("Law" :chi[2](df%~~%.(round(res$coef[1],0)))),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topleft",border="gray",c("empirical Distr","chisq Law"),pch=c("*",""),col=c("black","red"),lwd=c("",2),lty=c(0,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="Beta")
      {
res  <- Ajdbeta(X,starts = list(shape1 =1, shape2=1))
res1 <-ks.test(X,"pbeta",res$coef[1],res$coef[2])
n <- length(X)
y <- sort(X)
f <- numeric()
for ( i in 1:n)
     {
     f[i] <- (i/n)
     }
plot(y,f,type="p",pch="*",las=1,xlab = expression(X),ylab="Frequence")
curve(pbeta(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(expression("Empirical Distribution for the Random Variable X "),line=2.8,adj=0.5,cex=1,col="black")
mtext(bquote("Law" :beta( list(shape1==.(round(res$coef[1],3)),shape2==.(round(res$coef[2],3)))  )  ),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topleft",border="gray",c("empirical Distr","beta Law"),pch=c("*",""),col=c("black","red"),lwd=c("",2),lty=c(0,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="fisher")
      {
res  <- Ajdf(X,starts = list(df1=1, df2=1))
res1 <-ks.test(X,"pf",res$coef[1],res$coef[2])
n <- length(X)
y <- sort(X)
f <- numeric()
for ( i in 1:n)
     {
     f[i] <- (i/n)
     }
plot(y,f,type="p",pch="*",las=1,xlab = expression(X),ylab="Frequence")
curve(pf(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(expression("Empirical Distribution for the Random Variable X "),line=2.8,adj=0.5,cex=1,col="black")
mtext(bquote("Law" :F( list(df1%~~%.(round(res$coef[1],0)),df2%~~%.(round(res$coef[2],0)))  )  ),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topleft",border="gray",c("empirical Distr","fisher Law"),pch=c("*",""),col=c("black","red"),lwd=c("",2),lty=c(0,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="student")
      {
res  <- Ajdt(X,starts = list(df=1))
res1 <-ks.test(X,"pt",res$coef[1])
n <- length(X)
y <- sort(X)
f <- numeric()
for ( i in 1:n)
     {
     f[i] <- (i/n)
     }
plot(y,f,type="p",pch="*",las=1,xlab = expression(X),ylab="Frequence")
curve(pt(x,res$coef[1]), col = 2, lwd = 2, add = TRUE)
mtext(expression("Empirical Distribution for the Random Variable X "),line=2.8,adj=0.5,cex=1,col="black")
mtext(bquote("Law" :St( list(df%~~%.(round(res$coef[1],0))))  ),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topleft",border="gray",c("empirical Distr","student Law"),pch=c("*",""),col=c("black","red"),lwd=c("",2),lty=c(0,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="weibull")
      {
res  <- Ajdweibull(X,starts = list(shape =1, scale=1))
res1 <-ks.test(X,"pweibull",res$coef[1],res$coef[2])
n <- length(X)
y <- sort(X)
f <- numeric()
for ( i in 1:n)
     {
     f[i] <- (i/n)
     }
plot(y,f,type="p",pch="*",las=1,xlab = expression(X),ylab="Frequence")
curve(pweibull(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(expression("Empirical Distribution for the Random Variable X "),line=2.8,adj=0.5,cex=1,col="black")
mtext(bquote("Law" :W( list(shape==.(round(res$coef[1],3)),scale==.(round(res$coef[2],3)))  )  ),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topleft",border="gray",c("empirical Distr","weibull Law"),pch=c("*",""),col=c("black","red"),lwd=c("",2),lty=c(0,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="Normlog")
      {
res  <- Ajdlognorm(X,starts = list(meanlog =1, sdlog=1))
res1 <-ks.test(X,"plnorm",res$coef[1],res$coef[2])
n <- length(X)
y <- sort(X)
f <- numeric()
for ( i in 1:n)
     {
     f[i] <- (i/n)
     }
plot(y,f,type="p",pch="*",las=1,xlab = expression(X),ylab="Frequence")
curve(plnorm(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(expression("Empirical Distribution for the Random Variable X "),line=2.8,adj=0.5,cex=1,col="black")
mtext(bquote("Law" :logN( list(mean==.(round(res$coef[1],3)),sd==.(round(res$coef[2],3)))  )  ),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topleft",border="gray",c("empirical Distr","Log Normal Law"),pch=c("*",""),col=c("black","red"),lwd=c("",2),lty=c(0,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}


if (Law=="Norm")
      {
res  <- Ajdnorm(X,starts = list(mean =1, sd=1))
res1 <-ks.test(X,"pnorm",res$coef[1],res$coef[2])
n <- length(X)
y <- sort(X)
f <- numeric()
for ( i in 1:n)
     {
     f[i] <- (i/n)
     }
plot(y,f,type="p",pch="*",las=1,xlab = expression(X),ylab="Frequence")
curve(pnorm(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(expression("Empirical Distribution for the Random Variable X "),line=2.8,adj=0.5,cex=1,col="black")
mtext(bquote("Law" :N( list(mean==.(round(res$coef[1],3)),sd==.(round(res$coef[2],3)))  )  ),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topleft",border="gray",c("empirical Distr","Normal Law"),pch=c("*",""),col=c("black","red"),lwd=c("",2),lty=c(0,1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}
}

.hist_general <-
function(Data,Breaks,
             Law=c("exp","GAmma","chisq","Beta","fisher",
                   "student","weibull","Normlog","Norm"))
           {
if (length(dim(Data)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

X <- Data
x <- seq(min(X),max(X),by=0.001)			
if (Breaks=='Sturges')   {Breaks <- nclass.Sturges(X)}
if (Breaks=='scott')     {Breaks <- nclass.scott(X)}
if (Breaks=='FD')        {Breaks <- nclass.FD(X)}

Law <- match.arg(Law)

if (Law=="exp")
      {
res  <- Ajdexp(X,starts = list(lambda = 1))
res1 <-ks.test(X,"dexp",res$coef[1])
hist(X,breaks = Breaks,probability = TRUE,col="light blue",border="dark blue",
     xlab = expression(X),las=1,main="")
box()
mtext(expression("Histogram for the Random Variable X"),line=2.5,adj=0.5,cex=1,col="black")
mtext(bquote(nclass==.(Breaks)),adj=0.5,line=0.2,cex=1,col="blue")
curve(dexp(x,res$coef[1]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Law" :exp(lambda==.(round(res$coef[1],3)))),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("exp Law"),col=c("red"),lwd=c(2),lty=c(1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="GAmma")
      {
res  <- Ajdgamma(X,starts = list(shape =1, rate=1))
res1 <-ks.test(X,"dgamma",res$coef[1],res$coef[2])
hist(X,breaks = Breaks,probability = TRUE,col="light blue",border="dark blue",
     xlab = expression(X),las=1,main="")
box()
mtext(expression("Histogram for the Random Variable X"),line=2.5,adj=0.5,cex=1,col="black")
mtext(bquote(nclass==.(Breaks)),adj=0.5,line=0.2,cex=1,col="blue")
curve(dgamma(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Law" :Gamma( list(shape==.(round(res$coef[1],3)),rate==.(round(res$coef[2],3)))  )  ),adj=0,line=0.2,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("gamma Law"),col=c("red"),lwd=c(2),lty=c(1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="chisq")
      {
res  <- Ajdchisq(X,starts = list(df = 1))
res1 <-ks.test(X,"dchisq",res$coef[1])
hist(X,breaks = Breaks,probability = TRUE,col="light blue",border="dark blue",
     xlab = expression(X),las=1,main="")
box()
mtext(expression("Histogram for the Random Variable X"),line=2.5,adj=0.5,cex=1,col="black")
mtext(bquote(nclass==.(Breaks)),adj=0.5,line=0.2,cex=1,col="blue")
curve(dchisq(x,res$coef[1]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Law" :chi[2](df%~~%.(round(res$coef[1],0)))),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("chisq Law"),col=c("red"),lwd=c(2),lty=c(1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="Beta")
      {
res  <- Ajdbeta(X,starts = list(shape1 =1, shape2=1))
res1 <-ks.test(X,"dbeta",res$coef[1],res$coef[2])
hist(X,breaks = Breaks,probability = TRUE,col="light blue",border="dark blue",
     xlab = expression(X),las=1,main="")
box()
mtext(expression("Histogram for the Random Variable X"),line=2.5,adj=0.5,cex=1,col="black")
mtext(bquote(nclass==.(Breaks)),adj=0.65,line=0.2,cex=1,col="blue")
curve(dbeta(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Law" :beta( list(shape1==.(round(res$coef[1],3)),shape2==.(round(res$coef[2],3)))  )  ),adj=0,line=0.2,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("beta Law"),col=c("red"),lwd=c(2),lty=c(1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="fisher")
      {
res  <- Ajdf(X,starts = list(df1=1, df2=1))
res1 <-ks.test(X,"df",res$coef[1],res$coef[2])
hist(X,breaks = Breaks,probability = TRUE,col="light blue",border="dark blue",
     xlab = expression(X),las=1,main="")
box()
mtext(expression("Histogram for the Random Variable X"),line=2.5,adj=0.5,cex=1,col="black")
mtext(bquote(nclass==.(Breaks)),adj=0.65,line=0.2,cex=1,col="blue")
curve(df(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Law" :F( list(df1%~~%.(round(res$coef[1],0)),df2%~~%.(round(res$coef[2],0)))  )  ),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("fisher Law"),col=c("red"),lwd=c(2),lty=c(1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="student")
      {
res  <- Ajdt(X,starts = list(df=1))
res1 <-ks.test(X,"dt",res$coef[1])
hist(X,breaks = Breaks,probability = TRUE,col="light blue",border="dark blue",
     xlab = expression(X),las=1,main="")
box()
mtext(expression("Histogram for the Random Variable X"),line=2.5,adj=0.5,cex=1,col="black")
mtext(bquote(nclass==.(Breaks)),adj=0.5,line=0.2,cex=1,col="blue")
curve(dt(x,res$coef[1]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Law" :St( list(df%~~%.(round(res$coef[1],0))))  ),adj=0,line=0.2,cex=1,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("student Law"),col=c("red"),lwd=c(2),lty=c(1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="weibull")
      {
res  <- Ajdweibull(X,starts = list(shape =1, scale=1))
res1 <-ks.test(X,"dweibull",res$coef[1],res$coef[2])
hist(X,breaks = Breaks,probability = TRUE,col="light blue",border="dark blue",
     xlab = expression(X),las=1,main="")
box()
mtext(expression("Histogram for the Random Variable X"),line=2.5,adj=0.5,cex=1,col="black")
mtext(bquote(nclass==.(Breaks)),adj=0.65,line=0.2,cex=1,col="blue")
curve(dweibull(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Law" :W( list(shape==.(round(res$coef[1],3)),scale==.(round(res$coef[2],3)))  )  ),adj=0,line=0.2,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("weibull Law"),col=c("red"),lwd=c(2),lty=c(1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="Normlog")
      {
res  <- Ajdlognorm(X,starts = list(meanlog =1, sdlog=1))
res1 <-ks.test(X,"dlnorm",res$coef[1],res$coef[2])
hist(X,breaks = Breaks,probability = TRUE,col="light blue",border="dark blue",
     xlab = expression(X),las=1,main="")
box()
mtext(expression("Histogram for the Random Variable X"),line=2.5,adj=0.5,cex=1,col="black")
mtext(bquote(nclass==.(Breaks)),adj=0.65,line=0.2,cex=1,col="blue")
curve(dlnorm(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Law" :logN( list(mean==.(round(res$coef[1],3)),sd==.(round(res$coef[2],3)))  )  ),adj=0,line=0.2,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("Log Normal Law"),col=c("red"),lwd=c(2),lty=c(1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}

if (Law=="Norm")
      {
res  <- Ajdnorm(X,starts = list(mean =1, sd=1))
res1 <-ks.test(X,"dnorm",res$coef[1],res$coef[2])
hist(X,breaks = Breaks,probability = TRUE,col="light blue",border="dark blue",
     xlab = expression(X),las=1,main="")
box()
mtext(expression("Histogram for the Random Variable X"),line=2.5,adj=0.5,cex=1,col="black")
mtext(bquote(nclass==.(Breaks)),adj=0.65,line=0.2,cex=1,col="blue")
curve(dnorm(x,res$coef[1],res$coef[2]), col = 2, lwd = 2, add = TRUE)
mtext(bquote("Law" :N( list(mean==.(round(res$coef[1],3)),sd==.(round(res$coef[2],3)))  )  ),adj=0,line=0.2,cex=0.7,col="red")
mtext(bquote(AIC ==.(round(res$AIC[1],3) )),adj=1,line=2,cex=0.8,col="blue")
mtext(bquote(p.value==.(round(res1$p.value,3) )),adj=1,line=1,cex=0.8,col="blue")
mtext(bquote(D.statistic==.(round(res1$statistic,3) )),adj=1,line=0.3,cex=0.8,col="blue")
legend("topright",border="gray",c("Normal Law"),col=c("red"),lwd=c(2),lty=c(1),cex=0.7)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
return(invisible(res))
return(invisible(res1))
}
}

.Ajdexp <-
function(X,starts = list(lambda = 1), leve = 0.95)
       {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

cd <- function(x,lambda,log = TRUE)
{
lik <- dexp(x,lambda  , log = TRUE )
if(!log )
lik <- exp(lik)
lik
}

lik <-function(lambda) 
       {
n  <- length(X)
-sum(cd(x=X[2:n],lambda,log = TRUE ))
       }
res <- mle(lik, start = starts)
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.test_ks_dexp <-
function(X,lambda)
             {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

{return(print(ks.test(X,"pexp",lambda)))}
}

.Ajdgamma <-
function(X,starts = list(shape =1, rate=1), leve = 0.95)
       {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

cd <- function(x,shape, rate,log = TRUE)
{
lik <- dgamma(x,shape, rate ,log = TRUE)
if(!log )
lik <- exp(lik)
lik
}
lik <-function(shape, rate) 
       {
n  <- length(X)
-sum(cd(x=X[2:n],shape, rate,log = TRUE ))
       }
res <- mle(lik, start = starts)
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.test_ks_dgamma <-
function(X,shape, rate)
             {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

{return(print(ks.test(X,"pgamma",shape, rate)))}
}

.Ajdchisq <-
function(X,starts = list(df=1), leve = 0.95)
       {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

cd <- function(x,df,log = TRUE)
{
lik <- dchisq(x,df, ncp=0 ,log = TRUE)
if(!log )
lik <- exp(lik)
lik
}
lik <-function(df) 
       {
n  <- length(X)
-sum(cd(x=X[2:n],df,log = TRUE ))
       }
res <- mle(lik, start = starts)
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.test_ks_dchisq <-
function(X,df)
             {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

{return(print(ks.test(X,"pchisq",df)))}
}

.Ajdbeta <-
function(X,starts = list(shape1=1, shape2=1), leve = 0.95)
       {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

cd <- function(x,shape1, shape2,log = TRUE)
{
lik <- dbeta(x,shape1, shape2, ncp=0 ,log = TRUE)
if(!log )
lik <- exp(lik)
lik
}
lik <-function(shape1, shape2) 
       {
n  <- length(X)
-sum(cd(x=X[2:n],shape1, shape2,log = TRUE ))
       }
res <- mle(lik, start = starts)
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.test_ks_dbeta <-
function(X,shape1, shape2)
             {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

{return(print(ks.test(X,"pbeta",shape1, shape2)))}
}

.Ajdf <-
function(X,starts = list(df1=1, df2=1), leve = 0.95)
       {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

cd <- function(x,df1, df2,log = TRUE)
{
lik <- df(x,df1, df2,log = TRUE)
if(!log )
lik <- exp(lik)
lik
}
lik <-function(df1, df2) 
       {
n  <- length(X)
-sum(cd(x=X[2:n],df1, df2,log = TRUE ))
       }
res <- mle(lik, start = starts)
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.test_ks_df <-
function(X,df1, df2)
             {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

{return(print(ks.test(X,"pf",df1, df2)))}
}

.Ajdt <-
function(X,starts = list(df=1), leve = 0.95)
       {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

cd <- function(x,df,log = TRUE)
{
lik <- dt(x,df,log = TRUE)
if(!log )
lik <- exp(lik)
lik
}
lik <-function(df) 
       {
n  <- length(X)
-sum(cd(x=X[2:n],df,log = TRUE ))
       }
res <- mle(lik, start = starts)
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.test_ks_dt <-
function(X,df)
             {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

{return(print(ks.test(X,"pt",df)))}
}

.Ajdweibull <-
function(X,starts = list(shape =1, scale=1), leve = 0.95)
       {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

cd <- function(x,shape, scale,log = TRUE)
{
lik <- dweibull(x,shape, scale ,log = TRUE)
if(!log )
lik <- exp(lik)
lik
}
lik <-function(shape, scale) 
       {
n  <- length(X)
-sum(cd(x=X[2:n],shape, scale,log = TRUE ))
       }
res <- mle(lik, start = starts)
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.test_ks_dweibull <-
function(X,shape, scale)
             {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

{return(print(ks.test(X,"pweibull",shape, scale)))}
}

.Ajdnorm <-
function(X,starts = list(mean =1, sd=1), leve = 0.95)
       {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

cd <- function(x,mean, sd,log = TRUE)
{
lik <- dnorm(x,mean, sd ,log = TRUE)
if(!log )
lik <- exp(lik)
lik
}
lik <-function(mean, sd) 
       {
n  <- length(X)
-sum(cd(x=X[2:n],mean, sd,log = TRUE ))
       }
res <- mle(lik, start = starts)
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.test_ks_dnorm <-
function(X,mean, sd)
             {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

{return(print(ks.test(X,"pnorm",mean, sd)))}
}

.Ajdlognorm <-
function(X,starts = list(meanlog =1, sdlog=1), leve = 0.95)
       {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

cd <- function(x,meanlog, sdlog,log = TRUE)
{
lik <- dlnorm(x,meanlog, sdlog ,log = TRUE)
if(!log )
lik <- exp(lik)
lik
}
lik <-function(meanlog, sdlog) 
       {
n  <- length(X)
-sum(cd(x=X[2:n],meanlog, sdlog,log = TRUE ))
       }
res <- mle(lik, start = starts)
{return(print(list(summary = summary(res), coef = coef(res), 
        AIC = AIC(res), vcov = vcov(res), confint = confint(res, 
        level = leve))))}
}

.test_ks_dlognorm <-
function(X,meanlog, sdlog)
             {
if (length(dim(X)) > 0) 
            stop(tkmessageBox(title="Error",message=paste( "data must be numeric vector" ),icon="error"))

{return(print(ks.test(X,"plnorm",meanlog, sdlog)))}
}

.tho_1 <-
function(N,t0,T,R0,v,K,Sigma,
                    Methods = c("Euler", "Milstein", "MilsteinS", 
                                "Ito-Taylor", "Heun", "RK3"), ...)
{

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))


if( R0 < 0 || v <= 0)
            stop(tkmessageBox(title="Error",message=paste( "R0 > 0 and v > 0" ),icon="error"))

if( R0 <= v)
            stop(tkmessageBox(title="Error",message=paste( "R0 > v" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( Sigma < 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K < 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

Methods <- match.arg(Methods)

if(Methods=="Euler")
                  {
Eul <- function(N,T,t0,R0,v,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]}      
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-Eul(N,T,t0,R0,v,K,Sigma)
return(invisible(X))
}

if(Methods=="Milstein")
                  {
Mil <- function(N,T,t0,R0,v,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

DSx  <- D(diff,"x")
A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)

Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
            R[i] = R[i-1]+ A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]+ 
                   0.5 *S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1])^2 -Dt)}     
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-Mil(N,T,t0,R0,v,K,Sigma)
return(invisible(X))
}

if(Methods=="MilsteinS")
                  {
MilS <- function(N,T,t0,R0,v,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
             R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1] +
                    0.5 *S(t[i-1],R[i-1]) * Sx(t[i-1],R[i-1])*(D[i-1]^2-Dt)+ 
                    Dt^(3 /2)*(0.5 *A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1]) +
                    0.5 *Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])+
                    0.25 *(S(t[i -1] ,R[i -1])^2) * Sxx(t[i -1] ,R[i -1]))*D[i -1]+ 
                   (Dt^2) * (0.5*A(t[i -1],R[i -1])*Ax(t[i-1],R[i-1])+
                    0.25 *Axx(t[i-1],R[i-1])*(S(t[i-1],R[i-1])^2))}     
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-MilS(N,T,t0,R0,v,K,Sigma)
return(invisible(X))
}

if(Methods=="Ito-Taylor")
                  {
ITO <- function(N,T,t0,R0,v,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
DZ= rnorm(N,0,sqrt((1/3)*Dt^3))
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
                 R[i]=R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]+
                      0.5*S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1]^2)-Dt)+
                      Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])*DZ[i-1]+0.5*(A(t[i-1],R[i-1])*Ax(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Axx(t[i-1],R[i-1]))*(Dt^2)+(A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Sxx(t[i-1],R[i-1]))*(D[i-1]*Dt-DZ[i-1])+
                      0.5*S(t[i-1],R[i-1])*(S(t[i-1],R[i-1])*Sxx(t[i-1],R[i-1])+
                     (Sx(t[i-1],R[i-1])^2))*((1/3)*(D[i-1]^2)-Dt)*D[i-1]}     
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-ITO(N,T,t0,R0,v,K,Sigma)
return(invisible(X))
}

if(Methods=="Heun")
                  {
Heu <- function(N,T,t0,R0,v,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
         Y[i-1]= R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]
         R[i]  = R[i-1]+0.5*Dt*(A(t[i-1],R[i-1])+A(t[i-1],Y[i-1]))+
                 0.5*(S(t[i-1],R[i-1])+S(t[i-1],Y[i-1]))*D[i-1]}     
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-Heu(N,T,t0,R0,v,K,Sigma)
return(invisible(X))
}

if(Methods=="RK3")
                  {
RK <- function(N,T,t0,R0,v,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
              Y    <- numeric()
              Z    <- numeric()
              Y[i-1]=R[i-1]+0.5*Dt*A(t[i-1],R[i-1])+S(t[i-1],R[i-1])*D[i-1]
              Z[i-1]=R[i-1]-A(t[i-1],R[i-1])*Dt+2*Dt*A(t[i-1]+0.5*Dt,Y[i-1])+
                     (2*S(t[i-1]+0.5*Dt,Y[i-1])-S(t[i-1],R[i-1]))*D[i-1]
              R[i] = R[i-1]+(Dt/6)*(A(t[i-1],R[i-1])+4*A(t[i-1]+0.5*Dt,Y[i-1])+A(t[i-1]+Dt,Z[i-1]))+
                     (1/6)*(S(t[i-1],R[i-1])+4*S(t[i-1]+0.5*Dt,Y[i-1])+S(t[i-1]+Dt,Z[i-1]))*D[i-1]}     
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-RK(N,T,t0,R0,v,K,Sigma)
return(invisible(X))
}
}

.tho_2 <-
function(N,t0,T,R0,v,K,s,Sigma,
                    Methods = c("Euler", "Milstein", "MilsteinS", 
                                "Ito-Taylor", "Heun", "RK3"), ...)
{

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))


if( R0 < 0 || v <= 0)
            stop(tkmessageBox(title="Error",message=paste( "R0 > 0 and v > 0" ),icon="error"))

if( R0 <= v)
            stop(tkmessageBox(title="Error",message=paste( "R0 > v" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( Sigma < 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K < 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

if( s <= 1 )
            stop(tkmessageBox(title="Error",message=paste( "s > 1" ),icon="error"))

Methods <- match.arg(Methods)

if(Methods=="Euler")
                  {
Eul <- function(N,T,t0,R0,v,K,s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]}      
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-Eul(N,T,t0,R0,v,K,s,Sigma)
return(invisible(X))
}

if(Methods=="Milstein")
                  {
Mil <- function(N,T,t0,R0,v,K,s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

DSx  <- D(diff,"x")
A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)

Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
            R[i] = R[i-1]+ A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]+ 
                   0.5 *S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1])^2 -Dt)}     
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-Mil(N,T,t0,R0,v,K,s,Sigma)
return(invisible(X))
}

if(Methods=="MilsteinS")
                  {
MilS <- function(N,T,t0,R0,v,K,s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
             R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1] +
                    0.5 *S(t[i-1],R[i-1]) * Sx(t[i-1],R[i-1])*(D[i-1]^2-Dt)+ 
                    Dt^(3 /2)*(0.5 *A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1]) +
                    0.5 *Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])+
                    0.25 *(S(t[i -1] ,R[i -1])^2) * Sxx(t[i -1] ,R[i -1]))*D[i -1]+ 
                   (Dt^2) * (0.5*A(t[i -1],R[i -1])*Ax(t[i-1],R[i-1])+
                    0.25 *Axx(t[i-1],R[i-1])*(S(t[i-1],R[i-1])^2))}     
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-MilS(N,T,t0,R0,v,K,s,Sigma)
return(invisible(X))
}

if(Methods=="Ito-Taylor")
                  {
ITO <- function(N,T,t0,R0,v,K,s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
DZ= rnorm(N,0,sqrt((1/3)*Dt^3))
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
                 R[i]=R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]+
                      0.5*S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1]^2)-Dt)+
                      Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])*DZ[i-1]+0.5*(A(t[i-1],R[i-1])*Ax(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Axx(t[i-1],R[i-1]))*(Dt^2)+(A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Sxx(t[i-1],R[i-1]))*(D[i-1]*Dt-DZ[i-1])+
                      0.5*S(t[i-1],R[i-1])*(S(t[i-1],R[i-1])*Sxx(t[i-1],R[i-1])+
                     (Sx(t[i-1],R[i-1])^2))*((1/3)*(D[i-1]^2)-Dt)*D[i-1]}     
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-ITO(N,T,t0,R0,v,K,s,Sigma)
return(invisible(X))
}

if(Methods=="Heun")
                  {
Heu <- function(N,T,t0,R0,v,K,s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
         Y[i-1]= R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]
         R[i]  = R[i-1]+0.5*Dt*(A(t[i-1],R[i-1])+A(t[i-1],Y[i-1]))+
                 0.5*(S(t[i-1],R[i-1])+S(t[i-1],Y[i-1]))*D[i-1]}     
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-Heu(N,T,t0,R0,v,K,s,Sigma)
return(invisible(X))
}

if(Methods=="RK3")
                  {
RK <- function(N,T,t0,R0,v,K,s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

Dt = (T-t0)/N
t <- seq (t0 ,T, length =N+1)
u = runif(N,0,1)
o = rep(1,N)
o [ which(u < 0.5) ] = -1
w = cumsum(c(0,o))*sqrt((T-t0)/N)
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
              Y    <- numeric()
              Z    <- numeric()
              Y[i-1]=R[i-1]+0.5*Dt*A(t[i-1],R[i-1])+S(t[i-1],R[i-1])*D[i-1]
              Z[i-1]=R[i-1]-A(t[i-1],R[i-1])*Dt+2*Dt*A(t[i-1]+0.5*Dt,Y[i-1])+
                     (2*S(t[i-1]+0.5*Dt,Y[i-1])-S(t[i-1],R[i-1]))*D[i-1]
              R[i] = R[i-1]+(Dt/6)*(A(t[i-1],R[i-1])+4*A(t[i-1]+0.5*Dt,Y[i-1])+A(t[i-1]+Dt,Z[i-1]))+
                     (1/6)*(S(t[i-1],R[i-1])+4*S(t[i-1]+0.5*Dt,Y[i-1])+S(t[i-1]+Dt,Z[i-1]))*D[i-1]}     
if (length(which(R <= v)) > 0 )  {xx  <- t[min(which(R <= v))]} else {xx <- NA}
}
X <-RK(N,T,t0,R0,v,K,s,Sigma)
return(invisible(X))
}
}

.Sim_tho02diff <-
function(N,t0,Dt,T=1,X1_0,X2_0,Y1_0,Y2_0,v,K,m,Sigma)
              {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if( v <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "v > 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( Sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

if( m <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "m > 0" ),icon="error"))  
  

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
wx1 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx1     <- diff(wx1)

wx2 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx2     <- diff(wx2)

wy1 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy1     <- diff(wy1)

wy2 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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

if (length(which(D <= v)) > 0 )  {xx  <- t[min(which(D <= v))]} else {xx <- NA}
xx
}

.RadialP_1flow <-
function(N,M,t0,Dt,T=1,R0,K,Sigma,Output=FALSE,
                   Methods = c("Euler", "Milstein", "MilsteinS", 
                                "Ito-Taylor", "Heun", "RK3"), ...)
{
##
##if( t0 >= T || t0 < 0 ) 
##            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
##
##if( N <= 1 )   
##            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
##
##if (M < 1)
##            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))
##
##if ( Dt <= 0 )
##            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))
##
##
##if( R0 < 0 )
##            stop(tkmessageBox(title="Error",message=paste( "R0 > 0" ),icon="error"))
##
##if( 2 * K <=  Sigma^2 )
##            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))
##
##if( Sigma < 0 )
##           stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))
##
##if( K < 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

Methods <- match.arg(Methods)

if(Methods=="Euler")
                {
Eul <- function(N,T=1,Dt,t0,R0,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]}      
R
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Eul,T=T,t0=t0,R0=R0,Dt=Dt,K=K,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Euler scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_1_Euler <- data.frame(time,R)
if (M >=2) {RadialP_1_Euler  <- data.frame(time,R,R.mean)}
showData(RadialP_1_Euler, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_Euler, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_Euler)
}

if(Methods=="Milstein")
                  {
Mil <- function(N,T=1,Dt,t0,R0,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

DSx  <- D(diff,"x")
A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
            R[i] = R[i-1]+ A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]+ 
                   0.5 *S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1])^2 -Dt)}     
R
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Mil,T=T,t0=t0,R0=R0,Dt=Dt,K=K,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Milstein scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_1_Milstein <- data.frame(time,R)
if (M >=2) {RadialP_1_Milstein  <- data.frame(time,R,R.mean)}
showData(RadialP_1_Milstein, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_Milstein, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_Milstein)
}

if(Methods=="MilsteinS")
                  {
MilS <- function(N,T=1,Dt,t0,R0,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
             R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1] +
                    0.5 *S(t[i-1],R[i-1]) * Sx(t[i-1],R[i-1])*(D[i-1]^2-Dt)+ 
                    Dt^(3 /2)*(0.5 *A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1]) +
                    0.5 *Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])+
                    0.25 *(S(t[i -1] ,R[i -1])^2) * Sxx(t[i -1] ,R[i -1]))*D[i -1]+ 
                   (Dt^2) * (0.5*A(t[i -1],R[i -1])*Ax(t[i-1],R[i-1])+
                    0.25 *Axx(t[i-1],R[i-1])*(S(t[i-1],R[i-1])^2))}     
R
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),MilS,T=T,t0=t0,R0=R0,Dt=Dt,K=K,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Second Milstein scheme")),line=3.3,adj=1,cex=0.6,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_1_MilsteinS <- data.frame(time,R)
if (M >=2) {RadialP_1_MilsteinS  <- data.frame(time,R,R.mean)}
showData(RadialP_1_MilsteinS, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_MilsteinS, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_MilsteinS)
}

if(Methods=="Ito-Taylor")
                  {
ITO <- function(N,T=1,Dt,t0,R0,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
DZ= rnorm(N,0,sqrt((1/3)*Dt^3))
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
                 R[i]=R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]+
                      0.5*S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1]^2)-Dt)+
                      Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])*DZ[i-1]+0.5*(A(t[i-1],R[i-1])*Ax(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Axx(t[i-1],R[i-1]))*(Dt^2)+(A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Sxx(t[i-1],R[i-1]))*(D[i-1]*Dt-DZ[i-1])+
                      0.5*S(t[i-1],R[i-1])*(S(t[i-1],R[i-1])*Sxx(t[i-1],R[i-1])+
                     (Sx(t[i-1],R[i-1])^2))*((1/3)*(D[i-1]^2)-Dt)*D[i-1]}     
R
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),ITO,T=T,t0=t0,R0=R0,Dt=Dt,K=K,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Strong Ito-Taylor \nScheme Order 1.5")),line=3,adj=1,cex=0.6,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_1_ITY <- data.frame(time,R)
if (M >=2) {RadialP_1_ITY  <- data.frame(time,R,R.mean)}
showData(RadialP_1_ITY, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_ITY, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_ITY)
}

if(Methods=="Heun")
                  {
Heu <- function(N,T=1,Dt,t0,R0,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
         Y[i-1]= R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]
         R[i]  = R[i-1]+0.5*Dt*(A(t[i-1],R[i-1])+A(t[i-1],Y[i-1]))+
                 0.5*(S(t[i-1],R[i-1])+S(t[i-1],Y[i-1]))*D[i-1]}     
R
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Heu,T=T,t0=t0,R0=R0,Dt=Dt,K=K,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Heun scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_1_Heun <- data.frame(time,R)
if (M >=2) {RadialP_1_Heun  <- data.frame(time,R,R.mean)}
showData(RadialP_1_Heun, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_Heun, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_Heun)
}

if(Methods=="RK3")
                  {
RK <- function(N,T=1,Dt,t0,R0,K,Sigma)
     {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
              Y    <- numeric()
              Z    <- numeric()
              Y[i-1]=R[i-1]+0.5*Dt*A(t[i-1],R[i-1])+S(t[i-1],R[i-1])*D[i-1]
              Z[i-1]=R[i-1]-A(t[i-1],R[i-1])*Dt+2*Dt*A(t[i-1]+0.5*Dt,Y[i-1])+
                     (2*S(t[i-1]+0.5*Dt,Y[i-1])-S(t[i-1],R[i-1]))*D[i-1]
              R[i] = R[i-1]+(Dt/6)*(A(t[i-1],R[i-1])+4*A(t[i-1]+0.5*Dt,Y[i-1])+A(t[i-1]+Dt,Z[i-1]))+
                     (1/6)*(S(t[i-1],R[i-1])+4*S(t[i-1]+0.5*Dt,Y[i-1])+S(t[i-1]+Dt,Z[i-1]))*D[i-1]
}     
R
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),RK,T=T,t0=t0,R0=R0,Dt=Dt,K=K,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Runge-Kutta \nscheme Order3")),line=2.8,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_1_RK <- data.frame(time,R)
if (M >=2) {RadialP_1_RK  <- data.frame(time,R,R.mean)}
showData(RadialP_1_RK, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_RK, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_RK)
}
}

.RadialP_1 <-
function(N,t0,Dt,T=1,R0,K,Sigma,Output=FALSE,
                    Methods = c("Euler", "Milstein", "MilsteinS", 
                                "Ito-Taylor", "Heun", "RK3"), ...)
{

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))


if( R0 < 0 )
            stop(tkmessageBox(title="Error",message=paste( "R0 > 0" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( Sigma < 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K < 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

Methods <- match.arg(Methods)

if(Methods=="Euler"){
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]}      
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Euler scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_1_Euler <- data.frame(time,R)
showData(RadialP_1_Euler, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_Euler, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_Euler)
}


if(Methods=="Milstein")
                  {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

DSx  <- D(diff,"x")
A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
            R[i] = R[i-1]+ A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]+ 
                   0.5 *S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1])^2 -Dt)}     
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Milstein scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_1_Milstein <- data.frame(time,R)
showData(RadialP_1_Milstein, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_Milstein, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_Milstein)
}

if(Methods=="MilsteinS")
                  {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
             R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1] +
                    0.5 *S(t[i-1],R[i-1]) * Sx(t[i-1],R[i-1])*(D[i-1]^2-Dt)+ 
                    Dt^(3 /2)*(0.5 *A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1]) +
                    0.5 *Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])+
                    0.25 *(S(t[i -1] ,R[i -1])^2) * Sxx(t[i -1] ,R[i -1]))*D[i -1]+ 
                   (Dt^2) * (0.5*A(t[i -1],R[i -1])*Ax(t[i-1],R[i-1])+
                    0.25 *Axx(t[i-1],R[i-1])*(S(t[i-1],R[i-1])^2))}     
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Second Milstein scheme")),line=3.3,adj=1,cex=0.6,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_1_MilsteinS <- data.frame(time,R)
showData(RadialP_1_MilsteinS, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_MilsteinS, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_MilsteinS)
}

if(Methods=="Ito-Taylor")
                  {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
DZ= rnorm(N,0,sqrt((1/3)*Dt^3))
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
                 R[i]=R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]+
                      0.5*S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1]^2)-Dt)+
                      Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])*DZ[i-1]+0.5*(A(t[i-1],R[i-1])*Ax(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Axx(t[i-1],R[i-1]))*(Dt^2)+(A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Sxx(t[i-1],R[i-1]))*(D[i-1]*Dt-DZ[i-1])+
                      0.5*S(t[i-1],R[i-1])*(S(t[i-1],R[i-1])*Sxx(t[i-1],R[i-1])+
                     (Sx(t[i-1],R[i-1])^2))*((1/3)*(D[i-1]^2)-Dt)*D[i-1]}     
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Strong Ito-Taylor \nScheme Order 1.5")),line=3,adj=1,cex=0.6,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_1_ITY <- data.frame(time,R)
showData(RadialP_1_ITY, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_ITY, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_ITY)
}

if(Methods=="Heun")
                  {

drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
         Y[i-1]= R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]
         R[i]  = R[i-1]+0.5*Dt*(A(t[i-1],R[i-1])+A(t[i-1],Y[i-1]))+
                 0.5*(S(t[i-1],R[i-1])+S(t[i-1],Y[i-1]))*D[i-1]}     
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Heun scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_1_Heun <- data.frame(time,R)
showData(RadialP_1_Heun, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_Heun, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_Heun)
}

if(Methods=="RK3")
                  {
drif <- expression(  (0.5*Sigma^2 - K )/ x )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
              Y    <- numeric()
              Z    <- numeric()
              Y[i-1]=R[i-1]+0.5*Dt*A(t[i-1],R[i-1])+S(t[i-1],R[i-1])*D[i-1]
              Z[i-1]=R[i-1]-A(t[i-1],R[i-1])*Dt+2*Dt*A(t[i-1]+0.5*Dt,Y[i-1])+
                     (2*S(t[i-1]+0.5*Dt,Y[i-1])-S(t[i-1],R[i-1]))*D[i-1]
              R[i] = R[i-1]+(Dt/6)*(A(t[i-1],R[i-1])+4*A(t[i-1]+0.5*Dt,Y[i-1])+A(t[i-1]+Dt,Z[i-1]))+
                     (1/6)*(S(t[i-1],R[i-1])+4*S(t[i-1]+0.5*Dt,Y[i-1])+S(t[i-1]+Dt,Z[i-1]))*D[i-1]}  
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s==1,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Runge-Kutta \nscheme Order3")),line=2.8,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac((sigma^2/2) -K,R[t])*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== 1),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_1_RK <- data.frame(time,R)
showData(RadialP_1_RK, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_1_RK, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_1_RK)
}
}

.RadialP2D_1 <-
function(N,t0,Dt,T=1,X0,Y0,v,K,Sigma,Output=FALSE)
       {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if( X0 == 0 & Y0 == 0 )
            stop(tkmessageBox(title="Error",message=paste( "X0 =! 0 or Y0 =! 0" ),icon="error"))

if( v <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "v > 0" ),icon="error"))

if ((X0^2) + (Y0^2) < v^2 || (X0^2) + (Y0^2) == v^2)
            stop(tkmessageBox(title="Error",message=paste( "X0^2 + Y0 ^2 > v^2 > 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( Sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))


drifx     <- expression( (-K*x) / (x^2+y^2) )
drify     <- expression( (-K*y) / (x^2+y^2) )
diff      <- expression( Sigma ) 

Ax    <- function(t,x,y)  eval(drifx)
Ay    <- function(t,x,y)  eval(drify)
S     <- function(t,x,y)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}

Dt= (T-t0)/N 
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx   <- diff(wx)

wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)

X    <- numeric()
Y    <- numeric()
X[1] <- X0
Y[1] <- Y0
for (i in 2:(N+1)){ 
        X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + S(t[i-1],X[i-1],Y[i-1])*Dx[i-1]
        Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + S(t[i-1],X[i-1],Y[i-1])*Dy[i-1]
                   } 
R = sqrt((X^2) + (Y^2))
n = which(R <= v)
if (length(n) > 0){
X    <- X[c(seq(1,min(n),by=1))]
Y    <- Y[c(seq(1,min(n),by=1))]
time <- t[c(seq(1,min(n),by=1))]}else{X <- X ; Y <- Y ; time <- t}

plot(X,Y,las=1,type="n",xlim=c(pmin(0,min(X)),max(X)),ylim=c(pmin(0,min(Y)),max(Y)),xlab=expression(X[t]),ylab=expression(Y[t]))
theta = seq(0,2*pi,length=1000)
a = v *cos(theta)
b = v *sin(theta)
points(a,b,type="l",col="gray")
polygon(a, b, col="gray")

if (X0 >= 0 || Y0 >= 0){
arrows(0,0,v*cos(pi/4),v*sin(pi/4), col= "black",code = 2,length = 0.1, angle = 20,lwd=2)
segments(0,0,v*cos(pi/4) ,v*sin(pi/4),col= "black",lwd=2)
points(0,0,type="p",pch=20,col="black",cex=1.6)
text((1/2)*v*cos(pi/4),(1/1.5)*v*sin(pi/4), "v", col="black", adj=c(-.1,-.1))}

if (X0 < 0 || Y0 < 0){
arrows(0,0,v*cos((7*pi)/4),v*sin((7*pi)/4), col= "black",code = 2,length = 0.1, angle = 20,lwd=2)
segments(0,0,v*cos((7*pi)/4) ,v*sin((7*pi)/4),col= "black",lwd=2)
points(0,0,type="p",pch=20,col="black",cex=1.6)
text((1/2)*v*cos((7*pi)/4),(1/1.5)*v*sin((7*pi)/4), "v", col="black", adj=c(-.1,-.1))}

mtext(expression("Simulation 2-Dimensional Attractive Model":bolditalic(M[(list(s==1,sigma))])),line=3.1,adj=0.8,cex=0.85,col="black")
mtext(bquote(dX[t]== frac(-K*X[t],X[t]^{2}+Y[t]^{2})*dt + sigma *d*W[t]^1),cex=0.7,adj=0,line=2.1,col="blue")
mtext(bquote(dY[t]== frac(-K*Y[t],X[t]^{2}+Y[t]^{2})*dt + sigma *d*W[t]^2),cex=0.7,adj=0,line=0.2,col="blue")

abline(h=0, v=0, col = "gray60",lwd=2)
mtext(bquote(X[t[0]]==.(X0)),line=1.4,adj=0.85,cex=0.9,col="green4")
mtext(bquote(Y[t[0]]==.(Y0)),line=0.4,adj=0.85,cex=0.9,col="green4")
mtext(bquote(K== .(K)),line=1.6,adj=1,cex=0.9,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.6,adj=1,cex=0.9,col="blue")
mtext(bquote(v== .(v)),line=0.7,adj=0.7,cex=0.9,col="black")

points(X0,Y0,type="p",pch=20,col="green4",cex=1.8)
text(X0,Y0, expression((list(X[t[0]],Y[t[0]]))), col="green4", adj=c(.5,-.2),cex = 0.8)

if (length(n) > 0 ){
for (i in 1:min(n)){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="blue",lwd=2)}}else{
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="blue",lwd=2)}}

if (length(n) > 0 ){points(X[min(n)],Y[min(n)],type="p",col="red",cex=1.2,pch="*")
                    text(X[min(n)],Y[min(n)], expression(tau[v]^(1)), col=2, adj=c(-.1,-.1),cex = 1.2)
                    mtext(bquote(tau[v]^(1)== .(t[min(n)])),line=0.5,adj=0.45,cex=1,col="red")}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
Result <- data.frame(time,X,Y)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.RadialP2D_1PC <-
function(N,R0,t0,T,ThetaMax,K,sigma,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( R0 <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "R0 > 0" ),icon="error"))


if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if( 2 * K <=  sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > sigma^2" ),icon="error"))

if( sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "sigma > 0" ),icon="error"))

if( K <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))



drif     <- expression((((sigma^2)/2)-K)/r)
diff     <- expression( sigma ) 

DAx  <- D(drif,"r")
DAxx <- D(D(drif,"r"),"r")
DSx  <- D(diff,"r")
DSxx <- D(D(diff,"r"),"r")
A    <- function(t,r)  eval(drif)
Ax   <- function(t,r)  eval(DAx)
Axx  <- function(t,r)  eval(DAxx)
S    <- function(t,r)  eval(diff)
Sx   <- function(t,r)  eval(DSx)
Sxx  <- function(t,r)  eval(DSxx)

t = seq(t0,T,length=N+1)
Dt= (T-t0)/N 
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
                  R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1] +
                         0.5 *S(t[i-1],R[i-1]) * Sx(t[i-1],R[i-1])*(D[i-1]^2-Dt)+ 
                         Dt^(3 /2)*(0.5 *A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1]) +
                         0.5 *Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])+
                         0.25 *(S(t[i -1] ,R[i -1])^2) * Sxx(t[i -1] ,R[i -1]))*D[i -1]+ 
                         (Dt^2) * (0.5*A(t[i -1],R[i -1])*Ax(t[i-1],R[i-1])+
                         0.25 *Axx(t[i-1],R[i-1])*(S(t[i-1],R[i-1])^2))
                   }      
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
theta = seq(0,ThetaMax,length=length(R))
X <- R * cos(theta)
Y <- R * sin(theta)
plot(X,Y,type="n",xlim=c(min(X),max(X)),ylim=c(min(Y),max(Y)),xlab="",ylab="",las=1)
mtext(expression("Simulation Two-Dimensional Attractive Model \nModels M(1,sigma) in Polar coordinates"),line=2.4,adj=0.85,cex=0.9,col="black")
mtext(bquote(dX[t]== frac(-K*X[t],X[t]^2+Y[t]^2)*dt + sigma *d*W[t]^1),cex=0.7,adj=0,line=2.1,col="blue")
mtext(bquote(dY[t]== frac(-K*Y[t],X[t]^2+Y[t]^2)*dt + sigma *d*W[t]^2),cex=0.7,adj=0,line=0,col="blue")
abline(h=0, v=0, col = "gray60",lwd=2)
mtext(bquote(R[t[0]]==.(R0)),line=1.5,adj=0.75,cex=0.9,col="blue")
mtext(bquote(theta[T]== .(round(ThetaMax,2))),line=0.6,adj=0.75,cex=0.9,col="blue")
mtext(bquote(K== .(K)),line=1.8,adj=1,cex=0.9,col="blue")
mtext(bquote(sigma== .(sigma)),line=0.8,adj=1,cex=0.9,col="blue")
mtext(paste("Polar coordinates"),side = 1, line = 2, adj = 0, cex = .8,col="blue")
points(R0,0,type="p",pch=20,col="blue",cex=1.8)
text(R0,0, expression((list(R[t[0]],theta[t[0]]))), col="blue", adj=c(.5,-.2),cex = 1)
for (i in 1:length(R)){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="blue",lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
Result <- data.frame(time,theta,R)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.RadialP3D_1 <-
function(N,t0,Dt,T=1,X0,Y0,Z0,v,K,Sigma,Output=FALSE)
    {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( X0 == 0 & Y0 == 0 & Z0 == 0)
            stop(tkmessageBox(title="Error",message=paste( "X0 =! 0 or Y0 =! 0 or Z0 =!" ),icon="error"))

if( v <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "v > 0" ),icon="error"))

if ((X0^2) + (Y0^2) + (Z0^2) < v^2 || (X0^2) + (Y0^2) + (Z0^2) == v^2)
            stop(tkmessageBox(title="Error",message=paste( "X0^2 + Y0^2 + Z0^2 > v^2" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( Sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

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
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)

wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)

wz = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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


Result <- data.frame(time,X,Y,Z)
if (length(n) > 0 ){points3d(X[min(n)],Y[min(n)],Z[min(n)],col="blue",size=8)
                    text3d(X[min(n)],Y[min(n)],Z[min(n)],texts=c("FPT = "),adj=c(0.5,-0.8),color = c("blue"),cex=1.2,family=c("serif"))
                    text3d(X[min(n)],Y[min(n)],Z[min(n)],texts=c(t[min(n)]),adj=c(-0.9,-0.8),color = c("blue"),cex=1.2,family=c("serif"))
                    text3d(V,V,V,c("FPT : First Passage Time"),adj=c(0.5,-0.25),cex=1.2,col="blue",family=c("serif"))
                    }
showData(Result , placement='+200-200', font = "vourier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "RadialP3D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
attach(Result)
}

.tho_M1 <-
function(N,M,t0,T,R0,v,K,sigma,Output = FALSE,
                   Methods = c("Euler", "Milstein", "MilsteinS", 
                               "Ito-Taylor", "Heun", "RK3"), ...)
       {
FPT <- numeric()
i = 1 
while( i <= M) { FPT[i] <- .tho_1(N,t0,T,R0,v,K,sigma,Methods)
                if ( !is.na(FPT[i]) ) {i = i +1} 
                       }
thoM1 <- data.frame(FPT)
showData(thoM1, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(thoM1, file = "FPT.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }                                    
attach(thoM1)
}


.RadialP_2flow <-
function(N,M,t0,Dt,T=1,R0,K,s,Sigma,Output=FALSE,
                    Methods = c("Euler", "Milstein", "MilsteinS", 
                                "Ito-Taylor", "Heun", "RK3"), ...)
{

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be >= 1" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))


if( R0 < 0 )
            stop(tkmessageBox(title="Error",message=paste( "R0 > 0" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( Sigma < 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K < 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

if( s <= 1 )
            stop(tkmessageBox(title="Error",message=paste( "s >= 2" ),icon="error"))


Methods <- match.arg(Methods)

if(Methods=="Euler")
                  {
Eul <- function(N,T=1,Dt,t0,R0,K,s,Sigma)
     {
drif <- expression( ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]}      
R
   }
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Eul,T=T,t0=t0,R0=R0,Dt=Dt,K=K,s=s,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Euler scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_2_Euler <- data.frame(time,R)
if (M >=2) {RadialP_2_Euler  <- data.frame(time,R,R.mean)}
showData(RadialP_2_Euler, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_Euler, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_2_Euler)
}

if(Methods=="Milstein")
                  {
Mil <- function(N,T=1,Dt,t0,R0,K,s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

DSx  <- D(diff,"x")
A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
            R[i] = R[i-1]+ A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]+ 
                   0.5 *S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1])^2 -Dt)}     
R
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Mil,T=T,t0=t0,R0=R0,Dt=Dt,K=K,s=s,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Milstein scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_2_Milstein <- data.frame(time,R)
if (M >=2) {RadialP_2_Milstein  <- data.frame(time,R,R.mean)}
showData(RadialP_2_Milstein, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_Milstein, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_2_Milstein)
}

if(Methods=="MilsteinS")
                  {
MilS <- function(N,T=1,Dt,t0,R0,K,s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
             R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1] +
                    0.5 *S(t[i-1],R[i-1]) * Sx(t[i-1],R[i-1])*(D[i-1]^2-Dt)+ 
                    Dt^(3 /2)*(0.5 *A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1]) +
                    0.5 *Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])+
                    0.25 *(S(t[i -1] ,R[i -1])^2) * Sxx(t[i -1] ,R[i -1]))*D[i -1]+ 
                   (Dt^2) * (0.5*A(t[i -1],R[i -1])*Ax(t[i-1],R[i-1])+
                    0.25 *Axx(t[i-1],R[i-1])*(S(t[i-1],R[i-1])^2))}     
R
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),MilS,T=T,t0=t0,R0=R0,Dt=Dt,K=K,s=s,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Second Milstein scheme")),line=3.3,adj=1,cex=0.6,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_2_MilsteinS <- data.frame(time,R)
if (M >=2) {RadialP_2_MilsteinS  <- data.frame(time,R,R.mean)}
showData(RadialP_2_MilsteinS, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_MilsteinS, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
(RadialP_2_MilsteinS)
}

if(Methods=="Ito-Taylor")
                  {
ITO <- function(N,T=1,Dt,t0,R0,K,s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
DZ= rnorm(N,0,sqrt((1/3)*Dt^3))
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
                 R[i]=R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]+
                      0.5*S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1]^2)-Dt)+
                      Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])*DZ[i-1]+0.5*(A(t[i-1],R[i-1])*Ax(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Axx(t[i-1],R[i-1]))*(Dt^2)+(A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Sxx(t[i-1],R[i-1]))*(D[i-1]*Dt-DZ[i-1])+
                      0.5*S(t[i-1],R[i-1])*(S(t[i-1],R[i-1])*Sxx(t[i-1],R[i-1])+
                     (Sx(t[i-1],R[i-1])^2))*((1/3)*(D[i-1]^2)-Dt)*D[i-1]}     
R
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),ITO,T=T,t0=t0,R0=R0,Dt=Dt,K=K,s=s,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Strong Ito-Taylor \nScheme Order 1.5")),line=3,adj=1,cex=0.6,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_2_ITY <- data.frame(time,R)
if (M >=2) {RadialP_2_ITY  <- data.frame(time,R,R.mean)}
showData(RadialP_2_ITY, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_ITY, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_2_ITY)
}

if(Methods=="Heun")
                  {
Heu <- function(N,T=1,Dt,t0,R0,K,s=s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
         Y[i-1]= R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]
         R[i]  = R[i-1]+0.5*Dt*(A(t[i-1],R[i-1])+A(t[i-1],Y[i-1]))+
                 0.5*(S(t[i-1],R[i-1])+S(t[i-1],Y[i-1]))*D[i-1]}     
R
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),Heu,T=T,t0=t0,R0=R0,Dt=Dt,K=K,s=s,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Heun scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_2_Heun <- data.frame(time,R)
if (M >=2) {RadialP_2_Heun  <- data.frame(time,R,R.mean)}
showData(RadialP_2_Heun, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_Heun, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_2_Heun)
}

if(Methods=="RK3")
                  {
RK <- function(N,T=1,Dt,t0,R0,K,s=s,Sigma)
     {
drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
              Y    <- numeric()
              Z    <- numeric()
              Y[i-1]=R[i-1]+0.5*Dt*A(t[i-1],R[i-1])+S(t[i-1],R[i-1])*D[i-1]
              Z[i-1]=R[i-1]-A(t[i-1],R[i-1])*Dt+2*Dt*A(t[i-1]+0.5*Dt,Y[i-1])+
                     (2*S(t[i-1]+0.5*Dt,Y[i-1])-S(t[i-1],R[i-1]))*D[i-1]
              R[i] = R[i-1]+(Dt/6)*(A(t[i-1],R[i-1])+4*A(t[i-1]+0.5*Dt,Y[i-1])+A(t[i-1]+Dt,Z[i-1]))+
                     (1/6)*(S(t[i-1],R[i-1])+4*S(t[i-1]+0.5*Dt,Y[i-1])+S(t[i-1]+Dt,Z[i-1]))*D[i-1]
}     
R
}
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
Q = sapply(rep(N,length=M),RK,T=T,t0=t0,R0=R0,Dt=Dt,K=K,s=s,Sigma=Sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(R[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Runge-Kutta \nscheme Order3")),line=2.8,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topright",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
R.mean <- Q.mean
R <- Q
RadialP_2_RK <- data.frame(time,R)
if (M >=2) {RadialP_2_RK  <- data.frame(time,R,R.mean)}
showData(RadialP_2_RK, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_RK, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_2_RK)
}
}


.RadialP_2 <-
function(N,t0,Dt,T=1,R0,K,s,Sigma,Output=FALSE,
                    Methods = c("Euler", "Milstein", "MilsteinS", 
                                "Ito-Taylor", "Heun", "RK3"), ...)
{

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))


if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))


if( R0 < 0 )
            stop(tkmessageBox(title="Error",message=paste( "R0 > 0" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( Sigma < 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K < 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

if( s <= 1 )
            stop(tkmessageBox(title="Error",message=paste( "s >= 2" ),icon="error"))

Methods <- match.arg(Methods)

if(Methods=="Euler")
                  {
drif <- expression( ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]}      
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Euler scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_2_Euler <- data.frame(time,R)
showData(RadialP_2_Euler, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_Euler, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_2_Euler)
}

if(Methods=="Milstein")
                  {

drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

DSx  <- D(diff,"x")
A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
            R[i] = R[i-1]+ A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1]+ 
                   0.5 *S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1])^2 -Dt)}     
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Milstein scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_2_Milstein <- data.frame(time,R)
showData(RadialP_2_Milstein, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_Milstein, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_2_Milstein)
}

if(Methods=="MilsteinS")
                  {

drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
             R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1] +
                    0.5 *S(t[i-1],R[i-1]) * Sx(t[i-1],R[i-1])*(D[i-1]^2-Dt)+ 
                    Dt^(3 /2)*(0.5 *A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1]) +
                    0.5 *Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])+
                    0.25 *(S(t[i -1] ,R[i -1])^2) * Sxx(t[i -1] ,R[i -1]))*D[i -1]+ 
                   (Dt^2) * (0.5*A(t[i -1],R[i -1])*Ax(t[i-1],R[i-1])+
                    0.25 *Axx(t[i-1],R[i-1])*(S(t[i-1],R[i-1])^2))}     
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Second Milstein scheme")),line=3.3,adj=1,cex=0.6,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_2_MilsteinS <- data.frame(time,R)
showData(RadialP_2_MilsteinS, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_MilsteinS, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
(RadialP_2_MilsteinS)
}

if(Methods=="Ito-Taylor")
                  {

drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

DAx  <- D(drif,"x")
DAxx <- D(D(drif,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drif)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
S    <- function(t,x)  eval(diff)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
DZ= rnorm(N,0,sqrt((1/3)*Dt^3))
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
                 R[i]=R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]+
                      0.5*S(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])*((D[i-1]^2)-Dt)+
                      Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])*DZ[i-1]+0.5*(A(t[i-1],R[i-1])*Ax(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Axx(t[i-1],R[i-1]))*(Dt^2)+(A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1])+
                      0.5*(S(t[i-1],R[i-1])^2)*Sxx(t[i-1],R[i-1]))*(D[i-1]*Dt-DZ[i-1])+
                      0.5*S(t[i-1],R[i-1])*(S(t[i-1],R[i-1])*Sxx(t[i-1],R[i-1])+
                     (Sx(t[i-1],R[i-1])^2))*((1/3)*(D[i-1]^2)-Dt)*D[i-1]}     
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Strong Ito-Taylor \nScheme Order 1.5")),line=3,adj=1,cex=0.6,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_2_ITY <- data.frame(time,R)
showData(RadialP_2_ITY, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_ITY, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_2_ITY)
}

if(Methods=="Heun")
                  {

drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
         Y[i-1]= R[i-1]+A(t[i-1],R[i-1])*Dt+S(t[i-1],R[i-1])*D[i-1]
         R[i]  = R[i-1]+0.5*Dt*(A(t[i-1],R[i-1])+A(t[i-1],Y[i-1]))+
                 0.5*(S(t[i-1],R[i-1])+S(t[i-1],Y[i-1]))*D[i-1]}     
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l")
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Heun scheme")),line=3.3,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_2_Heun <- data.frame(time,R)
showData(RadialP_2_Heun, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_Heun, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_2_Heun)
}

if(Methods=="RK3")
                  {

drif <- expression(  ( (0.5*Sigma^2)*(x^(s-1)) - K ) / (x^s) )
diff <- expression( Sigma )

A    <- function(t,x)  eval(drif)
S    <- function(t,x)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
Y    <- numeric()
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
              Y    <- numeric()
              Z    <- numeric()
              Y[i-1]=R[i-1]+0.5*Dt*A(t[i-1],R[i-1])+S(t[i-1],R[i-1])*D[i-1]
              Z[i-1]=R[i-1]-A(t[i-1],R[i-1])*Dt+2*Dt*A(t[i-1]+0.5*Dt,Y[i-1])+
                     (2*S(t[i-1]+0.5*Dt,Y[i-1])-S(t[i-1],R[i-1]))*D[i-1]
              R[i] = R[i-1]+(Dt/6)*(A(t[i-1],R[i-1])+4*A(t[i-1]+0.5*Dt,Y[i-1])+A(t[i-1]+Dt,Z[i-1]))+
                     (1/6)*(S(t[i-1],R[i-1])+4*S(t[i-1]+0.5*Dt,Y[i-1])+S(t[i-1]+Dt,Z[i-1]))*D[i-1]}    
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
plot(t,R,type="n",ylab=expression(R[t]),xlab="time",las=1)
points(t,R,type="l") 
mtext(bquote("Attractive Model": bolditalic(M[(list(s >= 2,sigma))])),line=2.8,adj=0.5,cex=1.1,col="black")
mtext(bquote(bolditalic("Runge-Kutta \nscheme Order3")),line=2.8,adj=1,cex=0.7,col="red")
mtext(bquote(dR[t]== frac(frac(sigma^2,2) *R[t]^(s-1) -K,R[t]^s)*dt + sigma *d*tilde(W)[t]),cex=1,adj=0,line=0.2,col="red")
mtext(bquote(S== .(s)),line=1.8,adj=0.75,cex=0.8,col="blue")
mtext(bquote(K== .(K)),line=1,adj=0.75,cex=0.8,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.2,adj=0.75,cex=0.8,col="blue")
mtext(bquote(R[0]==.(R0)),line=1.7,adj=1,cex=0.8,col="blue")
mtext(bquote(Delta*t==.(Dt)),line=0.2,cex=0.8,adj=1,col="blue")
mtext(bquote(T==.(T)),line=1,cex=0.8,adj=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
RadialP_2_RK <- data.frame(time,R)
showData(RadialP_2_RK, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(RadialP_2_RK, file = "RadialP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(RadialP_2_RK)
}
}

.RadialP2D_2 <-
function(N,t0,Dt,T=1,X0,Y0,v,K,s,Sigma,Output=FALSE)
       {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if( X0 == 0 & Y0 == 0 )
            stop(tkmessageBox(title="Error",message=paste( "X0 =! 0 or Y0 =! 0" ),icon="error"))

if( v <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "v > 0" ),icon="error"))

if ((X0^2) + (Y0^2) < v^2 || (X0^2) + (Y0^2) == v^2)
            stop(tkmessageBox(title="Error",message=paste( "X0^2 + Y0 ^2 > v^2 > 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( s <= 1 )
            stop(tkmessageBox(title="Error",message=paste( "s > 1" ),icon="error"))

if( Sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))


drifx     <- expression( (-K*x) / (sqrt( x^2 + y^2))^(s+1) )
drify     <- expression( (-K*y) / (sqrt( x^2 + y^2))^(s+1) )
diff      <- expression( Sigma ) 

Ax    <- function(t,x,y)  eval(drifx)
Ay    <- function(t,x,y)  eval(drify)
S     <- function(t,x,y)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}

Dt= (T-t0)/N 
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx   <- diff(wx)

wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)

X    <- numeric()
Y    <- numeric()
X[1] <- X0
Y[1] <- Y0
for (i in 2:(N+1)){ 
        X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + S(t[i-1],X[i-1],Y[i-1])*Dx[i-1]
        Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + S(t[i-1],X[i-1],Y[i-1])*Dy[i-1]
                   } 
R = sqrt((X^2) + (Y^2))
n = which(R <= v)
if (length(n) > 0){
X    <- X[c(seq(1,min(n),by=1))]
Y    <- Y[c(seq(1,min(n),by=1))]
time <- t[c(seq(1,min(n),by=1))]}else{X <- X ; Y <- Y ; time <- t}

plot(X,Y,las=1,type="n",xlim=c(pmin(0,min(X)),max(X)),ylim=c(pmin(0,min(Y)),max(Y)),xlab=expression(X[t]),ylab=expression(Y[t]))
theta = seq(0,2*pi,length=1000)
a = v *cos(theta)
b = v *sin(theta)
points(a,b,type="l",col="gray")
polygon(a, b, col="gray")

if (X0 >= 0 || Y0 >= 0){
arrows(0,0,v*cos(pi/4),v*sin(pi/4), col= "black",code = 2,length = 0.1, angle = 20,lwd=2)
segments(0,0,v*cos(pi/4) ,v*sin(pi/4),col= "black",lwd=2)
points(0,0,type="p",pch=20,col="black",cex=1.6)
text((1/2)*v*cos(pi/4),(1/1.5)*v*sin(pi/4), "v", col="black", adj=c(-.1,-.1))}

if (X0 < 0 || Y0 < 0){
arrows(0,0,v*cos((7*pi)/4),v*sin((7*pi)/4), col= "black",code = 2,length = 0.1, angle = 20,lwd=2)
segments(0,0,v*cos((7*pi)/4) ,v*sin((7*pi)/4),col= "black",lwd=2)
points(0,0,type="p",pch=20,col="black",cex=1.6)
text((1/2)*v*cos((7*pi)/4),(1/1.5)*v*sin((7*pi)/4), "v", col="black", adj=c(-.1,-.1))}

mtext(expression("Simulation 2-Dimensional Attractive Model":bolditalic(M[(list(s>=2,sigma))])),line=3.1,adj=0.93,cex=0.8,col="black")
mtext(bquote(dX[t]== frac(-K*X[t],(sqrt(X[t]^2+Y[t]^2))^{S+1})*dt + sigma *d*W[t]^1),cex=0.7,adj=0,line=2.1,col="blue")
mtext(bquote(dY[t]== frac(-K*Y[t],(sqrt(X[t]^2+Y[t]^2))^{S+1})*dt + sigma *d*W[t]^2),cex=0.7,adj=0,line=0.2,col="blue")

abline(h=0, v=0, col = "gray60",lwd=2)
mtext(bquote(X[t[0]]==.(X0)),line=1.4,adj=0.85,cex=0.9,col="green4")
mtext(bquote(Y[t[0]]==.(Y0)),line=0.4,adj=0.85,cex=0.9,col="green4")
mtext(bquote(K== .(K)),line=1.6,adj=1,cex=0.9,col="blue")
mtext(bquote(sigma== .(Sigma)),line=0.6,adj=1,cex=0.9,col="blue")
mtext(bquote(s== .(s)),line=1.6,adj=0.7,cex=0.9,col="blue")
mtext(bquote(v== .(v)),line=0.7,adj=0.7,cex=0.9,col="black")

points(X0,Y0,type="p",pch=20,col="green4",cex=1.8)
text(X0,Y0, expression((list(X[t[0]],Y[t[0]]))), col="green4", adj=c(.5,-.2),cex = 0.8)

if (length(n) > 0 ){
for (i in 1:min(n)){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="blue",lwd=2)}}else{
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="blue",lwd=2)}}

if (length(n) > 0 ){points(X[min(n)],Y[min(n)],type="p",col="red",cex=1.2,pch="*")
                    text(X[min(n)],Y[min(n)], expression(tau[v]^(s)), col=2, adj=c(-.1,-.1),cex = 1.2)
                    mtext(bquote(tau[v]^(s)== .(t[min(n)])),line=0.5,adj=0.45,cex=1,col="red")}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
Result <- data.frame(time,X,Y)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "Models2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.RadialP2D_2PC <-
function(N,R0,t0,T,ThetaMax,K,s,sigma,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( R0 <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "R0 > 0" ),icon="error"))


if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if( 2 * K <=  sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > sigma^2" ),icon="error"))

if( sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "sigma > 0" ),icon="error"))

if( K <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

if( s <= 1 )
            stop(tkmessageBox(title="Error",message=paste( "s > 1" ),icon="error"))



drif     <- expression(( (0.5*sigma^2)*(r^(s-1)) - K ) / (r^s) )
diff     <- expression( sigma ) 

DAx  <- D(drif,"r")
DAxx <- D(D(drif,"r"),"r")
DSx  <- D(diff,"r")
DSxx <- D(D(diff,"r"),"r")
A    <- function(t,r)  eval(drif)
Ax   <- function(t,r)  eval(DAx)
Axx  <- function(t,r)  eval(DAxx)
S    <- function(t,r)  eval(diff)
Sx   <- function(t,r)  eval(DSx)
Sxx  <- function(t,r)  eval(DSxx)

t = seq(t0,T,length=N+1)
Dt= (T-t0)/N 
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
D    <- diff(w)
R    <- numeric()
R[1] <- R0
for (i in 2:(N+1)){
                  R[i] = R[i-1] + A(t[i-1],R[i-1])*Dt + S(t[i-1],R[i-1])*D[i-1] +
                         0.5 *S(t[i-1],R[i-1]) * Sx(t[i-1],R[i-1])*(D[i-1]^2-Dt)+ 
                         Dt^(3 /2)*(0.5 *A(t[i-1],R[i-1])*Sx(t[i-1],R[i-1]) +
                         0.5 *Ax(t[i-1],R[i-1])*S(t[i-1],R[i-1])+
                         0.25 *(S(t[i -1] ,R[i -1])^2) * Sxx(t[i -1] ,R[i -1]))*D[i -1]+ 
                         (Dt^2) * (0.5*A(t[i -1],R[i -1])*Ax(t[i-1],R[i-1])+
                         0.25 *Axx(t[i-1],R[i-1])*(S(t[i-1],R[i-1])^2))
                   }      
NN <- which(R <= 0)
if(length(NN)>0){nn <- min(NN)}
if(length(NN)> 0){R <- R[seq(1,nn,by=1)]}
if(length(NN)> 0){t <- t[seq(1,nn,by=1)]}
theta = seq(0,ThetaMax,length=length(R))
X <- R * cos(theta)
Y <- R * sin(theta)
plot(X,Y,type="n",xlim=c(min(X),max(X)),ylim=c(min(Y),max(Y)),xlab="",ylab="",las=1)
mtext(expression("Simulation Two-Dimensional Attractive Model \nModels M(S>=2,sigma) in Polar coordinates"),line=2.4,adj=0.9,cex=0.8,col="black")
mtext(bquote(dX[t]== frac(-K*X[t],(sqrt(X[t]^2+Y[t]^2))^{S+1})*dt + sigma *d*W[t]^1),cex=0.7,adj=0,line=2.1,col="blue")
mtext(bquote(dY[t]== frac(-K*Y[t],(sqrt(X[t]^2+Y[t]^2))^{S+1})*dt + sigma *d*W[t]^2),cex=0.7,adj=0,line=0,col="blue")
abline(h=0, v=0, col = "gray60",lwd=2)
mtext(bquote(R[t[0]]==.(R0)),line=1.5,adj=0.75,cex=0.9,col="blue")
mtext(bquote(theta[T]== .(round(ThetaMax,2))),line=0.6,adj=0.75,cex=0.9,col="blue")
mtext(bquote(K== .(K)),line=1.8,adj=1,cex=0.9,col="blue")
mtext(bquote(sigma== .(sigma)),line=1,adj=1,cex=0.9,col="blue")
mtext(bquote(s== .(s)),line=0.2,adj=1,cex=0.9,col="blue")
mtext(paste("Polar coordinates"),side = 1, line = 2, adj = 0, cex = .8,col="blue")
points(R0,0,type="p",pch=20,col="blue",cex=1.8)
text(R0,0, expression((list(R[t[0]],theta[t[0]]))), col="blue", adj=c(.5,-.2),cex = 1)
for (i in 1:length(R)){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="blue",lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
Result <- data.frame(time,theta,R)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "Models2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.RadialP3D_2 <-
function(N,t0,Dt,T=1,X0,Y0,Z0,v,K,s,Sigma,Output=FALSE)
    {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( X0 == 0 & Y0 == 0 & Z0 == 0)
            stop(tkmessageBox(title="Error",message=paste( "X0 =! 0 or Y0 =! 0 or Z0 =!" ),icon="error"))

if( v <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "v > 0" ),icon="error"))

if ((X0^2) + (Y0^2) + (Z0^2) < v^2 || (X0^2) + (Y0^2) + (Z0^2) == v^2)
            stop(tkmessageBox(title="Error",message=paste( "X0^2 + Y0^2 + Z0^2 > v^2" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( s <= 1 )
            stop(tkmessageBox(title="Error",message=paste( "s > 1" ),icon="error"))

if( Sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

drifx     <- expression( (-K*x) / (sqrt(x^2 + y^2 + z^2))^(s+1) )
drify     <- expression( (-K*y) / (sqrt(x^2 + y^2 + z^2))^(s+1) )
drifz     <- expression( (-K*z) / (sqrt(x^2 + y^2 + z^2))^(s+1) )
diff      <- expression( Sigma ) 

Ax    <- function(t,x,y,z)  eval(drifx)
Ay    <- function(t,x,y,z)  eval(drify)
Az    <- function(t,x,y,z)  eval(drifz)
S     <- function(t,x,y,z)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt= (T-t0)/N 
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)

wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)

wz = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
title3d(family=c("serif"),main="Simulation Three-Dimensional Attractive Model M(S>=2,Sigma)",color = c("black"),cex=1.2)

if (length(n) > 0){
for (i in 1:min(n)) {lines3d(c(G[i,1],G[i+1,1]),c(G[i,2],G[i+1,2]),c(G[i,3],G[i+1,3]),col="red",from ="lines",lwd=2)}}else
{for (i in 1:N) {lines3d(c(G[i,1],G[i+1,1]),c(G[i,2],G[i+1,2]),c(G[i,3],G[i+1,3]),col="red",from ="lines",lwd=2)}}


Result <- data.frame(time,X,Y,Z)
if (length(n) > 0 ){points3d(X[min(n)],Y[min(n)],Z[min(n)],col="blue",size=8)
                    text3d(X[min(n)],Y[min(n)],Z[min(n)],texts=c("FPT = "),adj=c(0.5,-0.8),color = c("blue"),cex=1.2,family=c("serif"))
                    text3d(X[min(n)],Y[min(n)],Z[min(n)],texts=c(t[min(n)]),adj=c(-0.9,-0.8),color = c("blue"),cex=1.2,family=c("serif"))
                    text3d(V,V,V,c("FPT : First Passage Time"),adj=c(0.5,-0.25),cex=1.2,col="blue",family=c("serif"))
                    }
showData(Result , placement='+200-200', font = "vourier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "Models2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
attach(Result)
}

.tho_M2 <-
function(N,M,t0,T,R0,v,K,s,Sigma,Output = FALSE,
                   Methods = c("Euler", "Milstein", "MilsteinS", 
                               "Ito-Taylor", "Heun", "RK3"), ...)
       {
FPTT <- numeric()
i = 1 
while( i <= M) { FPTT[i] <- .tho_2(N,t0,T,R0,v,K,s,Sigma,Methods)
                if ( !is.na(FPTT[i]) ) {i = i +1} 
                       }
thoM2 <- data.frame(FPTT)
showData(thoM2, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(thoM2, file = "FPTModels2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }                                    
attach(thoM2)
}


.TwoDiffAtra2D <-
function(N,t0,Dt,T=1,X1_0,X2_0,Y1_0,Y2_0,v,K,m,Sigma,Output=FALSE)
            {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if( v <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "v > 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( Sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

if( m <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "m > 0" ),icon="error"))

Q = sqrt((X1_0-Y1_0)^2 + (X2_0-Y2_0)^2)

if( Q <= 0.01 )
            stop(tkmessageBox(title="Error",message=paste( "D1 = sqrt((X1_0-Y1_0)^2 + (X2_0-Y2_0)^2) > 0" ),icon="error"))


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
wx1 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx1     <- diff(wx1)

wx2 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx2     <- diff(wx2)

wy1 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy1     <- diff(wy1)

wy2 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
lines(c(Y1[i],Y1[i+1]),c(Y2[i],Y2[i+1]),type="l",col="blue",lwd=1)
lines(c(X1[i],X1[i+1]),c(X2[i],X2[i+1]),type="l",col="red",lwd=1)
                 }
n <- which(D <= v)
if (length(n) > 0) {thoV1V2 <- t[min(n)]
mtext(bquote( bolditalic( tau[group("||",D[t],"||")<=v]^(m)*(list(V[t]^(1),V[t]^(2)))==.(round(thoV1V2,4)) )  ),cex=0.9,adj=0.5,line=0.8,col="black")
##points(X1[length(D)],X2[length(D)],pch=19,col="black",cex=1.1)
##legend("topleft",bg="gray85",border="gray",expression(group("||",D[t],"||")<=v),col=c("black"),pch=19,cex=0.7)
points(c(X1[length(D)],Y1[length(D)]),c(X2[length(D)],Y2[length(D)]),type="l",col="black",cex=1.1,lwd=2)
legend("topleft",border="gray",expression(group("||",D[t],"||")<=v),col=c("black"),lty=1,cex=0.7,lwd=2)
}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
Result <- data.frame(t,X1,X2,Y1,Y2,D)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "2Diffforattraction.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.TwoDiffAtra3D <-
function(N,t0,Dt,T=1,X1_0,X2_0,X3_0,Y1_0,Y2_0,Y3_0,v,K,m,Sigma,Output=FALSE)
              {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if( v <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "v > 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if( 2 * K <=  Sigma^2 )
            stop(tkmessageBox(title="Error",message=paste( "2*K > Sigma^2" ),icon="error"))

if( Sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if( K <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "K > 0" ),icon="error"))

if( m <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "m > 0" ),icon="error"))

Q = sqrt((X1_0-Y1_0)^2 + (X2_0-Y2_0)^2 + (X3_0-Y3_0)^2)

if( Q <= 0.01 )
            stop(tkmessageBox(title="Error",message=paste( "D1 = sqrt((X1_0-Y1_0)^2 + (X2_0-Y2_0)^2 + (X3_0-Y3_0)^2) > 0" ),icon="error"))
			
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
wx1 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx1     <- diff(wx1)

wx2 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx2     <- diff(wx2)

wx3 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx3     <- diff(wx3)

wy1 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy1     <- diff(wy1)

wy2 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy2     <- diff(wy2)

wy3 = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
text3d(X1_0,X2_0,X3_0,c("(X1,X2,X3)"),adj=c(0.5,-0.25),cex=0.8,family=c("serif"))
points3d(Gx[1,],color = c("red"),size=6)
text3d(Y1_0,Y2_0,Y3_0,c("(Y1,Y2,Y3)"),adj=c(0.5,1),cex=0.8,family=c("serif"))
points3d(Gy[1,],color = c("blue"),size=6)
title3d(family=c("serif"),main="3-Dimensional Attractive Model for 2-Diffusion Processes",color = c("black"),cex=1.2)

if (length(n) > 0){
for (i in 1:min(n)) {lines3d(c(Gy[i,1],Gy[i+1,1]),c(Gy[i,2],Gy[i+1,2]),c(Gy[i,3],Gy[i+1,3]),col="blue",from ="lines",lwd=2)
                     lines3d(c(Gx[i,1],Gx[i+1,1]),c(Gx[i,2],Gx[i+1,2]),c(Gx[i,3],Gx[i+1,3]),col="red",from ="lines",lwd=2)}}else
{for (i in 1:N) {lines3d(c(Gy[i,1],Gy[i+1,1]),c(Gy[i,2],Gy[i+1,2]),c(Gy[i,3],Gy[i+1,3]),col="blue",from ="lines",lwd=2)
                 lines3d(c(Gx[i,1],Gx[i+1,1]),c(Gx[i,2],Gx[i+1,2]),c(Gx[i,3],Gx[i+1,3]),col="red",from ="lines",lwd=2)}}


if (length(n) > 0 ){points3d(X1[min(n)],X2[min(n)],X3[min(n)],col="green4",size=8)
                    text3d(X1[min(n)],X2[min(n)],X3[min(n)],texts=c("FPT = "),adj=c(0.5,-0.8),color = c("green4"),cex=1.2,family=c("serif"))
                    text3d(X1[min(n)],X2[min(n)],X3[min(n)],texts=c(t[min(n)]),adj=c(-0.9,-0.8),color = c("green4"),cex=1.2,family=c("serif"))
                    text3d(V,V,V,c("FPT : First Passage Time"),adj=c(0.5,-0.25),cex=1.2,col="green4",family=c("serif"))
                    }
Result <- data.frame(t,X1,X2,X3,Y1,Y2,Y3,D)
showData(Result , placement='+200-200', font = "vourier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "3Diffforattraction.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
attach(Result)
}

.tho_02diff <-
function(N,M,t0,Dt,T=1,X1_0,X2_0,Y1_0,Y2_0,v,K,m,Sigma,Output=FALSE)
       {
FPT <- numeric()
i = 1 
while( i <= M) { FPT[i] <-.Sim_tho02diff(N,t0,Dt,T=1,X1_0,X2_0,Y1_0,Y2_0,v,K,m,Sigma)
                if ( !is.na(FPT[i]) ) {i = i +1} 
                       }
tho02diff <- data.frame(FPT)
showData(tho02diff, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(tho02diff, file = "FPT2Diffforattraction.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(tho02diff)
}

.BMRW2D <- function(N,t0,T,x0,y0,Sigma,Step=FALSE,Output=FALSE)
       {
t = seq(t0,T,length=N+1)
dt = (T-t0)/N
u = runif(N,0,1)
o = rep(1,N)
o[ which( u < 0.5) ]= -1
w1 = cumsum(c(0,o))*sqrt(dt)
Dx    <- diff(w1)
u1 = runif(N,0,1)
o1 = rep(1,N)
o1[ which( u1 < 0.5) ]= -1
w2 = cumsum(c(0,o1))*sqrt(dt)
Dy    <- diff(w2)
drifx     <- expression( 0 )
drify     <- expression( 0 )
diff      <- expression( Sigma )
Ax    <- function(t,x,y)  eval(drifx)
Ay    <- function(t,x,y)  eval(drify)
S     <- function(t,x,y)  eval(diff) 
X <- numeric()
Y <- numeric()
X[1] <- x0
Y[1] <- y0
for (i in 2:(N+1)){     
        X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*dt + S(t[i-1],X[i-1],Y[i-1])*Dx[i-1]
        Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*dt + S(t[i-1],X[i-1],Y[i-1])*Dy[i-1]
                   } 
plot(X,Y,las=1,lwd=3,type="n",xlab=expression(W[t]),ylab=expression(W[t]))
mtext("Brownian motion in 2D plane",line=2,cex=1.2) 
mtext("by a Random Walk",line=0.5,cex=1.2,col="red")
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
if(Step==FALSE){points(X,Y,type="l",lwd=1,col="black")}
if(Step==TRUE){for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",
                col="black",lwd=1)}}
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
BMRW2_D <- data.frame(t,X,Y)
showData(BMRW2_D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(BMRW2_D, file = "BMRW2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(BMRW2_D)
}

.BMN2D <- function(N,t0,T,x0,y0,Sigma,Step=FALSE,Output=FALSE)
       {
t = seq(t0,T,length=N+1)
dt = (T-t0)/N
w1 = cumsum(rnorm(N+1,mean=0,sd=sqrt(dt)))
w2 = cumsum(rnorm(N+1,mean=0,sd=sqrt(dt)))
Dx    <- diff(w1)
Dy    <- diff(w2)
drifx     <- expression( 0 )
drify     <- expression( 0 )
diff      <- expression( Sigma )
Ax    <- function(t,x,y)  eval(drifx)
Ay    <- function(t,x,y)  eval(drify)
S     <- function(t,x,y)  eval(diff) 
X <- numeric()
Y <- numeric()
X[1] <- x0
Y[1] <- y0
for (i in 2:(N+1)){     
        X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*dt + S(t[i-1],X[i-1],Y[i-1])*Dx[i-1]
        Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*dt + S(t[i-1],X[i-1],Y[i-1])*Dy[i-1]
                   } 
plot(X,Y,las=1,lwd=3,type="n",xlab=expression(W[t]),ylab=expression(W[t]))
mtext("Brownian motion in 2D plane",line=2,cex=1.2) 
mtext("by normal law",line=0.5,cex=1.2,col="red")
points(x0,y0,type="p",pch=20,col="red2",cex=1.4)
if(Step==FALSE){points(X,Y,type="l",lwd=1,col="black")}
if(Step==TRUE){for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",
                col="black",lwd=1)}}
legend("topleft",border="gray",c("(X0,Y0)"),pch=c(20),col=c("red2"))
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
BMN2_D <- data.frame(t,X,Y)
showData(BMN2_D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(BMN2_D, file = "BMN2D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(BMN2_D)
}

.BMRW3D <- function(N,t0,T,X0,Y0,Z0,Sigma,Output=FALSE)
       {
drifx     <- expression( 0 )
drify     <- expression( 0 )
drifz     <- expression( 0 )
diff      <- expression( Sigma ) 

Ax    <- function(t,x,y,z)  eval(drifx)
Ay    <- function(t,x,y,z)  eval(drify)
Az    <- function(t,x,y,z)  eval(drifz)
S     <- function(t,x,y,z)  eval(diff)
t <- seq (t0 ,T, length =N+1)
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
G <- data.frame(X,Y,Z)
V = 1
a <- c(0,V,0,0)
b <- c(0,0,V,0)
c <- c(0,0,0,V)
labels <- c("Origin", "X", "Y", "Z")
i <- c(1,2,1,3,1,4)
open3d()
segments3d(a[i],b[i],c[i],color = c("black"),lwd= 2.0)
text3d(a,b,c,labels,adj=0.5,col="red",cex=1.2,family=c("serif"))
text3d(X0,Y0,Z0,c("(X0,Y0,Z0)"),adj=c(0.5,-0.25),cex=1.2,family=c("serif"))
points3d(G[1,],color = c("blue"),size=6)
title3d(family=c("serif"),main="Simulation Three-Dimensional for Brownian Motion by a Random Walk",color = c("black"),cex=1.2)
for (i in 1:N) {lines3d(c(G[i,1],G[i+1,1]),c(G[i,2],G[i+1,2]),c(G[i,3],G[i+1,3]),col="red",from ="lines",lwd=2)}
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)

BMRW3_D <- data.frame(t,X,Y,Z)
showData(BMRW3_D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(BMRW3_D, file = "BMRW3D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(BMRW3_D)
}

.BMN3D <- function(N,t0,T,X0,Y0,Z0,Sigma,Output=FALSE)
       {
drifx     <- expression( 0 )
drify     <- expression( 0 )
drifz     <- expression( 0 )
diff      <- expression( Sigma ) 

Ax    <- function(t,x,y,z)  eval(drifx)
Ay    <- function(t,x,y,z)  eval(drify)
Az    <- function(t,x,y,z)  eval(drifz)
S     <- function(t,x,y,z)  eval(diff)
t <- seq (t0 ,T, length =N+1)
Dt= (T-t0)/N 
wx = cumsum(rnorm(N+1,mean=0,sd=sqrt(Dt)))
Dx    <- diff(wx)
wy = cumsum(rnorm(N+1,mean=0,sd=sqrt(Dt)))
Dy    <- diff(wy)
wz = cumsum(rnorm(N+1,mean=0,sd=sqrt(Dt)))
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
G <- data.frame(X,Y,Z)
V = 1
a <- c(0,V,0,0)
b <- c(0,0,V,0)
c <- c(0,0,0,V)
labels <- c("Origin", "X", "Y", "Z")
i <- c(1,2,1,3,1,4)
open3d()
segments3d(a[i],b[i],c[i],color = c("black"),lwd= 2.0)
text3d(a,b,c,labels,adj=0.5,col="red",cex=1.2,family=c("serif"))
text3d(X0,Y0,Z0,c("(X0,Y0,Z0)"),adj=c(0.5,-0.25),cex=1.2,family=c("serif"))
points3d(G[1,],color = c("blue"),size=6)
title3d(family=c("serif"),main="Simulation Three-Dimensional for Brownian Motion by Normal law",color = c("black"),cex=1.2)
for (i in 1:N) {lines3d(c(G[i,1],G[i+1,1]),c(G[i,2],G[i+1,2]),c(G[i,3],G[i+1,3]),col="red",from ="lines",lwd=2)}
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
BMN3_D <- data.frame(t,X,Y,Z)
showData(BMN3_D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(BMN3_D, file = "BMN3D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(BMN3_D)
}

.STS3D <-
function(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty) || !is.expression(driftz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y,Z)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy) || !is.expression(diffz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y,Z)" ),icon="error"))

DDAx   <- D(driftx,"x")
DDDAx  <- D(D(driftx,"x"),"x")
DDSx   <- D(diffx,"x")
DDDSx  <- D(D(diffx,"x"),"x")
Ax     <- function(t,x,y,z)  eval(driftx)
DAxx   <- function(t,x,y,z)  eval(DDAx)
DAxxx  <- function(t,x,y,z)  eval(DDDAx)
Sx     <- function(t,x,y,z)  eval(diffx)
DSx    <- function(t,x,y,z)  eval(DDSx)
DSxx   <- function(t,x,y,z)  eval(DDDSx)

DDAy   <- D(drifty,"y")
DDDAy  <- D(D(drifty,"y"),"y")
DDSy   <- D(diffy,"y")
DDDSy  <- D(D(diffy,"y"),"y")
Ay     <- function(t,x,y,z)  eval(drifty)
DAyy   <- function(t,x,y,z)  eval(DDAy)
DAyyy  <- function(t,x,y,z)  eval(DDDAy)
Sy     <- function(t,x,y,z)  eval(diffy)
DSy    <- function(t,x,y,z)  eval(DDSy)
DSyy   <- function(t,x,y,z)  eval(DDDSy)

DDAz   <- D(driftz,"z")
DDDAz  <- D(D(driftz,"z"),"z")
DDSz   <- D(diffz,"z")
DDDSz  <- D(D(diffz,"z"),"z")
Az     <- function(t,x,y,z)  eval(driftz)
DAzz   <- function(t,x,y,z)  eval(DDAz)
DAzzz  <- function(t,x,y,z)  eval(DDDAz)
Sz     <- function(t,x,y,z)  eval(diffz)
DSz    <- function(t,x,y,z)  eval(DDSz)
DSzz   <- function(t,x,y,z)  eval(DDDSz)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}
Dt = (T-t0)/N
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
wz = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dz    <- diff(wz)
DZx= rnorm(N,0,sqrt((1/3)*Dt^3))
DZy= rnorm(N,0,sqrt((1/3)*Dt^3))
DZz= rnorm(N,0,sqrt((1/3)*Dt^3))
X    <- numeric()
Y    <- numeric()
Z    <- numeric()
X[1] <- x0
Y[1] <- y0
Z[1] <- z0
for (i in 2:(N+1)){
X[i]=X[i-1]+Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dx[i-1]+
     0.5*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSx(t[i-1],X[i-1],Y[i-1],Z[i-1])*((Dx[i-1]^2)-Dt)+
     DAxx(t[i-1],X[i-1],Y[i-1],Z[i-1])*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*DZx[i-1]+0.5*(Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*DAxx(t[i-1],X[i-1],Y[i-1],Z[i-1])+
     0.5*(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*DAxxx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*(Dt^2)+(Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSx(t[i-1],X[i-1],Y[i-1],Z[i-1])+
     0.5*(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*DSxx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*(Dx[i-1]*Dt-DZx[i-1])+
     0.5*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSxx(t[i-1],X[i-1],Y[i-1],Z[i-1])+
     (DSx(t[i-1],X[i-1],Y[i-1],Z[i-1])^2))*((1/3)*(Dx[i-1]^2)-Dt)*Dx[i-1]
Y[i]=Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dy[i-1]+
     0.5*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSy(t[i-1],X[i-1],Y[i-1],Z[i-1])*((Dy[i-1]^2)-Dt)+
     DAyy(t[i-1],X[i-1],Y[i-1],Z[i-1])*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*DZy[i-1]+0.5*(Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*DAyy(t[i-1],X[i-1],Y[i-1],Z[i-1])+
     0.5*(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*DAyyy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*(Dt^2)+(Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSy(t[i-1],X[i-1],Y[i-1],Z[i-1])+
     0.5*(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*DSyy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*(Dy[i-1]*Dt-DZy[i-1])+
     0.5*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSyy(t[i-1],X[i-1],Y[i-1],Z[i-1])+
     (DSx(t[i-1],X[i-1],Y[i-1],Z[i-1])^2))*((1/3)*(Dy[i-1]^2)-Dt)*Dy[i-1]
Z[i]=Z[i-1]+Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dz[i-1]+
     0.5*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSz(t[i-1],X[i-1],Y[i-1],Z[i-1])*((Dz[i-1]^2)-Dt)+
     DAzz(t[i-1],X[i-1],Y[i-1],Z[i-1])*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*DZz[i-1]+0.5*(Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*DAzz(t[i-1],X[i-1],Y[i-1],Z[i-1])+
     0.5*(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*DAzzz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*(Dt^2)+(Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSz(t[i-1],X[i-1],Y[i-1],Z[i-1])+
     0.5*(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])^2)*DSzz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*(Dz[i-1]*Dt-DZz[i-1])+
     0.5*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSzz(t[i-1],X[i-1],Y[i-1],Z[i-1])+
     (DSz(t[i-1],X[i-1],Y[i-1],Z[i-1])^2))*((1/3)*(Dz[i-1]^2)-Dt)*Dz[i-1]
                  } 
G <- data.frame(X,Y,Z)
if(Step==FALSE){
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
lines3d(G[,1],G[,2],G[,3],col="black",from ="lines",lwd=2)
title3d(family=c("serif"),main="Strong Taylor Scheme Order 1.5 : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}

if(Step==TRUE){
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
title3d(family=c("serif"),main="Strong Taylor Scheme Order 1.5 : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
Diff3D <- data.frame(t,X,Y,Z)
showData(Diff3D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Diff3D, file = "SYSDiff3D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Diff3D)
}

.RK33D <-
function(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty) || !is.expression(driftz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y,Z)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy) || !is.expression(diffz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y,Z)" ),icon="error"))

Ax    <- function(t,x,y,z)  eval(driftx)
Ay    <- function(t,x,y,z)  eval(drifty)
Az    <- function(t,x,y,z)  eval(driftz)
Sx    <- function(t,x,y,z)  eval(diffx)
Sy    <- function(t,x,y,z)  eval(diffy)
Sz    <- function(t,x,y,z)  eval(diffz)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
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
XX   <- numeric()
YY   <- numeric()
ZZ   <- numeric()
XXX  <- numeric()
YYY  <- numeric()
ZZZ  <- numeric()
X[1] <- x0
Y[1] <- y0
Z[1] <- z0
for (i in 2:(N+1)){
              XX[i-1]=X[i-1]+0.5*Dt*Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dx[i-1]
              YY[i-1]=Y[i-1]+0.5*Dt*Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dy[i-1]
              ZZ[i-1]=Z[i-1]+0.5*Dt*Az(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dz[i-1]
              XXX[i-1]=X[i-1]-Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+2*Dt*Ax(t[i-1]+0.5*Dt,XX[i-1],Y[i-1],Z[i-1])+
                     (2*Sx(t[i-1]+0.5*Dt,XX[i-1],Y[i-1],Z[i-1])-Sx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dx[i-1]
              YYY[i-1]=Y[i-1]-Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+2*Dt*Ay(t[i-1]+0.5*Dt,X[i-1],YY[i-1],Z[i-1])+
                     (2*Sy(t[i-1]+0.5*Dt,X[i-1],YY[i-1],Z[i-1])-Sy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dy[i-1]
              ZZZ[i-1]=Z[i-1]-Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+2*Dt*Az(t[i-1]+0.5*Dt,X[i-1],Y[i-1],ZZ[i-1])+
                     (2*Sz(t[i-1]+0.5*Dt,X[i-1],Y[i-1],ZZ[i-1])-Sz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dz[i-1]
              X[i] = X[i-1]+(Dt/6)*(Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Ax(t[i-1]+0.5*Dt,XX[i-1],Y[i-1],Z[i-1])+Ax(t[i-1]+Dt,XXX[i-1],Y[i-1],Z[i-1]))+
                     (1/6)*(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Sx(t[i-1]+0.5*Dt,XX[i-1],Y[i-1],Z[i-1])+Sx(t[i-1]+Dt,XXX[i-1],Y[i-1],Z[i-1]))*Dx[i-1]
              Y[i] = Y[i-1]+(Dt/6)*(Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Ay(t[i-1]+0.5*Dt,X[i-1],YY[i-1],Z[i-1])+Ay(t[i-1]+Dt,X[i-1],YYY[i-1],Z[i-1]))+
                     (1/6)*(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Sy(t[i-1]+0.5*Dt,X[i-1],YY[i-1],Z[i-1])+Sy(t[i-1]+Dt,X[i-1],YYY[i-1],Z[i-1]))*Dy[i-1]
              Z[i] = Z[i-1]+(Dt/6)*(Az(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Az(t[i-1]+0.5*Dt,X[i-1],Y[i-1],ZZ[i-1])+Az(t[i-1]+Dt,X[i-1],Y[i-1],ZZZ[i-1]))+
                     (1/6)*(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])+4*Sz(t[i-1]+0.5*Dt,X[i-1],Y[i-1],ZZ[i-1])+Sz(t[i-1]+Dt,X[i-1],Y[i-1],ZZZ[i-1]))*Dz[i-1]
                  } 
G <- data.frame(X,Y,Z)
if(Step==FALSE){
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
lines3d(G[,1],G[,2],G[,3],col="black",from ="lines",lwd=2)
title3d(family=c("serif"),main="Runge-Kutta scheme Order3 : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}

if(Step==TRUE){
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
title3d(family=c("serif"),main="Runge-Kutta scheme Order3 : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
Diff3D <- data.frame(t,X,Y,Z)
showData(Diff3D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Diff3D, file = "SYSDiff3D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Diff3D)
}

.Heun3D <-
function(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty) || !is.expression(driftz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y,Z)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy) || !is.expression(diffz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y,Z)" ),icon="error"))

Ax    <- function(t,x,y,z)  eval(driftx)
Ay    <- function(t,x,y,z)  eval(drifty)
Az    <- function(t,x,y,z)  eval(driftz)
Sx    <- function(t,x,y,z)  eval(diffx)
Sy    <- function(t,x,y,z)  eval(diffy)
Sz    <- function(t,x,y,z)  eval(diffz)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
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
XX   <- numeric()
YY   <- numeric()
ZZ   <- numeric()
X[1] <- x0
Y[1] <- y0
Z[1] <- z0
for (i in 2:(N+1)){
         XX[i-1]= X[i-1]+Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dx[i-1]
         YY[i-1]= Y[i-1]+Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dy[i-1]
         ZZ[i-1]= Z[i-1]+Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt+Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dz[i-1]
         X[i]   = X[i-1]+0.5*Dt*(Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])+Ax(t[i-1],XX[i-1],Y[i-1],Z[i-1]))+
                 0.5*(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sx(t[i-1],XX[i-1],Y[i-1],Z[i-1]))*Dx[i-1]
         Y[i]   = Y[i-1]+0.5*Dt*(Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])+Ay(t[i-1],X[i-1],YY[i-1],Z[i-1]))+
                 0.5*(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sy(t[i-1],X[i-1],YY[i-1],Z[i-1]))*Dy[i-1]
         Z[i]   = Z[i-1]+0.5*Dt*(Az(t[i-1],X[i-1],Y[i-1],Z[i-1])+Az(t[i-1],X[i-1],Y[i-1],ZZ[i-1]))+
                 0.5*(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])+Sz(t[i-1],X[i-1],Y[i-1],ZZ[i-1]))*Dz[i-1]
                  } 
G <- data.frame(X,Y,Z)
if(Step==FALSE){
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
lines3d(G[,1],G[,2],G[,3],col="black",from ="lines",lwd=2)
title3d(family=c("serif"),main="Heun scheme : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}

if(Step==TRUE){
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
title3d(family=c("serif"),main="Heun scheme : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
Diff3D <- data.frame(t,X,Y,Z)
showData(Diff3D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Diff3D, file = "SYSDiff3D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Diff3D)
}

.MilsteinS3D  <-
function(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty) || !is.expression(driftz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y,Z)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy) || !is.expression(diffz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y,Z)" ),icon="error"))

DDAx   <- D(driftx,"x")
DDDAx  <- D(D(driftx,"x"),"x")
DDSx   <- D(diffx,"x")
DDDSx  <- D(D(diffx,"x"),"x")
Ax     <- function(t,x,y,z)  eval(driftx)
DAxx   <- function(t,x,y,z)  eval(DDAx)
DAxxx  <- function(t,x,y,z)  eval(DDDAx)
Sx     <- function(t,x,y,z)  eval(diffx)
DSx    <- function(t,x,y,z)  eval(DDSx)
DSxx   <- function(t,x,y,z)  eval(DDDSx)

DDAy   <- D(drifty,"y")
DDDAy  <- D(D(drifty,"y"),"y")
DDSy   <- D(diffy,"y")
DDDSy  <- D(D(diffy,"y"),"y")
Ay     <- function(t,x,y,z)  eval(drifty)
DAyy   <- function(t,x,y,z)  eval(DDAy)
DAyyy  <- function(t,x,y,z)  eval(DDDAy)
Sy     <- function(t,x,y,z)  eval(diffy)
DSy    <- function(t,x,y,z)  eval(DDSy)
DSyy   <- function(t,x,y,z)  eval(DDDSy)

DDAz   <- D(driftz,"z")
DDDAz  <- D(D(driftz,"z"),"z")
DDSz   <- D(diffz,"z")
DDDSz  <- D(D(diffz,"z"),"z")
Az     <- function(t,x,y,z)  eval(driftz)
DAzz   <- function(t,x,y,z)  eval(DDAz)
DAzzz  <- function(t,x,y,z)  eval(DDDAz)
Sz     <- function(t,x,y,z)  eval(diffz)
DSz    <- function(t,x,y,z)  eval(DDSz)
DSzz   <- function(t,x,y,z)  eval(DDDSz)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
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
X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dx[i-1] +
       0.5 *Sx(t[i-1],X[i-1],Y[i-1],Z[i-1]) * DSx(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Dx[i-1]^2-Dt)+ 
       Dt^(3/2)*(0.5 *Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSx(t[i-1],X[i-1],Y[i-1],Z[i-1]) +
       0.5 *DAxx(t[i-1],X[i-1],Y[i-1],Z[i-1])*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])+
       0.25 *(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])^2) * DSxx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dx[i -1]+ 
       (Dt^2) * (0.5*Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*DAxx(t[i-1],X[i-1],Y[i-1],Z[i-1])+
       0.25 *DAxxx(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])^2))

Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dy[i-1] +
       0.5 *Sy(t[i-1],X[i-1],Y[i-1],Z[i-1]) * DSy(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Dy[i-1]^2-Dt)+ 
       Dt^(3/2)*(0.5 *Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSy(t[i-1],X[i-1],Y[i-1],Z[i-1]) +
       0.5 *DAyy(t[i-1],X[i-1],Y[i-1],Z[i-1])*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])+
       0.25 *(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])^2) * DSyy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dy[i -1]+ 
       (Dt^2) * (0.5*Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*DAyy(t[i-1],X[i-1],Y[i-1],Z[i-1])+
       0.25 *DAyyy(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])^2))

Z[i] = Z[i-1] + Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dz[i-1] +
       0.5 *Sz(t[i-1],X[i-1],Y[i-1],Z[i-1]) * DSz(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Dz[i-1]^2-Dt)+ 
       Dt^(3/2)*(0.5 *Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSz(t[i-1],X[i-1],Y[i-1],Z[i-1]) +
       0.5 *DAzz(t[i-1],X[i-1],Y[i-1],Z[i-1])*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])+
       0.25 *(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])^2) * DSzz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dz[i -1]+ 
       (Dt^2) * (0.5*Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*DAzz(t[i-1],X[i-1],Y[i-1],Z[i-1])+
       0.25 *DAzzz(t[i-1],X[i-1],Y[i-1],Z[i-1])*(Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])^2))
                  } 
G <- data.frame(X,Y,Z)
if(Step==FALSE){
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
lines3d(G[,1],G[,2],G[,3],col="black",from ="lines",lwd=2)
title3d(family=c("serif"),main="Second Milstein scheme : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria") ,side = 1, line = 4, adj = 0.5, cex = .66)
}

if(Step==TRUE){
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
title3d(family=c("serif"),main="Second Milstein scheme : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
Diff3D <- data.frame(t,X,Y,Z)
showData(Diff3D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Diff3D, file = "SYSDiff3D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Diff3D)
}

.Milstein3D  <-
function(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty) || !is.expression(driftz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y,Z)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy) || !is.expression(diffz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y,Z)" ),icon="error"))

DSxx   <- D(diffx,"x")
DSyy   <- D(diffy,"y")
DSzz   <- D(diffz,"z")			
			
Ax    <- function(t,x,y,z)  eval(driftx)
Ay    <- function(t,x,y,z)  eval(drifty)
Az    <- function(t,x,y,z)  eval(driftz)
Sx    <- function(t,x,y,z)  eval(diffx)
Sy    <- function(t,x,y,z)  eval(diffy)
Sz    <- function(t,x,y,z)  eval(diffz)
DSx   <- function(t,x,y,z)  eval(DSxx)
DSy   <- function(t,x,y,z)  eval(DSyy)
DSz   <- function(t,x,y,z)  eval(DSzz)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
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
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dx[i-1]+
           0.5 *Sx(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSx(t[i-1],X[i-1],Y[i-1],Z[i-1])*((Dx[i-1])^2 -Dt)
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dy[i-1]+
           0.5 *Sy(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSy(t[i-1],X[i-1],Y[i-1],Z[i-1])*((Dy[i-1])^2 -Dt) 
    Z[i] = Z[i-1] + Az(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dt + Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*Dz[i-1]+
           0.5 *Sz(t[i-1],X[i-1],Y[i-1],Z[i-1])*DSz(t[i-1],X[i-1],Y[i-1],Z[i-1])*((Dz[i-1])^2 -Dt) 
                  }
G <- data.frame(X,Y,Z)
if(Step==FALSE){
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
lines3d(G[,1],G[,2],G[,3],col="black",from ="lines",lwd=2)
title3d(family=c("serif"),main="Milstein scheme : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}

if(Step==TRUE){
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
title3d(family=c("serif"),main="Milstein scheme : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
Diff3D <- data.frame(t,X,Y,Z)
showData(Diff3D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Diff3D, file = "SYSDiff3D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Diff3D)
}

.Euler3D <-
function(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty) || !is.expression(driftz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y,Z)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy) || !is.expression(diffz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y,Z)" ),icon="error"))

Ax    <- function(t,x,y,z)  eval(driftx)
Ay    <- function(t,x,y,z)  eval(drifty)
Az    <- function(t,x,y,z)  eval(driftz)
Sx    <- function(t,x,y,z)  eval(diffx)
Sy    <- function(t,x,y,z)  eval(diffy)
Sz    <- function(t,x,y,z)  eval(diffz)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
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
if(Step==FALSE){
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
lines3d(G[,1],G[,2],G[,3],col="black",from ="lines",lwd=2)
title3d(family=c("serif"),main="Euler scheme : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria") ,side = 1, line = 4, adj = 0.5, cex = .66)
}

if(Step==TRUE){
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
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
Diff3D <- data.frame(t,X,Y,Z)
showData(Diff3D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Diff3D, file = "SYSDiff3D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Diff3D)
}

.PredCorr3D <-
function(N,T=1,t0,x0,y0,z0,Dt,alpha=0.5,mu=0.5,driftx,drifty,driftz,diffx,diffy,diffz,Step=FALSE,Output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))
			
if ( alpha > 1 || alpha < 0 )
            stop(tkmessageBox(title="Error",message=paste( "0 <= alpha <= 1" ),icon="error"))

if ( mu > 1 || mu < 0 )
            stop(tkmessageBox(title="Error",message=paste( "0 <= mu <= 1" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty) || !is.expression(driftz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y,Z)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy) || !is.expression(diffz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y,Z)" ),icon="error"))

DSx  <- D(diffx,"x")
Ax    <- function(t,x,y,z)  eval(driftx)
Sx    <- function(t,x,y,z)  eval(diffx)
Sxx   <- function(t,x,y,z)  eval(DSx)
SSx   <- function(t,x,y,z)  eval(driftx) - mu * eval(diffx) * eval(DSx)

DSy  <- D(diffy,"y")
Ay    <- function(t,x,y,z)  eval(drifty)
Sy    <- function(t,x,y,z)  eval(diffy)
Syy   <- function(t,x,y,z)  eval(DSy)
SSy   <- function(t,x,y,z)  eval(drifty) - mu * eval(diffy) * eval(DSy)

DSz  <- D(diffz,"z")
Az    <- function(t,x,y,z)  eval(driftz)
Sz    <- function(t,x,y,z)  eval(diffz)
Szz   <- function(t,x,y,z)  eval(DSz)
SSz   <- function(t,x,y,z)  eval(driftz) - mu * eval(diffz) * eval(DSz)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} 
          else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
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
XX   <- numeric()
YY   <- numeric()
ZZ   <- numeric()
X[1] <- XX[1] <- x0
Y[1] <- YY[1] <- y0
Z[1] <- ZZ[1] <- z0
for (i in 2:(N+1)){
XX[i] = XX[i-1] + Ax(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*Dt + Sx(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*Dx[i-1]
YY[i] = YY[i-1] + Ay(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*Dt + Sy(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*Dy[i-1]
ZZ[i] = ZZ[i-1] + Az(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*Dt + Sz(t[i-1],XX[i-1],YY[i-1],ZZ[i-1])*Dz[i-1]
}

for (i in 2:(N+1)){
X[i] = X[i-1] +(alpha*SSx(t[i],XX[i],YY[i],ZZ[i])+(1-alpha)*SSx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dt+
       (mu*Sx(t[i],XX[i],YY[i],ZZ[i])+(1-mu)*Sx(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dx[i-1]
Y[i] = Y[i-1] +(alpha*SSy(t[i],XX[i],YY[i],ZZ[i])+(1-alpha)*SSy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dt+
       (mu*Sy(t[i],XX[i],YY[i],ZZ[i])+(1-mu)*Sy(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dy[i-1]
Z[i] = Z[i-1] +(alpha*SSz(t[i],XX[i],YY[i],ZZ[i])+(1-alpha)*SSz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dt+
       (mu*Sz(t[i],XX[i],YY[i],ZZ[i])+(1-mu)*Sz(t[i-1],X[i-1],Y[i-1],Z[i-1]))*Dz[i-1]
} 
G <- data.frame(X,Y,Z)
if(Step==FALSE){
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
lines3d(G[,1],G[,2],G[,3],col="black",from ="lines",lwd=2)
title3d(family=c("serif"),main="Predictor-Corrector Method : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria") ,side = 1, line = 4, adj = 0.5, cex = .66)
}

if(Step==TRUE){
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
title3d(family=c("serif"),main="Predictor-Corrector Method : Simulation SDE Three-Dimensional",color = c("black"),cex=1.2)
title3d(family=c("serif"),font=4,sub='Copyright 2012, USTHB. Algeria',color = c("blue"),cex=0.8)
par( mar =c(3 ,3 ,3 ,1))
par(mfrow=c(3,1))
plot(t,X,type="l",xlab=expression(time),ylab=expression(X[t]^1),las=1,col="red")
mtext(bquote(dX[t]^1== a[1](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[1](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^1),cex=0.8,adj=0,line=0.1,col="red")
plot(t,Y,type="l",xlab=expression(time),ylab=expression(X[t]^2),las=1,col="blue")
mtext(bquote(dX[t]^2== a[2](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[2](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^2),cex=0.8,adj=0,line=0.1,col="blue")
plot(t,Z,type="l",xlab=expression(time),ylab=expression(X[t]^3),las=1,col="green4")
mtext(bquote(dX[t]^3== a[3](t,X[t]^1,X[t]^2,X[t]^3)*dt + sigma[3](t,X[t]^1,X[t]^2,X[t]^3) *d*W[t]^3),cex=0.8,adj=0,line=0.1,col="green4")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
Diff3D <- data.frame(t,X,Y,Z)
showData(Diff3D, placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Diff3D, file = "SYSDiff3D.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                 }
attach(Diff3D)
}

.snssde3D <-
function(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step=FALSE,Output=FALSE,Methods=c("SchEuler","SchMilstein",
                   "SchMilsteinS","SchTaylor","SchHeun",
                   "SchRK3"),...)
        {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if ( Dt <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Dt > 0" ),icon="error"))

if(!is.expression(driftx) || !is.expression(drifty) || !is.expression(driftz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X,Y,Z)" ),icon="error"))

if(!is.expression(diffx) || !is.expression(diffy) || !is.expression(diffz))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X,Y,Z)" ),icon="error"))

Methods <- match.arg(Methods)

if ( Methods=="SchEuler" )    {R <- .Euler3D(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step,Output)}
if ( Methods=="SchMilstein")  {R <- .Milstein3D(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step,Output)}
if ( Methods=="SchMilsteinS") {R <- .MilsteinS3D(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step,Output)}
if ( Methods=="SchTaylor")    {R <- .STS3D(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step,Output)}
if ( Methods=="SchHeun")      {R <- .Heun3D(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step,Output)}
if ( Methods=="SchRK3")       {R <- .RK33D(N,T=1,t0,x0,y0,z0,Dt,driftx,drifty,driftz,diffx,diffy,diffz,Step,Output)}
      }
	  
.dconShoji <- function(x, t, x0, t0, drift, diff, Output=FALSE)
           {
if(!is.expression(drift))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(diff))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

if( length(x) <= 1)
            stop(tkmessageBox(title="Error",message=paste( "x vector of quantiles, example x<-seq(-3,3,by=0.1)." ),icon="error"))

if(t[1] < 0)
            stop(tkmessageBox(title="Error",message=paste( "time >= 0, value or vector of the times." ),icon="error"))

DAt  <- D(drift,"t")
DAx  <- D(drift,"x")
DAxx <- D(D(drift,"x"),"x")
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
At   <- function(t,x)  eval(DAt)
for(i in 1:length(t)){
plot(x,dnorm(x, mean=(x0 + A(t[i],x)*(exp(Ax(t[i],x)*t[i])-1)/Ax(t[i],x) + (S(t[i],x)^2 * Axx(t[i],x)/2 + At(t[i],x))*(exp(Ax(t[i],x)*t[i]) -1 -Ax(t[i],x)*t[i])/Ax(t[i],x)^2) ,sd=sqrt(S(t[i],x)^2*(exp(2*Ax(t[i],x)*t[i])-1)/(2*Ax(t[i],x)))),
     type="l",xlab="x",ylab=expression(bold(f(list(t,y)/x))),las=1)
mtext(bquote("Evolution Conditional Density at time":.(round(t[i],2))),line=2.5,cex=1.2,adj=0.5)
mtext("Shoji-Ozaki method",line=1,cex=1.2,adj=0.5)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
                     }
f_x <- dnorm(x, mean=(x0 + A(t[length(t)],x)*(exp(Ax(t[length(t)],x)*t[length(t)])-1)/Ax(t[length(t)],x) + (S(t[length(t)],x)^2 * Axx(t[length(t)],x)/2 + At(t[length(t)],x))*(exp(Ax(t[length(t)],x)*t[length(t)]) -1 -Ax(t[length(t)],x)*t[length(t)])/Ax(t[length(t)],x)^2) ,sd=sqrt(S(t[length(t)],x)^2*(exp(2*Ax(t[length(t)],x)*t[length(t)])-1)/(2*Ax(t[length(t)],x))))
Result <- data.frame(x,f_x)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "CD.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.dconKessler <- function(x, t, x0, t0, drift, diff, Output=FALSE)
          {
if(!is.expression(drift))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(diff))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

if( length(x) <= 1)
            stop(tkmessageBox(title="Error",message=paste( "x vector of quantiles, example x<-seq(-3,3,by=0.1)." ),icon="error"))

if(t[1] < 0)
            stop(tkmessageBox(title="Error",message=paste( "time >= 0, value or vector of the times." ),icon="error"))

DAx  <- D(drift,"x")
DAxx <- D(D(drift,"x"),"x")
DSx  <- D(diff,"x")
DSxx <- D(D(diff,"x"),"x")
A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)
Ax   <- function(t,x)  eval(DAx)
Axx  <- function(t,x)  eval(DAxx)
Sx   <- function(t,x)  eval(DSx)
Sxx  <- function(t,x)  eval(DSxx)

for (i in 1:length(t)){
plot(x,dnorm(x, mean = (x0 + A(t0, x0) * t[i] + (A(t0, x0) * Ax(t0, x0) + 0.5 * (S(t0, x0)^2 * Axx(t0, x0))) * (t[i]^2)/2),
       sd = sqrt((x0^2 + (2 * A(t0, x0) * x0 + Sxx(t0, x0)^2) * t[i] + (2 * A(t0, x0) * (Ax(t0, x0) *
                  x0 + A(t0, x0) + S(t0, x0) * Sx(t0, x0)) + S(t0, x0)^2 * (Axx(t0, x0) * x0 + 2 * Ax(t0, x0) + Sx(t0, x0)^2 +
                  S(t0, x0) * Sxx(t0, x0))) * (t[i]^2)/2 - ((x0 + A(t0, x0) * t[i] + (A(t0, x0) * Ax(t0, x0) + 0.5 * (S(t0, x0)^2 * Axx(t0, x0))) * (t[i]^2)/2))^2))),
                  type="l",xlab="x",ylab=expression(bold(f(list(t,y)/x))),las=1)
mtext(bquote("Evolution Conditional Density at time":.(round(t[i],2))),line=2.5,cex=1.2,adj=0.5)
mtext("Kessler method",line=1,cex=1.2,adj=0.5)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)                     
}
f_x <- dnorm(x, mean = (x0 + A(t0, x0) * t[length(t)] + (A(t0, x0) * Ax(t0, x0) + 0.5 * (S(t0, x0)^2 * Axx(t0, x0))) * (t[length(t)]^2)/2),
                sd = sqrt((x0^2 + (2 * A(t0, x0) * x0 + Sxx(t0, x0)^2) * t[length(t)] + (2 * A(t0, x0) * (Ax(t0, x0) *
                  x0 + A(t0, x0) + S(t0, x0) * Sx(t0, x0)) + S(t0, x0)^2 * (Axx(t0, x0) * x0 + 2 * Ax(t0, x0) + Sx(t0, x0)^2 +
                  S(t0, x0) * Sxx(t0, x0))) * (t[length(t)]^2)/2 - ((x0 + A(t0, x0) * t[length(t)] + (A(t0, x0) * Ax(t0, x0) + 0.5 * (S(t0, x0)^2 * Axx(t0, x0))) * (t[length(t)]^2)/2))^2)))
Result <- data.frame(x,f_x)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "CD.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.dconEuler <- function(x, t, x0,t0, drift, diff, Output=FALSE)
                    {
if(!is.expression(drift))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(diff))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

if( length(x) <= 1)
            stop(tkmessageBox(title="Error",message=paste( "x vector of quantiles, example x<-seq(-3,3,by=0.1)." ),icon="error"))

if(t[1] < 0)
            stop(tkmessageBox(title="Error",message=paste( "time >= 0, value or vector of the times." ),icon="error"))

A    <- function(t,x)  eval(drift)
S    <- function(t,x)  eval(diff)
for(i in 1:length(t)){
plot(x,dnorm (x, mean = x0 - A(t[i],x)*t[i], sd= sqrt(t[i])*S(t[i],x)),type="l",xlab="x",ylab=expression(bold(f(list(t,y)/x))),las=1)
mtext(bquote("Evolution Conditional Density at time":.(round(t[i],2))),line=2.5,cex=1.2,adj=0.5)
mtext("Euler method",line=1,cex=1.2,adj=0.5)
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
                     }
f_x <- dnorm (x, mean = x0 - A(t[length(t)],x)*t[length(t)], sd= sqrt(t[length(t)])*S(t[length(t)],x))
Result <- data.frame(x,f_x)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "CD.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}
	  
.Appdcon <- function(x, t, x0,t0, drift, diff, Output=FALSE, Methods=c("Euler","Shoji-Ozaki","Kessler"),...)
         {

if(!is.expression(drift))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of drift must be expressions f(t,X)" ),icon="error"))

if(!is.expression(diff))
            stop(tkmessageBox(title="Error",message=paste( "The coefficients of diffusion must be expressions f(t,X)" ),icon="error"))

if( length(x) <= 1)
            stop(tkmessageBox(title="Error",message=paste( "x vector of quantiles, example x<-seq(-3,3,by=0.1)." ),icon="error"))

if(t[1] < 0)
            stop(tkmessageBox(title="Error",message=paste( "time >= 0, value or vector of the times." ),icon="error"))

Methods <- match.arg(Methods)

if ( Methods=="Euler" )     {R <- .dconEuler(x, t, x0,t0, drift, diff, Output)}
if ( Methods=="Shoji-Ozaki"){R <- .dconShoji(x, t, x0, t0, drift, diff, Output)}
if ( Methods=="Kessler")    {R <- .dconKessler(x, t, x0, t0, drift, diff, Output)}
}


.Sharosc <- function(N,T,x0,v0,lambda,omega,sigma,Step=FALSE,Output=FALSE)
        {

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if( T <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "T > 0" ),icon="error"))

if( lambda < 0) 
            stop(tkmessageBox(title="Error",message=paste( "lambda >= 0" ),icon="error"))

if( omega <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "omega > 0" ),icon="error"))

if( sigma < 0) 
            stop(tkmessageBox(title="Error",message=paste( "sigma >= 0" ),icon="error"))

diffx  <- expression(0)
diffy  <- expression(sigma)
driftx <- expression(y)
drifty <- expression(-(2*lambda*y + omega^2 * x))
DSxx   <- D(diffx,"x")
DSyy   <- D(diffy,"y")
Ax     <- function(t,x,y)  eval(driftx)
Ay     <- function(t,x,y)  eval(drifty)
Sx     <- function(t,x,y)  eval(diffx)
DSx    <- function(t,x,y)  eval(DSxx)
Sy     <- function(t,x,y)  eval(diffy)
DSy    <- function(t,x,y)  eval(DSyy)
Dt = T/N
t <- seq(0,T,length=N+1)
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X     <- numeric()
Y     <- numeric()
X[1]  <- x0
Y[1]  <- v0
for (i in 2:(N+1)){
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]+
           0.5 *Sx(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1])*((Dx[i-1])^2 -Dt)
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]+
           0.5 *Sy(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1])*((Dy[i-1])^2 -Dt) 
                  } 
if(Step==FALSE){
plot(X,(Y/omega),type="l",axes = FALSE,xlab=expression(x[t]*(mm)),ylab=expression(V/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
points(x0,v0/omega,pch=20,col="red")
mtext(expression("The phase portrait of stochastic harmonic oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second+2*lambda*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(lambda==.(lambda)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="l",axes = FALSE,ylab=expression(x[t]*(mm)),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
points(0,x0,pch=20,col="red")
mtext(expression("Temporal evolution of stochastic harmonic oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second+2*lambda*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(lambda==.(lambda)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
V <- Y/omega
plot(X,(Y/omega),type="n",axes = FALSE,xlab=expression(x[t]*(mm)),ylab=expression(V/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
mtext(expression("The phase portrait of stochastic harmonic oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second+2*lambda*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(lambda==.(lambda)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
points(x0,v0/omega,pch=20,col="red")
for (i in 1:N){lines(c(X[i],X[i+1]),c(V[i],V[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="n",axes = FALSE,ylab=expression(x[t]*(mm)),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
mtext(expression("Temporal evolution of stochastic harmonic oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second+2*lambda*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(lambda==.(lambda)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
points(0,x0,pch=20,col="red")
for (i in 1:N){lines(c(t[i],t[i+1]),c(X[i],X[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
V  <- Y
Result <- data.frame(t,X,V/omega)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "SHO.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.Spendu <- function(N,T,theta0,theta1,lambda,omega,sigma,Step=FALSE,Output=FALSE)
        {
if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
if( T <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "T > 0" ),icon="error"))
			
if( theta0 >= pi || theta0 <= -pi ) 
            stop(tkmessageBox(title="Error",message=paste( "-pi < theta0 < pi" ),icon="error"))
if( lambda < 0) 
            stop(tkmessageBox(title="Error",message=paste( "lambda >= 0" ),icon="error"))
if( omega <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "omega > 0" ),icon="error"))
if( sigma < 0) 
            stop(tkmessageBox(title="Error",message=paste( "sigma >= 0" ),icon="error"))
diffx  <- expression(0)
diffy  <- expression(sigma)
driftx <- expression(y)
drifty <- expression(-(2*lambda*y + omega^2 * sin(x)))
DSxx   <- D(diffx,"x")
DSyy   <- D(diffy,"y")
Ax     <- function(t,x,y)  eval(driftx)
Ay     <- function(t,x,y)  eval(drifty)
Sx     <- function(t,x,y)  eval(diffx)
DSx    <- function(t,x,y)  eval(DSxx)
Sy     <- function(t,x,y)  eval(diffy)
DSy    <- function(t,x,y)  eval(DSyy)
Dt = T/N
t <- seq(0,T,length=N+1)
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X     <- numeric()
Y     <- numeric()
X[1]  <- theta0
Y[1]  <- theta1
for (i in 2:(N+1)){
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]+
           0.5 *Sx(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1])*((Dx[i-1])^2 -Dt)
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]+
           0.5 *Sy(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1])*((Dy[i-1])^2 -Dt) 
                  } 
if(Step==FALSE){
plot(X,(Y/omega),type="l",axes = FALSE,xlab=expression(theta*(rad)),ylab=expression(theta*minute/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
points(theta0,theta1/omega,pch=20,col="red")
mtext(expression("The phase portrait of stochastic pendulum"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(theta[t]*second+2*lambda*theta[t]*minute+omega^2*sin(theta[t])==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(theta[0]==.(theta0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(theta[0]*minute==.(theta1)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(lambda==.(lambda)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="l",axes = FALSE,ylab=expression(theta[t]*(rad)),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
points(0,theta0,pch=20,col="red")
mtext(expression("Temporal evolution of stochastic pendulum"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(theta[t]*second+2*lambda*theta[t]*minute+omega^2*sin(theta[t])==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(theta[0]==.(theta0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(theta[0]*minute==.(theta1)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(lambda==.(lambda)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
V = Y/omega
plot(X,(Y/omega),type="n",axes = FALSE,xlab=expression(theta*(rad)),ylab=expression(theta*minute/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
mtext(expression("The phase portrait of stochastic pendulum"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(theta[t]*second+2*lambda*theta[t]*minute+omega^2*sin(theta[t])==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(theta[0]==.(theta0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(theta[0]*minute==.(theta1)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(lambda==.(lambda)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
points(theta0,theta1/omega,pch=20,col="red")
for (i in 1:N){lines(c(X[i],X[i+1]),c(V[i],V[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="n",axes = FALSE,ylab=expression(theta[t]*(rad)),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
mtext(expression("Temporal evolution of stochastic pendulum"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(theta[t]*second+2*lambda*theta[t]*minute+omega^2*sin(theta[t])==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(theta[0]==.(theta0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(theta[0]*minute==.(theta1)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(lambda==.(lambda)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
points(0,theta0,pch=20,col="red")
for (i in 1:N){lines(c(t[i],t[i+1]),c(X[i],X[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
theta1  <- Y
theta   <- X
Result <- data.frame(t,theta,theta1/omega)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.Svandp <- function(N,T,x0,v0,a,b,omega,sigma,Step=FALSE,Output=FALSE)
        {
if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
if( T <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "T > 0" ),icon="error"))
if( a < 0) 
            stop(tkmessageBox(title="Error",message=paste( "a >= 0" ),icon="error"))
if( b <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "b > 0" ),icon="error"))
if( omega < 0) 
            stop(tkmessageBox(title="Error",message=paste( "omega >= 0" ),icon="error"))
if( sigma < 0) 
            stop(tkmessageBox(title="Error",message=paste( "sigma >= 0" ),icon="error"))			
diffx  <- expression(0)
diffy  <- expression(sigma)
driftx <- expression(y)
drifty <- expression(-(a*y*(b^-2 * x^2 -1)+omega^2 * x))
DSxx   <- D(diffx,"x")
DSyy   <- D(diffy,"y")
Ax     <- function(t,x,y)  eval(driftx)
Ay     <- function(t,x,y)  eval(drifty)
Sx     <- function(t,x,y)  eval(diffx)
DSx    <- function(t,x,y)  eval(DSxx)
Sy     <- function(t,x,y)  eval(diffy)
DSy    <- function(t,x,y)  eval(DSyy)
Dt = T/N
t <- seq(0,T,length=N+1)
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X     <- numeric()
Y     <- numeric()
X[1]  <- x0
Y[1]  <- v0
for (i in 2:(N+1)){
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]+
           0.5 *Sx(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1])*((Dx[i-1])^2 -Dt)
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]+
           0.5 *Sy(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1])*((Dy[i-1])^2 -Dt) 
                  } 
if(Step==FALSE){
plot(X,(Y/omega),type="l",axes = FALSE,xlab=expression(x[t]),ylab=expression(x[t]*minute/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
points(x0,v0/omega,pch=20,col="red")
mtext(expression("The phase portrait of stochastic Van Der Pol oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second+a*x[t]*minute*(x[t]^2 / b - 1)+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=2.9,cex=1,adj=1,col="blue")
mtext(bquote(b==.(b)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="l",axes = FALSE,ylab=expression(x[t]),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
points(0,x0,pch=20,col="red")
mtext(expression("Temporal evolution of stochastic Van Der Pol oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second+a*x[t]*minute*(x[t]^2 / b - 1)+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=2.9,cex=1,adj=1,col="blue")
mtext(bquote(b==.(b)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
V <- Y/omega
plot(X,(Y/omega),type="n",axes = FALSE,xlab=expression(x[t]),ylab=expression(x[t]*minute/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
mtext(expression("The phase portrait of stochastic Van Der Pol oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second+a*x[t]*minute*(x[t]^2 / b - 1)+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=2.9,cex=1,adj=1,col="blue")
mtext(bquote(b==.(b)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
points(x0,v0/omega,pch=20,col="red")
for (i in 1:N){lines(c(X[i],X[i+1]),c(V[i],V[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="n",axes = FALSE,ylab=expression(x[t]),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
mtext(expression("Temporal evolution of stochastic Van Der Pol oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second+a*x[t]*minute*(x[t]^2 / b - 1)+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=2.9,cex=1,adj=1,col="blue")
mtext(bquote(b==.(b)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
points(0,x0,pch=20,col="red")
for (i in 1:N){lines(c(t[i],t[i+1]),c(X[i],X[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
X1  <- Y
Result <- data.frame(t,X,X1/omega)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "SVDPO.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.Srayle <- function(N,T,x0,v0,a,omega,sigma,Step=FALSE,Output=FALSE)
        {
if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
if( T <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "T > 0" ),icon="error"))
if( a <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "a > 0" ),icon="error"))
if( omega <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "omega > 0" ),icon="error"))
if( sigma < 0) 
            stop(tkmessageBox(title="Error",message=paste( "sigma >= 0" ),icon="error"))
diffx  <- expression(0)
diffy  <- expression(sigma)
driftx <- expression(y)
drifty <- expression( a*(1-y^2)*y - x*omega^2 )
DSxx   <- D(diffx,"x")
DSyy   <- D(diffy,"y")
Ax     <- function(t,x,y)  eval(driftx)
Ay     <- function(t,x,y)  eval(drifty)
Sx     <- function(t,x,y)  eval(diffx)
DSx    <- function(t,x,y)  eval(DSxx)
Sy     <- function(t,x,y)  eval(diffy)
DSy    <- function(t,x,y)  eval(DSyy)
Dt = T/N
t <- seq(0,T,length=N+1)
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X     <- numeric()
Y     <- numeric()
X[1]  <- x0
Y[1]  <- v0
for (i in 2:(N+1)){
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]+
           0.5 *Sx(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1])*((Dx[i-1])^2 -Dt)
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]+
           0.5 *Sy(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1])*((Dy[i-1])^2 -Dt) 
                  } 
if(Step==FALSE){
plot(X,(Y/omega),type="l",axes = FALSE,xlab=expression(X[t]),ylab=expression(V/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
points(x0,v0/omega,pch=20,col="red")
mtext(expression("The phase portrait of stochastic rayleigh oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second-a*(1-x[t]*minute^2)*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="l",axes = FALSE,ylab=expression(x[t]),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
points(0,x0,pch=20,col="red")
mtext(expression("Temporal evolution of rayleigh oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second-a*(1-x[t]*minute^2)*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
V = Y/omega
plot(X,(Y/omega),type="n",axes = FALSE,xlab=expression(x[t]),ylab=expression(V/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
mtext(expression("The phase portrait of stochastic rayleigh oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second-a*(1-x[t]*minute^2)*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
points(x0,v0/omega,pch=20,col="red")
for (i in 1:N){lines(c(X[i],X[i+1]),c(V[i],V[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="n",axes = FALSE,ylab=expression(x[t]),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
mtext(expression("Temporal evolution of stochastic rayleigh oscillator"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second-a*(1-x[t]*minute^2)*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
points(0,x0,pch=20,col="red")
for (i in 1:N){lines(c(t[i],t[i+1]),c(X[i],X[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
v0  <- Y
theta   <- X
Result <- data.frame(t,theta,v0/omega)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "SRO.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.SSCPP <- function(N,T,theta0,theta1,a,b,omega,sigma,K0=1,Prd=6,Step=FALSE,Output=FALSE)
        {
if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
if( Prd <= 0 )   
            stop(tkmessageBox(title="Error",message=paste( " Period > 0" ),icon="error"))
if( K0 <= 0 )   
            stop(tkmessageBox(title="Error",message=paste( " K0 > 0" ),icon="error"))
if( T <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "T > 0" ),icon="error"))	
if( theta0 >= pi || theta0 <= -pi ) 
            stop(tkmessageBox(title="Error",message=paste( "-pi < theta0 < pi" ),icon="error"))
if( a < 0) 
            stop(tkmessageBox(title="Error",message=paste( "a >= 0" ),icon="error"))
if( b < 0) 
            stop(tkmessageBox(title="Error",message=paste( "b >= 0" ),icon="error"))			
if( omega < 0) 
            stop(tkmessageBox(title="Error",message=paste( "omega >= 0" ),icon="error"))
if( sigma < 0) 
            stop(tkmessageBox(title="Error",message=paste( "sigma >= 0" ),icon="error"))
diffx  <- expression(0)
diffy  <- expression(sigma)
driftx <- expression(y)
drifty <- expression(-(a*y+b+omega*sin(x)))
DSxx   <- D(diffx,"x")
DSyy   <- D(diffy,"y")
Ax     <- function(t,x,y)  eval(driftx)
Ay     <- function(t,x,y)  eval(drifty)
Sx     <- function(t,x,y)  eval(diffx)
DSx    <- function(t,x,y)  eval(DSxx)
Sy     <- function(t,x,y)  eval(diffy)
DSy    <- function(t,x,y)  eval(DSyy)
Dt = T/N
t <- seq(0,T,length=N+1)
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X     <- numeric()
Y     <- numeric()
X[1]  <- theta0
Y[1]  <- theta1
for (i in 2:(N+1)){
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]+
           0.5 *Sx(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1])*((Dx[i-1])^2 -Dt)
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]+
           0.5 *Sy(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1])*((Dy[i-1])^2 -Dt) 
                  } 
if(Step==FALSE){
plot(X,(Y/omega),type="l",axes = FALSE,xlab=expression(theta*(rad)),ylab=expression(theta*minute/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
points(theta0,theta1/omega,pch=20,col="red")
mtext(expression("The phase portrait of stochastic system with a cylindric phase plane"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(theta[t]*second+a*theta[t]*minute+b+omega^2*sin(theta[t])==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(theta[0]==.(theta0)),line=1.9,adj=0.78,cex=1,col="blue")
mtext(bquote(theta[0]*minute==.(theta1)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=0.1,cex=1,adj=0.78,col="blue")
mtext(bquote(b==.(b)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="l",axes = FALSE,ylab=expression(theta[t]*(rad)),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
points(0,theta0,pch=20,col="red")
mtext(expression("Temporal evolution of stochastic system with a cylindric phase plane"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(theta[t]*second+a*theta[t]*minute+b+omega^2*sin(theta[t])==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(theta[0]==.(theta0)),line=1.9,adj=0.78,cex=1,col="blue")
mtext(bquote(theta[0]*minute==.(theta1)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=0.1,cex=1,adj=0.78,col="blue")
mtext(bquote(b==.(b)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
x <- seq(-Prd*pi, Prd*pi, length.out = 100)
y <- seq(-Prd*pi*0.7, Prd*pi*0.7, length.out = 100)
rotsinc <- function(x,y){
     sinc <- function(x) { y <- x ; y[is.na(y)] <- 1; y }
     exp( -(a*(y^2 + 2*b*x - 2*omega^2*cos(x)))/(2*pi*K0) ) 
 }
sinc.exp <- expression(italic(p[s]) == exp*bgroup("(",frac(-a*(y^2 + 2*b*x -2*omega^2*cos(x)),2*pi*K[0]),")"))
z <- outer(x, y, rotsinc)
oldpar <- par(bg = "white")
persp(x, y, z, theta = -40, phi = 20, expand = 0.5, col = "lightblue",
       ltheta = 180, shade = 0.75, ticktype = "detailed",
       xlab = "X", ylab = "Y", zlab = "Z")
title(sub="The Fokker-Planck equation\nSystem with a Cylindric Phase Plane\n\n\n")
title(main = sinc.exp)
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
V = Y/omega
plot(X,(Y/omega),type="n",axes = FALSE,xlab=expression(theta*(rad)),ylab=expression(theta*minute/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
mtext(expression("The phase portrait of stochastic system with a cylindric phase plane"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(theta[t]*second+a*theta[t]*minute+b+omega^2*sin(theta[t])==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(theta[0]==.(theta0)),line=1.9,adj=0.78,cex=1,col="blue")
mtext(bquote(theta[0]*minute==.(theta1)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=0.1,cex=1,adj=0.78,col="blue")
mtext(bquote(b==.(b)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
points(theta0,theta1/omega,pch=20,col="red")
for (i in 1:N){lines(c(X[i],X[i+1]),c(V[i],V[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="n",axes = FALSE,ylab=expression(theta[t]*(rad)),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
mtext(expression("Temporal evolution of stochastic system with a cylindric phase plane"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(theta[t]*second+a*theta[t]*minute+b+omega^2*sin(theta[t])==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(theta[0]==.(theta0)),line=1.9,adj=0.78,cex=1,col="blue")
mtext(bquote(theta[0]*minute==.(theta1)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=0.1,cex=1,adj=0.78,col="blue")
mtext(bquote(b==.(b)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
points(0,theta0,pch=20,col="red")
for (i in 1:N){lines(c(t[i],t[i+1]),c(X[i],X[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
x <- seq(-Prd*pi, Prd*pi, length.out = 100)
y <- seq(-Prd*pi*0.7, Prd*pi*0.7, length.out = 100)
rotsinc <- function(x,y){
     sinc <- function(x) { y <- x ; y[is.na(y)] <- 1; y }
     exp( -(a*(y^2 + 2*b*x - 2*omega^2*cos(x)))/(2*pi*K0) ) 
 }
sinc.exp <- expression(italic(p[s]) == exp*bgroup("(",frac(-a*(y^2 + 2*b*x -2*omega^2*cos(x)),2*pi*K[0]),")"))
z <- outer(x, y, rotsinc)
oldpar <- par(bg = "white")
persp(x, y, z, theta = -40, phi = 20, expand = 0.5, col = "lightblue",
       ltheta = 180, shade = 0.75, ticktype = "detailed",
       xlab = "X", ylab = "Y", zlab = "Z")
title(sub="The Fokker-Planck equation\nSystem with a Cylindric Phase Plane\n\n\n")
title(main = sinc.exp)
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
theta1  <- Y
theta   <- X
Result <- data.frame(t,theta,theta1/omega)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "SCPP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.Sosadd <- function(N,T,x0,v0,a,omega,sigma,K0=1,Step=FALSE,Output=FALSE)
        {
if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
if( K0 <= 0 )   
            stop(tkmessageBox(title="Error",message=paste( " K0 > 0" ),icon="error"))
if( T <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "T > 0" ),icon="error"))
if( a < 0) 
            stop(tkmessageBox(title="Error",message=paste( "a >= 0" ),icon="error"))
if( omega <= 0) 
            stop(tkmessageBox(title="Error",message=paste( "omega > 0" ),icon="error"))
if( sigma < 0) 
            stop(tkmessageBox(title="Error",message=paste( "sigma >= 0" ),icon="error"))
diffx  <- expression(0)
diffy  <- expression(sigma)
driftx <- expression(y)
drifty <- expression(a*(1- x^2 - y^2)*y-omega^2*x)
DSxx   <- D(diffx,"x")
DSyy   <- D(diffy,"y")
Ax     <- function(t,x,y)  eval(driftx)
Ay     <- function(t,x,y)  eval(drifty)
Sx     <- function(t,x,y)  eval(diffx)
DSx    <- function(t,x,y)  eval(DSxx)
Sy     <- function(t,x,y)  eval(diffy)
DSy    <- function(t,x,y)  eval(DSyy)
Dt = T/N
t <- seq(0,T,length=N+1)
wx = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dx    <- diff(wx)
wy = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
Dy    <- diff(wy)
X     <- numeric()
Y     <- numeric()
X[1]  <- x0
Y[1]  <- v0
for (i in 2:(N+1)){
    X[i] = X[i-1] + Ax(t[i-1],X[i-1],Y[i-1])*Dt + Sx(t[i-1],X[i-1],Y[i-1])*Dx[i-1]+
           0.5 *Sx(t[i-1],X[i-1],Y[i-1])*DSx(t[i-1],X[i-1],Y[i-1])*((Dx[i-1])^2 -Dt)
    Y[i] = Y[i-1] + Ay(t[i-1],X[i-1],Y[i-1])*Dt + Sy(t[i-1],X[i-1],Y[i-1])*Dy[i-1]+
           0.5 *Sy(t[i-1],X[i-1],Y[i-1])*DSy(t[i-1],X[i-1],Y[i-1])*((Dy[i-1])^2 -Dt) 
                  } 
if(Step==FALSE){
plot(X,(Y/omega),type="l",axes = FALSE,xlab=expression(x[t]),ylab=expression(V/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
points(x0,v0/omega,pch=20,col="red")
mtext(expression("The phase portrait of stochastic oscillator with additive noise"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second-a*(1-x[t]^2-x[t]*minute^2)*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="l",axes = FALSE,ylab=expression(x[t]),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
points(0,x0,pch=20,col="red")
mtext(expression("Temporal evolution of stochastic oscillator with additive noise"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second-a*(1-x[t]^2-x[t]*minute^2)*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
x <- seq(-2*pi, 2*pi, length.out = 100)
y <- seq(-2*pi, 2*pi,length.out=100)
rotsinc <- function(x,y){
     sinc <- function(x) { y <- x ; y[is.na(y)] <- 1; y }
     exp( -(a*(x^2 + y^2)^2)/(4*pi*K0) ) 
 }
sinc.exp <- expression(italic(p[s]) == exp*bgroup("(",frac(-a*(x^2+y^2)^2,2*pi*K[0]),")"))
z <- outer(x, y, rotsinc)
oldpar <- par(bg = "white")
persp(x, y, z, theta = -40, phi = 20, expand = 0.5, col = "lightblue",
       ltheta = 180, shade = 0.75, ticktype = "detailed",
       xlab = "X", ylab = "Y", zlab = "Z")
title(sub="The Fokker-Planck equation\nStochastic Oscillator with Additive Noise\n\n\n")
title(main = sinc.exp)
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
if(Step==TRUE){
V <- Y/omega
plot(X,(Y/omega),type="n",axes = FALSE,xlab=expression(x[t]),ylab=expression(V/omega))
box()
axis(1, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
axis(2, at = round(seq(min(Y)/omega,max(Y)/omega,length=10),0), labels = TRUE,las=1)
mtext(expression("The phase portrait of stochastic oscillator with additive noise"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second-a*(1-x[t]^2-x[t]*minute^2)*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
points(x0,v0/omega,pch=20,col="red")
for (i in 1:N){lines(c(X[i],X[i+1]),c(V[i],V[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
plot(seq(0,T,length=N+1),X,type="n",axes = FALSE,ylab=expression(x[t]),xlab="t (s)",las=1)
box()
axis(1,labels=TRUE)
axis(2, at = round(seq(min(X),max(X),length=10),0), labels = TRUE,las=1)
mtext(expression("Temporal evolution of stochastic oscillator with additive noise"),line=2.7,adj=0.5,cex=1,col="black")
mtext(expression(x[t]*second-a*(1-x[t]^2-x[t]*minute^2)*x[t]*minute+omega^2*x[t]==epsilon[t]),line=1.3,adj=0,cex=1,col="red")
mtext(expression(bold(E)(epsilon[t]*epsilon[t+h])==sigma*delta*(h)),line=0.05,adj=0,cex=1,col="red")
mtext(bquote(x[0]==.(x0)),line=1,adj=0.78,cex=1,col="blue")
mtext(bquote(v[0]==.(v0)),line=0.1,adj=0.78,cex=1,col="blue")
mtext(bquote(a==.(a)),line=1.9,cex=1,adj=1,col="blue")
mtext(bquote(omega==.(omega)),line=1.0,cex=1,adj=1,col="blue")
mtext(bquote(sigma==.(sigma)),line=0.1,cex=1,adj=1,col="blue")
points(0,x0,pch=20,col="red")
for (i in 1:N){lines(c(t[i],t[i+1]),c(X[i],X[i+1]),type="l",col="black",lwd=1)}
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
windows()
x <- seq(-2*pi, 2*pi, length.out = 100)
y <- seq(-2*pi, 2*pi,length.out=100)
rotsinc <- function(x,y){
     sinc <- function(x) { y <- x ; y[is.na(y)] <- 1; y }
     exp( -(a*(x^2 + y^2)^2)/(4*pi*K0) ) 
 }
sinc.exp <- expression(italic(p[s]) == exp*bgroup("(",frac(-a*(x^2+y^2)^2,2*pi*K[0]),")"))
z <- outer(x, y, rotsinc)
oldpar <- par(bg = "white")
persp(x, y, z, theta = -40, phi = 20, expand = 0.5, col = "lightblue",
       ltheta = 180, shade = 0.75, ticktype = "detailed",
       xlab = "X", ylab = "Y", zlab = "Z")
title(sub="The Fokker-Planck equation\nStochastic Oscillator with Additive Noise\n\n\n")
title(main = sinc.exp)
mtext(paste("Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
V  <- Y
Result <- data.frame(t,X,V/omega)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "SOAN.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
}

.Stweibull <-
function(N,t0,x0,T,shape,scale,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( shape <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "shape > 0" ),icon="error"))
if( scale <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "scale > 0" ),icon="error"))

temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rweibull(u,shape,scale)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :Weibull(.(shape),.(scale))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stexp <-
function(N,t0,x0,T,rate,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( rate <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "rate > 0" ),icon="error"))
temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rexp(u,rate)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :Exp(.(rate))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stchisq <-
function(N,t0,x0,T,df,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( df <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "df > 0, degrees of freedom (non-negative, but can be non-integer)." ),icon="error"))
temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rchisq(u,df)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :chi[.(df)]^2),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stbeta <-
function(N,t0,x0,T,shape1, shape2,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( shape1 <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "shape1 > 0." ),icon="error"))
if( shape2 <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "shape2 > 0." ),icon="error"))
temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rbeta(u,shape1, shape2)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :B(.(shape1),.(shape2))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stcauchy <-
function(N,t0,x0,T,location, scale,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( scale <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "scale > 0." ),icon="error"))
temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rcauchy(u,location, scale)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :Cauchy(.(location),.(scale))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stlnorm <-
function(N,t0,x0,T,meanlog, sdlog,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( sdlog <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "sdlog > 0." ),icon="error"))

temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rlnorm(u, meanlog, sdlog)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :logN(.(meanlog),.(sdlog))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stlnorm3 <-
function(N,t0,x0,T,meanlog, sdlog,thres,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( sdlog <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "sdlog > 0." ),icon="error"))

temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rlnorm(u,meanlog, sdlog)+thres
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :logN[3](.(meanlog),.(sdlog),.(thres))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stgamma3 <-
function(N,t0,x0,T,shape,rate,thres,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( shape <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "shape > 0." ),icon="error"))
if( rate <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "rate > 0." ),icon="error"))
temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rgamma(u,shape,rate)+thres
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :Gamma[3](.(shape),.(rate),.(thres))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stlgamma3 <-
function(N,t0,x0,T,shape,rate,thres,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( shape <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "shape > 0." ),icon="error"))
if( rate <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "rate > 0." ),icon="error"))
temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- exp(rgamma(u,shape,rate)+thres)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :log*Gamma[3](.(shape),.(rate),.(thres))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stweibull3 <-
function(N,t0,x0,T,shape,scale,thres,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( shape <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "shape > 0" ),icon="error"))
if( scale <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "scale > 0" ),icon="error"))

temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- thres + rweibull(u,shape,scale)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :Weibull[3](.(shape),.(scale),.(thres))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stlogis <-
function(N,t0,x0,T,location, scale,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( scale <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "scale > 0." ),icon="error"))
temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rlogis(u,location, scale)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :Logistic(.(location),.(scale))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stllogis <-
function(N,t0,x0,T,location, scale,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( scale <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "scale > 0." ),icon="error"))
temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- exp(rlogis(u,location, scale))
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :Log*Logistic(.(location),.(scale))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stllogis3 <-
function(N,t0,x0,T,location, scale,thres,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( scale <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "scale > 0." ),icon="error"))
temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- exp(rlogis(u,location, scale))+thres
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :Log*Logistic[3](.(location),.(scale),.(thres))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stgumbel <-
function(N,t0,x0,T,location, scale,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( scale <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "scale > 0." ),icon="error"))

qgumbel <- function(p,scale=1,location=0,lower.tail=TRUE,log.p=FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	xF <- location-scale*log(-log(p))
	return(xF)
}
rgumbel <- function(n,scale=1,location=0)
{
	qgumbel(runif(n),scale,location)
}

temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rgumbel(u,scale,location)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :Gumbel(.(location),.(scale))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}

.Stgp <-
function(N,t0,x0,T,shape,scale,output=FALSE)
       {
if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))
if( shape <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "shape > 0" ),icon="error"))
if( scale <= 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "scale > 0" ),icon="error"))

qgp <-
function(p,shape=1,scale=1,lower.tail=TRUE,log.p=FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	xF <- scale/shape*(1-(1-p)^shape)
	return(xF)
}
rgp <- function(n,shape=1,scale=1)
{
 qgp(runif(n),shape,scale)
}
temps <- seq(t0,T,length=N+1)
dt <- (T-t0)/N
u <- runif(N)
x <- rgp(u,shape,scale)
y <- vector()
for (i in 1:N){if ( u[i] <= 0.5)
                   y[i] = -x[i] 
                   else
                   y[i] = x[i]}
procg <- c(x0,y*sqrt(dt))
plot(temps,procg,las=1,type="n",xlab="time",ylab=expression(X[t]),main=bquote("Stochastic process law" :G*Pareto(.(shape),.(scale))),font.main=2)
points(temps,procg,type="l",col="black",lwd=1)
mtext(paste("  Copyright 2012, USTHB. Algeria"),
      side = 1, line = 4, adj = 0.5, cex = .66)
X <- procg
time <- temps
Result <- data.frame(time,X)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "SP.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
                }
attach(Result)
}


.SLVM <-
function(N,t0,T,x0,y0,a,b,c,d,sigma,Step=FALSE,Output=FALSE)
       {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if (x0 <= 0 || y0 <= 0)
            stop(tkmessageBox(title="Error",message=paste( "x0,y0 > 0" ),icon="error"))

if( N <= 1 )   
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))

if (a <= 0 || b <= 0 || c <= 0 || d <= 0)
            stop(tkmessageBox(title="Error",message=paste( "a,b,c,d > 0" ),icon="error"))

if (sigma < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "sigma >= 0" ),icon="error"))

Dt <- (T-t0)/N

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
if(Step==FALSE){
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
windows()
par( mar =c(3 ,3 ,2,1))
par(mfrow=c(2,1))
plot(seq(t0,T,length=N+1),X,type="l",ylab=expression(x[t]),xlab="t",las=1,col="red")
points(t0,x0,pch=20,col="red")
mtext(expression(x[t]*minute-a*x[t]+b*x[t]*y[t]==epsilon[t]),line=0.1,adj=0,cex=1,col="red")
plot(seq(t0,T,length=N+1),Y,type="l",ylab=expression(y[t]),xlab="t",las=1,col="blue")
points(t0,y0,pch=20,col="blue")
mtext(expression(y[t]*minute-c*x[t]*y[t]+d*y[t]==epsilon[t]),line=0.1,adj=0,cex=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}

if(Step==TRUE){
plot(X,Y,type="n",xlab=expression(X[t]),ylab=expression(Y[t]),las=1)
points(x0,y0,pch=20,col="red")
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="black",lwd=1)}
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
windows()
par( mar =c(3 ,3 ,2,1))
par(mfrow=c(2,1))
plot(seq(t0,T,length=N+1),X,type="l",ylab=expression(x[t]),xlab="t",las=1,col="red")
points(t0,x0,pch=20,col="red")
mtext(expression(x[t]*minute-a*x[t]+b*x[t]*y[t]==epsilon[t]),line=0.1,adj=0,cex=1,col="red")
plot(seq(t0,T,length=N+1),Y,type="l",ylab=expression(y[t]),xlab="t",las=1,col="blue")
points(t0,y0,pch=20,col="blue")
mtext(expression(y[t]*minute-c*x[t]*y[t]+d*y[t]==epsilon[t]),line=0.1,adj=0,cex=1,col="blue")
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
}
time <- seq(t0,T,length=N+1)
Result <- data.frame(time,X,Y)
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(Output==TRUE){
write.table(Result, file = "SLVM.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
    }

.FBD <-
function(N,M,t0,T,x0,mu,sigma,output=FALSE)
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

FB <- function(N,T,t0,x0,mu,sigma)
   {
Dt <- (T-t0)/N
a <- expression(mu*x)
s <- expression(sigma*sqrt(x))
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
Q = sapply(rep(N,length=M),FB,T=T,t0=t0,x0=x0,mu,sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("Feller Branching Diffusion",adj=0.5,line=2.5,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==.(mu)*X[t]*dt+.(sigma)*sqrt(X[t])*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "FBD.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
    }

.WFD <-
function(N,M,t0,T,x0,gamma1,gamma2,sigma,output=FALSE)
       {

if( t0 >= T || t0 < 0 ) 
            stop(tkmessageBox(title="Error",message=paste( "T > t0 >= 0" ),icon="error"))

if ( x0 <= 0 || x0 >= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " 0 < x0 < 1" ),icon="error"))

if( N <= 1 ) 
            stop(tkmessageBox(title="Error",message=paste( " N must be very large N >>> 0" ),icon="error"))
 
if (M < 1)
            stop(tkmessageBox(title="Error",message=paste( "M must be  >= 1" ),icon="error"))

if (sigma <= 0 )
            stop(tkmessageBox(title="Error",message=paste( "Sigma > 0" ),icon="error"))

if ( gamma1 < 0 || gamma2 < 0)
            stop(tkmessageBox(title="Error",message=paste( "gamma1, gamma2 >= 0" ),icon="error"))

WF <- function(N,T,t0,x0,gamma1,gamma2,sigma)
   {
Dt <- (T-t0)/N
a <- expression(-gamma1*x+gamma2*(1-x))
s <- expression(sigma*sqrt(x*(1-x)))
DSx  <- D(s,"x")
A    <- function(t,x)  eval(a)
S    <- function(t,x)  eval(s)
Sx   <- function(t,x)  eval(DSx)
t = seq(t0,T,length=N+1)
w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
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
Q = sapply(rep(N,length=M),WF,T=T,t0=t0,x0=x0,gamma1,gamma2,sigma)
Q.mean <- apply(Q,1,mean)
r1 <- min(apply(Q,2,min))
r2 <- max(apply(Q,2,max))
plot(t,Q[,1],type="n",ylab=expression(X[t]),ylim=c(r1,r2),xlab="time",las=1)
mtext("Wright-Fisher Diffusion",adj=0.5,line=2.5,cex=1.2)
for (i in 1:M){points(t,Q[,i],type="l",col="black",lwd=1)}
mtext(bquote( dX[t]==(-.(gamma1)*X[t]+.(gamma2)*(1-X[t]))*dt+.(sigma)*sqrt(X[t]*(1-X[t]))*dW[t] ),cex=1.2,adj=0.5,line=0.25,col="red")
mtext(bquote(x[.(0)]==.(x0)),line=0.1,cex=0.9,adj=1,col="red")
mtext(bquote(t[0]==.(t0)),line=0.9,cex=0.9,adj=1,col="red")
mtext(bquote(Delta*t==.(Dt)),line=0.4,cex=1,adj=0,col="red")
if (M >=2){lines(t,Q.mean,lwd=2,col="red")
legend("topleft",border="gray",c("Average trajectory"),lty=c(1),col=c("red"),lwd=2)}
mtext(paste("  Copyright 2012, USTHB. Algeria"),side = 1, line = 4, adj = 0.5, cex = .66)
time <- t
X.mean <- Q.mean
X <- Q
Result <- data.frame(time,X)
if (M >=2) {Result <- data.frame(time,Q,X.mean)}
showData(Result , placement='+200-200', font = "Courier 11", body.textcolor = "black")
if(output==TRUE){
write.table(Result, file = "WFD.csv", sep = ";", col.names = TRUE,row.names = FALSE,
            qmethod = "double")
}
attach(Result)
    }

