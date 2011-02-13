## Sat Feb 05 16:24:36 2011
## Guidoum Arsalane (PG-PS/USTHB)


N=1000;t0=0;Dt=0.001;T=1;X0=2;Y0=1;v=0.2;K=3;Sigma=0.2;

drifx     <- expression( (-K*x) / (x^2+y^2) )
drify     <- expression( (-K*y) / (x^2+y^2) )
diff      <- expression( Sigma ) 

Ax    <- function(t,x,y)  eval(drifx)
Ay    <- function(t,x,y)  eval(drify)
S     <- function(t,x,y)  eval(diff)

if(missing(Dt)){t <- seq (t0 ,T, length =N+1)} else {t <- c(t0 ,t0+ cumsum(rep(Dt,N)))
                T <- t[N +1]}

Dt= (T-t0)/N 
ux = runif(N,0,1)
ox = rep(1,N)
ox [ which(ux < 0.5) ] = -1
wx = cumsum(c(0,ox))*sqrt((T-t0)/N)
Dx   <- diff(wx)

uy = runif(N,0,1)
oy = rep(1,N)
oy [ which(uy < 0.5) ] = -1
wy = cumsum(c(0,oy))*sqrt((T-t0)/N)
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
for (i in 1:min(n)){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="blue",lwd=2,panel.frist=grid(col="gray"))}}else{
for (i in 1:N){lines(c(X[i],X[i+1]),c(Y[i],Y[i+1]),type="l",col="blue",lwd=2,panel.frist=grid(col="gray"))}}

if (length(n) > 0 ){points(X[min(n)],Y[min(n)],type="p",col="red",cex=1.2,pch="*")
                    text(X[min(n)],Y[min(n)], expression(tau[v]^(1)), col=2, adj=c(-.1,-.1),cex = 1.2)
                    mtext(bquote(tau[v]^(1)== .(t[min(n)])),line=0.5,adj=0.45,cex=1,col="red")}
mtext(paste("USTHB,Faculty of Mathematics,Department of Probabilities and Statistics,Algeria Sat Feb 05 16:24:36 2011"),side = 1, line = 4, adj = 0.5, cex = .66)

