## 22:45 20/03/2010
## Guidoum Arsalane(PG-PS/USTHB)
## Calculate empirical covariance of the Brownian Motion (cov(t,s)= C*min(t,s))
## N Size of process
## M Numbers of the trajectories simulated
## T final time
## C Constant (if C=1 ==> Standard Brownian Motion)

 
N = 100
M = 300
T = 10
C = 1
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
