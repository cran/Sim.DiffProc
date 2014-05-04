library(Sim.DiffProc)


## Itô sde

f <- expression( -3*(1+x) )
g <- expression( 0.5*x )
res <- fptsde1d(drift=f,diffusion=g,x0=1,c=0,M=100,N=1000)
res
summary(res)
bconfint(res,level=0.95)
moment(res,order=c(2,3,4,5))
plot(density(res$tau[!is.na(res$tau)]))

## Stratonovich sde

f <- expression( -3*(1+x) )
g <- expression( 0.5*x )
res <- fptsde1d(drift=f,diffusion=g,x0=1,c=0,M=100,N=1000,type="str")
res
summary(res)
bconfint(res,level=0.95)
moment(res,order=c(2,3,4,5))
plot(density(res$tau[!is.na(res$tau)]))
