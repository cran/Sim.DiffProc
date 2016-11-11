library(Sim.DiffProc)


## Example 1: 

f <- expression( 0.5*x*t )
g <- expression( sqrt(1+x^2) )
St <- expression(-0.5*sqrt(t)+exp(t^2))
mod2 <- snssde1d(drift=f,diffusion=g,x0=2,M=40,type="srt")
fptmod2 <- rfptsde1d(mod2,boundary=St)
summary(fptmod2)
plot(dfptsde1d(mod2,boundary=St))

## Example 2: 

fx <- expression(4*(-1-x)*y, 4*(1-y)*x)
gx <- rep(expression(0.2),2)

St <- expression(-3+5*t)

mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(2,-2),M=50)
fptmod2d <- rfptsde2d(mod2d,boundary=St)
fptmod2d
summary(fptmod2d)
plot(dfptsde2d(mod2d,boundary=St))