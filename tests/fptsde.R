library(Sim.DiffProc)


## Example 1: 

f <- expression( 0.5*x*t )
g <- expression( sqrt(1+x^2) )
St <- expression(-0.5*sqrt(t)+exp(t^2))
mod2 <- snssde1d(drift=f,diffusion=g,x0=2,M=40,type="srt")
fptmod2 <- fptsde1d(mod2,boundary=St)
summary(fptmod2)
plot(density(fptmod2$fpt[!is.na(fptmod2$fpt)]),main="Kernel Density of a First-Passage-Time")

## Example 2: 

fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)
fz <- expression(4*(1-z)*y)
gz <- expression(0.2)

St <- expression(-3+5*t)

mod3d <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,
                x0=2,y0=-2,z0=0,M=50)
fptmod3d <- fptsde3d(mod3d,boundary=St)
fptmod3d
summary(fptmod3d)
plot(fptmod3d,union=TRUE)
dev.new()
plot(fptmod3d,union=FALSE)