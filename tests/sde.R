library(Sim.DiffProc)


## Itô sde 3-dim

fx <- expression(y)
gx <- expression(z)
fy <- expression(0)
gy <- expression(1)
fz <- expression(0)
gz <- expression(1)

res <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,N=1000)
res
Sim <- res$XYZ
dev.new()
plot(res)
dev.new()
plot(res,plot.type="single")
dev.new()
plot3D(res,display="persp")
plot3D(res,display="rgl") 

## Stratonovich sde 3-dim

fx <- expression(2*(3-x))
gx <- expression(y+z)
fy <- expression(2*(3-y))
gy <- expression(x+z)
fz <- expression(2*(3-z))
gz <- expression(x+y)

res <- snssde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,
                ,N=1000,type="str")
res
dev.new()
plot(res,plot.type="single")
dev.new()
plot(res)
dev.new()
plot3D(res,display="rgl") 
plot3D(res,display="persp")
