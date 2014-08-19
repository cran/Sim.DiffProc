library(Sim.DiffProc)


### 2-dim Ito bridge sde

fx <- expression(4*(-2-x))
gx <- expression(y)
fy <- expression(4*(-2-y))
gy <- expression(x)

res <- bridgesde2d(x0=c(1,1),y=c(1,1),driftx=fx,diffx=gx,drifty=fy,diffy=gy)
res
plot(res)
dev.new()
plot(res,plot.type="single")
dev.new()
plot2d(res,type="n")
points2d(res,col=rgb(0,100,0,50,maxColorValue=255), pch=16)

