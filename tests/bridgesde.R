library(Sim.DiffProc)


### 1-dim Ito bridge sde

fx <- expression(4*(-2-x))
gx <- expression(y)
fy <- expression(4*(-2-y))
gy <- expression(x)

XY <- bridgesde2d(x0=c(1,1),y=c(1,1),driftx=fx,diffx=gx,drifty=fy,diffy=gy)
XY
summary(XY)
plot(XY)
dev.new()
plot(XY,plot.type="single")
dev.new()
plot2d(XY,type="n")
points2d(XY,col=rgb(0,100,0,50,maxColorValue=255), pch=16)
