## ----setup, echo = F, message = F, results = 'hide',screenshot.force=FALSE----
library(Sim.DiffProc)
library(knitr)
knitr::opts_chunk$set(comment="",prompt=TRUE, fig.show='hold',warning=FALSE, message=FALSE)
options(prompt="R> ",scipen=16,digits=5,warning=FALSE, message=FALSE,width = 80)

## -----------------------------------------------------------------------------
mu=1;sigma=0.5;theta=2
x0=0;y0=0;init=c(x0,y0)
f <- expression(1/mu*(theta-x), x)  
g <- expression(sqrt(sigma),0)
mod2d <- snssde2d(drift=f,diffusion=g,M=500,Dt=0.015,x0=c(x=0,y=0))
## true values of first and second moment at time 10
Ex <- function(t) theta+(x0-theta)*exp(-t/mu)
Vx <- function(t) 0.5*sigma*mu *(1-exp(-2*(t/mu)))
Ey <- function(t) y0+theta*t+(x0-theta)*mu*(1-exp(-t/mu))
Vy <- function(t) sigma*mu^3*((t/mu)-2*(1-exp(-t/mu))+0.5*(1-exp(-2*(t/mu))))
covxy <- function(t) 0.5*sigma*mu^2 *(1-2*exp(-t/mu)+exp(-2*(t/mu)))
tvalue = list(m1=Ex(15),m2=Ey(15),S1=Vx(15),S2=Vy(15),C12=covxy(15))
## function of the statistic(s) of interest.
sde.fun2d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),mean(d$y),var(d$x),var(d$y),cov(d$x,d$y)))
}
## Parallel Monte-Carlo of 'OUI' at time 10
mcm.mod2d = MCM.sde(mod2d,statistic=sde.fun2d,time=15,R=10,exact=tvalue,parallel="snow",ncpus=2)
mcm.mod2d$MC


## -----------------------------------------------------------------------------
TEX.sde(object = mcm.mod2d, booktabs = TRUE, align = "r", caption ="LaTeX 
          table for Monte Carlo results generated by `TEX.sde()` method.")

## ----echo=FALSE---------------------------------------------------------------
kable(mcm.mod2d$MC, format = "html",booktabs = TRUE,align = "r", caption ="LaTeX 
          table for Monte Carlo results generated by `TEX.sde()` method.")

## -----------------------------------------------------------------------------
mem.oui <- MEM.sde(drift = f, diffusion = g)
mem.oui

## -----------------------------------------------------------------------------
TEX.sde(object = mem.oui)


## -----------------------------------------------------------------------------
f <- expression((alpha*x *(1 - x / beta)- delta * x^2 * y / (kappa + x^2)),
                (gamma * x^2 * y / (kappa + x^2) - mu * y^2)) 
g <- expression(sqrt(sigma1)*x*(1-y), abs(sigma2)*y*(1-x))  
TEX.sde(object=c(drift = f, diffusion = g))

