library(Sim.DiffProc)

## Converting Sim.DiffProc Objects to LaTeX

# Example 1

f <- expression(-mu1 * x) 
g <- expression(mu2 * sqrt(x)) 
TEX.sde(object = c(drift = f, diffusion = g))

# Example 2

f <- expression(mu1*cos(mu2+z),mu1*sin(mu2+z),0) 
g <- expression(sigma,sigma,alpha) 
TEX.sde(object = c(drift = f, diffusion = g))

## LaTeX mathematic for object of class 'MEM.sde'
## Copy and paste the following output in your LaTeX file

# Example 3

mem.mod3d <- MEM.sde(drift = f, diffusion = g)
TEX.sde(object = mem.mod3d)
