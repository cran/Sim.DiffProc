

\name{rfptsde2d}
\alias{rfptsde2d}
\alias{rfptsde2d.default}
\alias{dfptsde2d}
\alias{dfptsde2d.default}
\alias{print.dfptsde2d}
\alias{plot.dfptsde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density and random generation for first passage time in 2-Dim SDE's}
\description{Density and random generation for the joint or marginal for first-passage-time (f.p.t) in 2-dim stochastic differential equations.}
\usage{
rfptsde2d(object, \dots)
dfptsde2d(object, \dots)
\method{rfptsde2d}{default}(object, boundary, \dots)
\method{dfptsde2d}{default}(object, boundary, pdf=c("Joint","Marginal"), \dots)
\method{plot}{dfptsde2d}(x,display=c("persp","rgl","image","contour"),
                         hist=FALSE, \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde2d}}.}
  \item{boundary}{an \code{\link[base]{expression}} of a constant or time-dependent boundary.} 
  \item{pdf}{probability density function \code{Joint} or \code{Marginal}.}    
  \item{x}{an object inheriting from class \code{dfptsde2d}.}
  \item{display}{display plots.} 
  \item{hist}{if \code{hist=TRUE} plot histogram.} 
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{rfptsde1d} returns a random variable \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)))} "first passage time", is defined as : 
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \geq S(t) \},\quad if \quad X(t_{0}) < S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) >= S(t)}, if X(t0) < S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \geq S(t) \},\quad if \quad Y(t_{0}) < S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) >= S(t)}, if Y(t0) < S(t0)}
and:
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \leq S(t) \},\quad if \quad X(t_{0}) > S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) <= S(t)}, if X(t0) > S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \leq S(t) \},\quad if \quad Y(t_{0}) > S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) <= S(t)}, if Y(t0) > S(t0)}

And \code{dfptsde2d} returns a kernel density approximation for \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)))} "first passage time".

with \eqn{S(t)}{S(t)} is through a continuous boundary (barrier).
}
\value{
\code{dfptsde2d} {gives the kernel density approximation for fpt.}
\code{rfptsde2d} {generates random of fpt.}
}
\references{
  Argyrakisa, P. and G.H. Weiss (2006).
  A first-passage time problem for many random walkers.
  \emph{Physica A}. \bold{363}, 343--347.
  
  Aytug H., G. J. Koehler (2000). 
  New stopping criterion for genetic algorithms. 
  \emph{European Journal of Operational Research}, \bold{126}, 662--674.
  
  Boukhetala, K. (1996)
  Modelling and simulation of a dispersion pollutant with attractive centre.
  ed by Computational Mechanics Publications, Southampton ,U.K and Computational Mechanics Inc, Boston, USA, 245--252.
  
  Boukhetala, K. (1998a). 
  Estimation of the first passage time distribution for a simulated diffusion process.
  \emph{Maghreb Math.Rev}, \bold{7}(1), 1--25.
  
  Boukhetala, K. (1998b). 
  Kernel density of the exit time in a simulated diffusion. 
  \emph{les Annales Maghrebines De L ingenieur}, \bold{12}, 587--589.
  
  Ding, M. and G. Rangarajan. (2004). 
  First Passage Time Problem: A Fokker-Planck Approach. 
  \emph{New Directions in Statistical Physics}. ed by L. T. Wille. Springer. 31--46.
  
   Roman, R.P., Serrano, J. J., Torres, F. (2008).
   First-passage-time location function: Application to determine first-passage-time densities in diffusion processes.
   \emph{Computational Statistics and Data Analysis}. \bold{52}, 4132--4146.
  
   Roman, R.P., Serrano, J. J., Torres, F. (2012). 
   An R package for an efficient approximation of first-passage-time densities for diffusion processes based on the FPTL function. 
   \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.
    
  Gardiner, C. W. (1997). 
  \emph{Handbook of Stochastic Methods}.
  Springer-Verlag, New York.
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link[Sim.DiffProc]{rfptsde1d}} for simulation fpt in sde 1-dim. \code{\link{rfptsde3d}} for simulation fpt in sde 3-dim. 

\code{\link[fptdApprox]{FPTL}} for computes values of the first passage time location (FPTL) function, and \code{\link[fptdApprox]{Approx.fpt.density}} 
for approximate first-passage-time (f.p.t.) density in package "fptdApprox". 
}
\examples{

## dX(t) = 5*(-1-Y(t))*X(t) * dt + 0.5 * dW1(t)          
## dY(t) = 5*(-1-X(t))*Y(t) * dt + 0.5 * dW2(t)
## x0 = 2, y0 = -2, and barrier -3+5*t.
## W1(t) and W2(t) two independent Brownian motion
set.seed(1234)

# SDE's 2d
fx <- expression(5*(-1-y)*x , 5*(-1-x)*y)
gx <- expression(0.5 , 0.5)
mod2d <- snssde2d(drift=fx,diffusion=gx,x0=c(2,-2),M=50)

# boundary

St <- expression(-1+5*t)

# random fpt

out <- rfptsde2d(mod2d,boundary=St)
summary(out)

# Marginal density 

denM <- dfptsde2d(mod2d,boundary=St,pdf="M")
denM
plot(denM)

# Joint density

denJ <- dfptsde2d(mod2d,boundary=St,pdf="J")
denJ
plot(denJ)
plot(denJ,display="image")
plot(denJ,display="image",drawpoints=TRUE,cex=0.5,pch=19,col.pt='green')
plot(denJ,display="contour")
plot(denJ,display="contour",color.palette=colorRampPalette(c('white','green','blue','red')))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
\keyword{mts}
