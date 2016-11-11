

\name{rfptsde3d}
\alias{rfptsde3d}
\alias{rfptsde3d.default}
\alias{dfptsde3d}
\alias{dfptsde3d.default}
\alias{print.dfptsde3d}
\alias{plot.dfptsde3d}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density and random generation for first passage time in 3-Dim SDE's}
\description{Density and random generation for the marginal for first-passage-time (f.p.t) in 3-dim stochastic differential equations.}
\usage{
rfptsde3d(object, \dots)
dfptsde3d(object, ...)
\method{rfptsde3d}{default}(object, boundary, \dots)
\method{dfptsde3d}{default}(object, boundary, ...)
\method{plot}{dfptsde3d}(x,hist=FALSE,...)					   
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde3d}}.}
  \item{boundary}{an \code{\link[base]{expression}} of a constant or time-dependent boundary.}
  \item{x}{an object inheriting from class \code{dfptsde3d}.}
  \item{hist}{if \code{hist=TRUE} plot histogram.}   
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{rfptsde3d} returns a random variable \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))},\tau_{(Z(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)),tau(Z(t),S(t)))} "first passage time", is defined as : 
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \geq S(t) \},\quad if \quad X(t_{0}) < S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) >= S(t)}, if X(t0) < S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \geq S(t) \},\quad if \quad Y(t_{0}) < S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) >= S(t)}, if Y(t0) < S(t0)}
\deqn{\tau_{(Z(t),S(t))} = \{ t \geq 0 ; Z_{t} \geq S(t) \},\quad if \quad Z(t_{0}) < S(t_{0})}{tau(Z(t),S(t))={t>=0; Z(t) >= S(t)}, if Z(t0) < S(t0)}
and:
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \leq S(t) \},\quad if \quad X(t_{0}) > S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) <= S(t)}, if X(t0) > S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \leq S(t) \},\quad if \quad Y(t_{0}) > S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) <= S(t)}, if Y(t0) > S(t0)}
\deqn{\tau_{(Z(t),S(t))} = \{ t \geq 0 ; Z_{t} \leq S(t) \},\quad if \quad Z(t_{0}) > S(t_{0})}{tau(Z(t),S(t))={t>=0; Z(t) <= S(t)}, if Z(t0) > S(t0)}

And \code{dfptsde3d} returns a marginal kernel density approximation for \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))},\tau_{(Z(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)),tau(Z(t),S(t)))} "first passage time".

with \eqn{S(t)}{S(t)} is through a continuous boundary (barrier).
}
\value{
\code{dfptsde3d} {gives the marginal kernel density approximation for fpt.}
\code{rfptsde3d} {generates random of fpt.}
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
\code{\link{rfptsde1d}} for simulation fpt in sde 1-dim. \code{\link{rfptsde2d}} for simulation fpt in sde 2-dim.

\code{\link[fptdApprox]{FPTL}} for computes values of the first passage time location (FPTL) function, and \code{\link[fptdApprox]{Approx.fpt.density}} 
for approximate first-passage-time (f.p.t.) density in package "fptdApprox".  
}
\examples{

## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 * dW1(t) 
## dY(t) = 4*(1-Y(t)) *X(t) dt + 0.2 * dW2(t) 
## dZ(t) = 4*(1-Z(t)) *Y(t) dt + 0.2 * dW3(t) 
## x0 = 0, y0 = -2, z0 = 0, and barrier -3+5*t.       
## W1(t), W2(t) and W3(t) three independent Brownian motion      
set.seed(1234)

# SDE's 3d

fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,x0=c(2,-2,0),M=50)

# boundary 
St <- expression(-3+5*t)

# random

out <- rfptsde3d(mod3d,boundary=St)
out
summary(out)

# Marginal density

denM <- dfptsde3d(mod3d,boundary=St)
denM
plot(denM)


## for Joint density see package 'sm' or 'ks'

## library(sm)
## sm.density(out,display="rgl")

## library(ks)
## fhat <- kde(x=out)
## plot(fhat, drawpoints=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
\keyword{mts}
