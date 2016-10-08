

\name{fptsde2d}
\alias{fptsde2d}
\alias{fptsde2d.default}
\alias{summary.fptsde2d}
\alias{plot.fptsde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{First Passage Time in 2-Dim SDE}
\description{The (S3) generic function \code{fptsde2d} for simulate first-passage-time (f.p.t) in 2-dim stochastic differential equations.}
\usage{
fptsde2d(object, \dots)
\method{fptsde2d}{default}(object, boundary, \dots)
\method{summary}{fptsde2d}(object, digits=5, \dots)
\method{plot}{fptsde2d}(x, \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde2d}}.}
  \item{boundary}{an \code{\link[base]{expression}} of a constant or time-dependent boundary.}
  \item{digits}{integer, used for number formatting.}    
  \item{x}{an object inheriting from class \code{\link{fptsde2d}}.}
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{fptsde1d} returns a random variable \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)))} "first passage time", is defined as : 
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \geq S(t) \},\quad if \quad X(t_{0}) < S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) >= S(t)}, if X(t0) < S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \geq S(t) \},\quad if \quad Y(t_{0}) < S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) >= S(t)}, if Y(t0) < S(t0)}
and:
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \leq S(t) \},\quad if \quad X(t_{0}) > S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) <= S(t)}, if X(t0) > S(t0)}
\deqn{\tau_{(Y(t),S(t))} = \{ t \geq 0 ; Y_{t} \leq S(t) \},\quad if \quad Y(t_{0}) > S(t_{0})}{tau(Y(t),S(t))={t>=0; Y(t) <= S(t)}, if Y(t0) > S(t0)}
with \eqn{S(t)}{S(t)} is through a continuous boundary (barrier).
}
\value{
\code{fptsde2d} returns an object inheriting from \code{\link[base]{class}} \code{"fptsde2d"}.
\item{fptx, fpty }{a vector of couple 'fpt' \eqn{(\tau_{(X(t),S(t))},\tau_{(Y(t),S(t))})}{(tau(X(t),S(t)),tau(Y(t),S(t)))}.}
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
\code{\link[Sim.DiffProc]{fptsde1d}} for simulation fpt in sde 1-dim. \code{\link{fptsde3d}} for simulation fpt in sde 3-dim. 

\code{\link[fptdApprox]{FPTL}} for computes values of the first passage time location (FPTL) function, and \code{\link[fptdApprox]{Approx.fpt.density}} 
for approximate first-passage-time (f.p.t.) density in package "fptdApprox". 
}
\examples{

## dX(t) = 5*(-1-Y(t))*X(t) * dt + 0.5 * dW1(t)          
## dY(t) = 5*(-1-X(t))*Y(t) * dt + 0.5 * dW2(t)
## x0 = 2, y0 = -2, and barrier -3+5*t.
## W1(t) and W2(t) two independent Brownian motion
set.seed(1234)

fx <- expression(5*(-1-y)*x)
gx <- expression(0.5)
fy <- expression(5*(-1-x)*y)
gy <- expression(0.5)
St <- expression(-1+5*t)

mod2d <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,x0=2,y0=-2,M=50)
fptmod2d <- fptsde2d(mod2d,boundary=St)
fptmod2d
summary(fptmod2d)
plot(fptmod2d )

##

X1 <- cbind(fptmod2d$fptx,fptmod2d$fpty)
## library(sm)
## sm.density(X1,display="persp")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
\keyword{mts}
