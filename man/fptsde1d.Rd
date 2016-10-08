

\name{fptsde1d}
\alias{fptsde1d}
\alias{fptsde1d.default}
\alias{summary.fptsde1d}
\alias{plot.fptsde1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{First Passage Time in 1-Dim SDE}
\description{The (S3) generic function \code{fptsde1d} for simulate first-passage-time (f.p.t) in 1-dim stochastic differential equations.}
\usage{
fptsde1d(object, \dots)
\method{fptsde1d}{default}(object, boundary, \dots)
\method{summary}{fptsde1d}(object, digits=5, \dots)
\method{plot}{fptsde1d}(x, \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde1d}}.}
  \item{boundary}{an \code{\link[base]{expression}} of a constant or time-dependent boundary.}
  \item{digits}{integer, used for number formatting.}    
  \item{x}{an object inheriting from class \code{\link{fptsde1d}}.}
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{fptsde1d} returns a random variable \eqn{\tau_{(X(t),S(t))}}{tau(X(t),S(t))} "first passage time", is defined as : 
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \geq S(t) \},\quad if \quad X(t_{0}) < S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) >= S(t)}, if X(t0) < S(t0)}
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \leq S(t) \},\quad if \quad X(t_{0}) > S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) <= S(t)}, if X(t0) > S(t0)}
with \eqn{S(t)}{S(t)} is through a continuous boundary (barrier).
}
\value{
\code{fptsde1d} returns an object inheriting from \code{\link[base]{class}} \code{"fptsde1d"}.
\item{fpt }{numeric vector of fpt (first passage time).}
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
\code{\link[Sim.DiffProc]{fptsde2d}} and \code{\link[Sim.DiffProc]{fptsde3d}} simulation fpt for 2 and 3-dim SDE.

\code{\link[fptdApprox]{FPTL}} for computes values of the first passage time location (FPTL) function, and \code{\link[fptdApprox]{Approx.fpt.density}} 
for approximate first-passage-time (f.p.t.) density in package "fptdApprox". 
}
\examples{

## Example 1: Ito SDE
## dX(t) = -4*X(t) *dt + 0.5*dW(t)
## S(t) = 0 (constant boundary)
set.seed(1234)

f <- expression( -4*x )
g <- expression( 0.5 )
St <- expression(0) 
mod <- snssde1d(drift=f,diffusion=g,x0=2,M=50)
fptmod <- fptsde1d(mod, boundary=St)
fptmod
summary(fptmod)
plot(fptmod)
plot(density(fptmod$fpt[!is.na(fptmod$fpt)]),main="Kernel Density of a First-Passage-Time")

## Example 2: Ito SDE
## X(t) Brownian motion
## S(t) = 0.3+0.2*t (time-dependent boundary)
set.seed(1234)

f <- expression( 0 )
g <- expression( 1 )
St <- expression(0.5-0.5*t) 
mod1 <- snssde1d(drift=f,diffusion=g,M=40)
fptmod1 <- fptsde1d(mod1,boundary=St)
summary(fptmod1)
plot(fptmod1)
plot(density(fptmod1$fpt[!is.na(fptmod1$fpt)]),main="Kernel Density of a First-Passage-Time")

## Example 3: Stratonovich SDE
## dX(t) = 0.5*X(t)*t *dt + sqrt(1+X(t)^2) o dW(t)
## S(t) = -0.5*sqrt(t) + exp(t^2) (time-dependent boundary)
set.seed(1234)

f <- expression( 0.5*x*t )
g <- expression( sqrt(1+x^2) )
St <- expression(-0.5*sqrt(t)+exp(t^2))
mod2 <- snssde1d(drift=f,diffusion=g,x0=2,M=40,type="srt")
fptmod2 <- fptsde1d(mod2,boundary=St)
summary(fptmod2)
plot(density(fptmod2$fpt[!is.na(fptmod2$fpt)]),main="Kernel Density of a First-Passage-Time")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
