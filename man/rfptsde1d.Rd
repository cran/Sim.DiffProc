

\name{rfptsde1d}
\alias{rfptsde1d}
\alias{rfptsde1d.default}
\alias{dfptsde1d}
\alias{dfptsde1d.default}
\alias{print.dfptsde1d}
\alias{plot.dfptsde1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density and random generation for first passage time in 1-Dim SDE}
\description{Density and random generation for first-passage-time (f.p.t) in 1-dim stochastic differential equations.}
\usage{
rfptsde1d(object, \dots)
dfptsde1d(object, \dots)
\method{rfptsde1d}{default}(object, boundary, \dots)
\method{dfptsde1d}{default}(object, boundary, \dots)
\method{plot}{dfptsde1d}(x, \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde1d}}.}
  \item{boundary}{an \code{\link[base]{expression}} of a constant or time-dependent boundary.}  
  \item{x}{an object inheriting from class \code{dfptsde1d}.}
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{rfptsde1d} returns a random variable \eqn{\tau_{(X(t),S(t))}}{tau(X(t),S(t))} "first passage time", is defined as : 
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \geq S(t) \},\quad if \quad X(t_{0}) < S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) >= S(t)}, if X(t0) < S(t0)}
\deqn{\tau_{(X(t),S(t))} = \{ t \geq 0 ; X_{t} \leq S(t) \},\quad if \quad X(t_{0}) > S(t_{0})}{tau(X(t),S(t))={t>=0; X(t) <= S(t)}, if X(t0) > S(t0)}

And \code{dfptsde1d} returns a kernel density approximation for \eqn{\tau_{(X(t),S(t))}}{tau(X(t),S(t))} "first passage time".

with \eqn{S(t)}{S(t)} is through a continuous boundary (barrier).
}
\value{
\code{dfptsde1d} {gives the density estimate of fpt.}
\code{rfptsde1d} {generates random of fpt.}
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
\code{\link[Sim.DiffProc]{rfptsde2d}} and \code{\link[Sim.DiffProc]{rfptsde3d}} simulation fpt for 2 and 3-dim SDE.

\code{\link[fptdApprox]{FPTL}} for computes values of the first passage time location (FPTL) function, and \code{\link[fptdApprox]{Approx.fpt.density}} 
for approximate first-passage-time (f.p.t.) density in package "fptdApprox". 

\code{\link[DiffusionRgqd]{GQD.TIpassage}} for compute the First Passage Time Density of a GQD With Time Inhomogeneous Coefficients.

}
\examples{

## Example 1: Ito SDE
## dX(t) = -4*X(t) *dt + 0.5*dW(t)
## S(t) = 0 (constant boundary)
set.seed(1234)

# SDE 1d
f <- expression( -4*x )
g <- expression( 0.5 )
mod <- snssde1d(drift=f,diffusion=g,x0=2,M=50)

# boundary
St <- expression(0) 
# random
out <- rfptsde1d(mod, boundary=St)
summary(out)
# density approximate
den <- dfptsde1d(mod, boundary=St)
plot(den)

## Example 2: Stratonovich SDE
## dX(t) = 0.5*X(t)*t *dt + sqrt(1+X(t)^2) o dW(t)
## S(t) = -0.5*sqrt(t) + exp(t^2) (time-dependent boundary)
set.seed(1234)

# SDE 1d
f <- expression( 0.5*x*t )
g <- expression( sqrt(1+x^2) )
mod2 <- snssde1d(drift=f,diffusion=g,x0=2,M=40,type="srt")

# boundary
St <- expression(-0.5*sqrt(t)+exp(t^2))
# random
out2 <- rfptsde1d(mod2,boundary=St)
summary(out2)
# density approximate
plot(dfptsde1d(mod2,boundary=St,bw='ucv'))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
