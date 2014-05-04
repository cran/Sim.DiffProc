\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{fptsde1d}
\alias{fptsde1d}
\alias{fptsde1d.default}
\alias{summary.fptsde1d}
\alias{mean.fptsde1d}
\alias{median.fptsde1d}
\alias{quantile.fptsde1d}
\alias{kurtosis.fptsde1d}
\alias{skewness.fptsde1d}
\alias{moment.fptsde1d}
\alias{bconfint.fptsde1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{First Passage Time in 1-Dim SDE}
\description{The (S3) generic function \code{fptsde1d} for simulate first-passage-time (f.p.t) in 1-dim stochastic differential equations.}
\usage{
fptsde1d(N, \dots)
\method{fptsde1d}{default}(N = 100, M = 10, x0 = 0, t0 = 0, T = 1, Dt, c = 0.5, 
   drift, diffusion, alpha = 0.5, mu = 0.5, type = c("ito", "str"), 
   method = c("euler", "milstein", "predcorr", "smilstein", "taylor", 
   "heun", "rk1", "rk2", "rk3"), \dots)
						  
						  
\method{summary}{fptsde1d}(object, \dots)
\method{mean}{fptsde1d}(x, \dots)
\method{median}{fptsde1d}(x, \dots)
\method{quantile}{fptsde1d}(x, \dots)
\method{kurtosis}{fptsde1d}(x, \dots)
\method{skewness}{fptsde1d}(x, \dots)
\method{moment}{fptsde1d}(x, order = 2, \dots)
\method{bconfint}{fptsde1d}(x, level=0.95, \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{size of sde.}
  \item{M}{size of fpt.}
  \item{x0}{initial value of the process at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{c}{boundary or barrier (threshold).}
  \item{drift}{drift coefficient: an \code{\link[base]{expression}} of two variables \code{t} and \code{x}.}
  \item{diffusion}{diffusion coefficient: an \code{\link[base]{expression}} of two variables \code{t} and \code{x}.}
  \item{alpha}{weight \code{alpha} of the predictor-corrector scheme; the default \code{alpha = 0.5}.}
  \item{mu}{weight \code{mu} of the predictor-corrector scheme; the default \code{mu = 0.5}.}
  \item{type}{sde of the type Ito or Stratonovich.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde1d}}.}
  \item{x, object}{an object inheriting from class \code{"fptsde1d"}.}
  \item{order}{order of moment.}
  \item{level}{the confidence level required.}    
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{fptsde1d} returns a random variable \eqn{\tau_{c}}{tau(c)} "first passage time", defined by : 
\deqn{\tau_{c} = \{ t \geq 0 ; X_{t} \geq c \},\quad if \quad (x_{0} < c)}{tau(c)={t>=0; X(t) >= c}, if (x0 < c)}
\deqn{\tau_{c} = \{ t \geq 0 ; X_{t} \leq c \},\quad if \quad (x_{0} > c)}{tau(c)={t>=0; X(t) <= c}, if (x0 > c)}
with \eqn{c}{c} is a fixed boundary or barrier.
}
\value{
\code{fptsde1d} returns an object inheriting from \code{\link[base]{class}} \code{"fptsde1d"}.
\item{tau }{numeric vector of fpt.}
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
   
   Roman, R.P., Serrano, J. J., Torres, F. (2013). 
   \CRANpkg{fptdApprox}: Approximation of first-passage-time densities for diffusion processes. 
   \emph{\R package version 2.0}. 
  
  Gardiner, C. W. (1997). 
  \emph{Handbook of Stochastic Methods}.
  Springer-Verlag, New York.
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link[Sim.DiffProc]{fptsde2d}} and \code{\link[Sim.DiffProc]{fptsde3d}} simulation fpt for sde 2 and 3-dim.

\code{\link[fptdApprox]{FPTL}} in package \CRANpkg{fptdApprox} for computes values of the first passage time location. 
}
\examples{

## dX(t) = -3*(1+X(t)) *dt + 0.5*X(t) *dW(t)
## tau(c=0) = inf{t>=0, X(t) <= 0}

f <- expression( -3*(1+x) )
g <- expression( 0.5*x )
res <- fptsde1d(drift=f,diffusion=g,x0=1,c=0,M=100,N=1000)
res
summary(res)
bconfint(res,level=0.95)
moment(res,order=c(2,3,4,5))
plot(density(res$tau[!is.na(res$tau)]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
