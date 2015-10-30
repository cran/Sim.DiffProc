\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{bridgesde1d}
\alias{bridgesde1d}
\alias{bridgesde1d.default}
\alias{print.bridgesde1d}
\alias{time.bridgesde1d}
\alias{mean.bridgesde1d}
\alias{median.bridgesde1d}
\alias{quantile.bridgesde1d}
\alias{kurtosis.bridgesde1d}
\alias{skewness.bridgesde1d}
\alias{moment.bridgesde1d}
\alias{bconfint.bridgesde1d}
\alias{plot.bridgesde1d}
\alias{points.bridgesde1d}
\alias{lines.bridgesde1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 1-Dim Diffusion Bridge}
\description{The (S3) generic function \code{bridgesde1d} for simulation of 1-dim diffusion bridge.}
\usage{
bridgesde1d(N, \dots)
\method{bridgesde1d}{default}(N = 1000, M=1, x0 = 0, y = 0, t0 = 0, T = 1, Dt, 
   drift, diffusion, alpha = 0.5, mu = 0.5, type = c("ito", "str"), 
   method = c("euler", "milstein", "predcorr", "smilstein", "taylor", 
   "heun", "rk1", "rk2", "rk3"), \dots)
							  
							  
\method{time}{bridgesde1d}(x, \dots)
\method{mean}{bridgesde1d}(x, \dots)
\method{median}{bridgesde1d}(x, \dots)
\method{quantile}{bridgesde1d}(x, \dots)
\method{kurtosis}{bridgesde1d}(x, \dots)
\method{skewness}{bridgesde1d}(x, \dots)
\method{moment}{bridgesde1d}(x, order = 2, \dots)
\method{bconfint}{bridgesde1d}(x, level=0.95, \dots)
\method{plot}{bridgesde1d}(x, \dots)
\method{lines}{bridgesde1d}(x, \dots)
\method{points}{bridgesde1d}(x, \dots)							  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories.}
  \item{x0}{initial value of the process at time \code{t0}.}
  \item{y}{terminal value of the process at time \code{T}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{drift}{drift coefficient: an \code{\link[base]{expression}} of two variables \code{t} and \code{x}.}
  \item{diffusion}{diffusion coefficient: an \code{\link[base]{expression}} of two variables \code{t} and \code{x}.}
  \item{alpha, mu}{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation diffusion bridge of Ito type, else \code{type="str"} simulation diffusion bridge of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde1d}}.}
  \item{x}{an object inheriting from class \code{"bridgesde1d"}.}
  \item{order}{order of moment.}
  \item{level}{the confidence level required.} 
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{bridgesde1d} returns a trajectory of the diffusion bridge starting at \code{x} at time \code{t0} and ending at \code{y} at time \code{T}.

The methods of approximation are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong 
and the weak (orders of) convergence. The \code{method} of simulation can be one among: \code{Euler-Maruyama Order 0.5}, \code{Milstein Order 1}, \code{Milstein Second-Order}, 
\code{Predictor-Corrector method}, \code{Ito-Taylor Order 1.5}, \code{Heun Order 2} and \code{Runge-Kutta Order 1, 2 and 3}.

For more details see \code{vignette("SDEs")}.
}
\value{
\code{bridgesde1d} returns an object inheriting from \code{\link[base]{class}} \code{"bridgesde1d"}.
\item{X }{an invisible \code{\link[stats]{ts}} object.}
\item{drift }{drift coefficient.}
\item{diffusion }{diffusion coefficient.}
\item{C}{numbre of crossing realized.}
\item{type }{type of sde.}
\item{method }{the numerical method used.}
}
\references{
   Bladt, M. and Sorensen, M. (2007).
   Simple simulation of diffusion bridges with application to likelihood inference for diffusions.
   \emph{Working Paper, University of Copenhagen}. 

   Iacus, S.M. (2008). 
   \emph{Simulation and inference for stochastic differential equations: with R examples}. 
   Springer-Verlag, New York
}
\author{
A.C. Guidoum, K. Boukhetala.
}
\seealso{
\code{\link{bridgesde2d}} and \code{\link{bridgesde3d}} for 2 and 3-dim. 

\code{\link[sde]{DBridge}} in package \CRANpkg{sde}. 
}
\examples{
## Ito Bridge sde
## dX(t) = 2*(1-X(t)) *dt + dW(t)
## x0 = 0 at time t0=0 , and y = 1 at time T=1

f <- expression( 2*(1-x) )
g <- expression( 1 )
X <- bridgesde1d(drift=f,diffusion=g,y=1,N=1000,M=100)
X
plot(X,plot.type="single")
lines(time(X),mean(X),col=2)
}

\keyword{sde}
\keyword{ts}
\keyword{mts}
