\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{bridgesde3d}
\alias{bridgesde3d}
\alias{bridgesde3d.default}
\alias{print.bridgesde3d}
\alias{time.bridgesde3d}
\alias{mean.bridgesde3d}
\alias{median.bridgesde3d}
\alias{quantile.bridgesde3d}
\alias{kurtosis.bridgesde3d}
\alias{skewness.bridgesde3d}
\alias{moment.bridgesde3d}
\alias{bconfint.bridgesde3d}
\alias{plot.bridgesde3d}
\alias{points.bridgesde3d}
\alias{lines.bridgesde3d}
\alias{plot3D.bridgesde3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 3-Dim Diffusion Bridge}
\description{The (S3) generic function \code{bridgesde3d} for simulation of 3-dim diffusion bridge.}
\usage{
bridgesde3d(N, \dots)
\method{bridgesde3d}{default}(N=1000,M=1, x0=c(0,0, 0), y=c(1,-1, 2), 
   t0 = 0, T = 1, Dt, driftx, diffx, drifty, diffy, driftz, diffz, 
   alpha = 0.5, mu = 0.5, type = c("ito", "str"), method = c("euler", 
   "milstein","predcorr", "smilstein", "taylor", "heun", "rk1", "rk2", 
   "rk3"), \dots)

\method{time}{bridgesde3d}(x, \dots)
\method{mean}{bridgesde3d}(x, \dots)
\method{median}{bridgesde3d}(x, \dots)
\method{quantile}{bridgesde3d}(x, \dots)
\method{kurtosis}{bridgesde3d}(x, \dots)
\method{skewness}{bridgesde3d}(x, \dots)
\method{moment}{bridgesde3d}(x, order = 2, \dots)
\method{bconfint}{bridgesde3d}(x, level=0.95, \dots)
\method{plot}{bridgesde3d}(x, \dots)
\method{lines}{bridgesde3d}(x, \dots)
\method{points}{bridgesde3d}(x, \dots)	
\method{plot3D}{bridgesde3d}(x, display = c("persp","rgl"), \dots)						  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories.}
  \item{x0}{initial value (numeric vector of length 3) of the process \eqn{X_t}{X(t)}, \eqn{Y_t}{Y(t)} and \eqn{Z_t}{Z(t)} at time \eqn{t_0}{t0}.}
  \item{y}{terminal value (numeric vector of length 3) of the process \eqn{X_t}{X(t)}, \eqn{Y_t}{Y(t)} and \eqn{Z_t}{Z(t)} at time \eqn{T}{T}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{driftx}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}.}
  \item{diffx}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}.}
  \item{drifty}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Y_t}{Y(t)}.}
  \item{diffy}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Y_t}{Y(t)}.}
  \item{driftz}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Z_t}{Z(t)}.}
  \item{diffz}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Z_t}{Z(t)}.}
  \item{alpha}{weight \code{alpha} of the predictor-corrector scheme; the default \code{alpha = 0.5}.}
  \item{mu}{weight \code{mu} of the predictor-corrector scheme; the default \code{mu = 0.5}.}
  \item{type}{if \code{type="ito"} simulation diffusion bridge of Ito type, else \code{type="str"} simulation diffusion bridge of Stratonovich type; the default \code{type="ito"}.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde3d}}.}
  \item{x}{an object inheriting from class \code{"bridgesde3d"}.}
  \item{order}{order of moment.}
  \item{level}{the confidence level required.} 
  \item{display}{\code{"persp"} perspective or \code{"rgl"} plots.} 
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{bridgesde3d} returns a \code{mts} of the diffusion bridge starting at \code{x} at time \code{t0} and ending at \code{y} at time \code{T}.

The methods of approximation are classified according to their different properties. Mainly two criteria of optimality are used in the literature: the strong 
and the weak (orders of) convergence. The \code{method} of simulation can be one among: \code{Euler-Maruyama Order 0.5}, \code{Milstein Order 1}, \code{Milstein Second-Order}, 
\code{Predictor-Corrector method}, \code{Ito-Taylor Order 1.5}, \code{Heun Order 2} and \code{Runge-Kutta Order 1, 2 and 3}. 

For more details see \code{vignette("SDEs")}.
}
\value{
\code{bridgesde3d} returns an object inheriting from \code{\link[base]{class}} \code{"bridgesde3d"}.
\item{X, Y, Z }{an invisible \code{mts} (3-dim) object (X(t),Y(t),Z(t)).}
\item{driftx, drifty, driftz }{drift coefficient of X(t), Y(t) and Z(t).}
\item{diffx, diffy, diffz }{diffusion coefficient of X(t), Y(t) and Z(t).}
\item{Cx, Cy, Cz}{numbre of crossing realized of X(t) (Y(t))(Z(t)).}
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
\code{\link{bridgesde1d}} for simulation of 1-dim diffusion bridge. \code{\link[sde]{DBridge}} in package \CRANpkg{sde}.
}
\examples{
## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 * dW1(t) ; x01 = 0 and y01 = 0
## dY(t) = 4*(1-Y(t)) *X(t) dt + 0.2 * dW2(t) ; x02 = -1 and y02 = -2
## dZ(t) = 4*(1-Z(t)) *Y(t) dt + 0.2 * dW3(t) ; x03 = 0.5 and y03 = 0.5       
## W1(t), W2(t) and W3(t) three independent Brownian motion
set.seed(1234)

fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)
fz <- expression(4*(1-z)*y)
gz <- expression(0.2)

res <- bridgesde3d(x0=c(0,-1,0.5),y=c(0,-2,0.5),driftx=fx,diffx=gx,
                drifty=fy,diffy=gy,driftz=fz,diffz=gz,M=20)
res
plot(res,union=TRUE)
dev.new()
plot3D(res,display = "persp",main="3-dim bridge sde")
}
\keyword{sde}
\keyword{ts}
\keyword{mts}
