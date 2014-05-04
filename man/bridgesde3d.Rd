\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{bridgesde3d}
\alias{bridgesde3d}
\alias{bridgesde3d.default}
\alias{print.bridgesde3d}
\alias{time.bridgesde3d}
\alias{summary.bridgesde3d}
\alias{plot.bridgesde3d}
\alias{points.bridgesde3d}
\alias{lines.bridgesde3d}
\alias{plot3D.bridgesde3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of 3-Dim Diffusion Bridge}
\description{The (S3) generic function \code{bridgesde3d} for simulation of 3-dim diffusion bridge.}
\usage{
bridgesde3d(N, \dots)
\method{bridgesde3d}{default}(N = 1000, x0 = c(0, 0, 0), y = c(1, -1, 2), t0 = 0, T = 1, Dt, 
   driftx, diffx, drifty, diffy, driftz, diffz, alpha = 0.5, mu = 0.5, 
   type = c("ito", "str"), method = c("euler", "milstein", "predcorr", 
   "smilstein", "taylor", "heun", "rk1", "rk2","rk3"), \dots)
							  
\method{summary}{bridgesde3d}(object, \dots)							  
\method{time}{bridgesde3d}(x, \dots)
\method{plot3D}{bridgesde3d}(x, display = c("persp","rgl"), \dots)
\method{plot}{bridgesde3d}(x, \dots)
\method{lines}{bridgesde3d}(x, \dots)
\method{points}{bridgesde3d}(x, \dots)								  
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{number of simulation steps.}
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
  \item{x, object}{an object inheriting from class \code{"bridgesde3d"}.}
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
\item{XYZ }{an invisible \code{\link[stats]{ts}} object \eqn{(X_{t},Y_{t},Z_{t})}{(X(t),Y(t),Z(t))}.}
}
\references{
   Bladt, M. and Sorensen, M. (2007).
   Simple simulation of diffusion bridges with application to likelihood inference for diffusions.
   \emph{Working Paper, University of Copenhagen}. 
   Available at \url{http://www.math.ku.dk/~michael/diffusionbridgepreprint.pdf}

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
## dX(t) = Y(t)*Z(t) * dt +  dW1(t)          
## dY(t) = X(t)*Z(t) * dt +  dW2(t)
## dZ(t) = X(t)*Y(t) * dt +  dW3(t)
## x01 = 1 , y01 = 1
## x02 = 1 , y02 = 1 
## x03 = 1 , y03 = 1
## W1(t), W2(t) and W3(t) three independent Brownian motion

fx <- expression(y*z)
gx <- expression(1)
fy <- expression(x*z)
gy <- expression(1)
fz <- expression(x*y)
gz <- expression(1)

XYZ <- bridgesde3d(x0=c(1,1,1),y=c(1,1,1),driftx=fx,diffx=gx,drifty=fy,
                   diffy=gy,driftz=fz,diffz=gz)
XYZ
summary(XYZ)
plot(XYZ)
dev.new()
plot(XYZ,plot.type="single")
dev.new()
plot3D(XYZ,display = "persp",main="3-dim bridge sde")
plot3D(XYZ,display = "rgl")
}
\keyword{sde}
\keyword{ts}
\keyword{mts}
