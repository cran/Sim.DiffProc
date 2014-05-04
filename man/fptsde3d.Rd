\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{fptsde3d}
\alias{fptsde3d}
\alias{fptsde3d.default}
\alias{summary.fptsde3d}
\alias{mean.fptsde3d}
\alias{median.fptsde3d}
\alias{quantile.fptsde3d}
\alias{kurtosis.fptsde3d}
\alias{skewness.fptsde3d}
\alias{moment.fptsde3d}
\alias{bconfint.fptsde3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{First Passage Time in 3-Dim SDE}
\description{The (S3) generic function \code{fptsde3d} for simulate first-passage-time (f.p.t) in 3-dim stochastic differential equations.}
\usage{
fptsde3d(N, \dots)
\method{fptsde3d}{default}(N = 100, M = 10, x0 = 0, y0 = 0, z0 = 0, t0 = 0, T = 1, Dt, c = 0.5, 
   driftx, diffx, drifty, diffy, driftz, diffz, alpha = 0.5, mu = 0.5, 
   type = c("ito", "str"), method = c("euler", "milstein", "predcorr", 
   "smilstein", "taylor", "heun", "rk1", "rk2", "rk3"), \dots)
						   
						   
\method{summary}{fptsde3d}(object, \dots)	
\method{mean}{fptsde3d}(x, \dots)
\method{median}{fptsde3d}(x, \dots)
\method{quantile}{fptsde3d}(x, \dots)
\method{kurtosis}{fptsde3d}(x, \dots)
\method{skewness}{fptsde3d}(x, \dots)
\method{moment}{fptsde3d}(x, order = 2, \dots)
\method{bconfint}{fptsde3d}(x, level=0.95, \dots)							   
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{size of sde.}
  \item{M}{size of fpt.}
  \item{x0}{initial value of the process \eqn{X_{t}}{X(t)} at time \code{t0}.}
  \item{y0}{initial value of the process \eqn{Y_{t}}{Y(t)} at time \code{t0}.}
  \item{z0}{initial value of the process \eqn{Z_{t}}{Z(t)} at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{c}{boundary or barrier (threshold).}
  \item{driftx}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}.}
  \item{diffx}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}.}
  \item{drifty}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Y_t}{Y(t)}.}
  \item{diffy}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Y_t}{Y(t)}.}
  \item{driftz}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Z_t}{Z(t)}.}
  \item{diffz}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{Z_t}{Z(t)}.}
  \item{alpha}{weight \code{alpha} of the predictor-corrector scheme; the default \code{alpha = 0.5}.}
  \item{mu}{weight \code{mu} of the predictor-corrector scheme; the default \code{mu = 0.5}.}
  \item{type}{sde of the type Ito or Stratonovich.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde3d}}.}
  \item{x, object}{an object inheriting from class \code{"fptsde3d"}.}
  \item{order}{order of moment.}
  \item{level}{the confidence level required.}    
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{fptsde3d} returns a random variables \eqn{\tau_{c}(x)}{tau_x(c)}, \eqn{\tau_{c}(y)}{tau_y(c)} and \eqn{\tau_{c}(z)}{tau_z(c)}
"first passage time" for \eqn{(X(t),Y(t),Z(t))}{(X(t),Y(t),Z(t))}, defined by : 
\deqn{\tau_{c}(x) = \{ t \geq 0 ; X_{t} \geq c \},\quad if \quad (x_{0} < c)}{tau_x(c)={t>=0; X(t) >= c}, if (x0 < c)}
\deqn{\tau_{c}(y) = \{ t \geq 0 ; Y_{t} \geq c \},\quad if \quad (y_{0} < c)}{tau_y(c)={t>=0; Y(t) >= c}, if (y0 < c)}
\deqn{\tau_{c}(z) = \{ t \geq 0 ; Z_{t} \geq c \},\quad if \quad (z_{0} < c)}{tau_z(c)={t>=0; Z(t) >= c}, if (z0 < c)}
and 
\deqn{\tau_{c}(x) = \{ t \geq 0 ; X_{t} \leq c \},\quad if \quad (x_{0} > c)}{tau_x(c)={t>=0; X(t) <= c}, if (x0 > c)}
\deqn{\tau_{c}(y) = \{ t \geq 0 ; Y_{t} \leq c \},\quad if \quad (y_{0} > c)}{tau_y(c)={t>=0; Y(t) <= c}, if (y0 > c)}
\deqn{\tau_{c}(z) = \{ t \geq 0 ; Z_{t} \leq c \},\quad if \quad (z_{0} > c)}{tau_z(c)={t>=0; Z(t) <= c}, if (z0 > c)}
with \eqn{c}{c} is a fixed boundary or barrier.
}
\value{
\code{fptsde3d} returns an object inheriting from \code{\link[base]{class}} \code{"fptsde3d"}.
\item{tau_x, tau_y, tau_z }{a vector of triplet 'fpt' \eqn{(\tau_{c}(x),\tau_{c}(y),\tau_{c}(z))}{tau_x(c),tau_y(c),tau_z(c)}.}
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
\code{\link{fptsde1d}} for simulation fpt in sde 1-dim. \code{\link[fptdApprox]{FPTL}} in package \CRANpkg{fptdApprox} for computes values of the first passage time location. 
}
\examples{

## Example :
## dX(t) = W2(t) dt + W3(t) * dW1(t) 
## dY(t) = dW2(t) and dZ(t) = dW3(t)        
## x0 = y0 = = z0 = 0, and barrier c = 0.5.
## W1(t), W2(t) and W3(t) three independent Brownian motion

fx <- expression(y)
gx <- expression(z)
fy <- expression(0)
gy <- expression(1)
fz <- expression(0)
gz <- expression(1)

res3 <- fptsde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,N=1000,M=50,c=0.5)
res3
summary(res3)
bconfint(res3,level=0.95)
moment(res3,order=c(2,3,4,5))
X <- cbind(res3$tau_x,res3$tau_y,res3$tau_z)
## library(sm)
## sm.density(X,display="rgl")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fpt}
\keyword{sde}
\keyword{ts}
\keyword{mts}
