\name{CEV}
\alias{CEV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Constant Elasticity of Variance (CEV) Models (by Milstein Scheme)
}
\description{
Simulation constant elasticity of variance models by milstein scheme.
}
\usage{
CEV(N, M, t0, T, x0, mu, sigma, gamma, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{M}{
number of trajectories.
}
  \item{t0}{
initial time.
}
  \item{T}{
final time.
}
  \item{x0}{
initial value of the process at time \code{t0}.
}
  \item{mu}{
constant (\code{mu * X(t) :drift coefficient}). 
}
  \item{sigma}{
constant positive (\code{sigma * X(t)^gamma :diffusion coefficient}).
}
  \item{gamma}{
constant positive (\code{sigma * X(t)^gamma :diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
The Constant Elasticity of Variance \code{(CEV)} model also derives directly from the linear drift class, the discretization \code{dt = (T-t0)/N}. 

The stochastic differential equation of \code{CEV} is  : \emph{\deqn{dX(t) = mu * X(t)* dt + sigma * X(t)^gamma *dW(t)}} with \code{mu * X(t) :drift coefficient} and \code{sigma * X(t)^gamma :diffusion coefficient}, \code{W(t)} is Wiener process.

This process is quite useful in modeling a skewed implied volatility. In particular,for \code{gamma < 1}, the skewness is negative, and for \code{gamma > 1} the skewness is positive. For \code{gamma = 1}, the \code{CEV} process is a particular version of the \code{geometric Brownian motion}.  
}
\value{
data.frame(time,x) and plot of process. 
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{CIR}} Cox-Ingersoll-Ross Models, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{CKLS}} Chan-Karolyi-Longstaff-Sanders Models, \code{\link{DWP}} Double-Well Potential Model, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{INFSR}} Inverse of Feller s Square Root models, \code{\link{JDP}} Jacobi Diffusion Process, \code{\link{PDP}} Pearson Diffusions Process, \code{\link{ROU}} Radial Ornstein-Uhlenbeck Process, \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
## Constant Elasticity of Variance Models
## dX(t) = 0.3 *X(t) *dt + 2 * X(t)^1.2 * dW(t)
## One trajectorie
CEV(N=1000,M=1,t0=0,T=1,x0=0.1,mu=0.3,sigma=2,gamma=1.2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}
