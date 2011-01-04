\name{CIR}
\alias{CIR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Cox-Ingersoll-Ross (CIR) Square Root Diffusion Models (by Milstein Scheme)
}
\description{
Simulation cox-ingersoll-ross models by milstein scheme.
}
\usage{
CIR(N, M, t0, T, x0, theta, r, sigma, output = FALSE)
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
  \item{theta}{
constant positive ( \code{(r - theta * X(t)) :drift coefficient}). 
}
  \item{r}{
constant positive ( \code{(r - theta * X(t)) :drift coefficient}). 
}
  \item{sigma}{
constant positive ( \code{sigma * sqrt(X(t)) :diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
Another interesting family of parametric models is that of the Cox-Ingersoll-Ross process. This model was introduced by Feller as a model for population growth and became quite popular in finance after Cox, Ingersoll, and Ross proposed it to model short-term interest rates. It was recently adopted to model nitrous oxide emission from soil by Pedersen and to model the evolutionary rate variation across sites in molecular evolution.

 
The discretization \code{dt = (T-t0)/N}, and the stochastic differential equation of CIR is  : \emph{\deqn{dX(t) = (r - theta *X(t))*dt + sigma *sqrt(X(t)) *dW(t)}}
With \code{(r - theta *X(t)) :drift coefficient} and \code{sigma*sqrt(X(t)) :diffusion coefficient}, \code{W(t)} is Wiener process.

Constraints: \code{2*r > sigma^2}. 
}
\value{
data.frame(time,x) and plot of process.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{CEV}} Constant Elasticity of Variance Models, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{CKLS}} Chan-Karolyi-Longstaff-Sanders Models, \code{\link{DWP}} Double-Well Potential Model, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{INFSR}} Inverse of Feller s Square Root models, \code{\link{JDP}} Jacobi Diffusion Process, \code{\link{PDP}} Pearson Diffusions Process, \code{\link{ROU}} Radial Ornstein-Uhlenbeck Process, \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
## Cox-Ingersoll-Ross Models 
## dX(t) = (0.1 - 0.2 *X(t)) *dt + 0.05 * sqrt(X(t)) * dW(t)
## One trajectorie
CIR(N=1000,M=1,t0=0,T=1,x0=0.2,theta=0.2,r=0.1,sigma=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}
