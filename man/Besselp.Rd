\name{Besselp}
\alias{Besselp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Bessel process (by Milstein Scheme)
}
\description{
Simulation Besselp process by milstein scheme.
}
\usage{
Besselp(N, M, t0, T, x0, alpha, output = FALSE)
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
  \item{alpha}{
constant positive \code{alpha >=2}.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
The stochastic differential equation of \code{Bessel process} is  : \emph{\deqn{dX(t) = (alpha-1)/(2*X(t)) * dt + dW(t)}} with \code{(alpha-1)/(2*X(t)) :drift coefficient} and \code{1 :diffusion coefficient}, \code{W(t)} is Wiener process, and the discretization \code{dt = (T-t0)/N}. 

Constraints: \code{alpha >= 2} \code{and x0 =! 0}.
}
\value{
data.frame(time,x) and plot of process. 
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{CEV}} Constant Elasticity of Variance Models, \code{\link{CIR}} Cox-Ingersoll-Ross Models, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{CKLS}} Chan-Karolyi-Longstaff-Sanders Models, \code{\link{DWP}} Double-Well Potential Model, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{INFSR}} Inverse of Feller s Square Root models, \code{\link{JDP}} Jacobi Diffusion Process, \code{\link{PDP}} Pearson Diffusions Process, \code{\link{ROU}} Radial Ornstein-Uhlenbeck Process, \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
## Bessel Process
## alpha = 4
## dX(t) = 3/(2*x) * dt + dW(t)
## One trajectorie
Besselp(N=1000,M=1,t0=0,T=100,x0=1,alpha=4,output=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
