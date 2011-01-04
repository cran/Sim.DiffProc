\name{ROU}
\alias{ROU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Radial Ornstein-Uhlenbeck Process (by Milstein Scheme)
}
\description{
Simulation the radial ornstein-uhlenbeck process by milstein scheme.
}
\usage{
ROU(N, M, t0, T, x0, theta, output = FALSE)
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
constant positive.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
The radial Ornstein-Uhlenbeck process is the solution to the stochastic differential equation :\emph{\deqn{dX(t) = (theta * X(t)^-(1)  - X(t))  * dt + dW(t)}}
With \code{(theta * X(t)^-1  - X(t)) :drift coefficient} and \code{1 :diffusion coefficient}, the discretization \code{dt = (T-t0)/N}, \code{W(t)} is Wiener process.
}
\value{
data.frame(time,x) and plot of process. 
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{CEV}} Constant Elasticity of Variance Models, \code{\link{CIR}} Cox-Ingersoll-Ross Models, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{CKLS}} Chan-Karolyi-Longstaff-Sanders Models, \code{\link{DWP}} Double-Well Potential Model, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{INFSR}} Inverse of Feller s Square Root models, \code{\link{JDP}} Jacobi Diffusion Process, \code{\link{PDP}} Pearson Diffusions Process,  \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
## Radial Ornstein-Uhlenbeck 
## dX(t) = (0.05*X(t)^(-1) - X(t)) *dt + dW(t)
## One trajectorie
ROU(N=1000,M=1,T=1,t0=0,x0=1,theta=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}
