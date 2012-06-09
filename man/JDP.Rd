\name{JDP}
\alias{JDP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating The Jacobi Diffusion Process (by Milstein Scheme)
}
\description{
Simulation the jacobi diffusion process by milstein scheme.
}
\usage{
JDP(N, M, t0, T, x0, theta, output = FALSE)
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
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
The Jacobi diffusion process is the solution to the stochastic differential equation : \emph{\deqn{dX(t) = -theta * (X(t) - 0.5)* dt + sqrt( theta*X(t)*(1-X(t))) * dW(t)}}
With \code{-theta * (X(t) - 0.5) :drift coefficient} and \code{sqrt( theta*X(t)*(1-X(t))) :diffusion coefficient}, \code{W(t)} is Wiener process, discretization \code{dt = (T-t0)/N}.
 
For \code{theta > 0}. It has an invariant distribution that is \code{uniform on [0,1]}.
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{CEV}} Constant Elasticity of Variance Models, \code{\link{CIR}} Cox-Ingersoll-Ross Models, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{CKLS}} Chan-Karolyi-Longstaff-Sanders Models, \code{\link{DWP}} Double-Well Potential Model, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{INFSR}} Inverse of Feller s Square Root models,  \code{\link{PDP}} Pearson Diffusions Process, \code{\link{ROU}} Radial Ornstein-Uhlenbeck Process, \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## Jacobi Diffusion Process 
## dX(t) = -0.05 * (X(t)-0.5)* dt + sqrt(0.05*X(t)*(1-X(t))) * dW(t),  
## One trajectorie
 JDP(N=1000,M=1,T=100,t0=0,x0=0,theta=0.05)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}