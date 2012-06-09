\name{INFSR}
\alias{INFSR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Ahn and Gao model or Inverse of Feller Square Root Models (by Milstein Scheme) 
}
\description{
Simulation the inverse of feller square root model by milstein scheme.
}
\usage{
INFSR(N, M, t0, T, x0, theta, r, sigma, output = FALSE)
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
constant ( \code{X(t)*(theta-(sigma^3-theta*r)*X(t)) :drift coefficient}). 
}
  \item{r}{
constant ( \code{X(t)*(theta-(sigma^3-theta*r)*X(t)) :drift coefficient}).
}
  \item{sigma}{
constant positive ( \code{sigma * X(t)^(3/2) :diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
A process \code{X} satisfying : \emph{\deqn{dX(t) = X(t)*(theta-(sigma^3-theta*r)*X(t)) * dt + sigma * X(t)^(3/2) * dW(t)}}
With \code{X(t)*(theta-(sigma^3-theta*r)*X(t)) :drift coefficient} and \code{sigma * X(t)^(3/2) :diffusion coefficient}, \code{W(t)} is Wiener process, discretization \code{dt = (T-t0)/N}.
 
The conditional distribution of this process is related to that of the Cox-Ingersoll-Ross \code{(CIR)} model.
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{CEV}} Constant Elasticity of Variance Models, \code{\link{CIR}} Cox-Ingersoll-Ross Models, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{CKLS}} Chan-Karolyi-Longstaff-Sanders Models, \code{\link{DWP}} Double-Well Potential Model, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models,   \code{\link{JDP}} Jacobi Diffusion Process, \code{\link{PDP}} Pearson Diffusions Process, \code{\link{ROU}} Radial Ornstein-Uhlenbeck Process, \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## Inverse of Feller Square Root Models 
## dX(t) = X(t)*(0.5-(1^3-0.5*0.5)*X(t)) * dt + 1 * X(t)^(3/2) * dW(t)
## One trajectorie
 INFSR(N=1000,M=1,T=50,t0=0,x0=0.5,theta=0.5,r=0.5,sigma=1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}

