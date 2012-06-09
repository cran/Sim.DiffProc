\name{CIRhy}
\alias{CIRhy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating The modified CIR and hyperbolic Process (by Milstein Scheme)}
\description{
Simulation the modified CIR and hyperbolic process by milstein scheme.
}
\usage{
CIRhy(N, M, t0, T, x0, r, sigma, output = FALSE)
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
  \item{r}{
constant ( \code{-r * X(t) :drift coefficient}). 
}
  \item{sigma}{
constant positive ( \code{sigma * sqrt(1+X(t)^2) :diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
The stochastic differential equation of the modified CIR is  : \emph{\deqn{dX(t) = -r *X(t)*dt + sigma *sqrt(1 + X(t)^2) *dW(t)}} 
With \code{-r*X(t) :drift coefficient} and \code{sigma*sqrt(1+X(t)^2) :diffusion coefficient}, \code{W(t)} is Wiener process, the discretization \code{dt = (T-t0)/N}.
 
Constraints: \code{r + (sigma^2)/2 > 0} (this is needed to make the process positive recurrent).
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{CEV}} Constant Elasticity of Variance Models, \code{\link{CIR}} Cox-Ingersoll-Ross Models , \code{\link{CKLS}} Chan-Karolyi-Longstaff-Sanders Models, \code{\link{DWP}} Double-Well Potential Model, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{INFSR}} Inverse of Feller s Square Root models, \code{\link{JDP}} Jacobi Diffusion Process, \code{\link{PDP}} Pearson Diffusions Process, \code{\link{ROU}} Radial Ornstein-Uhlenbeck Process, \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## The modified CIR and hyperbolic Process  
## dX(t) = - 0.3 *X(t) *dt + 0.9 * sqrt(1+X(t)^2) * dW(t)
## One trajectorie
 CIRhy(N=1000,M=1,T=1,t0=0,x0=1,r=0.3,sigma=0.9)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}
