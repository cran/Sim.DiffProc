\name{DWP}
\alias{DWP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Double-Well Potential Model (by Milstein Scheme)
}
\description{
Simulation double-well potential model by milstein scheme.
}
\usage{
DWP(N, M, t0, T, x0, output = FALSE)
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
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
This model is interesting because of the fact that its density has a bimodal shape. 

The process satisfies the stochastic differential equation : \emph{\deqn{dX(t) = ( X(t) - X(t)^3) * dt + dW(t)}} With \code{(X(t) - X(t)^3) :drift coefficient} and \code{1 is diffusion coefficient}, \code{W(t)} is Wiener process,and the discretization \code{dt = (T-t0)/N}.
 
This model is challenging in the sense that the Milstein approximation.  
}
\value{
data.frame(time,x) and plot of process.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{CEV}} Constant Elasticity of Variance Models, \code{\link{CIR}} Cox-Ingersoll-Ross Models, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{CKLS}} Chan-Karolyi-Longstaff-Sanders Models, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{INFSR}} Inverse of Feller s Square Root models, \code{\link{JDP}} Jacobi Diffusion Process, \code{\link{PDP}} Pearson Diffusions Process, \code{\link{ROU}} Radial Ornstein-Uhlenbeck Process, \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
## Double-Well Potential Model 
## dX(t) = (X(t) - X(t)^3) * dt + dW(t)
## One trajectorie
DWP(N=1000,M=1,T=1,t0=0,x0=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}
