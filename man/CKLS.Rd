\name{CKLS}
\alias{CKLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating The Chan-Karolyi-Longstaff-Sanders (CKLS) family of models (by Milstein Scheme)
}
\description{
Simulation the chan-karolyi-longstaff-sanders models by milstein scheme.
}
\usage{
CKLS(N, M, t0, T, x0, r, theta, sigma, gamma, output = FALSE)
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
constant ( \code{(r + theta *X(t)) :drift coefficient}). 
}
  \item{theta}{
constant ( \code{(r + theta *X(t)) :drift coefficient}). 
}
  \item{sigma}{
constant positive ( \code{sigma * X(t)^gamma :diffusion coefficient}).
}
  \item{gamma}{
constant positive ( \code{sigma * X(t)^gamma :diffusion coefficient}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
The Chan-Karolyi-Longstaff-Sanders \code{(CKLS)} family of models is a class of parametric stochastic differential equations widely used in many finance applications, in particular to model interest rates or asset prices. 

The CKLS process solves the stochastic differential equation : \emph{\deqn{dX(t) = (r + theta * X(t))*dt + sigma *X(t)^gamma * dW(t)}}
With \code{(r + theta * X(t)) :drift coefficient} and \code{sigma* X(t)^gamma :diffusion coefficient}, \code{W(t)} is Wiener process, the discretization \code{dt = (T-t0)/N}.
 
This \code{CKLS} model is a further extension of the Cox-Ingersoll-Ross model and hence embeds all previous models. 

The \code{CKLS} model does not admit an explicit transition density unless \code{r = 0} or \code{gamma = 0.5}. It takes values in \code{(0, + lnf)} if \code{r},\code{theta > 0}, and \code{gamma > 0.5}. In all cases, \code{sigma} is assumed to be positive.   
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{CEV}} Constant Elasticity of Variance Models, \code{\link{CIR}} Cox-Ingersoll-Ross Models, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{DWP}} Double-Well Potential Model, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{INFSR}} Inverse of Feller s Square Root models, \code{\link{JDP}} Jacobi Diffusion Process, \code{\link{PDP}} Pearson Diffusions Process, \code{\link{ROU}} Radial Ornstein-Uhlenbeck Process, \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## Chan-Karolyi-Longstaff-Sanders Models 
## dX(t) = (0.3 + 0.01 *X(t)) *dt + 0.1 * X(t)^0.2 * dW(t)
## One trajectorie
 CKLS(N=1000,M=1,T=1,t0=0,x0=1,r=0.3,theta=0.01,sigma=0.1,gamma=   0.2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}
