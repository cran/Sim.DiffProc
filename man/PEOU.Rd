\name{PEOU}
\alias{PEOU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Estimation of Ornstein-Uhlenbeck Model (Exact likelihood inference)
}
\description{
Parametric estimation of Ornstein-Uhlenbeck Model.
}
\usage{
PEOU(X, delta, starts = list(r= 1, sigma= 1), leve = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of the observed time-series values.
}
  \item{delta}{
the fraction of the sampling period between successive observations.
}
  \item{starts}{
named list. Initial values for optimizer.
}
  \item{leve}{
the confidence level required.
}
}
\details{
This process solves the stochastic differential equation : \emph{\deqn{dX(t) = - r * X(t) * dt + sigma * dW(t)}} 
It is ergodic for \code{r > 0}. We have also shown its exact conditional and stationary densities. In particular, the conditional density \code{p(t,.|x)} is the density of a \code{Gaussian law} with \code{mean = x0 * exp(-r*t)} and \code{variance = ((sigma^2)/(2*r))*(1-exp(-2*r*t))}.

\code{R} has the \code{[dqpr]norm} functions to evaluate the density, the quantiles, and the cumulative distribution or generate pseudo random numbers from the normal distribution.
}
\value{
\item{coef}{Coefficients extracted from the model.}
\item{AIC}{A numeric value with the corresponding AIC.}
\item{vcov}{A matrix of the estimated covariances between the parameter estimates in the linear or non-linear predictor of the model.}
\item{confint}{A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2.}
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{PEABM}} Parametric Estimation of Arithmetic Brownian Motion, \code{\link{PEOUexp}} Explicit Estimators of Ornstein-Uhlenbeck Model, \code{\link{PEOUG}} Parametric Estimation of Hull-White/Vasicek Models, \code{\link{PEBS}} Parametric Estimation of model Black-Scholes.
}
\examples{
## Parametric estimation of Ornstein-Uhlenbeck Model.
## t0 = 0 ,T = 10
data(DATA1)
res <- PEOU(DATA1,delta=0.01,starts=list(r=2,sigma=1),leve = 0.90)
res
OUF(N=1000,M=10,t0=0,T=10,x0=40,r=0.1979284,sigma=3.972637)
points(seq(0,10,length=length(DATA1)),DATA1,type="l",lwd=3,col="blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Parametric Estimation }
