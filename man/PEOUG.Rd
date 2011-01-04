\name{PEOUG}
\alias{PEOUG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Estimation of Hull-White/Vasicek (HWV) Gaussian Diffusion Models(Exact likelihood inference)
}
\description{
Parametric estimation of Hull-White/Vasicek Model.
}
\usage{
PEOUG(X, delta, starts = list(r= 1, theta= 1, sigma= 1), leve = 0.95)
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
the Vasicek or Ornstein-Uhlenbeck model solves the stochastic differential equation : \emph{\deqn{dX(t) = r * (theta -  X(t)) * dt + sigma * dW(t)}}
It is ergodic for \code{r > 0}. We have also shown its exact conditional and stationary densities. In particular, the conditional density \code{p(t,.|x)} is the density of a \code{Gaussian law} with \code{mean = theta+(x0-theta)*exp(-r*t)} and \code{variance = (sigma^2/(2*r))*(1-exp(-2*r*t))}. 

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
\code{\link{PEABM}} Parametric Estimation of Arithmetic Brownian Motion, \code{\link{PEOUexp}} Explicit Estimators of Ornstein-Uhlenbeck Model, \code{\link{PEOU}} Parametric Estimation of Ornstein-Uhlenbeck Model, \code{\link{PEBS}} Parametric Estimation of model Black-Scholes.
}
\examples{
## example 1
## t0 = 0 ,T = 10
data(DATA1)
res <- PEOUG(DATA1,delta=0.01,starts=list(r=2,theta=0,sigma=1))
res
HWVF(N=1000,M=10,t0=0,T=10,x0=40,r=0.9979465,theta=16.49602,sigma=3.963486)
points(seq(0,10,length=length(DATA1)),DATA1,type="l",lwd=3,col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Parametric Estimation }
