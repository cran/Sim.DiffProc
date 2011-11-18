\name{PEBS}
\alias{PEBS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Estimation of Model Black-Scholes (Exact likelihood inference)
}
\description{
Parametric estimation of model Black-Scholes.
}
\usage{
PEBS(X, delta, starts = list(theta= 1, sigma= 1), leve = 0.95)
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
The Black and Scholes, or geometric Brownian motion model solves the stochastic differential equation: \emph{\deqn{dX(t) = theta * X(t)* dt + sigma * X(t)* dW(t)}} 
The conditional density function \code{p(t,.|x)} is \code{log-normal} with \code{mean = x * exp(theta*t)} and \code{variance =  x^2 * exp(2*theta*t)*(exp(sigma^2 *t ) -1 )}.

\code{R} has the \code{[dqpr]lnorm} functions to evaluate the density, the quantiles, and the cumulative distribution or generate pseudo random numbers from the lognormal distribution.
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
\code{\link{PEABM}} Parametric Estimation of Arithmetic Brownian Motion, \code{\link{PEOU}} Parametric Estimation of Ornstein-Uhlenbeck Model, \code{\link{PEOUexp}} Explicit Estimators of Ornstein-Uhlenbeck Model, \code{\link{PEOUG}} Parametric Estimation of Hull-White/Vasicek Models.
}
\examples{
## Parametric estimation of model Black-Scholes.
## t0 = 0 ,T = 1
data(DATA2)
res <- PEBS(DATA2,delta=0.001,starts=list(theta=2,sigma=1))
res
GBMF(N=1000,M=10,T=1,t0=0,x0=DATA2[1],theta=res$coef[1],sigma=res$coef[2])
points(seq(0,1,length=length(DATA2)),DATA2,type="l",lwd=3,col="blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Parametric Estimation }
