\name{PEABM}
\alias{PEABM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Estimation of Arithmetic Brownian Motion(Exact likelihood inference)
}
\description{
Parametric estimation of Arithmetic Brownian Motion
}
\usage{
PEABM(X, delta, starts = list(theta= 1, sigma= 1), leve = 0.95)
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
This process solves the stochastic differential equation : \emph{\deqn{dX(t) = theta * dt + sigma * dW(t)}} 
The conditional density \code{p(t,.|x)} is the density of a \code{Gaussian law} with \code{mean = x0 + theta * t} and \code{variance = sigma^2 * t}.

\code{R} has the \code{[dqpr]norm} functions to evaluate the density, the quantiles, and the cumulative distribution or generate pseudo random numbers from the normal distribution.
}
\value{
\item{coef}{Coefficients extracted from the model.}
\item{AIC}{A numeric value with the corresponding AIC.}
\item{vcov}{A matrix of the estimated covariances between the parameter estimates in the linear or non-linear predictor of the model.}
\item{confint}{A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2.}
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{PEOU}} Parametric Estimation of Ornstein-Uhlenbeck Model, \code{\link{PEOUexp}} Explicit Estimators of Ornstein-Uhlenbeck Model, \code{\link{PEOUG}} Parametric Estimation of Hull-White/Vasicek Models, \code{\link{PEBS}} Parametric Estimation of model Black-Scholes.
}
\examples{
\donttest{
## Parametric estimation of Arithmetic Brownian Motion.
## t0 = 0 ,T = 100
 data(DATA3)
 res <- PEABM(DATA3,delta=0.1,starts=list(theta=1,sigma=1),leve = 0.95)
 res
 ABMF(N=1000,M=10,t0=0,T=100,x0=DATA3[1],theta=res$coef[1],sigma=res$coef[2])
 points(seq(0,100,length=length(DATA3)),DATA3,type="l",lwd=3,col="blue")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Parametric Estimation }
