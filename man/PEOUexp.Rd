\name{PEOUexp}
\alias{PEOUexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Estimation of Ornstein-Uhlenbeck Model (Explicit Estimators)
}
\description{
Explicit estimators of Ornstein-Uhlenbeck Model.
}
\usage{
PEOUexp(X, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of the observed time-series values.
}
  \item{delta}{
the fraction of the sampling period between successive observations.
}
}
\details{
This process solves the stochastic differential equation : \emph{\deqn{dX(t) = - r * X(t) * dt + sigma * dW(t)}}
It is ergodic for \code{r > 0}. 

We have also shown its exact conditional and stationary densities. In particular, the conditional density \code{p(t,.|x)} is the density of a \code{Gaussian law} with \code{mean = x0 * exp(-r*t)} and \code{variance = ((sigma^2)/(2*r))*(1-exp(-2*r*t))}, the maximum likelihood estimator of \code{r} is available in explicit form and takes the form : \emph{\deqn{r = -(1/dt) * log(sum(X(t)*X(t-1)) / sum( X(t-1)^2))}}
which is defined only if  \code{sum(X(t)*X(t-1)) > 0}, this estimator is consistent and asymptotically Gaussian.

The maximum likelihood estimator of : \emph{\deqn{sigma^2 = (2*r)/(N*(1-exp(-2*dt*r))) * sum(X(t)-X(t-1)**exp(-dt*r))^2}}
}
\value{
\item{r}{Estimator of speed of reversion.}
\item{sigma}{Estimator of volatility.}
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{PEABM}} Parametric Estimation of Arithmetic Brownian Motion, \code{\link{PEOU}} Parametric Estimation of Ornstein-Uhlenbeck Model, \code{\link{PEOUG}} Parametric Estimation of Hull-White/Vasicek Models, \code{\link{PEBS}} Parametric Estimation of model Black-Scholes.
}
\examples{
## t0 = 0 ,T = 10
data(DATA1)
res <- PEOUexp(DATA1,delt=0.01)
res
OUF(N=1000,M=10,t0=0,T=10,x0=DATA1[1],r=res$r,sigma=res$sigma)
points(seq(0,10,length=length(DATA1)),DATA1,type="l",lwd=3,col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Parametric Estimation }