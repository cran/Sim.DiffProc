

\name{rsde1d}
\alias{rsde1d}
\alias{rsde1d.default}
\alias{summary.rsde1d}
\alias{mean.rsde1d}
\alias{median.rsde1d}
\alias{quantile.rsde1d}
\alias{kurtosis.rsde1d}
\alias{skewness.rsde1d}
\alias{moment.rsde1d}
\alias{bconfint.rsde1d}
\alias{plot.rsde1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random Number Generators for 1-Dim SDE}
\description{The (S3) generic function \code{rsde1d} for simulate random number generators to generate 1-dim sde.}
\usage{
rsde1d(N, \dots)
\method{rsde1d}{default}(N = 1000, M = 100, x0 = 0, t0 = 0, T = 1, Dt, tau = 0.5, 
   drift, diffusion, alpha = 0.5, mu = 0.5,type = c("ito", "str"), 
   method = c("euler", "milstein", "predcorr", "smilstein", "taylor",
   "heun", "rk1", "rk2", "rk3"), \dots)
   
\method{summary}{rsde1d}(object, \dots)
\method{mean}{rsde1d}(x, \dots)
\method{median}{rsde1d}(x, \dots)
\method{quantile}{rsde1d}(x, \dots)
\method{kurtosis}{rsde1d}(x, \dots)
\method{skewness}{rsde1d}(x, \dots)
\method{moment}{rsde1d}(x, order = 2, \dots)
\method{bconfint}{rsde1d}(x, level=0.95, \dots)
\method{plot}{rsde1d}(x, \dots)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{size of sde.}
  \item{M}{number of random numbers to be geneated.}
  \item{x0}{initial value of the process at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{tau}{time between \code{t0} and \code{T}. Random number generated at \code{time=tau}.}
  \item{drift}{drift coefficient: an \code{\link[base]{expression}} of two variables \code{t} and \code{x}.}
  \item{diffusion}{diffusion coefficient: an \code{\link[base]{expression}} of two variables \code{t} and \code{x}.}
  \item{alpha, mu}{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{sde of  the type Ito or Stratonovich.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde1d}}.}
  \item{x, object}{an object inheriting from class \code{"rsde1d"}.}
  \item{order}{order of moment.}
  \item{level}{the confidence level required.}  
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{rsde1d} returns a random variable \eqn{x_{\tau}}{x(tau)} realize at time \eqn{t=\tau}{t=tau} defined by : 
\deqn{x_{ \tau } = \{ t \geq 0 ; x = X_{ \tau } \}}{x(tau) = { t>=0 ; x = X(tau) }}
with \eqn{ \tau }{tau} is a fixed time between \code{t0} and \code{T}.
}
\value{
\code{rsde1d} returns an object inheriting from \code{\link[base]{class}} \code{"rsde1d"}.
\item{x }{a vector of random numbers of 1-dim sde realize at time \eqn{t=\tau}{t=tau}.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link{rsde2d}} and \code{\link{rsde3d}} simulation RNs in sde 2 and 3-dim.

\code{\link[yuima]{rng}} random number generators in "yuima" package. 

\code{\link[sde]{rcBS}}, \code{\link[sde]{rcCIR}}, \code{\link[sde]{rcOU}} and \code{\link[sde]{rsOU}}  in package "sde". 
}
\examples{

## Example 1: Ito sde
## dX(t) = 2*(3-X(t)) *dt + dW(t)
set.seed(1234)

f <- expression( 4*(2-x) )
g <- expression( 0.2 )
res <- rsde1d(drift=f,diffusion=g,tau=1.75,T=2)
res
summary(res)
plot(res,pos=7,cex=1)
dev.new()
plot(density(res$x))

## Example 2: Stratonovich sde
## dX(t) = (-2*(X(t)<=0)+2*(X(t)>=0)) *dt + 0.5 o dW(t)
set.seed(1234)

f <- expression(-2*(x<=0)+2*(x>=0))
g <- expression(0.5)
res1 <- rsde1d(drift=f,diffusion=g,tau=0.95123,type="str")
res1
summary(res1)
plot(res1,pos=3,cex=1)
dev.new()
plot(density(res1$x))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
\keyword{random generators}
