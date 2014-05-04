\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{st.int}
\alias{st.int}
\alias{st.int.default}
\alias{summary.st.int}
\alias{print.st.int}
\alias{time.st.int}
\alias{mean.st.int}
\alias{median.st.int}
\alias{quantile.st.int}
\alias{kurtosis.st.int}
\alias{skewness.st.int}
\alias{moment.st.int}
\alias{bconfint.st.int}
\alias{plot.st.int}
\alias{points.st.int}
\alias{lines.st.int}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stochastic Integrals}
\description{
The (S3) generic function \code{st.int} of simulation of stochastic integrals of Ito or Stratonovich type.
}
\usage{
st.int(expr, \dots)
\method{st.int}{default}(expr, lower = 0, upper = 1, M = 1, subdivisions = 1000L, 
               type = c("ito", "str"), \dots)

\method{summary}{st.int}(object, \dots)
\method{time}{st.int}(x, \dots)
\method{mean}{st.int}(x, \dots)
\method{median}{st.int}(x, \dots)
\method{quantile}{st.int}(x, \dots)
\method{kurtosis}{st.int}(x, \dots)
\method{skewness}{st.int}(x, \dots)
\method{moment}{st.int}(x, order = 2, \dots)
\method{bconfint}{st.int}(x, level=0.95, \dots)
\method{plot}{st.int}(x, \dots)
\method{lines}{st.int}(x, \dots)
\method{points}{st.int}(x, \dots)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{expr}{an \code{\link[base]{expression}} of two variables \code{t} (time) and \code{w} (\code{w}: standard Brownian motion).}
  \item{lower, upper}{the lower and upper end points of the interval to be integrate.}
  \item{M}{number of trajectories.}
  \item{subdivisions}{the maximum number of subintervals.}
  \item{type}{Ito or Stratonovich integration.}
  \item{x, object}{an object inheriting from class \code{"st.int"}.}
  \item{order}{order of moment.}
  \item{level}{the confidence level required.}    
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{st.int} returns a \code{\link[stats]{ts}} x of length \code{N+1}; i.e. simulation of stochastic integrals 
of Ito or Stratonovich type.

The Ito interpretation is:
\deqn{\int_{t_{0}}^{t} f(s) dW_{s} = \lim_{N \rightarrow \infty} \sum_{i=1}^{N} f(t_{i-1})(W_{t_{i}}-W_{t_{i-1}})}{int(f(s)*dw(s),t0,T) = sum(f(t(i-1)) * (W(t(i)) - W(t(i-1))),i=1,\dots,N)}
The Stratonovich interpretation is:
\deqn{\int_{t_{0}}^{t} f(s) \circ dW_{s} = \lim_{N \rightarrow \infty} \sum_{i=1}^{N} f\left(\frac{t_{i}+t_{i-1}}{2}\right)(W_{t_{i}}-W_{t_{i-1}})}{int(f(s) o dw(s),t0,T) = sum(f((t(i)+t(i-1))/2) * (W(t(i)) - W(t(i-1))),i=1,\dots,N)}

For more details see \code{vignette("SDEs")}.
}
\value{
\code{st.int} returns an object inheriting from \code{\link[base]{class}} \code{"st.int"}.
\item{X }{the final simulation of the integral, an invisible \code{\link[stats]{ts}} object.}
\item{fun }{function to be integrated.}
\item{type }{type of stochastic integral.}
\item{subdivisions }{the number of subintervals produced in the subdivision process.}
}
\references{
  Ito, K. (1944). 
  Stochastic integral. 
  \emph{Proc. Jap. Acad, Tokyo}, \bold{20}, 19--529.
  
  Kloeden, P.E, and Platen, E. (1995).
  \emph{Numerical Solution of Stochastic Differential Equations}. 
  Springer-Verlag, New York.

  Oksendal, B. (2000).
  \emph{Stochastic Differential Equations: An Introduction with Applications}.
  5th edn. Springer-Verlag, Berlin.
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link{snssde1d}}, \code{\link{snssde2d}} and \code{\link{snssde3d}} for 1,2 and 3-dim sde.
}
\examples{

## Example 1: Ito integral
## f(t,w(t)) = int(exp(w(t) - 0.5*t) * dw(s)) with t in [0,1]

fexpr <- expression( exp(w-0.5*t) )
res <- st.int(fexpr,type="ito",M=10,lower=0,upper=1)
res
## res$X 
summary(res)
## Display
plot(res,plot.type="single")
lines(time(res),mean(res),col=2,lwd=2)
lines(time(res),bconfint(res,level=0.95)[,1],col=4,lwd=2)
lines(time(res),bconfint(res,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95," confidence")),
        inset = .01,col=c(2,4),lwd=2,cex=0.8)

## Example 2: Stratonovich integral
## f(t,w(t)) = int(w(s)  o dw(s)) with t in [0,1]

fexpr <- expression( w )
res1 <- st.int(fexpr,type="str",M=10,lower=0,upper=1)
res1
## res1$X 
summary(res1)
## Display
plot(res1,plot.type="single")
lines(time(res1),mean(res1),col=2,lwd=2)
lines(time(res1),bconfint(res1,level=0.95)[,1],col=4,lwd=2)
lines(time(res1),bconfint(res1,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95," confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
