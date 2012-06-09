\name{BMStra}
\alias{BMStra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stratonovitch Integral [1]
}
\description{
Simulation of the Stratonovitch \code{integral(W(s) o dW(s),0,t)}. 
}
\usage{
BMStra(N, T, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{T}{
final time.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
Stratonovitch integral as defined : \emph{\deqn{integral(f(t) o dW(s),0,t) = lim(sum(0.5*(f(t[i])+f(t[i+1]))*(W(t[i+1])-W(t[i]))))}}
calculus for Stratonovitch integral with \code{w(0) = 0}: \emph{\deqn{integral(W(s) o dW(s),0,t) = 0.5 * W(t)^2}}

The discretization \code{dt = T/N}, and  \code{W(t)} is Wiener process.
}
\value{
data frame(time,Stra) and plot of the Stratonovitch integral.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{BMStraC}} Stratonovitch Integral [2], \code{\link{BMStraP}} Stratonovitch Integral [3], \code{\link{BMStraT}} Stratonovitch Integral [4]. 
}
\examples{
\donttest{
 BMStra(N=1000, T=1, output = FALSE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic integral }
