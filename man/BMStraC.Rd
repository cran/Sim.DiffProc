\name{BMStraC}
\alias{BMStraC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stratonovitch Integral [2]
}
\description{
Simulation of the Stratonovitch \code{integral(alpha o dW(s),0,t)}. 
}
\usage{
BMStraC(N, T, alpha, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{T}{
final time.
}
  \item{alpha}{
constant.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
Stratonovitch integral as defined : \emph{\deqn{integral(f(t) o dW(s),0,t) = lim(sum(0.5*(f(t[i])+f(t[i+1]))*(W(t[i+1])-W(t[i]))))}}
calculus for Stratonovitch integral with \code{w(0) = 0}: \emph{\deqn{integral(alpha o dW(s),0,t) = alpha * W(t)}}

The discretization \code{dt = T/N}, and  \code{W(t)} is Wiener process.
}
\value{
data frame(time,Stra) and plot of the Stratonovitch integral.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{BMStra}} Stratonovitch Integral [1], \code{\link{BMStraP}} Stratonovitch Integral [3], \code{\link{BMStraT}} Stratonovitch Integral [4].
 }
\examples{
\donttest{
 BMStraC(N=1000, T=1, alpha = 2,output = FALSE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic integral }
