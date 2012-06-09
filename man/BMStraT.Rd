\name{BMStraT}
\alias{BMStraT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stratonovitch Integral [4]
}
\description{
Simulation of the Stratonovitch \code{integral(s o dW(s),0,t)}.
}
\usage{
BMStraT(N, T, output = FALSE)
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
calculus for Stratonovitch integral with \code{w(0) = 0}: \emph{\deqn{integral(s o dW(s),0,t) = lim(sum(0.5*(t[i]*(W(t[i+1])-W(t[i]))+t[i+1]*(W(t[i+1])-W(t[i])))))}}

The discretization \code{dt = T/N}, and  \code{W(t)} is Wiener process.
}
\value{
data frame(time,Stra) and plot of the Stratonovitch integral.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{BMStra}} Stratonovitch Integral [1], \code{\link{BMStraC}} Stratonovitch Integral [2], \code{\link{BMStraC}} Stratonovitch Integral [3].
}
\examples{
\donttest{
 BMStraT(N=1000, T=1,output = FALSE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic integral }
