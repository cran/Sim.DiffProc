\name{BMStraP}
\alias{BMStraP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stratonovitch Integral [3]
}
\description{
Simulation of the Stratonovitch \code{integral(W(s)^n o dW(s),0,t)}.
}
\usage{
BMStraP(N, T, power, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{T}{
final time.
}
  \item{power}{
constant.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
Stratonovitch integral as defined : \emph{\deqn{integral(f(t) o dW(s),0,t) = lim(sum(0.5*(f(t[i])+f(t[i+1]))*(W(t[i+1])-W(t[i]))))}}
calculus for Stratonovitch integral with \code{w(0) = 0}: \emph{\deqn{integral(W(s)^n o dW(s),0,t) = lim(sum(0.5*(W(t[i])^(n-1)+W(t[i+1])^(n-1))*(W(t[i+1])^2 - W(t[i])^2)))}}

The discretization \code{dt = T/N}, and  \code{W(t)} is Wiener process.
}
\value{
data frame(time,Stra) and plot of the Stratonovitch integral.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{BMStra}} Stratonovitch Integral [1], \code{\link{BMStraC}} Stratonovitch Integral [2], \code{\link{BMStraT}} Stratonovitch Integral [4].
}
\examples{
##
BMStraP(N=1000, T=1, power = 2,output = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic integral }
