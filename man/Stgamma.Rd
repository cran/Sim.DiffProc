\name{Stgamma}
\alias{Stgamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Stochastic Process The Gamma Distribution
}
\description{
Simulation stochastic process by a gamma distribution.
}
\usage{
Stgamma(N, t0, T, alpha, beta, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{t0}{
initial time.
}
  \item{T}{
final time.
}
  \item{alpha}{
constant positive.
}
  \item{beta}{
an alternative way to specify the scale.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\note{
\itemize{\item Time step of the simulation \code{T/N}.}
}
\seealso{
\code{\link{Stbeta}},\code{\link{Stweibull}},\code{\link{Stexp}},\code{\link{Stchisq}},\code{\link{Stcauchy}},\code{\link{Stlnorm}},\code{\link{Stlnorm3}},
\code{\link{Stgamma3}},\code{\link{Stlgamma3}},\code{\link{Stweibull3}},\code{\link{Stlogis}},\code{\link{Stllogis}},\code{\link{Stllogis3}},\code{\link{Stgumbel}},
\code{\link{Stgp}},\code{\link{SRW}},\code{\link{WNG}},\code{\link{Stst}}.
}
\examples{
\donttest{
## Stochastic Process The Gamma Distribution
 Stgamma(N=1000,t0=0,T=5,alpha=1,beta=1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{Stochastic processes}
