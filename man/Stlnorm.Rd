\name{Stlnorm}
\alias{Stlnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Stochastic Process The Log Normal Distribution 
}
\description{
Simulation stochastic process by a log normal distribution.
}
\usage{
Stlnorm(N, t0, x0, T, meanlog, sdlog, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{N}{
size of process.
}
  \item{t0}{
initial time.
}
  \item{x0}{
initial value at time \code{t0}.
}
  \item{T}{
final time.
}
  \item{meanlog}{
mean of the distribution on the log scale.
}
  \item{sdlog}{
standard deviation of the distribution on the log scale.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Guidoum Arsalane.
}
\note{
\itemize{\item Time step of the simulation \code{(T-t0)/N}.}
}
\seealso{
\code{\link{Stgamma}},\code{\link{Stweibull}},\code{\link{Stchisq}},\code{\link{Stcauchy}},\code{\link{Stbeta}},\code{\link{Stllogis3}},\code{\link{Stlnorm3}},
\code{\link{Stexp}},\code{\link{Stlgamma3}},\code{\link{Stweibull3}},\code{\link{Stlogis}},\code{\link{Stlgamma3}},\code{\link{Stllogis}},\code{\link{Stgp}},
\code{\link{Stgamma3}},\code{\link{SRW}},\code{\link{WNG}},\code{\link{Stst}}.
}
\examples{
\donttest{
 Stlnorm(N=1000,t0=0,x0=0,T=1,meanlog=1,sdlog=1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{Stochastic processes}

