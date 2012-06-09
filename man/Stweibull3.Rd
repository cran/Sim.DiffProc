\name{Stweibull3}
\alias{Stweibull3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Stochastic Process The Three-Parameter Weibull Distribution 
}
\description{
Simulation stochastic process by a three-parameter weibull distribution.
}
\usage{
Stweibull3(N, t0, x0, T, shape, scale, thres, output = FALSE)
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
  \item{shape}{
shape parameter.
}
  \item{scale}{
scale parameter.
}
  \item{thres}{
thres parameter.
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
\code{\link{Stgamma}},\code{\link{Stlogis}},\code{\link{Stchisq}},\code{\link{Stcauchy}},\code{\link{Stbeta}},\code{\link{Stllogis3}},\code{\link{Stlnorm}},
\code{\link{Stexp}},\code{\link{Stlgamma3}},\code{\link{Stweibull}},\code{\link{Stlnorm3}},\code{\link{Stlgamma3}},\code{\link{Stllogis}},\code{\link{Stgp}},
\code{\link{Stgamma3}},\code{\link{SRW}},\code{\link{WNG}},\code{\link{Stst}}.
}
\examples{
\donttest{
 Stweibull3(N=1000,t0=0,x0=0,T=1,shape=1,scale=2,thres=-1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{Stochastic processes}
