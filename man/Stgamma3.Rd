\name{Stgamma3}
\alias{Stgamma3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Stochastic Process The Three-Parameter Gamma Distribution 
}
\description{
Simulation stochastic process by a three-parameter gamma distribution (also known as Pearson type III distribution).
}
\usage{
Stgamma3(N, t0, x0, T, shape, rate, thres, output = FALSE)
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
  \item{rate}{
rate parameter.
}
  \item{thres}{
threshold or shift parameter.
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
\code{\link{Stgamma}},\code{\link{Stweibull}},\code{\link{Stchisq}},\code{\link{Stcauchy}},\code{\link{Stbeta}},\code{\link{Stlnorm}},\code{\link{Stlnorm3}},
\code{\link{Stexp}},\code{\link{Stlgamma3}},\code{\link{Stweibull3}},\code{\link{Stlogis}},\code{\link{Stllogis}},\code{\link{Stllogis3}},\code{\link{Stgumbel}},
\code{\link{Stgp}},\code{\link{SRW}},\code{\link{WNG}},\code{\link{Stst}}.
}
\examples{
\donttest{
 Stgamma3(N=1000,t0=0,x0=0,T=1,shape=0.2,rate=0.1,thres=-5)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{Stochastic processes}
