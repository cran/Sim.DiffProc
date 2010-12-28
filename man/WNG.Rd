\name{WNG}
\alias{WNG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating White Noise Gaussian 
}
\description{
Simulation white noise gaussian. 
}
\usage{
WNG(N, t0, T, m, sigma2, output = FALSE)
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
  \item{m}{
mean.
}
  \item{sigma2}{
variance.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\value{
data.frame(time,x) and plot of process. 
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\examples{
## White Noise Gaussian 
WNG(N=1000,t0=0,T=1,m=0,sigma2=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
