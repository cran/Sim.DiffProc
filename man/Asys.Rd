\name{Asys}
\alias{Asys}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evolution a Telegraphic Process in Time
}
\description{
Simulation the evolution of the telegraphic process (the availability of a system).
}
\usage{
Asys(lambda, mu, t, T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
the rate so that the system functions. 
}
  \item{mu}{
the rate so that the system is broken down.
}
  \item{t}{
calculate the matrix of transition \code{p(t)} has at the time \code{t}.
}
  \item{T}{
final time of evolution the process \code{[0,T]}.
}
}
\details{
Calculate the matrix of transition \code{p(t)} at time \code{t}, the space states of the telegraphic process is \code{(0,1)} with \code{0 : the system is broken down} and \code{1 : the system functions}, the initial distribution at time \code{t = 0} of the process is \code{p(t=0)=(1,0)} or \code{p(t=0)=(0,1)}.
}
\value{
matrix \code{p(t)} at time \code{t}, and plot of evolution the process.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{Telegproc}} simulation a telegraphic process. 
}
\examples{
## evolution a telegraphic process in time [0 , 5] 
## calculate the matrix of transition p(t = 10) 
Asys(0.5,0.5,10,5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Simulation }
