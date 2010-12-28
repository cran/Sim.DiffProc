\name{BMItoP}
\alias{BMItoP}
\alias{ItovsStraP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Properties of the stochastic integral and Ito Process [4]
}
\description{
Simulation of the Ito \code{integral(W(s)^n*dW(s),0,t)}. 
}
\usage{
BMItoP(N, T, power, output = FALSE)
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
However the Ito integral also has the peculiar property, amongst others, that : \emph{\deqn{integral(W(s)^n*dW(s),0,t) = W(t)^(n+1) /(n+1) - (n/2) * integral(W(s)^n-1 *ds,0,t)}}
from classical calculus for Ito integral with \code{w(0) = 0}.

The  follows from the algebraic rearrangement : \emph{\deqn{integral(W(s)^n*dW(s),0,t) = sum ( W(t)^n*(W(t+1)-W(t)),0,t)}}
}
\value{
data frame(time,Ito,sum.Ito) and plot of the Ito integral.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{BMIto1}} simulation of the Ito integral[1], \code{\link{BMIto2}} simulation of the Ito integral[2], \code{\link{BMItoC}} properties of the stochastic integral and Ito processes[3], \code{\link{BMItoT}} properties of the stochastic integral and Ito processes[5].
}
\examples{
## if power = 1
## integral(W(s) * dW(s),0,t) = W(t)^2/2 - 1/2 * t
BMItoP(N=1000,T=1,power =1)
## if power = 2
## integral(W(s)^2 * dW(s),0,t) = W(t)^3/3 - 2/2 * integral(W(s)*ds,0,t)
BMItoP(N=1000,T=1,power =2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Stochastic integral }
