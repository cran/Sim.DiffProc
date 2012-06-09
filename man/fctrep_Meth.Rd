\name{fctrep_Meth}
\alias{fctrep_Meth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating the Empirical Distribution of Random Variable X 
}
\description{
Calculating your empirical distribution of random variable X.
}
\usage{
fctrep_Meth(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of the observed values.
}
}
\details{
calculating the empirical distribution \code{F[i]= (1/n)*Sum(V[i])} with \code{V[i]= 1 if x[i] <= X else V[i] = 0}. 
}
\value{
Plot the empirical distribution.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{hist_meth}} Histograms,\code{\link{Kern_meth}} Kernel Density.
}
\examples{
 X <- rexp(1000,2)
 Y <- rgamma(1000,1,2)
 Z <- rweibull(1000,1,1)
 G <- rnorm(1000,mean(X),sd(X))
 par(mfrow=c(2,2))
 fctrep_Meth(X)
 fctrep_Meth(Y)
 fctrep_Meth(Z)
 fctrep_Meth(G)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
