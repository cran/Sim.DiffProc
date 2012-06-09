\name{Kern_meth}
\alias{Kern_meth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel Density of Random Variable X 
}
\description{
kernel density estimates. Its default method does so with the given kernel and bandwidth for univariate observations. 
}
\usage{
Kern_meth(X, bw, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of the observed values.
}
  \item{bw}{
the smoothing bandwidth to be used. The kernels are scaled such that this is the standard deviation of the smoothing kernel.
bw=c('Irt','scott','Ucv','Bcv','SJ') or manual, see details \code{bw.nrd0}
}
  \item{k}{
a character string giving the smoothing kernel to be used. This must be one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine"
}
}
\details{
see details \code{plot.density}
}
\value{
plot.density for your data.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\note{
\itemize{\item bw='Irt'   ===> bw= bw.nrd0(X), implements a rule-of-thumb for choosing the bandwidth of a Gaussian kernel density estimator.\item bw='scott' ===> bw= bw.nrd(X) ,is the more common variation given by Scott.\item bw='Ucv'   ===> bw= bw.ucv(X) , implement unbiased cross-validation.\item bw='Bcv'   ===> bw= bw.bcv(X) , implement biased cross-validation.\item bw='SJ'    ===> bw= bw.SJ(X)  , implements the methods of Sheather & Jones.}
}
\seealso{
\code{\link{hist_meth}} Histograms,\code{\link{fctrep_Meth}} Empirical Distribution.
}
\examples{
\donttest{
## Example 1
## fixed bw with different kernel
 X <- rbeta(1000,1,2)
 par(mfrow=c(2,2))
 Kern_meth(X, bw='Ucv', k="rectangular")
 Kern_meth(X, bw='Ucv',k="triangular")
 Kern_meth(X, bw='Ucv',k="epanechnikov")
 Kern_meth(X, bw='Ucv',k="cosine")

## Example 2
## fixed kernel with different bw
 Y <- rlnorm(1000)
 par(mfrow=c(2,2))
 Kern_meth(Y, bw='Irt', k="epanechnikov")
 Kern_meth(Y, bw='Ucv',k="epanechnikov")
 Kern_meth(Y, bw='scott',k="epanechnikov")
 Kern_meth(Y, bw=0.4,k="epanechnikov")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
