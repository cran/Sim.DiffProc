\name{Kern_general}
\alias{Kern_general}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjustment the Density of Random Variable by Kernel Methods
}
\description{
kernel density estimates. Its default method does so with the given kernel and bandwidth for univariate observations, and adjusted your density with distributions.
}
\usage{
Kern_general(Data, bw, k, Law = c("exp", "GAmma", "chisq", "Beta", 
            "fisher", "student", "weibull", "Normlog", "Norm"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
a numeric vector of the observed values.
}
  \item{bw}{
the smoothing bandwidth to be used. The kernels are scaled such that this is the standard deviation of the smoothing kernel.
bw=c('Irt','scott','Ucv','Bcv','SJ') or manual, see details \code{bw.nrd0}
}
  \item{k}{
a character string giving the smoothing kernel to be used. This must be one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine"
}
  \item{Law}{
distribution function with Adjusted. see details \code{Distributions} (R >= 2.12.1)
}
}
\details{
see details \code{density}
}
\value{
plot.density estimated with Adjustment.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
\itemize{\item bw='Irt'   ===> bw= bw.nrd0(X), implements a rule-of-thumb for choosing the bandwidth of a Gaussian kernel density estimator.\item bw='scott' ===> bw= bw.nrd(X) ,is the more common variation given by Scott.\item bw='Ucv'   ===> bw= bw.ucv(X) , implement unbiased cross-validation.\item bw='Bcv'   ===> bw= bw.bcv(X) , implement biased cross-validation.\item bw='SJ'    ===> bw= bw.SJ(X)  , implements the methods of Sheather & Jones.\item Choose your best distribution with minimum AIC.}
}
\seealso{
\code{\link{fctgeneral}} empirical distribution,\code{\link{hist_general}} Histograms Methods.
}
\examples{
##
X <- rexp(1000,1)
par(mfrow=c(2,2))
Kern_general(Data=X, bw='Irt', k="gaussian", Law = c("exp"))
Kern_general(Data=X, bw='scott', k="gaussian", Law = c("exp"))
Kern_general(Data=X, bw='Ucv', k="gaussian", Law = c("exp"))
Kern_general(Data=X, bw=0.3, k="gaussian", Law = c("exp"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
