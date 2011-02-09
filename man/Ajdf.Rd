\name{Ajdf}
\alias{Ajdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjustment By F Distribution
}
\description{
Adjusted your sample by the F law, estimated these parameters using the method of maximum likelihood, and
calculating the Akaike information criterion for one or several fitted model objects for which a log-likelihood value can be obtained, 
according to the formula -2*log-likelihood + k*npar, where npar represents the number of parameters in the fitted model, and k = 2 for the usual AIC,
and computes confidence intervals for one or more parameters in a fitted model (Law).
}
\usage{
Ajdf(X, starts = list(df1 = 1, df2 = 1), leve = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of the observed values.
}
  \item{starts}{
named list. Initial values for optimizer.
}
  \item{leve}{
the confidence level required.
}
}
\details{
The \code{optim} optimizer is used to find the minimum of the negative log-likelihood. An approximate covariance matrix for the parameters is obtained by inverting the Hessian matrix at the optimum.

For more detail consulted \code{mle},\code{confint},\code{AIC}.

\code{R} has the \code{[dqpr]f} functions to evaluate the density, the quantiles, and the cumulative distribution or generate pseudo random numbers from the F distribution.
}
\value{
\item{coef}{Coefficients extracted from the model.}
\item{AIC}{A numeric value with the corresponding AIC.}
\item{vcov}{A matrix of the estimated covariances between the parameter estimates in the linear or non-linear predictor of the model.}
\item{confint}{A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2.
}
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{Ajdchisq}} Adjustment By Chi-Squared Distribution,\code{\link{Ajdexp}} Adjustment By Exponential Distribution,
\code{\link{Ajdgamma}} Adjustment By Gamma Distribution,\code{\link{Ajdlognorm}} Adjustment By Log Normal Distribution,
\code{\link{Ajdnorm}} Adjustment By Normal Distribution,\code{\link{Ajdt}} Adjustment By Student t Distribution,
\code{\link{Ajdweibull}} Adjustment By Weibull Distribution,\code{\link{Ajdbeta}} Adjustment By Beta Distribution.
}
\examples{
X <- rf(100,df1=5,df2=5)
Ajdf(X, starts = list(df1 = 1, df2 = 1), leve = 0.95)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
