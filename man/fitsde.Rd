

\name{fitsde}
\alias{fitsde}
\alias{fitsde.default}
\alias{summary.fitsde}
\alias{print.fitsde}
\alias{vcov.fitsde}
\alias{AIC.fitsde}
\alias{BIC.fitsde}
\alias{logLik.fitsde}
\alias{coef.fitsde}
\alias{confint.fitsde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum Pseudo-Likelihood Estimation of 1-D SDE}
\description{
The (S3) generic function \code{"fitsde"} of estimate drift and diffusion parameters by the method of maximum pseudo-likelihood 
of the 1-dim stochastic differential equation. 
}
\usage{
fitsde(data, \dots)
\method{fitsde}{default}(data, drift, diffusion, start = list(), pmle = c("euler","kessler", 
   "ozaki", "shoji"), optim.method = "L-BFGS-B",
   lower = -Inf, upper = Inf, \dots)
   
\method{summary}{fitsde}(object, \dots)
\method{coef}{fitsde}(object, \dots)
\method{vcov}{fitsde}(object, \dots)
\method{logLik}{fitsde}(object, \dots)
\method{AIC}{fitsde}(object, \dots)
\method{BIC}{fitsde}(object, \dots)
\method{confint}{fitsde}(object,parm, level=0.95, \dots)   
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a univariate time series (\code{\link{ts}} class).}
  \item{drift}{drift coefficient: an \code{\link[base]{expression}} of two variables \code{t}, \code{x} and \code{theta} a vector of parameters of sde. See Examples.}
  \item{diffusion}{diffusion coefficient: an \code{\link[base]{expression}} of two variables \code{t}, \code{x} and \code{theta} a vector of parameters of sde. See Examples.}
  \item{start}{named list of starting values for optimizer. See Examples.}
  \item{pmle}{a \code{\link{character}} string specifying the method; can be either:
    \code{"euler"} (Euler pseudo-likelihood),
    \code{"ozaki"} (Ozaki pseudo-likelihood),
    \code{"shoji"} (Shoji pseudo-likelihood), and
	\code{"kessler"} (Kessler pseudo-likelihood).}
  \item{optim.method}{the \code{method} for \code{\link{optim}}.}
  \item{lower, upper}{bounds on the variables for the \code{"Brent"} or \code{"L-BFGS-B"} method.}
  \item{object}{an object inheriting from class \code{"fitsde"}.}
  \item{parm}{a specification of which parameters are to be given confidence intervals, either a vector of names (example \code{parm='theta1'}). If missing, all parameters are considered.}
  \item{level}{the confidence level required.}    
  \item{\dots}{potentially further arguments to pass to \code{\link{optim}}.}
}
\details{
The function \code{fitsde} returns a pseudo-likelihood estimators of the drift and diffusion parameters in 1-dim stochastic 
differential equation. The \code{\link{optim}} optimizer is used to find the maximum of the negative log pseudo-likelihood. An 
approximate covariance matrix for the parameters is obtained by inverting the Hessian matrix at the optimum.

The \code{pmle} of pseudo-likelihood can be one among:\code{"euler"}: Euler pseudo-likelihood), \code{"ozaki"}: Ozaki pseudo-likelihood,
\code{"shoji"}: Shoji pseudo-likelihood, and \code{"kessler"}: Kessler pseudo-likelihood.

An overview of this package, see \code{\link[utils]{browseVignettes}('Sim.DiffProc')} for more informations.
}
\value{
\code{fitsde} returns an object inheriting from \code{\link[base]{class}} \code{"fitsde"}.
}
\references{
 %  Elerian, O. (1998). 
 %  A note on the existence of a closed form conditional density for the Milstein scheme. 
 %  \emph{Working Paper, Nuffield College, Oxford University}. 
 %  Available at \url{http://www.nuff.ox.ac.uk/economics/papers/}
   
   Kessler, M. (1997).  
   Estimation of an ergodic diffusion from discrete observations. 
   \emph{Scand. J. Statist.}, 24, 211-229.

   Iacus, S.M. (2008). 
   \emph{Simulation and inference for stochastic differential equations: with R examples}. 
   Springer-Verlag, New York.

   Iacus, S.M. (2009). 
   sde: Simulation and Inference for Stochastic Differential Equations. 
   \emph{\R package version 2.0.10}. 
   
   Iacus, S.M. and all. (2014).
   The yuima Project: A Computational Framework for Simulation and Inference of Stochastic Differential Equations.
  \emph{Journal of Statistical Software}, \bold{57}(4).
   
   Ozaki, T. (1992). 
   A bridge between nonlinear time series models and nonlinear stochastic dynamical systems: A local linearization approach. 
   \emph{Statistica Sinica}, 2, 25-83.
   
   Shoji, L., Ozaki, T. (1998). 
   Estimation for nonlinear stochastic differential equations by a local linearization method. 
   \emph{Stochastic Analysis and Applications}, 16, 733-752.
   
  Dacunha, D.C. and Florens, D.Z. (1986). 
  Estimation of the Coefficients of a Diffusion from Discrete Observations. 
  \emph{Stochastics}. 19, 263--284.  

  Dohnal, G. (1987).
  On estimating the diffusion coefficient.
  \emph{J. Appl.Prob.}, 24, 105--114.

  Genon, V.C. (1990).
  Maximum constrast estimation for diffusion processes from discrete observation.
  \emph{Statistics}, 21, 99--116.
   
   Nicolau, J. (2004).
   Introduction to the estimation of stochastic differential equations based on discrete observations.
   \emph{Autumn School and International Conference, Stochastic Finance}.
  
   Ait-Sahalia, Y. (1999). 
   Transition densities for interest rate and other nonlinear diffusions.
   \emph{The Journal of Finance}, 54, 1361--1395.
   
   Ait-Sahalia, Y. (2002). 
   Maximum likelihood estimation of discretely sampled diffusions: a closed-form approximation approach. 
   \emph{Econometrica}. 70, 223--262.
   
   B.L.S. Prakasa Rao. (1999).
   \emph{Statistical Inference for Diffusion Type Processes}.
   Arnold, London and Oxford University press, New York.

   Kutoyants, Y.A. (2004).
   \emph{Statistical Inference for Ergodic Diffusion Processes}.
   Springer, London.
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{dcEuler}, \code{dcElerian}, \code{dcOzaki}, \code{dcShoji}, 
\code{dcKessler} and \code{dcSim} for approximated conditional law of a diffusion process. \code{gmm} estimator of the generalized method of moments by Hansen, and \code{HPloglik} these functions are useful 
to calculate approximated maximum likelihood estimators when the transition density of the process is not known, in package "sde". 

\code{qmle} in package "yuima" calculate quasi-likelihood and ML estimator of least squares estimator.

}
\examples{
#####  Example 1:

## Modele GBM (BS)
## dX(t) = theta1 * X(t) * dt + theta2 * x * dW(t)
## Simulation of data
set.seed(1234)

X <- GBM(N =1000,theta=4,sigma=1)
## Estimation: true theta=c(4,1)
fx <- expression(theta[1]*x)
gx <- expression(theta[2]*x)

fres <- fitsde(data=X,drift=fx,diffusion=gx,start = list(theta1=1,theta2=1),
              lower=c(0,0))
fres			   
summary(fres)
coef(fres)
logLik(fres)
AIC(fres)
BIC(fres)
vcov(fres)
confint(fres,level=0.95)
		  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fit}
\keyword{ts}
\keyword{sde}
