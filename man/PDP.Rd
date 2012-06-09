\name{PDP}
\alias{PDP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Pearson Diffusions Process (by Milstein Scheme)
}
\description{
Simulation the pearson diffusions process by milstein scheme.
}
\usage{
PDP(N, M, t0, T, x0, theta, mu, a, b, c, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{M}{
number of trajectories.
}
  \item{t0}{
initial time.
}
  \item{T}{
final time.
}
  \item{x0}{
initial value of the process at time \code{t0}.
}
  \item{theta}{
constant positive.
}
  \item{mu}{
constant.
}
  \item{a}{
constant.
}
  \item{b}{
constant.
}
  \item{c}{
constant.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
A class that further generalizes the Ornstein-Uhlenbeck and Cox-Ingersoll-Ross processes is the class of Pearson diffusion, the pearson diffusions process is the solution to the stochastic differential equation : \emph{\deqn{dX(t) = -theta *(X(t)-mu)*dt + sqrt( 2*theta*(a*X(t)^2 + b *X(t)+ c))* dW(t)}}
With \code{-theta *(X(t)-mu) :drift coefficient} and \code{sqrt( 2*theta*(a*X(t)^2 + b *X(t)+ c)) :diffusion coefficient}, \code{W(t)} is Wiener process, discretization \code{dt = (T-t0)/N}. 

With \code{theta > 0} and \code{a}, \code{b}, and \code{c} such that the diffusion coefficient is well-defined i.e., the square root can be extracted for all the values of the state space of \code{X(t)}.

\enumerate{\item When the \code{diffusion coefficient = sqrt(2*theta*c)} i.e, \code{(a=0,b=0)}, we recover the Ornstein-Uhlenbeck process.
\item For \code{diffusion coefficient = sqrt(2*theta*X(t))} and \code{0 < mu <= 1} i.e, \code{(a=0,b=1,c=0)}, we obtain the Cox-Ingersoll-Ross process, and if \code{mu > 1} the invariant distribution is a \code{Gamma law} with scale parameter \code{1} and shape parameter \code{mu}.
\item For \code{a > 0} and \code{diffusion coefficient = sqrt(2*theta*a*(X(t)^2+1))} i.e, \code{(b=0,c=a)}, the invariant distribution always exists on the real line, and for \code{mu = 0} the invariant distribution is a scaled \code{t distribution} with \code{v=(1+a^(-1))} degrees of freedom and scale parameter \code{v^(-0.5)}, while for \code{mu =! 0} the distribution is a form of skewed \code{t distribution} that is called Pearson type \code{IV} distribution.
\item For \code{a > 0}, \code{mu > 0}, and \code{diffusion coefficient = sqrt(2*theta*a*X(t)^2)} i.e, \code{(b=0,c=0)}, the distribution is defined on the positive half line and it is an \code{inverse Gamma distribution} with shape parameter \code{1 + a^-1} and scale parameter \code{a/mu}.
\item For \code{a > 0}, \code{mu >= a}, and \code{diffusion coefficient = sqrt(2*theta*a*X(t)*(X(t)+1))} i.e, \code{(b=a,c=0)}, the invariant distribution is the scaled \code{F distribution} with \code{(2*mu)/a} and \code{(2/a)+2} degrees of freedom and scale parameter \code{mu / (a+1)}. For \code{0 < mu < 1}, some reflecting conditions on the boundaries are also needed.
\item If \code{a < 0} and \code{mu > 0} are such that \code{min(mu,1-mu) >=  -a} and  \code{diffusion coefficient = sqrt(2*theta*a*X(t)*(X(t)-1))} i.e, \code{(b=-a,c=0)}, the invariant distribution exists on the interval \code{[0,1]} and is a \code{Beta distribution} with parameters \code{-mu/a} and \code{(mu-1)/a}.
}  
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{CEV}} Constant Elasticity of Variance Models, \code{\link{CIR}} Cox-Ingersoll-Ross Models, \code{\link{CIRhy}} modified CIR and hyperbolic Process, \code{\link{CKLS}} Chan-Karolyi-Longstaff-Sanders Models, \code{\link{DWP}} Double-Well Potential Model, \code{\link{GBM}} Model of Black-Scholes, \code{\link{HWV}} Hull-White/Vasicek Models, \code{\link{INFSR}} Inverse of Feller s Square Root models,  \code{\link{JDP}} Jacobi Diffusion Process, \code{\link{ROU}} Radial Ornstein-Uhlenbeck Process, \code{\link{diffBridge}} Diffusion Bridge Models, \code{\link{snssde}} Simulation Numerical Solution of SDE.
}
\examples{
\donttest{
## example 1
## theta = 5, mu = 10, (a=0,b=0,c=0.5)
## dX(t) = -5 *(X(t)-10)*dt + sqrt( 2*5*0.5)* dW(t)
 PDP(N=1000,M=1,T=1,t0=0,x0=1,theta=5,mu=10,a=0,b=0,c=0.5)

## example 2
## theta = 0.1, mu = 0.25, (a=0,b=1,c=0)
## dX(t) = -0.1 *(X(t)-0.25)*dt + sqrt( 2*0.1*X(t))* dW(t)
 PDP(N=1000,M=1,T=1,t0=0,x0=1,theta=0.1,mu=0.25,a=0,b=1,c=0)

## example 3
## theta = 0.1, mu = 1, (a=2,b=0,c=2)
## dX(t) = -0.1*(X(t)-1)*dt + sqrt( 2*0.1*(2*X(t)^2+2))* dW(t)
 PDP(N=1000,M=1,T=1,t0=0,x0=1,theta=0.1,mu=1,a=2,b=0,c=2)

## example 4
## theta = 0.1, mu = 1, (a=2,b=0,c=0)
## dX(t) = -0.1*(X(t)-1)*dt + sqrt( 2*0.1*2*X(t)^2)* dW(t)
 PDP(N=1000,M=1,T=1,t0=0,x0=1,theta=0.1,mu=1,a=2,b=0,c=0)

## example 5
## theta = 0.1, mu = 3, (a=2,b=2,c=0)
## dX(t) = -0.1*(X(t)-3)*dt + sqrt( 2*0.1*(2*X(t)^2+2*X(t)))* dW(t)
 PDP(N=1000,M=1,T=1,t0=0,x0=0.1,theta=0.1,mu=3,a=2,b=2,c=0)

## example 6
## theta = 0.1, mu = 0.5, (a=-1,b=1,c=0)
## dX(t) = -0.1*(X(t)-0.5)*dt + sqrt( 2*0.1*(-X(t)^2+X(t)))* dW(t)
 PDP(N=1000,M=1,T=1,t0=0,x0=0.1,theta=0.1,mu=0.5,a=-1,b=1,c=0)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
\keyword{ Stochastic Differential Equation }
\keyword{Numerical Solution of Stochastic Differential Equation}
