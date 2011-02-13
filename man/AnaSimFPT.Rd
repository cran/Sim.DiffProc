\name{AnaSimFPT}
\alias{AnaSimFPT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation The First Passage Time FPT For A Simulated Diffusion Process
}
\description{
Simulation M-samples of the first passage time (FPT) by a simulated diffusion process with a fixed the threshold v.
}
\usage{
AnaSimFPT(N, M, t0, Dt, T = 1, X0, v, drift, diff, 
          ELRENA=c("No","Yes","Mean","Median"), 
          Output = FALSE, Methods = c("Euler", "Milstein", 
          "MilsteinS","Ito-Taylor", "Heun", "RK3"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of the diffusion process.
}
  \item{M}{
size of the FPT.
}
  \item{t0}{
initial time.
}
  \item{Dt}{
time step of the simulation \code{(discretization)}.
}
  \item{T}{
final time.
}
  \item{X0}{
initial value of the process at time \code{t0}.
}
  \item{v}{
threshold (Risk).
}
  \item{drift}{
drift coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{diff}{
diffusion coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{ELRENA}{
if \code{ELRENA = "No"} not eliminate \code{NA} (Not Available),and if \code{ELRENA="Yes"} eliminate \code{NA} (Not Available),
or replace NA  by : \code{mean(FPT)} ,\code{median(FPT)}.
}
  \item{Output}{
if \code{Output = TRUE} write a \code{Output} to an Excel 2007.
}
  \item{Methods}{
method of simulation ,see details \code{\link{snssde}}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
The stochastic differential equation of  is  : \emph{\deqn{dX(t) = a(t,X(t)) *dt + b(t,X(t)) *dW(t)}} with \code{a(t,X(t)) :drift coefficient} and \code{b(t,X(t)) :diffusion coefficient}, \code{W(t)} is Wiener process.

We take interest in the random variable \code{tau} "first passage time", is defined by : \emph{\deqn{tau = inf(t>=0 \ X(t) <= v   OR  X(t) >= v)}} with \code{v is the threshold}.

For more detail consulted \code{References}.
}
\value{
Random variable tau "FPT".
}
\references{
\enumerate{\item K.Boukhetala, Estimation of the first passage time distribution for a simulated diffusion process, Maghreb Math.Rev, Vol.7, No 1, Jun 1998, pp. 1-25. \item K.Boukhetala, Simulation study of a dispersion about an attractive centre. In proceedings of 11th Symposium Computational Statistics, edited by R.Dutter and W.Grossman, Wien , Austria, 1994, pp. 128-130. \item K.Boukhetala,Modelling and simulation of a dispersion pollutant with attractive centre, Edited by Computational Mechanics Publications, Southampton ,U.K and Computational Mechanics Inc, Boston, USA, pp. 245-252.\item K.Boukhetala, Kernel density of the exit time in a simulated diffusion, les Annales Maghrebines De L ingenieur, Vol , 12, N Hors Serie. Novembre 1998, Tome II, pp 587-589. 
}
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\note{
Time of Calculating

The Ornstein-Uhlenbeck Process (example)
drift <- expression(-5*x)
diff <- expression(1)

system.time(AnaSimFPT(N=1000, M=30, t0=0, Dt=0.001, T = 1, X0=10,
             v=0.05,drift, diff, ELRENA ="No", Output = FALSE))
					  
utilisateur     systeme      ecoule 

       1.89        0.55        2.62

system.time(AnaSimFPT(N=1000, M=100, t0=0, Dt=0.001, T = 1, X0=10, 
            v=0.05,drift, diff, ELRENA ="No", Output = FALSE))
					  
utilisateur     systeme      ecoule

       5.74        1.64        7.78

system.time(AnaSimFPT(N=1000, M=500, t0=0, Dt=0.001, T = 1, X0=10,
            v=0.05,drift, diff, ELRENA ="Mean", Output = FALSE))
					  
utilisateur     systeme      ecoule

      26.07        7.78       37.93

system.time(AnaSimFPT(N=1000, M=500, t0=0, Dt=0.001, T = 1, X0=10,
            v=0.05,drift, diff, ELRENA ="Mean", 
            Output = FALSE,Methods="RK3"))
					  
utilisateur     systeme      ecoule 

     125.64        8.90      150.85
}
\seealso{
\code{\link{AnaSimX}} Simulation M-Samples of Random Variable X(v[t]) For A Simulated Diffusion Process, \code{\link{tho_M1}} Simulation The FPT For Attractive Model(S = 1,Sigma), \code{\link{tho_M1}} Simulation The FPT For Attractive Model(S >= 2,Sigma), \code{\link{tho_02diff}} Simulation FPT For Attractive Model for 2-Diffusion Processes.  
}
\examples{
## Example 1
## tau = inf(t>=0 \ X(t) <= v 
## Ornstein-Uhlenbeck Process  or Gaussian Diffusion Models
## v = 0.05 
drift <- expression(5*(-2-x))
diff <- expression(1)
AnaSimFPT(N=1000, M=30, t0=0, Dt=0.001, T = 1, X0=10, v=0.05, drift, 
          diff,ELRENA ="No", Output = FALSE)
summary(tau)
hist(tau)
plot(density(tau,kernel ="gaussian"),col="red")
## v = -0.05
AnaSimFPT(N=1000, M=30, t0=0, Dt=0.001, T = 1, X0=10, v=-0.05, drift, 
          diff,ELRENA ="No", Output = FALSE)
summary(tau)
hist(tau)
plot(density(tau,kernel ="gaussian"),col="red")
## Attention
## v = -3
## AnaSimFPT(N=1000, M=30, t0=0, Dt=0.001, T = 1, X0=10, v=-3, drift,
##           diff,ELRENA ="No", Output = FALSE)

## Example 2
## tau = inf(t>=0 \ X(t) >= v )
## v = 1
drift <- expression(2*(3-x))
diff <- expression(0.1)
AnaSimFPT(N=1000, M=30, t0=0, Dt=0.001, T = 1, X0=-5, v=1, drift, 
          diff,ELRENA ="No", Output = FALSE)
summary(tau)
hist(tau)
plot(density(tau,kernel ="gaussian"),col="red")
## v = 3 
AnaSimFPT(N=1000, M=30, t0=0, Dt=0.01, T = 1, X0=-5, v=3, drift, 
          diff,ELRENA ="No", Output = FALSE)
summary(tau)
hist(tau)
plot(density(tau,kernel ="gaussian"),col="red")
## v = 3.1
##AnaSimFPT(N=1000, M=30, t0=0, Dt=0.01, T = 1, X0=-5, v=3.1, drift, 
##          diff,ELRENA ="No", Output = FALSE)
## Remplaced NA by mean(tau) or median(tau)
##AnaSimFPT(N=1000, M=30, t0=0, Dt=0.01, T = 1, X0=-5, v=3.1, drift, 
##          diff,ELRENA ="Yes", Output = FALSE)
##AnaSimFPT(N=1000, M=30, t0=0, Dt=0.01, T = 1, X0=-5, v=3.1, drift, 
##          diff,ELRENA ="Mean", Output = FALSE)
##AnaSimFPT(N=1000, M=30, t0=0, Dt=0.01, T = 1, X0=-5, v=3.1, drift, 
##          diff,ELRENA ="Median", Output = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Attractive models }
\keyword{ Stochastic Differential Equation }
\keyword{ Statistical Analysis }
