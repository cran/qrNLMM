\name{QRNLMM}
\alias{QRNLMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile Regression for Nonlinear Mixed-Effects Models
}
\description{
Performs a quantile regression for a NLMEM using the Stochastic-Approximation of the EM Algorithm (SAEM) for an unique or a set of quantiles.
}
\usage{
QRNLMM(y,x,groups,initial,exprNL,covar=NA,p=0.5,precision=0.0001,MaxIter=500,
       M=20,cp=0.25,beta=NA,sigma=NA,Psi=NA,show.convergence=TRUE,CI=95,verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the response vector of dimension \eqn{N} where \eqn{N} is the total of observations.
}
  \item{x}{
vector of longitudinal (repeated measures) covariate of dimension \eqn{N}. For example: Time, location, etc.  
}
  \item{groups}{
factor of dimension \eqn{N} specifying the partitions of the data over which the random effects vary.
}
  \item{initial}{an numeric vector, or list of initial estimates
   for the fixed effects. It must be provide adequately (see \code{details} section) in order to ensure a proper convergence.
   }
  \item{exprNL}{
expression containing the proposed nonlinear function. It can be of class \code{character} or \code{expression}. It must have a defined structure defined in the \code{details} section in order to be correctly read by the derivate R function \code{deriv}.
}
  \item{covar}{
a matrix of dimension \eqn{N \times r} where \eqn{r} represents the number of covariates.
}
 \item{p}{
unique quantile or a set of quantiles related to the quantile regression.
}
  \item{precision}{
the convergence maximum error.
}
  \item{MaxIter}{
the maximum number of iterations of the SAEM algorithm. Default = 500.
}
  \item{M}{
Number of Monte Carlo simulations used by the SAEM Algorithm. Default = 20. For more accuracy we suggest to use \code{M=20}.
}
  \item{cp}{
cut point \eqn{(0 \le cp \le 1)} which determines
the percentage of initial iterations with no memory.
}
  \item{beta}{
fixed effects vector of initial parameters, if desired. 
}
  \item{sigma}{
dispersion initial parameter for the error term, if desired.
}
  \item{Psi}{
Variance-covariance random effects matrix of initial parameters, if desired.
}
  \item{show.convergence}{
if \code{TRUE}, it will show a graphical summary for the convergence of the estimates of all parameters for each quantile in order to assess the convergence.
}
  \item{CI}{
Confidence to be used for the Confidence Interval when a grid of quantiles is provided. Default=95.
}
  \item{verbose}{
if \code{TRUE}, an output summary is printed.
}
}
\details{
This algorithm performs the SAEM algorithm proposed by Delyon et al. (1999), a stochastic version of the usual EM Algorithm deriving exact maximum likelihood estimates of the fixed-effects and variance components. Covariates are allowed, the longitudinal (repeated measures) coded \code{x} and a set of covariates \code{covar}.
 
\eqn{About initial values:} Estimation for fixed effects parameters envolves a Newton-Raphson step. In adition, NL models are highly sensitive to initial values. So, we suggest to set of intial values quite good, this based in the parameter interpretation of the proposed NL function.

\eqn{About the nonlinear expression:} For the NL expression \code{exprNL} just the variables \code{x}, \code{covar}, \code{fixed} and \code{random} can be defined. Both \code{x} and \code{covar} represent the covariates defined above. The fixed effects must be declared as \code{fixed[1]}, \code{fixed[2]},..., \code{fixed[d]} representing the first, second and \eqn{d}th fixed effect. Exactly the same for the random effects and covariates where the term \code{fixed} should be replace for \code{random} and \code{covar} respectively.

For instance, if we use the exponential nonlinear function with two parameters, each parameter represented by a fixed and a random effect, this will be defined by

\deqn{y_{ij} = (\beta_1 + b_1)\exp^{-(\beta_2 + b_2)x_{ij}}}

and the \code{exprNL} should be a character or and expression defined by

\code{exprNL = "(fixed[1]+random[1])*exp(-(fixed[2]+random[2])*x)"}

or

\code{exprNL = expression((fixed[1]+random[1])*exp(-(fixed[2]+random[2])*x))}.


If we are interested in adding two covariates in order to explain on of the parameters, the covariates \code{covar[1]} and \code{covar[2]} must be included in the model. For example, for the nonlinear function

\deqn{y_{ij} = (\beta_1 + \beta_3*covar1_{ij} + b_1)\exp^{-(\beta_2 + \beta_4*
covar2_{ij} + b_2)x_{ij}}}

the \code{exprNL} should be

\code{exprNL = "(fixed[1]+fixed[3]*covar[1]+random[1])*exp(-(fixed[2]+fixed[4]*covar[2]+random[2])*x)"}

or

\code{exprNL = expression((fixed[1]+fixed[3]*covar[1]+random[1])*exp(-(fixed[2]+
fixed[4]*covar[2]+random[2])*x))}.

Note that the mathematical function \code{exp} was used. For derivating the \code{deriv} R function recognizes in the \code{exprNL} expression the arithmetic operators +, -, *, / and ^, and the single-variable functions \code{exp}, \code{log}, \code{sin}, \code{cos}, \code{tan}, \code{sinh}, \code{cosh}, \code{sqrt}, \code{pnorm}, \code{dnorm}, \code{asin}, \code{acos}, \code{atan}, \code{gamma}, \code{lgamma}, \code{digamma} and \code{trigamma}, as well as \code{psigamma} for one or two arguments (but derivative only with respect to the first).


\eqn{General details:} When a grid of quantiles is provided, a graphical summary with point estimates and Confidence Intervals for model parameters is shown and also a graphical summary for the convergence of these estimates (for each quantile), if \code{show.convergence=TRUE}.

If the convergence graphical summary shows that convergence has not be attained, it's suggested to increase the total number of iterations \code{MaxIter}.

About the cut point parameter \code{cp}, a number between 0 and 1 \eqn{(0 \le cp \le 1)} will assure an initial convergence in distribution to a solution neighborhood for the first \code{cp}*\code{MaxIter} iterations and an almost sure convergence for the rest of the iterations. If you do not know how SAEM algorithm works, these parameters SHOULD NOT be changed.

This program uses progress bars that will close when the algorithm ends. They must not be closed before, if not, the algorithm will stop.
}
\value{
The function returns a list with two objects
\item{conv}{A two elements list with the matrices \code{teta} and \code{se} containing the point estimates and standard error estimate for all parameters along all iterations.}

The second element of the list is \code{res}, a list of 13 elements detailed as

\item{p}{quantile(s) fitted.}
\item{iter}{number of iterations.}
\item{criteria}{attained criteria value.}
\item{nlmodel}{the proposed nonlinear function.}
\item{beta}{fixed effects estimates.}
\item{weights}{random effects weights (\eqn{b_i}).}
\item{sigma}{scale parameter estimate for the error term.}
\item{Psi}{Random effects variance-covariance estimate matrix.}
\item{SE}{Standard Error estimates.}
\item{table}{Table containing the inference for the fixed effects parameters.}
\item{loglik}{Log-likelihood value.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{HQ}{Hannan-Quinn information criterion.}
\item{fitted.values}{vector containing the fitted values}
\item{residuals}{vector containing the residuals.}
\item{time}{processing time.}
}
\references{
Galarza, C.E., Castro, L.M., Louzada, F. & Lachos, V. (2020) Quantile regression for nonlinear mixed effects models: a likelihood based perspective. Stat Papers 61, 1281-1307.
\doi{10.1007/s00362-018-0988-y}

Delyon, B., Lavielle, M. & Moulines, E. (1999). Convergence of a stochastic approximation version of the EM algorithm. Annals of Statistics, pages 94-128.
}
\author{
Christian E. Galarza <\email{chedgala@espol.edu.ec}> and 
Victor H. Lachos <\email{hlachos@uconn.edu}>
}
\note{
If a grid of quantiles was provided, the result is a list of the same dimension where each element corresponds to each quantile as detailed above.
}

\seealso{
\code{\link{Soybean}}, \code{\link{HIV}}, \code{\link[lqr]{lqr}} , \code{\link{group.plots}}
}
\examples{
\dontrun{
#Using the Soybean data

data(Soybean)
attach(Soybean)

#################################
#A full model (no covariate)

y     = weight          #response
x     = Time            #time

#Expression for the three parameter logistic curve

exprNL = expression((fixed[1]+random[1])/
                    (1 + exp(((fixed[2]+random[2])- x)/(fixed[3]+random[3]))))

#Initial values for fixed effects
initial = c(max(y),0.6*max(y),0.73*max(y))

#A median regression (by default)
median_reg = QRNLMM(y,x,Plot,initial,exprNL)

#Assesing the fit

fxd     = median_reg$res$beta
nlmodel = median_reg$res$nlmodel
weights = median_reg$res$weights
seqc    = seq(min(x),max(x),length.out = 500)

group.plot(x = Time,y = weight,groups = Plot,type="l",
           main="Soybean profiles",xlab="time (days)",
           ylab="mean leaf weight (gr)",col="gray")

for(i in 1:nlevels(Plot)){
  lines(seqc,nlmodel(x = seqc,fixed = fxd,random = weights[i,]),lty=2)
}

# median population curve
lines(seqc,nlmodel(x = seqc,fixed = fxd,random = rep(0,ncol(weights))),
      lwd=3,col="blue")         

#########################################
#A model for compairing the two genotypes

y     = weight          #response
x     = Time            #time
covar = as.numeric(Variety)-1    #factor genotype (0=Forrest, 1=Plan Introduction)

#Expression for the three parameter logistic curve with a covariate in the asymp growth
exprNL = expression((fixed[1]+(fixed[4]*covar[1])+random[1])/
                      (1 + exp(((fixed[2]+random[2])- x)/(fixed[3]+random[3]))))

#Initial values for fixed effects
initial = c(max(y),0.6*max(y),0.73*max(y),3)

# A quantile regression for the three quartiles (just 200 iterations)
box_reg = QRNLMM(y,x,Plot,initial,exprNL,covar,p=c(0.25,0.50,0.75),MaxIter = 200)

# Assing the fit for the first quartile Q1 curve
fxd_q1    = box_reg[[1]]$res$beta
fxd_q2    = box_reg[[2]]$res$beta #median
fxd_q3    = box_reg[[3]]$res$beta
nlmodel = box_reg[[1]]$res$nlmodel
seqc    = seq(min(x),max(x),length.out = 500)


group.plot(x = Time[Variety=="P"],y = weight[Variety=="P"],
           groups = Plot[Variety=="P"],type="l",col="lightblue",
           main="Soybean profiles by genotype",xlab="time (days)",
           ylab="mean leaf weight (gr)")

group.lines(x = Time[Variety=="F"],y = weight[Variety=="F"],
            groups = Plot[Variety=="F"],col="gray")

# Add the three quantile lines
lines(seqc,nlmodel(x = seqc, fixed = fxd_q1, random = rep(0, 3), covar = 1),
      lwd = 2, col = "blue", lty = "dashed") # Q1, dashed
lines(seqc, nlmodel(x = seqc, fixed = fxd_q2, random = rep(0, 3), covar = 1),
      lwd = 2, col = "blue")             # Median, solid
lines(seqc, nlmodel(x = seqc, fixed = fxd_q3, random = rep(0, 3), covar = 1),
      lwd = 2, col = "blue", lty = "dashed") # Q3, dashed

lines(seqc,nlmodel(x = seqc,fixed = fxd_q1,random = rep(0,3),covar=0),
      lwd=2,col="black",lty="dashed") #q1
lines(seqc,nlmodel(x = seqc,fixed = fxd_q2,random = rep(0,3),covar=0),
      lwd=2,col="black")             # Median
lines(seqc,nlmodel(x = seqc,fixed = fxd_q3,random = rep(0,3),covar=0),
      lwd=2,col="black",lty="dashed") #q3


#########################################
#A simple output example

---------------------------------------------------
Quantile Regression for Nonlinear Mixed Model
---------------------------------------------------
Quantile = 0.5
Subjects = 48 ; Observations = 412

- Nonlinear function

function(x,fixed,random,covar=NA){
  resp = (fixed[1] + random[1])/(1 + exp(((fixed[2] +
          random[2]) - x)/(fixed[3] + random[3])))
  return(resp)}

-----------
Estimates
-----------
- Fixed effects

       Estimate Std. Error   z value Pr(>|z|)
beta 1 18.80029    0.53098  35.40704        0
beta 2 54.47930    0.29571 184.23015        0
beta 3  8.25797    0.09198  89.78489        0

sigma = 0.31569

Random effects Variance-Covariance Matrix matrix
         b1       b2      b3
b1 24.36687 12.27297 3.24721
b2 12.27297 15.15890 3.09129
b3  3.24721  3.09129 0.67193

------------------------
Model selection criteria
------------------------
        Loglik      AIC      BIC       HQ
Value -622.899 1265.798 1306.008 1281.703

-------
Details
-------
Convergence reached? = FALSE
Iterations = 300 / 300
Criteria = 0.00058
MC sample = 20
Cut point = 0.25
Processing time = 22.83885 mins
}
}

\keyword{SAEM}
\keyword{quantile}
\keyword{quantile regression}
\keyword{nonlinear mixed models}
\keyword{ALD}
