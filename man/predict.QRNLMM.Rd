\name{predict.QRNLMM}
\alias{predict.QRNLMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict method for Quantile Regression for Nonlinear Mixed-Effects (QRNLMM) fits
}

\description{
Takes a fitted object produced by \code{QRNLMM()} and produces predictions given a new set of values for the model covariates.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{

\method{predict}{QRNLMM}(object,x = NULL,groups = NULL,covar = NULL, y = NULL,MC = 1000,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted QRNLMM object as produced by QRNLMM().
}
  \item{x}{
vector of longitudinal (repeated measures) covariate of dimension \eqn{N}. For example: Time, location, etc.  
}
  \item{groups}{
factor of dimension \eqn{N} specifying the partitions of the data over which the random effects vary.
}
  \item{covar}{
a matrix of dimension \eqn{N \times r} where \eqn{r} represents the number of covariates.
}
  \item{y}{
the response vector of dimension \eqn{N} where \eqn{N} is the total of observations. Optional. See details.
}
\item{MC}{
number of MC replicates for the computation of new individual values (only when \code{y} is provided). By default \code{MC = 1000}. See details.
}

\item{\dots}{
additional arguments affecting the predictions produced.
}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A \code{data.frame} containing the predicted values, one column per quantile.
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\details{

Prediction for QRNLMM objects can be performed under three different scenarios:
\enumerate{
\item \code{predict(object)}: if no newdata is provided, fitted values for the original dataset is returned. Please refer to the \code{fitted.values} value in \code{\link{QRNLMM}}.

\item \code{predict(object,x,groups,covar = NULL)}: if new data is provided, but only the independent variables (no response), population curves are provided. If no covariates are provided, the predicted curves will be the same.

\item \code{predict(object,x,groups,covar = NULL,y)} if the response values are provided, a Metropolis-Hastings algorithm (with \code{MC} replicates and \code{thin = 5}) is performed in order to compute the random-effects for new subjects. The method is based on Galarza et.al. (2020).
}

}


\note{
For scenario 3, results may vary a little each time. For more precision, please increase \code{MC}.
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

\seealso{
\code{\link{QRNLMM}},\code{\link{group.plots}},\code{\link{group.lines}},\code{\link{Soybean}}, \code{\link{HIV}}, \code{\link[lqr]{lqr}}
}

\examples{
  \dontrun{
    
#A model for comparing the two genotypes  (with covariates)

data(Soybean)
attach(Soybean)

y     = weight          #response
x     = Time            #time
covar = c(Variety)-1    #factor genotype (0=Forrest, 1=Plan Intro)

#Expression for the three parameter logistic curve with a covariate

exprNL = expression((fixed[1]+(fixed[4]*covar[1])+random[1])/
                      (1 + exp(((fixed[2]+random[2])- x)/
                      (fixed[3]+random[3]))))

#Initial values for fixed effects
initial = c(max(y),0.6*max(y),0.73*max(y),3)

# A quantile regression for percentiles p = c(0.05,0.50,0.95)

#Take your sit and some popcorn

results = qrNLMM::QRNLMM(
  y = y,
  x = x,
  groups = Plot,
  initial = initial,
  exprNL = exprNL,
  covar = covar,
  p= c(0.05,0.50,0.95),# quantiles to estimate
  MaxIter = 50,M = 15, # to accelerate
  verbose = FALSE      # show no output
)

######################################################################
# Predicting
######################################################################

# now we select two random subjects from original data
set.seed(19)
index = Plot \%in\% sample(Plot,size = 2)
index2 = c(1,diff(as.numeric(Plot)))>0 & index

# 1. Original dataset
#####################

prediction = predict(object = results)
head(prediction) # fitted values

if(TRUE){
  group.plot(x = Time[index],
             y = weight[index],
             groups = Plot[index],
             type="b",
             main="Soybean profiles",
             xlab="time (days)",
             ylab="mean leaf weight (gr)",
             col= ifelse(covar[index2],"gray70","gray90"),
             ylim = range(prediction[index,]),
             lty = 1
  )
  
  legend("bottomright",
         legend = c("Forrest","Plan Intro"),
         bty = "n",col = c(4,2),lty = 1)
  
  # predictions for these two plots
  
  group.lines(x = Time[index], # percentile 5
              y = prediction[index,1],
              groups = Plot[index],
              type = "l",
              col=ifelse(covar[index2],"red","blue"),
              lty = 2
  )
  
  group.lines(x = Time[index], # median
              y = prediction[index,2],
              groups = Plot[index],
              type = "l",
              col="black",
              lty = 2)
  
  group.lines(x = Time[index], # percentile 95
              y = prediction[index,3],
              groups = Plot[index],
              type = "l",
              col=ifelse(covar[index2],"red","blue"),
              lty = 2)
  
  legend("topleft",
         legend = c("p = 5","p = 50","p = 95"),
         col = c(4,1,4),lty = c(2,2,2),bty = "n")
}

# 2. New covariates with no response
####################################

# For the two randomly selected plots (index == TRUE)

# newdata
newdata = data.frame(new.groups = Plot[index],
                new.x = Time[index],
                new.covar = covar[index])

newdata
attach(newdata)

prediction2 = predict(object = results,
                      groups = new.groups,
                      x = new.x,
                      covar = new.covar)

# population curves

if(TRUE){
  group.plot(x = new.x, # percentile 5
             y = prediction2[,1],
             groups = new.groups,
             type = "l",
             col=ifelse(covar[index2],"red","blue"), 
             lty = 2,
             main="Soybean profiles",
             xlab="time (days)",
             ylab="mean leaf weight (gr)",
             ylim = range(prediction2)
  )
  
  legend("bottomright",
         legend = c("Forrest","Plan Intro"),
         bty = "n",col = c(4,2),lty = 1)
  
  # predictions for these two plots
  
  group.lines(x = new.x, # median
              y = prediction2[,2],
              groups = new.groups,
              type = "l",
              col="black",
              lty = 1)
  
  group.lines(x = new.x, # percentile 95
              y = prediction2[,3],
              groups = new.groups,
              type = "l",
              col=ifelse(covar[index2],"red","blue"),
              lty = 2)
  
  legend("topleft",
         legend = c("p = 5","p = 50","p = 95"),
         col = c(4,1,4),lty = c(2,1,2),bty = "n")
  
  segments(x0 = new.x[new.covar==1],
           y0 = prediction2[new.covar==1,1],
           y1 = prediction2[new.covar==1,3],
           lty=2,col = "red")
  
  segments(x0 = new.x[new.covar==0],
           y0 = prediction2[new.covar==0,1],
           y1 = prediction2[new.covar==0,3],
           lty=2,col = "blue")
}

# 3. New covariates + response
####################################

# newdata
newdata2 = data.frame(new.groups = Plot[index],
                     new.x = Time[index],
                     new.covar = covar[index],
                     new.y = weight[index])

newdata2
attach(newdata2)

prediction2 = predict(object = results,
                      groups = new.groups,
                      x = new.x,
                      covar = new.covar,
                      y = new.y)

# individual curves (random-effects to be computed)

if(TRUE){
  group.plot(x = Time[index],
             y = weight[index],
             groups = Plot[index],
             type="b",
             main="Soybean profiles",
             xlab="time (days)",
             ylab="mean leaf weight (gr)",
             col= ifelse(covar[index2],"gray70","gray90"),
             ylim = range(prediction[index,]),
             lty = 1
  )
  
  legend("bottomright",
         legend = c("Forrest","Plan Intro"),
         bty = "n",col = c(4,2),lty = 1)
  
  # predictions for these two plots
  
  group.lines(x = new.x, # percentile 5
              y = prediction2[,1],
              groups = new.groups,
              type = "l",
              col=ifelse(covar[index2],"red","blue"),
              lty = 2)
  
  group.lines(x = new.x, # median
              y = prediction2[,2],
              groups = new.groups,
              type = "l",
              col="black",
              lty = 1)
  
  group.lines(x = new.x, # percentile 95
              y = prediction2[,3],
              groups = new.groups,
              type = "l",
              col=ifelse(covar[index2],"red","blue"),
              lty = 2)
  
  legend("topleft",
         legend = c("p = 5","p = 50","p = 95"),
         col = c(4,1,4),lty = c(2,1,2),bty = "n")
}
  }
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SAEM}
\keyword{prediction}
\keyword{quantile}
\keyword{mixed-effects}
\keyword{quantile regression}
\keyword{nonlinear mixed models}
