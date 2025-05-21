\name{Soybean}
\alias{Soybean}
\title{Growth of soybean plants}
\description{
  The \code{Soybean} data frame has 412 rows and 5 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Plot}{
      a factor giving a unique identifier for each plot.
    }
    \item{Variety}{
      a factor indicating the variety; Forrest (F) or Plant
      Introduction #416937 (P).
    }
    \item{Year}{
      a factor indicating the year the plot was planted.
    }
    \item{Time}{
      a numeric vector giving the time the sample was taken (days
      after planting).
    }
    \item{weight}{
      a numeric vector giving the average leaf weight per plant (g).
    }
  }
}
\details{
  These data are described in
  Davidian and Giltinan (1995, 1.1.3, p.7) as ``Data from an
  experiment to compare growth patterns of two genotypes of soybeans:
  Plant Introduction #416937 (P), an experimental strain, and Forrest
  (F), a commercial variety.'' In order to fit the Nonlinear data we sugest to use the three parameter logistic model as in Pinheiro & Bates (1995).
}
\source{
  Pinheiro, J. C. and Bates, D. M. (2000), \emph{Mixed-Effects Models in S
    and S-PLUS}, Springer, New York.  (Appendix A.27)

  Davidian, M. and Giltinan, D. M. (1995), \emph{Nonlinear Models for
    Repeated Measurement Data},  Chapman and Hall, London.
}
\examples{
\dontrun{
data(Soybean)
attach(Soybean)

#################################
#A full model (no covariate)

y     = weight          #response
x     = Time            #time

#Expression for the three parameter logistic curve

exprNL = expression((fixed[1]+random[1])/(1 + exp(((fixed[2]+random[2])- x)/(fixed[3]+random[3]))))

#Initial values for fixed effects
initial = c(max(y),0.6*max(y),0.73*max(y))

#A median regression (by default)
median_reg = QRNLMM(y,x,Plot,initial,exprNL)

#Assing the fit

fxd     = median_reg$res$beta
nlmodel = median_reg$res$nlmodel
seqc    = seq(min(x),max(x),length.out = 500)

group.plot(x = Time,y = weight,groups = Plot,type="l",
           main="Soybean profiles",xlab="time (days)",
           ylab="mean leaf weight (gr)",col="gray")

lines(seqc,nlmodel(x = seqc,fixed = fxd,random = rep(0,3)),
      lwd=2,col="blue")             

#Histogram for residuals
hist(median_reg$res$residuals,breaks = 20)

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

# A quantile regression for the three quartiles ()
box_reg = QRNLMM(y,x,Plot,initial,exprNL,covar,p=c(0.25,0.50,0.75),MaxIter = 20)

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


}
}
\keyword{datasets}
