\name{qrNLMM-package}
\alias{qrNLMM-package}
\docType{package}
\title{
Package for Quantile Regression for Linear Mixed-Effects Model
}
\description{
This package contains a principal function that performs a quantile regression for a Nonlinear Mixed-Effects Model using the Stochastic-Approximation of the EM Algorithm (SAEM) for an unique or a set of quantiles.

Exploiting the nice hierarchical representation of the ALD, our classical approach follows the Stochastic Approximation of the EM(SAEM) algorithm for deriving exact maximum likelihood estimates of the fixed-effects and variance components.
}
\details{
\tabular{ll}{
Package: \tab qrNLMM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-01-30\cr
License: \tab What license is it under?\cr
}
}
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and 
Victor H. Lachos <\email{hlachos@ime.unicamp.br}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\references{
Galarza, C.M., Castro, LM. Louzada, F. and V. H. Lachos (2016) Quantile Regression for Nonlinear Mixed Effects Models: A Likelihood Based Perspective. Technical Report 8, Universidade Estadual de Campinas. 
\url{http://www.ime.unicamp.br/sites/default/files/rp08-2016.pdf}

Yu, K. & Moyeed, R. (2001). Bayesian quantile regression. Statistics & Probability Letters, 54(4), 437-447.

Yu, K., & Zhang, J. (2005). A three-parameter asymmetric Laplace distribution and its extension. Communications in Statistics-Theory and Methods, 34(9-10), 1867-1879.
}

\keyword{package}
\keyword{quantile}
\keyword{ALD}
\seealso{\code{\link{Soybean}}, \code{\link{HIV}}, \code{\link{QRNLMM}},
\code{\link[lqr]{lqr}}, \code{\link[qrLMM]{QRLMM}}, \code{\link{group.plots}}
}
\examples{
#See examples for the QRNLMM function linked above.
}