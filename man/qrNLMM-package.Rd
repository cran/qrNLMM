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
\author{Christian E. Galarza <\email{chedgala@espol.edu.ec}> and 
Victor H. Lachos <\email{hlachos@ime.unicamp.br}>

Maintainer: Christian E. Galarza <\email{chedgala@espol.edu.ec}>
}
\references{
Galarza, C.E., Castro, L.M., Louzada, F. & Lachos, V. (2020) Quantile regression for nonlinear mixed effects models: a likelihood based perspective. Stat Papers 61, 1281-1307.
\doi{10.1007/s00362-018-0988-y}

Yu, K. & Moyeed, R. (2001). Bayesian quantile regression. Statistics & Probability Letters, 54(4), 437-447.

Yu, K., & Zhang, J. (2005). A three-parameter asymmetric Laplace distribution and its extension. Communications in Statistics-Theory and Methods, 34(9-10), 1867-1879.
}

\keyword{package}
\keyword{quantile}
\keyword{ALD}
\seealso{\code{\link{Soybean}}, \code{\link{HIV}}, \code{\link{QRNLMM}},
\code{\link[lqr]{lqr}} , \code{\link{group.plots}}
}
\examples{
#See examples for the QRNLMM function linked above.
}
