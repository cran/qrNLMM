\name{HIV}
\alias{HIV}
\docType{data}
\title{
HIV viral load study
}
\description{
The data set belongs to a clinical trial (ACTG 315) studied in previous researches by Wu (2002) and Lachos et al. (2013). In this study, we analyze the HIV viral load of 46 HIV-1 infected patients under antiretroviral treatment (protease inhibitor and reverse transcriptase inhibitor drugs). The viral load and some other covariates were mesured several times days after the start of treatment
been 4 and 10 the minimum and maximum number of measures per patient respectively.
}

\usage{data(HIV)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{patid}}{a numeric vector indicating the patient register number.}
      \item{\code{ind}}{a numeric vector indicating the number patient on which the
      measurement was made. It represents the subject number in the study.}
    \item{\code{day}}{time in days.}
    \item{\code{cd4}}{cd4 count in cells/\eqn{mm^{3}}.}
    \item{\code{lgviral}}{viral load in log10 scale.}
    \item{\code{cd8}}{cd8 count in cells/\eqn{mm^{3}}.}
  }
}
\details{
In order to fit the nonlinear data we sugest to use the Nonlinear model proposed by Wu (2002) and also used by Lachos et al. (2013).}

\source{
Wu, L. (2002). A joint model for nonlinear mixed-effects models with censoring and covariates measured with error, with application to aids studies. Journal of the American Statistical association,
97(460), 955-964.

Lachos, V. H., Castro, L. M. & Dey, D. K. (2013). Bayesian inference in nonlinear mixed-effects
models using normal independent distributions. Computational Statistics & Data Analysis, 64,
237-252.
}

\keyword{datasets}