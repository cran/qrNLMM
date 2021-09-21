\name{group.plots}
\alias{group.plots}
\alias{group.plot}
\alias{group.lines}
\alias{group.points}

\title{Plot function for grouped data
}
\description{Functions for plotting a profiles plot for grouped data.
}
\usage{
group.plot(x,y,groups,...)
group.lines(x,y,groups,...)
group.points(x,y,groups,...)
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
  \item{...}{
additional graphical arguments passed to \code{\link[graphics]{matplot}}. See \code{\link[graphics]{par}}.
}
}
\author{
Christian E. Galarza <\email{chedgala@espol.edu.ec}> and 
Victor H. Lachos <\email{hlachos@uconn.edu}>
}
\seealso{\code{\link{Soybean}}, \code{\link{HIV}}, \code{\link{QRNLMM}},
\code{\link[qrLMM]{QRLMM}}
}

\examples{
\dontrun{
#A full profile plot for Soybean data

data(Soybean)
attach(Soybean)

group.plot(x = Time,y = weight,groups = Plot,type="b",
             main="Soybean profiles",xlab="time (days)",
             ylab="mean leaf weight (gr)")

#Profile plot by genotype

group.plot(x = Time[Variety=="P"],y = weight[Variety=="P"],
             groups = Plot[Variety=="P"],type="l",col="blue",
             main="Soybean profiles by genotype",xlab="time (days)",
             ylab="mean leaf weight (gr)")

group.lines(x = Time[Variety=="F"],y = weight[Variety=="F"],
              groups = Plot[Variety=="F"],col="black")
}
}

\keyword{plot grouped data}
\keyword{plot}
\keyword{nested data}
\keyword{mixed models}
