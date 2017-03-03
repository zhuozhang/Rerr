\name{Dosimetry}
\alias{Dose}
\alias{annualDoseToCumulativeDose}
\alias{cumulativeDoseToAnnualDose}
\alias{hsumaAnnualDose}
\alias{hsumaCumulativeDose}

\title{
Functions related to manipulating yearly dose data
}

\description{

  \code{annualDoseToCumulativeDose} Convert annual dose to cumulative
dose with given id vector. The years associated with dose records for a
given id should be in increasing order.
  
  \code{cumulativeDoseToAnnualDose} Convert cumulative dose to annual
dose with given id vector. The years associated with dose records for a
given id should be in increasing order.

  \code{hsumaAnnualDose} Generate dose realizations with annual dose
 records using hierachical SUMA model (with only multiplicative errors).


  \code{hsumaCumulativeDose} Generate dose realizations with cumulative
dose records using hierachical SUMA model (with only multiplicative
errors).

}

\usage{

annualDoseToCumulativeDose(dose, doseyr.id)

cumulativeDoseToAnnualDose(dose, doseyr.id)

hsumaAnnualDose(seed, uM, sM, pM, dose.mean, doseyr.id)

hsumaCumulativeDose(seed, uM, sM, pM, dose.mean, doseyr.id)

}

\arguments{
\item{dose}{A vector describing year exposure, either annually or
cumulative}

\item{doseyr.id}{A vector describing the id information of the dose
vector.  Both dose and doseyr.id should be in increasing order with
respect to doseyr.id}

\item{uM, sM, pM}{Multiplicative error, unshared, shared and partially
shared}

\item{dose.mean}{Average dose vecotr}

\item{seed}{Random seed to ensure reproducibility.}
}

\details{
Function \code{cumulativeDoseToAnnualDose} and
\code{annualDoseToCumulativeDose} help with converting cumulative data to
annual data and the other way around.

Function \code{hsumaAnnualDose} and \code{hsumaCumulativeDose} generates
annual and cumulative dose replications using the SUMA model.
}

\references{
}

\keyword{SUMA}

\seealso{
}

\examples{

n <- 5000
n.year <- 10
id <- rep(1:n, each = n.year)
z.annual <- rchisq(n * n.year, 1)
z.cumulative <- annualDoseToCumulativeDose(z.annual, id)

z.annual.back <- cumulativeDoseToAnnualDose(z.cumulative, id)

all.equal(z.annual, z.annual.back)

um <- 0.5
sm <- 0.2

z.reps.annual <- hsumaAnnualDose(0, um, sm, 0, z.annual, id)
z.reps.cumulative <- hsumaCumulativeDose(0, um, sm, 0, z.cumulative, id)

}