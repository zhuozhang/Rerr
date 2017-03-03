\name{Rerr-package}
\alias{Rerr-package}
\alias{Rerr}
\docType{package}
\title{
Excessive Relative Risk (ERR) model fitting and inference
}
\description{
A package that fits and make inference of a ERR model with adjustment
for measurement errors described by a Monte Carlo dosimetry system.
}
\details{
  
Rerr is a package that fits excess relative risk (ERR) models in
survival analysis, while allowing for measurement error in exposures,
with or without effect modifications.  Currently this package only allow
one of the exposure variables to be measured with error.

To learn to use this package, check the Rerr Vignette, or see the
  example provided below.

}
\author{
Zhuo Zhang, Daniel Stram
Maintainer: Zhuo Zhang 
}
\references{
  TODO
}
\keyword{ Rerr }
\seealso{
  \code{\link{simulateERR}}
  \code{\link{parseERR}}
  \code{\link{fitERR}}
  \code{\link{inferERR}}
}
\examples{

## a simple simulation study

library(Rerr)

n <- 5000                              # size of the cohort
n.year <- 10                           # years of follow-up
n.rep <- 400                           # number of dose realizations to generate
id <- rep(1:n, each = n.year)          # id for py table
start.age <- sample(20:50, n, TRUE)    # age of enrollment of each person
age <- unlist(lapply(start.age, function(x)
    x + 1:n.year - 1))                 # age for py table
sex <- rep(rbinom(n, 1, 0.5),
           each = n.year)              # sex for py table
time <- rep(1:n.year - 1, n)           # time variable for py table
censor.time <- rep(n.year, n * n.year) # censor time for py table
x1 <- rchisq(n * n.year, 1)            # exposure 1 for py table
x2 <- rchisq(n * n.year, 1)            # exposure 2 for py table

sm.var <- 0.2                          # variance of shared multiplicative error
sm.log.se <- log(1 + sm.var)           # standard error on log scale
sm.log.mu <- - 1/2 * sm.var            # mean on log scale
um.var <- 0.5                          # variance of unshared multiplicative error
um.log.se <- log(1 + um.var)           # standard error on log scale
um.log.mu <- - 1/2 * um.var            # mean on log scale
x1.reps <- do.call(cbind, lapply(1:n.rep, function(i) {
    x1 * exp(rnorm(1, sm.log.mu, sm.log.se)) *
        exp(rnorm(n * n.year, um.log.mu, um.log.se))}))

rep.index <- 1:n.rep
rep.index <- rep.index[-1]
x1.mean <- apply(x1.reps[, rep.index], 1, mean)
data <- data.frame(sex = sex,
                   age = age,
                   x1 = x1.reps[, 1],
                   x2 = x2,
                   id = id,
                   time = time,
                   censor.time = censor.time)
data[1:10, ]

model <- case ~ I(log(age/40)^2) + sex +
    ERR(U(x1) | sex) + ERR(x2) +
    offset(py)
theta <- c(-5, log(2), 1, 2, 0.5, 1)

data.new <- simulateERR(model, data, theta,
                        id = id,
                        time = time,
                        censor.time = censor.time)
number.case <- sum(data.new$case)

number.case

data.new[1:10,]

data.new$x1 <- x1.mean[data.new$index]
data.err <- parseERR(model, data.new)
xm <- fitERR(data.err)
xm$conv
xm$theta
theta.CI <- inferERR(data.err, xm$theta,
                     x1.reps[data.new$index, rep.index],
                     naive = TRUE)
theta.CI
}
