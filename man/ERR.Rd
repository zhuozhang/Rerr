\name{ERR}
\alias{err}
\alias{simulateERR}
\alias{parseERR}
%% \alias{msiERR}
\alias{fitERR}
%% \alias{varERR}
\alias{inferERR}

\title{
Functions related to Excessive Relative Risk (ERR) model fitting and inference
}

\description{
  \code{simulateERR} simulates survival data based on ERR model.

  \code{parseERR} parses a ERR model with data, and gives a structure
  that contains information for fitting the model and making inference.
  
  \code{fitERR} fits a ERR model based on the parsed structure returned
  by \code{parseERR}.

  \code{inferERR} makes inference about a ERR model based on the parsed
  structure returned by \code{parseERR}.
}

\usage{
  simulateERR(formula, data, theta, id, time, censor.time, seed)

  parseERR(formula, data, env=parent.frame())

  fitERR(err.struct, theta.start, n.iter=50, epsilon=1e-4, step.depth=5,
  index)

  inferERR(err.struct, theta, dose.reps, alpha=0.05, score.approx=FALSE,
  naive=FALSE)
}

\arguments{

  parseERR()
  
  \item{formula}{The formula object will be used to describe the model.
  Excessive Relative Risk terms will be enclosed by 'ERR()'. Inside
  'ERR()', the dose effect and effect modifiers are separated by
  '|', with dose effect on the left side. Error-prone dose effect will
  be enclosed by 'U()'.}

  \item{data}{data.frame that contains needed data in the formula}

  \item{env}{environment where the formula will be evaluated}

  fitERR()
  
  \item{err.struct}{ERR model structure returned by parseERR}

  \item{theta.start}{parameter values to start with for fitting}

  \item{n.iter}{maximum number of iterations}

  \item{epsilon}{criteria for determining convergence}

  \item{step.depth}{number of steps to reduce the step size during fitting}

  \item{index}{index of boolean vector to determine which parameters are
  free to be estimated. All the other parameters will be kept fixed}

  inferERR()
  
  \item{theta}{estimated parameter values used for constructing
  confidence interval} 

  \item{dose.reps}{Dose realizations for error-prone covariates. For
  multiple error-prone covariates, a list should be provided}

  \item{alpha}{Significance level for constructing confidence intervals}

  \item{score.approx}{Use approximation of variance of parameters as in
    score test to construct confidence interval}

  \item{naive}{Give naive confidence interval}

  simulateERR()

  \item{id}{id information to identify individuals in the data.frame}

  \item{time}{time points for each record in the data.frame}

  \item{censor.time}{censor time for each individual. It can have the
  same length as 'id', or it can only contain censor time for each
  unique 'id'}

  \item{seed}{random seed for generating the survival data}
}


\details{
  These 4 functions are used for simulating survival data, fitting
  survival data, and making inference of model parameters, all based on
  excess relative risk models (ERR).

  simulateERR takes a data.frame containing survival data and a formula
  object that specifying the hazard function to calculate the event
  time. The returned data.frame object contains outcome information
  formulated as in the original data.frame in a person-time
  fashion.  Person-time spent in each period is also provided.  One
  extra column keeps track of the follow-up index.

  parseERR takes the survival data in tabulated Poisson format and a ERR
  model formula, and parse the information to a ERR object for further
  analysis.

  fitERR uses Fisher's scoring to calculate the maximum likelhood
  estimator with the ERR object returned by parseERR, and returns the
  MLE and convergence information.

  inferERR uses the ERR object returned by parseERR and estimates of the
  paramters to calculate confidence intervals of these estimates.
  
}
\references{
  For construction of Poisson format data for survival analysis, refer
  to Laird and Oliver 1981.
}

\keyword{parseERR, fitERR, inferERR, simulateERR, ERR}
%% \seealso{
%%   TODO
%% Optional links to other man pages, e.g.
%% \code{\link[<pkg>:<pkg>-package]{<pkg>}} 
%% }

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

sm.var <- 0.2                # variance of shared multiplicative error
sm.log.se <- log(1 + sm.var) # standard error on log scale
sm.log.mu <- - 1/2 * sm.var  # mean on log scale
um.var <- 0.5                # variance of unshared multiplicative error
um.log.se <- log(1 + um.var) # standard error on log scale
um.log.mu <- - 1/2 * um.var  # mean on log scale
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
