\name{ERR}
\alias{err}
\alias{parseERR}
%% \alias{msiERR}
\alias{fitERR}
%% \alias{varERR}
\alias{inferERR}
\alias{simulateERR}

\title{
Functions related to Excessive Relative Risk (ERR) model fitting and inference
}

\description{
  \name{parseERR} parses a ERR model with data, and gives a structure
  that contains information for fitting the model and making inference.
  
  \name{fitERR} fits a ERR model based on the parsed structure returned
  by \name{parseERR}.

  \name{inferERR} makes inference about a ERR model based on the parsed
  structure returned by \name{parseERR}.

  \name{simulateERR} simulates survival data based on ERR model.
}

\usage{
  parseERR(formula, data, env=parent.frame())

  fitERR(err.struct, theta.start, n.iter=50, epsilon=1e-4, step.depth=5,
  index)

  inferERR(err.struct, theta, dose.reps, alpha=0.05, score.approx=FALSE,
  naive=FALSE)

  simulateERR(formula, data, theta, id, time, censor.time, seed)
}

\arguments{

  parseERR()
  
  \item{formula}{The formula object will be used to describe the model.
  Excessive Relative Risk terms will be enclosed by 'ERR(\cdot)'. Inside
  'ERR(\cdot)', the dose effect and effect modifiers are separated by
  '|', with dose effect on the left side. Error-prone dose effect will
  be enclosed by 'U(\cdot)'.}

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
  TODO
}
\references{
  TODO
}
\keyword{parseERR, fitERR, inferERR, ERR}
\seealso{
  TODO
%% Optional links to other man pages, e.g.
%% \code{\link[<pkg>:<pkg>-package]{<pkg>}} 
}
\examples{
  %% TODO

model.gen <- cell.cases ~ covariate.age + I(covariate.age^2) + sex +
    offset(cell.pyr) + ERR(U(dose.internal.true) | covariate.age + sex)
    + ERR(dose.external)

model.est <- cell.cases ~ covariate.age + I(covariate.age^2) + sex +
    offset(cell.pyr) + ERR(U(dose.internal) | covariate.age + sex) +
    ERR(dose.external)

data <- simulateERR(model.gen, data.prep, theta, id, age, censor.time, seed)

ERR <- parseERR(model.est, data=data)

xm <- fitERR(ERR)

CI <- inferERR(ERR, xm$theta, dose.reps, score.approx=TRUE, naive=TRUE)

}
