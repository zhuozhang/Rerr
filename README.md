# Rerr

## Installation

Rerr package can be installed using [install.packages](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/install.packages.html) in R.
Download the zip file, and unzip.  Then use the following command to install.
```r
install.packages(path_to_folder, repos=NULL, type="source")
```

## Example Code

```r
library(Rerr)
## define cohort parameter and generate covariates
n <- 1000
n.years <- 5
theta <- c(-6, 1, 2)
id <- rep(1:n, each=n.years)
time <- rep(1:n.years, n)
censor.time <- rep(n.years, n * n.years)
sex <- rbinom(n, 1, 0.5)
Z <- rchisq(n * n.years, 1)

## generate dose replications
n.reps <- 400
sM <- 0.1 # variance of shared multiplicative error
uM <- 0.1 # variance of unshared multiplicative error
sM.sigma2 <- log(1 + sM)
sM.mu <- -1/2 * sM.sigma2
uM.sigma2 <- log(1 + uM)
uM.mu <- -1/2 * uM.sigma2
X.reps <- do.call(cbind, lapply(1:n.reps, function(i) {
          sM.error <- exp(rnorm(1, sM.mu, sqrt(sM.sigma2)))
          uM.error <- exp(rnorm(n, uM.mu, sqrt(uM.sigma2)))
          return(Z * sM.error * uM.error)
    	  }))

## define ERR models
model.true <- cell.case ~ sex + ERR(U(dose)) + offset(cell.pyr)
model <- cell.case ~ sex + ERR(U(dose.mean)) + offset(cell.pyr)

dose <- X.reps[,1]
dose.reps <- X.reps[, -1]
dose.mean <- apply(dose.reps, 1, mean)
data <- data.frame(sex=sex, dose=dose, id=id, time=time, censor.time=censor.time, dose.mean=dose.mean)

## simulate survival data into tabulated Poisson format
data.survival <- simulateERR(model.true, data, theta, id, time, censor.time)

## parameter estimation and inference
err <- parseERR(model, data.survival)
theta.hat <- fitERR(err)$theta
CI <- inferERR(err, theta.hat, dose.reps[data.survival$index,], naive=TRUE)
```
