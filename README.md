# Rerr

Due to the limitation of Markdown engine on Github, equations and
outputs cannot be displayed. A full version of this document can be
viewed
[here](http://htmlpreview.github.io/?https://github.com/zhuozhang/Rerr/blob/master/vignettes/Rerr.html).

## Introduction

__Rerr__ is a package that fits excess relative risk (ERR) models in
survival analysis, while allowing for measurement error in exposures,
with or without effect modifications.  They are widely used in
Radiation epidemiology.

Currently this package only allow one of the exposure variables to be
measured with error.  This document will walk you through the package
from installation to analysis of a simple survival dataset, simulated
also through the package.

## Installation

The package can be downloaded from Github at
[here](https://github.com/zhuozhang/Rerr).  After unzipping the
downloaded file, assume the full path to the folder of the package is
**path_to_foler**.  Then type the following command in R console to
install the package (replace `path_to_folder` with the real path):

```{r eval=FALSE}
install.packages(path_to_folder, repos=NULL, type="source")
```

After installation, a full version of this document with outputs and
mathematical notations can also be found in `R` by running
```{r}
browseVignettes("Rerr")
```

## Quick Start

### Model specification

The package allows user to specify the ERR model using `formula`
object in R.  The simulation example provided in this document is
based on a hypothetical ERR model, formulated as follows (ommiting the
paramters) 

	h(t) = exp(1 + log^2(age/40) + sex)(1 + x1 exp(sex) + x2)

where `x1`, `x2` are some time-dependent or time-independent
exposure variables.  The corresponding `formula` specification is

```{r eval=FALSE}
case ~ I(log(age/60)^2) + sex + ERR( U(x1) | sex) + ERR(x2) + offset(pt)
```

As can be seen, 

*  variables associated with baseline risks are specified as in
   regular linear regressions;
*  the relative risk part associated with `x1`, `x2` are enclosed by
   `ERR()`;
*  dose effect modifiers, if any, are separated from exposure by `|`;
*  the exposure variable that is associated with measurement error is
   further enclosed by `U()`;

There are two more variables introduced in the model specification,
`case` and `offset(pt)`. `offset(pt)` specify that `pt` will
represents the amount of person-time in each period, and `case`
represents the outcome indicator variable.  These two variables are
standard variables used in Poisson table used to analyze survival
data.  For details, refer to Holford 1980, Laird and Olivier 1981.

### Data preparation

Now we simulate a cohort consists of 5,000 people.  For simplicity we
assume that each person is followed for 10 years.

```{r}
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

```

To generate a Monte Carlo dosimetry system (MCDS), we use the SUMA
model proposed by Stram and Kopecky 2003, with only multiplicative
errors.  We assume the variance is 0.2 for shared error and 0.5 for
unshared errors. Here we generated 400 dose replications, as defined
by `n.rep`.

```{r}
sm.var <- 0.2                # variance of shared multiplicative error
sm.log.se <- log(1 + sm.var) # standard error on log scale
sm.log.mu <- - 1/2 * sm.var  # mean on log scale
um.var <- 0.5                # variance of unshared multiplicative error
um.log.se <- log(1 + um.var) # standard error on log scale
um.log.mu <- - 1/2 * um.var  # mean on log scale
x1.reps <- do.call(cbind, lapply(1:n.rep, function(i) {
    x1 * exp(rnorm(1, sm.log.mu, sm.log.se)) *
        exp(rnorm(n * n.year, um.log.mu, um.log.se))}))
```

We use the first dose replication for `x1` as the true dose for
simulation of the survival data.

```{r}
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

```
Shown above is the structure of the `data.frame` object constructed
for the survival data simulation.  

### Simulation and Inference

We use function __simulateERR__ to generate the outcome data using the
specified ERR model.  We also specified the parameters for the model,
following the order of the variables in the `formula` object for the
model.

```{r}

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
```

We can see that there are `r number.case` cases. There is a `index`
column in the generated `data.frame` object, `data.new`.  This column
indexes the rows of the original `data` that are kept. The reason of
having such column is whenever an event occurs, the hypothetical
exposure in the 'future' should not exist. This column is necessary
for aligning the dose replications to the generated `data.frame`.

We fit the model using the average dose in the remaining dose
replications, and then calculate the confidence intervals of the
parameters.  We use three functions, __parseERR__, __fitERR__, and
__inferERR__.  

-  `parseERR` parses the data and the model into a ERR object, which
   contains essential information for fitting and making inference;
-  `fitERR` use the parsed ERR object to fit the ERR model to get the
   MLE;
-  `inferERR` takes the parsed ERR object and the MLE to calculate the
   confidence intervals.

```{r}
data.new$x1 <- x1.mean[data.new$index]
data.err <- parseERR(model, data.new)
xm <- fitERR(data.err)
xm$conv
xm$theta
theta.CI <- inferERR(data.err, xm$theta,
                     x1.reps[data.new$index, rep.index],
                     naive = TRUE)
theta.CI
```

This concludes a simple simulation. As we can see here, the naive
confidence interval is narrower for `x1` since it does not account for
the variance introduced by shared error component in `x1`.

