## generate piecewise exponential random variables
rpexp <- function (n, rate=1, t=0) {
    if (length(t) != length(rate))
        stop("Different length for rate and t")
    index <- rep(TRUE, n)
    result <- rep(0, n)
    for (i in seq_along(rate)) {
        result[index] <- result[index] + rexp(n, rate[i])[index]
        if (i == length(rate))
            break
        index[result <= t[i+1]] <- FALSE
        result[index] <- t[i+1]
        if (!any(index))
            break
    }
    return (result)
}

## Convert annual dose to cumulative dose with given id vector. The
## years associated with dose records for a given id should be in
## increasing order.
annualDoseToCumulativeDose <- function(dose, doseyr.id) {
    unique.id <- unique(doseyr.id)
    id.start <- match(unique.id, doseyr.id)
    id.end <- length(doseyr.id) + 1 - match(unique.id, rev(doseyr.id))
    if (is.matrix(dose)) {
        return(apply(dose, 2, function(x) {
            unlist(lapply(seq_along(unique.id), function(i){ cumsum(x[id.start[i]:id.end[i]]) }))
        }))
    } else {
        return(unlist(lapply(seq_along(unique.id), function(i){ cumsum(dose[id.start[i]:id.end[i]]) })))
    }
}

## Convert cumulative dose to annual dose with given id vector. The
## years associated with dose records for a given id should be in
## increasing order.
cumulativeDoseToAnnualDose <- function(dose, doseyr.id) {
    unique.id <- unique(doseyr.id)
    id.start <- match(unique.id, doseyr.id)
    id.end <- length(doseyr.id) + 1 - match(unique.id, rev(doseyr.id))
    if (is.matrix(dose)) {
        return(apply(dose, 2, function(x) {
            unlist(lapply(seq_along(unique.id), function(i){ c(x[id.start[i]], diff(x[id.start[i]:id.end[i]])) }))
        }))
    } else {
        return(unlist(lapply(seq_along(unique.id), function(i){ c(dose[id.start[i]], diff(dose[id.start[i]:id.end[i]])) })))
    }
}

## Generate dose realizations with annual dose records using
## Hierachical SUMA model (with only multiplicative errors). 
hsumaAnnualDose <- function(seed, uM, sM, pM, dose.mean, doseyr.id) {
    n <- length(dose.mean)
    sM.sigma2 <- log(1 + sM) # log standard error unshared error
    sM.mu <- -1/2 * sM.sigma2  # log mean
    uM.sigma2 <- log(1 + uM) # log standard error unshared error
    uM.mu <- -1/2 * uM.sigma2  # log mean
    pM.sigma2 <- log(1 + pM) # log standard error unshared error
    pM.mu <- -1/2 * pM.sigma2  # log mean
    unique.id <- unique(doseyr.id)
    id.start <- match(unique.id, doseyr.id)
    id.end <- length(doseyr.id) + 1 - match(unique.id, rev(doseyr.id))
    if (length(seed) > 1) {
        return(do.call(cbind, lapply(seed, function(seed.i) {
            set.seed(seed.i)
            sM.error <- rep(exp(rnorm(1, sM.mu, sqrt(sM.sigma2))), n)
            uM.error <- exp(rnorm(n, uM.mu, sqrt(uM.sigma2)))
            pM.error.unique <- exp(rnorm(length(unique.id), pM.mu, pM.sigma2))
            pM.error <- pM.error.unique[unlist(lapply(seq_along(unique.id),
                                                      function(i) rep(i, id.end[i]-id.start[i]+1)))]
            return(dose.mean * sM.error * uM.error * pM.error)
        })))
    } else {
        set.seed(seed)
        sM.error <- rep(exp(rnorm(1, sM.mu, sqrt(sM.sigma2))), n)
        uM.error <- exp(rnorm(n, uM.mu, sqrt(uM.sigma2)))
        pM.error.unique <- exp(rnorm(length(unique.id), pM.mu, pM.sigma2))
        pM.error <- pM.error.unique[unlist(lapply(seq_along(unique.id),
                                                  function(i) rep(i, id.end[i]-id.start[i]+1)))]
        return(dose.mean * sM.error * uM.error * pM.error)
    }
}

## Generate dose realizations with cumulative dose records using
## Hierachical SUMA model (with only multiplicative errors). 
hsumaCumulativeDose <- function(seed, uM, sM, pM, dose.mean, doseyr.id) {
    annualDoseToCumulativeDose(
        hsumaAnnualDose(seed, uM, sM, pM,
                        cumulativeDoseToAnnualDose(dose.mean, doseyr.id),
                        doseyr.id), doseyr.id)
}

## Wrapper function to generate cumulative dose records using
## Hierachical SUMA model
genHSUMA <- function(out.file, param, dose, id, n.rep, seed) {
    if (missing(seed)) seed <- runif(n.rep) * .Machine$integer.max
    hsuma.rep <- hsumaCumulativeDose(
        seed,
        param[1], # unshared eror
        param[2], # shared error
        param[3], # partially shared error
        dose,
        id)
    if (is.null(out.file)) {
        return(hsuma.rep)
    } else {
        saveRDS(hsuma.rep, out.file)
        rm(hsuma.rep)
        gc()
    }
}

