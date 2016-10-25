## Given the arithmetic mean and standard deviation, calculate the
## mean and standard deviation of a log normal random variable in
## log scale.
##
logNormal2Normal <- function(nMean, nSD){
    lNvar <- log(1 + nSD^2/nMean^2) # log variance
    lNmean <- log(nMean) - (1/2)*lNvar # log mean
    return(c(lNmean, sqrt(lNvar)))
}

## Extract function components in formula 'x' given by argument
## 'name'. Return either the index of the component or a list
## seperating the desired components and other components
extractFunc <- function(x, name, returnIndex=FALSE) {
    x.terms <- terms(x)
    if (attr(x.terms, "response") == 1)
        l <- as.list(attr(x.terms, "variables"))[-c(1:2)]
    else
        l <- as.list(attr(x.terms, "variables"))[-1]
    index <- grep(paste0(name,"\\(.*\\)"), l)
    if (returnIndex)
        return(index)
    else 
        return (list(l[index], l[setdiff(1:length(l), index)]))
}

## Strip the function name of the component
stripFunc <- function(x, name) {
    if (class(x) == "call") {
        if (x[[1]] == name) {
            return (x[[2]])
        }
        for (i in seq_along(x)) {
            x[[i]] <- stripFunc(x[[i]], name)
        }
        return (x)
    }
    return (x)
}

## Given a call, extract the 'U' component of it
processU <- function(x) {
    x.formula <- eval(call("~", call("+", 0, x)))
    u.list <- extractFunc(x.formula, "U")
    u.list.index <- extractFunc(x.formula, "U", TRUE)
    invisible(lapply(u.list[[1]], function(i) { # sanity check
        if (length(i) == 1) {
            stop("valid input for U() is needed")
        }
        i <- i[[2]]
        if (length(i) > 1) {
            stop("no operations should be performed inside U()")
        }
    }))
    return (list(stripFunc(x, quote(U)), u.list.index))
}

## Transform a call of summed components into a list of components
sumCallToList <- function(x) {
    if (is.call(x) && x[[1]]==quote(`+`)) {
        return (c(sumCallToList(x[[2]]), sumCallToList(x[[3]])))
    }
    return (as.character(x))
}

## Return the length of the call of summed components
sumCallLength <- function(x) {
    if (length(x) == 0) return (0)
    if (class(x) == "name") return (1)
    if (class(x) == "call") return (sumCallLength(x[[2]]) + sumCallLength(x[[3]]))
}

## Extract ERR components of a formula, and return a parsed structure
## representing the ERR model
extractERR <- function(x) {
    if (class(x) != "formula") {
        stop("the input needs to be a formula")
    }
    err.call.list <- extractFunc(x, "ERR")[[1]]
    err.list <- lapply(err.call.list, function(i) {
        if (length(i) == 1) {
            stop("valid input for ERR() is needed")
        }
        i <- i[[2]]
        if (class(i) == "call" && i[[1]] == quote(`|`)) # we have effect modifier here
            result <- list(i[[2]], i[[3]])
        else ## no effect modifier
            result <- list(i, list())
        ## check if there is any operator on the dose 
    })
    err.list <- lapply(err.list, function(i) {
        list(processU(i[[1]]), i[[2]])
    })
    ## sanity check: no variable is used twice for dose effect
    dose.union <- do.call(c, lapply(err.list, function(i) sumCallToList(i[[1]][[1]])))
    dose.duplicated <- unique(dose.union[duplicated(dose.union)])
    if (length(dose.duplicated) > 0) {
        stop("dose effects are duplicated: ", paste(dose.duplicated, collapse=", "))
    }
    ## now we will construct design matrix based on this
    return (err.list)
}

## Parse err model with data and return an intermediate structure to
## be used in further fitting procedure
parseERR <- function(formula, data, env=parent.frame()){
    call.env <- env
    base.list <- extractFunc(formula, "ERR")[[2]]
    base.model <- paste0("~", paste(as.character(base.list), collapse="+"))
    formula.terms <- terms(formula)
    base.n <- length(base.list) - any(grepl("^offset\\(", base.list))
    if (!attr(formula.terms, "intercept"))
        base.model <- paste0(base.model, "-1")
    else
        base.n <- base.n + 1
    if (length(formula) == 3)
        base.model <- paste0(formula[[2]], base.model)
    base.model <- as.formula(base.model, env=call.env)
    err.list <- extractERR(formula)
    theta.n.list <- c(base.n, unlist(lapply(err.list, function(i)
        sumCallLength(i[[1]][[1]]) + sumCallLength(i[[2]]))))
    theta.n <- sum(theta.n.list)
    theta.list <- lapply(seq_along(theta.n.list), (function() {
        tmp <- c(0, cumsum(theta.n.list))
        function(i) {
            (tmp[i]+1):tmp[i+1]
        }
    })())
    U.index <- lapply(seq_along(err.list), function(i)
        err.list[[i]][[1]][[2]]+cumsum(theta.n.list)[i])
    ## construct design matrix
    mf.call <- match.call()
    mf.base.call <- mf.call
    mf.base.call[[2]] <- base.model
    mf.base.call[[1]] <- quote(stats::model.frame)
    mf.base <- eval(mf.base.call, call.env)
    mt.base <- attr(mf.base, "terms")
    y <- model.response(mf.base)
    offset <- model.offset(mf.base)
    mat.base <- model.matrix(mt.base, mf.base)
    n <- nrow(mat.base)
    err.mat.list <- rapply(err.list, function(x) {
        mf.err <- mf.call
        mf.err[[2]] <- as.formula(paste0("~-1+", deparse(x)), call.env)
        mf.err[[1]] <- quote(stats::model.matrix)
        mf.err$object <- mf.err$formula
        mf.err$formula <- NULL
        mf.err <- eval(mf.err, call.env)
        return (mf.err)
    }, c("name", "call"), how="replace")
    theta.name <- c(colnames(mat.base), unlist(lapply(err.mat.list, function(x) {
        name.dose <- colnames(x[[1]][[1]])
        name.modifier <- colnames(x[[2]])
        name.dose[x[[1]][[2]]] <- paste0("U(", name.dose[x[[1]][[2]]], ")")
        if (is.null(name.modifier))
            return (name.dose)
        else 
            return (c(name.dose, paste0(paste(name.dose, collapse="+"), "|", name.modifier)))
    })))
    return (list(call=mf.call, y=y, offset=offset,
                 err.list=err.list, U.index=U.index,
                 model.base=mat.base, model.err=err.mat.list,
                 theta.n=theta.n, theta.list=theta.list,
                 theta.name=theta.name))
}

## Helper function to calculate expected value, score, and information
## based on parsed ERR structure and theta, the parameter values.
msiERR <- function(err.struct, theta) {
    return(with(err.struct, {
        mu <- exp(model.base %*% theta[theta.list[[1]]])
        if (!is.null(offset)) mu <- mu * offset
        Q.mat <- matrix(NA, nrow(model.base), theta.n)
        G.mat.list <- list()
        dose.effect.list <- list()
        for (i in seq_along(err.list)) {
            dose.call <- err.list[[i]]
            dose.pos <- theta.list[[i+1]][1] - 1 +
                1:sumCallLength(dose.call[[1]][[1]])
            dose.theta <- theta[dose.pos]
            dose.mat <- model.err[[i]][[1]][[1]]
            dose.effect <- dose.mat %*% dose.theta
            dose.effect.list[[i]] <- dose.effect
            G.mat.list[[i]] <- c(mu)
            Q.mat[, dose.pos] <- dose.mat
            if (length(dose.call[[2]])) { # with dose effect modifiers
                modifier.pos <- dose.pos[length(dose.pos)] +
                    1:sumCallLength(dose.call[[2]])
                modifier.theta <- theta[modifier.pos]
                modifier.mat <- model.err[[i]][[2]]
                modifier.effect <- modifier.mat %*% modifier.theta
                dose.effect.list[[i]] <-
                    dose.effect.list[[i]] * exp(modifier.effect)
                G.mat.list[[i]] <- c(mu * exp(modifier.effect))
                Q.mat[, theta.list[[i+1]]] <-
                    cbind(dose.mat * c(exp(modifier.effect)),
                          modifier.mat * c(dose.effect * exp(modifier.effect)))
            }
        }
        if (length(dose.effect.list)) {
            dose.err <- 1 + Reduce("+", dose.effect.list)
            mu <- mu * dose.err
            Q.mat <- Q.mat / c(dose.err)
        }
        Q.mat[, 1:ncol(model.base)] <- model.base
        colnames(Q.mat) <- theta.name
        result <- list(mu=mu, Q.mat=Q.mat, G.mat.list=G.mat.list)
        if (!is.null(y)) {
            score <- t(Q.mat) %*% (y - mu) ## score vector
            expinf <- t(c(mu) * Q.mat) %*% Q.mat ## information matrix
            result$score <- score
            result$expinf <- expinf
        }
        return (result)
    }))
}

## Use parsed ERR structure to fit the model. If not given,
## 'theta.start' will be set in default to all 0's. 'index' is a
## vector indexing the free parameters to be fitted.
fitERR <- function(err.struct, theta.start, n.iter=50, epsilon=1e-4,
                   step.depth=5, index) {
    theta.n <- err.struct$theta.n
    y <- err.struct$y
    if (missing(theta.start))
        theta.start <- rep(0, theta.n)
    names(theta.start) <- err.struct$theta.name
    if (missing(index))
        index <- rep(TRUE, theta.n)
    if (length(index) != theta.n)
        stop("Invalid index given")
    msi.start <- msiERR(err.struct, theta.start)
    if (any(msi.start$mu < 0))
        stop("Invalid initial value for theta: ", paste(theta.start, ", "))
    logL.start <- sum(dpois(y, msi.start$mu, log=TRUE))
    for (iter in 1:n.iter) {
        score <- msi.start$score[index]
        expinf <- msi.start$expinf[index, index]
        step.size <- 1
        step.depth.current <- step.depth
        if (kappa(expinf) > 1e6)
            expinf <- expinf + diag(1, sum(index))
        inc <- solve(expinf, score)
        theta <- theta.start
        theta[index] <- theta.start[index] + inc * step.size
        while (step.depth.current > 0) {
            msi <- msiERR(err.struct, theta)
            if (any(msi$mu < 0))
                logL <- -Inf
            else
                logL <- sum(dpois(y, msi$mu, log=TRUE))
            if (logL <= logL.start) {
                step.size <- step.size / 2
                step.depth.current <- step.depth.current - 1
                theta[index] <- theta.start[index] + inc * step.size
            } else break
        }
        if (logL < logL.start - epsilon) { # divergence detected
            cat(paste0("Failed to reach convergence. ",
                       "More step depth needed.\n"))
            return (list(theta=theta.start, iter=iter, conv=FALSE))
        }
        if (logL - logL.start < epsilon) { # convergence reached
            return (list(theta=theta, iter=iter, conv=TRUE))
        }
        theta.start <- theta
        msi.start <- msi
        logL.start <- logL
    }
    cat(paste0("Failed to reach convergence. ",
               "More iterations needed.\n"))
    return (list(theta=theta.start, iter=iter, conv=FALSE))
}

## Calculate variance of each parameters based on parsed ERR
## structure, given theta, and covariance information of dose
## replications
varERR <- function(err.struct, theta, dose.reps) {
    W.list <- list()
    msi <- msiERR(err.struct, theta)
    for (i in seq_along(err.struct$U.index)) {
        if (length(err.struct$U.index[[i]])) {
            M <- (msi$Q.mat * msi$G.mat.list[[i]])
            index <- length(W.list) + 1:length(err.struct$U.index[[i]])
            W.list[index] <- lapply(dose.reps[index], function(x) t(M) %*% x)
        }
    }
    cov.W.list <- lapply(W.list, function(x) cov(t(x)))
    U.index <- unlist(err.struct$U.index)
    nU.index <- setdiff(1:err.struct$theta.n, U.index)
    cov.list <- list()
    inv.expinf <- solve(msi$expinf)
    for (i in seq_along(U.index)) {
        cov.list[[i]] <- solve(msi$expinf, cov.W.list[[i]]) %*% inv.expinf
    }
    theta.CI <- matrix(NA, err.struct$theta.n, 2)
    rownames(theta.CI) <- err.struct$theta.name
    nU.var <- inv.expinf + # calculate CI for error-free covariates
        Reduce(`+`, lapply(seq_along(U.index), function(x)
            cov.list[[x]] * theta[U.index[x]]^2))
    var.mat <- matrix(NA, nrow=err.struct$theta.n, ncol=2)
    for (i in nU.index) { 
        v.1 <- nU.var[i,i]
        var.mat[i,1] <- v.1
    }
    for (i in seq_along(U.index)) { # calculate CI for error-prone covariates
        v.2 <- cov.list[[i]][U.index[i], U.index[i]]
        v.1 <- nU.var[U.index[i], U.index[i]] - theta[U.index[i]]^2 * v.2
        var.mat[U.index[i],] <- c(v.1, v.2)
    }
    return (var.mat)
}


## Calulate confidence intervals based on parsed ERR structure and
## estimated theta. If desired, score-type CI's and naive CI's can be
## also be obtained.
inferERR <- function(err.struct, theta, dose.reps, alpha=0.05,
                     score.approx=FALSE, naive=FALSE) {
    if (class(dose.reps) != "list")
        dose.reps <- list(dose.reps)
    if (length(unlist(err.struct$U.index)) != length(dose.reps))
        stop("number of dosimetry reps does not match model")
    theta.CI <- matrix(NA, err.struct$theta.n, 2)
    rownames(theta.CI) <- err.struct$theta.name
    v.mat <- varERR(err.struct, theta, dose.reps)
    for (i in seq_along(theta)) {
        theta.CI[i,] <- getCI(theta[i], v.mat[i,1], v.mat[i,2], alpha=alpha)
    }
    result <- list(CI=theta.CI)
    if (score.approx) {
        theta.CI.score <- matrix(NA, err.struct$theta.n, 2)
        rownames(theta.CI.score) <- err.struct$theta.name
        index <- rep(TRUE, length(theta))
        ## define poins to be evaluated relative to original confidence interval
        theta.CI.cutoff.lower <- c(0.75, 1, 1.25)
        theta.CI.cutoff.upper <- c(0.75, 1, 1.25)
        for (i in seq_along(theta)) {
            theta.points.lower <- (theta.CI[i,1] - theta[i]) * theta.CI.cutoff.lower + theta[i]
            theta.points.upper <- (theta.CI[i,2] - theta[i]) * theta.CI.cutoff.lower + theta[i]
            theta.err.tmp <- function(x) {
                theta[i] <- x
                index[i] <- FALSE
                if (any(msiERR(err.struct, theta)$mu < 0)) {
                    return (NA)
                } else {
                    theta <- fitERR(err.struct, theta, index=index)$theta
                    return (varERR(err.struct, theta, dose.reps)[i,])
                }
            }
            theta.var.lower <- do.call(rbind, lapply(theta.points.lower, theta.err.tmp))
            theta.var.upper <- do.call(rbind, lapply(theta.points.upper, theta.err.tmp))
            theta.points <- c(theta.points.lower, theta[i], theta.points.upper)
            theta.var <- rbind(theta.var.lower, v.mat[i,], theta.var.upper)
            valid.index <- !is.na(theta.var[,1])
            theta.points <- theta.points[valid.index]
            theta.var <- theta.var[valid.index,]
            v.1 <- splinefun(theta.points, theta.var[,1], method="monoH.FC")
            if (is.na(v.mat[i,2]))
                v.2 <- NA
            else
                v.2 <- splinefun(theta.points, theta.var[,2], method="monoH.FC")
            theta.CI.score[i,] <- getCI(theta[i], v.1, v.2, alpha)
        }
        result$CI.score <- theta.CI.score
    }
    if (naive) {
        theta.CI.naive <- matrix(NA, err.struct$theta.n, 2)
        rownames(theta.CI.naive) <- err.struct$theta.name
        z.crit <- qnorm(0.975)
        msi <- msiERR(err.struct, theta)
        theta.sd <- sqrt(diag(solve(msi$expinf)))
        for (i in seq_along(theta)) {
            theta.CI.naive[i,] <- c(theta[i] - z.crit * theta.sd[i],
                                    theta[i] + z.crit * theta.sd[i])
        }
        result$CI.naive <- theta.CI.naive
    }
    return (result)
}

## Calculate confidence interval for individual paramter given
## variance of normal and logNormal components.
getCI <- function(theta, v.1, v.2, alpha, eps=1e-8) {
    if (is.function(v.1)) {
        if (!is.function(v.2)) {
            theta.var <- v.1(theta)
            theta.lower.func <- function(b)
                pnorm((theta - b) / sqrt(max(v.1(b), eps))) - (1 - alpha/2)
            theta.upper.func <- function(b)
                pnorm((theta - b) / sqrt(max(v.1(b), eps))) - alpha/2
        } else {
            theta.var <- v.1(theta) + theta^2 * v.2(theta)
            cdf.approx <- function(b) {
                v.1 <- max(v.1(b), eps)
                v.2 <- max(v.2(b), eps)
                logN <- logNormal2Normal(1, sqrt(v.2))
                if (v.1 < 1e-6) {
                    if (b > 0)
                        return (plnorm(theta / b, logN[1], logN[2]))
                    else
                        return (plnorm(theta / b, logN[1], logN[2],
                                       lower.tail=FALSE))
                } else {
                    if (b > 0)
                        return (integrate(function(y){
                            plnorm((theta - y) / b, logN[1], logN[2]) *
                                dnorm(y, 0, sqrt(v.1))}, -Inf, Inf)[[1]])
                    else
                        return (integrate(function(y){
                            plnorm((theta - y) / b, logN[1], logN[2],
                                   lower.tail=FALSE) *
                                dnorm(y, 0, sqrt(v.1))}, -Inf, Inf)[[1]])
                }
            }
            theta.lower.func <- function(b)
                cdf.approx(b) - ( 1- alpha/2)
            theta.upper.func <- function(b)
                cdf.approx(b) - alpha/2
        }
        theta.bound.inc <- 10 * sqrt(theta.var)
        theta.lower.bound <- c(theta, theta - theta.bound.inc)
        theta.upper.bound <- c(theta, theta + theta.bound.inc)
        for (i in 1:20) {
            if (theta.lower.func(theta.lower.bound[2]) < 0) {
                theta.lower.bound <- theta.lower.bound - theta.bound.inc
            } else {
                break
            }
        }
        for (j in 1:20) {
            if (theta.upper.func(theta.upper.bound[2]) > 0) {
                theta.upper.bound <- theta.upper.bound + theta.bound.inc
            } else {
                break
            }
        }
        if (i == 20 || j == 20) {
            print("Adaptive confidence interval cannot be computed")
            return (NA)
        }
        theta.lower <- uniroot(theta.lower.func, theta.lower.bound)$root
        theta.upper <- uniroot(theta.upper.func, theta.upper.bound)$root
        return (c(theta.lower, theta.upper))
    } else {
        if (is.na(v.2)) {
            z.crit <- qnorm(alpha/2, lower.tail=FALSE)
            return(c(theta - z.crit * sqrt(v.1), theta + z.crit * sqrt(v.1)))
        } else {
            logN <- logNormal2Normal(1, sqrt(v.2))
            cdf.approx <- function(b) {
                if (v.1 < 1e-6) {
                    if (b > 0)
                        return (plnorm(theta / b, logN[1], logN[2]))
                    else
                        return (plnorm(theta / b, logN[1], logN[2],
                                       lower.tail=FALSE))
                } else {
                    if (b > 0)
                        return (integrate(function(y){
                            plnorm((theta - y) / b, logN[1], logN[2]) *
                                dnorm(y, 0, sqrt(v.1))}, -Inf, Inf)[[1]])
                    else
                        return (integrate(function(y){
                            plnorm((theta - y) / b, logN[1], logN[2],
                                   lower.tail=FALSE) *
                                dnorm(y, 0, sqrt(v.1))}, -Inf, Inf)[[1]])
                }}
            if ((theta - qnorm(1-alpha/2, 0, sqrt(v.1)))<0) {
                theta.lower.bound <- (theta - qnorm(1 - alpha / 3, 0, sqrt(v.1))) /
                    qlnorm(alpha / 3, logN[1], logN[2])
            } else {
                theta.lower.bound <- (theta - qnorm(1 - alpha / 3, 0, sqrt(v.1))) /
                    qlnorm(1 - alpha / 3, logN[1], logN[2])
            }
            if (theta - qnorm(alpha / 2, 0, sqrt(v.1)) < 0) {
                theta.upper.bound <- (theta - qnorm(alpha / 3, 0, sqrt(v.1))) /
                    qlnorm(1 - alpha/3, logN[1], logN[2])
            } else {
                theta.upper.bound <- (theta - qnorm(alpha / 3, 0, sqrt(v.1))) /
                    qlnorm(alpha / 3, logN[1], logN[2])
            }
            theta.lower <- uniroot(function(b) cdf.approx(b) - 1 + alpha / 2,
                                   c(theta.lower.bound, theta))[[1]]
            theta.upper <- uniroot(function(b) cdf.approx(b) - alpha / 2,
                                   c(theta, theta.upper.bound))[[1]]
            return (c(theta.lower, theta.upper))
        }
    }
}

## Simulate survival data in Poisson data format. 'formula' contains
## the model for generating the survival dat, and data is the
## associated data frame.  We also need to know the variable name for
## 'id', 'time', and 'censor.time' to set piecewise exponential
## intervals and censor time.
simulateERR <- function(formula, data, theta, id, time, censor.time, seed) {
    if (!missing(seed)) set.seed(seed)
    response.name <- as.character(formula[[2]])
    offset.name <- as.character(
        stripFunc(extractFunc(formula, "offset")[[1]][[1]], "offset"))
    intercept <- attr(terms(formula), "intercept")
    formula.component <- paste(unlist(
        extractFunc(formula, "offset")[[2]]), collapse="+")
    if (!intercept) formula.component <- paste0(formula[[2]], "-1")
    formula <- as.formula(paste0("~", formula.component),
                          env=environment(formula))
    pERR <- match.call()
    pERR <- pERR[c(1, match(c("formula", "data"), names(pERR)))]
    pERR[[1]] <- quote(parseERR)
    pERR[[2]] <- formula
    err.struct <- eval(pERR, parent.frame())
    mu <- msiERR(err.struct, theta)$mu
    ##
    err.call <- match.call()
    err.itc.index <- match(c("id", "time", "censor.time"), names(err.call))
    err.temp <- as.environment(data)
    parent.env(err.temp) <- parent.frame()
    id <- eval(err.call[[err.itc.index[1]]], err.temp)
    time <- eval(err.call[[err.itc.index[2]]], err.temp)
    censor.time <- eval(err.call[[err.itc.index[3]]], err.temp)
    unique.id <- unique(id)
    if (sum(diff(id)!=0)+1 != length(unique.id))
        stop("id should be grouped")
    n <- length(unique.id)
    unique.id.start <- match(unique.id, id)
    unique.id.end <- length(id)+1 - match(unique.id, rev(id))
    if (length(censor.time) == length(id))
        censor.time <- censor.time[unique.id.start]
    result <- lapply(1:n, function(i) {
        start <- unique.id.start[i]
        end <- unique.id.end[i]
        index <- start:end
        cuts <- time[index] - time[start]
        event.time <- rpexp(1, rate=mu[index], t=cuts)
        censor.ind <- FALSE
        if (event.time > censor.time[i]) {
            event.time <- censor.time[i]
            censor.ind <- TRUE
        }
        cell.n <- sum((event.time > cuts))
        cell.case <- rep(0, cell.n)
        if (!censor.ind) cell.case[cell.n] <- 1
        cell.pyr <- diff(cuts)[seq(1, by=1, length.out=cell.n)]
        last.index <- sum(cuts < event.time)
        cell.pyr[cell.n] <- event.time - cuts[last.index]
        return (list(event.time, cell.case, cell.pyr,
        (start:end)[seq(1, by=1, length.out=cell.n)]))
    })
    event.time <- unlist(lapply(result, `[`, 1))
    cell.case <- unlist(lapply(result, `[`, 2))
    cell.pyr <- unlist(lapply(result, `[`, 3))
    cell.id.index <- unlist(lapply(result, `[`, 4))
    data <- cbind(cell.case, cell.pyr, data[cell.id.index,], index=cell.id.index)
    colnames(data)[1:2] <- c(response.name, offset.name)
    return (data)
}

inferERR.null <- function(err.struct, theta, index) {
    if (length(index) > 1)
        stop("index should only contain one number")
    if (theta[index] != 0)
        stop("theta provided does not match index")
    msi.null <- msiERR(err.struct, theta)
    ## score.statistic <- msi.null$score[index]^2 * 
    ##     solve(msi.null$expinf)[index, index]
    ## p <- pchisq(score.statistic, 1, lower.tail=FALSE)
    score.statistic <- msi.null$score[index] * sqrt(solve(msi.null$expinf)[index, index])
    return(score.statistic)
}
