# Copyright 2016 Patrick O. Perry
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# improve_tol = 1e-4 suggested by Nocedal & Wright (c1, p. 38)
# curvature_tol = 0.9 suggested by Nocedal & Wright (c2, p. 39)

poisgamma_control <- function(epsilon = 1e-8, maxit = 25, qr_tol = 1e-7,
                              improve_tol = 1e-4, curvature_tol = 0.9,
                              linesearch_maxit = 20, trace = FALSE)
{
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(qr_tol) || qr_tol <= 0)
        stop("value of 'qr_tol' must be > 0")
    if (!is.numeric(improve_tol) || improve_tol <= 0)
        stop("value of 'improve_tol' must be > 0")
    if (!is.numeric(curvature_tol) || curvature_tol <= 0)
        stop("value of 'curvature_tol' must be > 0")
    if (!is.numeric(linesearch_maxit) || linesearch_maxit <= 0)
        stop("maximum number of linesearch iterations must be > 0")

    list(epsilon = epsilon, maxit = maxit, qr_tol = qr_tol,
         improve_tol = improve_tol, curvature_tol = curvature_tol,
         linesearch_maxit = linesearch_maxit, trace = trace)
}


poisgamma_eval <- function(coefficients, x, y, weights, offset, lhetero,
                           control)
{
    ret <- list(indomain=FALSE) # default return value
    
    # model parameters
    eta <- drop(offset + x %*% coefficients)
    a <- exp(lhetero)
    k <- exp(-lhetero)

    mu <- exp(eta)
    if (!(all(is.finite(mu)) & all(mu > 0))) {
        return(ret)
    }

    log_a_mu <- eta + lhetero
    a_mu <- if (a == 0) 0 else exp(log_a_mu)
    a_y <- a * y

    # qr decomposition
    wt <- weights * (mu + a_mu * y) / (1 + a_mu)^2
    xw <- x * sqrt(wt)
    qr <- qr(xw, LAPACK=TRUE)
    R <- qr.R(qr)
    Q <- qr.Q(qr)
    rownames(R) <- colnames(R)

    # deviance
    dev <- if (a == 0) {
        -2 * sum(weights * (y * eta - mu))
    } else -2 * sum(weights * (y * (log_a_mu - log1p(a_mu)) - log1p(a_mu) / a))

    # residuals, socre
    residuals <- ((y - mu) / mu) * ((1 + a_mu) / (1 + a_y))
    score <- drop(t(xw) %*% (sqrt(wt) * residuals))

    # domain check
    indomain <- is.finite(dev) && all(is.finite(score))

    list(eta = eta, mu = mu, residuals = residuals, R = R, rank = qr$rank,
         qr = qr, weights = wt, prior_weights = weights,
         deviance = dev, score = score, indomain = indomain)
}


poisgamma_lhetero_deriv <- function(eval, x, y, weights, lhetero)
{
    a <- exp(lhetero)
    ymax <- max(y)
    mu <- eval$mu
    a_mu <- a * mu
    a_y <- a * y
    y_mu <- y * mu

    j <- seq(0, max(0, ymax - 1))

    da_tab <- c(0, cumsum(a * j / (1 + a * j)))
    da <- sum(weights * (da_tab[y+1]
                         + (if (a == 0) mu else log1p(a_mu) / a)
                         - (1 + a_y)/(1 + a_mu) * mu))

    #da2_tab <- c(0, cumsum((a * j / (1 + a * j))^2))
    #da2 <- da - sum(weights * (da2_tab[y+1]
    #                           + 2 * (if (a == 0) mu else log1p(a_mu) / a)
    #                           - 2 * mu / (1 + a_mu)
    #                           - (1 + a_y) /(1 + a_mu)^2 * a_mu * mu))

    da2_tab <- c(0, cumsum(a * j / (1 + a * j)^2))
    da2 <- sum(weights * (da2_tab[y+1]
                          - (if (a == 0) mu else log1p(a_mu) / a)
                          + (1 - a_y) / (1 + a_mu) * mu
                          + (1 + a_y) /(1 + a_mu)^2 * a_mu * mu))

    dadb <- drop(t(x) %*% (-weights * a_mu * (y - mu)/(1 + a_mu)^2))

    d1 <- da
    Rinv_dadb <- backsolve(eval$R, transpose=TRUE, dadb[eval$qr$pivot])
    d2 <- da2 - sum(Rinv_dadb^2)

    list(deriv1 = d1, deriv2 = d2)
}


poisgamma_fit <- function(x, y, weights = rep(1, nobs), start = NULL,
                          etastart = NULL, mustart = NULL,
                          offset = rep(0, nobs), lhetero = 0,
                          control = list(...), intercept = TRUE,
                          singular_ok = TRUE, ...)
{
    # control
    control <- do.call("poisgamma_control", control)

    # design matrix, dimensions
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if(is.matrix(y)) rownames(y) else names(y)
    nobs <- NROW(y)
    nvar <- ncol(x)

    # weights, offset
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)

    # determine valid range of eta values
    etamax <- .Machine$double.xmax
    etamin <- -(etamax)

    # initial parameters
    if (is.null(mustart)) {
        mustart <- y + 0.5
    }

    eta <- offset
    if (nobs > 0L) {
        mu <- exp(offset)
    } else {
        mu <- numeric()
    }
    wt <- weights * mu

    qr <- qr(x * sqrt(wt), tol=control$qr_tol)
    rank <- qr$rank
    pivot <- qr$pivot

    # compute initial parameters
    if (!is.null(start)) {
        if (length(start) != nvar) {
            stop(gettextf(paste("length of 'start' should equal %d",
                                "and correspond to initial coefs for %s"),
                          nvar, paste(deparse(xnames), collapse=", ")),
                 domain=NA)
        }
    } else {
        if (is.null(etastart)) {
            if (length(mustart) > 0L) {
                etastart <- log(mustart)
            } else {
                etastart <- numeric()
            }
        }
        start <- qr.coef(qr, sqrt(weights) * (etastart - offset))
        start[is.na(start)] <- 0
    }


    # check for rank-deficiency
    xorig <- x
    if (rank < nvar) {
        if (!singular_ok)
            stop("singular fit encountered")
        xdrop <- x[,pivot[(rank+1L):nvar],drop=FALSE]
        x <- x[,pivot[seq_len(rank)],drop=FALSE]
        start <- start[pivot[seq_len(rank)]]
    }

    # constant
    a <- exp(lhetero)
    k <- exp(-lhetero)
    ymax <- max(y)
    if (is.finite(k)) {
        if (ymax > 0) {
            log_table <- c(0, -lhetero, log(seq(k + 1, length.out=ymax - 1)))
            const_table <- cumsum(log_table)
            const <- -2 * sum(weights * const_table[y + 1])
        } else {
            const <- 0
        }
    } else if (is.finite(lhetero)) {
        const <- 2 * sum(weights * y) * lhetero
    } else {
        const <- 0
    }

    if (rank == 0L) { # empty model
        coefficients <- rep(NA, nvar)
        names(coefficients) <- xnames
        R <- matrix(NA, 0, 0)
        effects <- qr.qty(qr, sqrt(wt) * y)
        if (!is.null(xnames))
            names(effects) <- character(nobs)

        log_a_mu <- eta + lhetero
        a_mu <- if (a == 0) 0 else exp(log_a_mu)
        a_y <- a * y

        dev <- if (a == 0) {
            -2 * sum(weights * (y * eta - mu))
        } else {
            -2 * sum(weights * (y * (log_a_mu - log1p(a_mu)) - log1p(a_mu) / a))
        }
        dev <- dev + const
        residuals <- ((y - mu) / mu) * ((1 + a_mu) / (1 + a_y))

        # TODO: add lhetero_deriv1, lhetero_deriv2

        return(list(coefficients = coefficients,
             residuals = residuals,
             fitted_values = mu, effects = effects,
             R = R, qr = qr, rank = rank,
             lhetero = lhetero, linear_predictors = eta,
             deviance = dev,
             iter = 0L, eval = 0L, weights = wt, prior_weights = weights,
             df_residual = nobs, df_null = nobs, y = y, converged = TRUE,
             boundary = FALSE))
    }

    objective <- function(coef)
        poisgamma_eval(coef, x, y, weights, offset, lhetero, control)
    
    coef0 <- start
    obj0 <- objective(coef0)

    if (!obj0$indomain)
        stop("cannot find valid starting values: please specify some",
             call. = FALSE)

    conv <- FALSE
    eval <- 1
    ftol <- control$improve_tol
    gtol <- control$curvature_tol

    for (iter in seq_len(control$maxit)) {
        if (control$trace)
            cat("Deviance = ", obj0$deviance, " Iterations - ", iter, "\n",
                sep = "")

        eta0 <- obj0$eta
        val0 <- 0.5 * (obj0$deviance)
        grad0 <- -(obj0$score)
        search <- numeric(rank)
        search[obj0$qr$pivot] <-
            backsolve(obj0$R, backsolve(obj0$R, transpose=TRUE,
                                        obj0$score[obj0$qr$pivot]))
        deriv0 <- sum(search * grad0)
        search_eta <- drop(x %*% search)

        # Test for convergence
        if ((deriv0)^2 <= 2 * control$epsilon) {
            conv <- TRUE
            break
        }

        # Fall back to gradient descent if Hessian is ill-conditioned
        if (deriv0 >= 0 || .kappa_tri(obj0$R, LINPACK=FALSE) >= 1e8) {
            search <- -grad0
            deriv0 <- sum(search * grad0)
            search_eta <- drop(x %*% search)
        }

        # determine maximum step size to ensure
        #   |eta[i] - eta0[i]| < 10 * (|eta0[i]| + 1)
        #   for all i
        step_max <- 10 * min((abs(eta0) + 1) / abs(search_eta))
        step_max <- min(step_max, .Machine$double.xmax)

        # determine minimum step size to ensure
        #   |eta[i] - eta0[i]| > eps * (|eta0[i]| + 1)
        #   for at least one i
        step_min <- .Machine$double.eps * min((abs(eta0) + 1) / abs(search_eta))
        step_min <- max(step_min, .Machine$double.xmin)

        # determine initial step, shrinking step.max if necessary
        if (step_min <= 1.0 && 1.0 <= step_max) {
            step0 <- 1.0
        } else {
            step0 <- step_min + 0.5 * (step_max - step_min)
        }
        repeat {
            coef <- coef0 + step0 * search
            obj <- objective(coef)
            eval <- eval + 1
            if (obj$indomain)
                break

            step_max <- step0
            if (step0 < 0.01) {
                step0 <- 2^(0.5 * log2(step_min) + 0.5 * log2(step_max))
            } else {
                step0 <- step_min + 0.5 * (step_max - step_min)
            }
            stopifnot(step0 > step_min)
        }

        # perform line search
        lsctrl <- linesearch_control(value_tol = ftol, deriv_tol = gtol,
                                     step_min = step_min, step_max = step_max)
        ls <- linesearch(val0, deriv0, step0, control = lsctrl)

        for (lsiter in seq_len(control$linesearch_maxit)) {
            val <- 0.5 * (obj$deviance)
            grad <- -(obj$score)
            deriv <- sum(search * grad)

            ls <- update(ls, val, deriv)
            if (ls$converged)
                break

            if (control$trace)
                  cat("New step size (", ls$step, ");",
                      " current deviance = ", obj$deviance, "\n", sep = "")
            coef <- coef0 + ls$step * search
            obj <- objective(coef)
            eval <- eval + 1
            stopifnot(obj$indomain)
        }

        if (!ls$converged) {
            warning("poisgamma_fit: line search failed to converge")
            break
        }

        coef0 <- coef
        obj0 <- obj
        rm(coef, obj)
    }

    eta <- obj0$eta
    mu <- obj0$mu
    dev <- obj0$deviance + const
    wt <- obj0$weights

    if (rank < nvar) {
        coefficients <- rep(NA, nvar)
        coefficients[pivot[1L:rank]] <- coef0
        names(coefficients) <- xnames
    } else {
        coefficients <- coef0
    }


    # qr.  This is tricky; we can't just call qr(sqrt(wt) * xorig), because
    # we need control of the pivoting
    qr1 <- obj0$qr
    ## qr1 <- obj0$qr.modified
    if (rank < nvar) {
        useLAPACK <- attr(qr1, "useLAPACK")
        i1 <- seq_len(rank)
        i2 <- rank + seq_len(nobs - rank)
        j1 <- seq_len(rank)
        j2 <- (rank+1L):nvar

        x2 <- qr.qty(qr1, sqrt(wt) * xdrop)
        ## q0 <- qr.Q(obj0$qr)
        ## H <- obj0$H
        ## x2 <- qr.qty(qr1, q0 %*% (H %*% (t(q0) %*% (sqrt(wt) * xdrop))))
        x21 <- x2[i1,,drop=FALSE]
        x22 <- x2[i2,,drop=FALSE]

        # LAPACK fails if nrow(x22) == 0
        LAPACK <- !is.null(useLAPACK) && useLAPACK
        if (LAPACK && nrow(x22) == 0L) {
            qr22 <- structure(list(qr = x22, rank = 0L, qraux = numeric(),
                                   pivot = seq_len(ncol(x22))),
                              useLAPACK=TRUE, class="qr")
        } else {
            qr22 <- qr(x22, LAPACK=LAPACK)
        }

        qr <- list()

        qr$qr <- matrix(0, nobs, nvar)
        rownames(qr$qr) <- rownames(qr1$qr)
        colnames(qr$qr) <- c(colnames(qr1$qr), colnames(qr22$qr))
        qr$qr[,j1] <- qr1$qr
        qr$qr[i1,j2] <- x21[,qr22$pivot]
        qr$qr[i2,j2] <- qr22$qr

        qr$rank <- rank
        qr$qraux <- c(qr1$qraux, qr22$qraux)
        qr$pivot <- c(pivot[j1][qr1$pivot], pivot[j2][qr22$pivot])

        class(qr) <- "qr"
        if (!is.null(useLAPACK))
            attr(qr, "useLAPACK") <- useLAPACK ## Important!
    } else {
        qr <- qr1
        qr$pivot <- pivot[qr1$pivot]
    }

    # effects
    effects <- qr.qty(qr, sqrt(wt) * y)
    if (!is.null(xnames))
        names(effects) <- c(xnames[qr$pivot[1L:rank]], rep("", nobs - rank))

    # df
    n_ok <- nobs - sum(weights == 0)
    nulldf <- n_ok - as.integer(intercept)
    resdf <- n_ok - rank

    # derivatives with respect to heterogeneity
    deriv <- poisgamma_lhetero_deriv(obj0, x, y, weights, lhetero)

    list(coefficients = coefficients, residuals = obj0$residuals,
         fitted_values = mu, effects = effects,
         R = obj0$R, qr = qr, rank = rank,
         lhetero = lhetero, linear_predictors = eta,
         deviance = dev,
         iter = iter, eval = eval, weights = wt, prior_weights = weights,
         df_residual = resdf, df_null = nulldf, y = y, converged = conv,
         lhetero_score = deriv$deriv1,
         lhetero_info = -(deriv$deriv2),
         boundary = FALSE)
}
