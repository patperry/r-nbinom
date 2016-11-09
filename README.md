
Overview
========

This package implements negative binomial regression. The mean of the
distribution is parametrized as

    log(mu) = x %*% coefficients + offset

The "heterogeneity" parameter is the inverse of the gamma shape parameter,
so that the variance is given by

    var = mu + exp(lhetero) * mu^2

Here, `lhetero` is the logarithm of the heterogeneity.


Demonstration
=============

Generate negative binomial data using the `rnbinom` function from the
`stats` package:

    set.seed(0)
    nobs <- 10000

    hetero <- 3
    rate <- 0.5
    exposure <- runif(nobs, 0, 50)

    y <- rnbinom(nobs, mu = rate * exposure, size = 1/hetero)
    

Fit a negative binomial model:

    (fit <- nbinom(y ~ 1, offset = log(exposure)))


Evaluate the fitted log-likelihood on the data:

    hetero_est <- exp(fit$lhetero)
    mu_est <- predict(fit, type="response")
    lik <- dnbinom(y, mu=mu_est, size=1/hetero_est, log=TRUE)


Evaluate the fitted likelihood on new data:

    y1 <- c(5, 10, 15)
    mu1 <- predict(fit, newdata=data.frame(exposure=c(10, 20, 30)),
                   "response")
    lik1 <- dnbinom(y, mu=mu1, size=1/hetero_est, log=TRUE)
