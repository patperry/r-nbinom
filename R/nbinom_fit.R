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


nbinom_fit <- function(x, y, weights = rep(1, nobs), start = NULL,
                       etastart = NULL, mustart = NULL,
                       offset = rep(0, nobs), control = list(...),
                       intercept = TRUE,
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

    fit0 <- poisgamma_fit(x=x, y=y, weights=weights, start=start,
                          etastart=etastart, mustart=mustart, offset=offset,
                          lhetero=0, control=control,
                          intercept=intercept, singular_ok=singular_ok)
    start <- fit0$coefficients

    f <- function(lhetero) {
        fit <- poisgamma_fit(x=x, y=y, weights=weights, start=start, offset=offset,
                             lhetero=lhetero, control=control,
                             intercept=intercept, singular_ok=singular_ok)
        res <- fit$deviance
        attr(res, "fit") <- fit
        res
    }

    opt <- optim(0, f, method="Brent", lower=-100, upper=+100)
    lhetero <- opt$par
    fit <- attr(opt$value, "fit")

    list(coefficients = fit$coefficients, lhetero = lhetero, 
         deviance = fit$deviance,
         rank = fit$rank, pivot=fit$qr$pivot, R = fit$R,
         lhetero_score = fit$lhetero_score,
         lhetero_info = fit$lhetero_info)
}
