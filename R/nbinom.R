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

nbinom_control <- function(fit_control = list(...), ...)
{
    fit_control <- do.call("poisgamma_control", fit_control)
}


nbinom <- function(formula, data, weights, subset, na.action,
                   start = NULL, etastart, mustart, offset,
                   control = list(), model = TRUE, x = FALSE,
                   y = TRUE, contrasts = NULL)
{
    # call
    call <- match.call()

    # data
    if (missing(data))
        data <- environment(formula)

    # model frame
    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    # terms
    mt <- attr(mf, "terms")

    # response
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    
    # predictors
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)

    # weights
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")

    # offset
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }

    # starting values
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")

    fit <- eval(call("nbinom_fit", x = X, y = Y,
                     weights = weights, start = start, etastart = etastart,
                     mustart = mustart, offset = offset,
                     control = control,
                     intercept = attr(mt, "intercept") > 0L))
    fit$contrasts <- attr(X, "contrasts")


    xlevels <- .getXlevels(mt, mf)
    fit$xlevels <- xlevels

    fit$call <- call
    fit$control <- control
    fit$terms <- mt
    if (model)
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x)
        fit$x <- X
    if (!y)
        fit$y <- NULL

    class(fit) <- "nbinom"

    fit
}


coef.nbinom <- function(object, ...)
{
    object$coefficients
}


print.nbinom <- function(x, digits = max(3L, getOption("digits") - 3L),
                        signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients")
        if (is.character(co <- x$contrasts)) 
            cat("  [contrasts: ", apply(cbind(names(co), co), 
                1L, paste, collapse = "="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nHeterogeneity: ", format(exp(x$lhetero), digits = digits), "\n", sep="")
    cat("\n")
    invisible(x)
}


predict.nbinom <- function (object, newdata = NULL, type = c("link", "response", 
    "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL, 
    na.action = na.pass, ...)
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    class(object) <- c(class(object), "lm") # suppress warnings to calls of predict.lm
    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, link = object$linear_predictors, 
                response = object$fitted_values, terms = predict.lm(object, 
                  se.fit = se.fit, scale = 1, type = "terms", 
                  terms = terms))
            if (!is.null(na.act)) 
                pred <- napredict(na.act, pred)
        }
        else {
            pred <- predict.lm(object, newdata, se.fit, scale = 1, 
                type = ifelse(type == "link", "response", type), 
                terms = terms, na.action = na.action)
            switch(type, response = {
                pred <- exp(pred)
            }, link = , terms = )
        }
    }
    else {
        stop("not implemented")
    }
    pred
}
