#model functions

#' Exponentially weighted Poisson regression model
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family choice of "ewp2" or "ewp3"
#' @param data a data frame containing the variables in the model.
#' @param verbose logical, defaults to TRUE; print model fitting progress
#' @param method string, passed to optim, defaults to 'BFGS'
#' @param hessian logical, defaults to TRUE; calculate Hessian?
#' @param autoscale logical, defaults to TRUE; automatically scale model parameters inside the optimisation routine based on initial estimates from a Poisson regression.
#' @param maxiter numeric, maximum number of iterations for optim
#' @param sum_limit numeric, defaults to 3*maximum count; upper limit for the sum used for the normalizing factor.
#'
#' @return an ewp model
#' @importFrom stats .getXlevels coef delete.response glm.fit model.frame model.matrix model.response na.omit na.pass optim optimHess poisson terms
#' @export
#'
ewp_reg <- function(formula, family = 'ewp3', data, verbose = TRUE, method = 'Nelder-Mead', hessian = TRUE, autoscale = TRUE, maxiter = 5000, sum_limit = round(max(Y)*3)){
  cl <- match.call()
  mt <- terms(formula, data = data)
  #if(missing(data)) data <- environment(formula)
  #mf <- match.call(expand.dots = FALSE)
  #m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  #mf <- mf[c(1, m)]
  #mf$drop.unused.levels <- TRUE

  ## call model.frame()
  #mf[[1]] <- as.name("model.frame")
  #mf <- eval(mf, parent.frame())
  mf <- model.frame(formula, data, drop.unused.levels = TRUE)

  ## model matrix, response
  mm <- model.matrix(formula, mf)
  Y <- model.response(mf, "numeric")

  if(any(Y >= 30)) ("Counts >= 30 detected. The likelihood estimation procedure is not currently set up to deal with this.")
  if(any(Y > 20)) warning("Counts > 20 detected. The likelihood estimation procedure is not currently set up to deal with counts in excess of 30. Results may be misleading if lambda >= 25 and beta2 < 1.")


  # set a warning message if the sum_limit < 15. Occurs for small values of max count (<5), or when sum_limit is set manually.
  if (sum_limit<15){
    warning("sum_limit < 15 detected. A sum_limit of 3 times the maximum count value is recommended, or of
            at least 15, in the case of a small maximum count.")
  }

  #get start values for lambda linpred from simple poisson regression
  start_values <- coef(glm.fit(x = mm, y = Y, family = poisson()))
  #estimate relative effect sizes for optim - assumes dispersion parameter is approx 1!
  if (autoscale) {
     parscale_est = abs(c(start_values, beta1 = 1, beta2 = 1)/start_values[1])
  } else {
     parscale_est = rep(1, length(start_values)+2)
    }
  #add dispersion parameter start values
  start_values = c(start_values, beta1 = 0, beta2 = 0)
  if(verbose){
    cat('start values are: \n')
    print(start_values)
  }

  pllik3 <- function(par, mm, Y){
    lambda = exp(mm %*% par[1:ncol(mm)])
    beta1 = unname(par['beta1'])
    beta2 = unname(par['beta2'])
    #ll = numeric(nrow(mm))
    #for (i in 1:nrow(mm)){
    #  ll[i] = log(dewp3_cpp(Y[i],lambda[i],beta1,beta2))
    #}
    #return(-1*sum(ll))
    return(pllik3_part_cpp(Y, lambda, beta1, beta2, sum_limit))
  }

  resultp3 <- optim(par = start_values,
                    fn = pllik3, mm = mm, Y = Y,
                    method = method,
                    hessian = FALSE,
                    control = list(trace = verbose,
                                   REPORT=4*verbose,
                                   ndeps=rep(1e-5, ncol(mm)+2),
                                   parscale = parscale_est,
                                   maxit = maxiter))

  if(hessian){
    if(verbose) cat('\nCalculating Hessian. This may take a while.\n')
    resultp3$hessian <- optimHess(resultp3$par, fn = pllik3, mm = mm, Y = Y)
    #estimate vcov
    vc = solve(resultp3$hessian)
  } else {
    vc <- resultp3$hessian <- matrix(NA_real_, nrow = ncol(mm) + 2, ncol = ncol(mm) + 2)
  }



  ## fitted and residuals
  Yhat = exp(mm %*% resultp3$par[1:ncol(mm)])
  res <- (Y - Yhat)

  #output structure
  out <- list(
    coefficients = resultp3$par,
    vcov = vc,
    se = sqrt(diag(vc)),
    optim = resultp3,
    loglik = -resultp3$value,
    residuals = res,
    fitted.values = as.vector(Yhat),
    terms = mt,
    call = cl,
    levels = .getXlevels(mt, mf),
    start = start_values,
    n = length(Y),
    df.residual = length(Y) - ncol(mm) - 2,
    converged = resultp3$convergence < 1,
    formula = formula,
    sum_limit=sum_limit,
    dist= 'ewp3'
  )
  class(out) <- "ewp"
  return(out)
}

#' Extract coefficients
#'
#' @param object an object of class ewp
#' @param ... ignored
#'
#' @return a vector of coefficient values. Beware that the lambda parameters are on the log-link scale, whereas the betas are estimated using an identity link.
#' @importFrom stats coef
#' @export
#'
coef.ewp <- function(object, ...) {
  object$coefficients
}

#' Extract estimated variance-covariance matrix
#'
#' @param object an object of class ewp
#' @param ... ignored
#'
#' @return a matrix
#' @importFrom stats vcov
#' @export
#'
vcov.ewp <- function(object, ...) {
  object$vcov
}

#' Extract log likelihood
#'
#' @param object an object of class ewp
#' @param ... ignored
#'
#' @return a numeric
#' @importFrom stats logLik
#' @export
#'
logLik.ewp <- function(object, ...) {
  structure(object$loglik, df = object$n - object$df.residual, nobs = object$n, class = "logLik")
}

#' Extract fitted values
#'
#' @param object an object of class ewp
#' @param ... ignored
#'
#' @return a vector of fitted values on the response scale
#' @importFrom stats fitted
#' @export
#'
fitted.ewp <- function(object, ...) {
  object$fitted.values
}


#' Print ewp model object
#'
#' @param x ewp model object
#' @param digits digits to print
#' @param ... ignored
#'
#' @export
#'
print.ewp <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    dist <- x$dist
    fixed <- FALSE

  #cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
      } else {
    cat(paste("Coefficients (", dist, " with log link on lambda):\n", sep = ""))
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    }

  invisible(x)
}

#' Model summary
#'
#' @param object ewp model fit
#' @param ... ignored
#'
#' @importFrom stats pnorm
#' @importFrom utils tail
#' @export
#'
summary.ewp <- function(object,...)
{
  ## deviance residuals
  #object$residuals <- residuals(object, type = "deviance")

  ## compute z statistics
  cf <- object$coefficients
  se <- sqrt(diag(object$vcov))
  k <- length(cf)


  zstat <- cf/se
  pval <- 2*pnorm(-abs(zstat))
  cf <- cbind(cf, se, zstat, pval)
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients <- cf

  ## number of iterations
  object$iterations <- tail(na.omit(object$optim$count), 1)

  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.ewp"
  object
}

#' Print ewp model summary
#'
#' @param x ewp model summary
#' @param digits number of digits to print
#' @param ... additional arguments to printCoefmat()
#'
#' @return printout of the summary object
#' @importFrom stats printCoefmat
#' @export
#'
print.summary.ewp <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  #cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {


      dist <- x$dist
      fixed <- FALSE
      npar_lambda <- nrow(x$coefficients) - 2


    cat("Deviance residuals:\n")
    #print(structure(quantile(x$residuals),
    #                names = c("Min", "1Q", "Median", "3Q", "Max")), digits = digits, ...)

    cat(paste("\nlambda coefficients (", dist, " with log link):\n", sep = ""))
    printCoefmat(x$coefficients[1:npar_lambda, , drop = FALSE], digits = digits, ...)
    cat(paste("\ndispersion coefficients:\n"))
    printCoefmat(x$coefficients[(npar_lambda+1):(npar_lambda+2),], digits = digits, ...)

    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations, "\n"))
    cat("Log-likelihood:", formatC(x$loglik, digits = digits), "on", x$n - x$df.residual, "Df\n")
  }

  invisible(x)
}


#' Predict from fitted model
#'
#' @param object ewp model object
#' @param newdata optional data.frame
#' @param type character; default="response", no other type implemented
#' @param na.action  defaults to na.pass()
#' @param ... ignored
#'
#' @return a vector of predictions
#' @importFrom stats predict
#' @importFrom stats weighted.mean
#' @export
#'
predict.ewp <- function(object, newdata, type = c("response"),
                              na.action = na.pass, ...)
{
  type <- match.arg(type)

  ## if no new data supplied
  if(missing(newdata)) {
    if(type != "response") {
      stop('Unknown prediction type')
    } else {
      x <- seq(0,object[["sum_limit"]], by= 1)

      pred_ewp <- vector()
      for (i in 1:length(object$fitted.values)){
        pmf_ewp <- dewp3(x, object$fitted.values[i], object$coefficients[["beta1"]],object$coefficients[["beta2"]])
        pred_ewp[i] <- weighted.mean(x,w=(pmf_ewp))
      }
      return(pred_ewp)
    }
  } else {
    mf <- model.frame(delete.response(object$terms), newdata, na.action = na.action, xlev = object$levels)
    X <- model.matrix(delete.response(object$terms), mf)#, contrasts = object$contrasts)
    offset <- if(!is.null(off.num <- attr(object$terms, "offset")))
      eval(attr(object$terms, "variables")[[off.num + 1]], newdata)
    else if(!is.null(object$offset)) eval(object$call$offset, newdata)
    if(is.null(offset)) offset <- rep(0, NROW(X))
  }

  rval <- exp(X %*% object$coefficients[1:ncol(X)])[,1]

  x <- seq(0,object[["sum_limit"]], by= 1)

  pred_ewp <- vector()
  for (i in 1:nrow(newdata)){
    pmf_ewp <- dewp3(x, rval[i], object$coefficients[["beta1"]],object$coefficients[["beta2"]])
    pred_ewp[i] <- weighted.mean(x, w=(pmf_ewp))
  }

  return(pred_ewp)
}


# ## covariances
# vc <- if(hessian) {
#   tryCatch(-solve(as.matrix(fit$hessian)),
#            error = function(e) {
#              warning(e$message, call = FALSE)
#              matrix(NA_real_, nrow = k + (dist == "negbin"), ncol = k + (dist == "negbin"))
#            })
# } else {
#   matrix(NA_real_, nrow = k + (dist == "negbin"), ncol = k + (dist == "negbin"))
# }
# if(dist == "negbin") {
#   SE.logtheta <- as.vector(sqrt(diag(vc)[k + 1]))
#   vc <- vc[-(k+1), -(k+1), drop = FALSE]
# } else {
#   SE.logtheta <- NULL
# }
# colnames(vc) <- rownames(vc) <- colnames(X)


#' simulate from fitted model
#'
#' @param object ewp model object
#' @param nsim number of response vectors to simulate. Defaults to 1.
#' @param ... ignored
#'
#' @return a data frame with `nsim` columns.
#'
#' @importFrom stats simulate
#' @export
#'
simulate.ewp <- function(object, nsim=1, ...){
  ftd <- fitted(object)
  n <- length(ftd)
  ntot <- n * nsim

  ncoef <- length(coef(object))
  nm <- names(ftd)
  val <- matrix(NA_integer_, nrow = n, ncol = nsim)
  for (i in 1:nsim){
    val[,i] <- vapply(ftd, function(x)rewp3(n = 1, lambda = x, beta1 = coef(object)[ncoef-1], beta2 = coef(object)[ncoef]), numeric(1))
  }

  val <- as.data.frame(val)
  names(val) <- paste0("sim_", seq_len(nsim))
  if (!is.null(nm))
    row.names(val) <- nm
  val
}
