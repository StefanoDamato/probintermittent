#' Gamma-Poisson dynamic model
#'
#' Bayesian model with Gamma prior and Poisson likelihood
#'
#' Short description #TODO.
#'
#' @param data A time series or a numeric vector or 1-dimensional array.
#' @param h The number of steps ahead used as forecasting horizon.
#' @param levels The confidence levels for predictive intervals.
#' @param holdout Whether to remove the last h values from the in-sample part.
#' @param cumulative Whether to just return the sum over the forecast horizon.
#' @param side The side of the confidence intervals.
#' @param nsim The amount of predictive sample paths.
#' @param paths Whether to return sample trajectories.
#'
#' @return The function returns a model of a class "#TODO", which contains:
#' \describe{
#'   \item{model}{The model name.}
#'   \item{call}{The arguments of the function call.}
#'   \item{params}{Internal parameters for the estimation of the model.}
#'   \item{forecast}{The predicted mean.}
#'   \item{probzero}{The predicted probability of observing zero.}
#'   \item{levels_upper}{The upper confidence levels of the interval.}
#'   \item{levels_lower}{The lower confidence levels of the interval.}
#'   \item{upperCI}{The upper bounds of the confidence region.}
#'   \item{lowerCI}{The lower bounds of the confidence region.}
#'   \item{samples}{The predictive sample paths.}
#' }
#'
#' @references
#' \itemize{
#' \item paper to cite # TODO
#' }
#'
#' @author Name Surname, \email{name@domain.#TODO}
#'
#' @keywords some words #TODO
#'
#' @seealso \code{\link[package]{function}, \link[package2]{function2}} #TODO
#'
#' @examples
#' y <- c(rpois(30,0.9), rpois(20,0.3))
#' test <- gampois(y)
#'
#' @importFrom nloptr nloptr
#' @importFrom stats rpois rgamma ts
#'
#' @family models
#' @export
gampois = function(data, h=10, levels=0.9, holdout=FALSE, cumulative=FALSE,
                   side=c("upper", "both", "lower"), nsim=10000, paths=FALSE){

    # Check the arguments
    if (h <= 0 | (h!=as.integer(h))) stop("h should be a positive integer");
    if (any(levels<0) | any(levels>1)) stop("levels should be in [0, 1]");
    if (!is.logical(holdout)) stop("holdout should be a boolean");
    if (!is.logical(cumulative)) stop("cumulative should be a boolean");
    side <- match.arg(side);
    if (nsim <= 0) stop("nsim should be a positive integer");
    if (!is.logical(paths)) stop("paths should be a boolean");

    # Save the arguments of the function in a list
    call <- list("data" = data, "h" = h, "levels" = levels, "holdout" = holdout,
                 "cumulative" = cumulative, "side" = side, "nsim" = nsim,
                 "paths" = paths);

    # Pre-process the data to get in-sample, holdout and infos
    data_list <- init_data(data, holdout, h);
    y_insample <- data_list$y_insample;
    y_holdout <- data_list$y_houldout;
    freq <- data_list$freq;
    pred_start <- data_list$pred_start;
    L <- data_list$L;

    # Define the negative log-likelihood of the data
    nll <- function(y_insample, a0, b0, w){
      gammaparams <- gammaDynamic(y_insample, a0, b0, w);
      a <- gammaparams$a;
      b <- gammaparams$b;
      return(-mean(lchoose(a+y_insample-1, y_insample) +
                     a*log(b) - (a+y_insample)*log(1+b)));
    }

    # Optimise the parameters
    opt <- nloptr(
      x0 = c(1, 1, 0.8),
      eval_f = function(x) nll(y_insample, x[1], x[2], x[3]),
      lb = rep(0, 3) + epsilon,
      ub = c(Inf, Inf, 1) - epsilon,
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 200)
    );

    # Get the optimal values
    x <- opt$solution;
    a0 <- x[1];
    b0 <- x[2];
    w <- x[3];

    # Run the dynamic process
    gammaparams <- gammaDynamic(y_insample, a0, b0, w);
    a <- gammaparams$a;
    b <- gammaparams$b;

    # Set the last states for the smoothing
    forecast_samples <- matrix(NA, nsim, h);
    a_state <- rep(w*a[L] + y_insample[L], nsim);
    b_state <- rep(w*b[L] + 1, nsim);

    # Hierachically generate sample paths
    for (i in 1:h){
      lambda_state <- rgamma(nsim, a_state, b_state);
      y_new <- rpois(nsim, lambda_state);
      forecast_samples[,i] <- y_new;
      a_state <- w*a_state + y_new;
      b_state <- w*b_state + 1;
    }

    # Save relevant parameters
    params = list("a0" = a0, "b0" = b0, "w" = w);

    # Store the results in some arrays
    forecast_list <- get_forecast(forecast_samples, levels, cumulative, side, h);
    mean_forecast <- forecast_list$mean_forecast;
    prob_zero <- forecast_list$prob_zero;
    levels_upper <- forecast_list$levels_upper;
    levels_lower <- forecast_list$levels_lower;
    CI_upper <- forecast_list$CI_upper;
    CI_lower <- forecast_list$CI_lower;

    # Save the results as time series
    meanforecast <- ts(mean_forecast, start=pred_start ,frequency=freq);
    probzero <- ts(prob_zero, start=pred_start, frequency=freq);
    upperCI <- ts(CI_upper, start=pred_start, frequency=freq, names=levels_upper);
    lowerCI <- ts(CI_lower, start=pred_start, frequency=freq, names=levels_lower);
    if (paths==TRUE){
      samples <- ts(t(forecast_samples), start=pred_start, frequency=freq,
                    names=paste0('Path ', 1:nsim));
    } else{
      samples <- NULL;
    }

    # Instantiate the predictions in a class
    model <- list("model" = "GamPois", "call" = call, "params" = params,
                  "meanforecast" = meanforecast, "probzero" = probzero,
                  "levels_upper" = levels_upper, levels_lower = "levels_lower",
                  "upperCI" = upperCI, "lowerCI" = lowerCI, "samples" = samples);
    return(structure(model, class="#TODO"))
}





