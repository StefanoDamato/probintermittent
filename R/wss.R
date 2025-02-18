#' WSS
#'
#' Markov-Chain-based bootstrapping model
#'
#' The model separates occurrence and the demand. The occurrence is estimated
#' assuming a Markov Chain whose transition matrix is estimated from past data.
#' The demand is resampled from past data with Gaussian jittering.
#'
#' @param data A time series or a numeric vector or 1-dimensional array.
#' @param h The number of steps ahead used as forecasting horizon.
#' @param levels The confidence levels for predictive intervals.
#' @param holdout Whether to remove the last h values from the in-sample part.
#' @param cumulative Whether to just return the sum over the forecast horizon.
#' @param side The side of the confidence intervals.
#' @param nsim The amount of predictive sample paths.
#'
#' @return Function returns a model of a class "#TODO", which contains:
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
#' y <- c(rpois(50,0.3),rpois(50,0.8))
#' test <- wss(y)
#'
#' @importFrom stats rnorm ts
#'
#' @family models
#' @export
wss = function(data, h=10, levels=0.9, holdout=FALSE, cumulative=FALSE,
               side=c("upper", "both", "lower"), nsim=10000){

    # Check the arguments
    if (h <= 0 | (h!=as.integer(h))) stop("h should be a positive integer");
    if (any(levels<0) | any(levels>1)) stop("levels should be in [0, 1]");
    if (!is.logical(holdout)) stop("holdout should be a boolean");
    if (!is.logical(cumulative)) stop("cumulative should be a boolean");
    side <- match.arg(side);
    if (nsim <= 0) stop("nsim should be a positive integer");

    # Save the arguments of the function in a list
    call <- list("data" = data, "h" = h, "levels" = levels, "holdout" = holdout,
                 "cumulative" = cumulative, "side" = side, "nsim" = nsim);

    # Pre-process the data to get in-sample, holdout and infos
    data_list <- init_data(data, holdout, h);
    y_insample <- data_list$y_insample;
    y_holdout <- data_list$y_houldout;
    freq <- data_list$freq;
    pred_start <- data_list$pred_start;

    # Compute Croston's decomposition
    decomp <- crostonsdecomp(y_insample);
    occurrence <- decomp$occurrence;
    demand <- decomp$demand;

    # Estimate the the transition matrix for the Markov Chain
    occ_diff <- 2*occurrence[2:L] - occurrence[1:(L-1)]
    P <- matrix(c(sum(occ_diff == 0), sum(occ_diff == -1),
                  sum(occ_diff == 2), sum(occ_diff == 1)), 2, 2);
    P <- P/rowSums(P);

    # Set the chain in the last observed state
    forecast_samples <- matrix(NA, nsim, h);
    occ_state <- rep(occ_state[L], nsim);

    # Sample sequentially the next h occurrences
    for(i in 1:h){
      P_x1 <- P[,2][occ_state+1];
      occ_state <- as.integer(runif(nsim) <= P_x1);
      forecast_samples[,i] <- occ_state;
    }

    # Sample the demand from past data with jittering
    to_sample <- forecast_samples==1;
    new_demand <- sample(demand, sum(to_sample), replace=TRUE);
    jittered <- 1 + round(new_demand + rnorm(sum(to_sample))*sqrt(new_demand));
    jittered[jittered<=0] <- new_demand[jittered<=0];
    forecast_samples[to_sample] <- jittered;

    # Save relevant parameters
    params <- list("P" = P);

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

    # Instantiate the predictions in a class
    model <- list("model" = "WSS", "call" = call, "params" = params,
                  "meanforecast" = meanforecast, "probzero" = probzero,
                  "levels_upper" = levels_upper, levels_lower = "levels_lower",
                  "upperCI" = upperCI, "lowerCI" = lowerCI);
    return(structure(model, class="#TODO"))
}





