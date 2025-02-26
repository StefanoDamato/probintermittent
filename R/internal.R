#' Initialize the data.
#'
#' This function extrapolates in-sample, holdout, frequency, starting time and
#' length from a time series. It is intended for internal use only and is not exported.
#'
#' @param data A time series or a numeric vector or 1-dimensional array.
#' @param holdout Whether to remove the last h values from the in-sample part.
#' @param h The number of steps ahead used as forecasting horizon.
#' @return A list with the following components:
#' \describe{
#'   \item{y_insample}{The in-sample part of the data.}
#'   \item{y_holdout}{The holdout part of the data.}
#'   \item{freq}{The time frequency.}
#'   \item{pred_start}{The starting time for the predictions.}
#'   \item{L}{The length of the in-sample part.}
#' }
#' @details This function is used within other package functions.
#'
#' @keywords internal
init_data <- function(data, holdout, h){

    y <- as.vector(data);
    L <- length(y) - h*holdout;
    if (L <= 0) stop("Not enough in-sample observations");

    y_insample <- y[1:L];
    y_holdout <- ifelse(holdout, y[L+1:h], rep(NA, h));
    freq <- frequency(data);
    pred_start <- time(data)[L] + deltat(data);

    return(list("y_insample" = y_insample, "y_holdout" = y_holdout,
                "freq" = freq, "pred_start" = pred_start, "L" = L));
}

#' Produce the forecasts from the samples.
#'
#' This functions generates mean forecast with confidence intervals
#' and predicts the probabilities of zero starting from samples.
#' It is intended for internal use only and is not exported.
#'
#' @param forecast_samples An array containing predictive sample paths.
#' @param levels The confidence levels for predictive intervals.
#' @param cumulative Whether to just return the sum over the forecast horizon.
#' @param side The side of the confidence intervals.
#' @param h The number of steps ahead used as forecasting horizon.
#' @return A list with the following components:
#' \describe{
#'   \item{y_forecast}{The predicted mean.}
#'   \item{prob_zero}{The predicted probability of observing zero.}
#'   \item{levels_upper}{The upper confidence levels of the interval.}
#'   \item{levels_lower}{The lower confidence levels of the interval.}
#'   \item{CI_upper}{The upper bounds of the confidence region.}
#'   \item{CI_lower}{The lower bounds of the confidence region.}
#' }
#' @details This function is used within other package functions.
#'
#' @keywords internal
get_forecast = function(forecast_samples, levels, cumulative, side, h){

    levels_upper <- levels_lower <- rep(NA, length(levels));

    h_steps <- ifelse(cumulative, 1, h);

    if (side == "both"){
        levels_upper[] <- (1+levels)/2;
        levels_lower[] <- (1-levels)/2;
    }
    else if (side == "upper"){
        levels_upper[] <- levels;
        levels_lower[] <- 0;
    }
    else if( side == "lower"){
        levels_upper[] <- 1;
        levels_lower[] <- 1-levels;
    }

    CI_upper <- CI_lower <- matrix(NA, h_steps, length(levels));
    y_forecast <- prob_zero <- rep(NA, h_steps);

    if (cumulative){
        forecast_samples[] <- rowSums(forecast_samples);
        y_forecast[] <- mean(forecast_samples);
        prob_zero[] <- mean(1*(forecast_samples == 0));
        CI_upper[] <- quantile(forecast_samples, levels_upper);
        CI_lower[] <- quantile(forecast_samples, levels_lower);
    }
    else{
        y_forecast[] <- colMeans(forecast_samples);
        prob_zero[] <- colMeans(1*(forecast_samples == 0));
        CI_upper[] <- t(apply(forecast_samples, 2, quantile, levels_upper));
        CI_lower[] <- t(apply(forecast_samples, 2, quantile, levels_lower));
    }

    CI_upper[,levels_upper>=1] <- Inf;
    CI_lower[,levels_lower<=0] <- 0;

    return(list("y_forecast" = y_forecast, "prob_zero" = prob_zero,
                "levels_upper" = levels_upper, "levels_lower" = levels_lower,
                "CI_upper" = CI_upper, "CI_lower" = CI_lower));

}
