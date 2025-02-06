#' Quantile loss
#'
#' A proper scoring rule for probabilistic forecasts.
#'
#' @param actuals The realisation of the data.
#' @param qforecast The predicted quantile.
#' @param level The quantile level of the quantile.
#' @return A vector containing tha value of the loss.
#'
#' @export
quantile_loss <- function(actuals, qforecast, level){
  return(2*abs(actuals-qforecast)*ifelse(actuals>qforecast, level, 1-level));
}

#' Quantile loss scaled by its in-sample LOO value.
#'
#' A proper and scale-independent scoring rule for probabilistic forecasts.
#'
#' @param actuals The realisation of the data.
#' @param qforecast The predicted quantile.
#' @param level The quantile level of the quantile.
#' @param y_insample In-sample values of the time series.
#' @return A vector containing tha value of the loss.
#'
#' @export
quantile_loss_scaled <- function(actuals, qforecast, level, y_insample){

  y_insample <- sort(y_insample);
  K <- ceiling((length(y_insample)-1)*level);
  scale = c((y_insample[K+1]-y_insample[1:K])*level,
            (y_insample[(K+1):length(y_insample)]-y[K])*(1-level));

  return(quantile_loss(actuals, qforecast, level)/mean(scale));
}
