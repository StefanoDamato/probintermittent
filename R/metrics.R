#' Quantile loss
#'
#' A proper scoring rule for probabilistic forecasts.
#'
#' @param actuals The realisation of the data.
#' @param qforecast The predicted quantile.
#' @param level The quantile level of the quantile.
#' @return A vector containing the value of the loss.
#'
#' @family metrics
#' @export
quantile_loss <- function(actuals, qforecast, level){
  return(2*abs(actuals-qforecast)*ifelse(actuals>qforecast, level, 1-level));
}

#' Interval score
#'
#' A proper scoring rule for probabilistic forecasts.
#'
#' @param actuals The realisation of the data.
#' @param qforecast The predicted quantile.
#' @param level The quantile level of the quantile.
#' @return A vector containing the value of the loss.
#'
#' @family metrics
#' @export
interval_score <- function(actuals, qforecast, level){
  return(qforecast + (1/1-levels)*ifelse(actuals>qforecast, actuals-qforecast, 0));
}

#' Quantile loss scaled by its in-sample LOO value.
#'
#' A proper and scale-independent scoring rule for probabilistic forecasts.
#'
#' @param actuals The realisation of the data.
#' @param qforecast The predicted quantile.
#' @param level The quantile level of the quantile.
#' @param y_insample In-sample values of the time series.
#' @return A vector containing the value of the loss.
#'
#' @family metrics
#' @export
scaled_quantile_loss <- function(actuals, qforecast, level, y_insample){

  y_insample <- sort(y_insample);
  K <- ceiling((length(y_insample)-1)*level);
  scale = c((y_insample[K+1]-y_insample[1:K])*(1-level),
            (y_insample[(K+1):length(y_insample)]-y_insample[K])*level);

  return(quantile_loss(actuals, qforecast, level)/(2*mean(scale)));
}

#' Interval score scaled by its in-sample LOO value.
#'
#' A proper and scale-independent scoring rule for probabilistic forecasts.
#'
#' @param actuals The realisation of the data.
#' @param qforecast The predicted quantile.
#' @param level The quantile level of the quantile.
#' @param y_insample In-sample values of the time series.
#' @return A vector containing the value of the loss.
#'
#' @family metrics
#' @export
scaled_interval_score <- function(actuals, qforecast, level, y_insample){

  y_insample <- sort(y_insample);
  K <- ceiling((length(y_insample)-1)*level);
  scale = c(rep(y_insample[K+1], K), rep(y_insample[K], length(y_insample)-K) +
              (y_insample[(K+1):length(y_insample)]-y_insample[K])*(1/1-level));

  return(interval_score(actuals, qforecast, level)/(mean(scale)));
}

#' The Brier score.
#'
#' A quadratic loss on a binary classification problem.
#'
#' @param actuals The realisation of the data.
#' @param prob_zero The predicted probability of zero.
#' @return A vector containing the value of the loss.
#'
#' @family metrics
#' @export
brier <- function(actuals, prob_zero){
  return((as.integer(actuals==0) - prob_zero)^2);
}

#' The binary log-loss.
#'
#' The negative log-likelihood of a Bernoulli random variable.
#'
#' @param actuals The realisation of the data.
#' @param prob_zero The predicted probability of zero.
#' @return A vector containing the value of the loss.
#'
#' @family metrics
#' @export
log_loss <- function(actuals, prob_zero){
  return(log(ifelse(as.integer(actuals==0), prob_zero, 1-prob_zero)));
}

#' Forecast coverage.
#'
#' The proportion of forecasts that are greater or equal than the true value.
#'
#' @param actuals The realisation of the data.
#' @param qforecasts The predicted quantile.
#' @return A vector containing the value of the loss.
#'
#' @family metrics
#' @export
coverage <- function(actuals, qforecast){
  return(as.integer(qforecast >= actuals));
}

#' The log-score.
#'
#' The log-likelihood of the forecast, evaluated on the data.
#'
#' @param actuals The realisation of the data.
#' @param samples The samples from the forecast distribution.
#' @return A vector containing the value of the loss.
#'
#' @family metrics
#' @export
log_score <- function(actuals, samples){
  likelihood <- rep(NA, length(actuals));
  for (i in 1:h){
    likelihood[i] <-  mean(samples[i,] == actuals[i]);
  }
  return(log(likelihood));
}


#' The binary log-loss.
#'
#' The negative log-likelihood of a Bernoulli random variable.
#'
#' @param actuals The realisation of the data.
#' @param prob_zero The predicted probability of zero.
#' @return A vector containing the value of the loss.
#'
#' @family metrics
#' @export
log_loss <- function(actuals, prob_zero){
  return(log(ifelse(as.integer(actuals==0), prob_zero, 1-prob_zero)));
}

#' Forecast coverage.
#'
#' The proportion of forecasts that are greater or equal than the true value.
#'
#' @param actuals The realisation of the data.
#' @param qforecasts The predicted quantile.
#' @return A vector containing the value of the loss.
#'
#' @family metrics
#' @export
coverage <- function(actuals, qforecast){
  return(as.integer(qforecast >= actuals));
}

#' t-test with FDR correction on test.
#'
#' Selection of the best method based on the mean (paired t-test).
#'
#' @param scores A dataframe where scores are stored as columns.
#' @param alpha The level of the test.
#' @return A vector containing neames of methods equivalent to the best.
#'
#' @family metrics
#' @export
test_FDR <- function(scores, alpha) {

  best_method <- names(which.min(colMeans(scores, na.rm = TRUE)));
  methods <- setdiff(colnames(scores), best_method);

  pvals <- sapply(methods,
                  function(method) t.test(scores[[best_method]],
                                          scores[[method]],
                                          paired = TRUE)$p.value)

  n <- length(pvals)
  pvals <- sort(pvals)
  thresholds <- ((1:n)/n)* alpha;
  reject <- pvals <= thresholds;

  if (any(reject)){
    if (all(reject)){
      equivalent_methods <- c()
    } else{
      equivalent_methods <-names(pvals[pvals > pvals[reject]]);
    }
  } else{
  equivalent_methods <- names(pvals)
  }

  return(list('best_method' = best_method,
              'equivalent_methods' = rev(equivalent_methods)))
}

