#' Static parametric distribution.
#'
#' Maximise the marginal likelihood of the time series using distributions as Poisson,
#' negative binomial or their hurdle-shifted versions.
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
#' @param rmzeros Whether to remove the initial zeros from the series.
#' @param distr The distribution choice to be fit.
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
#' @author Stefano Damato, \email{stefano.damato@idsia.ch}
#'
#' @keywords some words #TODO
#'
#' @seealso \code{\link[package]{function}, \link[package2]{function2}} #TODO
#'
#' @examples
#' y <- c(rpois(30,0.9), rpois(20,0.3))
#' staticd <- staticd(y)
#'
#' @importFrom stats dnbinom rnbinom dpois rpois runif
#'
#' @family models
#' @export
staticd = function(data, h=10, levels=0.9, holdout=FALSE, cumulative=FALSE,
                   side=c("upper", "both", "lower"), nsim=10000, paths=FALSE,
                   rmzeros=FALSE, distr=c("auto", "pois", "hsp",
                                          "nbinom", "hsnb", "mixture")){

  # Check the arguments
  if (h <= 0 | (h!=as.integer(h))) stop("h should be a positive integer");
  if (any(levels<0) | any(levels>1)) stop("levels should be in [0, 1]");
  if (!is.logical(holdout)) stop("holdout should be a boolean");
  if (!is.logical(cumulative)) stop("cumulative should be a boolean");
  side <- match.arg(side);
  if (nsim <= 0) stop("nsim should be a positive integer");
  if (!is.logical(paths)) stop("paths should be a boolean");
  if (!is.logical(rmzeros)) stop("rmzeros should be a boolean");
  distr <- match.arg(distr);

  # Save the arguments of the function in a list
  call <- list("data" = data, "h" = h, "levels" = levels, "holdout" = holdout,
               "cumulative" = cumulative, "side" = side, "nsim" = nsim,
               "paths" = paths, "rmzeros" = rmzeros, "distr" = distr);

  # Pre-process the data to get in-sample, holdout and infos
  data_list <- init_data(data, holdout, h);
  y_insample <- data_list$y_insample;
  y_holdout <- data_list$y_houldout;
  freq <- data_list$freq;
  pred_start <- data_list$pred_start;
  L <- data_list$L;

  # Compute Croston's decomposition
  decomp <- crostonsdecomp(y_insample);
  occurrence <- decomp$occurrence;
  demand <- decomp$demand;
  shifted_demand <- demand-1;

  # Remove initial zeros
  if (rmzeros == TRUE){
    new_start <- which(y_insample > 0)[1]
    y_insample <- y_insample[new_start:L];
    occurrence <- occurrence[new_start:L];
  }

  #Select the distributions to be evaluated
  if (distr == "auto" || distr == "mixture"){
    to_eval <- c("nbinom", "pois", "hsnb", "hsp");
  } else{
    to_eval <- distr;
  }
  aic <- mles <- c();

  # Fit distribution with computable MLE
  if ("pois" %in% to_eval){
    lambda <- mean(y_insample);
    loglik <- sum(dpois(y_insample, lambda, log=TRUE))
    npar <- 1;
    aic["pois"] <- -2*loglik + 2*npar;
    mles["pois_lambda"] <- lambda;
  }
  if ("hsp" %in% to_eval){
    lambda <- mean(shifted_demand);
    pzero <- mean(1-occurrence);
    loglik <- sum(dhsp(y_insample, pzero, lambda, log=TRUE));
    npar <- 2;
    aic["zip"] <- -2*loglik + 2*npar;
    mles[c("zip_lambda", "zip_pzero")] <- c(lambda, pzero);
  }

  # Fit distributions with no closed form expression for ML estimators
  if ("nbinom" %in% to_eval){
    nbinom_fitted <- fit_nbinom(y_insample);
    size <- nbinom_fitted['size'];
    prob <- nbinom_fitted['prob'];
    loglik <- sum(dnbinom(y_insample, size, prob, log=TRUE));
    npar <- 2;
    aic["nbinom"] <- -2*loglik + 2*npar;
    mles[c("nbinom_size", "nbinom_prob")] <- c(size, prob);
  }
  if ("zinb" %in% to_eval){
    pzero <- mean(1-occurrence);
    nbinom_fitted <- fit_nbinom(shifted_demand);
    size <- nbinom_fitted['size'];
    prob <- nbinom_fitted['prob'];
    loglik <- sum(dhsnb(y_insample, pzero, size, prob, log=TRUE));
    npar <- 3;
    aic["zinb"] <- -2*loglik + 2*npar;
    mles[c("zinb_pzero", "zinb_size", "zinb_prob")] <- c(pzero, size, prob);
  }

  # Select the best_model, and replicate its result
  if (distr == "mixture"){
    pred_distr <- to_eval;
  } else{
    pred_distr <- names(which.min(aic));
  }
  iid_sim <- h*nsim/length(pred_distr);
  samples <- c()

  # Generate samples and replicate them for h steps
  if ("pois" %in% pred_distr){
    samples <- c(samples, rpois(iid_sim, mles["pois_lambda"]));
  }
  if ("zip" %in% pred_distr){
    samples <- c(samples, rhsp(iid_sim, mles["zip_pzero"], mles["zip_lambda"]));
  }
  if ("nbinom" %in% pred_distr){
    samples <- c(samples, rnbinom(iid_sim, mles["nbinom_size"],
                                 mles["nbinom_prob"]));
  }
  if ("zinb" %in% pred_distr){
    samples <- c(samples, rhsnb(iid_sim, mles["zinb_pzero"], mles["zinb_size"],
                                   mles["zinb_prob"]));
  }
  if (distr == "mixture"){
    samples <- samples[sample.int(length(samples))];
  }

  forecast_samples <- matrix(samples, ncol=h);


  # Save relevant parameters
  params = list("aic" = aic, "mles" = mles, "pred_distr" = pred_distr);

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
  model <- list("model" = "StaticD", "call" = call, "params" = params,
                "meanforecast" = meanforecast, "probzero" = probzero,
                "levels_upper" = levels_upper, levels_lower = "levels_lower",
                "upperCI" = upperCI, "lowerCI" = lowerCI, "samples" = samples);
  return(structure(model, class="#TODO"))
}
