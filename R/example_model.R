#' Big title #TODO
#'
#' One-sentence description #TODO
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
#' test <- example_model(y) #TODO
#'
#' @importFrom package function
example_model = function(data, h=10, levels=0.9, holdout=FALSE,
                         cumulative=FALSE, side=c("upper", "both", "lower"),
                         nsim=10000){

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
    list2env(data_list, parent.frame());

    # Whatever happens in the model, generate an array of samples
    params = list();
    forecast_samples = matrix(rpois(nsim*h, 1), nsim, h);

    # Store the results in some arrays
    forecast_list <- get_forecast(forecast_samples, levels, cumulative, side, h);
    list2env(forecast_list, parent.frame());

    # Save the results as time series
    meanforecast <- ts(mean_forecast, start=pred_start ,frequency=freq);
    probzero <- ts(prob_zero, start=pred_start, frequency=freq);
    upperCI <- ts(CI_upper, start=pred_start, frequency=freq, names=levels_upper);
    lowerCI <- ts(CI_lower, start=pred_start, frequency=freq, names=levels_lower);

    # Instantiate the predictions in a class
    model <- list("model" = "#TODO", "call" = call, "params" = params,
                  "meanforecast" = meanforecast, "probzero" = probzero,
                  "levels_upper" = levels_upper, levels_lower = "levels_lower",
                  "upperCI" = upperCI, "lowerCI" = lowerCI);
    return(structure(model, class="#TODO"))
}





