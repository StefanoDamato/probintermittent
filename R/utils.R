#' Croston's decomposition
#'
#' Divide a time series into occurrence and demand.
#'
#' @param y A vector of data.
#' @return The function return a list of 3 series:
#' \describe{
#'   \item{occurrence}{The binarised series between positive and negative values.}
#'   \item{demand}{The subseries of positive values.}
#'   \item{intervals}{A sequence of time intervals between positive observations.}
#' }
#'
#' @export
crostonsdecomp <- function(y){

  if (all(y>0)) stop("The time series is smooth.");
  if (all(y==0)) stop("The time series is all zero.");

  occurrence <- ifelse(y>0, 1, 0);
  d_times <- which(y > 0);
  demand <- y[d_times];
  intervals <- diff(c(0,d_times));

  return(list("occurrence" = occurrence, "demand" = demand,
              "intervals" = intervals));
}
