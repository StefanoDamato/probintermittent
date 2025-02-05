#quantile loss
quantile_loss <- function(actuals, qforecast, level){
  return(2*abs(actuals-qforecast)*ifelse(actuals>qforecast, level, 1-level));
}

#ql scaled in sample loo
quantile_loss_scaled <- function(actuals, qforecast, level, y_insample){

  y_insample <- sort(y_insample);
  K <- ceiling((length(y_insample)-1)*level);
  scale = c(quantile_loss(y_insample[1:K], y_insample[K+1], level),
            quantile_loss(y_insample[(K+1):length(y_insample)], y[K], level));

  return(quantile_loss(actuals, qforecast, level)/mean(scale));
}

x = c(1,2,3,4)
quantile(ecdf(x), 0.5001, type=1)
