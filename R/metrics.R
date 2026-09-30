pinball_loss <- function(actual, predicted, tau = 0.5) {
  error <- actual - predicted
  mean(pmax(tau * error, (tau - 1) * error))
}


pseudo_r2_quantile <- function(actual, predicted, tau = 0.5) {
  model_loss <- pinball_loss(actual, predicted, tau)

  null_prediction <- rep(
    as.numeric(stats::quantile(actual, probs = tau, na.rm = TRUE)),
    length(actual)
  )
  null_loss <- pinball_loss(actual, null_prediction, tau)

  1 - model_loss / null_loss
}
