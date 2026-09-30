# Bandwidth selection for GW-QXGB
#
# Uses leave-one-location-out prediction and golden-section search.

cv_gw_qxgb <- function(
  formula,
  cordx,
  cordy,
  min_h,
  max_h,
  kernel = c("bisquare", "gaussian", "global"),
  data,
  tau = 0.5,
  model_type = c("quantile", "mean"),
  dist_km = FALSE,
  params = list(
    eta = 0.05,
    max_depth = 6,
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 3
  ),
  nrounds = 200,
  tolerance = 1
) {
  kernel <- match.arg(kernel)
  model_type <- match.arg(model_type)

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required.")
  }

  model_frame <- model.frame(formula, data = data)
  y <- model_frame[, 1]
  X <- model.matrix(formula, data = data)[, -1, drop = FALSE]
  n <- nrow(X)

  if (model_type == "quantile") {
    params$objective <- "reg:quantileerror"
    params$quantile_alpha <- tau
  } else {
    params$objective <- "reg:squarederror"
    params$quantile_alpha <- NULL
  }

  score_bandwidth <- function(h) {
    pred <- numeric(n)

    for (i in seq_len(n)) {
      d <- spatial_distance(cordx, cordy, i, dist_km)
      w <- spatial_kernel_weights(d, h, kernel)

      # Leave the target location out of its local training set.
      w[i] <- 0
      idx <- which(w > 0)

      dtrain <- xgboost::xgb.DMatrix(
        data = X[idx, , drop = FALSE],
        label = y[idx],
        weight = w[idx]
      )

      model <- xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = nrounds,
        verbose = 0
      )

      pred[i] <- predict(
        model,
        newdata = X[i, , drop = FALSE]
      )
    }

    error <- y - pred

    if (model_type == "quantile") {
      return(sum(error * (tau - (error < 0))))
    }

    sum(error^2)
  }

  ratio <- (sqrt(5) - 1) / 2
  left <- min_h
  right <- max_h

  x1 <- ratio * left + (1 - ratio) * right
  x2 <- (1 - ratio) * left + ratio * right
  f1 <- score_bandwidth(x1)
  f2 <- score_bandwidth(x2)

  while ((right - left) > tolerance) {
    if (f1 < f2) {
      right <- x2
      x2 <- x1
      f2 <- f1
      x1 <- ratio * left + (1 - ratio) * right
      f1 <- score_bandwidth(x1)
    } else {
      left <- x1
      x1 <- x2
      f1 <- f2
      x2 <- (1 - ratio) * left + ratio * right
      f2 <- score_bandwidth(x2)
    }
  }

  if (f1 <= f2) left else right
}
