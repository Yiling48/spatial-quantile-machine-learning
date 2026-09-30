# Geographically Weighted Quantile XGBoost
#
# Cleaned portfolio version of the thesis implementation.
# The function fits a local XGBoost model around each observation using
# spatial kernel weights.

gw_qxgb <- function(
  formula,
  cordx,
  cordy,
  h,
  kernel = c("bisquare", "gaussian", "global"),
  data,
  tau = 0.5,
  dist_km = FALSE,
  params = list(
    objective = "reg:quantileerror",
    quantile_alpha = 0.5,
    eta = 0.05,
    max_depth = 6,
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 3
  ),
  nrounds = 200
) {
  kernel <- match.arg(kernel)

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required.")
  }

  model_frame <- model.frame(formula, data = data)
  y <- model_frame[, 1]
  X <- model.matrix(formula, data = data)[, -1, drop = FALSE]

  n <- nrow(X)
  if (length(cordx) != n || length(cordy) != n) {
    stop("Coordinate lengths must match the number of observations.")
  }

  params$objective <- "reg:quantileerror"
  params$quantile_alpha <- tau

  yhat <- numeric(n)

  for (i in seq_len(n)) {
    d <- spatial_distance(
      cordx = cordx,
      cordy = cordy,
      i = i,
      dist_km = dist_km
    )

    weights <- spatial_kernel_weights(
      d = d,
      h = h,
      kernel = kernel
    )

    idx <- which(weights > 0)
    X_local <- X[idx, , drop = FALSE]
    y_local <- y[idx]
    w_local <- weights[idx]

    dtrain <- xgboost::xgb.DMatrix(
      data = X_local,
      label = y_local,
      weight = w_local
    )

    model <- xgboost::xgb.train(
      params = params,
      data = dtrain,
      nrounds = nrounds,
      verbose = 0
    )

    yhat[i] <- predict(
      model,
      newdata = X[i, , drop = FALSE]
    )
  }

  list(
    y = y,
    yhat = yhat,
    residuals = y - yhat,
    tau = tau,
    bandwidth = h,
    kernel = kernel,
    params = params
  )
}


spatial_distance <- function(cordx, cordy, i, dist_km = FALSE) {
  if (dist_km) {
    if (!requireNamespace("geosphere", quietly = TRUE)) {
      stop("Package 'geosphere' is required when dist_km = TRUE.")
    }

    return(
      geosphere::distCosine(
        cbind(cordx, cordy),
        cbind(cordx[i], cordy[i]),
        r = 6371
      )
    )
  }

  sqrt((cordx - cordx[i])^2 + (cordy - cordy[i])^2)
}


spatial_kernel_weights <- function(
  d,
  h,
  kernel = c("bisquare", "gaussian", "global")
) {
  kernel <- match.arg(kernel)
  n <- length(d)

  if (kernel == "global") {
    return(rep(1, n))
  }

  h_index <- max(1, min(n, round(h)))
  distance_bandwidth <- sort(d)[h_index]

  if (distance_bandwidth <= 0) {
    stop("Bandwidth distance must be positive.")
  }

  if (kernel == "gaussian") {
    return(exp(-0.5 * (d / distance_bandwidth)^2))
  }

  # Adaptive bisquare kernel
  w <- numeric(n)
  inside <- which(d <= distance_bandwidth)
  w[inside] <- (1 - (d[inside] / distance_bandwidth)^2)^2
  w
}
