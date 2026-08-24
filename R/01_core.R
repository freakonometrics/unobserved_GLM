# Core functions for audit-calibrated fair GLMs.
# Version 5.0 -- base R only.

PROXY_FAIR_GLM_CORE_VERSION <- "5.1"

`%||%` <- function(x, y) if (is.null(x)) y else x

clip_prob <- function(x, eps = 1e-6) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

log1pexp <- function(x) {
  x <- as.numeric(x)
  out <- numeric(length(x))
  positive <- x > 0
  out[positive] <- x[positive] + log1p(exp(-x[positive]))
  out[!positive] <- log1p(exp(x[!positive]))
  out
}

log_loss <- function(y, p, eps = 1e-12) {
  y <- as.numeric(y)
  p <- clip_prob(p, eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

auc_rank <- function(y, score) {
  y <- as.integer(y)
  score <- as.numeric(score)
  keep <- is.finite(score) & !is.na(y)
  y <- y[keep]
  score <- score[keep]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

brier_score <- function(y, p) {
  mean((as.numeric(y) - clip_prob(p))^2)
}

expected_calibration_error <- function(y, p, n_bins = 10L) {
  y <- as.numeric(y)
  p <- clip_prob(p)
  breaks <- unique(stats::quantile(
    p,
    probs = seq(0, 1, length.out = n_bins + 1L),
    na.rm = TRUE,
    type = 8
  ))
  if (length(breaks) < 3L) return(abs(mean(y) - mean(p)))
  groups <- cut(p, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  pieces <- split(seq_along(p), groups)
  sum(vapply(pieces, function(idx) {
    length(idx) / length(p) * abs(mean(y[idx]) - mean(p[idx]))
  }, numeric(1L)))
}

signed_group_gap <- function(score, s, subset = rep(TRUE, length(score))) {
  score <- as.numeric(score)
  s <- as.integer(s)
  keep <- as.logical(subset) & is.finite(score) & !is.na(s)
  n1 <- sum(keep & s == 1L)
  n0 <- sum(keep & s == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  mean(score[keep & s == 1L]) - mean(score[keep & s == 0L])
}

fairness_metrics <- function(y, p, s, threshold = 0.5) {
  y <- as.integer(y)
  p <- as.numeric(p)
  s <- as.integer(s)
  decision <- as.integer(p >= threshold)

  dp_score_signed <- signed_group_gap(p, s)
  eo_score_signed <- signed_group_gap(p, s, subset = y == 1L)
  dp_decision_signed <- signed_group_gap(decision, s)
  eo_decision_signed <- signed_group_gap(decision, s, subset = y == 1L)

  data.frame(
    log_loss = log_loss(y, p),
    auc = auc_rank(y, p),
    accuracy = mean(decision == y),
    dp_score_signed = dp_score_signed,
    eo_score_signed = eo_score_signed,
    dp_score = abs(dp_score_signed),
    eo_score = abs(eo_score_signed),
    dp_decision_signed = dp_decision_signed,
    eo_decision_signed = eo_decision_signed,
    dp_decision = abs(dp_decision_signed),
    eo_decision = abs(eo_decision_signed)
  )
}

stratified_split <- function(
    y,
    s,
    proportions = c(train = 0.6, validation = 0.2, test = 0.2),
    seed = 1L) {
  stopifnot(abs(sum(proportions) - 1) < 1e-10)
  if (!all(c("train", "validation", "test") %in% names(proportions))) {
    stop("proportions must be named train, validation, and test.")
  }
  set.seed(seed)
  strata <- interaction(y, s, drop = TRUE)
  split_label <- rep(NA_character_, length(y))
  for (level in levels(strata)) {
    idx <- sample(which(strata == level))
    n <- length(idx)
    n_train <- floor(proportions[["train"]] * n)
    n_validation <- floor(proportions[["validation"]] * n)
    if (n_train > 0L) split_label[idx[seq_len(n_train)]] <- "train"
    if (n_validation > 0L) {
      split_label[idx[n_train + seq_len(n_validation)]] <- "validation"
    }
    used <- n_train + n_validation
    if (used < n) split_label[idx[(used + 1L):n]] <- "test"
  }
  factor(split_label, levels = c("train", "validation", "test"))
}

standardize_from_train <- function(X, train_index, intercept = TRUE) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  intercept_column <- if (intercept && !is.null(colnames(X))) {
    which(colnames(X) == "(Intercept)")
  } else {
    integer(0L)
  }
  scale_columns <- setdiff(seq_len(ncol(X)), intercept_column)
  center <- colMeans(X[train_index, scale_columns, drop = FALSE])
  scale <- apply(X[train_index, scale_columns, drop = FALSE], 2L, stats::sd)
  scale[!is.finite(scale) | scale < 1e-12] <- 1
  X[, scale_columns] <- sweep(X[, scale_columns, drop = FALSE], 2L, center, "-")
  X[, scale_columns] <- sweep(X[, scale_columns, drop = FALSE], 2L, scale, "/")
  list(
    X = X,
    center = center,
    scale = scale,
    scaled_columns = colnames(X)[scale_columns]
  )
}

ridge_logistic_fit <- function(
    X,
    y,
    ridge = 1e-3,
    weights = NULL,
    beta_start = NULL,
    maxit = 100L,
    tol = 1e-8) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  if (length(y) != n) stop("X and y have incompatible dimensions.")

  if (is.null(weights)) weights <- rep(1, n)
  weights <- as.numeric(weights)
  if (length(weights) != n || any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be finite and nonnegative, with one value per row.")
  }
  if (mean(weights) <= 0) stop("At least one weight must be positive.")
  weights <- weights / mean(weights)

  beta <- if (is.null(beta_start)) rep(0, p) else as.numeric(beta_start)
  if (length(beta) != p) stop("beta_start has the wrong length.")
  if (is.null(beta_start)) {
    beta[1L] <- qlogis(clip_prob(stats::weighted.mean(y, weights), 1e-5))
  }

  penalty <- rep(ridge, p)
  intercept_index <- if (!is.null(colnames(X)) && "(Intercept)" %in% colnames(X)) {
    which(colnames(X) == "(Intercept)")[1L]
  } else {
    1L
  }
  penalty[intercept_index] <- 0

  objective <- function(b) {
    eta <- drop(X %*% b)
    mean(weights * (log1pexp(eta) - y * eta)) + 0.5 * sum(penalty * b^2)
  }

  converged <- FALSE
  iteration <- 0L
  for (iteration in seq_len(maxit)) {
    eta <- drop(X %*% beta)
    probability <- clip_prob(plogis(eta), 1e-8)
    gradient <- drop(crossprod(X, weights * (probability - y))) / n + penalty * beta
    working_weight <- weights * probability * (1 - probability)
    hessian <- crossprod(X, X * working_weight) / n + diag(penalty + 1e-10, p)
    step <- tryCatch(
      solve(hessian, gradient),
      error = function(e) qr.solve(hessian, gradient, tol = 1e-10)
    )

    old_objective <- objective(beta)
    step_size <- 1
    repeat {
      candidate <- beta - step_size * step
      new_objective <- objective(candidate)
      armijo <- old_objective - 1e-4 * step_size * sum(gradient * step)
      if (is.finite(new_objective) && new_objective <= armijo) break
      step_size <- step_size / 2
      if (step_size < 2^-20) {
        candidate <- beta
        break
      }
    }

    if (max(abs(candidate - beta)) < tol * (1 + max(abs(beta)))) {
      beta <- candidate
      converged <- TRUE
      break
    }
    beta <- candidate
  }

  eta <- drop(X %*% beta)
  probability <- clip_prob(plogis(eta), 1e-8)
  gradient <- drop(crossprod(X, weights * (probability - y))) / n + penalty * beta

  list(
    coefficients = beta,
    converged = converged,
    iterations = iteration,
    objective = objective(beta),
    gradient_norm = sqrt(sum(gradient^2))
  )
}

predict_ridge_logistic <- function(object, newx) {
  clip_prob(plogis(drop(as.matrix(newx) %*% object$coefficients)), 1e-5)
}

make_folds <- function(n, K = 5L, seed = 1L) {
  if (K < 2L) stop("K must be at least 2.")
  set.seed(seed)
  sample(rep(seq_len(K), length.out = n))
}

crossfit_sensitive_proxy <- function(
    X_proxy,
    s,
    r,
    K = 5L,
    ridge = 1e-2,
    seed = 1L,
    X_new = NULL,
    fallback_probability = 0.5) {
  X_proxy <- as.matrix(X_proxy)
  storage.mode(X_proxy) <- "double"
  s <- as.integer(s)
  r <- as.integer(r)
  n <- nrow(X_proxy)
  if (length(s) != n || length(r) != n) stop("X_proxy, s, and r are incompatible.")
  if (!all(r %in% c(0L, 1L))) stop("r must be binary.")

  folds <- make_folds(n, K = K, seed = seed)
  q <- rep(NA_real_, n)
  minimum_fit_n <- max(20L, ncol(X_proxy) + 2L)

  for (k in seq_len(K)) {
    fit_index <- which(folds != k & r == 1L)
    prediction_index <- which(folds == k)

    # The fallback must not use sensitive labels from the held-out fold.
    # This preserves the unitwise exclusion condition required for exact
    # design unbiasedness under an independent Bernoulli audit.
    if (length(fit_index) < minimum_fit_n || length(unique(s[fit_index])) < 2L) {
      fold_prevalence <- if (length(fit_index) > 0L) {
        mean(s[fit_index])
      } else {
        fallback_probability
      }
      q[prediction_index] <- clip_prob(fold_prevalence, 1e-4)
    } else {
      fit <- ridge_logistic_fit(
        X_proxy[fit_index, , drop = FALSE],
        s[fit_index],
        ridge = ridge
      )
      q[prediction_index] <- predict_ridge_logistic(
        fit,
        X_proxy[prediction_index, , drop = FALSE]
      )
    }
  }

  audited_index <- which(r == 1L)
  if (length(audited_index) < minimum_fit_n ||
      length(unique(s[audited_index])) < 2L) {
    full_prevalence <- if (length(audited_index) > 0L) {
      mean(s[audited_index])
    } else {
      fallback_probability
    }
    full_fit <- NULL
    q_new <- if (is.null(X_new)) NULL else {
      rep(clip_prob(full_prevalence, 1e-4), nrow(X_new))
    }
  } else {
    full_fit <- ridge_logistic_fit(
      X_proxy[audited_index, , drop = FALSE],
      s[audited_index],
      ridge = ridge
    )
    q_new <- if (is.null(X_new)) NULL else {
      predict_ridge_logistic(full_fit, X_new)
    }
  }

  list(
    q = clip_prob(q, 1e-4),
    q_new = if (is.null(q_new)) NULL else clip_prob(q_new, 1e-4),
    full_fit = full_fit,
    folds = folds,
    audited_n = length(audited_index)
  )
}

fit_sensitive_proxy_from_audit <- function(
    X_proxy,
    s,
    r,
    ridge = 1e-2,
    X_new = X_proxy,
    fallback_probability = 0.5) {
  audited_index <- which(as.integer(r) == 1L)
  minimum_fit_n <- max(20L, ncol(X_proxy) + 2L)
  if (length(audited_index) < minimum_fit_n ||
      length(unique(s[audited_index])) < 2L) {
    prevalence <- if (length(audited_index) > 0L) {
      mean(s[audited_index])
    } else {
      fallback_probability
    }
    return(list(
      q = rep(clip_prob(prevalence, 1e-4), nrow(X_new)),
      fit = NULL,
      audited_n = length(audited_index)
    ))
  }
  fit <- ridge_logistic_fit(
    X_proxy[audited_index, , drop = FALSE],
    s[audited_index],
    ridge = ridge
  )
  list(
    q = predict_ridge_logistic(fit, X_new),
    fit = fit,
    audited_n = length(audited_index)
  )
}

pseudo_sensitive <- function(s, r, pi, q = 0) {
  s <- as.numeric(s)
  r <- as.numeric(r)
  n <- length(s)
  if (length(r) != n) stop("s and r have incompatible lengths.")
  if (length(pi) == 1L) pi <- rep(pi, n)
  if (length(q) == 1L) q <- rep(q, n)
  if (length(pi) != n || length(q) != n) stop("pi and q must be scalar or length n.")
  if (any(!is.finite(pi)) || any(pi <= 0) || any(pi > 1)) {
    stop("All audit probabilities must lie in (0, 1].")
  }
  q + r / pi * (s - q)
}

estimate_score_moment <- function(score, z) {
  score <- as.numeric(score)
  z <- as.numeric(z)
  if (length(score) != length(z)) stop("score and z have incompatible lengths.")
  mean((z - mean(z)) * score)
}

estimate_score_dp <- function(score, z, prevalence_floor = 0.02) {
  alpha <- mean(as.numeric(z))
  alpha_stable <- pmin(pmax(alpha, prevalence_floor), 1 - prevalence_floor)
  signed_dp <- estimate_score_moment(score, z) / (alpha_stable * (1 - alpha_stable))
  c(signed = signed_dp, absolute = abs(signed_dp), prevalence = alpha)
}

exact_design_variance <- function(score, s, q, pi) {
  score <- as.numeric(score)
  s <- as.numeric(s)
  q <- as.numeric(q)
  n <- length(score)
  if (length(pi) == 1L) pi <- rep(pi, n)
  centered_score <- score - mean(score)
  sum((1 - pi) / pi * (s - q)^2 * centered_score^2) / n^2
}

estimate_design_variance <- function(score, s, r, q, pi) {
  score <- as.numeric(score)
  s <- as.numeric(s)
  r <- as.numeric(r)
  q <- as.numeric(q)
  n <- length(score)
  if (length(pi) == 1L) pi <- rep(pi, n)
  centered_score <- score - mean(score)
  sum(r * (1 - pi) / pi^2 * (s - q)^2 * centered_score^2) / n^2
}

fairness_direction <- function(X, z) {
  X <- as.matrix(X)
  centered_z <- as.numeric(z) - mean(z)
  drop(crossprod(X, centered_z)) / nrow(X)
}

fair_logistic_objective <- function(beta, X, y, cvec, lambda, ridge) {
  eta <- drop(X %*% beta)
  penalty <- rep(ridge, length(beta))
  intercept_index <- if (!is.null(colnames(X)) && "(Intercept)" %in% colnames(X)) {
    which(colnames(X) == "(Intercept)")[1L]
  } else {
    1L
  }
  penalty[intercept_index] <- 0
  mean(log1pexp(eta) - y * eta) +
    0.5 * sum(penalty * beta^2) +
    lambda * sum(cvec * beta)^2
}

fit_fair_logistic <- function(
    X,
    y,
    cvec,
    lambda = 0,
    ridge = 1e-4,
    beta_start = NULL,
    maxit = 120L,
    tol = 1e-8) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  cvec <- as.numeric(cvec)
  if (length(cvec) != p) stop("cvec has the wrong length.")

  beta <- if (is.null(beta_start)) rep(0, p) else as.numeric(beta_start)
  if (is.null(beta_start)) beta[1L] <- qlogis(clip_prob(mean(y), 1e-5))

  penalty <- rep(ridge, p)
  intercept_index <- if (!is.null(colnames(X)) && "(Intercept)" %in% colnames(X)) {
    which(colnames(X) == "(Intercept)")[1L]
  } else {
    1L
  }
  penalty[intercept_index] <- 0

  converged <- FALSE
  iteration <- 0L
  for (iteration in seq_len(maxit)) {
    eta <- drop(X %*% beta)
    probability <- clip_prob(plogis(eta), 1e-8)
    surrogate_moment <- sum(cvec * beta)
    gradient <- drop(crossprod(X, probability - y)) / n +
      penalty * beta + 2 * lambda * surrogate_moment * cvec
    working_weight <- probability * (1 - probability)
    hessian <- crossprod(X, X * working_weight) / n +
      diag(penalty + 1e-10, p) + 2 * lambda * tcrossprod(cvec)
    step <- tryCatch(
      solve(hessian, gradient),
      error = function(e) qr.solve(hessian, gradient, tol = 1e-10)
    )

    old_objective <- fair_logistic_objective(beta, X, y, cvec, lambda, ridge)
    step_size <- 1
    repeat {
      candidate <- beta - step_size * step
      new_objective <- fair_logistic_objective(candidate, X, y, cvec, lambda, ridge)
      armijo <- old_objective - 1e-4 * step_size * sum(gradient * step)
      if (is.finite(new_objective) && new_objective <= armijo) break
      step_size <- step_size / 2
      if (step_size < 2^-20) {
        candidate <- beta
        break
      }
    }

    if (max(abs(candidate - beta)) < tol * (1 + max(abs(beta)))) {
      beta <- candidate
      converged <- TRUE
      break
    }
    beta <- candidate
  }

  eta <- drop(X %*% beta)
  probability <- clip_prob(plogis(eta), 1e-8)
  surrogate_moment <- sum(cvec * beta)
  gradient <- drop(crossprod(X, probability - y)) / n +
    penalty * beta + 2 * lambda * surrogate_moment * cvec

  list(
    coefficients = beta,
    converged = converged,
    convergence_code = if (converged) 0L else 1L,
    iterations = iteration,
    lambda = lambda,
    surrogate_moment = surrogate_moment,
    objective = fair_logistic_objective(beta, X, y, cvec, lambda, ridge),
    gradient_norm = sqrt(sum(gradient^2)),
    path_type = "convex"
  )
}

fit_fair_path <- function(X, y, z, lambdas, ridge = 1e-4, beta_start = NULL) {
  cvec <- fairness_direction(X, z)
  lambdas <- sort(unique(as.numeric(lambdas)))
  fits <- vector("list", length(lambdas))
  beta <- beta_start
  for (j in seq_along(lambdas)) {
    fits[[j]] <- fit_fair_logistic(
      X = X,
      y = y,
      cvec = cvec,
      lambda = lambdas[j],
      ridge = ridge,
      beta_start = beta
    )
    beta <- fits[[j]]$coefficients
  }
  names(fits) <- format(lambdas, scientific = TRUE)
  list(
    fits = fits,
    lambdas = lambdas,
    cvec = cvec,
    z = as.numeric(z),
    path_type = "convex"
  )
}

exact_fair_logistic_objective <- function(beta, X, y, z, lambda, ridge) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  centered_z <- as.numeric(z) - mean(z)
  eta <- drop(X %*% beta)
  probability <- clip_prob(plogis(eta), 1e-8)
  penalty <- rep(ridge, length(beta))
  intercept_index <- if (!is.null(colnames(X)) && "(Intercept)" %in% colnames(X)) {
    which(colnames(X) == "(Intercept)")[1L]
  } else {
    1L
  }
  penalty[intercept_index] <- 0
  moment <- mean(centered_z * probability)
  mean(log1pexp(eta) - y * eta) +
    0.5 * sum(penalty * beta^2) +
    lambda * moment^2
}

exact_fair_logistic_gradient <- function(beta, X, y, z, lambda, ridge) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  centered_z <- as.numeric(z) - mean(z)
  n <- nrow(X)
  eta <- drop(X %*% beta)
  probability <- clip_prob(plogis(eta), 1e-8)
  penalty <- rep(ridge, length(beta))
  intercept_index <- if (!is.null(colnames(X)) && "(Intercept)" %in% colnames(X)) {
    which(colnames(X) == "(Intercept)")[1L]
  } else {
    1L
  }
  penalty[intercept_index] <- 0
  moment <- mean(centered_z * probability)
  gradient_moment <- drop(crossprod(
    X,
    centered_z * probability * (1 - probability)
  )) / n
  drop(crossprod(X, probability - y)) / n +
    penalty * beta + 2 * lambda * moment * gradient_moment
}

fit_exact_fair_logistic <- function(
    X,
    y,
    z,
    lambda = 0,
    ridge = 1e-4,
    beta_start = NULL,
    maxit = 300L,
    reltol = 1e-9,
    gradient_tolerance = 1e-5) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.numeric(y)
  z <- as.numeric(z)
  p <- ncol(X)
  if (is.null(beta_start)) {
    beta_start <- rep(0, p)
    beta_start[1L] <- qlogis(clip_prob(mean(y), 1e-5))
  }

  run_optim <- function(start) {
    stats::optim(
      par = as.numeric(start),
      fn = exact_fair_logistic_objective,
      gr = exact_fair_logistic_gradient,
      X = X,
      y = y,
      z = z,
      lambda = lambda,
      ridge = ridge,
      method = "BFGS",
      control = list(maxit = maxit, reltol = reltol)
    )
  }

  candidate <- run_optim(beta_start)
  gradient <- exact_fair_logistic_gradient(
    candidate$par, X, y, z, lambda, ridge
  )
  gradient_norm <- sqrt(sum(gradient^2))

  # A zero-vector fallback is only used when the warm-start solution reports
  # poor convergence. The better objective is retained.
  if (candidate$convergence != 0L || !is.finite(gradient_norm) ||
      gradient_norm > gradient_tolerance) {
    fallback_start <- rep(0, p)
    fallback_start[1L] <- qlogis(clip_prob(mean(y), 1e-5))
    fallback <- run_optim(fallback_start)
    if (is.finite(fallback$value) && fallback$value < candidate$value) {
      candidate <- fallback
      gradient <- exact_fair_logistic_gradient(
        candidate$par, X, y, z, lambda, ridge
      )
      gradient_norm <- sqrt(sum(gradient^2))
    }
  }

  eta <- drop(X %*% candidate$par)
  probability <- clip_prob(plogis(eta), 1e-8)
  centered_z <- z - mean(z)

  list(
    coefficients = candidate$par,
    converged = identical(candidate$convergence, 0L) &&
      is.finite(gradient_norm) && gradient_norm <= 10 * gradient_tolerance,
    convergence_code = candidate$convergence,
    iterations = unname(candidate$counts[["function"]]),
    lambda = lambda,
    exact_score_moment = mean(centered_z * probability),
    objective = candidate$value,
    gradient_norm = gradient_norm,
    path_type = "exact"
  )
}

fit_exact_fair_path <- function(X, y, z, lambdas, ridge = 1e-4, beta_start = NULL) {
  lambdas <- sort(unique(as.numeric(lambdas)))
  fits <- vector("list", length(lambdas))
  beta <- beta_start
  for (j in seq_along(lambdas)) {
    fits[[j]] <- fit_exact_fair_logistic(
      X = X,
      y = y,
      z = z,
      lambda = lambdas[j],
      ridge = ridge,
      beta_start = beta
    )
    beta <- fits[[j]]$coefficients
  }
  names(fits) <- format(lambdas, scientific = TRUE)
  list(
    fits = fits,
    lambdas = lambdas,
    z = as.numeric(z),
    path_type = "exact"
  )
}

fit_fair_path_by_type <- function(
    path_type,
    X,
    y,
    z,
    lambdas,
    ridge = 1e-4,
    beta_start = NULL) {
  path_type <- match.arg(path_type, c("exact", "convex"))
  if (path_type == "exact") {
    fit_exact_fair_path(X, y, z, lambdas, ridge, beta_start)
  } else {
    fit_fair_path(X, y, z, lambdas, ridge, beta_start)
  }
}

predict_fair_path <- function(path, X) {
  X <- as.matrix(X)
  predictions <- vapply(path$fits, function(fit) {
    clip_prob(plogis(drop(X %*% fit$coefficients)), 1e-8)
  }, numeric(nrow(X)))
  if (is.vector(predictions)) predictions <- matrix(predictions, ncol = 1L)
  colnames(predictions) <- format(path$lambdas, scientific = TRUE)
  predictions
}

path_diagnostics <- function(path, X, y, z) {
  predictions <- predict_fair_path(path, X)
  rows <- lapply(seq_along(path$lambdas), function(j) {
    beta <- path$fits[[j]]$coefficients
    eta <- drop(as.matrix(X) %*% beta)
    dp <- estimate_score_dp(predictions[, j], z)
    data.frame(
      lambda = path$lambdas[j],
      log_loss = log_loss(y, predictions[, j]),
      estimated_moment_signed = estimate_score_moment(predictions[, j], z),
      estimated_moment = abs(estimate_score_moment(predictions[, j], z)),
      estimated_dp_signed = unname(dp[["signed"]]),
      estimated_dp = unname(dp[["absolute"]]),
      prevalence_estimate = unname(dp[["prevalence"]]),
      surrogate_moment_signed = sum(fairness_direction(X, z) * beta),
      surrogate_moment = abs(sum(fairness_direction(X, z) * beta)),
      logit_spread = mean((eta - mean(eta))^2),
      converged = isTRUE(path$fits[[j]]$converged),
      gradient_norm = path$fits[[j]]$gradient_norm %||% NA_real_
    )
  })
  do.call(rbind, rows)
}

select_lambda_by_relative_moment_reduction <- function(
    path,
    X_select,
    y_select,
    z_select,
    reduction_fraction = 0.50) {
  if (!is.finite(reduction_fraction) || reduction_fraction < 0 ||
      reduction_fraction > 1) {
    stop("reduction_fraction must lie in [0, 1].")
  }
  table <- path_diagnostics(path, X_select, y_select, z_select)
  zero_candidates <- which(abs(table$lambda) <= sqrt(.Machine$double.eps))
  baseline_index <- if (length(zero_candidates) > 0L) {
    zero_candidates[1L]
  } else {
    which.min(table$lambda)
  }
  baseline <- table$estimated_moment[baseline_index]
  target <- reduction_fraction * baseline

  if (!is.finite(baseline) || baseline <= 1e-12) {
    chosen <- baseline_index
    feasible <- TRUE
  } else {
    feasible_indices <- which(
      is.finite(table$estimated_moment) &
        table$estimated_moment <= target + 1e-12
    )
    if (length(feasible_indices) > 0L) {
      chosen <- feasible_indices[which.min(table$lambda[feasible_indices])]
      feasible <- TRUE
    } else {
      best <- min(table$estimated_moment, na.rm = TRUE)
      candidates <- which(abs(table$estimated_moment - best) <= 1e-12)
      chosen <- candidates[which.min(table$lambda[candidates])]
      feasible <- FALSE
    }
  }

  list(
    index = chosen,
    lambda = table$lambda[chosen],
    table = table,
    baseline = baseline,
    target = target,
    reduction_fraction = reduction_fraction,
    feasible = feasible
  )
}

summarize_replicates <- function(df, group_vars, value_vars) {
  key <- interaction(df[group_vars], drop = TRUE, lex.order = TRUE)
  pieces <- split(df, key)
  output <- lapply(pieces, function(d) {
    base <- d[1L, group_vars, drop = FALSE]
    for (variable in value_vars) {
      x <- d[[variable]]
      finite <- is.finite(x)
      count <- sum(finite)
      base[[paste0(variable, "_mean")]] <- if (count > 0L) mean(x[finite]) else NA_real_
      base[[paste0(variable, "_sd")]] <- if (count > 1L) stats::sd(x[finite]) else NA_real_
      base[[paste0(variable, "_se")]] <- if (count > 1L) {
        stats::sd(x[finite]) / sqrt(count)
      } else {
        NA_real_
      }
      base[[paste0(variable, "_q025")]] <- if (count > 0L) {
        unname(stats::quantile(x[finite], 0.025, type = 8))
      } else {
        NA_real_
      }
      base[[paste0(variable, "_q975")]] <- if (count > 0L) {
        unname(stats::quantile(x[finite], 0.975, type = 8))
      } else {
        NA_real_
      }
    }
    base$n_rep <- nrow(d)
    base
  })
  rownames(output) <- NULL
  do.call(rbind, output)
}

degrade_proxy <- function(q, severity) {
  if (any(severity < 0 | severity > 1)) stop("severity must lie in [0, 1].")
  q <- clip_prob(q)
  clip_prob((1 - severity) * q + severity * (1 - q), 1e-4)
}

allocate_audit_probabilities <- function(
    variance_score,
    budget_rate,
    pi_min = NULL,
    pi_max = 1) {
  variance_score <- pmax(as.numeric(variance_score), 0)
  n <- length(variance_score)
  if (!is.finite(budget_rate) || budget_rate <= 0 || budget_rate > 1) {
    stop("budget_rate must lie in (0, 1].")
  }
  if (is.null(pi_min)) pi_min <- min(1e-3, budget_rate / 20)
  if (pi_min <= 0 || pi_min > budget_rate) {
    stop("pi_min must be positive and no larger than budget_rate.")
  }
  target <- n * budget_rate
  if (target < n * pi_min - 1e-10 || target > n * pi_max + 1e-10) {
    stop("The requested budget is incompatible with pi_min and pi_max.")
  }

  root_score <- sqrt(variance_score + 1e-16)
  if (max(root_score) <= 1e-12) return(rep(budget_rate, n))

  total_probability <- function(multiplier) {
    sum(pmin(pi_max, pmax(pi_min, multiplier * root_score)))
  }

  lower <- 0
  upper <- 1
  while (total_probability(upper) < target) upper <- upper * 2
  multiplier <- stats::uniroot(
    function(value) total_probability(value) - target,
    interval = c(lower, upper),
    tol = 1e-12
  )$root
  probability <- pmin(pi_max, pmax(pi_min, multiplier * root_score))

  # Correct tiny numerical budget discrepancies without violating bounds.
  discrepancy <- target - sum(probability)
  if (abs(discrepancy) > 1e-8) {
    free <- which(probability > pi_min + 1e-12 & probability < pi_max - 1e-12)
    if (length(free) > 0L) probability[free] <- probability[free] + discrepancy / length(free)
    probability <- pmin(pi_max, pmax(pi_min, probability))
  }
  probability
}

draw_audit <- function(pi, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stats::rbinom(length(pi), size = 1L, prob = pi)
}

proxy_diagnostics <- function(s, q) {
  data.frame(
    auc = auc_rank(s, q),
    brier = brier_score(s, q),
    ece = expected_calibration_error(s, q),
    mean_probability = mean(q),
    prevalence = mean(s)
  )
}

# -----------------------------------------------------------------------------
# Version 5.1 additions
# -----------------------------------------------------------------------------

# Select the smallest lambda that attains a relative reduction in either the
# exact probability-score moment or the convex logit surrogate moment.
select_lambda_by_relative_diagnostic_reduction <- function(
    path,
    X_select,
    y_select,
    z_select,
    reduction_fraction = 0.50,
    selection_metric = c("exact", "surrogate")) {
  selection_metric <- match.arg(selection_metric)
  if (!is.finite(reduction_fraction) || reduction_fraction < 0 ||
      reduction_fraction > 1) {
    stop("reduction_fraction must lie in [0, 1].")
  }

  table <- path_diagnostics(path, X_select, y_select, z_select)
  metric_column <- if (selection_metric == "exact") {
    "estimated_moment"
  } else {
    "surrogate_moment"
  }
  metric_values <- table[[metric_column]]

  zero_candidates <- which(abs(table$lambda) <= sqrt(.Machine$double.eps))
  baseline_index <- if (length(zero_candidates) > 0L) {
    zero_candidates[1L]
  } else {
    which.min(table$lambda)
  }
  baseline <- metric_values[baseline_index]
  target <- reduction_fraction * baseline

  if (!is.finite(baseline) || baseline <= 1e-12) {
    chosen <- baseline_index
    feasible <- TRUE
  } else {
    feasible_indices <- which(
      is.finite(metric_values) & metric_values <= target + 1e-12
    )
    if (length(feasible_indices) > 0L) {
      chosen <- feasible_indices[which.min(table$lambda[feasible_indices])]
      feasible <- TRUE
    } else {
      best <- min(metric_values, na.rm = TRUE)
      candidates <- which(abs(metric_values - best) <= 1e-12)
      chosen <- candidates[which.min(table$lambda[candidates])]
      feasible <- FALSE
    }
  }

  list(
    index = chosen,
    lambda = table$lambda[chosen],
    table = table,
    baseline = baseline,
    target = target,
    reduction_fraction = reduction_fraction,
    feasible = feasible,
    selection_metric = selection_metric,
    metric_column = metric_column
  )
}

# Backward-compatible wrapper: the main application selects on the exact
# probability-score moment.
select_lambda_by_relative_moment_reduction <- function(
    path,
    X_select,
    y_select,
    z_select,
    reduction_fraction = 0.50) {
  select_lambda_by_relative_diagnostic_reduction(
    path = path,
    X_select = X_select,
    y_select = y_select,
    z_select = z_select,
    reduction_fraction = reduction_fraction,
    selection_metric = "exact"
  )
}

# One-sided uncertainty-aware selector for an independent validation audit.
# The criterion is U_j = |M_j| + z_alpha * se(M_j), and the selected lambda is
# the smallest one satisfying U_j <= reduction_fraction * U_0. This is a
# conservative stability diagnostic, not a formal simultaneous confidence
# guarantee over the full path.
select_lambda_by_relative_moment_ucb <- function(
    path,
    X_select,
    y_select,
    z_select,
    s_select,
    r_select,
    pi_select,
    q_select,
    reduction_fraction = 0.50,
    z_value = 1.645) {
  if (!is.finite(z_value) || z_value < 0) stop("z_value must be nonnegative.")
  table <- path_diagnostics(path, X_select, y_select, z_select)
  predictions <- predict_fair_path(path, X_select)

  moment_se <- vapply(seq_along(path$lambdas), function(j) {
    variance <- estimate_design_variance(
      score = predictions[, j],
      s = s_select,
      r = r_select,
      q = q_select,
      pi = pi_select
    )
    sqrt(max(variance, 0))
  }, numeric(1L))

  table$estimated_moment_se <- moment_se
  table$estimated_moment_ucb <- table$estimated_moment + z_value * moment_se

  zero_candidates <- which(abs(table$lambda) <= sqrt(.Machine$double.eps))
  baseline_index <- if (length(zero_candidates) > 0L) {
    zero_candidates[1L]
  } else {
    which.min(table$lambda)
  }
  baseline <- table$estimated_moment[baseline_index]
  baseline_ucb <- table$estimated_moment_ucb[baseline_index]
  target_ucb <- reduction_fraction * baseline_ucb

  if (!is.finite(baseline_ucb) || baseline_ucb <= 1e-12) {
    chosen <- baseline_index
    feasible <- TRUE
  } else {
    feasible_indices <- which(
      is.finite(table$estimated_moment_ucb) &
        table$estimated_moment_ucb <= target_ucb + 1e-12
    )
    if (length(feasible_indices) > 0L) {
      chosen <- feasible_indices[which.min(table$lambda[feasible_indices])]
      feasible <- TRUE
    } else {
      best <- min(table$estimated_moment_ucb, na.rm = TRUE)
      candidates <- which(abs(table$estimated_moment_ucb - best) <= 1e-12)
      chosen <- candidates[which.min(table$lambda[candidates])]
      feasible <- FALSE
    }
  }

  list(
    index = chosen,
    lambda = table$lambda[chosen],
    table = table,
    baseline = baseline,
    baseline_ucb = baseline_ucb,
    target = reduction_fraction * baseline,
    target_ucb = target_ucb,
    reduction_fraction = reduction_fraction,
    feasible = feasible,
    z_value = z_value,
    selection_metric = "exact_ucb"
  )
}

# Fixed-size pilot sample, useful when pilot cost must be counted exactly.
draw_fixed_size_audit <- function(n, m, seed = NULL) {
  n <- as.integer(n)
  m <- as.integer(round(m))
  if (n <= 0L || m < 0L || m > n) stop("Require 0 <= m <= n.")
  if (!is.null(seed)) set.seed(seed)
  r <- integer(n)
  if (m > 0L) r[sample.int(n, size = m, replace = FALSE)] <- 1L
  r
}

# Two-phase estimator: pilot units are observed exactly; among non-pilot units,
# residuals are corrected using a second-stage Bernoulli audit.
two_phase_pseudo_sensitive <- function(s, pilot_r, main_r, main_pi, q) {
  s <- as.numeric(s)
  pilot_r <- as.integer(pilot_r)
  main_r <- as.integer(main_r)
  q <- as.numeric(q)
  n <- length(s)
  if (length(pilot_r) != n || length(main_r) != n || length(q) != n) {
    stop("s, pilot_r, main_r, and q must have the same length.")
  }
  if (length(main_pi) == 1L) main_pi <- rep(main_pi, n)
  if (length(main_pi) != n) stop("main_pi must be scalar or length n.")

  pilot <- pilot_r == 1L
  nonpilot <- !pilot
  if (any(main_r[pilot] != 0L)) stop("main_r must be zero for pilot units.")
  if (any(main_pi[nonpilot] <= 0 | main_pi[nonpilot] > 1)) {
    stop("Second-stage probabilities for non-pilot units must lie in (0,1].")
  }

  z <- q
  z[pilot] <- s[pilot]
  z[nonpilot] <- q[nonpilot] +
    main_r[nonpilot] / main_pi[nonpilot] * (s[nonpilot] - q[nonpilot])
  z
}

two_phase_exact_design_variance <- function(
    score, s, pilot_r, q, main_pi) {
  score <- as.numeric(score)
  s <- as.numeric(s)
  pilot_r <- as.integer(pilot_r)
  q <- as.numeric(q)
  n <- length(score)
  if (length(main_pi) == 1L) main_pi <- rep(main_pi, n)
  nonpilot <- pilot_r == 0L
  centered_score <- score - mean(score)
  sum(
    (1 - main_pi[nonpilot]) / main_pi[nonpilot] *
      (s[nonpilot] - q[nonpilot])^2 * centered_score[nonpilot]^2
  ) / n^2
}

two_phase_estimate_design_variance <- function(
    score, s, pilot_r, main_r, q, main_pi) {
  score <- as.numeric(score)
  s <- as.numeric(s)
  pilot_r <- as.integer(pilot_r)
  main_r <- as.integer(main_r)
  q <- as.numeric(q)
  n <- length(score)
  if (length(main_pi) == 1L) main_pi <- rep(main_pi, n)
  nonpilot <- pilot_r == 0L
  centered_score <- score - mean(score)
  sum(
    main_r[nonpilot] * (1 - main_pi[nonpilot]) / main_pi[nonpilot]^2 *
      (s[nonpilot] - q[nonpilot])^2 * centered_score[nonpilot]^2
  ) / n^2
}

# Delta-method diagnostic for the Monte Carlo standard error of the sample
# variance, allowing non-Gaussian estimator draws through their fourth moment.
sample_variance_standard_error <- function(x) {
  x <- as.numeric(x[is.finite(x)])
  n <- length(x)
  if (n < 4L) return(NA_real_)
  centered <- x - mean(x)
  s2 <- stats::var(x)
  mu4 <- mean(centered^4)
  variance_of_s2 <- (mu4 - ((n - 3) / (n - 1)) * s2^2) / n
  sqrt(max(variance_of_s2, 0))
}

