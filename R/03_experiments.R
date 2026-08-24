# Reproduction workflows for Audit-Certified Fair GLMs.
# This file consolidates the numerical scripts used in the manuscript.
# It assumes R/01_core.R and R/02_data.R have already been sourced.

credit_config <- function(
    quick = TRUE,
    sensitive = "sex",
    path_type = "exact",
    audit_rates = c(0.01, 0.02, 0.05, 0.10, 0.20),
    n_rep_moment = NULL,
    n_rep_learning = NULL,
    k_proxy = 5L,
    ridge_proxy = 1e-2,
    ridge_outcome = 1e-4,
    lambda_grid = NULL,
    reduction_fraction = 0.50,
    seed_split = 20260801L) {

  sensitive <- match.arg(tolower(sensitive), c("sex", "age"))
  path_type <- match.arg(tolower(path_type), c("exact", "convex"))

  if (is.null(n_rep_moment)) {
    n_rep_moment <- if (quick) 20L else 300L
  }
  if (is.null(n_rep_learning)) {
    n_rep_learning <- if (quick) 3L else 40L
  }
  if (is.null(lambda_grid)) {
    lambda_grid <- c(0, 10^seq(-3, 4, length.out = if (quick) 15L else 31L))
  }

  list(
    quick = quick,
    sensitive = sensitive,
    path_type = path_type,
    audit_rates = audit_rates,
    n_rep_moment = as.integer(n_rep_moment),
    n_rep_learning = as.integer(n_rep_learning),
    k_proxy = as.integer(k_proxy),
    ridge_proxy = ridge_proxy,
    ridge_outcome = ridge_outcome,
    lambda_grid = lambda_grid,
    reduction_fraction = reduction_fraction,
    seed_split = as.integer(seed_split)
  )
}


run_smoke_simulation <- function(
    n = 2500L,
    audit_rate = 0.05,
    seed = 20260802L) {

  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  X <- stats::model.matrix(~ x1 + x2 + I(x1^2))
  s_probability <- plogis(-0.3 + 0.7 * x1 - 0.4 * x2 + 1.1 * x1^2)
  s <- stats::rbinom(n, 1L, s_probability)
  y <- stats::rbinom(
    n, 1L, plogis(-0.2 + 0.8 * x1 - 0.5 * x2 + 0.5 * s)
  )
  score <- plogis(drop(X %*% c(-0.1, 0.7, -0.4, 0.5)))
  truth <- estimate_score_moment(score, s)

  r <- draw_audit(rep(audit_rate, n), seed = seed + 1L)
  proxy <- crossfit_sensitive_proxy(
    X_proxy = X,
    s = s,
    r = r,
    K = 5L,
    ridge = 1e-2,
    seed = seed + 2L
  )

  z_ipw <- pseudo_sensitive(s, r, audit_rate, q = 0)
  z_aug <- pseudo_sensitive(s, r, audit_rate, q = proxy$q)

  data.frame(
    estimator = c("Truth", "Plug-in", "IPW", "Augmented"),
    moment = c(
      truth,
      estimate_score_moment(score, proxy$q),
      estimate_score_moment(score, z_ipw),
      estimate_score_moment(score, z_aug)
    )
  )
}


run_credit_main <- function(
    root = ".",
    config = credit_config(),
    save_results = TRUE,
    run_tag = NULL) {

  data <- load_credit_application_data(
    root = root,
    sensitive = config$sensitive,
    seed_split = config$seed_split
  )

  X_train <- data$X$train
  X_validation <- data$X$validation
  X_test <- data$X$test
  y_train <- data$y_split$train
  y_validation <- data$y_split$validation
  y_test <- data$y_split$test
  s_train <- data$s_split$train
  s_validation <- data$s_split$validation
  s_test <- data$s_split$test

  if (is.null(run_tag)) {
    run_tag <- paste(
      "credit", config$sensitive, config$path_type,
      if (config$quick) "quick" else "full",
      sep = "_"
    )
  }
  result_dir <- file.path(root, "results", run_tag)
  if (save_results) dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  unconstrained_fit <- fit_fair_logistic(
    X_train,
    y_train,
    cvec = rep(0, ncol(X_train)),
    lambda = 0,
    ridge = config$ridge_outcome
  )
  p_train_fixed <- clip_prob(plogis(drop(
    X_train %*% unconstrained_fit$coefficients
  )))
  p_test_unconstrained <- clip_prob(plogis(drop(
    X_test %*% unconstrained_fit$coefficients
  )))
  base_test_metrics <- fairness_metrics(
    y_test, p_test_unconstrained, s_test
  )
  true_train_moment <- estimate_score_moment(
    p_train_fixed, s_train
  )

  # -------------------------------------------------------------------------
  # Experiment 1: fixed-score fairness-moment estimation
  # -------------------------------------------------------------------------
  moment_rows <- vector(
    "list",
    length(config$audit_rates) * config$n_rep_moment * 7L
  )
  moment_id <- 0L

  for (rate_index in seq_along(config$audit_rates)) {
    audit_rate <- config$audit_rates[rate_index]

    for (replicate_index in seq_len(config$n_rep_moment)) {
      seed <- 100000L + 1000L * rate_index + replicate_index
      r_train <- draw_audit(
        rep(audit_rate, length(s_train)),
        seed = seed
      )

      rich <- crossfit_sensitive_proxy(
        data$X_proxy_rich$train,
        s_train,
        r_train,
        K = config$k_proxy,
        ridge = config$ridge_proxy,
        seed = seed
      )
      misspecified <- crossfit_sensitive_proxy(
        data$X_proxy_misspecified$train,
        s_train,
        r_train,
        K = config$k_proxy,
        ridge = config$ridge_proxy,
        seed = seed + 17L
      )

      estimators <- list(
        True = s_train,
        SoftProxy = rich$q,
        HardProxy = as.integer(rich$q >= 0.5),
        IPW = pseudo_sensitive(s_train, r_train, audit_rate, q = 0),
        Augmented = pseudo_sensitive(
          s_train, r_train, audit_rate, q = rich$q
        ),
        SoftProxyMisspecified = misspecified$q,
        AugmentedMisspecified = pseudo_sensitive(
          s_train, r_train, audit_rate, q = misspecified$q
        )
      )

      for (method in names(estimators)) {
        estimate <- estimate_score_moment(
          p_train_fixed, estimators[[method]]
        )
        moment_id <- moment_id + 1L
        moment_rows[[moment_id]] <- data.frame(
          audit_rate = audit_rate,
          replicate = replicate_index,
          method = method,
          audited_n = sum(r_train),
          estimate = estimate,
          truth = true_train_moment,
          error = estimate - true_train_moment,
          squared_error = (estimate - true_train_moment)^2
        )
      }
    }
  }

  moment <- do.call(rbind, moment_rows[seq_len(moment_id)])
  moment_summary <- summarize_replicates(
    moment,
    group_vars = c("audit_rate", "method"),
    value_vars = c("estimate", "error", "squared_error")
  )
  moment_summary$rmse <- sqrt(moment_summary$squared_error_mean)

  # -------------------------------------------------------------------------
  # Experiment 2: fair-GLM learning
  # -------------------------------------------------------------------------
  oracle_path <- fit_fair_path_by_type(
    config$path_type,
    X_train,
    y_train,
    z = s_train,
    lambdas = config$lambda_grid,
    ridge = config$ridge_outcome,
    beta_start = unconstrained_fit$coefficients
  )
  oracle_choice <- select_lambda_by_relative_moment_reduction(
    oracle_path,
    X_train,
    y_train,
    z_select = s_train,
    reduction_fraction = config$reduction_fraction
  )

  learning_rows <- list()
  frontier_rows <- list()
  learning_id <- 0L
  frontier_id <- 0L

  for (rate_index in seq_along(config$audit_rates)) {
    audit_rate <- config$audit_rates[rate_index]

    for (replicate_index in seq_len(config$n_rep_learning)) {
      seed <- 200000L + 1000L * rate_index + replicate_index
      r_train <- draw_audit(
        rep(audit_rate, length(s_train)),
        seed = seed
      )
      r_validation <- draw_audit(
        rep(audit_rate, length(s_validation)),
        seed = seed + 1L
      )

      rich <- crossfit_sensitive_proxy(
        data$X_proxy_rich$train,
        s_train,
        r_train,
        K = config$k_proxy,
        ridge = config$ridge_proxy,
        seed = seed,
        X_new = data$X_proxy_rich$validation
      )
      misspecified <- crossfit_sensitive_proxy(
        data$X_proxy_misspecified$train,
        s_train,
        r_train,
        K = config$k_proxy,
        ridge = config$ridge_proxy,
        seed = seed + 17L,
        X_new = data$X_proxy_misspecified$validation
      )

      methods <- list(
        Unconstrained = list(
          train = rich$q,
          validation = rich$q_new,
          fixed_path = list(
            fits = list(unconstrained_fit),
            lambdas = 0,
            path_type = "unconstrained"
          )
        ),
        Oracle = list(
          train = s_train,
          validation = s_validation,
          fixed_path = oracle_path,
          fixed_choice = oracle_choice
        ),
        SoftProxy = list(
          train = rich$q,
          validation = rich$q_new
        ),
        IPW = list(
          train = pseudo_sensitive(
            s_train, r_train, audit_rate, q = 0
          ),
          validation = pseudo_sensitive(
            s_validation, r_validation, audit_rate, q = 0
          )
        ),
        Augmented = list(
          train = pseudo_sensitive(
            s_train, r_train, audit_rate, q = rich$q
          ),
          validation = pseudo_sensitive(
            s_validation, r_validation, audit_rate, q = rich$q_new
          )
        ),
        AugmentedMisspecified = list(
          train = pseudo_sensitive(
            s_train, r_train, audit_rate, q = misspecified$q
          ),
          validation = pseudo_sensitive(
            s_validation, r_validation, audit_rate,
            q = misspecified$q_new
          )
        )
      )

      for (method in names(methods)) {
        information <- methods[[method]]

        if (!is.null(information$fixed_path)) {
          path <- information$fixed_path
        } else {
          path <- fit_fair_path_by_type(
            config$path_type,
            X_train,
            y_train,
            z = information$train,
            lambdas = config$lambda_grid,
            ridge = config$ridge_outcome,
            beta_start = unconstrained_fit$coefficients
          )
        }

        if (!is.null(information$fixed_choice)) {
          choice <- information$fixed_choice
        } else if (method == "Unconstrained") {
          diagnostic <- path_diagnostics(
            path, X_train, y_train, information$train
          )
          choice <- list(
            index = 1L,
            lambda = 0,
            table = diagnostic,
            baseline = diagnostic$estimated_moment[1L],
            target = diagnostic$estimated_moment[1L],
            feasible = TRUE
          )
        } else {
          choice <- select_lambda_by_relative_moment_reduction(
            path,
            X_train,
            y_train,
            z_select = information$train,
            reduction_fraction = config$reduction_fraction
          )
        }

        validation_diagnostic <- path_diagnostics(
          path, X_validation, y_validation, information$validation
        )
        test_predictions <- predict_fair_path(path, X_test)

        for (path_index in seq_along(path$lambdas)) {
          metrics <- fairness_metrics(
            y_test,
            test_predictions[, path_index],
            s_test
          )
          frontier_id <- frontier_id + 1L
          frontier_rows[[frontier_id]] <- data.frame(
            audit_rate = audit_rate,
            replicate = replicate_index,
            method = method,
            path_type = config$path_type,
            lambda = path$lambdas[path_index],
            selected = path_index == choice$index,
            audited_train_n = sum(r_train),
            audited_validation_n = sum(r_validation),
            training_estimated_moment =
              choice$table$estimated_moment[path_index],
            validation_estimated_moment =
              validation_diagnostic$estimated_moment[path_index],
            converged = choice$table$converged[path_index],
            gradient_norm = choice$table$gradient_norm[path_index],
            metrics
          )
        }

        selected_metrics <- fairness_metrics(
          y_test,
          test_predictions[, choice$index],
          s_test
        )
        learning_id <- learning_id + 1L
        learning_rows[[learning_id]] <- data.frame(
          audit_rate = audit_rate,
          replicate = replicate_index,
          method = method,
          path_type = config$path_type,
          lambda = choice$lambda,
          audited_train_n = sum(r_train),
          audited_validation_n = sum(r_validation),
          selection_estimated_moment =
            choice$table$estimated_moment[choice$index],
          selection_baseline_moment = choice$baseline,
          selection_target_moment = choice$target,
          selection_feasible = choice$feasible,
          validation_estimated_moment =
            validation_diagnostic$estimated_moment[choice$index],
          selected_converged =
            choice$table$converged[choice$index],
          selected_gradient_norm =
            choice$table$gradient_norm[choice$index],
          selected_metrics
        )
      }
    }
  }

  learning <- do.call(rbind, learning_rows)
  frontier <- do.call(rbind, frontier_rows)
  learning_summary <- summarize_replicates(
    learning,
    group_vars = c("audit_rate", "method", "path_type"),
    value_vars = c(
      "log_loss", "auc", "accuracy",
      "dp_score_signed", "eo_score_signed", "dp_score", "eo_score",
      "dp_decision_signed", "eo_decision_signed",
      "dp_decision", "eo_decision",
      "lambda", "selection_estimated_moment",
      "selection_baseline_moment", "selection_target_moment",
      "validation_estimated_moment"
    )
  )

  output <- list(
    config = config,
    data = data,
    unconstrained_fit = unconstrained_fit,
    base_test_metrics = base_test_metrics,
    true_train_moment = true_train_moment,
    moment = moment,
    moment_summary = moment_summary,
    learning = learning,
    learning_summary = learning_summary,
    frontier = frontier
  )

  if (save_results) {
    utils::write.csv(
      moment,
      file.path(result_dir, "credit_moment_replicates.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      moment_summary,
      file.path(result_dir, "credit_moment_summary.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      learning,
      file.path(result_dir, "credit_learning_selected.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      learning_summary,
      file.path(result_dir, "credit_learning_summary.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      frontier,
      file.path(result_dir, "credit_learning_frontiers.csv"),
      row.names = FALSE
    )
    saveRDS(
      list(
        config = config,
        sensitive_description = data$sensitive_description,
        n = nrow(data$credit),
        n_train = nrow(X_train),
        n_validation = nrow(X_validation),
        n_test = nrow(X_test),
        sensitive_prevalence = mean(data$s),
        favorable_outcome_prevalence = mean(data$y),
        core_version = PROXY_FAIR_GLM_CORE_VERSION,
        R_version = R.version.string
      ),
      file.path(result_dir, "credit_metadata.rds")
    )
  }

  output
}


plot_credit_main <- function(result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to make figures.")
  }

  moment <- result$moment
  learning <- result$learning
  frontier <- result$frontier

  moment$audit_percent <- 100 * moment$audit_rate
  learning$audit_percent <- 100 * learning$audit_rate
  frontier$audit_percent <- 100 * frontier$audit_rate

  rmse <- stats::aggregate(
    squared_error ~ audit_percent + method,
    data = subset(moment, method != "True"),
    FUN = mean
  )
  rmse$rmse <- sqrt(rmse$squared_error)
  rmse <- subset(
    rmse,
    method %in% c(
      "Augmented", "AugmentedMisspecified",
      "IPW", "SoftProxy", "SoftProxyMisspecified"
    )
  )

  p_rmse <- ggplot2::ggplot(
    rmse,
    ggplot2::aes(
      x = audit_percent, y = rmse,
      group = method, linetype = method, shape = method
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Audit rate (%)",
      y = "RMSE of the fairness moment"
    ) +
    ggplot2::theme_minimal()

  selected <- stats::aggregate(
    cbind(log_loss, dp_score, eo_score) ~ audit_percent + method,
    data = learning,
    FUN = mean
  )
  selected <- subset(
    selected,
    method %in% c(
      "Augmented", "AugmentedMisspecified",
      "IPW", "Oracle", "SoftProxy"
    )
  )

  p_dp <- ggplot2::ggplot(
    selected,
    ggplot2::aes(
      x = audit_percent, y = dp_score,
      group = method, linetype = method, shape = method
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Audit rate (%)",
      y = "Test score demographic-parity gap"
    ) +
    ggplot2::theme_minimal()

  representative_rate <- if (5 %in% frontier$audit_percent) {
    5
  } else {
    sort(unique(frontier$audit_percent))[1L]
  }
  frontier_5 <- subset(
    frontier,
    audit_percent == representative_rate &
      method %in% c(
        "Augmented", "AugmentedMisspecified",
        "IPW", "Oracle", "SoftProxy"
      )
  )
  frontier_summary <- stats::aggregate(
    cbind(log_loss, dp_score) ~ method + lambda,
    data = frontier_5,
    FUN = mean
  )
  frontier_summary <- frontier_summary[
    order(frontier_summary$method, frontier_summary$lambda),
  ]

  p_frontier <- ggplot2::ggplot(
    frontier_summary,
    ggplot2::aes(
      x = dp_score, y = log_loss,
      group = method, linetype = method, shape = method
    )
  ) +
    ggplot2::geom_path() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Test score demographic-parity gap",
      y = "Test log loss"
    ) +
    ggplot2::theme_minimal()

  list(
    rmse = p_rmse,
    selected_dp = p_dp,
    frontier = p_frontier,
    selected_table = selected
  )
}


run_proxy_stress <- function(
    root = ".",
    quick = TRUE,
    sensitive = "sex",
    audit_rates = c(0.01, 0.02, 0.05, 0.10, 0.20),
    n_repetitions = if (quick) 25L else 300L) {

  data <- load_credit_application_data(root, sensitive = sensitive)
  X_train <- data$X$train
  y_train <- data$y_split$train
  s_train <- data$s_split$train
  n_train <- length(s_train)

  outcome_fit <- fit_fair_logistic(
    X_train, y_train,
    cvec = rep(0, ncol(X_train)),
    lambda = 0, ridge = 1e-4
  )
  score <- clip_prob(plogis(drop(
    X_train %*% outcome_fit$coefficients
  )))
  truth <- estimate_score_moment(score, s_train)

  X_intercept <- matrix(1, nrow = n_train, ncol = 1L)
  colnames(X_intercept) <- "(Intercept)"

  rows <- list()
  id <- 0L

  for (rate_index in seq_along(audit_rates)) {
    audit_rate <- audit_rates[rate_index]
    for (replicate_index in seq_len(n_repetitions)) {
      seed <- 310000L + 1000L * rate_index + replicate_index
      r <- draw_audit(rep(audit_rate, n_train), seed)

      rich <- crossfit_sensitive_proxy(
        data$X_proxy_rich$train, s_train, r,
        K = 5L, ridge = 1e-2, seed = seed
      )$q
      reduced <- crossfit_sensitive_proxy(
        data$X_proxy_misspecified$train, s_train, r,
        K = 5L, ridge = 1e-2, seed = seed + 17L
      )$q
      intercept_only <- crossfit_sensitive_proxy(
        X_intercept, s_train, r,
        K = 5L, ridge = 0, seed = seed + 29L,
        fallback_probability = 0.5
      )$q

      cases <- list(
        IPW = rep(0, n_train),
        InterceptOnly = intercept_only,
        ConstantHalf = rep(0.5, n_train),
        FittedRich = rich,
        FittedMisspecified = reduced,
        InvertedFitted = 1 - rich,
        ConstantOne = rep(1, n_train),
        OracleAdversarial = 1 - s_train
      )

      for (proxy_case in names(cases)) {
        q <- cases[[proxy_case]]
        z <- pseudo_sensitive(s_train, r, audit_rate, q)
        estimate <- estimate_score_moment(score, z)
        id <- id + 1L
        rows[[id]] <- data.frame(
          audit_rate = audit_rate,
          replicate = replicate_index,
          proxy_case = proxy_case,
          estimate = estimate,
          error = estimate - truth,
          squared_error = (estimate - truth)^2
        )
      }
    }
  }

  results <- do.call(rbind, rows)
  variance <- stats::aggregate(
    estimate ~ audit_rate + proxy_case,
    data = results,
    FUN = stats::var
  )
  names(variance)[3L] <- "empirical_variance"
  ipw <- subset(
    variance, proxy_case == "IPW",
    select = c("audit_rate", "empirical_variance")
  )
  names(ipw)[2L] <- "ipw_variance"
  variance <- merge(variance, ipw, by = "audit_rate")
  variance$variance_ratio_to_ipw <-
    variance$empirical_variance / variance$ipw_variance

  list(results = results, variance = variance)
}


run_variance_validation <- function(
    root = ".",
    quick = TRUE,
    sensitive = "sex",
    pilot_rate = 0.20,
    audit_rates = c(0.01, 0.02, 0.05, 0.10, 0.20),
    reps_1pct = if (quick) 200L else 5000L,
    reps_other = if (quick) 100L else 1500L) {

  data <- load_credit_application_data(root, sensitive = sensitive)
  X_train <- data$X$train
  y_train <- data$y_split$train
  s_train <- data$s_split$train

  outcome_fit <- fit_fair_logistic(
    X_train, y_train,
    cvec = rep(0, ncol(X_train)),
    lambda = 0, ridge = 1e-4
  )
  score <- clip_prob(plogis(drop(
    X_train %*% outcome_fit$coefficients
  )))
  truth <- estimate_score_moment(score, s_train)

  pilot_r <- draw_audit(
    rep(pilot_rate, length(s_train)),
    seed = 510001L
  )
  pilot_proxy <- fit_sensitive_proxy_from_audit(
    data$X_proxy_rich$train,
    s_train,
    pilot_r,
    ridge = 1e-2,
    X_new = data$X_proxy_rich$train
  )$q

  rows <- list()
  id <- 0L

  for (rate_index in seq_along(audit_rates)) {
    audit_rate <- audit_rates[rate_index]
    n_repetitions <- if (rate_index == 1L) reps_1pct else reps_other

    exact_variances <- c(
      Augmented = exact_design_variance(
        score, s_train, pilot_proxy, audit_rate
      ),
      IPW = exact_design_variance(
        score, s_train, rep(0, length(s_train)), audit_rate
      )
    )

    for (replicate_index in seq_len(n_repetitions)) {
      r <- draw_audit(
        rep(audit_rate, length(s_train)),
        510000L + 100000L * rate_index + replicate_index
      )

      q_methods <- list(
        Augmented = pilot_proxy,
        IPW = rep(0, length(s_train))
      )

      for (method in names(q_methods)) {
        q <- q_methods[[method]]
        z <- pseudo_sensitive(s_train, r, audit_rate, q)
        estimate <- estimate_score_moment(score, z)
        vhat <- estimate_design_variance(
          score, s_train, r, q, audit_rate
        )
        id <- id + 1L
        rows[[id]] <- data.frame(
          audit_rate = audit_rate,
          replicate = replicate_index,
          method = method,
          estimate = estimate,
          error = estimate - truth,
          estimated_variance = vhat,
          exact_design_variance = exact_variances[[method]]
        )
      }
    }
  }

  results <- do.call(rbind, rows)
  pieces <- split(
    results,
    interaction(results$audit_rate, results$method, drop = TRUE)
  )
  summary <- do.call(rbind, lapply(pieces, function(d) {
    empirical_variance <- stats::var(d$estimate)
    data.frame(
      audit_rate = d$audit_rate[1L],
      method = d$method[1L],
      empirical_variance = empirical_variance,
      empirical_variance_se =
        sample_variance_standard_error(d$estimate),
      empirical_variance_ci_lower = max(
        empirical_variance -
          1.96 * sample_variance_standard_error(d$estimate),
        .Machine$double.eps
      ),
      empirical_variance_ci_upper =
        empirical_variance +
        1.96 * sample_variance_standard_error(d$estimate),
      mean_estimated_variance = mean(d$estimated_variance),
      exact_design_variance = unique(d$exact_design_variance),
      exact_to_empirical_ratio =
        unique(d$exact_design_variance) / empirical_variance,
      bias = mean(d$error),
      rmse = sqrt(mean(d$error^2)),
      n_rep = nrow(d)
    )
  }))
  rownames(summary) <- NULL

  list(results = results, summary = summary)
}


run_exact_vs_convex <- function(
    root = ".",
    quick = TRUE,
    sensitive = "sex",
    audit_rate = 0.05,
    n_repetitions = if (quick) 3L else 50L,
    reduction_fraction = 0.50) {

  data <- load_credit_application_data(root, sensitive = sensitive)
  X_train <- data$X$train
  X_test <- data$X$test
  y_train <- data$y_split$train
  y_test <- data$y_split$test
  s_train <- data$s_split$train
  s_test <- data$s_split$test

  lambda_grid <- c(
    0, 10^seq(-3, 4, length.out = if (quick) 15L else 31L)
  )
  unconstrained <- fit_fair_logistic(
    X_train, y_train,
    cvec = rep(0, ncol(X_train)),
    lambda = 0, ridge = 1e-4
  )

  oracle_paths <- lapply(c("exact", "convex"), function(path_type) {
    fit_fair_path_by_type(
      path_type,
      X_train, y_train, s_train,
      lambda_grid, 1e-4,
      unconstrained$coefficients
    )
  })
  names(oracle_paths) <- c("exact", "convex")

  selected_rows <- list()
  frontier_rows <- list()
  selected_id <- 0L
  frontier_id <- 0L

  for (replicate_index in seq_len(n_repetitions)) {
    seed <- 410000L + replicate_index
    r <- draw_audit(rep(audit_rate, length(s_train)), seed)
    proxy <- crossfit_sensitive_proxy(
      data$X_proxy_rich$train,
      s_train, r,
      K = 5L, ridge = 1e-2, seed = seed
    )$q

    methods <- list(
      Oracle = s_train,
      SoftProxy = proxy,
      IPW = pseudo_sensitive(s_train, r, audit_rate, q = 0),
      Augmented = pseudo_sensitive(
        s_train, r, audit_rate, q = proxy
      )
    )

    for (method in names(methods)) {
      z <- methods[[method]]

      for (path_type in c("exact", "convex")) {
        path <- if (method == "Oracle") {
          oracle_paths[[path_type]]
        } else {
          fit_fair_path_by_type(
            path_type, X_train, y_train, z,
            lambda_grid, 1e-4,
            unconstrained$coefficients
          )
        }

        diagnostics <- path_diagnostics(path, X_train, y_train, z)
        predictions <- predict_fair_path(path, X_test)

        zero_index <- which.min(abs(diagnostics$lambda))
        exact_baseline <- diagnostics$estimated_moment[zero_index]
        surrogate_baseline <- diagnostics$surrogate_moment[zero_index]

        diagnostics$relative_exact_moment <- if (
          is.finite(exact_baseline) && exact_baseline > 1e-12
        ) {
          diagnostics$estimated_moment / exact_baseline
        } else {
          0
        }
        diagnostics$relative_surrogate_moment <- if (
          is.finite(surrogate_baseline) && surrogate_baseline > 1e-12
        ) {
          diagnostics$surrogate_moment / surrogate_baseline
        } else {
          0
        }

        specs <- if (path_type == "convex") {
          c(ExactScoreMoment = "exact", SurrogateMoment = "surrogate")
        } else {
          c(ExactScoreMoment = "exact")
        }

        choices <- lapply(specs, function(metric) {
          select_lambda_by_relative_diagnostic_reduction(
            path, X_train, y_train, z,
            reduction_fraction, metric
          )
        })

        for (j in seq_along(path$lambdas)) {
          metrics <- fairness_metrics(
            y_test, predictions[, j], s_test
          )
          frontier_id <- frontier_id + 1L
          frontier_rows[[frontier_id]] <- data.frame(
            replicate = replicate_index,
            audit_rate = audit_rate,
            method = method,
            path_type = path_type,
            lambda = path$lambdas[j],
            training_exact_moment = diagnostics$estimated_moment[j],
            training_surrogate_moment = diagnostics$surrogate_moment[j],
            relative_training_exact_moment =
              diagnostics$relative_exact_moment[j],
            relative_training_surrogate_moment =
              diagnostics$relative_surrogate_moment[j],
            logit_spread = diagnostics$logit_spread[j],
            converged = diagnostics$converged[j],
            gradient_norm = diagnostics$gradient_norm[j],
            metrics
          )
        }

        for (target in names(choices)) {
          choice <- choices[[target]]
          metrics <- fairness_metrics(
            y_test, predictions[, choice$index], s_test
          )
          selected_id <- selected_id + 1L
          selected_rows[[selected_id]] <- data.frame(
            replicate = replicate_index,
            audit_rate = audit_rate,
            method = method,
            path_type = path_type,
            selection_target = target,
            selection_metric = choice$selection_metric,
            lambda = choice$lambda,
            feasible = choice$feasible,
            training_exact_moment =
              diagnostics$estimated_moment[choice$index],
            training_surrogate_moment =
              diagnostics$surrogate_moment[choice$index],
            relative_training_exact_moment =
              diagnostics$relative_exact_moment[choice$index],
            relative_training_surrogate_moment =
              diagnostics$relative_surrogate_moment[choice$index],
            metrics
          )
        }
      }
    }
  }

  selected <- do.call(rbind, selected_rows)
  frontiers <- do.call(rbind, frontier_rows)

  selected$feasible_numeric <- as.numeric(selected$feasible)
  summary <- summarize_replicates(
    selected,
    group_vars = c(
      "method", "path_type",
      "selection_target", "selection_metric"
    ),
    value_vars = c(
      "log_loss", "dp_score", "eo_score",
      "lambda", "training_exact_moment",
      "training_surrogate_moment",
      "relative_training_exact_moment",
      "relative_training_surrogate_moment",
      "feasible_numeric"
    )
  )

  list(
    selected = selected,
    frontiers = frontiers,
    summary = summary,
    audit_rate = audit_rate,
    lambda_grid = lambda_grid
  )
}

run_neyman_total_budget <- function(
    root = ".",
    quick = TRUE,
    sensitive = "sex",
    total_audit_rates = if (quick) c(0.05, 0.10) else c(0.05, 0.10, 0.20),
    pilot_fractions = if (quick) 0.20 else c(0.10, 0.20, 0.40),
    n_repetitions = if (quick) 10L else 200L) {

  data <- load_credit_application_data(root, sensitive = sensitive)
  X_train <- data$X$train
  y_train <- data$y_split$train
  s_train <- data$s_split$train
  n_train <- length(s_train)

  outcome_fit <- fit_fair_logistic(
    X_train, y_train,
    cvec = rep(0, ncol(X_train)),
    lambda = 0, ridge = 1e-4
  )
  score <- clip_prob(plogis(drop(
    X_train %*% outcome_fit$coefficients
  )))
  centered_score <- score - mean(score)
  truth <- estimate_score_moment(score, s_train)

  rows <- list()
  id <- 0L

  for (total_index in seq_along(total_audit_rates)) {
    total_rate <- total_audit_rates[total_index]
    total_n <- max(2L, round(total_rate * n_train))

    for (fraction_index in seq_along(pilot_fractions)) {
      pilot_fraction <- pilot_fractions[fraction_index]
      pilot_n <- max(
        1L,
        min(total_n - 1L, round(pilot_fraction * total_n))
      )
      main_budget_n <- total_n - pilot_n

      for (replicate_index in seq_len(n_repetitions)) {
        seed_base <- 620000L +
          100000L * total_index +
          10000L * fraction_index +
          replicate_index

        pilot_r <- draw_fixed_size_audit(
          n_train, pilot_n, seed = seed_base
        )
        nonpilot <- which(pilot_r == 0L)
        main_rate_nonpilot <- main_budget_n / length(nonpilot)

        pilot_fit <- fit_sensitive_proxy_from_audit(
          data$X_proxy_rich$train,
          s_train,
          pilot_r,
          ridge = 1e-2,
          X_new = data$X_proxy_rich$train
        )
        q <- pilot_fit$q

        feasible_score <- q[nonpilot] * (1 - q[nonpilot]) *
          centered_score[nonpilot]^2
        oracle_score <- (s_train[nonpilot] - q[nonpilot])^2 *
          centered_score[nonpilot]^2

        designs <- list(
          TwoPhaseUniform = rep(
            main_rate_nonpilot, length(nonpilot)
          ),
          TwoPhasePilotNeyman = allocate_audit_probabilities(
            feasible_score,
            budget_rate = main_rate_nonpilot,
            pi_min = main_rate_nonpilot / 20
          ),
          TwoPhaseOracleNeyman = allocate_audit_probabilities(
            oracle_score,
            budget_rate = main_rate_nonpilot,
            pi_min = main_rate_nonpilot / 20
          )
        )

        for (design_index in seq_along(designs)) {
          design <- names(designs)[design_index]
          pi_nonpilot <- designs[[design]]
          pi_full <- rep(1, n_train)
          pi_full[nonpilot] <- pi_nonpilot
          main_r <- integer(n_train)
          main_r[nonpilot] <- draw_audit(
            pi_nonpilot,
            seed_base + 1000L * design_index
          )

          z <- two_phase_pseudo_sensitive(
            s_train, pilot_r, main_r, pi_full, q
          )
          estimate <- estimate_score_moment(score, z)

          id <- id + 1L
          rows[[id]] <- data.frame(
            total_audit_rate = total_rate,
            pilot_fraction = pilot_fraction,
            replicate = replicate_index,
            design = design,
            estimate = estimate,
            error = estimate - truth,
            squared_error = (estimate - truth)^2
          )
        }
      }
    }
  }

  results <- do.call(rbind, rows)
  pieces <- split(
    results,
    interaction(
      results$total_audit_rate,
      results$pilot_fraction,
      results$design,
      drop = TRUE
    )
  )
  summary <- do.call(rbind, lapply(pieces, function(d) {
    data.frame(
      total_audit_rate = d$total_audit_rate[1L],
      pilot_fraction = d$pilot_fraction[1L],
      design = d$design[1L],
      empirical_variance = stats::var(d$estimate),
      rmse = sqrt(mean(d$squared_error)),
      bias = mean(d$error),
      n_rep = nrow(d)
    )
  }))
  rownames(summary) <- NULL

  uniform <- subset(
    summary,
    design == "TwoPhaseUniform",
    select = c(
      "total_audit_rate", "pilot_fraction",
      "empirical_variance", "rmse"
    )
  )
  names(uniform)[3:4] <- c("uniform_variance", "uniform_rmse")
  summary <- merge(
    summary, uniform,
    by = c("total_audit_rate", "pilot_fraction")
  )
  summary$variance_ratio_to_uniform <-
    summary$empirical_variance / summary$uniform_variance
  summary$rmse_ratio_to_uniform <-
    summary$rmse / summary$uniform_rmse

  list(results = results, summary = summary)
}


run_selection_robustness <- function(
    root = ".",
    quick = TRUE,
    sensitive = "sex",
    training_audit_rate = 0.05,
    validation_audit_rates = c(0.05, 0.10, 0.20),
    n_repetitions = if (quick) 3L else 40L,
    reduction_fraction = 0.50,
    ucb_z_value = 1.645) {

  data <- load_credit_application_data(root, sensitive = sensitive)
  X_train <- data$X$train
  X_validation <- data$X$validation
  X_test <- data$X$test
  y_train <- data$y_split$train
  y_validation <- data$y_split$validation
  y_test <- data$y_split$test
  s_train <- data$s_split$train
  s_validation <- data$s_split$validation
  s_test <- data$s_split$test

  lambda_grid <- c(
    0, 10^seq(-3, 4, length.out = if (quick) 15L else 31L)
  )
  unconstrained <- fit_fair_logistic(
    X_train, y_train,
    cvec = rep(0, ncol(X_train)),
    lambda = 0, ridge = 1e-4
  )

  rows <- list()
  id <- 0L

  append_row <- function(
      replicate_index, method, rule, validation_rate,
      choice, predictions_test) {
    metrics <- fairness_metrics(
      y_test, predictions_test[, choice$index], s_test
    )
    data.frame(
      replicate = replicate_index,
      method = method,
      selection_rule = rule,
      validation_audit_rate = validation_rate,
      lambda = choice$lambda,
      feasible = choice$feasible,
      metrics
    )
  }

  for (replicate_index in seq_len(n_repetitions)) {
    seed <- 710000L + replicate_index
    r_train <- draw_audit(
      rep(training_audit_rate, length(s_train)), seed
    )

    proxy <- crossfit_sensitive_proxy(
      data$X_proxy_rich$train,
      s_train, r_train,
      K = 5L, ridge = 1e-2, seed = seed,
      X_new = data$X_proxy_rich$validation
    )

    methods <- list(
      SoftProxy = proxy$q,
      IPW = pseudo_sensitive(
        s_train, r_train, training_audit_rate, q = 0
      ),
      Augmented = pseudo_sensitive(
        s_train, r_train, training_audit_rate, q = proxy$q
      )
    )

    for (method in names(methods)) {
      z_train <- methods[[method]]
      path <- fit_exact_fair_path(
        X_train, y_train, z_train,
        lambda_grid, 1e-4,
        unconstrained$coefficients
      )
      predictions_test <- predict_fair_path(path, X_test)

      choice_training <- select_lambda_by_relative_moment_reduction(
        path, X_train, y_train, z_train, reduction_fraction
      )
      id <- id + 1L
      rows[[id]] <- append_row(
        replicate_index, method, "TrainingAudit",
        NA_real_, choice_training, predictions_test
      )

      choice_full <- select_lambda_by_relative_moment_reduction(
        path, X_validation, y_validation,
        s_validation, reduction_fraction
      )
      id <- id + 1L
      rows[[id]] <- append_row(
        replicate_index, method, "FullDemographicsValidation",
        1, choice_full, predictions_test
      )

      if (method %in% c("IPW", "Augmented")) {
        for (rate_index in seq_along(validation_audit_rates)) {
          validation_rate <- validation_audit_rates[rate_index]
          r_validation <- draw_audit(
            rep(validation_rate, length(s_validation)),
            seed + 1000L * rate_index
          )
          q_validation <- if (method == "IPW") {
            rep(0, length(s_validation))
          } else {
            proxy$q_new
          }
          z_validation <- pseudo_sensitive(
            s_validation, r_validation,
            validation_rate, q_validation
          )

          choice_point <- select_lambda_by_relative_moment_reduction(
            path, X_validation, y_validation,
            z_validation, reduction_fraction
          )
          id <- id + 1L
          rows[[id]] <- append_row(
            replicate_index, method, "ValidationAuditPoint",
            validation_rate, choice_point, predictions_test
          )

          choice_ucb <- select_lambda_by_relative_moment_ucb(
            path,
            X_validation,
            y_validation,
            z_validation,
            s_validation,
            r_validation,
            validation_rate,
            q_validation,
            reduction_fraction,
            ucb_z_value
          )
          id <- id + 1L
          rows[[id]] <- append_row(
            replicate_index, method, "ValidationAuditUCB",
            validation_rate, choice_ucb, predictions_test
          )
        }
      }
    }
  }

  results <- do.call(rbind, rows)
  results$feasible_numeric <- as.numeric(results$feasible)
  results$validation_rate_key <- ifelse(
    is.na(results$validation_audit_rate),
    -1,
    results$validation_audit_rate
  )
  summary <- summarize_replicates(
    results,
    group_vars = c(
      "method", "selection_rule", "validation_rate_key"
    ),
    value_vars = c(
      "lambda", "log_loss", "dp_score",
      "eo_score", "feasible_numeric"
    )
  )
  names(summary)[names(summary) == "validation_rate_key"] <-
    "validation_audit_rate"
  summary$validation_audit_rate[
    summary$validation_audit_rate < 0
  ] <- NA_real_

  list(results = results, summary = summary)
}


# -----------------------------------------------------------------------------
# Additional v5.1 experiment recovered from the original scripts
# -----------------------------------------------------------------------------

run_neyman_external_pilot <- function(
    root = ".",
    quick = TRUE,
    sensitive = "sex",
    pilot_rates = if (quick) c(0.05, 0.20) else c(0.01, 0.02, 0.05, 0.20),
    audit_rates = c(0.01, 0.02, 0.05, 0.10, 0.20),
    n_repetitions = if (quick) 50L else 1000L,
    save_results = FALSE,
    run_tag = NULL) {

  data <- load_credit_application_data(root, sensitive = sensitive)
  X_train <- data$X$train
  y_train <- data$y_split$train
  s_train <- data$s_split$train
  n_train <- length(s_train)

  outcome_fit <- fit_fair_logistic(
    X_train,
    y_train,
    cvec = rep(0, ncol(X_train)),
    lambda = 0,
    ridge = 1e-4
  )
  score <- clip_prob(plogis(drop(
    X_train %*% outcome_fit$coefficients
  )))
  centered_score <- score - mean(score)
  truth <- estimate_score_moment(score, s_train)

  rows <- list()
  allocation_rows <- list()
  pilot_rows <- list()
  row_id <- 0L
  allocation_id <- 0L
  pilot_id <- 0L

  for (pilot_index in seq_along(pilot_rates)) {
    pilot_rate <- pilot_rates[pilot_index]
    pilot_n <- max(1L, round(pilot_rate * n_train))
    pilot_r <- draw_fixed_size_audit(
      n_train,
      pilot_n,
      seed = 610001L + 1000L * pilot_index
    )
    pilot_fit <- fit_sensitive_proxy_from_audit(
      data$X_proxy_rich$train,
      s_train,
      pilot_r,
      ridge = 1e-2,
      X_new = data$X_proxy_rich$train
    )
    pilot_proxy <- pilot_fit$q
    diagnostics <- proxy_diagnostics(s_train, pilot_proxy)

    pilot_id <- pilot_id + 1L
    pilot_rows[[pilot_id]] <- data.frame(
      pilot_rate = pilot_rate,
      pilot_n = sum(pilot_r),
      proxy_auc = diagnostics$auc,
      proxy_brier = diagnostics$brier,
      proxy_ece = diagnostics$ece
    )

    feasible_variance_score <-
      pilot_proxy * (1 - pilot_proxy) * centered_score^2
    oracle_variance_score <-
      (s_train - pilot_proxy)^2 * centered_score^2

    for (rate_index in seq_along(audit_rates)) {
      audit_rate <- audit_rates[rate_index]
      designs <- list(
        Uniform = rep(audit_rate, n_train),
        PilotNeyman = allocate_audit_probabilities(
          feasible_variance_score,
          budget_rate = audit_rate,
          pi_min = audit_rate / 20
        ),
        OracleNeyman = allocate_audit_probabilities(
          oracle_variance_score,
          budget_rate = audit_rate,
          pi_min = audit_rate / 20
        )
      )

      for (design in names(designs)) {
        pi <- designs[[design]]
        allocation_id <- allocation_id + 1L
        allocation_rows[[allocation_id]] <- data.frame(
          pilot_rate = pilot_rate,
          pilot_n = sum(pilot_r),
          audit_rate = audit_rate,
          design = design,
          mean_pi = mean(pi),
          min_pi = min(pi),
          q10_pi = unname(stats::quantile(pi, 0.10, type = 8)),
          median_pi = stats::median(pi),
          q90_pi = unname(stats::quantile(pi, 0.90, type = 8)),
          max_pi = max(pi),
          exact_design_variance = exact_design_variance(
            score, s_train, pilot_proxy, pi
          )
        )
      }

      for (replicate_index in seq_len(n_repetitions)) {
        for (design_index in seq_along(designs)) {
          design <- names(designs)[design_index]
          pi <- designs[[design]]
          seed <- 610000L + 1000000L * pilot_index +
            10000L * rate_index +
            1000L * design_index +
            replicate_index
          r <- draw_audit(pi, seed)
          z <- pseudo_sensitive(s_train, r, pi, pilot_proxy)
          estimate <- estimate_score_moment(score, z)

          row_id <- row_id + 1L
          rows[[row_id]] <- data.frame(
            pilot_rate = pilot_rate,
            pilot_n = sum(pilot_r),
            audit_rate = audit_rate,
            replicate = replicate_index,
            design = design,
            audited_n = sum(r),
            estimate = estimate,
            error = estimate - truth,
            squared_error = (estimate - truth)^2,
            estimated_variance = estimate_design_variance(
              score, s_train, r, pilot_proxy, pi
            )
          )
        }
      }
    }
  }

  results <- do.call(rbind, rows)
  allocations <- do.call(rbind, allocation_rows)
  pilots <- do.call(rbind, pilot_rows)

  pieces <- split(
    results,
    interaction(
      results$pilot_rate,
      results$audit_rate,
      results$design,
      drop = TRUE
    )
  )
  summary <- do.call(rbind, lapply(pieces, function(d) {
    data.frame(
      pilot_rate = d$pilot_rate[1L],
      pilot_n = d$pilot_n[1L],
      audit_rate = d$audit_rate[1L],
      design = d$design[1L],
      mean_audited_n = mean(d$audited_n),
      bias = mean(d$error),
      empirical_variance = stats::var(d$estimate),
      mean_estimated_variance = mean(d$estimated_variance),
      rmse = sqrt(mean(d$squared_error)),
      n_rep = nrow(d)
    )
  }))
  rownames(summary) <- NULL

  uniform <- subset(
    summary,
    design == "Uniform",
    select = c(
      "pilot_rate", "audit_rate",
      "empirical_variance", "rmse"
    )
  )
  names(uniform)[3:4] <- c("uniform_variance", "uniform_rmse")
  summary <- merge(
    summary,
    uniform,
    by = c("pilot_rate", "audit_rate")
  )
  summary$variance_ratio_to_uniform <-
    summary$empirical_variance / summary$uniform_variance
  summary$rmse_ratio_to_uniform <-
    summary$rmse / summary$uniform_rmse

  output <- list(
    results = results,
    summary = summary,
    allocations = allocations,
    pilots = pilots,
    truth = truth
  )

  if (save_results) {
    if (is.null(run_tag)) {
      run_tag <- paste0(
        "neyman_external_pilot_v51_",
        sensitive, "_",
        if (quick) "quick" else "full"
      )
    }
    result_dir <- file.path(root, "results", run_tag)
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
      results,
      file.path(result_dir, "neyman_external_replicates.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      summary,
      file.path(result_dir, "neyman_external_summary.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      allocations,
      file.path(result_dir, "neyman_external_allocations.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      pilots,
      file.path(result_dir, "neyman_external_pilots.csv"),
      row.names = FALSE
    )
  }

  output
}


# -----------------------------------------------------------------------------
# Manuscript figures
# -----------------------------------------------------------------------------

save_pdf_plot <- function(plot, filename, width = 6.6, height = 3.9) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to save figures.")
  }
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf
  )
  invisible(filename)
}


save_main_manuscript_figures <- function(
    result,
    figures_dir = "figures",
    sensitive = "sex") {

  plots <- plot_credit_main(result)

  files <- c(
    moment_rmse = file.path(
      figures_dir,
      paste0("credit_moment_rmse_", sensitive, ".pdf")
    ),
    selected_dp = file.path(
      figures_dir,
      paste0("credit_selected_dp_", sensitive, ".pdf")
    ),
    frontier = file.path(
      figures_dir,
      paste0("credit_frontier_", sensitive, ".pdf")
    )
  )

  save_pdf_plot(plots$rmse, files[["moment_rmse"]])
  save_pdf_plot(plots$selected_dp, files[["selected_dp"]])
  save_pdf_plot(
    plots$frontier,
    files[["frontier"]],
    width = 6.6,
    height = 4.15
  )

  files
}


save_proxy_stress_figure <- function(
    stress,
    filename = file.path(
      "figures", "proxy_stress_variance_ratio.pdf"
    )) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to save figures.")
  }

  d <- stress$variance
  d$audit_percent <- 100 * d$audit_rate

  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = proxy_case,
      y = variance_ratio_to_ipw,
      group = factor(audit_percent),
      linetype = factor(audit_percent),
      shape = factor(audit_percent)
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Auxiliary proxy",
      y = "Variance ratio: augmented / IPW",
      linetype = "Audit rate (%)",
      shape = "Audit rate (%)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 28, hjust = 1
      )
    )

  save_pdf_plot(p, filename, width = 7.5, height = 4.4)
  invisible(filename)
}


# =============================================================================
# Exact recovered source/target proxy-provenance simulation
# =============================================================================

run_proxy_provenance <- function(
    root = ".",
    quick = TRUE,
    seed = 918273L,
    out_dir = file.path(root, "simulation_outputs_sameX"),
    figures_dir = file.path(root, "figures")) {

  QUICK <- isTRUE(quick)
  SEED <- as.integer(seed)

  if (QUICK) {
    N_SOURCE <- 30000L
    N_TARGET <- 8000L
    N_REP <- 100L
  } else {
    N_SOURCE <- 120000L
    N_TARGET <- 30000L
    N_REP <- 500L
  }

  AUDIT_RATES <- c(0.01, 0.02, 0.05, 0.10, 0.20)

  K_FOLDS <- 5L
  RIDGE_EXTERNAL <- 1e-2
  RIDGE_TARGET <- 1e-2
  RIDGE_RECAL <- 1e-3
  SMOOTHING <- 0.5

  MIS_ALPHA <- 0.90
  MIS_GAMMA <- 1.35

  P_SOURCE <- 0.28
  P_TARGET_TRANSPORTED <- P_SOURCE
  P_TARGET_SHIFTED <- 0.40

  N_SURNAME <- 80L
  N_GEO <- 24L
  NAME_STRENGTH <- 1.20
  GEO_STRENGTH <- 0.95

  SOURCE_X1_S <- 0.85
  SOURCE_X2_S <- -0.55
  SOURCE_X3_X1 <- 0.35
  SOURCE_X3_S <- 0.20

  SHIFT_X1_S <- -0.20
  SHIFT_X2_S <- 0.15
  SHIFT_X3_X1 <- 0.20
  SHIFT_X3_S <- 0.10

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # =============================================================================
  # 2. Generic helpers
  # =============================================================================
  
  expit <- function(x) plogis(x)
  
  logit_clip <- function(p, eps = 1e-6) {
    qlogis(pmin(pmax(p, eps), 1 - eps))
  }
  
  normalize_prob <- function(x) {
    x <- pmax(as.numeric(x), 1e-14)
    x / sum(x)
  }
  
  sample_cat <- function(prob, n) {
    if (n <= 0L) return(integer(0))
    sample.int(length(prob), size = n, replace = TRUE, prob = prob)
  }
  
  make_binary_category_tables <- function(k, strength, seed) {
    set.seed(seed)
    base <- rgamma(k, shape = 2.5, rate = 1)
    z <- rnorm(k)
  
    p0 <- normalize_prob(base * exp(-0.5 * strength * z))
    p1 <- normalize_prob(base * exp( 0.5 * strength * z))
  
    list(p0 = p0, p1 = p1)
  }
  
  
  # =============================================================================
  # 3. Population generator
  # =============================================================================
  
  generate_population <- function(
      n,
      p_s,
      surname_tab,
      geo_tab,
      seed,
      x1_s = SOURCE_X1_S,
      x2_s = SOURCE_X2_S,
      x3_x1 = SOURCE_X3_X1,
      x3_s = SOURCE_X3_S
  ) {
    set.seed(seed)
  
    s <- rbinom(n, 1L, p_s)
  
    surname <- integer(n)
    geo <- integer(n)
  
    i0 <- which(s == 0L)
    i1 <- which(s == 1L)
  
    surname[i0] <- sample_cat(surname_tab$p0, length(i0))
    surname[i1] <- sample_cat(surname_tab$p1, length(i1))
  
    geo[i0] <- sample_cat(geo_tab$p0, length(i0))
    geo[i1] <- sample_cat(geo_tab$p1, length(i1))
  
    # Ordinary nonsensitive covariates.
    # Their S-association can differ in the feature_shift target scenario.
    x1 <- x1_s * s + rnorm(n)
    x2 <- x2_s * s + rnorm(n)
    x3 <- x3_x1 * x1 + x3_s * s + rnorm(n)
  
    data.frame(
      S = s,
      surname = surname,
      geo = geo,
      x1 = x1,
      x2 = x2,
      x3 = x3
    )
  }
  
  
  # =============================================================================
  # 4. External BISG-like model: surname + geography only
  # =============================================================================
  
  estimate_bisg_like <- function(source, n_surname, n_geo, smoothing = 0.5) {
    prior <- mean(source$S)
  
    tab_name <- matrix(0, nrow = 2L, ncol = n_surname)
    tab_geo <- matrix(0, nrow = 2L, ncol = n_geo)
  
    for (sval in 0:1) {
      ds <- source[source$S == sval, , drop = FALSE]
  
      cn <- tabulate(ds$surname, nbins = n_surname) + smoothing
      cg <- tabulate(ds$geo, nbins = n_geo) + smoothing
  
      tab_name[sval + 1L, ] <- cn / sum(cn)
      tab_geo[sval + 1L, ] <- cg / sum(cg)
    }
  
    list(
      prior = prior,
      name = tab_name,
      geo = tab_geo
    )
  }
  
  predict_bisg_like <- function(dat, model) {
    p1 <- model$prior *
      model$name[2L, dat$surname] *
      model$geo[2L, dat$geo]
  
    p0 <- (1 - model$prior) *
      model$name[1L, dat$surname] *
      model$geo[1L, dat$geo]
  
    pmin(pmax(p1 / (p1 + p0), 1e-6), 1 - 1e-6)
  }
  
  
  # =============================================================================
  # 5. SAME-X design matrix
  # =============================================================================
  #
  # Crucial design choice:
  # External same-X and Target-audit use EXACTLY THE SAME covariates and the same
  # source-defined preprocessing. This makes provenance / sample size the main
  # difference between them.
  
  proxy_design_matrix_raw <- function(dat, bisg_model) {
    name_llr <- log(
      bisg_model$name[2L, dat$surname] /
        bisg_model$name[1L, dat$surname]
    )
  
    geo_llr <- log(
      bisg_model$geo[2L, dat$geo] /
        bisg_model$geo[1L, dat$geo]
    )
  
    cbind(
      intercept = 1,
      name_llr = name_llr,
      geo_llr = geo_llr,
      x1 = dat$x1,
      x2 = dat$x2,
      x3 = dat$x3
    )
  }
  
  fit_scaler <- function(X) {
    X <- as.matrix(X)
  
    if (ncol(X) <= 1L) {
      return(list(mu = numeric(0), sd = numeric(0)))
    }
  
    mu <- colMeans(X[, -1L, drop = FALSE])
    sdv <- apply(X[, -1L, drop = FALSE], 2, sd)
  
    sdv[!is.finite(sdv) | sdv < 1e-8] <- 1
  
    list(mu = mu, sd = sdv)
  }
  
  apply_scaler <- function(X, scaler) {
    X <- as.matrix(X)
  
    if (ncol(X) <= 1L) return(X)
  
    X[, -1L] <- sweep(
      sweep(X[, -1L, drop = FALSE], 2, scaler$mu, "-"),
      2,
      scaler$sd,
      "/"
    )
  
    X
  }
  
  
  # =============================================================================
  # 6. Logistic fits
  # =============================================================================
  
  ridge_logistic_fit <- function(X, y, lambda = 1e-2, maxit = 250L) {
    X <- as.matrix(X)
    y <- as.numeric(y)
  
    p <- ncol(X)
  
    # Robust fallback for tiny audit folds.
    if (length(y) < max(12L, p + 3L) || length(unique(y)) < 2L) {
      pr <- (sum(y) + 0.5) / (length(y) + 1)
      beta <- rep(0, p)
      beta[1L] <- qlogis(pr)
      return(beta)
    }
  
    fn <- function(beta) {
      eta <- drop(X %*% beta)
  
      nll <- sum(
        ifelse(
          eta >= 0,
          log1p(exp(-eta)) + (1 - y) * eta,
          log1p(exp(eta)) - y * eta
        )
      )
  
      pen <- 0.5 * lambda * sum(beta[-1L]^2)
  
      nll + pen
    }
  
    gr <- function(beta) {
      eta <- drop(X %*% beta)
      p_hat <- plogis(eta)
  
      g <- drop(crossprod(X, p_hat - y))
      g[-1L] <- g[-1L] + lambda * beta[-1L]
  
      g
    }
  
    init <- rep(0, p)
    init[1L] <- qlogis((sum(y) + 0.5) / (length(y) + 1))
  
    fit <- optim(
      init,
      fn,
      gr,
      method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-9)
    )
  
    fit$par
  }
  
  predict_logistic_beta <- function(X, beta) {
    pmin(
      pmax(plogis(drop(as.matrix(X) %*% beta)), 1e-6),
      1 - 1e-6
    )
  }
  
  
  # =============================================================================
  # 7. Honest target-audit fitting / recalibration
  # =============================================================================
  
  make_foldid <- function(n, k, seed) {
    set.seed(seed)
    sample(rep(seq_len(k), length.out = n))
  }
  
  crossfit_target_proxy <- function(X, S, R, foldid, lambda = 1e-2) {
    n <- length(S)
    q <- rep(NA_real_, n)
  
    for (k in sort(unique(foldid))) {
      te <- which(foldid == k)
      tr <- which(foldid != k & R == 1L)
  
      beta <- ridge_logistic_fit(
        X[tr, , drop = FALSE],
        S[tr],
        lambda = lambda
      )
  
      q[te] <- predict_logistic_beta(
        X[te, , drop = FALSE],
        beta
      )
    }
  
    pmin(pmax(q, 1e-6), 1 - 1e-6)
  }
  
  crossfit_recalibration <- function(q_external, S, R, foldid, lambda = 1e-3) {
    X <- cbind(
      intercept = 1,
      logit_external = logit_clip(q_external)
    )
  
    n <- length(S)
    q <- rep(NA_real_, n)
  
    for (k in sort(unique(foldid))) {
      te <- which(foldid == k)
      tr <- which(foldid != k & R == 1L)
  
      beta <- ridge_logistic_fit(
        X[tr, , drop = FALSE],
        S[tr],
        lambda = lambda
      )
  
      q[te] <- predict_logistic_beta(
        X[te, , drop = FALSE],
        beta
      )
    }
  
    pmin(pmax(q, 1e-6), 1 - 1e-6)
  }
  
  distort_probability <- function(q, alpha, gamma) {
    plogis(alpha + gamma * logit_clip(q))
  }
  
  
  # =============================================================================
  # 8. Diagnostics
  # =============================================================================
  
  auc_rank <- function(y, score) {
    y <- as.integer(y)
  
    n1 <- sum(y == 1L)
    n0 <- sum(y == 0L)
  
    if (n1 == 0L || n0 == 0L) return(NA_real_)
  
    r <- rank(score, ties.method = "average")
  
    (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  
  proxy_quality <- function(S, q) {
    data.frame(
      brier = mean((S - q)^2),
      auc = auc_rank(S, q)
    )
  }
  
  
  # =============================================================================
  # 9. Fairness moments
  # =============================================================================
  
  fairness_moment <- function(S, f) {
    mean(S * f)
  }
  
  estimate_augmented_moment <- function(q, S, R, pi, f) {
    mean(
      (
        q +
          (R / pi) * (S - q)
      ) * f
    )
  }
  
  estimate_plugin_moment <- function(q, f) {
    mean(q * f)
  }
  
  analytic_design_variance <- function(q, S, pi, f) {
    mean(
      ((1 - pi) / pi) *
        (S - q)^2 *
        f^2
    ) / length(S)
  }
  
  
  # =============================================================================
  # 10. Summaries
  # =============================================================================
  
  summarize_replications <- function(rep_df) {
    keys <- interaction(
      rep_df$scenario,
      rep_df$audit_rate,
      rep_df$method,
      drop = TRUE
    )
  
    chunks <- split(rep_df, keys)
  
    out <- lapply(chunks, function(d) {
      data.frame(
        scenario = d$scenario[1],
        audit_rate = d$audit_rate[1],
        method = d$method[1],
        n_rep = nrow(d),
        mean_estimate = mean(d$estimate),
        true_moment = d$true_moment[1],
        bias = mean(d$estimate - d$true_moment),
        sd = sd(d$estimate),
        rmse = sqrt(mean((d$estimate - d$true_moment)^2)),
        mean_brier = if (all(is.na(d$brier))) NA_real_ else mean(d$brier, na.rm = TRUE),
        mean_auc = if (all(is.na(d$auc))) NA_real_ else mean(d$auc, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  
    do.call(rbind, out)
  }
  
  
  # =============================================================================
  # 11. Source population and external models
  # =============================================================================
  
  set.seed(SEED)
  
  surname_tab <- make_binary_category_tables(
    N_SURNAME,
    NAME_STRENGTH,
    SEED + 1L
  )
  
  geo_tab <- make_binary_category_tables(
    N_GEO,
    GEO_STRENGTH,
    SEED + 2L
  )
  
  source_pop <- generate_population(
    N_SOURCE,
    P_SOURCE,
    surname_tab,
    geo_tab,
    seed = SEED + 10L,
    x1_s = SOURCE_X1_S,
    x2_s = SOURCE_X2_S,
    x3_x1 = SOURCE_X3_X1,
    x3_s = SOURCE_X3_S
  )
  
  # 11a. Canonical external BISG-like model.
  external_bisg_model <- estimate_bisg_like(
    source_pop,
    N_SURNAME,
    N_GEO,
    smoothing = SMOOTHING
  )
  
  # 11b. SAME-X external model.
  #
  # The same source-derived design transform will also be used by Target-audit.
  X_source_raw <- proxy_design_matrix_raw(
    source_pop,
    external_bisg_model
  )
  
  source_scaler <- fit_scaler(X_source_raw)
  
  X_source_sameX <- apply_scaler(
    X_source_raw,
    source_scaler
  )
  
  beta_external_sameX <- ridge_logistic_fit(
    X_source_sameX,
    source_pop$S,
    lambda = RIDGE_EXTERNAL
  )
  
  
  # =============================================================================
  # 12. Target scenarios
  # =============================================================================
  
  scenario_specs <- data.frame(
    scenario = c(
      "transported",
      "prior_shift",
      "feature_shift"
    ),
    p_target = c(
      P_TARGET_TRANSPORTED,
      P_TARGET_SHIFTED,
      P_TARGET_SHIFTED
    ),
    x1_s = c(
      SOURCE_X1_S,
      SOURCE_X1_S,
      SHIFT_X1_S
    ),
    x2_s = c(
      SOURCE_X2_S,
      SOURCE_X2_S,
      SHIFT_X2_S
    ),
    x3_x1 = c(
      SOURCE_X3_X1,
      SOURCE_X3_X1,
      SHIFT_X3_X1
    ),
    x3_s = c(
      SOURCE_X3_S,
      SOURCE_X3_S,
      SHIFT_X3_S
    ),
    stringsAsFactors = FALSE
  )
  
  
  # =============================================================================
  # 13. Main repeated-audit experiment
  # =============================================================================
  
  all_reps <- list()
  quality_rows <- list()
  
  rep_counter <- 1L
  quality_counter <- 1L
  
  for (ss in seq_len(nrow(scenario_specs))) {
  
    scenario <- scenario_specs$scenario[ss]
  
    target <- generate_population(
      N_TARGET,
      scenario_specs$p_target[ss],
      surname_tab,
      geo_tab,
      seed = SEED + 100L + 1000L * ss,
      x1_s = scenario_specs$x1_s[ss],
      x2_s = scenario_specs$x2_s[ss],
      x3_x1 = scenario_specs$x3_x1[ss],
      x3_s = scenario_specs$x3_s[ss]
    )
  
    # EXACTLY THE SAME feature representation for:
    #   - External same-X
    #   - Target-audit
    X_target_raw <- proxy_design_matrix_raw(
      target,
      external_bisg_model
    )
  
    X_target_sameX <- apply_scaler(
      X_target_raw,
      source_scaler
    )
  
    # Fixed external SAME-X prediction.
    q_external_sameX <- predict_logistic_beta(
      X_target_sameX,
      beta_external_sameX
    )
  
    # Deliberately distorted SAME-X score.
    q_external_sameX_distorted <- distort_probability(
      q_external_sameX,
      MIS_ALPHA,
      MIS_GAMMA
    )
  
    # External BISG-like: name + geography only.
    q_external_bisg <- predict_bisg_like(
      target,
      external_bisg_model
    )
  
    # Fixed favorable-decision score.
    #
    # It is measurable with respect to the SAME-X feature set. Therefore, in the
    # transported scenario, a well estimated E[S | SAME-X] is the natural
    # plug-in benchmark. BISG-like sees only part of this information.
    score <- plogis(
      -0.35 +
        0.72 * target$x1 -
        0.52 * target$x2 +
        0.15 * target$x3 +
        0.24 * X_target_sameX[, "name_llr"] +
        0.18 * X_target_sameX[, "geo_llr"]
    )
  
    f <- score - mean(score)
  
    true_m <- fairness_moment(
      target$S,
      f
    )
  
    # Fixed proxy quality, independent of target audit.
    fixed_proxies <- list(
      "External same-X" = q_external_sameX,
      "External same-X distorted" = q_external_sameX_distorted,
      "External BISG-like" = q_external_bisg
    )
  
    for (nm in names(fixed_proxies)) {
      qq <- proxy_quality(
        target$S,
        fixed_proxies[[nm]]
      )
  
      quality_rows[[quality_counter]] <- data.frame(
        scenario = scenario,
        audit_rate = NA_real_,
        proxy = nm,
        brier = qq$brier,
        auc = qq$auc,
        stringsAsFactors = FALSE
      )
  
      quality_counter <- quality_counter + 1L
    }
  
    for (ar in AUDIT_RATES) {
  
      message(
        sprintf(
          "same-X provenance: %s, audit %.1f%%",
          scenario,
          100 * ar
        )
      )
  
      for (b in seq_len(N_REP)) {
  
        set.seed(
          SEED +
            100000L * ss +
            round(10000 * ar) +
            b
        )
  
        R <- rbinom(
          N_TARGET,
          1L,
          ar
        )
  
        foldid <- make_foldid(
          N_TARGET,
          K_FOLDS,
          seed =
            SEED +
            200000L * ss +
            round(10000 * ar) +
            b
        )
  
        # Locally learned target-audit model, SAME X as external same-X.
        q_target <- crossfit_target_proxy(
          X_target_sameX,
          target$S,
          R,
          foldid,
          lambda = RIDGE_TARGET
        )
  
        # Recalibrate the external same-X model on the small target audit.
        q_external_sameX_recal <- crossfit_recalibration(
          q_external_sameX,
          target$S,
          R,
          foldid,
          lambda = RIDGE_RECAL
        )
  
        proxies <- list(
          "Target-audit same-X" = q_target,
          "External same-X" = q_external_sameX,
          "External same-X distorted" = q_external_sameX_distorted,
          "External same-X recalibrated" = q_external_sameX_recal,
          "External BISG-like" = q_external_bisg
        )
  
        # --------------------------- IPW baseline -------------------------------
        est_ipw <- estimate_augmented_moment(
          rep(0, N_TARGET),
          target$S,
          R,
          ar,
          f
        )
  
        all_reps[[rep_counter]] <- data.frame(
          scenario = scenario,
          audit_rate = ar,
          replicate = b,
          method = "IPW",
          estimate = est_ipw,
          true_moment = true_m,
          brier = NA_real_,
          auc = NA_real_,
          stringsAsFactors = FALSE
        )
  
        rep_counter <- rep_counter + 1L
  
        # ------------------------ proxy-based methods ---------------------------
        for (nm in names(proxies)) {
  
          q <- proxies[[nm]]
  
          qq <- proxy_quality(
            target$S,
            q
          )
  
          # Augmented.
          est_aug <- estimate_augmented_moment(
            q,
            target$S,
            R,
            ar,
            f
          )
  
          all_reps[[rep_counter]] <- data.frame(
            scenario = scenario,
            audit_rate = ar,
            replicate = b,
            method = paste0("Augmented: ", nm),
            estimate = est_aug,
            true_moment = true_m,
            brier = qq$brier,
            auc = qq$auc,
            stringsAsFactors = FALSE
          )
  
          rep_counter <- rep_counter + 1L
  
          # Direct plug-in.
          est_plugin <- estimate_plugin_moment(
            q,
            f
          )
  
          all_reps[[rep_counter]] <- data.frame(
            scenario = scenario,
            audit_rate = ar,
            replicate = b,
            method = paste0("Plug-in: ", nm),
            estimate = est_plugin,
            true_moment = true_m,
            brier = qq$brier,
            auc = qq$auc,
            stringsAsFactors = FALSE
          )
  
          rep_counter <- rep_counter + 1L
        }
      }
    }
  }
  
  
  # =============================================================================
  # 14. Summaries and CSV outputs
  # =============================================================================
  
  rep_df <- do.call(
    rbind,
    all_reps
  )
  
  summary_df <- summarize_replications(
    rep_df
  )
  
  quality_df <- do.call(
    rbind,
    quality_rows
  )
  
  write.csv(
    rep_df,
    file.path(out_dir, "sameX_proxy_provenance_replicates.csv"),
    row.names = FALSE
  )
  
  write.csv(
    summary_df,
    file.path(out_dir, "sameX_proxy_provenance_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    quality_df,
    file.path(out_dir, "sameX_external_proxy_quality.csv"),
    row.names = FALSE
  )
  
  
  # =============================================================================
  # 15. Calibration-only sensitivity for the external SAME-X score
  # =============================================================================
  #
  # This is deliberately done in the transported scenario.
  # alpha shifts the logit intercept only:
  #
  #   q_alpha = expit(alpha + logit(q_external_sameX)).
  #
  # Hence AUC is EXACTLY unchanged. Calibration, plug-in bias, and augmented
  # variance may nevertheless change.
  
  cal_target <- generate_population(
    N_TARGET,
    P_TARGET_TRANSPORTED,
    surname_tab,
    geo_tab,
    seed = SEED + 999L,
    x1_s = SOURCE_X1_S,
    x2_s = SOURCE_X2_S,
    x3_x1 = SOURCE_X3_X1,
    x3_s = SOURCE_X3_S
  )
  
  X_cal_raw <- proxy_design_matrix_raw(
    cal_target,
    external_bisg_model
  )
  
  X_cal_sameX <- apply_scaler(
    X_cal_raw,
    source_scaler
  )
  
  q0_sameX <- predict_logistic_beta(
    X_cal_sameX,
    beta_external_sameX
  )
  
  score_cal <- plogis(
    -0.35 +
      0.72 * cal_target$x1 -
      0.52 * cal_target$x2 +
      0.15 * cal_target$x3 +
      0.24 * X_cal_sameX[, "name_llr"] +
      0.18 * X_cal_sameX[, "geo_llr"]
  )
  
  f_cal <- score_cal - mean(score_cal)
  
  true_cal <- fairness_moment(
    cal_target$S,
    f_cal
  )
  
  pi_cal <- 0.05
  
  var_ipw <- analytic_design_variance(
    rep(0, N_TARGET),
    cal_target$S,
    pi_cal,
    f_cal
  )
  
  alpha_grid <- seq(
    -2.0,
    2.0,
    by = 0.20
  )
  
  cal_sens <- lapply(
    alpha_grid,
    function(a) {
  
      q <- distort_probability(
        q0_sameX,
        alpha = a,
        gamma = 1
      )
  
      qq <- proxy_quality(
        cal_target$S,
        q
      )
  
      data.frame(
        alpha = a,
        brier = qq$brier,
        auc = qq$auc,
        plugin_bias =
          estimate_plugin_moment(q, f_cal) -
          true_cal,
        variance_ratio_to_ipw =
          analytic_design_variance(
            q,
            cal_target$S,
            pi_cal,
            f_cal
          ) / var_ipw,
        stringsAsFactors = FALSE
      )
    }
  )
  
  cal_sens <- do.call(
    rbind,
    cal_sens
  )
  
  write.csv(
    cal_sens,
    file.path(out_dir, "sameX_external_calibration_sensitivity.csv"),
    row.names = FALSE
  )
  
  
  # =============================================================================
  # 16. Plots
  # =============================================================================
  
  library(ggplot2)
  
  # --- 16a. RMSE: provenance comparison ----------------------------------------
  
  rmse_methods <- c(
    "IPW",
    "Augmented: Target-audit same-X",
    "Augmented: External same-X",
    "Augmented: External same-X recalibrated",
    "Augmented: External BISG-like"
  )
  
  rmse_df <- summary_df[
    summary_df$method %in% rmse_methods,
    ,
    drop = FALSE
  ]
  
  p_rmse <- ggplot(
    rmse_df,
    aes(
      x = 100 * audit_rate,
      y = rmse,
      color = method,
      shape = method
    )
  ) +
    geom_line() +
    geom_point(size = 2.1) +
    facet_wrap(
      ~ scenario,
      scales = "free_y",
      nrow = 1
    ) +
    labs(
      x = "Audit rate (%)",
      y = "RMSE of the fairness moment",
      color = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )
  
  ggsave(
    file.path(out_dir, "figure_sameX_provenance_rmse.pdf"),
    p_rmse,
    width = 11,
    height = 4.6
  )
  
  
  # --- 16b. External-vs-target RMSE ratio ---------------------------------------
  
  ratio_rows <- list()
  ratio_counter <- 1L
  
  for (sc in unique(summary_df$scenario)) {
    for (ar in AUDIT_RATES) {
  
      d_target <- summary_df[
        summary_df$scenario == sc &
          abs(summary_df$audit_rate - ar) < 1e-12 &
          summary_df$method == "Augmented: Target-audit same-X",
        ,
        drop = FALSE
      ]
  
      d_ext <- summary_df[
        summary_df$scenario == sc &
          abs(summary_df$audit_rate - ar) < 1e-12 &
          summary_df$method == "Augmented: External same-X",
        ,
        drop = FALSE
      ]
  
      if (nrow(d_target) == 1L && nrow(d_ext) == 1L) {
        ratio_rows[[ratio_counter]] <- data.frame(
          scenario = sc,
          audit_rate = ar,
          rmse_ratio_external_over_target =
            d_ext$rmse / d_target$rmse
        )
  
        ratio_counter <- ratio_counter + 1L
      }
    }
  }
  
  ratio_df <- do.call(
    rbind,
    ratio_rows
  )
  
  write.csv(
    ratio_df,
    file.path(out_dir, "sameX_external_vs_target_rmse_ratio.csv"),
    row.names = FALSE
  )
  
  p_ratio <- ggplot(
    ratio_df,
    aes(
      x = 100 * audit_rate,
      y = rmse_ratio_external_over_target,
      color = scenario,
      shape = scenario
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    geom_line() +
    geom_point(size = 2.1) +
    labs(
      x = "Audit rate (%)",
      y = "RMSE ratio: external same-X / target-audit same-X",
      color = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom"
    )
  
  ggsave(
    file.path(out_dir, "figure_sameX_external_vs_target_ratio.pdf"),
    p_ratio,
    width = 7.0,
    height = 4.4
  )
  
  
  # --- 16c. Plug-in vs augmented bias -------------------------------------------
  
  bias_methods <- c(
    "IPW",
    "Plug-in: External same-X",
    "Plug-in: External same-X distorted",
    "Plug-in: External same-X recalibrated",
    "Plug-in: External BISG-like",
    "Augmented: External same-X",
    "Augmented: External same-X distorted",
    "Augmented: External same-X recalibrated",
    "Augmented: External BISG-like"
  )
  
  bias_df <- summary_df[
    summary_df$method %in% bias_methods &
      summary_df$scenario %in% c("prior_shift", "feature_shift"),
    ,
    drop = FALSE
  ]
  
  p_bias <- ggplot(
    bias_df,
    aes(
      x = 100 * audit_rate,
      y = bias,
      color = method,
      shape = method
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4
    ) +
    geom_line() +
    geom_point(size = 1.8) +
    facet_wrap(
      ~ scenario,
      scales = "free_y"
    ) +
    labs(
      x = "Audit rate (%)",
      y = "Signed bias of the fairness moment",
      color = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )
  
  ggsave(
    file.path(out_dir, "figure_sameX_external_proxy_bias.pdf"),
    p_bias,
    width = 10.0,
    height = 5.4
  )
  
  
  # --- 16d. Proxy Brier by audit rate -------------------------------------------
  
  quality_methods <- c(
    "Augmented: Target-audit same-X",
    "Augmented: External same-X",
    "Augmented: External same-X recalibrated",
    "Augmented: External BISG-like"
  )
  
  quality_plot_df <- summary_df[
    summary_df$method %in% quality_methods,
    ,
    drop = FALSE
  ]
  
  p_quality <- ggplot(
    quality_plot_df,
    aes(
      x = 100 * audit_rate,
      y = mean_brier,
      color = method,
      shape = method
    )
  ) +
    geom_line() +
    geom_point(size = 2.0) +
    facet_wrap(
      ~ scenario,
      scales = "free_y",
      nrow = 1
    ) +
    labs(
      x = "Audit rate (%)",
      y = "Brier score of the auxiliary demographic model",
      color = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )
  
  ggsave(
    file.path(out_dir, "figure_sameX_proxy_brier.pdf"),
    p_quality,
    width = 11,
    height = 4.6
  )
  
  
  # --- 16e. Calibration sensitivity: variance ----------------------------------
  
  p_cal_var <- ggplot(
    cal_sens,
    aes(
      x = alpha,
      y = variance_ratio_to_ipw
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    geom_line() +
    geom_point(size = 1.8) +
    labs(
      x = "External same-X logit calibration shift (alpha)",
      y = "Analytic augmented variance / IPW variance (5% audit)"
    ) +
    theme_minimal(base_size = 11)
  
  ggsave(
    file.path(out_dir, "figure_sameX_calibration_variance.pdf"),
    p_cal_var,
    width = 6.6,
    height = 4.3
  )
  
  
  # --- 16f. Calibration sensitivity: plug-in bias -------------------------------
  
  p_cal_bias <- ggplot(
    cal_sens,
    aes(
      x = alpha,
      y = plugin_bias
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4
    ) +
    geom_line() +
    geom_point(size = 1.8) +
    labs(
      x = "External same-X logit calibration shift (alpha)",
      y = "Plug-in fairness-moment bias"
    ) +
    theme_minimal(base_size = 11)
  
  ggsave(
    file.path(out_dir, "figure_sameX_calibration_plugin_bias.pdf"),
    p_cal_bias,
    width = 6.6,
    height = 4.3
  )
  
  
  # =============================================================================
  # 17. Compact diagnostic tables
  # =============================================================================
  
  lookup <- function(scenario, audit_rate, method, field) {
    d <- summary_df[
      summary_df$scenario == scenario &
        abs(summary_df$audit_rate - audit_rate) < 1e-12 &
        summary_df$method == method,
      ,
      drop = FALSE
    ]
  
    if (nrow(d) != 1L) return(NA_real_)
  
    d[[field]][1]
  }
  
  diagnostic_rows <- list()
  diag_counter <- 1L
  
  for (sc in unique(summary_df$scenario)) {
    for (ar in AUDIT_RATES) {
  
      diagnostic_rows[[diag_counter]] <- data.frame(
        scenario = sc,
        audit_rate = ar,
        rmse_ipw = lookup(
          sc,
          ar,
          "IPW",
          "rmse"
        ),
        rmse_target = lookup(
          sc,
          ar,
          "Augmented: Target-audit same-X",
          "rmse"
        ),
        rmse_external_sameX = lookup(
          sc,
          ar,
          "Augmented: External same-X",
          "rmse"
        ),
        rmse_external_recal = lookup(
          sc,
          ar,
          "Augmented: External same-X recalibrated",
          "rmse"
        ),
        rmse_bisg = lookup(
          sc,
          ar,
          "Augmented: External BISG-like",
          "rmse"
        ),
        brier_target = lookup(
          sc,
          ar,
          "Augmented: Target-audit same-X",
          "mean_brier"
        ),
        brier_external_sameX = lookup(
          sc,
          ar,
          "Augmented: External same-X",
          "mean_brier"
        ),
        brier_recal = lookup(
          sc,
          ar,
          "Augmented: External same-X recalibrated",
          "mean_brier"
        ),
        brier_bisg = lookup(
          sc,
          ar,
          "Augmented: External BISG-like",
          "mean_brier"
        ),
        stringsAsFactors = FALSE
      )
  
      diag_counter <- diag_counter + 1L
    }
  }
  
  diagnostic_df <- do.call(
    rbind,
    diagnostic_rows
  )
  
  write.csv(
    diagnostic_df,
    file.path(out_dir, "sameX_key_diagnostics.csv"),
    row.names = FALSE
  )
  
  
  # =============================================================================
  # 18. Text diagnostic for quick inspection
  # =============================================================================
  
  fmt <- function(x, digits = 5L) {
    if (!is.finite(x)) return("NA")
    formatC(x, digits = digits, format = "f")
  }
  
  lines <- c(
    "SAME-X PROXY PROVENANCE DIAGNOSTIC",
    "==================================",
    "",
    sprintf(
      "N_source = %d; N_target = %d; repeated audits per cell = %d",
      N_SOURCE,
      N_TARGET,
      N_REP
    ),
    ""
  )
  
  for (sc in c("transported", "prior_shift", "feature_shift")) {
  
    lines <- c(
      lines,
      paste0("Scenario: ", sc),
      paste0(
        "  1% RMSE -- IPW: ",
        fmt(lookup(sc, 0.01, "IPW", "rmse")),
        "; target: ",
        fmt(lookup(sc, 0.01, "Augmented: Target-audit same-X", "rmse")),
        "; external same-X: ",
        fmt(lookup(sc, 0.01, "Augmented: External same-X", "rmse")),
        "; recalibrated: ",
        fmt(lookup(sc, 0.01, "Augmented: External same-X recalibrated", "rmse")),
        "; BISG-like: ",
        fmt(lookup(sc, 0.01, "Augmented: External BISG-like", "rmse"))
      ),
      paste0(
        " 20% RMSE -- IPW: ",
        fmt(lookup(sc, 0.20, "IPW", "rmse")),
        "; target: ",
        fmt(lookup(sc, 0.20, "Augmented: Target-audit same-X", "rmse")),
        "; external same-X: ",
        fmt(lookup(sc, 0.20, "Augmented: External same-X", "rmse")),
        "; recalibrated: ",
        fmt(lookup(sc, 0.20, "Augmented: External same-X recalibrated", "rmse")),
        "; BISG-like: ",
        fmt(lookup(sc, 0.20, "Augmented: External BISG-like", "rmse"))
      ),
      paste0(
        "  5% plug-in bias -- external same-X: ",
        fmt(lookup(sc, 0.05, "Plug-in: External same-X", "bias")),
        "; distorted: ",
        fmt(lookup(sc, 0.05, "Plug-in: External same-X distorted", "bias")),
        "; recalibrated: ",
        fmt(lookup(sc, 0.05, "Plug-in: External same-X recalibrated", "bias")),
        "; BISG-like: ",
        fmt(lookup(sc, 0.05, "Plug-in: External BISG-like", "bias"))
      ),
      ""
    )
  }
  
  writeLines(
    lines,
    file.path(out_dir, "sameX_quick_diagnostic.txt")
  )
  
  cat(
    paste(lines, collapse = "\n"),
    "\n"
  )
  
  
  # =============================================================================
  # 19. LaTeX paragraph and compact appendix table
  # =============================================================================
  
  ar <- 0.05
  
  body_txt <- sprintf(
    paste0(
      "At a 5\\%% target audit rate in the transported scenario, fairness-moment RMSE is ",
      "%.5f for IPW, %.5f for augmentation with a proxy learned from the target audit, ",
      "%.5f for augmentation with an external model trained on the same covariates, ",
      "and %.5f for the external BISG-like proxy. ",
      "Under target prior shift, direct substitution of the external same-covariate score has signed bias %.5f, ",
      "whereas its audit-augmented version has bias %.5f; after audit-based recalibration the augmented bias is %.5f. ",
      "The feature-shift scenario separately tests local adaptation when the conditional relation between ordinary ",
      "covariates and the protected attribute changes between source and target populations. ",
      "These experiments isolate external statistical strength from target-population adaptation while retaining ",
      "BISG as a canonical example of an externally constructed demographic score.\n"
    ),
    lookup("transported", ar, "IPW", "rmse"),
    lookup("transported", ar, "Augmented: Target-audit same-X", "rmse"),
    lookup("transported", ar, "Augmented: External same-X", "rmse"),
    lookup("transported", ar, "Augmented: External BISG-like", "rmse"),
    lookup("prior_shift", ar, "Plug-in: External same-X", "bias"),
    lookup("prior_shift", ar, "Augmented: External same-X", "bias"),
    lookup("prior_shift", ar, "Augmented: External same-X recalibrated", "bias")
  )
  
  writeLines(
    body_txt,
    file.path(out_dir, "sameX_proxy_provenance_body.tex")
  )
  
  
  # Compact table: 5% audit, all three target regimes.
  tab_methods <- c(
    "IPW",
    "Augmented: Target-audit same-X",
    "Augmented: External same-X",
    "Augmented: External same-X recalibrated",
    "Augmented: External BISG-like",
    "Plug-in: External same-X",
    "Plug-in: External BISG-like"
  )
  
  tab <- summary_df[
    abs(summary_df$audit_rate - 0.05) < 1e-12 &
      summary_df$method %in% tab_methods,
    c(
      "scenario",
      "method",
      "bias",
      "sd",
      "rmse",
      "mean_brier",
      "mean_auc"
    )
  ]
  
  con <- file(
    file.path(out_dir, "sameX_proxy_provenance_table.tex"),
    open = "wt"
  )
  
  writeLines("\\begin{table*}[t]", con)
  writeLines("\\centering\\small", con)
  writeLines(
    "\\caption{Proxy-provenance simulation at a 5\\% target audit rate. External same-X and target-audit proxies use identical predictors; BISG-like uses surname and geography only.}",
    con
  )
  writeLines("\\label{tab:sameX-proxy-provenance}", con)
  writeLines("\\begin{tabular}{llrrrr}", con)
  writeLines("\\toprule", con)
  writeLines("Scenario & Method & Bias & SD & RMSE & Brier \\\\", con)
  writeLines("\\midrule", con)
  
  for (i in seq_len(nrow(tab))) {
  
    br <- ifelse(
      is.finite(tab$mean_brier[i]),
      sprintf("%.4f", tab$mean_brier[i]),
      "--"
    )
  
    method_tex <- gsub(
      "_",
      "\\\\_",
      tab$method[i],
      fixed = TRUE
    )
  
    line <- sprintf(
      "%s & %s & %.5f & %.5f & %.5f & %s \\\\",
      tab$scenario[i],
      method_tex,
      tab$bias[i],
      tab$sd[i],
      tab$rmse[i],
      br
    )
  
    writeLines(line, con)
  }
  
  writeLines("\\bottomrule", con)
  writeLines("\\end{tabular}", con)
  writeLines("\\end{table*}", con)
  
  close(con)
  
  
  # =============================================================================
  # 20. Final message
  # =============================================================================
  
  message(
    "\nSimulation complete. Outputs written to:\n  ",
    out_dir
  )
  
  message(
    "N_source = ",
    N_SOURCE,
    "; N_target = ",
    N_TARGET,
    "; repeated audits per cell = ",
    N_REP
  )

  # Copy manuscript figure artifacts to the common figures/ directory.
  provenance_figures <- list.files(
    out_dir,
    pattern = "^figure_sameX_.*\\.pdf$",
    full.names = TRUE
  )
  if (length(provenance_figures) > 0L) {
    file.copy(
      provenance_figures,
      figures_dir,
      overwrite = TRUE
    )
  }

  invisible(list(
    replicates = rep_df,
    summary = summary_df,
    quality = quality_df,
    calibration = cal_sens,
    ratio = ratio_df,
    diagnostics = diagnostic_df,
    out_dir = normalizePath(out_dir, mustWork = FALSE),
    figures = basename(provenance_figures),
    n_source = N_SOURCE,
    n_target = N_TARGET,
    n_rep = N_REP
  ))
}


save_robustness_figures <- function(
    stress,
    variance_check,
    exact_convex,
    neyman_external,
    neyman_total,
    selection,
    figures_dir = "figures") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to save robustness figures.")
  }
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  files <- character(0)

  # 1. Proxy stress: variance ratio and bias.
  d <- stress$variance
  d$audit_percent <- 100 * d$audit_rate
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = proxy_case,
      y = variance_ratio_to_ipw,
      group = factor(audit_percent),
      linetype = factor(audit_percent),
      shape = factor(audit_percent)
    )
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Auxiliary proxy",
      y = "Variance ratio: augmented / IPW",
      linetype = "Audit rate (%)",
      shape = "Audit rate (%)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 28, hjust = 1)
    )
  f <- file.path(figures_dir, "figure_proxy_stress_variance_ratio_v51.pdf")
  save_pdf_plot(p, f, width = 7.5, height = 4.4)
  files <- c(files, f)

  bias <- stats::aggregate(
    error ~ audit_rate + proxy_case,
    data = stress$results,
    FUN = mean
  )
  bias$audit_percent <- 100 * bias$audit_rate
  p <- ggplot2::ggplot(
    bias,
    ggplot2::aes(
      x = proxy_case,
      y = error,
      group = factor(audit_percent),
      linetype = factor(audit_percent),
      shape = factor(audit_percent)
    )
  ) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Auxiliary proxy",
      y = "Signed bias of the fairness moment",
      linetype = "Audit rate (%)",
      shape = "Audit rate (%)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 28, hjust = 1)
    )
  f <- file.path(figures_dir, "figure_proxy_stress_bias_v51.pdf")
  save_pdf_plot(p, f, width = 7.5, height = 4.2)
  files <- c(files, f)

  # 2. Exact versus convex.
  d <- exact_convex$frontiers
  frontier_summary <- stats::aggregate(
    cbind(log_loss, dp_score) ~ method + path_type + lambda,
    data = d,
    FUN = mean
  )
  frontier_summary <- frontier_summary[
    order(
      frontier_summary$method,
      frontier_summary$path_type,
      frontier_summary$lambda
    ),
  ]
  p <- ggplot2::ggplot(
    frontier_summary,
    ggplot2::aes(
      x = dp_score,
      y = log_loss,
      linetype = path_type,
      shape = method,
      group = interaction(method, path_type)
    )
  ) +
    ggplot2::geom_path() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Test score demographic-parity gap",
      y = "Test log loss",
      linetype = "Objective",
      shape = NULL
    ) +
    ggplot2::theme_minimal()
  f <- file.path(figures_dir, "figure_exact_vs_convex_frontiers_v51.pdf")
  save_pdf_plot(p, f, width = 6.8, height = 4.2)
  files <- c(files, f)

  moment_long <- rbind(
    data.frame(
      method = d$method,
      path_type = d$path_type,
      lambda = d$lambda,
      moment_type = "Exact score moment",
      relative_moment = d$relative_training_exact_moment
    ),
    data.frame(
      method = d$method,
      path_type = d$path_type,
      lambda = d$lambda,
      moment_type = "Surrogate logit moment",
      relative_moment = d$relative_training_surrogate_moment
    )
  )
  positive <- moment_long$lambda[moment_long$lambda > 0]
  min_positive <- min(positive)
  moment_long$log10_lambda_plot <- ifelse(
    moment_long$lambda == 0,
    log10(min_positive) - 1,
    log10(moment_long$lambda)
  )
  moment_summary <- stats::aggregate(
    relative_moment ~
      method + path_type + log10_lambda_plot + moment_type,
    data = subset(
      moment_long,
      method %in% c("Oracle", "Augmented")
    ),
    FUN = mean
  )
  p <- ggplot2::ggplot(
    moment_summary,
    ggplot2::aes(
      x = log10_lambda_plot,
      y = relative_moment,
      linetype = moment_type,
      group = moment_type
    )
  ) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted") +
    ggplot2::geom_line() +
    ggplot2::facet_grid(method ~ path_type) +
    ggplot2::labs(
      x = "log10(lambda), with lambda=0 shown at left",
      y = "Moment relative to lambda = 0",
      linetype = NULL
    ) +
    ggplot2::theme_minimal()
  f <- file.path(figures_dir, "figure_exact_vs_convex_moment_paths_v51.pdf")
  save_pdf_plot(p, f, width = 7.2, height = 5.2)
  files <- c(files, f)

  # 3. Variance validation.
  v <- variance_check$summary
  v$audit_percent <- 100 * v$audit_rate
  long <- rbind(
    data.frame(
      audit_percent = v$audit_percent,
      method = v$method,
      variance_type = "Empirical",
      value = v$empirical_variance
    ),
    data.frame(
      audit_percent = v$audit_percent,
      method = v$method,
      variance_type = "Estimated",
      value = v$mean_estimated_variance
    ),
    data.frame(
      audit_percent = v$audit_percent,
      method = v$method,
      variance_type = "Exact design",
      value = v$exact_design_variance
    )
  )
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = audit_percent,
      y = value,
      linetype = variance_type,
      shape = variance_type,
      group = variance_type
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ method, scales = "free_y") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Audit rate (%)",
      y = "Variance (log scale)",
      linetype = NULL,
      shape = NULL
    ) +
    ggplot2::theme_minimal()
  f <- file.path(figures_dir, "figure_variance_validation_v51.pdf")
  save_pdf_plot(p, f, width = 6.8, height = 4.1)
  files <- c(files, f)

  # 4. External-pilot Neyman.
  n <- neyman_external$summary
  n$audit_percent <- 100 * n$audit_rate
  n$pilot_percent <- factor(100 * n$pilot_rate)
  p <- ggplot2::ggplot(
    subset(n, design != "Uniform"),
    ggplot2::aes(
      x = audit_percent,
      y = variance_ratio_to_uniform,
      linetype = design,
      shape = design,
      group = design
    )
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ pilot_percent) +
    ggplot2::labs(
      x = "Main-audit rate (%)",
      y = "Variance ratio to uniform main audit",
      linetype = NULL,
      shape = NULL
    ) +
    ggplot2::theme_minimal()
  f <- file.path(figures_dir, "figure_neyman_external_pilot_v51.pdf")
  save_pdf_plot(p, f, width = 7.2, height = 4.7)
  files <- c(files, f)

  # 5. Counted-budget Neyman.
  n <- neyman_total$summary
  n$total_percent <- 100 * n$total_audit_rate
  n$pilot_fraction_percent <- factor(100 * n$pilot_fraction)
  p <- ggplot2::ggplot(
    subset(
      n,
      design %in% c(
        "TwoPhasePilotNeyman",
        "TwoPhaseOracleNeyman"
      )
    ),
    ggplot2::aes(
      x = total_percent,
      y = variance_ratio_to_uniform,
      linetype = design,
      shape = design,
      group = design
    )
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ pilot_fraction_percent) +
    ggplot2::labs(
      x = "Total audit budget (%)",
      y = "Variance ratio to two-phase uniform allocation",
      linetype = NULL,
      shape = NULL
    ) +
    ggplot2::theme_minimal()
  f <- file.path(figures_dir, "figure_neyman_total_budget_v51.pdf")
  save_pdf_plot(p, f, width = 7.2, height = 4.7)
  files <- c(files, f)

  # 6. Validation-audit selection robustness.
  s <- subset(
    selection$summary,
    selection_rule %in% c(
      "ValidationAuditPoint",
      "ValidationAuditUCB"
    )
  )
  if (nrow(s) > 0L) {
    s$validation_percent <- 100 * s$validation_audit_rate

    p <- ggplot2::ggplot(
      s,
      ggplot2::aes(
        x = validation_percent,
        y = dp_score_mean,
        linetype = selection_rule,
        shape = selection_rule,
        group = selection_rule
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~ method) +
      ggplot2::labs(
        x = "Independent validation-audit rate (%)",
        y = "Test score demographic-parity gap",
        linetype = NULL,
        shape = NULL
      ) +
      ggplot2::theme_minimal()
    f <- file.path(figures_dir, "figure_selection_validation_dp_v51.pdf")
    save_pdf_plot(p, f, width = 7.0, height = 4.2)
    files <- c(files, f)

    p <- ggplot2::ggplot(
      s,
      ggplot2::aes(
        x = validation_percent,
        y = feasible_numeric_mean,
        linetype = selection_rule,
        shape = selection_rule,
        group = selection_rule
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~ method) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        x = "Independent validation-audit rate (%)",
        y = "Fraction attaining the 50% target",
        linetype = NULL,
        shape = NULL
      ) +
      ggplot2::theme_minimal()
    f <- file.path(figures_dir, "figure_selection_feasibility_v51.pdf")
    save_pdf_plot(p, f, width = 7.0, height = 4.2)
    files <- c(files, f)
  }

  invisible(files)
}
