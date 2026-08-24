# Data preparation for the UCI Default of Credit Card Clients data.

resolve_project_root <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE)))
  }

  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) {
    if (!is.null(frame$ofile)) as.character(frame$ofile) else NA_character_
  }, character(1L))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles) > 0L) {
    # When this helper is sourced by a top-level script, use the directory of
    # the top-level script rather than the R/ subdirectory.
    candidate <- dirname(normalizePath(tail(ofiles, 1L), mustWork = FALSE))
    if (basename(candidate) == "R") candidate <- dirname(candidate)
    return(candidate)
  }
  getwd()
}

load_credit_application_data <- function(
    root,
    sensitive = c("sex", "age"),
    seed_split = 20260801L) {
  sensitive <- match.arg(tolower(sensitive), c("sex", "age"))
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Install readxl before running: install.packages('readxl')")
  }

  data_dir <- file.path(root, "data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  zip_url <- paste0(
    "https://archive.ics.uci.edu/static/public/350/",
    "default%2Bof%2Bcredit%2Bcard%2Bclients.zip"
  )
  zip_path <- file.path(data_dir, "default_credit_clients.zip")
  xls_path <- file.path(data_dir, "default of credit card clients.xls")

  if (!file.exists(xls_path)) {
    if (!file.exists(zip_path)) {
      message("Downloading the UCI Credit Default dataset...")
      utils::download.file(zip_url, zip_path, mode = "wb", quiet = FALSE)
    }
    utils::unzip(zip_path, exdir = data_dir)
    candidates <- list.files(
      data_dir,
      pattern = "\\.xls$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(candidates) == 0L) {
      stop("The UCI archive did not contain an .xls file.")
    }
    xls_path <- candidates[1L]
  }

  credit <- as.data.frame(readxl::read_excel(xls_path, skip = 1L))
  names(credit) <- toupper(gsub("[^A-Za-z0-9]+", "_", names(credit)))
  names(credit) <- sub("_$", "", names(credit))
  target_candidates <- grep(
    "DEFAULT.*PAYMENT.*NEXT.*MONTH",
    names(credit),
    value = TRUE
  )
  if (length(target_candidates) != 1L) {
    stop(
      "Could not identify the default-payment outcome column. Columns are: ",
      paste(names(credit), collapse = ", ")
    )
  }
  target_name <- target_candidates[1L]
  y <- 1L - as.integer(credit[[target_name]])

  if (sensitive == "age") {
    age_cut <- stats::median(credit$AGE, na.rm = TRUE)
    s <- as.integer(credit$AGE >= age_cut)
    sensitive_description <- sprintf(
      "age at or above the sample median (%g years)",
      age_cut
    )
    excluded_sensitive_columns <- "AGE"
  } else {
    s <- as.integer(credit$SEX == 2L)
    sensitive_description <- "female sex"
    excluded_sensitive_columns <- "SEX"
  }

  signed_log1p <- function(x) sign(x) * log1p(abs(x))
  feature_data <- credit
  feature_data[[target_name]] <- NULL
  feature_data$ID <- NULL
  feature_data[excluded_sensitive_columns] <- NULL

  for (variable in intersect(c("SEX", "EDUCATION", "MARRIAGE"), names(feature_data))) {
    feature_data[[variable]] <- factor(feature_data[[variable]])
  }
  amount_variables <- grep(
    "^(LIMIT_BAL|BILL_AMT|PAY_AMT)",
    names(feature_data),
    value = TRUE
  )
  for (variable in amount_variables) {
    feature_data[[variable]] <- signed_log1p(as.numeric(feature_data[[variable]]))
  }

  X_all <- stats::model.matrix(~ ., data = feature_data)
  split <- stratified_split(y, s, seed = seed_split)
  index <- list(
    train = which(split == "train"),
    validation = which(split == "validation"),
    test = which(split == "test")
  )
  scaled <- standardize_from_train(X_all, index$train, intercept = TRUE)
  X_all <- scaled$X

  misspecified_names <- unique(c(
    "(Intercept)",
    grep("LIMIT_BAL", colnames(X_all), value = TRUE),
    grep("^PAY_0$|^PAY_2$", colnames(X_all), value = TRUE)
  ))
  misspecified_names <- intersect(misspecified_names, colnames(X_all))

  list(
    credit = credit,
    X_all = X_all,
    y = y,
    s = s,
    split = split,
    index = index,
    X = list(
      train = X_all[index$train, , drop = FALSE],
      validation = X_all[index$validation, , drop = FALSE],
      test = X_all[index$test, , drop = FALSE]
    ),
    y_split = list(
      train = y[index$train],
      validation = y[index$validation],
      test = y[index$test]
    ),
    s_split = list(
      train = s[index$train],
      validation = s[index$validation],
      test = s[index$test]
    ),
    X_proxy_rich = list(
      train = X_all[index$train, , drop = FALSE],
      validation = X_all[index$validation, , drop = FALSE],
      test = X_all[index$test, , drop = FALSE]
    ),
    X_proxy_misspecified = list(
      train = X_all[index$train, misspecified_names, drop = FALSE],
      validation = X_all[index$validation, misspecified_names, drop = FALSE],
      test = X_all[index$test, misspecified_names, drop = FALSE]
    ),
    sensitive = sensitive,
    sensitive_description = sensitive_description,
    sensitive_prevalence = mean(s),
    favorable_outcome_prevalence = mean(y),
    scaling = scaled,
    misspecified_names = misspecified_names,
    source_file = xls_path
  )
}
