competitive_models_cv <- function(data, max_cv = 10, seed = 12345) {

  # ============================================================
  # Storage objects (same structure as your script)
  # ============================================================
  test_pred_ordASDA <- list()

  ordinal_logistic_class_CV <- list()
  ordinal_logistic_prob_CV  <- list()

  multi_logistic_class_CV <- list()
  multi_logistic_prob_CV  <- list()

  rf_class_CV <- list()
  rf_prob_CV  <- list()

  ordfor_class_CV <- list()

  # ============================================================
  # Cross-validation loop
  # ============================================================
  for (cv in 1:max_cv) {

    cat("CV =", cv, "\n")

    # ----------------------------------------------------------
    # Train/test split (20% test)
    # ----------------------------------------------------------
    set.seed(seed + cv)
    kk <- sample(1:nrow(data), size = round(0.2 * nrow(data)))

    train_data <- data[-kk, ]
    test_data  <- data[kk, ]

    X <- as.matrix(train_data[, -1])
    Y <- train_data[, 1]


    # ==========================================================
    # 1. Ordinal Logistic Regression (Proportional Odds)
    # ==========================================================
    set.seed(seed + cv)

    fit_polr <- MASS::polr(
      Y ~ .,
      data = train_data,
      method = "logistic"
    )

    # Class predictions
    pred_class <- predict(fit_polr, newdata = test_data[, -1], type = "class")

    # Probability predictions
    pred_probs <- predict(fit_polr, newdata = test_data[, -1], type = "probs")

    ordinal_logistic_class_CV[[cv]] <- pred_class
    ordinal_logistic_prob_CV[[cv]]  <- pred_probs


    # ==========================================================
    # 2. Multinomial Logistic Regression (ignores ordering)
    # ==========================================================
    set.seed(seed + cv)

    fit_multinom <- nnet::multinom(
      Y ~ .,
      data = train_data,
      trace = FALSE
    )

    pred_probs_nom <- predict(fit_multinom,
                              newdata = test_data[, -1],
                              type = "probs")

    pred_class_nom <- predict(fit_multinom,
                              newdata = test_data[, -1],
                              type = "class")

    multi_logistic_class_CV[[cv]] <- pred_class_nom
    multi_logistic_prob_CV[[cv]]  <- pred_probs_nom


    # ==========================================================
    # 3. Random Forest (nonparametric baseline)
    # ==========================================================
    set.seed(seed + cv)

    fit_rf <- randomForest::randomForest(
      Y ~ .,
      data = train_data,
      ntree = 500,
      mtry = floor(sqrt(ncol(X)))
    )

    pred_class_rf <- predict(fit_rf,
                             newdata = test_data[, -1],
                             type = "class")

    pred_probs_rf <- predict(fit_rf,
                             newdata = test_data[, -1],
                             type = "prob")

    rf_class_CV[[cv]] <- pred_class_rf
    rf_prob_CV[[cv]]  <- pred_probs_rf


    # ==========================================================
    # 4. Ordinal Forest (order-aware tree ensemble)
    # ==========================================================
    set.seed(seed + cv)

    fit_ordfor <- ordinalForest::ordfor(
      depvar = "Y",
      data = train_data,
      nsets = 1000,
      ntreeperdiv = 100,
      ntreefinal = 500,
      nbest = 10
    )

    pred_ordfor <- predict(fit_ordfor, newdata = test_data[, -1])

    # Ordered class predictions
    pred_class_ordfor <- pred_ordfor$ypred

    ordfor_class_CV[[cv]] <- pred_class_ordfor
  }

  # ============================================================
  # Return all results
  # ============================================================
  return(list(
    ordASDA_pred      = test_pred_ordASDA,
    ordinal_logistic  = list(class = ordinal_logistic_class_CV,
                             prob  = ordinal_logistic_prob_CV),
    multinomial       = list(class = multi_logistic_class_CV,
                             prob  = multi_logistic_prob_CV),
    random_forest     = list(class = rf_class_CV,
                             prob  = rf_prob_CV),
    ordinal_forest    = list(class = ordfor_class_CV)
  ))
}
