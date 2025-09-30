# ------------------------------------------------------------
# c-TREX (exact) via 2p SOCPs in CVXR  —  requires CVXR >= 1.0
# ------------------------------------------------------------
# Objective:
#  minimize   ||y - X beta||_2^2 / (v^T X^T (y - X beta))  +  c * ||beta||_1
#  subject to v^T X^T (y - X beta) >= ||X^T (y - X beta)||_inf
# for each v in {±e_j}, j = 1..p; pick the best among 2p problems.

ctrex <- function(X, y, c = 0.5,
                  center = TRUE, scale = TRUE,
                  solver = c("ECOS","SCS"), solver_opts = list(),
                  verbose = FALSE) {
  if (!requireNamespace("CVXR", quietly = TRUE))
    stop("Please install the 'CVXR' package.")
  solver <- match.arg(solver)

  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  if (length(y) != n) stop("length(y) must equal nrow(X).")

  # standardize predictors and center y (like in q-TREX helper I gave earlier)
  x_means <- if (center) colMeans(X) else rep(0, p)
  Xc <- sweep(X, 2, x_means, "-")
  x_sds <- if (scale) apply(Xc, 2, function(v) sd(v) + 1e-12) else rep(1, p)
  Xs <- sweep(Xc, 2, x_sds, "/")
  y_mean <- if (center) mean(y) else 0
  ys <- y - y_mean

  # CVXR constants & helper
  Xc_cvx <- CVXR::as_matrix(Xs)
  y_cvx  <- CVXR::as_vector(ys)

  best_obj <- Inf
  best_beta <- rep(0, p)
  best_tag <- NULL
  sub_objs  <- matrix(NA_real_, nrow = p, ncol = 2,
                      dimnames = list(colnames(X), c("+","-")))

  # loop over j and s in {+1, -1}
  for (j in seq_len(p)) {
    for (sgn in c(1,-1)) {
      v <- numeric(p); v[j] <- sgn

      beta <- CVXR::Variable(p)
      r    <- y_cvx - Xc_cvx %*% beta
      s    <- t(Xc_cvx) %*% r

      denom <- CVXR::t(v) %*% s
      # c-TREX subproblem:
      # minimize quad_over_lin(r, denom) + c * ||beta||_1
      # subject to denom >= max(abs(s))  (enforces index/sign)
      obj <- CVXR::quad_over_lin(r, denom) + c * CVXR::norm1(beta)
      constr <- list( denom >= CVXR::max_elemwise(CVXR::abs(s)) )

      prob <- CVXR::Problem(CVXR::Minimize(obj), constr)

      result <- CVXR::solve(prob,
                            solver = solver,
                            verbose = verbose,
                            feastol = solver_opts$feastol,
                            reltol = solver_opts$reltol,
                            abstol = solver_opts$abstol,
                            max_iters = solver_opts$max_iters)

      sub_objs[j, ifelse(sgn>0,1,2)] <- result$value
      if (is.finite(result$value) && result$value < best_obj) {
        best_obj  <- result$value
        best_beta <- as.numeric(CVXR::result(beta))
        best_tag  <- list(j = j, sign = sgn)
      }
    }
  }

  # map coefficients back to original scale
  beta_orig <- best_beta / x_sds
  intercept <- if (center) y_mean - sum(x_means * beta_orig) else 0
  names(beta_orig) <- colnames(X)

  structure(list(
    coefficients = c(`(Intercept)` = intercept, beta_orig),
    beta_standardized = best_beta,
    objective = best_obj,
    chosen_index = best_tag,
    subproblem_objectives = sub_objs,
    center = center, scale = scale,
    x_means = x_means, x_sds = x_sds, y_mean = y_mean,
    call = match.call()
  ), class = "trex_fit")
}

# ------------------------------------------------------------
# B-TREX  (bootstrap TREX selection)
# ------------------------------------------------------------
# Repeats TREX (by default: exact c-TREX above; you can switch to "qtrex")
# on bootstrap resamples of rows, tallies selection frequencies, and
# returns a ranked list + an optional LS refit on the selected set.
#
# Selection rule:
#   - "majority": select vars with freq > 0.5
#   - "threshold": select vars with freq >= pi_thr (e.g., 0.6)
#   - "topk": select top-K by frequency
#
# Set refit = TRUE to least-squares refit coefficients on full data.

btrex <- function(X, y, B = 31,
                  engine = c("ctrex","qtrex"),
                  c = 0.5,
                  selection = c("majority","threshold","topk"),
                  pi_thr = 0.6, k = 10,
                  refit = TRUE,
                  seed = NULL,
                  ... # passed to ctrex() or trex()
                  ) {
  engine <- match.arg(engine)
  selection <- match.arg(selection)
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  if (!is.null(seed)) set.seed(seed)

  sel_counts <- numeric(p)
  coefs_list <- vector("list", B)

  for (b in seq_len(B)) {
    idx <- sample.int(n, replace = TRUE)
    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]

    fit <- switch(engine,
      ctrex = do.call(ctrex, c(list(X = Xb, y = yb, c = c), list(...))),
      qtrex = do.call(trex,  c(list(X = Xb, y = yb, c = c), list(...)))
    )
    beta <- fit$coefficients[-1]  # drop intercept
    sel <- as.integer(abs(beta) > 0)
    sel_counts <- sel_counts + sel
    coefs_list[[b]] <- beta
  }

  freq <- sel_counts / B
  names(freq) <- colnames(X)

  selected <- switch(selection,
    majority  = which(freq > 0.5),
    threshold = which(freq >= pi_thr),
    topk      = {
      ord <- order(freq, decreasing = TRUE)
      ord[seq_len(min(k, length(ord)))]
    }
  )

  beta_refit <- rep(0, p); names(beta_refit) <- colnames(X)
  intercept  <- 0
  if (refit && length(selected) > 0) {
    Xs <- X[, selected, drop = FALSE]
    # Least squares with intercept
    df <- data.frame(y = y, Xs)
    fit_ls <- lm(y ~ . , data = df)
    cf <- coef(fit_ls)
    intercept <- unname(cf[1])
    beta_refit[selected] <- unname(cf[-1])
  }

  # create a tidy ranking table
  ranking <- data.frame(
    feature = colnames(X),
    frequency = as.numeric(freq)
  )
  ranking <- ranking[order(ranking$frequency, decreasing = TRUE), , drop = FALSE]
  rownames(ranking) <- NULL

  list(
    frequencies = freq,
    selected_idx = selected,
    selected_names = colnames(X)[selected],
    refit = list(intercept = intercept, beta = beta_refit),
    ranking = ranking,
    bootstrap_betas = coefs_list,
    call = match.call()
  )
}

# ---------------------------------------------
# Minimal usage example (synthetic)
# ---------------------------------------------
# set.seed(1)
# n <- 120; p <- 50
# X <- matrix(rnorm(n*p), n, p)
# beta_true <- c(rep(2, 5), rep(0, p-5))
# y <- X %*% beta_true + rnorm(n)
#
# # Exact c-TREX
# fit_c <- ctrex(X, y, solver = "ECOS", solver_opts = list(reltol = 1e-7))
# fit_c$coefficients[1:10]
#
# # B-TREX using c-TREX engines + majority vote + LS refit
# bt <- btrex(X, y, B = 31, engine = "ctrex",
#             selection = "majority", refit = TRUE,
#             solver = "SCS", solver_opts = list(reltol = 1e-5))
# bt$selected_names
# head(bt$ranking, 10)
