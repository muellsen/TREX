# -----------------------------
# TREX (proximal-gradient solver)
# -----------------------------

trex <- function(X, y, c = 0.5,
                 center = TRUE, scale = TRUE,
                 maxit = 5000, tol = 1e-6,
                 step_init = 1e-2, step_min = 1e-8,
                 backtrack = 0.5, verbose = FALSE,
                 beta0 = NULL) {
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  if (length(y) != n) stop("length(y) must equal nrow(X).")

  # standardize (optional)
  x_means <- if (center) colMeans(X) else rep(0, p)
  Xc <- sweep(X, 2, x_means, FUN = "-")
  x_sds <- if (scale) apply(Xc, 2, function(v) sd(v) + 1e-12) else rep(1, p)
  Xs <- sweep(Xc, 2, x_sds, FUN = "/")
  y_mean <- if (center) mean(y) else 0
  ys <- y - y_mean

  # helpers --------------------------------------------------------
  soft_thresh <- function(z, t) sign(z) * pmax(abs(z) - t, 0)

  obj_val <- function(beta) {
    r <- ys - Xs %*% beta
    s <- crossprod(Xs, r)               # p-vector
    denom <- max(abs(s))
    if (denom <= 0) return(Inf)
    num <- sum(r * r)
    num / denom + c * sum(abs(beta))
  }

  grad_g <- function(beta) {
    # subgradient of g(beta) = ||r||^2 / ||X^T r||_inf
    r <- as.numeric(ys - Xs %*% beta)
    s <- as.numeric(crossprod(Xs, r))
    a <- max(abs(s))
    jstar <- which.max(abs(s))          # pick any maximizer
    tsgn <- sign(s[jstar]); if (tsgn == 0) tsgn <- 1
    h <- sum(r * r)                      # ||r||^2
    # ∇h = -2 X^T r
    grad_h <- -2 * crossprod(Xs, r)     # p-vector
    # subgradient of a = ||X^T r||_inf at j*:
    # ∂a/∂beta = - tsgn * X^T x_{j*}
    xj <- Xs[, jstar, drop = FALSE]
    grad_a <- - tsgn * crossprod(Xs, xj)  # p-vector
    # quotient rule: (∇h * a - h * ∇a) / a^2
    as.numeric( (grad_h * a - h * grad_a) / (a^2) )
  }

  # init
  beta <- if (is.null(beta0)) rep(0, p) else as.numeric(beta0)
  step <- step_init
  f_curr <- obj_val(beta)

  for (iter in 1:maxit) {
    g <- grad_g(beta)
    # backtracking line search on composite objective
    accepted <- FALSE
    while (!accepted && step >= step_min) {
      beta_proposed <- soft_thresh(beta - step * g, c * step)
      f_prop <- obj_val(beta_proposed)
      if (is.finite(f_prop) && f_prop <= f_curr + 1e-12) {
        accepted <- TRUE
      } else {
        step <- step * backtrack
      }
    }
    if (!accepted) break

    # update & check convergence
    if (max(abs(beta_proposed - beta)) < tol) {
      beta <- beta_proposed; f_curr <- f_prop
      if (verbose) message(sprintf("Converged in %d iters; obj=%.6g", iter, f_curr))
      break
    }
    beta <- beta_proposed
    f_curr <- f_prop
    if (verbose && iter %% 100 == 0) {
      message(sprintf("iter %d: obj=%.6g, step=%.2e", iter, f_curr, step))
    }
  }

  # map back to original scale
  beta_orig <- beta / x_sds
  # intercept on original scale (if centered)
  intercept <- if (center) y_mean - sum(x_means * beta_orig) else 0

  structure(list(
    coefficients = c(`(Intercept)` = intercept, beta_orig),
    beta_standardized = beta,
    center = center, scale = scale,
    x_means = x_means, x_sds = x_sds, y_mean = y_mean,
    objective = f_curr, iterations = iter, step = step,
    call = match.call()
  ), class = "trex_fit")
}

predict.trex_fit <- function(object, newx, ...) {
  b <- object$coefficients
  intercept <- b[1]; beta <- b[-1]
  as.numeric(intercept + as.matrix(newx) %*% beta)
}

# -----------------------------
# Example (small synthetic)
# -----------------------------
# set.seed(1)
# n <- 100; p <- 50
# X <- matrix(rnorm(n*p), n, p)
# beta_true <- c(rep(2, 5), rep(0, p-5))
# y <- X %*% beta_true + rnorm(n, sd = 1)
# fit <- trex(X, y, verbose = TRUE)
# fit$coefficients
