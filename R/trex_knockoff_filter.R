# --------------------------------------------
# TREX knockoff importance statistic (signed-max)
# --------------------------------------------
# X, X_k: n-by-p originals and knockoffs (same column order)
# y: response
# engine: "ctrex" (exact SOCP via CVXR) or "qtrex" (prox-gradient)
# ...: forwarded to ctrex()/trex() (e.g., solver options)

trex_knock_stat <- function(X, X_k, y,
                            engine = c("ctrex","qtrex"),
                            c = 0.5, ...) {
  engine <- match.arg(engine)
  X_aug <- cbind(X, X_k)

  # Fit TREX once on augmented design
  fit <- switch(engine,
    ctrex = ctrex(X_aug, y, c = c, ...),
    qtrex = trex( X_aug, y, c = c, ...)
  )
  beta_aug <- fit$coefficients[-1]  # drop intercept
  p <- ncol(X)
  beta_orig <- beta_aug[1:p]
  beta_k    <- beta_aug[(p+1):(2*p)]

  Z  <- abs(beta_orig)
  Zk <- abs(beta_k)

  # Signed-max knockoff statistic:
  # W_j = max(Z_j, Z~_j) * sign(Z_j - Z~_j)
  W <- pmax(Z, Zk) * sign(Z - Zk)
  as.numeric(W)
}

# --------------------------------------------
# Full knockoff+ TREX wrapper
# --------------------------------------------
# knockoffs: function to create knockoffs (default: second-order Gaussian)
# fdr: target FDR; offset=1 yields knockoff+ (recommended)
# engine: "ctrex" or "qtrex"
# Returns the usual knockoff result plus the fitted W and beta split.

trex_knockoff_filter <- function(X, y,
                                 knockoffs = knockoff::create.second_order,
                                 fdr = 0.1, offset = 1,
                                 engine = c("ctrex","qtrex"),
                                 c = 0.5,
                                 ...) {
  engine <- match.arg(engine)

  # define the statistic function with engine/c bound in
  k_stat <- function(X, X_k, y) trex_knock_stat(X, X_k, y, engine = engine, c = c, ...)

  result <- knockoff::knockoff.filter(
    X, y,
    knockoffs = knockoffs,
    statistic = k_stat,
    fdr = fdr,
    offset = offset
  )
  result
}

# --------------------------------------------
# Example (synthetic; model-X knockoffs)
# --------------------------------------------
# set.seed(123)
# n <- 120; p <- 60
# rho <- 0.4
# Sigma <- toeplitz(rho^(0:(p-1)))
# X <- matrix(rnorm(n*p), n) %*% chol(Sigma)
# beta_true <- c(rep(2, 8), rep(0, p-8))
# y <- as.numeric(X %*% beta_true + rnorm(n))
#
# # Model-X knockoffs with known (mu, Sigma)
# mu <- rep(0, p)
# Kfun <- function(X) knockoff::create.gaussian(X, mu, Sigma, method="asdp")
#
# # Exact c-TREX engine
# res_ctrex <- trex_knockoff_filter(X, y, knockoffs = Kfun,
#                                   engine = "ctrex", fdr = 0.1,
#                                   solver = "ECOS", solver_opts = list(reltol=1e-7))
# res_ctrex$selected
#
# # Fast q-TREX engine
# res_qtrex <- trex_knockoff_filter(X, y, knockoffs = Kfun,
#                                   engine = "qtrex", fdr = 0.1,
#                                   maxit = 3000, tol = 1e-6)
# res_qtrex$selected
