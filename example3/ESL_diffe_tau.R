
Huber_ESL_IRLS <- function(beta_init, tau, lambda, X, Y, S_n,
                           eps = 1e-6, max_iter = 100) {
  # 初始值
  beta <- beta_init
  n    <- NROW(X)
  iter <- 1L
  
  while (TRUE) {
    ## 1. 标准化残差 u_i = (Y - X beta) / S_n
    r <- as.vector((Y - X %*% beta) / S_n)
    
    ## 2. 权重 w(u)：
    ##    |u| <= tau 时，w = 1/2；
    ##    |u| >  tau 时，w = (1/2) * exp(-(u^2 - tau^2) / lambda)
    w <- rep(1/2, length(r))          # 默认：|u| <= tau 区间
    idx <- which(abs(r) > tau)
    if (length(idx) > 0L) {
      w[idx] <- exp(-(r[idx]^2 - tau^2) / lambda) / 2
    }
    
    ## 3. 构造加权矩阵 W
    W_mat <- diag(w, nrow = n)
    
    ## 4. 加权最小二乘更新 beta：
    ##    beta_new = (X^T W X)^{-1} X^T W Y
    beta_new <- solve(t(X) %*% W_mat %*% X,
                      t(X) %*% W_mat %*% Y)
    
    ## 5. 收敛判定
    if (norm(beta_new - beta, type = "2") <= eps || iter >= max_iter) {
      return(beta_new)
    }
    
    beta <- beta_new
    iter <- iter + 1L
  }
}














ESL_IRLS <- function(beta_init, tau, X, Y, S_n,
                     eps = 1e-6, max_iter = 100) {
  beta <- beta_init
  n <- NROW(X)
  iter <- 1
  
  while (TRUE) {
    # 标准化残差
    r <- as.vector((Y - X %*% beta) / S_n)
    
    ## 权重 w(u) = psi(u)/(2u) = (1/tau)*exp(-u^2 / tau)
    ## 对 u=0，极限是 1/tau，这里直接用同一公式即可
    w <- (1 / tau) * exp(- (r^2) / tau)
    
    W_mat <- diag(w, nrow = n)
    
    # 加权最小二乘更新 β
    beta_new <- solve(t(X) %*% W_mat %*% X,
                      t(X) %*% W_mat %*% Y)
    
    # 收敛判定
    if (norm(beta_new - beta, type = "2") <= eps || iter >= max_iter) {
      return(beta_new)
    }
    
    beta <- beta_new
    iter <- iter + 1
  }
}



Huber_IRLS <- function(beta_init, tau, X, Y, S_n,
                       eps = 1e-6, max_iter = 100) {
  beta <- beta_init
  n <- NROW(X)
  iter <- 1
  
  while (TRUE) {
    # 标准化残差
    r <- as.vector((Y - X %*% beta) / S_n)
    
    # Huber 权重：ψ(z)/ (2z)
    w <- ifelse(abs(r) <= tau, 
                0.5, 
                tau * sign(r) / (2 * r))
    
    W_mat <- diag(w, nrow = n)
    
    # 加权最小二乘更新 β
    beta_new <- solve(t(X) %*% W_mat %*% X,
                      t(X) %*% W_mat %*% Y)
    
    # 收敛判定
    if (norm(beta_new - beta, type = "2") <= eps || iter >= max_iter) {
      return(beta_new)
    }
    
    beta <- beta_new
    iter <- iter + 1
  }
}


Tukey_IRLS <- function(beta_init, tau, X, Y, S_n,
                       eps = 1e-6, max_iter = 100) {
  beta <- beta_init
  n <- NROW(X)
  iter <- 1
  
  while (TRUE) {
    # 标准化残差
    r <- as.vector((Y - X %*% beta) / S_n)
    
    # Tukey biweight 权重：ψ(z)/(2z)
    w <- ifelse(
      abs(r) <= tau,
      0.5 * (1 - (r / tau)^2)^2,
      0
    )
    
    W_mat <- diag(w, nrow = n)
    
    # 加权最小二乘更新 β
    beta_new <- solve(t(X) %*% W_mat %*% X,
                      t(X) %*% W_mat %*% Y)
    
    # 收敛判定
    if (norm(beta_new - beta, type = "2") <= eps || iter >= max_iter) {
      return(beta_new)
    }
    
    beta <- beta_new
    iter <- iter + 1
  }
}









## ============================================================
## 0) Utilities
## ============================================================

dens_at <- function(u, x0, bw = "nrd0", n_grid = 4096, kernel = "gaussian") {
  d <- density(u, bw = bw, n = n_grid, kernel = kernel)
  approx(d$x, d$y, xout = x0, rule = 2)$y
}

abs_normal_q <- function(c_seg, total_node, alpha = 1.0, power = 3) {
  n <- as.integer(total_node)
  c_seg <- as.integer(c_seg)
  
  if (c_seg < 1L) stop("c_seg must be >= 1")
  if (alpha <= 0 || alpha > 1) stop("alpha must be in (0,1]")
  if (power <= 1) stop("power must be > 1")
  
  n_max <- max(1L, floor(alpha * n))
  if (n_max < c_seg) stop("alpha * total_node is too small: n_max < c_seg")
  
  t <- seq(0, 1, length.out = c_seg)
  i_idx <- 1L + floor((n_max - 1L) * (1 - (1 - t)^power))
  
  p <- 0.5 + i_idx / (2 * (n + 1L))
  q <- qnorm(p)
  
  list(i = i_idx, q = q)
}

## ---- build tau_vec: slopes for segments ----
make_tau_vec <- function(c_seg, total_node, knot_vec) {
  c_seg <- as.integer(c_seg)
  total_node <- as.integer(total_node)
  
  tau_vec <- numeric(c_seg + 1L)
  tau0<-knot_vec[1]
  tau_vec[1L] <- tau0  # tau_0
  
  for (i in 1:c_seg) {
    u <- i / (total_node + 1)
    tau_vec[i + 1L] <- qnorm((u + 1) / 2) + tau0  # tau_i
  }
  tau_vec
}

## ---- build knot_vec of length c+1: (z0, z1, ..., zc) ----
make_knot_vec <- function(c, total_node, alpha2, power) {
  c <- as.integer(c)
  out <- abs_normal_q(c_seg = c+1, total_node = total_node,
                      alpha = alpha2, power = power)
  knot_vec <- out$q                      # length = c
  knot_vec
}


## ============================================================
## 1) psi function for your loss (tail psi = 0)
##    knot_vec = (z0,...,zc), length c+1
## ============================================================


psi_tau_c <- function(z, knot_vec,total_node) {
  c_seg <- length(knot_vec) - 1
  z_abs <- abs(z)
  
  ## 2. 斜率 τ_0,...,τ_c
  tau0 <- knot_vec[1]
  tau_vec <- numeric(c_seg + 1)
  tau_vec[1] <- tau0
  
  for (i in 1:c_seg) {
    u <- i / (total_node + 1)
    tau_vec[i + 1] <- qnorm((u + 1) / 2) + tau0  # τ_i = F_T^+^{-1}(u)+τ_0
  }
  
  psi <- numeric(length(z))
  
  ## (1) |z| <= |z|_[0]: psi(z) = z
  idx0 <- which(z_abs <= knot_vec[1])
  if (length(idx0) > 0) {
    psi[idx0] <- z[idx0]
  }
  
  ## (2) 中间各段: |z|_[i] < |z| <= |z|_[i+1], i = 0,...,c-1
  for (i in 0:(c_seg - 1)) {
    idx <- which(z_abs > knot_vec[i + 1] & z_abs <= knot_vec[i + 2])
    if (length(idx) > 0) {
      psi[idx] <- tau_vec[i + 1] * sign(z[idx])  # τ_i * sign(z)
    }
  }
  
  ## (3) 尾部 |z| > |z|_[c]：rho 为常数 => psi = 0（保持初始 0）
  # 不需要额外赋值
  
  return(psi)
}
## ============================================================
## 2) IRLS for beta (fixed S_n)
##    Standard weight: w = psi(r)/2r
## ============================================================

## X: n x m 矩阵（含截距），每行 x_i = (x_{i1},...,x_{im})
## Y: 长度 n 向量
## br_vec: 折点 |z|_[0],...,|z|_[c]，长度 = c_seg + 1

Modifid_huber_IRLS <- function(beta_init, X, Y, S_n, knot_vec,total_node,
                               eps = 1e-6, max_iter = 100) {
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  beta <- as.numeric(beta_init)
  n <- NROW(X)
  iter <- 1L
  
  while (TRUE) {
    ## 1. 标准化残差 r = (Y - X beta)/S_n
    r <- as.vector((Y - X %*% beta) / S_n)
    
    ## 2. ψ(r)
    psi_val <- psi_tau_c(r, knot_vec,total_node)
    
    ## 3. 权重 w = ψ(r) / r，r≈0 时极限 = 1
    w <- numeric(length(r))
    small <- abs(r) <= knot_vec[1]
    
    w[small]    <- 1/2        # 因为 psi(r) ~ r => psi(r)/r ~ 1
    w[!small]   <- psi_val[!small] / (2*r[!small])
    
    ## 构造对角权重矩阵
    W_mat <- diag(w, nrow = n)
    
    ## 4. 加权最小二乘更新 β
    XtWX <- t(X) %*% W_mat %*% X
    XtWy <- t(X) %*% W_mat %*% Y
    
    beta_new <- solve(XtWX, XtWy)
    
    ## 5. 收敛判定
    if (norm(beta_new - beta, type = "2") <= eps || iter >= max_iter) {
      return(beta_new)
    }
    
    beta <- beta_new
    iter <- iter + 1L
  }
}

## ============================================================
## 3) Compute Jhat = Bhat / Ahat^2 (with jump-corrected Ahat)
## ============================================================

compute_Jhat <- function(X, Y, beta_init, S_n,
                         knot_vec, total_node,
                         bw = "nrd0",
                         eps = 1e-6, max_iter = 200) {
  X <- as.matrix(X); Y <- as.numeric(Y)
  
  c_seg <- length(knot_vec) - 1L
  tau_vec <- make_tau_vec(c_seg = c_seg, total_node = total_node, knot_vec = knot_vec)
  
  beta_hat <- Modifid_huber_IRLS(
    beta_init = beta_init, X = X, Y = Y, S_n = S_n,
    knot_vec = knot_vec, total_node = total_node,
    eps = eps, max_iter = max_iter
  )
  
  u <- (Y - as.numeric(X %*% beta_hat)) / S_n  # standardized residuals
  z0 <- knot_vec[1L]
  zc <- knot_vec[c_seg + 1L]
  
  ## Bhat = E[psi^2] empirical
  psi_u <- psi_tau_c(u, knot_vec, total_node)
  Bhat <- mean(psi_u^2)
  
  ## Ahat = integral phi dpsi estimated by: A0 + jumps
  A0 <- mean(abs(u) < z0)  # naive part (E[psi'(U)] within core)
  
  ## slope_i = tau_i used in middle segments: i=0,...,c-1
  slope <- tau_vec[1:c_seg]  # tau_0,...,tau_{c-1}
  
  ## internal jumps at z_1,...,z_{c-1}
  if (c_seg >= 2L) {
    z_in <- knot_vec[2:c_seg]                 # z_1,...,z_{c-1}
    f_in <- dens_at(u, z_in, bw = bw)
    delta_slope <- slope[2:c_seg] - slope[1:(c_seg - 1L)]
    jump_in <- 2 * sum(delta_slope * f_in)
  } else {
    jump_in <- 0
  }
  
  ## tail jump at z_c: psi drops from tau_{c-1}*sign to 0
  f_zc <- dens_at(u, zc, bw = bw)
  jump_tail <- -2 * slope[c_seg] * f_zc
  
  Ahat <- A0 + jump_in + jump_tail
  Jhat <- Bhat / (Ahat^2)
  
  list( u = u, Ahat = Ahat, Bhat = Bhat, Jhat = Jhat)
}
## ============================================================
## 4) Grid search over tau_1 and (tau_1,lambda)
## ============================================================



make_folds <- function(n, K = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx <- sample.int(n)
  split(idx, cut(seq_along(idx), breaks = K, labels = FALSE))
}





tune_esl_cv <- function(new_X, Y, LMS_beta, Sn,
                        tau_grid = c(3, 6, 10, 15, 20),
                        K = 5, seed = NULL,
                        eps = 1e-3, max_iter = 500) {
  
  folds <- make_folds(n = nrow(new_X), K = K, seed = seed)
  
  fold_scores <- matrix(NA_real_, nrow = length(tau_grid), ncol = K,
                        dimnames = list(paste0("tau=", tau_grid), paste0("fold", 1:K)))
  
  cv_vals <- numeric(length(tau_grid))
  
  for (j in seq_along(tau_grid)) {
    tau <- tau_grid[j]
    fs  <- numeric(K)
    
    for (k in seq_len(K)) {
      val_idx <- folds[[k]]
      tr_idx  <- setdiff(seq_len(nrow(new_X)), val_idx)
      
      Xtr <- new_X[tr_idx, , drop = FALSE]
      Ytr <- Y[tr_idx]
      Xva <- new_X[val_idx, , drop = FALSE]
      Yva <- Y[val_idx]
      
      beta_hat <- ESL_IRLS(LMS_beta, tau = tau,
                           X = Xtr, Y = Ytr,
                           S_n = Sn, eps = eps,
                           max_iter = max_iter)
      
      pred <- as.vector(Xva %*% beta_hat)
      fs[k] <- median(abs(Yva - pred))
      fold_scores[j, k] <- fs[k]
    }
    
    cv_vals[j] <- median(fs)  # 聚合：K 折得分的中位数
  }
  
  best_j   <- which.min(cv_vals)
  best_tau <- tau_grid[best_j]
  
  beta_best <- ESL_IRLS(LMS_beta, tau = best_tau,
                        X = new_X, Y = Y,
                        S_n = Sn, eps = eps,
                        max_iter = max_iter)
  
  list(
    best_tau = best_tau,
    cv_table = data.frame(tau = tau_grid, cv_score = cv_vals),
    fold_scores = fold_scores,
    beta_hat = beta_best
  )
}

tune_hesl_cv <- function(new_X, Y, LMS_beta, Sn,
                         tau_grid = c(3, 6, 10, 15, 20),
                         lambda_grid = c(0.2, 0.5, 1.0),
                         K = 5, seed = NULL,
                         eps = 1e-3, max_iter = 500) {
  folds <- make_folds(n = nrow(new_X), K = K, seed = seed)
  
  grid <- expand.grid(tau = tau_grid, lambda = lambda_grid, KEEP.OUT.ATTRS = FALSE)
  cv_vals <- numeric(nrow(grid))
  
  cv_score_hesl <- function(tau, lambda) {
    fs <- numeric(K)
    for (k in seq_len(K)) {
      val_idx <- folds[[k]]
      tr_idx  <- setdiff(seq_len(nrow(new_X)), val_idx)
      
      Xtr <- new_X[tr_idx, , drop = FALSE]
      Ytr <- Y[tr_idx]
      Xva <- new_X[val_idx, , drop = FALSE]
      Yva <- Y[val_idx]
      
      beta_hat <- Huber_ESL_IRLS(LMS_beta, tau = tau, lambda = lambda,
                                 X = Xtr, Y = Ytr,
                                 S_n = Sn, eps = eps,
                                 max_iter = max_iter)
      
      pred <- as.vector(Xva %*% beta_hat)
      fs[k] <- median(abs(Yva - pred))
    }
    median(fs)
  }
  
  for (g in seq_len(nrow(grid))) {
    cv_vals[g] <- cv_score_hesl(grid$tau[g], grid$lambda[g])
  }
  
  best_g <- which.min(cv_vals)
  best_tau <- grid$tau[best_g]
  best_lambda <- grid$lambda[best_g]
  
  beta_best <- Huber_ESL_IRLS(LMS_beta, tau = best_tau, lambda = best_lambda,
                              X = new_X, Y = Y,
                              S_n = Sn, eps = eps,
                              max_iter = max_iter)
  
  out_tbl <- cbind(grid, cv_score = cv_vals)
  out_tbl <- out_tbl[order(out_tbl$cv_score, out_tbl$tau, out_tbl$lambda), ]
  
  list(
    best_tau = best_tau,
    best_lambda = best_lambda,
    cv_table = out_tbl,
    beta_hat = beta_best
  )
}

#0.0235


library(robustbase)
library(MASS)
library(VGAM)
library(L1pack)
library(adaHuber)
library(quantreg)

n <- 500
p <- 10

# (beta0, beta1, ..., betap) = (1,2,3,3,0,...,0)
beta <- c(1, 2, 3, 3, rep(0, p - 3))


n_rep <- 1000


## 各种方法的 MSE / bias 容器（可以以后改成矩阵）
diff_tau_esl_MSE_lis=c()
diff_tau_esl_bias_lis=c()


tau_parmater_lis=c(1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3)
for (tau_parmater in  tau_parmater_lis) {
  esl_MSE_lis=c()
  esl_bias_lis=c() 
  for (i in 1:1000) {
    set.seed(i)
    X <- matrix(rnorm(n * p), n, p)
    new_X <- cbind(1, X)
    err <- rnorm(n,0,1)
    
    Y <- new_X %*% beta + err
    
    outlier_frac <- 0.05
    m <- floor(outlier_frac * n)
    idx_out <- sample.int(n, m)
    Y[idx_out]<-rnorm(m, mean = 1, sd = 1)
    
    ## --- LMS ---
    LMS_estimator <- lqs(X, Y, intercept = TRUE, method = "lms")
    LMS_beta <- LMS_estimator$coefficients
    err_ <- Y - new_X %*% LMS_beta
    med_r   <- median(err_, na.rm = TRUE)
    med_dev <- median(abs(err_ - med_r), na.rm = TRUE)
    Sn <- 1.4826 * med_dev
    
    
    
    
    
    
    ESL_beta <- ESL_IRLS(LMS_beta, tau = tau_parmater,
                         X = new_X, Y = Y,
                         S_n = Sn, eps = 1e-3,
                         max_iter = 500)
    
    ESL_MSE  <- norm(ESL_beta - beta, "2")^2
    ESL_bias <- sum(abs(ESL_beta - beta))
    esl_MSE_lis  <- c(esl_MSE_lis, ESL_MSE)
    esl_bias_lis <- c(esl_bias_lis, ESL_bias)
    
    
    
    if (i %% 50 == 0) cat("i =", i, "\n")
  }
  print(mean(esl_MSE_lis))
  diff_tau_esl_MSE_lis=c(diff_tau_esl_MSE_lis,mean(esl_MSE_lis))
  diff_tau_esl_bias_lis=c(diff_tau_esl_bias_lis,mean(esl_bias_lis))
  
}

diff_tau_esl_MSE_lis
diff_tau_esl_bias_lis
df <- data.frame(
  tau = tau_parmater_lis,
  mse   = diff_tau_esl_MSE_lis,
  bias  = diff_tau_esl_bias_lis
)

write.csv(df, "C:/Users/apple/Desktop/modified_huber/diff_tau_esl.csv", row.names = FALSE)








