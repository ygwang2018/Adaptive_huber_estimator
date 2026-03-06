








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



library(adaHuber)
library(robustbase)
library(MASS)
library(VGAM)
library(L1pack)





dataset1<-read.csv("C:\\Users\\PC\\Desktop\\Polonomial_block_descen_algorithm\\application\\plasma_data.csv",header=F)
y1<-dataset1$V11
x1<-cbind(dataset1$V5,dataset1$V10)



B <- 100     # bootstrap 次数
n <- NROW(x1)

bootstrap_list <- vector("list", B)



for (b in 1:B) {
  set.seed(b)
  # Step 1: bootstrap 抽样 n 个样本（有放回）
  boot_idx <- sample.int(n, size = n, replace = TRUE)
  
  # Step 2: 在 bootstrap 样本内部按“位置”划分训练/测试（不去重）
  pos_train <- sample.int(n, size = 250, replace = FALSE)  # 抽 200 个位置
  train_idx_in_boot <- boot_idx[pos_train]
  test_idx_in_boot  <- boot_idx[-pos_train]  # 剩余位置，保留重复
  
  bootstrap_list[[b]] <- list(
    X_train = x1[train_idx_in_boot, ],
    y_train = y1[train_idx_in_boot],
    X_test  = x1[test_idx_in_boot, ],
    y_test  = y1[test_idx_in_boot]
  )
}


## ---- 每个估计量的 MAPE 向量 ----
LSE_MAPE      <- numeric(B)
LTS_MAPE      <- numeric(B)
Huber_MAPE    <- numeric(B)
MM_MAPE       <- numeric(B)
LMS_MAPE      <- numeric(B) 
Tukey_MAPE    <- numeric(B)
S_MAPE        <- numeric(B)
Proposed_MAPE <- numeric(B)
Huber_ESL_MAPE<- numeric(B)
ESL_MAPE<- numeric(B)
opti_c_lis     <- numeric(B)
opti_alph_lis  <- numeric(B)
optial_power_lis<- numeric(B)
optial_c_total_lis<- numeric(B)
## ---- 存储每折 β 向量 ----
LSE_beta_list      <- list()
LTS_beta_list      <- list()
LMS_beta_list      <- list() 
Huber_beta_list    <- list()
MM_beta_list       <- list()
Tukey_beta_list    <- list()
S_beta_list        <- list()
Proposed_beta_list <- list()
Huber_ESL_beta_list<- list()
ESL_beta_list<-list()



for (i in 1:B) {
  
  X_train <- bootstrap_list[[i]]$X_train
  y_train <- bootstrap_list[[i]]$y_train
  X_test  <- bootstrap_list[[i]]$X_test
  y_test  <- bootstrap_list[[i]]$y_test
  
  use_n <- NROW(X_train)
  use_p <- NCOL(X_train)
  
  ## 定义一个 data.frame 方便 formula 接口使用
  df_train <- data.frame(y = y_train, X_train)
  colnames(df_train) <- c("y", paste0("x", 1:use_p))
  df_test  <- data.frame(y = y_test, X_test)
  colnames(df_test) <- c("y", paste0("x", 1:use_p))
  ## ------------ LSE ------------
  beta_LSE <- solve(t(X_train) %*% X_train) %*% t(X_train) %*% y_train
  LSE_MAPE[i] <- median(abs(y_test - X_test %*% beta_LSE))
  LSE_beta_list[[i]] <- as.vector(beta_LSE)
  
  ## ------------ LTS ------------
  set.seed(5)
  lts_model <- ltsReg(y ~ . - 1, data = df_train)
  LTS_beta  <- coef(lts_model)
  LTS_MAPE[i] <- median(abs(y_test - as.matrix(X_test) %*% LTS_beta))
  LTS_beta_list[[i]] <- as.vector(LTS_beta)
  ## ------------ LMS ------------
  LMS_est  <- lqs(X_train, y_train,
                  intercept = FALSE,
                  method = "lms")
  LMS_beta <- LMS_est$coefficients
  LMS_MAPE[i] <- median(abs(y_test - X_test %*% LMS_beta))
  LMS_beta_list[[i]] <- as.vector(LMS_beta)
  ## ------------ Compute Sn (共用) ------------
  err_ <- y_train - X_train %*% LTS_beta
  med_r  <- median(err_)
  med_dev <- median(abs(err_ - med_r))
  Sn <- 1.4826 * med_dev
  y_mat <- matrix(y_train, ncol = 1)
  ## --- ESL ---
  
  
 
  tau_hat <- tune_esl_cv(X_train, y_train, LMS_beta, Sn,
                                     tau_grid = c(3, 6, 10, 15, 20),
                                     K = 5, seed = 1,
                                     eps = 1e-3, max_iter = 500)$best_tau 
  #print(tau_hat)
  # 用选出来的 tau 在全训练集上重拟合
  ESL_beta <- ESL_IRLS(LMS_beta, tau = tau_hat,
                       X = X_train, Y = y_train,
                       S_n = Sn, eps = 1e-3,
                       max_iter = 500)
  
  ESL_MAPE[i] <- median(abs(y_test - X_test %*%  ESL_beta))
  ESL_beta_list[[i]] <- as.vector(ESL_beta)
  
  ## --- Huber_ESL ---
  
  best_lambda_tau<-tune_hesl_cv(X_train, y_train, LMS_beta, Sn,
                                tau_grid = c(3, 6, 10, 15, 20),
                                lambda_grid = c(0.2, 0.5, 1.0),
                                K = 5, seed = 1,
                                eps = 1e-3, max_iter = 500)
  tau_hat2 <- best_lambda_tau$best_tau
  lambda_hat <-  best_lambda_tau$best_lambda 
  # 用选出来的 (tau,lambda) 在全训练集上重拟合
  Huber_ESL_beta <- Huber_ESL_IRLS(LMS_beta, tau = tau_hat2, lambda = lambda_hat,
                                   X = X_train, Y =y_train,
                                   S_n = Sn, eps = 1e-3,
                                   max_iter = 500)
  
  Huber_ESL_MAPE[i] <- median(abs(y_test - X_test %*%  Huber_ESL_beta))
  Huber_ESL_beta_list[[i]] <- as.vector(Huber_ESL_beta)
  
  
  ## ------------ Huber ------------
  huber_beta <- Huber_IRLS(
    beta_init = LMS_beta, tau = 1.345,
    X = X_train, Y = y_mat,
    S_n = Sn, eps = 1e-3,
    max_iter = 500
  )
  Huber_MAPE[i] <- median(abs(y_test - X_test %*% huber_beta))
  Huber_beta_list[[i]] <- as.vector(huber_beta)
  
  ## ------------ MM ------------
  MM_model <- rlm(y ~ . - 1, data = df_train, method = "MM")
  MM_beta  <- as.matrix(MM_model$coefficients)
  MM_MAPE[i] <- median(abs(y_test - X_test %*% MM_beta))
  MM_beta_list[[i]] <- as.vector(MM_beta)
  
  ## ------------ Tukey ------------
  tukey_beta <- Tukey_IRLS(
    beta_init = LMS_beta, tau = 4.685,
    X = X_train, Y = y_mat,
    S_n = Sn, eps = 1e-3,
    max_iter = 500
  )
  Tukey_MAPE[i] <- median(abs(y_test - X_test %*% tukey_beta))
  Tukey_beta_list[[i]] <- as.vector(tukey_beta)
  
  ## ------------ S estimator ------------
  m.lmS <- lmrob.S(x = X_train, y = y_train,
                   control = lmrob.control(nRes = 20),
                   trace.lev = 0)
  S_beta <- as.matrix(m.lmS$coefficients)
  S_MAPE[i] <- median(abs(y_test - X_test %*% S_beta))
  S_beta_list[[i]] <- as.vector(S_beta)
  
  ## -----modified huber-------
  
  grid_tbl<- rbind(
    data.frame(alpha = 1, c = floor(50/10),total_node=50,power=1.1),
    data.frame(alpha = 1, c = floor((2*50)/10),total_node=50,power=1.1),
    data.frame(alpha = 1, c = floor((3*50)/10),total_node=50,power=1.1),
    data.frame(alpha = 1, c = floor((4*50)/10),total_node=50,power=1.1),
    data.frame(alpha = 1, c = floor((5*50)/10),total_node=50,power=1.1),
    data.frame(alpha = 1, c = floor((6*50)/10),total_node=50,power=1.1),
    data.frame(alpha = 1, c = floor((7*50)/10),total_node=50,power=1.1),
    data.frame(alpha = 1, c = floor((8*50)/10),total_node=50,power=1.1),
    data.frame(alpha = 1, c = floor((9*50)/10),total_node=50,power=1.1),
    data.frame(alpha = 1, c =29999,total_node=30000,power=1.1)
  )
  
  
  
  for (k in seq_len(nrow(grid_tbl))) {
    
    a2 <- grid_tbl$alpha[k]
    c  <- as.integer(grid_tbl$c[k])
    power<-grid_tbl$power[k]
    total_node<-grid_tbl$total_node[k]
    knot_vec <- make_knot_vec(c = c, total_node = total_node,
                              alpha2 = a2, power = power)
    res <- compute_Jhat(
      X =  X_train, Y = y_train, beta_init = LMS_beta, S_n = Sn,
      knot_vec = knot_vec, total_node = total_node,
      bw = 1, eps = 10^-3, max_iter = 500
    )
    
    grid_tbl$Jhat[k] <- res$Jhat
    
  }
  
  
  best_idx <- which.min(grid_tbl$Jhat)
  optial_c=grid_tbl$c[best_idx]
  optial_alpha=grid_tbl$alpha[best_idx]
  optial_power=grid_tbl$power[best_idx]
  optial_c_total=grid_tbl$total_node[best_idx]
  opti_c_lis[i]<- optial_c
  opti_alph_lis[i]<-optial_alpha
  optial_power_lis[i]<-optial_power
  optial_c_total_lis[i]<-optial_c_total
  
  #66.07257   66.04092
  res <- abs_normal_q(c_seg = optial_c+1, total_node = optial_c_total, alpha =optial_alpha, power = optial_power)
  op_knot_vec <- res$q
  print(optial_c)
  
  
  modifide_huber_beta <- Modifid_huber_IRLS(
    beta_init = LMS_beta,
    X         = X_train,
    Y         = y_train,
    S_n       = Sn,
    knot_vec    = op_knot_vec,
    total_node=optial_c_total,
    eps       = 1e-3,
    max_iter  = 500
  )
  
  
  
  Proposed_MAPE[i] <- median(abs(y_test - X_test %*% modifide_huber_beta))
  Proposed_beta_list[[i]] <- as.vector(modifide_huber_beta)
}

## =========================================================
## 汇总结果：MAPE
## =========================================================

cat("LSE   MAPE mean, sd:\n");      print(c(mean(LSE_MAPE),      sd(LSE_MAPE)))
cat("LTS   MAPE mean, sd:\n");      print(c(mean(LTS_MAPE),      sd(LTS_MAPE)))
cat("Huber MAPE mean, sd:\n");      print(c(mean(Huber_MAPE),    sd(Huber_MAPE)))
cat("MM    MAPE mean, sd:\n");      print(c(mean(MM_MAPE),       sd(MM_MAPE)))
cat("Tukey MAPE mean, sd:\n");      print(c(mean(Tukey_MAPE),    sd(Tukey_MAPE)))
#cat("S     MAPE mean, sd:\n");      print(c(mean(S_MAPE),        sd(S_MAPE)))
cat("Prop. MAPE mean, sd:\n");      print(c(mean(Proposed_MAPE), sd(Proposed_MAPE)))
cat("LMS   MAPE mean, sd:\n");      print(c(mean(LMS_MAPE),      sd(LMS_MAPE))) 
cat("Huber ESL   MAPE mean, sd:\n");      print(c(mean(Huber_ESL_MAPE),      sd(Huber_ESL_MAPE))) 
cat("ESL  MAPE mean, sd:\n");      print(c(mean(ESL_MAPE),      sd(ESL_MAPE))) 
cat("\nSelected c mean:",     mean(opti_c_lis),     "\n")
cat("Selected truncation alpha mean:",  mean(opti_alph_lis), "\n")
cat("selected power parameter:",        mean(optial_power_lis),       "\n")
cat("selected grid size:",        mean(optial_c_total_lis),       "\n")


## =========================================================
## 汇总结果：各估计量 β 的均值和方差
## =========================================================

beta_summary <- function(beta_list) {
  beta_mat <- do.call(cbind, beta_list)   # p × K 矩阵
  list(
    mean = rowMeans(beta_mat),
    var  = apply(beta_mat, 1, var)
  )
}

LSE_beta_stat      <- beta_summary(LSE_beta_list)
LTS_beta_stat      <- beta_summary(LTS_beta_list)
LMS_beta_stat      <- beta_summary(LMS_beta_list)   
Huber_beta_stat    <- beta_summary(Huber_beta_list)
MM_beta_stat       <- beta_summary(MM_beta_list)
Tukey_beta_stat    <- beta_summary(Tukey_beta_list)
S_beta_stat        <- beta_summary(S_beta_list)
Proposed_beta_stat <- beta_summary(Proposed_beta_list)
