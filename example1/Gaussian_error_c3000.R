
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
    small <- abs(r) < 1e-8
    
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





library(robustbase)
library(MASS)
library(VGAM)
library(L1pack)
library(adaHuber)
library(quantreg)

n  <- 500
p  <- 1
beta <- c(1, 2)

n_rep <- 1000

total_Huber_eff_lis <- c(); total_Huber_bias_lis <- c()
total_Mohuber_eff_lis <- c(); total_Mohuber_bias_lis <- c()

tau_lis=seq(0.1,5,0.3)
for (tau in tau_lis){
  Huber_MSE_lis <- c(); Huber_bias_lis <- c()
  LSE_MSE_lis <- c(); LSE_bias_lis <- c()
  Mohuber_MSE_lis <- c(); Mohuber_bias_lis <- c()
  for (i in 1:n_rep) {
    set.seed(i)
    
    X <- matrix(rnorm(n * p), n, p)
    new_X <- cbind(1, X)
    err<-rnorm(n,0,1)
    
    #opera <- rbinom(n, 1, 0.6)
    #err <- c(rt(sum(opera == 1), 1.2),
    # rnorm(n - sum(opera == 1), 0, 2))
    
    Y <- new_X %*% beta + err
    
    ## --- LTS ---
    lts_model <- ltsReg(Y ~ X)
    LTS_beta  <- coef(lts_model)
    err_ <- Y - new_X %*% LTS_beta
    med_r   <- median(err_, na.rm = TRUE)
    med_dev <- median(abs(err_ - med_r), na.rm = TRUE)
    Sn <- 1.4826 * med_dev
    
    
    ## --- Huber (固定 tau=1.345) ---
    
    
    huber_beta <- Huber_IRLS(LTS_beta, tau = tau,
                             X = new_X, Y = Y,
                             S_n = Sn, eps = 1e-3,
                             max_iter = 500)
    huber_MSE  <- norm(huber_beta - beta, "2")^2
    huber_bias <- sum(abs(huber_beta - beta))
    Huber_MSE_lis  <- c(Huber_MSE_lis, huber_MSE)
    Huber_bias_lis <- c(Huber_bias_lis, huber_bias)
    
    
    
    
    ## -----modified huber-------
    
    res <- abs_normal_q(c_seg = 3000, total_node =50000, alpha =1, power = 1.1)$q
    c=3000
    knot_vec <- numeric(c + 1)
    knot_vec[1] <- tau                   # |z|_[0] = tau
    knot_vec[2:(c + 1)] <- tau + res[1:c]  # |z|_[i] = tau + q_i
    modifide_huber_beta <- Modifid_huber_IRLS(
      beta_init = LTS_beta,
      X         = new_X,
      Y         = Y,
      S_n       = Sn,
      knot_vec    = knot_vec,
      total_node = 50000,
      eps       = 1e-3,
      max_iter  = 500
    )
    modifide_huber_MSE  <- norm(beta - modifide_huber_beta, "2")^2
    modifide_huber_bias <- sum(abs(beta - modifide_huber_beta))
    Mohuber_MSE_lis <- c(Mohuber_MSE_lis,modifide_huber_MSE)
    Mohuber_bias_lis <- c(Mohuber_bias_lis,modifide_huber_bias)
    ## --- LSE ---
    beta_LSE <- solve(t(new_X) %*% new_X,
                      t(new_X) %*% Y)
    LSE_MSE  <- norm(beta - beta_LSE, "2")^2
    LSE_bias <- sum(abs(beta - beta_LSE))
    LSE_MSE_lis  <- c(LSE_MSE_lis, LSE_MSE)
    LSE_bias_lis <- c(LSE_bias_lis, LSE_bias)
    
    if (i %% 50 == 0) cat("i =", i, "\n")
  }
  total_Huber_eff_lis <- c(total_Huber_eff_lis,mean(LSE_MSE_lis) / mean(Huber_MSE_lis))
  total_Huber_bias_lis <- c(total_Huber_bias_lis,mean(Huber_bias_lis))
  
  total_Mohuber_eff_lis <- c(total_Mohuber_eff_lis,mean(LSE_MSE_lis) / mean(Mohuber_MSE_lis)) 
  total_Mohuber_bias_lis <- c(total_Mohuber_bias_lis,mean(Mohuber_bias_lis))
  
}
tau_lis <- seq(0.1, 5, 0.3)

df <- data.frame(
  tau = tau_lis,
  Huber_eff = total_Huber_eff_lis ,
  Huber_Bias = total_Huber_bias_lis,
  MoHuber_eff = total_Mohuber_eff_lis,
  MoHuber_Bias = total_Mohuber_bias_lis
)

write.csv(df, "C://Users//PC//Desktop//modified_huber//example2//eff_bias_summary_c3000.csv", row.names = FALSE)

