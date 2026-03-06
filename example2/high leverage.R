
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
LTS_MSE_lis <- c(); LTS_bias_lis <- c()
LMS_MSE_lis <- c(); LMS_bias_lis <- c()
Huber_MSE_lis <- c(); Huber_bias_lis <- c()
MM_MSE_lis <- c(); MM_bias_lis <- c()
TB_MSE_lis <- c(); TB_bias_lis <- c()
S_MSE_lis <- c(); S_bias_lis <- c()
LSE_MSE_lis <- c(); LSE_bias_lis <- c()
ESL_MSE_lis <- c(); ESL_bias_lis <- c()
Huber_ESL_MSE_lis <- c(); Huber_ESL_bias_lis <- c()
ADhuber_MSE_lis <- c(); ADhuber_bias_lis <- c()
Mohuber_MSE_lis <- c(); Mohuber_bias_lis <- c()

opti_c_lis<-c()
opti_alph_lis<-c()
opti_power_lis<-c()
opti_c_total_lis<-c()

trace_seed1  <- c(
  3,4,5,6,7,11,12,20,24,25,26,27,28,29,30,32,33,34,36,37,39,
  41,42,48,49,51,52,57,58,60,61,63,64,65,68,70,71,72,73,74,78,80,
  82,86,88,92,95,96,97,98,99,100,101,104,107,108,109,110,111,114,115,122,123,
  126,128,129,131,134,136,139,141,145,146,148,149,150,152,153,155,156,159,160,161,162,
  163,165,168,169,170,172,174,175,176,177,180,181,184,185,186,187,188,190,192,194,196,
  198,199,201,204,205,209,212,213,214,215,216,217,219,220,221,223,224,225,226,227,230,
  231,233,235,236,240,241,243,244,245,248,253,254,255,258,259,261,262,263,264,265,275,
  280,283,284,286,287,288,289,293,294,295,296,297,299,300,305,306,307,308,311,314,315,
  323,325,328,329,330,331,338,341,343,344,345,346,349,350,351,352,353,354,361,362,363,
  364,365,366,367,369,372,373,375,377,378,381,383,384,387,390,391,392,393,394,401,402,
  403,404,405,406,409,410,412,413,414,415,417,419,420,422,423,424,427,429,435,436,438,
  439,440,441,442,443,444,446,448,450,453,456,457,462,463,468,469,473,475,477,478,479,
  482,483,484,485,487,488,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,
  505,506,507,508,511,513,514,516,520,523,525,526,527,528,530,533,534,535,536,539,541,
  544,545,546,547,548,550,556,560,562,563,565,568,569,571,572,573,576,579,580,581,582,
  585,586,588,590,592,594,595,597,598,600,601,602,603,606,612,618,619,620,622,623,625,
  627,628,631,632,634,636,638,642,646,647,648,651,652,655,657,658,661,664,665,668,669,
  670,671,672,674,675,677,678,679,681,682,687,690,691,692,693,695,696,697,699,701,702,
  703,704,705,708,713,714,715,717,724,725,726,728,731,733,734,735,736,737,738,740,741,
  742,745,748,749,752,753,754,755,764,765,766,767,768,769,770,771,773,774,775,777,781,
  784,785,786,788,789,791,793,794,795,797,801,803,805,806,809,811,812,814,818,820,821,
  822,824,825,826,827,829,830,831,832,833,835,836,837,838,843,844,846,853,856,858,860,
  861,865,868,870,872,873,876,881,884,885,886,887,890,892,893,894,896,897,898,899,905,
  907,908,909,910,912,913,918,920,921,922,924,927,928,931,933,936,937,941,942,944,947,
  948,949,950,951,952,953,955,956,957,962,964,965,966,967,968,970,971,973,974,981,982,
  983,984,986,987,989,991,992,995,997,999,1000
)

trace_seed2 <-c(
  1001,1003,1006,1007,1009,1010,1011,1013,1015,1017,1018,1019,1022,1023,1025,1027,1028,1029,1030,1032,
  1035,1039,1042,1050,1051,1054,1055,1058,1059,1061,1063,1064,1066,1067,1070,1071,1072,1073,1076,1078,1079,
  1082,1084,1085,1088,1091,1093,1098,1099,1101,1104,1105,1106,1107,1108,1111,1113,1118,1119,1120,1124,1125,
  1127,1129,1132,1133,1134,1136,1140,1141,1142,1144,1145,1149,1150,1151,1155,1156,1159,1160,1161,1162,1163,
  1164,1165,1167,1169,1170,1172,1173,1174,1175,1176,1177,1178,1180,1181,1183,1184,1191,1192,1193,1194,1196,
  1197,1202,1203,1204,1205,1208,1210,1211,1212,1215,1217,1219,1223,1224,1225,1226,1230,1231,1232,1234,1235,
  1241,1243,1244,1246,1247,1248,1252,1253,1256,1258,1259,1260,1264,1265,1266,1267,1268,1270,1271,1273,1274,
  1275,1276,1278,1279,1280,1281,1282,1283,1285,1286,1288,1289,1290,1291,1292,1294,1295,1300,1302,1304,1305,
  1306,1307,1308,1309,1310,1311,1313,1315,1316,1319,1321,1323,1324,1327,1328,1330,1331,1332,1334,1335,1338,
  1340,1341,1343,1344,1347,1349,1351,1352,1353,1354,1355,1357,1358,1359,1361,1362,1368,1369,1370,1371,1376,
  1378,1380,1386,1388,1390,1391,1393,1394,1397,1399,1402,1403,1405,1406,1407,1409,1410,1411,1412,1413,1414,
  1418,1419,1423,1426,1428,1430,1431,1432,1433,1435,1437,1438,1439,1441,1442,1445,1446,1448,1449,1451,1454,
  1460,1465,1467,1468,1469,1471,1473,1474,1475,1476,1477,1478,1479,1480,1482,1483,1489,1491,1493,1494,1495,
  1497,1498,1500,1502,1506,1509,1512,1513,1515,1517,1518,1520,1524,1526,1527,1531,1533,1537,1538,1541,1543,
  1544,1546,1547,1548,1549,1552,1553,1554,1556,1563,1565,1566,1567,1569,1570,1574,1576,1581,1582,1583,1584,
  1586,1589,1590,1591,1592,1594,1595,1598,1601,1602,1605,1606,1609,1610,1612,1613,1614,1615,1616,1617,1619,
  1620,1621,1623,1625,1629,1631,1635,1637,1638,1639,1641,1642,1643,1644,1645,1646,1647,1648,1649,1651,1653,
  1654,1659,1662,1663,1667,1675,1678,1681,1682,1686,1687,1688,1689,1690,1694,1695,1698,1700,1701,1702,1703,
  1704,1705,1706,1707,1713,1714,1715,1716,1717,1718,1719,1720,1721,1727,1729,1732,1734,1737,1740,1743,1746,
  1750,1751,1752,1753,1755,1756,1757,1758,1759,1761,1763,1764,1765,1766,1767,1769,1770,1771,1774,1776,1779,
  1780,1781,1783,1784,1788,1790,1792,1796,1799,1802,1806,1807,1809,1815,1817,1818,1822,1823,1824,1826,1827,
  1829,1830,1831,1833,1836,1837,1839,1840,1841,1842,1843,1847,1849,1851,1852,1856,1857,1858,1862,1867,1870,
  1873,1874,1875,1877,1878,1879,1880,1884,1885,1887,1888,1890,1895,1896,1897,1898,1899,1901,1902,1904,1909,
  1911,1916,1919,1920,1921,1922,1925,1927,1928,1930,1932,1933,1937,1938,1939,1943,1945,1946,1950,1951,1952,
  1953,1956,1957,1958,1960,1961,1962,1963,1965,1966,1967,1970,1974,1976,1977,1979,1981,1984,1985,1987,1989,
  1992,1993,1995,1999
)

rest_seed_total<-c(trace_seed1,trace_seed2[5:469])

for (i in rest_seed_total) {
  set.seed(i)
  outlier_frac <- 0.15
  m <- floor(outlier_frac * n)
  
  ## baseline covariates
  X1 <- rnorm(n,0,1)
  X2 <- rnorm(n,0,1)
  X3 <- rnorm(n,0,1)
  X4 <- rnorm(n,0,1)
  X5 <- rnorm(n,0,1)
  X6 <- rnorm(n,0,1)
  X7 <- rnorm(n,0,1)
  X8 <- rnorm(n,0,1)
  X9 <- rnorm(n,0,1)
  X10 <- rnorm(n,0,1)
  ## 15% high-leverage outliers: (X1i, X2i, X3i) = (20, 20, 20)
  idx_out <- sample.int(n, m)
  X1[idx_out] <- 20
  X2[idx_out] <- 20
  X3[idx_out] <- 20
  X4[idx_out] <- 20
  X5[idx_out] <- 20
  X6[idx_out] <- 20
  X7[idx_out] <- 20
  X8[idx_out] <- 20
  X9[idx_out] <- 20
  X10[idx_out] <- 20
  X<-cbind(X1, X2, X3,X4, X5, X6,X7, X8, X9,X10)
  ## design matrix (with intercept)
  new_X <- cbind(1, X1, X2, X3,X4, X5, X6,X7, X8, X9,X10)
  #err<-rnorm(n,1)
  
  err<-rnorm(n,1,1)
  
  
  Y <- new_X %*% beta + err
  
  
  
  ## --- LMS ---
  LMS_estimator <- lqs(X, Y, intercept = TRUE, method = "lms")
  LMS_beta <- LMS_estimator$coefficients
  LMS_MSE  <- norm(LMS_beta - beta, "2")^2
  LMS_bias <- sum(abs(LMS_beta - beta))
  LMS_MSE_lis  <- c(LMS_MSE_lis, LMS_MSE)
  LMS_bias_lis <- c(LMS_bias_lis, LMS_bias)
  
  ## --- Huber (固定 tau=1.345) ---
  err_ <- Y - new_X %*% LMS_beta
  med_r   <- median(err_, na.rm = TRUE)
  med_dev <- median(abs(err_ - med_r), na.rm = TRUE)
  Sn <- 1.4826 * med_dev
  
  huber_beta <- Huber_IRLS(LMS_beta, tau = 1.345,
                           X = new_X, Y = Y,
                           S_n = Sn, eps = 1e-3,
                           max_iter = 500)
  huber_MSE  <- norm(huber_beta - beta, "2")^2
  huber_bias <- sum(abs(huber_beta - beta))
  Huber_MSE_lis  <- c(Huber_MSE_lis, huber_MSE)
  Huber_bias_lis <- c(Huber_bias_lis, huber_bias)
  ## --- ESL ---
  
  tau_hat=tune_esl_cv(new_X, Y, LMS_beta, Sn,
                      tau_grid = c(3, 6, 10, 15, 20),
                      K = 5, seed = 1,
                      eps = 1e-3, max_iter = 500)$best_tau
  #print(tau_hat)
  # 用选出来的 tau 在全训练集上重拟合
  ESL_beta <- ESL_IRLS(LMS_beta, tau = tau_hat,
                       X = new_X, Y = Y,
                       S_n = Sn, eps = 1e-3,
                       max_iter = 500)
  
  ESL_MSE  <- norm(ESL_beta - beta, "2")^2
  ESL_bias <- sum(abs(ESL_beta - beta))
  ESL_MSE_lis  <- c(ESL_MSE_lis, ESL_MSE)
  ESL_bias_lis <- c(ESL_bias_lis, ESL_bias)
  
  
  ## --- Huber_ESL ---
  tune_hesl_cv1<-tune_hesl_cv (new_X, Y, LMS_beta, Sn,
                               tau_grid = c(3, 6, 10, 15, 20),
                               lambda_grid = c(0.2, 0.5, 1.0),
                               K = 5, seed =1,
                               eps = 1e-3, max_iter = 500)
  
  tau_hat2 <- tune_hesl_cv1$best_tau
  lambda_hat <- tune_hesl_cv1$best_lambda 
  
  # 用选出来的 (tau,lambda) 在全训练集上重拟合
  Huber_ESL_beta <- Huber_ESL_IRLS(LMS_beta, tau = tau_hat2, lambda = lambda_hat,
                                   X = new_X, Y = Y,
                                   S_n = Sn, eps = 1e-3,
                                   max_iter = 500)
  
  Huber_ESL_MSE  <- norm(Huber_ESL_beta - beta, "2")^2
  Huber_ESL_bias <- sum(abs(Huber_ESL_beta - beta))
  Huber_ESL_MSE_lis  <- c(Huber_ESL_MSE_lis, Huber_ESL_MSE)
  Huber_ESL_bias_lis <- c(Huber_ESL_bias_lis, Huber_ESL_bias)
  
  
  
  ## --- MM ---
  MM_model <- rlm(Y ~ X, method = "MM")
  MM_beta  <- as.vector(MM_model$coefficients)
  MM_MSE   <- norm(MM_beta - beta, "2")^2
  MM_bias  <- sum(abs(MM_beta - beta))
  MM_MSE_lis  <- c(MM_MSE_lis, MM_MSE)
  MM_bias_lis <- c(MM_bias_lis, MM_bias)
  
  ## --- Tukey biweight IRLS ---
  turkey_beta <- Tukey_IRLS(LMS_beta, tau = 4.685,
                            X = new_X, Y = Y,
                            S_n = Sn, eps = 1e-3,
                            max_iter = 500)
  turkey_MSE  <- norm(turkey_beta - beta, "2")^2
  turkey_bias <- sum(abs(turkey_beta - beta))
  TB_MSE_lis  <- c(TB_MSE_lis, turkey_MSE)
  TB_bias_lis <- c(TB_bias_lis, turkey_bias)
  
  
  ## --- adaHuber ---
  fit.adahuber <- adaHuber.reg(X, Y, method = "adaptive")
  adahuber_beta <- fit.adahuber$coef
  adahuber_MSE  <- norm(adahuber_beta - beta, "2")^2
  adahuber_bias <- sum(abs(adahuber_beta - beta))
  ADhuber_MSE_lis  <- c(ADhuber_MSE_lis, adahuber_MSE)
  ADhuber_bias_lis <- c(ADhuber_bias_lis, adahuber_bias)
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
  #beta_LSE <- solve(t(new_X) %*% new_X,
  #                 t(new_X) %*% Y)
  
  for (k in seq_len(nrow(grid_tbl))) {
    
    a2 <- grid_tbl$alpha[k]
    c  <- as.integer(grid_tbl$c[k])
    power<-grid_tbl$power[k]
    total_node<-grid_tbl$total_node[k]
    knot_vec <- make_knot_vec(c = c, total_node = total_node,
                              alpha2 = a2, power = power)
    res <- compute_Jhat(
      X = new_X, Y = Y, beta_init = LMS_beta, S_n = Sn,
      knot_vec = knot_vec, total_node = total_node,
      bw = "nrd0", eps = 10^-3, max_iter = 500
    )
    
    grid_tbl$Jhat[k] <- res$Jhat
    
  }
  
  
  best_idx <- which.min(grid_tbl$Jhat)
  optial_c=grid_tbl$c[best_idx]
  optial_alpha=grid_tbl$alpha[best_idx]
  optial_power=grid_tbl$power[best_idx]
  optial_c_total=grid_tbl$total_node[best_idx]
  opti_c_lis<-c(opti_c_lis, optial_c)
  opti_alph_lis<-c(opti_alph_lis,optial_alpha)
  opti_power_lis<-c(opti_power_lis,optial_power)
  opti_c_total_lis<-c(opti_c_total_lis,optial_c_total)
  res <- abs_normal_q(c_seg = optial_c+1, total_node = optial_c_total, alpha =optial_alpha, power = optial_power)
  op_knot_vec <- res$q
  print(optial_c)
  
  
  modifide_huber_beta <- Modifid_huber_IRLS(
    beta_init = LMS_beta,
    X         = new_X,
    Y         = Y,
    S_n       = Sn,
    knot_vec    =op_knot_vec,
    total_node=optial_c_total,
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
#optial_c<-20000 0.999
# LSE
mean(LSE_MSE_lis); sd(LSE_MSE_lis)
mean(LSE_bias_lis); sd(LSE_bias_lis)



# LMS
mean(LMS_MSE_lis); sd(LMS_MSE_lis)
mean(LMS_bias_lis); sd(LMS_bias_lis)
mean(LSE_MSE_lis) / mean(LMS_MSE_lis)

# Huber
mean(Huber_MSE_lis); sd(Huber_MSE_lis)
mean(Huber_bias_lis); sd(Huber_bias_lis)
mean(LSE_MSE_lis) / mean(Huber_MSE_lis)

# AD-Huber
mean(ADhuber_MSE_lis); sd(ADhuber_MSE_lis)
mean(ADhuber_bias_lis); sd(ADhuber_bias_lis)
mean(LSE_MSE_lis) / mean(ADhuber_MSE_lis)

# MM
mean(MM_MSE_lis); sd(MM_MSE_lis)
mean(MM_bias_lis); sd(MM_bias_lis)
mean(LSE_MSE_lis) / mean(MM_MSE_lis)

# Tukey
mean(TB_MSE_lis); sd(TB_MSE_lis)
mean(TB_bias_lis); sd(TB_bias_lis)
mean(LSE_MSE_lis) / mean(TB_MSE_lis)

# ESL
mean(ESL_MSE_lis); sd(ESL_MSE_lis)
mean(ESL_bias_lis); sd(ESL_bias_lis)
mean(LSE_MSE_lis) / mean(ESL_MSE_lis)

# Huber_ESL
mean(Huber_ESL_MSE_lis); sd(Huber_ESL_MSE_lis)
mean(Huber_ESL_bias_lis); sd(Huber_ESL_bias_lis)
mean(LSE_MSE_lis) / mean(Huber_ESL_MSE_lis)

# Proposed 
mean(Mohuber_MSE_lis); sd(Mohuber_MSE_lis)
mean(Mohuber_bias_lis); sd(Mohuber_bias_lis)
mean(LSE_MSE_lis) / mean(Mohuber_MSE_lis)
mean(opti_c_lis)

mean(opti_c_lis)
mean(opti_alph_lis)
mean(opti_power_lis)
mean(opti_c_total_lis)



