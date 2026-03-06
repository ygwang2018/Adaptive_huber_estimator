
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

trace_seed1 <- c(
  1,2,6,8,10,11,13,16,17,20,21,23,24,25,26,27,30,31,32,33,34,35,36,37,39,41,43,
  44,45,47,48,50,51,52,53,55,57,58,59,60,61,64,66,74,76,77,78,79,80,84,85,86,88,89,
  93,95,96,98,100,101,102,103,104,107,108,110,111,112,113,119,120,121,123,124,130,131,133,136,138,140,145,
  146,149,150,151,154,157,158,159,161,162,164,165,171,175,176,179,180,182,183,186,187,192,195,196,201,202,205,
  206,208,210,211,212,213,214,216,217,218,221,222,223,227,229,233,236,238,239,240,242,243,244,245,246,247,248,
  250,251,252,253,254,256,259,260,264,265,266,272,273,274,275,277,278,279,282,283,284,285,289,291,299,301,302,
  307,308,309,310,311,312,314,315,316,317,319,322,323,324,325,327,329,330,331,332,335,338,342,343,344,345,346,
  347,348,351,356,357,358,360,361,363,364,366,367,369,370,373,375,376,378,380,383,384,385,387,390,391,397,399,
  407,408,409,411,413,415,416,418,421,422,423,424,425,427,428,431,433,434,435,436,437,441,442,443,444,448,449,
  451,452,453,455,457,458,460,464,469,471,473,474,475,476,481,484,485,487,491,493,494,497,498,499,501,502,505,
  507,509,511,512,513,514,516,518,520,521,523,527,528,532,533,535,538,539,540,542,543,547,548,549,550,551,552,
  554,555,559,563,570,571,575,577,578,579,580,582,584,585,586,590,595,596,598,600,603,608,609,610,622,624,626,
  627,631,633,635,636,638,639,640,641,642,645,647,649,651,652,654,655,656,657,658,660,661,662,669,670,673,674,
  675,677,679,680,681,682,683,686,690,691,692,693,694,702,705,708,709,710,719,722,724,726,727,729,732,733,738,
  739,741,742,744,746,748,750,751,754,755,756,759,760,762,763,764,765,766,767,768,769,770,771,773,776,777,778,
  782,786,787,788,789,791,793,794,795,796,798,799,802,804,805,806,807,809,811,812,813,815,817,818,819,820,822,
  825,828,830,835,837,840,841,844,845,846,847,849,852,855,857,858,859,861,862,864,865,866,867,868,869,873,874,
  875,877,879,880,881,885,887,890,891,892,896,897,902,904,905,906,907,908,909,910,911,913,914,915,918,919,921,
  923,926,927,928,929,930,932,933,937,938,939,940,941,942,943,946,950,952,954,955,957,958,960,965,966,968,969,
  970,973,975,976,981,983,984,991,994,999
)



trace_seed2 <-c(
  1001,1002,1004,1005,1007,1009,1010,1011,1013,1014,1015,1016,1020,1021,1024,1025,1028,1030,1031,1033,1039,
  1040,1041,1042,1043,1047,1050,1051,1052,1055,1057,1063,1064,1067,1070,1071,1072,1073,1077,1079,1081,1083,
  1087,1088,1089,1090,1092,1094,1100,1101,1103,1104,1106,1107,1108,1111,1115,1116,1118,1122,1124,1127,1128,
  1129,1130,1132,1136,1137,1138,1139,1140,1142,1143,1145,1149,1150,1154,1156,1157,1158,1159,1164,1167,1168,
  1171,1172,1173,1174,1175,1177,1178,1180,1181,1184,1185,1186,1187,1189,1190,1191,1195,1196,1199,1200,1201,
  1204,1205,1207,1210,1211,1213,1215,1216,1217,1220,1221,1224,1225,1226,1227,1228,1230,1233,1234,1235,1236,
  1238,1239,1240,1241,1244,1245,1246,1247,1248,1252,1253,1255,1256,1258,1260,1261,1262,1263,1264,1266,1269,
  1270,1272,1274,1275,1276,1277,1280,1283,1284,1286,1287,1288,1289,1292,1294,1295,1296,1297,1300,1301,1303,
  1306,1308,1309,1310,1311,1312,1315,1316,1317,1318,1320,1321,1322,1323,1324,1325,1326,1327,1330,1331,1333,
  1334,1336,1341,1342,1343,1347,1348,1350,1352,1354,1355,1356,1357,1362,1363,1365,1366,1369,1370,1372,1374,
  1375,1376,1379,1383,1385,1386,1388,1390,1393,1397,1398,1399,1402,1404,1405,1407,1409,1411,1412,1413,1419,
  1420,1422,1424,1425,1426,1429,1430,1438,1439,1443,1444,1449,1452,1454,1456,1459,1460,1462,1463,1464,1465,
  1466,1468,1469,1470,1471,1472,1475,1476,1477,1478,1480,1481,1484,1489,1492,1493,1495,1496,1497,1498,1499,
  1502,1505,1507,1509,1510,1511,1512,1514,1515,1518,1519,1520,1521,1522,1527,1528,1529,1530,1531,1536,1537,
  1538,1539,1540,1543,1544,1545,1549,1550,1555,1556,1560,1561,1562,1563,1564,1565,1566,1567,1569,1573,1574,
  1575,1576,1578,1581,1582,1584,1585,1586,1587,1591,1592,1593,1595,1600,1601,1603,1604,1606,1608,1610,1612,
  1613,1614,1615,1618,1620,1622,1623,1630,1638,1643,1645,1646,1647,1649,1651,1652,1656,1657,1662,1663,1665,
  1666,1669,1672,1673,1674,1676,1680,1681,1682,1683,1684,1686,1688,1692,1694,1696,1697,1698,1699,1701,1702,
  1705,1706,1707,1708,1709,1710,1711,1713,1714,1715,1716,1722,1724,1727,1729,1731,1732,1734,1737,1739,1740,
  1741,1743,1744,1745,1748,1749,1751,1753,1756,1757,1758,1759,1761,1763,1769,1770,1774,1777,1780,1782,1783,
  1785,1786,1790,1791,1793,1794,1796,1797,1799,1800,1802,1803,1806,1808,1810,1812,1814,1816,1817,1819,1821,
  1824,1826,1827,1829,1830,1833,1834,1835,1836,1838,1839,1845,1847,1850,1851,1852,1855,1857,1859,1861,1862,
  1863,1864,1865,1866,1868,1869,1870,1877,1878,1879,1880,1881,1882,1883,1884,1885,1886,1888,1889,1891,1892,
  1894,1895,1899,1901,1902,1903,1909,1911,1912,1913,1915,1916,1917,1918,1919,1921,1922,1923,1926,1928,1929,
  1933,1934,1935,1938,1946,1947,1948,1950,1951,1952,1953,1954,1955,1957,1960,1961,1963,1965,1966,1975,1976,
  1977,1978,1980,1981,1983,1984,1986,1991,1992,1993,1994,1996,1997,1999
)
rest_seed_total<-c(trace_seed1,trace_seed2[2:478])

for (i in rest_seed_total) {
  set.seed(i)
  X <- matrix(rnorm(n * p), n, p)
  new_X <- cbind(1, X)
  err<-rnorm(n,0,1)
  
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



