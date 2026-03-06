## ============================================================
## 1) Your abs_normal_q() (keep as-is)
## ============================================================
abs_normal_q <- function(c_seg, total_node, alpha = 0.95, power = 3) {
  n <- as.integer(total_node)
  c_seg <- as.integer(c_seg)
  
  if (c_seg < 1L) stop("c_seg must be >= 1")
  if (alpha <= 0 || alpha > 1) stop("alpha must be in (0,1]")
  if (power <= 1) stop("power must be > 1 (controls decreasing gaps)")
  
  n_max <- max(1L, floor(alpha * n))
  # if (n_max < c_seg) stop("alpha*n is too small: n_max < c_seg")  # you commented
  
  t <- seq(0, 1, length.out = c_seg)
  i <- 1L + floor((n_max - 1L) * (1 - (1 - t)^power))
  
  p <- 0.5 + i / (2 * (n + 1L))
  q <- qnorm(p)
  
  list(i = i, q = q)
}



## ============================================================
## 3) Build slopes tau_0,...,tau_c  (your code)
##    tau_0 = tau0, tau_i = tau0 + qnorm((u+1)/2), u=i/(total_node+1)
## ============================================================
make_tau_vec <- function(c_seg, total_node,knot_vec) {
  c_seg <- as.integer(c_seg)
  total_node <- as.integer(total_node)
  tau0=knot_vec[1]
  tau_vec <- numeric(c_seg + 1L)
  tau_vec[1L] <- tau0
  
  for (i in 1:c_seg) {
    u <- i / (total_node + 1)
    tau_vec[i + 1L] <- qnorm((u + 1) / 2) + tau0
  }
  tau_vec
}

calc_ARE_cauchy_modified <- function(c_seg, total_node,
                                     alpha = 1, power = 3) {
  c_seg <- as.integer(c_seg)
  total_node <- as.integer(total_node)
  if (c_seg < 1L) stop("c_seg must be >= 1")
  
  ## --- knots & taus: keep your construction ---
  z   <- abs_normal_q (c_seg+1, total_node, alpha = alpha, power =power)$q  # length c+1
  tau <- make_tau_vec(c_seg, total_node = total_node,z)                                # length c+1
  
  
  
  z0 <- z[1L]
  zc <- z[c_seg + 1L]
  
  ## ============================================================
  ## B(theta) under Z ~ t(1) (standard Cauchy)
  ## E[psi^2(Z)] = 2/pi (z0 - atan z0) + 2/pi sum_{i=0}^{c-1} tau_i^2 (atan z_{i+1} - atan z_i)
  ## indices: tau_i = tau[i+1], z_i = z[i+1]
  ## ============================================================
  term0 <- (2 / pi) * (z0 - atan(z0))
  
  term_seg <- (2 / pi) * sum(
    (tau[1:c_seg]^2) * (atan(z[2:(c_seg + 1L)]) - atan(z[1:c_seg]))
  )
  
  B <- term0 + term_seg
  
  ## ============================================================
  ## A(theta) under Z ~ t(1)
  ## A = 2/pi atan(z0) + 2/pi sum_{i=1}^{c-1} (tau_i - tau_{i-1})/(1+z_i^2) - 2/pi * tau_{c-1}/(1+z_c^2)
  ## Here z_i = z[i+1], tau_i = tau[i+1].
  ## Note: your construction sets tau0 = z0, so no jump term at z0 is needed.
  ## ============================================================
  A_naive <- (2 / pi) * atan(z0)  # P(|Z|<z0) under Cauchy
  
  if (c_seg >= 2L) {
    jump_mid <- diff(tau[1:c_seg])           # (tau_1-tau_0),...,(tau_{c-1}-tau_{c-2})
    A_mid <- (2 / pi) * sum(jump_mid / (1 + z[2:c_seg]^2))  # z_1,...,z_{c-1}
  } else {
    A_mid <- 0
  }
  
  A_tail <- -(2 / pi) * (tau[c_seg] / (1 + zc^2))  # - 2/pi * tau_{c-1}/(1+z_c^2)
  
  A_corr <- A_naive + A_mid + A_tail
  
  ## ============================================================
  ## Asymptotic variance factor and "efficiency"
  ## J = B/A^2 (smaller better), Eff = A^2/B (larger better)
  ## ============================================================
  J_naive <- B / (A_naive^2)
  J_corr  <- B / (A_corr^2)
  
  ARE_naive <- (A_naive^2) / B
  ARE_corr  <- (A_corr^2)  / B
  
  list(
    ARE_corr  = ARE_corr
  )
}



## =========================
## Example
## =========================
c_total_lis=c(20,30,50,60,70,80,100,200,500,700,1000,2000)
ARE_lis1<-c()
for (c_total in c_total_lis){
  res <- calc_ARE_cauchy_modified(
    c_seg = c_total/2,
    total_node = c_total,
    alpha = 1,     # truncation
    power = 3      # power
  )
  ARE_lis1<-c(ARE_lis1,res$ARE_corr)
}

ARE_lis1


power_lis=c(1.1,1.2,1.5,1.8,2,2.2,2.5,2.7,3,3.2,3.4,3.6)
ARE_lis2<-c()

for (power in power_lis){
  res <- calc_ARE_cauchy_modified(
    c_seg = 25,
    total_node = 50,
    alpha = 1,     # truncation
    power = power      # power
  )
  ARE_lis2<-c(ARE_lis2,res$ARE_corr)
}

ARE_lis2




c_lis=c(50/13,50*(2/13),50*(3/13),50*(4/13),50*(5/13),50*(6/13),50*(7/13),50*(8/13),50*(9/13),50*(10/13),50*(11/13),50*(12/13))


ARE_lis3<-c()

for (c in c_lis){
  res <- calc_ARE_cauchy_modified(
    c_seg = c,
    total_node =50,
    alpha = 1,     # truncation
    power = 1.1      # power
  )
  ARE_lis3<-c(ARE_lis3,res$ARE_corr)
}

ARE_lis3

trunca_lis=c(0.7,0.73,0.75,0.76,0.8,0.83,0.85,0.9,0.93,0.95,0.99,1)
ARE_lis4<-c()

for (truc_para in trunca_lis){
  res <- calc_ARE_cauchy_modified(
    c_seg = 49,
    total_node = 50,
    alpha = truc_para,     # truncation
    power = 1.1      # power
  )
  ARE_lis4<-c(ARE_lis4,res$ARE_corr)
}

ARE_lis4
