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



## =========================
## Compute:
##  B = E[psi^2(Z)]  (your closed form)
##  A_naive = E[psi'(Z)] = P(|Z|<z0)
##  A_corr  = ∫ φ(u) dψ(u) = A_naive + jump contributions
##  ARE = A^2 / B
## where Z~N(0,1)
## =========================
calc_ARE_normal_modified <- function(c_seg, total_node,
                                     alpha = 1, power = 3) {
  
  z   <- abs_normal_q (c_seg+1, total_node, alpha = alpha, power =power)$q  # length c+1
  tau <- make_tau_vec(c_seg, total_node = total_node,z)                                # length c+1
  
  z0 <- z[1]
  
  ## ---- B = E psi^2(Z) ----
  term0 <- -2 * z0 * dnorm(z0) + 2 * pnorm(z0) - 1
  
  ## segments i=0,...,c-1 use tau_i ; in R: tau[1:c]
  seg_prob <- pnorm(z[2:(c_seg + 1L)]) - pnorm(z[1:c_seg])
  term_seg <- 2 * sum((tau[1:c_seg]^2) * seg_prob)
  
  B <- term0 + term_seg
  
  ## ---- A_naive (what you currently use) ----
  A_naive <- 2 * pnorm(z0) - 1
  
  ## ---- A_corr: add jump contributions ----
  ## ψ(u)=u on |u|<=z0, then ψ(u)=tau_i sign(u) on (z_i,z_{i+1}], i=0..c-1,
  ## and ψ(u)=0 for |u|>z_c.
  ##
  ## Jump at z0 if tau0 != z0: Δψ(z0)=tau0 - z0
  #jump_z0 <- (tau[1] - z0)
  
  ## Jumps at z_i for i=1,...,c-1: Δψ(z_i)=tau_i - tau_{i-1}
  ## here tau_i is tau[i+1] in R; but for segments we use tau[1:c]
  jump_mid <- diff(tau[1:c_seg])          # length c-1, corresponds to z_1,...,z_{c-1}
  
  ## Jump at z_c: from tau_{c-1} to 0 => Δψ(z_c)=0 - tau_{c-1}
  jump_tail <- (0 - tau[c_seg])
  
  ## Combine symmetric ±z contributions => multiply by 2 * φ(z_k)
  A_corr <- A_naive +
    2 * sum(jump_mid * dnorm(z[2:c_seg])) +
    2 * jump_tail * dnorm(z[c_seg+1])
  
  ## ---- AREs ----
  ARE_corr  <- (A_corr^2)  / B
  
  list(
    ARE_corr  = ARE_corr
  )
}

## =========================
## Example
## =========================
c_total_lis=c(2000,4000,10000,20000,30000,40000,50000,60000)
ARE_lis1<-c()
for (c_total in c_total_lis){
  res <- calc_ARE_normal_modified(
    c_seg = c_total-1,
    total_node = c_total,
    alpha = 1,     # truncation
    power = 3      # power
)
  ARE_lis1<-c(ARE_lis1,res$ARE_corr)
}

ARE_lis1


power_lis=c(1.1,1.5,1.8,2,2.2,2.5,3,3.2)
ARE_lis2<-c()

for (power in power_lis){
  res <- calc_ARE_normal_modified(
    c_seg = 49999,
    total_node = 50000,
    alpha = 1,     # truncation
    power = power      # power
  )
  ARE_lis2<-c(ARE_lis2,res$ARE_corr)
}

c_lis=c(1000,2000,5000,8000,10000,20000,30000,49999)


ARE_lis3<-c()

for (c in c_lis){
  res <- calc_ARE_normal_modified(
    c_seg = c,
    total_node = 50000,
    alpha = 1,     #truncation
    power = 1.1      # power
  )
  ARE_lis3<-c(ARE_lis3,res$ARE_corr)
}

ARE_lis3

trunca_lis=c(0.8,0.83,0.85,0.9,0.93,0.95,0.99,1)
ARE_lis4<-c()

for (truc_para in trunca_lis){
  res <- calc_ARE_normal_modified(
    c_seg = 49999,
    total_node = 50000,
    alpha = truc_para,     # truncation
    power = 1.1      # power
  )
  ARE_lis4<-c(ARE_lis4,res$ARE_corr)
}

ARE_lis4
