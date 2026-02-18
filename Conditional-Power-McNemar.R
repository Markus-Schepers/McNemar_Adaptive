
lower <- function(n,alpha=0.05){
  if(pbeta(0.5,n,1) > alpha/2) return(NA)
  for(k in 1:n){
    cdf_val <- pbeta(0.5,n-k,k+1)
    if(cdf_val > alpha/2){
      return(k-1)
    }
  }
}

upper <- function(n,alpha=0.05){
  if(pbeta(0.5,1,n) < 1-alpha/2) return(NA)
  for(k in (n-1):0){
    cdf_val <- pbeta(0.5,n-k+1,k)
    if( 1- cdf_val > alpha/2){
      return(k+1)
    }
  }
}


# exact power conditioned on a certain number of discordant pairs nd
exact_power <- function(p_star,nd,alpha=0.05){
  if(nd==0) return(0)
  k1 <- lower(nd,alpha)
  k2 <- upper(nd,alpha)
  if(is.na(k1)) return(0)
  if(is.na(k2)) return(0)
  return(pbinom(k1,nd,p_star)+1-pbinom(k2-1,nd,p_star))
}

# exact power conditioned on number of discordant pairs nd1 in first stage, 
# using size of first stage and size of second stage
# using the probability for a discordant pair
conditional_power <- function(p_star,psi, nd1,nsecond=38,alpha=0.05){
  nd2_vals <- 0:nsecond
  sum(vapply(nd2_vals,function(nd2){ return(exact_power(p_star,(nd1+nd2),alpha) * dbinom(nd2,nsecond,psi)) },numeric(1))) 
}


min_nsecond <- function(p_star, psi, nd1, target_power = 0.8, 
                        alpha = 0.05, nsecond_max = 200) {
  # 1. Check feasibility at maximum nsecond
  cp_max <- conditional_power(p_star = p_star, nd1 = nd1,
                              nsecond = nsecond_max,
                              psi = psi, alpha = alpha)
  
  if(is.na(cp_max) || cp_max < target_power){
    return(NA)   # impossible to reach target power
  }
  
  # 2. Search for smallest n2 achieving the target
  for(n2 in 0:nsecond_max){
    cp <- conditional_power(p_star = p_star, nd1 = nd1,
                            nsecond = n2,
                            psi = psi, alpha = alpha)
    
    if(!is.na(cp) && cp >= target_power){
      return(n2)
    }
  }
  
  # Should never happen because feasibility was checked
  return(NA)
}




# 2x2 table
# rows: surrogate (1=disc, 0=conc)
# cols: primary     (1=disc, 0=conc)
estimate_theta <- function(n11,
                           n10,
                           n01,
                           n00,
                           pseudo = 0.5) {
  
  theta1 <- (n11 + pseudo) / (n11 + n10 + 2*pseudo)
  theta0 <- (n01 + pseudo) / (n01 + n00 + 2*pseudo)
  
  return(c(theta_d = theta1, theta_c = theta0))
}



# when using the surrogate endpoint, by the law of total probability condition on all the number of discordant pairs in the primary endpoint at interim analysis in terms of transitions from surrogate to primary endpoint, and use the conditional power for the primary endpoint
# p_star is the effect size (once the fraction of discordant pairs psi is fixed)
power_surrogate <- function(p_star, psi,
                            theta_d, theta_c,
                            nd1,
                            n1,
                            m1,
                            m_d_SP,
                            n2,
                            alpha) {
  
  # pmfs
  p_X <- dbinom(0:n2, size = n2, prob = psi)
  p_Y <- dbinom(0:m_d_SP, size = m_d_SP, prob = theta_d)
  p_Z <- dbinom(0:(m1 - n1 - m_d_SP),
                size = m1 - n1 - m_d_SP, prob = theta_c)
  
  # convolution: Y + Z
  p_YZ <- convolve(p_Y, rev(p_Z), type = "open")
  
  # convolution: X + Y + Z
  p_XYZ <- convolve(p_X, rev(p_YZ), type = "open")
  
  # support of X + Y + Z
  k_vals <- 0:(length(p_XYZ) - 1)
  
  # conditional power evaluated once per total
  cp_vals <- vapply(
    k_vals+nd1,
    function(k)
      exact_power(
        p_star = p_star,
        nd = k,
        alpha = alpha
      ),
    numeric(1)
  )
  sum(cp_vals * p_XYZ)
}





min_nsecond_surrogate <- function(
    p_star, psi,
    theta_d, theta_c,
    nd1,
    n1,
    m1,
    m_d_SP,
    target_power = 0.8,
    alpha = 0.05,
    nsecond_max = 200
) {
  
  # feasibility check at upper bound
  pp_max <- power_surrogate(
    p_star = p_star, psi=psi,
    theta_d = theta_d, theta_c = theta_c,
    nd1 = nd1,
    n1 = n1,
    m1 = m1,
    m_d_SP = m_d_SP,
    n2 = nsecond_max,
    alpha = alpha
  )
  
  if (is.na(pp_max) || pp_max < target_power)
    return(NA)
  
  # binary search
  lo <- 0
  hi <- nsecond_max
  
  while (lo < hi) {
    mid <- floor((lo + hi) / 2)
    
    pp_mid <- power_surrogate(
      p_star = p_star, psi = psi,
      theta_d = theta_d, theta_c = theta_c,
      nd1 = nd1,
      n1 = n1,
      m1 = m1,
      m_d_SP = m_d_SP,
      n2 = mid,
      alpha = alpha
    )
    
    if (!is.na(pp_mid) && pp_mid >= target_power) {
      hi <- mid
    } else {
      lo <- mid + 1
    }
  }
  
  return(lo)
}





