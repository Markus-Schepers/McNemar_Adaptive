
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





estimate_theta <- function(n_disc_surr,
                           n_disc_final,
                           n_disc_surr_final,
                           n_total,
                           pseudo = 0.5) {
  
  # 2x2 table
  # rows: surrogate (1=disc, 0=conc)
  # cols: final     (1=disc, 0=conc)
  
  n11 <- n_disc_surr_final
  n10 <- n_disc_surr - n_disc_surr_final
  n01 <- n_disc_final - n_disc_surr_final
  n00 <- n_total - n_disc_surr - n_disc_final + n_disc_surr_final
  
  theta1 <- (n11 + pseudo) / (n11 + n10 + 2*pseudo)
  theta0 <- (n01 + pseudo) / (n01 + n00 + 2*pseudo)
  
  return(c(theta1 = theta1, theta0 = theta0))
}
predict_nd1_from_surrogate <- function(nd1_surr, n1, theta1, theta0) {
  round(theta1 * nd1_surr + theta0 * (n1 - nd1_surr))
}


# when using the surrogate endpoint, re-estimate the sample size by estimating the number of discordant pairs in the final endpoint at interim analysis, and use the conditional power for the final endpoint
# p_star is the effect size (once the fraction of discordant pairs psi is fixed)
predicted_power_surrogate <- function(p_star,psi,
                                      nd1_surr,
                                      n1,
                                      theta1,
                                      theta0,
                                      nsecond,
                                      alpha) {
  
  nd1 <- predict_nd1_from_surrogate(nd1_surr, n1, theta1, theta0)
  
  return(conditional_power(
    p_star = p_star, psi = psi,
    nd1 = nd1,
    nsecond = nsecond,
    alpha = alpha
  ))
}


min_nsecond_surrogate <- function(
    p_star, psi,
    nd1_surr,
    n1,
    theta1,
    theta0,
    target_power = 0.8,
    alpha = 0.05,
    nsecond_max = 200) {
  # feasibility check
  pp_max <- predicted_power_surrogate(p_star,nd1_surr,n1,theta1 = theta1,theta0 = theta0, psi = psi,nsecond=nsecond_max,alpha = alpha)
  if(is.na(pp_max) || pp_max < target_power) return(NA)
  
  # search minimal n2
  for(n2 in 0:nsecond_max) {
    pp <- predicted_power_surrogate(
      p_star = p_star,
      nd1_surr = nd1_surr,
      n1 = n1,
      theta1 = theta1,theta0 = theta0,psi=psi,
      nsecond=n2,
      alpha = alpha)
    if(!is.na(pp) && pp >= target_power) return(n2)
  }
  
  return(NA)
}






