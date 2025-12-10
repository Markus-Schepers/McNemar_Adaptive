
library(dplyr)

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
exact_power <- function(p,nd,alpha=0.05){
  if(nd==0) return(0)
  k1 <- lower(nd,alpha)
  k2 <- upper(nd,alpha)
  if(is.na(k1)) return(0)
  if(is.na(k2)) return(0)
  return(pbinom(k1,nd,p)+1-pbinom(k2-1,nd,p))
}

# exact power conditioned on number of discordant pairs nd1 in first stage, 
# using size of first stage and size of second stage
# using the probability for a discordant pair
conditional_power <- function(p,nd1,nsecond=38,psi=0.2,alpha=0.05){
  nd2_vals <- 0:nsecond
  sum(vapply(nd2_vals,function(nd2){ return(exact_power(p,(nd1+nd2),alpha) * dbinom(nd2,nsecond,psi)) },numeric(1))) 
}


min_nsecond <- function(p, nd1, target_power = 0.8, psi = 0.2,
                        alpha = 0.05, nsecond_max = 200) {
  # 1. Check feasibility at maximum nsecond
  cp_max <- conditional_power(p = p, nd1 = nd1,
                              nsecond = nsecond_max,
                              psi = psi, alpha = alpha)
  
  if(is.na(cp_max) || cp_max < target_power){
    return(NA)   # impossible to reach target power
  }
  
  # 2. Search for smallest n2 achieving the target
  for(n2 in 0:nsecond_max){
    cp <- conditional_power(p = p, nd1 = nd1,
                            nsecond = n2,
                            psi = psi, alpha = alpha)
    
    if(!is.na(cp) && cp >= target_power){
      return(n2)
    }
  }
  
  # Should never happen because feasibility was checked
  return(NA)
}


estimate_transition_matrix <- function(surr, final, pseudo = 0.5) {
  # surr  = factor/1..4 category for surrogate endpoint at interim
  # final = factor/1..4 category for final endpoint (same patients)
  # pseudo = pseudocount to stabilize very small samples (Laplace smoothing)
  
  if(length(surr) != length(final)) stop("Lengths differ")
  
  # Create full 4x4 table
  tab <- table(factor(final, levels=1:4), factor(surr, levels=1:4))
  T_hat <- matrix(0, nrow=4, ncol=4)
  
  for(j in 1:4) {
    # allowed i are i >= j (lower triangular)
    allowed_i <- j:4
    sub_counts <- as.numeric(tab[allowed_i, j])
    
    # smoothing to avoid 0/0 divisions
    sub_probs <- (sub_counts + pseudo) / sum(sub_counts + pseudo)
    
    T_hat[allowed_i, j] <- sub_probs
  }
  
  return(T_hat)
}


conditional_power_surrogate <- function(q_vec, Tmat, nd1_surr,
                                        nsecond = 38,
                                        alpha = 0.05) {
  
  # obtain final endpoint probabilities
  p_vec <- as.numeric(Tmat %*% q_vec)
  
  # discordant probability under final endpoint
  psi_p <- p_vec[2] + p_vec[3]
  p_star <- p_vec[2] / psi_p 
  
  cp <- conditional_power(p = p_star,
                          nd1 = nd1_surr, # recompute?
                          nsecond = nsecond,
                          psi = psi_p,
                          alpha = alpha)
  
  return(cp)
}

min_nsecond_surrogate <- function(q_vec, Tmat, nd1_surr,
                                  target_power = 0.8,
                                  alpha = 0.05,
                                  nsecond_max = 200) {
  
  p_vec <- as.numeric(Tmat %*% q_vec)
  psi_p <- p_vec[2] + p_vec[3]
  
  # feasibility check
  cp_max <- conditional_power_surrogate(q_vec, Tmat,
                                        nd1_surr,
                                        nsecond = nsecond_max,
                                        alpha = alpha)
  if(is.na(cp_max) || cp_max < target_power) return(NA)
  
  # search minimal n2
  for(n2 in 0:nsecond_max) {
    cp <- conditional_power_surrogate(q_vec, Tmat,
                                      nd1_surr,
                                      nsecond = n2,
                                      alpha = alpha)
    if(!is.na(cp) && cp >= target_power) return(n2)
  }
  
  return(NA)
}


simulate_trial <- function(n1 = 38,                # stage-1 sample size
                           n_link = 25,            # patients with both surr+final observed at interim
                           q_true,                 # true surrogate cell probs
                           T_true,                 # true transition
                           alpha = 0.05,
                           target_power = 0.8,
                           B = 1000) {
  
  results <- data.frame(
    reject = logical(B),
    n2     = integer(B)
  )
  
  for(b in 1:B) {
    
    ## ---- Stage 1 surrogate data ----
    surr_stage1 <- sample(1:4, size = n1, replace = TRUE, prob = q_true)
    nd1_surr <- sum(surr_stage1 == 2 | surr_stage1 == 3)
    
    ## ---- Sub-sample for transition estimation ----
    surr_link <- sample(surr_stage1, size = n_link, replace = FALSE)
    # simulate corresponding final endpoint via true T
    final_link <- sapply(surr_link, function(j) {
      sample(1:4, size = 1, prob = T_true[, j])
    })
    
    ## ---- Estimate T_hat ----
    T_hat <- estimate_transition_matrix(surr = surr_link,
                                        final = final_link,
                                        pseudo = 0.5)
    
    ## ---- Compute stage-2 sample size ----
    n2 <- min_nsecond_surrogate(q_vec = q_true,
                                Tmat = T_hat,
                                nd1_surr = nd1_surr,
                                target_power = target_power,
                                alpha = alpha)
    results$n2[b] <- n2
    
    ## ---- Generate remaining surrogate & final data ----
    if(is.na(n2)) {
      results$reject[b] <- FALSE
      next
    }
    
    # final sample consists of n1+n2 patients; simulate final endpoint directly
    # using true p = T_true q_true
    p_true <- as.numeric(T_true %*% q_true)
    
    final_all <- sample(1:4, size = n1 + n2, replace = TRUE, prob = p_true)
    n12 <- sum(final_all == 2)
    n21 <- sum(final_all == 3)
    test_result <- binom.test(n12, n12 + n21, p = 0.5, alternative = "two.sided")$p.value
    results$reject[b] <- (test_result <= alpha)
    
  }
  
  return(results)
}








