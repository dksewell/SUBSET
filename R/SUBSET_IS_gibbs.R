#' SUBspace Shrinkage via Exponential Tilting for fixed 
#' projection matrix
#' 
#' SUBSET_IS_gibbs() will obtain posterior inference 
#' using an exponential tilted prior that shrinks the 
#' parameters towards a linear subspace.
#' 
#' SUBSET_IS_gibbs() will compute perform a 2-block Gibbs sampler, 
#' where the blocks correspond to the primary model parameters 
#' and the projection matrix parameter (only one parameter 
#' supported at this time).  The Gibbs algorithm performs importance 
#' sampling to obtain draws from the conditional posterior of theta 
#' using an exponential tilted prior that shrinks the 
#' parameters towards a linear subspace for a given projection 
#' matrix.
#' *Note that that the parameter space must be convex.* (It usually is.)
#' 
#' 
#' @param draws0 matrix. The posterior draws  
#' corresponding to the base prior.
#' @param pi0_samples Samples of the model parameters 
#' from the * base prior* (not including phi)
#' @param P_phi Either (1) a function taking in a scalar and returning 
#' a projection matrix for the model parameters; or (2) "power" 
#' in which case the projection matrix will correspond to 
#' span((1,1,...)',(1^x,2^x,...)'); or (3) "geometric" in 
#' which case the projection matrix will correspond to 
#' span((1,1,...)',(x^1,x^2,...)').
#' @param prior_phi list with named argument "r".  prior_phi$r(n) 
#' should provide n random draws from the prior on phi.
#' @param initial_phi_to_get_u best guess for phi to get a u 
#' (the hyperparameter that  shrinks the proposal distribution 
#' towards the linear subspace) with which to use during the 
#' importance sampling step of the 2-block Gibbs sampler
#' @param v_sequence vector of positive values. See details.
#' @param Z_approximation If set to be NULL, Z will be 
#' evaluated "exactly" using the pi0_samples.  Otherwise a list 
#' must be provided in order to approximate Z(phi) via 
#' natural cubic splines.  If the latter, the list must 
#' have named arguments df, the degrees of freedom (see 
#' ns()), and phi_range, a vector of length 2 giving the 
#' lower and upper bound of \eqn{\phi}.  phi_range is 
#' optional, and if left missing will be obtained from 
#' the range of prior_phi(1000).
#' @param verbose logical. Should any output be provided? 
#' @param min_ESS integer.  Minimum number of effective sample size to get 
#' u value (based on initial_phi_to_get_u). WARNING:
#' Don't make this too high, or you risk not representing the full parameter space well.
#' @param n_u_values Number of values to try for u, the hyperparameter that 
#' shrinks the proposal distribution towards the linear subspace (based on 
#' initial_phi_to_get_u).
#' 
#' @examples 

if(FALSE){
  library(SUBSET)
  source("~/SUBSET/R/SUBSET_IS_fixed.R")
# Toy example, with response rate as an increasing function of dose

set.seed(2023)

## Set true response rates
true_response_rates = 
  c(0.05,0.1,0.15,0.2,0.3,0.4,0.55,0.7)

## Simulate data
N = 50
y = rbinom(length(true_response_rates),N,true_response_rates)

## Draws from the base prior
## (Using Jeffrey's prior on p_1, p_2, and p_3)
ndraws = 1e4

pi0_samples = 
  matrix(rbeta(ndraws*length(true_response_rates),1/2,1/2),
         ndraws,
         length(true_response_rates))

## Get posterior draws corresponding to the base prior
draws0 = NULL
for(j in 1:length(true_response_rates)){
  draws0 = 
    cbind(draws0,
          rbeta(ndraws,1/2 + y[j], 1/2 + N - y[j])
    )
}
colnames(draws0) = paste("p",1:length(true_response_rates),sep="")

## Get posterior draws corresponding to the tilted prior
library(parallel)
cl = makeCluster(10)
draws_tilted = 
  SUBSET_IS_gibbs(draws0,
                  pi0_samples,
                  P_phi = "power",
                  prior_phi = list(d = function(x) dgamma(x,shape = 2,rate = 1),
                                   q = function(p) qgamma(p,shape = 2,rate = 1),
                                   r = function(n = 1) rgamma(n,shape = 2,rate = 1)),
                  phi_sequence = 10,
                  v_sequence = c(5,10,20,40),
                  cl = cl)
stopCluster(cl)

library(tidyverse);library(magrittr)
phi_distribution = 
  draws_tilted$phi_pmf %>% 
  as_tibble() %>% 
  mutate(phi_sequence = round(phi_sequence,3))
for(v in 1:length(draws_tilted$v_sequence)){
  temp = 
    prop.table(table(draws_tilted$samples[,length(true_response_rates) + 1,v]))
  temp = 
    tibble(phi_sequence = as.numeric(names(temp)),
           v1 = as.numeric(temp)) %>% 
    mutate(phi_sequence = round(phi_sequence,3))
  phi_distribution = 
    left_join(phi_distribution,
              temp,
              by = "phi_sequence")
}
names(phi_distribution) = 
  c("phi_sequence","prior",paste("v",as.character(draws_tilted$v_sequence),sep=""))
plot(v40 ~ phi_sequence,
     data = phi_distribution,
     type = 'l',
     lwd = 2,
     col = gray(seq(0.2,0.8,l=ncol(phi_distribution))[ncol(phi_distribution)]))
for(j in 1:(ncol(phi_distribution)-1)){
  lines(unlist(phi_distribution[,j]) ~ phi_distribution$phi_sequence,
       lwd = 2,
       col = gray(seq(0.2,0.8,l=ncol(phi_distribution))[j]))
}

P_phi = "power"
prior_phi = list(d = function(x) dgamma(x,shape = 2,rate = 1),
                 q = function(p) qgamma(p,shape = 2,rate = 1),
                 r = function(n = 1) rgamma(n,shape = 2,rate = 1))
phi_sequence = 10
v_sequence = c(5,10,20,40)
verbose = TRUE
min_ESS = 1/2 * nrow(draws0)
n_u_values = 100
n_draws = nrow(draws0)
}

SUBSET_IS_gibbs = function(draws0,
                           pi0_samples,
                           P_phi = c("power","geometric")[1],
                           prior_phi,
                           phi_sequence = 10, #either seq or integer for length of seq
                           v_sequence = seq(0.25,5,by = 0.25),
                           verbose = TRUE,
                           min_ESS = 1/2 * nrow(draws0),
                           n_u_values = 100,
                           n_draws = nrow(draws0),
                           cl){
  
  p = ncol(pi0_samples)
  v_len = length(v_sequence)
  acc_rate = numeric(v_len)
  
  # Check
  if(missing(prior_phi) & (class(P_phi) == "function")) stop("Must provide prior for custom projection function")
  
  # Create P function
  if(class(P_phi) != "function"){
    if(P_phi == "power"){
      P = function(x) Proj(cbind(1,c(1:p)^x))
      if(missing(prior_phi)) prior_phi = list(d = function(x) dgamma(x,shape = 2,rate = 2),
                                              q = function(p) qgamma(p,shape = 2,rate = 2),
                                              r = function(n = 1) rgamma(n,shape = 2,rate = 2))
    }
    if(P_phi == "geometric"){
      P = function(x) Proj(cbind(1,x^(-c(1:p))))
      if(missing(prior_phi)) prior_phi = function(x) list(d = function(x) dbeta(x,2,2),
                                                          q = function(p) qgamma(p,shape = 2,rate = 2),
                                                          r = function(n = 1) rbeta(n,2,2))
    }
  }else{
    P = P_phi
  }
  
  # Get pmf for phi_sequence
  if(length(phi_sequence) == 1){
    phi_sequence = 
      prior_phi$q(seq(0,1,l=phi_sequence + 2)[-c(1,phi_sequence+2)])
  }
  phi_pmf = prior_phi$d(phi_sequence)
  phi_pmf = phi_pmf / sum(phi_pmf)
  
  
  cat("\nGetting weights for IS proposal distributions\n")
  fixed = list()
  for(i in 1:length(phi_sequence)){
    fixed[[i]] = 
      SUBSET_IS_fixed(draws0,
                      P(phi_sequence[i]),
                      v_sequence = v_sequence,
                      min_ESS = min_ESS,
                      n_u_values = n_u_values,
                      verbose = FALSE)
  }
  names(fixed) = as.character(phi_sequence)
  
  # Create array to store samples
  if(is.null(colnames(draws0))){
    samples = 
      array(0.0,
            c(n_draws,p + 1,v_len),
            dimnames = list(NULL,
                            c(paste("v",1:p,sep=""),"phi"),
                            v_sequence))
  }else{
    samples = 
      array(0.0,
            c(n_draws,p + 1,v_len),
            dimnames = list(NULL,
                            c(colnames(draws0),"phi"),
                            v_sequence))
  }
  
  
  # Get exact (actually, this is MC) values of Z_{\nu,\phi}
  if(missing(cl)){
    Z_exact = function(I_m_P,V){
      mean(
        sapply(1:nrow(pi0_samples), 
               function(i){
                 exp(-0.5 * V * drop(crossprod(pi0_samples[i,],I_m_P %*% pi0_samples[i,])))
               })
      )
    }
    
  }else{
    clusterExport(cl,c("pi0_samples"),envir = environment())
    Z_exact = function(I_m_P,V){
      clusterExport(cl,c("I_m_P","V"),envir = environment())
      
      mean(
        parSapply(cl,
                  1:nrow(pi0_samples), 
                  function(i){
                    exp(-0.5 * V * drop(crossprod(pi0_samples[i,],I_m_P %*% pi0_samples[i,])))
                  })
      )
    }
  }
  cat("\nGetting normalizing constances for each phi and each v\n")
  Z_phi = matrix(0.0,length(phi_sequence),length(v_sequence),
                 dimnames = list(as.character(phi_sequence),
                                 as.character(v_sequence)))
  for(x in 1:nrow(Z_phi)){
    for(y in 1:ncol(Z_phi)){
      Z_phi[x,y] = 
        Z_exact(diag(p) - P(phi_sequence[x]),v_sequence[y])
    }
  }
  
  
  #Perform Gibbs sampler for each v 
  for(v in 1:v_len){
    if(verbose){
      cat("\n")
      cat(paste(rep("-",4 + nchar(as.character(v_sequence[v]))),collapse=""))
      cat("\n")
      cat(paste0("v = ",v_sequence[v]))
      cat("\n")
      cat(paste(rep("-",4 + nchar(as.character(v_sequence[v]))),collapse=""))
      cat("\n")
    }
    
    # Get initial draw
    ## ...of phi
    samples[1,p + 1,v] = sample(phi_sequence,1,prob = phi_pmf)
    ## ...of theta
    new_draw_index = 
      sample(1:nrow(draws0),1,
             prob = fixed[[as.character(samples[1,p + 1,v])]]$is_draws[[v]][,p+1])
    samples[1,1:p,v] = 
      fixed[[as.character(samples[1,p + 1,v])]]$is_draws[[v]][new_draw_index,1:p]
    
    # Set up necessary objects for gibbs sampler
    P_old = P(samples[1,p+1,v])
    Z_old = Z_phi[as.character(samples[1,p + 1,v]),v]
    
    cat("\nPerforming Gibbs Sampling\n")
    if(verbose) pb = txtProgressBar(0,n_draws,style=3)
    for(it in 2:n_draws){
      
      # Draw phi from prior, adjust with MH
      phi_proposal = sample(phi_sequence,1,prob = phi_pmf)
      P_new = P(phi_proposal)
      Z_new = Z_phi[as.character(phi_proposal),v]
      acc_prob = 
        exp(-0.5 * v_sequence[v] * 
              drop(crossprod(samples[it - 1,1:p,v],(P_old - P_new) %*% samples[it - 1,1:p,v])) ) * 
        Z_old / Z_new # No prior terms since this cancels in the acceptance prob., i.e., g = pi
      
      if(runif(1) < acc_prob){
        samples[it,p+1,v] = phi_proposal
        P_old = P_new
        Z_old = Z_new
        acc_rate[v] = acc_rate[v] + 1 / (n_draws - 1)
      }else{
        samples[it,p+1,v] = samples[it - 1,p+1,v]
      }
      
      # Draw theta via IS
      new_draw_index = 
        sample(1:nrow(draws0),1,
               prob = fixed[[as.character(samples[it,p + 1,v])]]$is_draws[[v]][,p+1])
      samples[it,1:p,v] = 
        fixed[[as.character(samples[it,p + 1,v])]]$is_draws[[v]][new_draw_index,1:p]
      
      
      if(verbose) setTxtProgressBar(pb,it)
    }
    
    
  }
  
  
  
  ret = list(samples = samples,
             acc_rate = acc_rate,
             draws0 = draws0,
             pi0_samples = pi0_samples,
             P_phi = P,
             phi_pmf = cbind(phi_sequence,phi_pmf),
             v_sequence = v_sequence)
  class(ret) = "subset_is_gibbs"
  return(ret)
}





