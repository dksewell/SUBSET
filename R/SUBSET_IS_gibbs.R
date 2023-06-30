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
true_response_rates = c(0.1,0.25,0.4)

## Simulate data
N = 50
y = rbinom(3,N,true_response_rates)

## Draws from the base prior
## (Using Jeffrey's prior on p_1, p_2, and p_3)
ndraws = 1e4

pi0_samples = 
  cbind(p1 = rbeta(ndraws,1/2,1/2),
        p2 = rbeta(ndraws,1/2,1/2),
        p3 = rbeta(ndraws,1/2,1/2)) 

## Get posterior draws corresponding to the base prior
draws0 = 
  cbind(p1 = rbeta(ndraws,1/2 + y[1], 1/2 + N - y[1]),
        p2 = rbeta(ndraws,1/2 + y[2], 1/2 + N - y[2]),
        p3 = rbeta(ndraws,1/2 + y[3], 1/2 + N - y[3])) 

## Get posterior draws corresponding to the tilted prior
library(parallel)
cl = makeCluster(10)
draws_tilted = 
  SUBSET_IS_gibbs(draws0,
                  pi0_samples,
                  P_phi = "power",
                  v_sequence = c(5,10,20,40),
                  Z_approximation = list(df = 10 + 1,
                                         phi_range = NULL),
                  cl = cl)
P_phi = "power"
v_sequence = seq(0.25,5,by = 0.25)
Z_approximation = list(df = 10 + 1,
                       phi_range = NULL)
verbose = TRUE
min_ESS = 1/2 * nrow(draws0)
n_u_values = 100
n_draws = nrow(draws0)
}

SUBSET_IS_gibbs = function(draws0,
                           pi0_samples,
                           P_phi = c("power","geometric")[1],
                           prior_phi,
                           initial_phi_to_get_u,
                           v_sequence = seq(0.25,5,by = 0.25),
                           Z_approximation = list(df = 10 + 1,
                                                  phi_range = NULL),
                           verbose = TRUE,
                           min_ESS = 1/2 * nrow(draws0),
                           n_u_values = 100,
                           n_draws,
                           cl){
  if(missing(n_draws)) n_draws = nrow(draws0)
  p = ncol(pi0_samples)
  v_len = length(v_sequence)
  acc_rate = numeric(v_len)
  prior_phi = list(d = function(x) dgamma(x,shape = 2,rate = 1),
                   r = function(n = 1) rgamma(n,shape = 2,rate = 1))
  
  # Check
  if(missing(prior_phi) & (class(P_phi) == "function")) stop("Must provide prior for custom projection function")
  
  # Create P function
  if(class(P_phi) != "function"){
    if(P_phi == "power"){
      P = function(x) Proj(cbind(1,c(1:p)^x))
      if(missing(prior_phi)) prior_phi = list(d = function(x) dgamma(x,shape = 2,rate = 2),
                                              r = function(n = 1) rgamma(n,shape = 2,rate = 2))
    }
    if(P_phi == "geometric"){
      P = function(x) Proj(cbind(1,x^(-c(1:p))))
      if(missing(prior_phi)) prior_phi = function(x) list(d = function(x) dbeta(x,2,2), #Keeping the d in case I later allow a different proposal than the prior
                                                          r = function(n = 1) rbeta(n,2,2))
    }
  }else{
    P = P_phi
  }
  
  if(missing(initial_phi_to_get_u)) initial_phi_to_get_u = mean(prior_phi$r(1e3))
  
  cat("\nGetting u for IS proposal distributions\n")
  fixed = 
    SUBSET_IS_fixed(draws0,P(initial_phi_to_get_u),
                    verbose = FALSE)
  u_values = 
    fixed$summary$u[which( (fixed$summary$v > 0) &
                             (fixed$summary$variable == fixed$summary$variable[1]) )]
  
  
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
  samples[1,,1] = c(colMeans(draws0),initial_phi_to_get_u)
  
  
  # Get exact (actually, this is MC) method for computing Z_{\nu,\phi}
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
    clusterExport(cl,c("pi0_samples","v_sequence","draws0","u_values","p"),envir = environment())
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
  
  if(class(Z_approximation) == "list"){
    if(is.null(Z_approximation$phi_range)){
      phi_prior_draws = prior_phi$r(1e4)
      Z_approximation$phi_range = range(phi_prior_draws)
    }
    phi_seq = seq(Z_approximation$phi_range[1],
                  Z_approximation$phi_range[2],
                  l = ifelse(is.null(Z_approximation$seq_length),
                             2*Z_approximation$df,
                             Z_approximation$seq_length))
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
    
    if(class(Z_approximation) == "list"){
      library(splines)
      if(!missing(cl))clusterEvalQ(cl,{library(splines);library(SUBSET)})
      cat("\nFitting spline interpolation for Z(phi)\n")
      # Fit Z to a sequence of values of phi
      Z_vals = 
        sapply(phi_seq,function(phi){
          Z_exact(diag(p) - P(phi),v_sequence[v])
        })
      ns_fit = 
        lm(Z ~ ns(phi,df = Z_approximation$df),
           data = data.frame(Z = Z_vals,
                             phi = phi_seq))
      if(verbose){
        par(mar = c(5,5,3,0.5))
        plot(fitted(ns_fit) ~ phi_seq, 
             type = 'l', 
             bty = 'l',
             lwd = 2,
             xlab = expression(phi),
             ylab = expression(Z[phi]),
             main = paste0("Spline interpolation of Z for v = ",v_sequence[v]),
             cex.axis = 1.5,
             cex.lab = 1.5,
             cex.main = 1.5)
        points(Z_vals ~ phi_seq, col = adjustcolor("tomato",0.5))
      }
      
      Z = function(x){
        predict(ns_fit,
                newdata = data.frame(phi = x))
      }
    }else{
      Z = function(phi){
        Z_exact(diag(p) - P(phi), v_sequence[v])
      }
    }
    
    
    if(v > 1) samples[1,,v] = samples[1,,v - 1]
    P_old = P(samples[1,p+1,v])
    Z_old = Z(samples[1,p + 1,v])
    cat("\nPerforming Gibbs Sampling\n")
    if(verbose) pb = txtProgressBar(0,n_draws,style=3)
    for(it in 2:n_draws){
      
      # Draw phi from prior, adjust with MH
      phi_proposal = prior_phi$r(1)
      P_new = P(phi_proposal)
      # P_orth_new = diag(p) - P_new
      Z_new = Z(phi_proposal)
      acc_prob = 
        exp(-0.5 * v_sequence[v] * 
              drop(crossprod(samples[it - 1,1:p,v],(P_old - P_new) %*% samples[it - 1,1:p,v])) ) * 
        Z_old / Z_new # No prior terms since this cancels in the acceptance prob., i.e., g = pi
      
      if(runif(1) < acc_prob){
        samples[it,p+1,v] = phi_proposal
        P_old = P_new
        # P_orth_old = P_orth_new
        Z_old = Z_new
        
        acc_rate[v] = acc_rate[v] + 1 / (n_draws - 1)
      }else{
        samples[it,p+1,v] = samples[it - 1,p+1,v]
      }
      
      # Draw theta via IS
      if(missing(cl)){
        w = 
          sapply(1:ndraws,
                 function(i){
                   exp(-0.5 * v_sequence[v] * u_values[v]^2 *
                         tcrossprod(draws0[i,],draws0[i,] %*% (diag(p) - P_old)))
                 })
      }else{
        clusterExport(cl,c("v","P_old"),envir = environment())
        w = 
          parSapply(cl,
                    1:ndraws,
                    function(i){
                      exp(-0.5 * v_sequence[v] * u_values[v]^2 *
                            tcrossprod(draws0[i,],draws0[i,] %*% (diag(p) - P_old)))
                      })
      }
      samples[it,1:p,v] = 
        draws0[sample(nrow(draws0),1,prob = w),]
      
      
      if(verbose) setTxtProgressBar(pb,it)
    }
    
    
  }
  
  
  
  
  
  
}





