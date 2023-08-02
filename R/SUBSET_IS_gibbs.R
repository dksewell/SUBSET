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
#' from the *base prior* (not including phi)
#' @param P_phi Either (1) a function taking in a scalar and returning 
#' a projection matrix for the model parameters; or (2) "power" 
#' in which case the projection matrix will correspond to 
#' span((1,1,...)',(1^x,2^x,...)'); or (3) "geometric" in 
#' which case the projection matrix will correspond to 
#' span((1,1,...)',(x^1,x^2,...)').
#' @param prior_phi list with named arguments "q" and "d".  prior_phi$q() 
#' should provide the quantile function, and prior_phi$d() should provide 
#' the density function.
#' @param v_sequence vector of non-negative values for exponential shrinkage.  Higher 
#' implies more shrinkage.
#' @param min_ESS integer.  Minimum number of effective sample size to get 
#' u value. WARNING:
#' Don't make this too high, or you risk not representing the full parameter space well.
#' @param n_u_values Number of values to try for u, the hyperparameter of the  
#' proposal distribution that shrinks towards the linear subspace
#' @param verbose logical. Should any output be provided? 
#' @param cl parallel socket cluster (see parallel::makeCluster()).
#' 
#' @return Object of class "subset_gibbs", "subset_is_gibbs", with the following structure:
#' \itemize{
#'    \item samples - 3-dim array.  First index is the samples, 
#' second index is the variable, and third index is the
#' value of v (the exponential tilting parameter)
#'    \item acc_rate: numeric giving the acceptance rate for phi
#'    \item draws0
#'    \item pi0_samples
#'    \item P_phi
#'    \item prior_phi
#'    \item phi_pmf - probability mass function for \phi
#'    \item v_sequence
#' }
#' 
#' @export 

SUBSET_IS_gibbs = function(draws0,
                           pi0_samples,
                           P_phi = c("power","geometric")[1],
                           prior_phi,
                           phi_sequence = 10, #either seq or integer for length of seq
                           v_sequence = seq(0.25,5,by = 0.25),
                           min_ESS = 1/2 * nrow(draws0),
                           n_u_values = 100,
                           n_draws = nrow(draws0),
                           verbose = TRUE,
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
      if(missing(prior_phi)) prior_phi = list(d = function(x) dbeta(x,2,2),
                                              q = function(p) qbeta(p,2,2),
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
  cat("\nGetting normalizing constants for each phi and each v\n")
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
             prior_phi = prior_phi,
             phi_pmf = cbind(phi_sequence,phi_pmf),
             v_sequence = v_sequence)
  class(ret) = c("subset_gibbs","subset_is_gibbs")
  return(ret)
}





