#' SUBspace Shrinkage via Exponential Tilting for estimated 
#' projection matrix
#' 
#' SUBSET_gibbs() will compute perform a 2-block Gibbs sampler, 
#' where the blocks correspond to the primary model parameters 
#' and the projection matrix parameter (only one parameter 
#' supported at this time).  The Gibbs algorithm relies on 
#' the asymptotic normal approximation to the conditional posterior 
#' using an exponential tilted prior that shrinks the 
#' parameters towards a linear subspace for a given projection 
#' matrix.
#' 
#' SUBSET takes a base prior \eqn{\pi_0} and tilts it towards
#' a linear subspace via the multiplicative term 
#' \deqn{e^{-\frac{v}{2}\theta'(I-P)\theta}.}
#' The user must perform the initial Bayesian inference for 
#' this base prior, and compute the posterior mean and covariance. 
#' *It is imperative that the necessary conditions for an asymptotic 
#' normal approximation of the posterior hold for the base prior.* 
#' For a projection matrix \eqn{P} parameterized by \eqn{\phi}, 
#' SUBSET_gibbs() will obtain joint draws of the model parameters and 
#' \eqn{\phi}. \code{v_sequence} is a user provided positive 
#' number or vector of positive numbers that control the shrinkage: 
#' higher values imply more shrinkage of the posterior while \eqn{v=0} 
#' corresponds to the base prior with no shrinkage.
#' 
#' @param Mean0 vector. The posterior mean corresponding
#' to the base prior (i.e., untilted prior).
#' @param Sigma0 matrix. The posterior covariance 
#' matrix corresponding to the base prior. NOTE: The
#' variable names are extracted from Sigma0, so use 
#' colnames(Sigma0) = ...
#' @param n_draws positive integer of how many posterior 
#' samples are desired.
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
#' @param cl (optional) cluster from the parallel package.
#' 
#' @return Object of class "subset_gibbs".  
#' * samples: 3-dim array.  First index is the samples, 
#' second index is the variable, and third index is the
#' value of v (the exponential tilting parameter)
#' 
#' * acc_rate: numeric giving the acceptance rate for phi
#' 
#' * Mean0: the posterior mean without tilting
#' 
#' * Sigma0: the posterior covariance without tilting
#' 
#' * pi0_samples: Samples of the model parameters 
#' from the *prior* (not including phi)
#' 
#' * P_phi: a function taking in a scalar and returning 
#' a projection matrix for the model parameters
#' 
#' * prior_phi: list with named argument "r".
#' 
#' * v_sequence: vector of exponential tilting parameters.


SUBSET_gibbs = function(Mean0,
                        Sigma0,
                        n_draws,
                        pi0_samples,
                        P_phi = c("power","geometric")[1],
                        prior_phi,
                        v_sequence = seq(0.25,5,by = 0.25),
                        Z_approximation = list(df = 10 + 1,
                                               phi_range = NULL),
                        verbose = TRUE,
                        cl){
  
  p = ncol(pi0_samples)
  v_len = length(v_sequence)
  acc_rate = numeric(v_len)
  
  Sigma0_inv = chol2inv(chol(Sigma0))
  
  
  if(missing(n_draws)) n_draws = nrow(pi0_samples)
  
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
  
  # Create array to store samples
  if(is.null(colnames(Sigma0))){
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
                            c(colnames(Sigma0),"phi"),
                            v_sequence))
  }
  samples[1,,1] = c(Mean0,prior_phi$r())
  
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
    clusterExport(cl,"pi0_samples",envir = environment())
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
      phi_prior_draws = prior_phi$r(1e3)
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
      
      # Draw phi from prior
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
      
      # Draw theta
      Sigma_inv = Sigma0_inv + v * (diag(p) - P_old)
      Sigma =  chol2inv(chol(Sigma_inv))
      Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
      samples[it,1:p,v] = 
        drop(
          rmvnorm(1,
                  mean = Mean,
                  sigma = Sigma)
        )
      
      if(verbose) setTxtProgressBar(pb,it)
    }
    
    
  }
  
  ret = list(samples = samples,
             acc_rate = acc_rate,
             Mean0 = Mean0,
             Sigma0 = Sigma0,
             pi0_samples = pi0_samples,
             P_phi = P,
             prior_phi = prior_phi,
             v_sequence = v_sequence)
  class(ret) = "subset_gibbs"
  return(ret)
}
