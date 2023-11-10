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
#' @param prior_draws matrix of prior draws (each row is a draw), 
#' corresponding to the base prior.
#' @param n_draws positive integer of how many posterior 
#' samples are desired.
#' @param P_phi Either (1) a function taking in a scalar and returning 
#' a projection matrix for the model parameters; or (2) "power" 
#' in which case the projection matrix will correspond to 
#' span((1,1,...)',(1^x,2^x,...)'); or (3) "geometric" in 
#' which case the projection matrix will correspond to 
#' span((1,1,...)',(x^1,x^2,...)').
#' @param prior_phi list with named argument "r".  prior_phi$r(n) 
#' should provide n random draws from the prior on phi.
#' @param nu shrinkage parameter.  If missing, the optimal nu via Bayes factors will be used
#' @param nu_max maximum value of nu to be considered.
#' @param Z_approximation If set to be NULL, Z will be 
#' evaluated "exactly" using the pi0_samples.  Otherwise a list 
#' must be provided in order to approximate Z(phi) via 
#' natural cubic splines.  If the latter, the list must 
#' have named arguments df, the degrees of freedom (see 
#' ns()); phi_range, a vector of length 2 giving the 
#' lower and upper bound of \eqn{\phi}; and optionally an integer named 
#' seq_length giving the number of training points along phi_range to
#' train the spline regression model (if missing, the default is twice the 
#' degrees of freedom).  phi_range is also
#' optional, and if left missing will be obtained from 
#' the range of prior_phi(1000).
#' @param verbose logical. Should any output be provided?
#' @param cl (optional) cluster from the parallel package.
#' 
#' @return Object of class "subset_gibbs", "subset_asymp_gibbs", 
#' with the following structure:
#' \itemize{
#'    \item samples - matrix of posterior draws from the 2-block Gibbs sampler. Each 
#'    row is a posterior draw.  The last column corresponds to phi.
#'    \item acc_rate: numeric giving the acceptance rate for phi
#'    \item Mean0
#'    \item Sigma0
#'    \item P_phi
#'    \item prior_phi
#' }
#' 
#' 
#' @import parallel
#' @import Matrix
#' @import graphics
#' @import mvtnorm
#' @import splines
#' @import stats
#' @import utils
#' @importFrom grDevices adjustcolor
#' @export
#' @exportClass subset_gibbs
#' @exportClass subset_asymp_gibbs


SUBSET_gibbs = function(Mean0,
                        Sigma0,
                        prior_draws,
                        n_draws,
                        P_phi = c("power","geometric")[1],
                        prior_phi,
                        nu,
                        nu_max,
                        Z_approximation = list(df = 10 + 1,
                                               phi_range = NULL),
                        verbose = TRUE,
                        cl){
  
  p = ncol(prior_draws)
  acc_rate = numeric(1)
  Mean0 = drop(Mean0)
  
  Sigma0_inv = chol2inv(chol(Sigma0))
  
  
  if(missing(n_draws)) n_draws = nrow(prior_draws)
  
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
      if(missing(prior_phi)) prior_phi = list(d = function(x) dbeta(x,2,2),
                                              r = function(n = 1) rbeta(n,2,2))
    }
  }else{
    P = P_phi
  }
  
  # Create array to store samples
  if(is.null(colnames(Sigma0))){
    samples = 
      matrix(0.0,
             n_draws,p + 1,
             dimnames = list(NULL,
                             c(paste("theta",1:p,sep=""),"phi")))
  }else{
    samples = 
      matrix(0.0,
             n_draws,p + 1,
             dimnames = list(NULL,
                             c(colnames(Sigma0),"phi")))
  }
  
  
  samples[1,] = c(Mean0,mean(prior_phi$r(250)))
  
  
  # Get nu
  Sigma0_logdet = determinant(Sigma0)$modulus
  if(missing(nu)){
    Proj_mat = P(samples[1,p+1])
    I_m_P = diag(p) - Proj_mat
    
    if(missing(cl)){
      wtilde_1 = 
        sapply(1:nrow(prior_draws),
               function(i){
                 exp(-0.5 * tcrossprod(prior_draws[i,],prior_draws[i,] %*% I_m_P))
               })
    }else{
      clusterExport(cl,c("prior_draws","I_m_P"),envir = environment())
      wtilde_1 = 
        parSapply(cl,
                  1:ndraws,
                  function(i){
                    exp(-0.5 * tcrossprod(prior_draws[i,],prior_draws[i,] %*% I_m_P))
                  })
    }
    non_nu_term = 
      -0.5 * Sigma0_logdet -
      0.5 * drop( crossprod(Mean0,Sigma0_inv %*% Mean0) )
    
    if(missing(cl)){
      get_BF = function(nu){
        Z_nu = mean(sapply(wtilde_1,function(x)x^nu))
        
        Sigma_inv = Sigma0_inv + nu * I_m_P
        Sigma =  chol2inv(chol(Sigma_inv))
        Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
        
        non_nu_term -
          0.5 * determinant(Sigma_inv)$modulus +
          0.5 * drop( crossprod(Mean,Sigma_inv %*% Mean) ) - 
          log(Z_nu)
      }
    }else{
      clusterExport(cl,
                    c("wtilde_1","Sigma0_inv","Mean0","non_nu_term"),
                    envir = environment())
      get_BF = function(nu){
        Z_nu = mean(parSapply(cl,wtilde_1,function(x)x^nu))
        
        Sigma_inv = Sigma0_inv + nu * I_m_P
        Sigma =  chol2inv(chol(Sigma_inv))
        Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
        
        non_nu_term -
          0.5 * determinant(Sigma_inv)$modulus +
          0.5 * drop( crossprod(Mean,Sigma_inv %*% Mean) ) - 
          log(Z_nu)
      }
    }
    lower_bound = 0
    nu = upper_bound = 5
    safety = 0
    while( ( abs(nu - upper_bound) / upper_bound < 1e-3) & (safety < 25)){
      upper_bound = 2 * upper_bound
      opt = optimize(get_BF,
                     interval = c(lower_bound,upper_bound),
                     maximum = TRUE)
      nu = opt$maximum
    }
  }
    
  
  
  # Get exact (actually, this is MC) method for computing Z_{\nu,\phi}
  if(missing(cl)){
    Z_exact = function(I_m_P){
      mean(
        sapply(1:nrow(prior_draws), 
               function(i){
                 exp(-0.5 * nu * drop(crossprod(prior_draws[i,],I_m_P %*% prior_draws[i,])))
               })
      )
    }
    
  }else{
    clusterExport(cl,c("prior_draws","nu"),envir = environment())
    Z_exact = function(I_m_P){
      clusterExport(cl,c("I_m_P"),envir = environment())
      
      mean(
        parSapply(cl,
                  1:nrow(prior_draws), 
                  function(i){
                    exp(-0.5 * nu * drop(crossprod(prior_draws[i,],I_m_P %*% prior_draws[i,])))
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
  if(class(Z_approximation) == "list"){
    cat("\nFitting spline interpolation for Z(phi)\n")
    # Fit Z to a sequence of values of phi
    Z_vals = 
      sapply(phi_seq,function(phi){
        Z_exact(diag(p) - P(phi))
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
           main = paste0("Spline interpolation of Z for v = ",nu),
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
  
  
  P_old = P(samples[1,p+1])
  Z_old = Z(samples[1,p + 1])
  cat("\nPerforming Gibbs Sampling\n")
  if(verbose) pb = txtProgressBar(0,n_draws,style=3)
  for(it in 2:n_draws){
    
    # Draw phi from prior
    phi_proposal = prior_phi$r(1)
    P_new = P(phi_proposal)
    # P_orth_new = diag(p) - P_new
    Z_new = Z(phi_proposal)
    acc_prob = 
      exp(-0.5 * nu * 
            drop(crossprod(samples[it - 1,1:p],(P_old - P_new) %*% samples[it - 1,1:p])) ) * 
      Z_old / Z_new # No prior terms since this cancels in the acceptance prob., i.e., g = pi
    
    if(runif(1) < acc_prob){
      samples[it,p+1] = phi_proposal
      P_old = P_new
      # P_orth_old = P_orth_new
      Z_old = Z_new
      
      acc_rate = acc_rate + 1 / (n_draws - 1)
    }else{
      samples[it,p+1] = samples[it - 1,p+1]
    }
    
    # Draw theta
    Sigma_inv = Sigma0_inv + nu * (diag(p) - P_old)
    Sigma =  chol2inv(chol(Sigma_inv))
    Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
    samples[it,1:p] = 
      drop(
        rmvnorm(1,
                mean = Mean,
                sigma = Sigma)
      )
    
    if(verbose) setTxtProgressBar(pb,it)
  }
  
  
  ret = list(samples = samples,
             acc_rate = acc_rate,
             Mean0 = Mean0,
             Sigma0 = Sigma0,
             P_phi = P,
             prior_phi = prior_phi,
             nu = nu,
             prior_draws = prior_draws)
  class(ret) = c("subset_gibbs", "subset_asymp_gibbs")
  return(ret)
}
