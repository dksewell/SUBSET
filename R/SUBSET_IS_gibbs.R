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
#' @param draws0 named list with names 'posterior' and 'prior'. 
#' Each ought to be a matrix of posterior or prior draws (each row is a draw),
#' all of which should correspond to the base prior.
#' @param P_phi Either (1) a function taking in a scalar and returning 
#' a projection matrix for the model parameters; or (2) "power" 
#' in which case the projection matrix will correspond to 
#' span((1,1,...)',(1^x,2^x,...)'); or (3) "geometric" in 
#' which case the projection matrix will correspond to 
#' span((1,1,...)',(x^1,x^2,...)').
#' @param prior_phi list with named arguments "q" and "d".  prior_phi$q() 
#' should provide the quantile function, and prior_phi$d() should provide 
#' the density function.
#' @param n_draws integer giving the number of 2-block Gibbs samples desired
#' @param phi_sequence Either integer giving the number of evenly spaced quantiles to 
#' discretize the space of phi (using prior_phi$q), or else a vector of values of phi
#' @param nu shrinkage parameter.  If missing, the optimal nu via Bayes factors will be used
#' @param nu_max maximum value of nu to be considered.
#' @param min_ESS integer.  Minimum number of effective sample size in considering values of nu. 
#' @param verbose logical. Should any output be provided? 
#' @param cl parallel socket cluster (see parallel::makeCluster()).
#' 
#' @return Object of class "subset_gibbs", "subset_is_gibbs", with the following structure:
#' \itemize{
#'    \item samples - matrix of posterior draws from the 2-block Gibbs sampler. Each 
#'    row is a posterior draw.  The last column corresponds to phi.
#'    \item acc_rate - numeric giving the acceptance rate for phi
#'    \item draws0 - original posterior and prior draws (under base prior)
#'    \item P_phi
#'    \item prior_phi
#'    \item phi_pmf - matrix with columns giving the values of phi, the prior probability, 
#'    and the posterior probability
#'    \item nu - either the value of nu supplied by user or the value selected by maximizing 
#'    the Bayes factor where phi is set to the prior mode.
#'    \item IS_ESS - the effective sample sizes for the importance samplers for each value of phi
#' }
#' 
#' @examples
#' # Simulate where truth lies far from subspace
#' set.seed(2023)
#' theta_true = c(1,1.5)
#' N = 100
#' ndraws = 1e4
#' draws0 = list()
#' # y ~ N(mu, I_2)
#' y = mvtnorm::rmvnorm(N,theta_true,diag(2))
#' # mu ~ N(0,2.5^2 * I_2)
#' prior_cov = 2.5^2 * diag(2)
#' # Draw from base prior
#' draws0$prior = 
#'   mvtnorm::rmvnorm(ndraws, numeric(2), prior_cov)
#' # Draw from posterior (under base prior)
#' post_cov = 
#'   chol2inv(chol( N * diag(2) + chol2inv(chol(prior_cov)) ))
#' post_mean = 
#'   post_cov %*% diag(2) %*% (N * colMeans(y))
#' draws0$posterior = 
#'   mvtnorm::rmvnorm(ndraws, post_mean, post_cov)
#' 
#' # Create projection matrix as a function of phi
#' P_phi = function(phi) Proj(c(1,phi))
#' # Create prior on phi
#' prior_phi = list(d = function(x) dgamma(x,shape = 50, rate = 50 - 1),
#'                  q = function(p) qgamma(p, shape = 50, rate = 50 - 1),
#'                  r = function(n) rgamma(n, shape = 50, rate = 50 - 1))
#' 
#' if(FALSE){
#'   # Optionally can do this in parallel
#'   library(parallel)
#'   cl = makeCluster(10)
#'   tilted_post_off_subspace = 
#'     SUBSET_IS_gibbs(draws0 = draws0,
#'                     P_phi = P_phi,
#'                     prior_phi = prior_phi,
#'                     n_draws = 1e4,
#'                     cl = cl)
#'   stopCluster(cl)
#' }else{
#'   tilted_post_off_subspace = 
#'     SUBSET_IS_gibbs(draws0 = draws0,
#'                     P_phi = P_phi,
#'                     prior_phi = prior_phi,
#'                     n_draws = 1e4)
#' }
#' tilted_post_off_subspace_summary = 
#'   summary(tilted_post_off_subspace)
#' sum((colMeans(draws0$posterior) - theta_true)^2) # MSE for posterior mean under base prior
#' sum((tilted_post_off_subspace_summary$Estimate[1:2] - theta_true)^2) # MSE for posterior mean under subset prior
#' 
#' 
#' # Now the truth lies near the subspace
#' set.seed(2023)
#' theta_true = c(1.2,1.25)
#' N = 100
#' ndraws = 1e4
#' draws0 = list()
#' # y ~ N(mu, I_2)
#' y = mvtnorm::rmvnorm(N,theta_true,diag(2))
#' # mu ~ N(0,2.5^2 * I_2)
#' prior_cov = 2.5^2 * diag(2)
#' # Draw from base prior
#' draws0$prior = 
#'   mvtnorm::rmvnorm(ndraws, numeric(2), prior_cov)
#' # Draw from posterior (under base prior)
#' post_cov = 
#'   chol2inv(chol( N * diag(2) + chol2inv(chol(prior_cov)) ))
#' post_mean = 
#'   post_cov %*% diag(2) %*% (N * colMeans(y))
#' draws0$posterior = 
#'   mvtnorm::rmvnorm(ndraws, post_mean, post_cov)
#' 
#' tilted_post_near_subspace = 
#'   SUBSET_IS_gibbs(draws0 = draws0,
#'                   P_phi = P_phi,
#'                   prior_phi = prior_phi,
#'                   n_draws = 1e4)
#' tilted_post_near_subspace_summary = 
#'   summary(tilted_post_near_subspace)
#' sum((colMeans(draws0$posterior) - theta_true)^2) # MSE for posterior mean under base prior
#' sum((tilted_post_near_subspace_summary$Estimate[1:2] - theta_true)^2) # MSE for posterior mean under subset prior
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
#' @exportClass subset_is_gibbs


SUBSET_IS_gibbs = function(draws0,
                           P_phi = c("power","geometric")[1],
                           prior_phi,
                           n_draws = nrow(draws0$posterior),
                           phi_sequence = 10,
                           nu,
                           nu_max,
                           min_ESS = 1/10 * nrow(draws0$posterior),
                           verbose = TRUE,
                           cl){
  
  p = ncol(draws0$prior)
  acc_rate = numeric(1)
  
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
  
  # Get I - P(phi)
  if(missing(cl)){
    I_m_P = 
      lapply(1:length(phi_sequence),
             function(i){
               diag(p) - P(phi_sequence[i])
             })
  }else{
    clusterEvalQ(cl,{
        Proj = function(x) 
        {
          proj_inner = function(y) {
            if ("matrix" %in% class(y)) {
              return(tcrossprod(y %*% qr.solve(crossprod(y)), y))
            }
            else {
              return(tcrossprod(y, y)/drop(crossprod(y)))
            }
          }
          if (isTRUE(class(x) == "list")) {
            return(Matrix::as.matrix(Matrix::bdiag(lapply(x, proj_inner))))
          }
          else {
            return(proj_inner(x))
          }
        }
      })
    clusterExport(cl,c("p","P","phi_sequence"),envir=environment())
    I_m_P = 
      parLapply(cl,
                1:length(phi_sequence),
                function(i){
                  diag(p) - P(phi_sequence[i])
                })
  }
  
  if(verbose) cat("\nGetting normalizing constants and weights for IS proposal distributions\n")
  
  
  # Compute w_k(1)
  if(missing(cl)){
    w_1 = 
      lapply(1:length(phi_sequence),
             function(phi){
               sapply(1:nrow(draws0$posterior),
                      function(i){
                        exp(-0.5 * tcrossprod(draws0$posterior[i,],draws0$posterior[i,] %*% I_m_P[[phi]]))
                      })
             })
    
    wtilde_1 = 
      lapply(1:length(phi_sequence),
             function(phi){
               sapply(1:nrow(draws0$prior),
                      function(i){
                        exp(-0.5 * tcrossprod(draws0$prior[i,],draws0$prior[i,] %*% I_m_P[[phi]]))
                      })
             })
  }else{
    clusterExport(cl,c("draws0","I_m_P"),envir = environment())
    w_1 = 
      lapply(1:length(phi_sequence),
             function(phi){
               parSapply(cl,
                         1:ndraws,
                         function(i){
                           exp(-0.5 * tcrossprod(draws0$posterior[i,],draws0$posterior[i,] %*% I_m_P[[phi]]))
                         })
               })
    wtilde_1 = 
      lapply(1:length(phi_sequence),
             function(phi){
               parSapply(cl,
                         1:ndraws,
                         function(i){
                           exp(-0.5 * tcrossprod(draws0$prior[i,],draws0$prior[i,] %*% I_m_P[[phi]]))
                         })
             })
  }
  
  # Create function that creates w_k(\nu)
  if(missing(cl)){
    get_w_k = function(nu,phi_index){
      w_1[[phi_index]]^nu
    }
    get_wtilde_k = function(nu,phi_index){
      wtilde_1[[phi_index]]^nu
    }
  }else{
    clusterExport(cl,c("w_1","wtilde_1"),envir = environment())
    get_w_k = function(nu,phi_index){
      parSapply(cl,w_1[[phi_index]],function(x)x^nu)
    }
    get_wtilde_k = function(nu,phi_index){
      parSapply(cl,wtilde_1[[phi_index]],function(x)x^nu)
    }
  }
  
  # Get nu (if missing) based on mode of prior
  if(missing(nu)){
    ## Find maximum \nu allowed by ESS
    find_max_nu = function(x,phi_index = which.max(phi_pmf)){
      w_k = get_w_k(x,phi_index)
      w_k_sum = sum(w_k)
      w_k = w_k / ifelse(w_k_sum == 0,1,w_k_sum)
      
      ESS = 1 / sum(w_k^2)
      
      (ESS - min_ESS)^2
    }
    if(missing(nu_max)){
      nu_lower_bound = 0
      nu_max = nu_upper_bound = 5
      safety = 0
      while( ( abs(nu_max - nu_upper_bound) / nu_upper_bound < 1e-3) & 
             (safety < 25) ){
        nu_upper_bound = 2 * nu_upper_bound
        nu_max = 
          optimize(find_max_nu,
                   interval = c(nu_lower_bound,
                                nu_upper_bound))$min
        nu_lower_bound = nu_upper_bound
        safety = safety + 1
      }
    }
    
    ## Get optimal nu according to bayes factor
    best_nu = function(nu){
      # Give log( BF(0,nu) )
      log(mean(get_wtilde_k(nu, phi_index = which.max(phi_sequence)))) - 
        log(mean(get_w_k(nu, phi_index = which.max(phi_sequence))))
    }
    nu = optimize(best_nu,
                  interval = c(0,nu_max))$min
  }
  # Now we have nu (if not already provided), and can proceed to algo 2
  
  
  ################
  
  # Get exact (actually, this is MC) values of Z_{\nu,\phi}
  if(!missing(cl)) clusterExport(cl,"nu",envir=environment())
  Z_phi = 
    sapply(1:length(phi_sequence), 
           function(phi_index){mean(get_wtilde_k(nu,phi_index))})
  
  ################
  
  # Get final IS weights
  w_k = 
    lapply(1:length(phi_sequence),get_w_k,nu = nu)
  w_k_norm = 
    lapply(w_k,function(x) x / sum(x))
  
  # Get ESS
  ESS =
    lapply(1:length(phi_sequence),function(phi_index)1 / sum(w_k_norm[[phi_index]]^2))
  
  # Now ready for 2-block Gibbs sampler
  
  
  ################
  
  if(verbose) cat("\nPerforming Gibbs Sampling\n")
  
  # Create array to store samples
  if(is.null(colnames(draws0$posterior))){
    samples = 
      matrix(0.0,
             n_draws,p + 1,
             dimnames = list(NULL,
                             c(paste("theta",1:p,sep="_"),"phi")))
  }else{
    samples = 
      matrix(0.0,
             n_draws,p + 1,
             dimnames = list(NULL,
                             c(colnames(draws0$posterior),"phi")))
  }
  
  # Get initial draw
  ## ...of phi
  samples[1,p + 1] = sample.int(length(phi_sequence),1,prob = phi_pmf)
  # ...of theta
  new_draw_index = 
    sample.int(nrow(draws0$posterior),1,
               prob = w_k_norm[[samples[1,p + 1]]])
  samples[1,1:p] = 
    draws0$posterior[new_draw_index,]
  
  # Set up necessary objects for gibbs sampler
  P_old = P(phi_sequence[ samples[1,p+1] ])
  Z_old = Z_phi[ samples[1,p+1] ]
  
  # Now perform sampling algorithm
  if(verbose) pb = txtProgressBar(0,n_draws,style=3)
  for(it in 2:n_draws){
    
    # Draw phi from prior, adjust with MH
    phi_proposal = sample.int(length(phi_sequence),1,prob = phi_pmf)
    P_new = P(phi_sequence[phi_proposal])
    Z_new = Z_phi[phi_proposal]
    acc_prob = 
      exp(-0.5 * nu * 
            drop(crossprod(samples[it - 1,1:p],(P_old - P_new) %*% samples[it - 1,1:p])) ) * 
      Z_old / Z_new # No prior terms since this cancels in the acceptance prob., i.e., g = pi
    
    if(runif(1) < acc_prob){
      samples[it,p+1] = phi_proposal
      P_old = P_new
      Z_old = Z_new
      acc_rate = acc_rate + 1 / (n_draws - 1)
    }else{
      samples[it,p + 1] = samples[it - 1,p + 1]
    }
    
    # Draw theta via IS
    new_draw_index = 
      sample.int(nrow(draws0$posterior),1,
                 prob = w_k[[ samples[it,p + 1] ]])
    samples[it,1:p] = 
      draws0$posterior[new_draw_index,]
    
    if(verbose) setTxtProgressBar(pb,it)
  }
  
  tab = 
    prop.table(
      table(c(1:length(phi_sequence),samples[,p + 1])) - 1
    )
  
  samples[,p + 1] = 
    phi_sequence[samples[,p + 1]]
  
  
  ################
  
  # Create output object
  
  ret = list(samples = samples,
             acc_rate = acc_rate,
             draws0 = draws0,
             P_phi = P,
             phi_pmf = cbind(phi = phi_sequence, 
                             prior_probability = phi_pmf,
                             posterior_probability = tab),
             prior_phi = prior_phi,
             nu = nu,
             IS_ESS = ESS) 
  class(ret) = c("subset_gibbs","subset_is_gibbs")
  return(ret)
}





