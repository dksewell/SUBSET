#' SUBspace Shrinkage via Exponential Tilting for fixed 
#' projection matrix
#' 
#' SUBSET_IS_fixed() will perform importance sampling 
#' using an exponential tilted prior that shrinks the 
#' parameters towards a linear subspace.
#' 
#' SUBSET takes a base prior \eqn{\pi_0} and tilts it towards
#' a linear subspace via the multiplicative term 
#' \deqn{e^{-\frac{v}{2}\theta'(I-P)\theta}.}
#' The user must perform the initial Bayesian inference for 
#' this base prior, and input posterior draws. 
#' Then, for a fixed projection matrix \eqn{P}, SUBSET_IS_fixed() will 
#' perform the importance sampling algorithm corresponding to 
#' the tilted prior.  \code{v_sequence} is a user provided positive 
#' number or vector of positive numbers that control the shrinkage: 
#' higher values imply more shrinkage of the posterior while \eqn{v=0} 
#' corresponds to the base prior with no shrinkage. 
#' *Note that that the parameter space must be convex.* (It usually is.)
#' 
#' 
#' @param draws0 named list with names 'posterior' and 'prior'. 
#' Each ought to be a matrix of posterior or prior draws (each row is a draw),
#' all of which should correspond to the base prior.
#' @param Proj_mat matrix.  Projection matrix for the 
#' desired linear subspace you wish to shrink towards.
#' @param CI_level numeric in (0,1). Level of the 
#' posterior credible interval.
#' @param nu shrinkage parameter.  If missing, the optimal nu via Bayes factors will be used
#' @param nu_max maximum value of nu to be considered.
#' @param min_ESS integer.  Minimum number of effective sample size in considering values of nu. 
#' @param verbose logical. Should any output be provided?
#' @param cl optional object of class c("SOCKcluster", "cluster") generated from parallel::makeCluster.
#' 
#' @return Object of class "subset_fixed", "subset_IS_fixed", with the following named elements:
#' \itemize{
#' \item summary, a data frame with the following structure:
#'    \itemize{
#'      \item variable - Variable names extracted from colnames(Sigma0)
#'      \item mean - Posterior mean
#'      \item [(1-CI_level)/2]$% - lower CI bound
#'      \item [1 - (1-CI_level)/2]$% - upper CI bound
#'    }
#' \item nu the SUBSET shrinkage parameter either fixed by the user or selected by Bayes factor
#' \item ESS the effective sample size
#' \item BF_favoring_nu The Bayes factor in favor of shrinkage vs. no shrinkage
#' \item proposal_draws Draws input by user (posterior draws under the base prior).
#' \item is_weights Importance sampling weights corresponding to the proposal draws under the 
#' SUBSET prior
#' \item base_prior_draws Draws input by the user (prior draws under the base prior).
#' }
#' 
#' @examples 
#' # Simulate where truth lies far from subspace
#' set.seed(2023)
#' theta_true = c(1,1.5)
#' N = 100
#' ndraws = 1e5
#' draws0 = list()
#' # y ~ N(mu, I_2)
#' y = mvtnorm::rmvnorm(N,theta_true,diag(2))
#' # mu ~ N(0,2.5^2 * I_2)
#' prior_cov = 2.5^2 * diag(2)
#' # Draw from base prior
#' draws0$prior = 
#'   mvtnorm::rmvnorm(ndraws, numeric(2), 2.5^2 * diag(2))
#' # Draw from posterior (under base prior)
#' post_cov = 
#'   chol2inv(chol( N * diag(2) + chol2inv(chol(prior_cov)) ))
#' post_mean = 
#'   post_cov %*% diag(2) %*% (N * colMeans(y))
#' draws0$posterior = 
#'   mvtnorm::rmvnorm(ndraws, post_mean, post_cov)
#' 
#' tilted_post_off_subspace = 
#'   SUBSET_IS_fixed(draws0, 
#'                   Proj_mat = Proj(c(1,1)))
#' tilted_post_off_subspace$nu # The shrinkage parameter selected by Bayes factor
#' tilted_post_off_subspace$BF_favoring_nu # The Bayes factor favoring nu = tilted_post$nu vs. nu = 0
#' tilted_post_off_subspace$ESS # Effective sample size
#' tilted_post_off_subspace$summary
#' sum((colMeans(draws0$posterior) - theta_true)^2) # MSE for posterior mean under base prior
#' sum((tilted_post_off_subspace$summary$posterior_mean - theta_true)^2) # MSE for posterior mean under subset prior
#' 
#' # Now the truth lies near on the subspace
#' set.seed(2023)
#' theta_true = c(1.2,1.25)
#' N = 100
#' ndraws = 1e5
#' draws0 = list()
#' # y ~ N(mu, I_2)
#' y = mvtnorm::rmvnorm(N,theta_true,diag(2))
#' # mu ~ N(0,2.5^2 * I_2)
#' prior_cov = 2.5^2 * diag(2)
#' # Draw from base prior
#' draws0$prior = 
#'   mvtnorm::rmvnorm(ndraws, numeric(2), 2.5^2 * diag(2))
#' # Draw from posterior (under base prior)
#' post_cov = 
#'   chol2inv(chol( N * diag(2) + chol2inv(chol(prior_cov)) ))
#' post_mean = 
#'   post_cov %*% diag(2) %*% (N * colMeans(y))
#' draws0$posterior = 
#'   mvtnorm::rmvnorm(ndraws, post_mean, post_cov)
#' 
#' tilted_post_near_subspace = 
#'   SUBSET_IS_fixed(draws0, 
#'                   Proj_mat = Proj(c(1,1)))
#' tilted_post_near_subspace$nu # The shrinkage parameter selected by Bayes factor
#' tilted_post_near_subspace$BF_favoring_nu # The Bayes factor favoring nu = tilted_post$nu vs. nu = 0
#' tilted_post_near_subspace$ESS # Effective sample size
#' tilted_post_near_subspace$summary
#' sum((colMeans(draws0$posterior) - theta_true)^2) # MSE for posterior mean under base prior
#' sum((tilted_post_near_subspace$summary$posterior_mean - theta_true)^2) # MSE for posterior mean under subset prior
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
#' @exportClass subset_fixed
#' @exportClass subset_IS_fixed



SUBSET_IS_fixed = function(draws0,
                           Proj_mat = Proj(matrix(1,ncol(draws0),1)),
                           CI_level = 0.95,
                           nu,
                           nu_max,
                           min_ESS = 1/10 * nrow(draws0$posterior),
                           verbose = TRUE,
                           cl){
  ndraws = sapply(draws0,nrow)
  CI_level = 1 - CI_level
  p = ncol(draws0$posterior)
  I_m_P = diag(p) - Proj_mat
  if(verbose)if(missing(nu))cat("\nNo nu provided.  Picking level of shrinkage via Bayes Factor.\n")
  
  if(!missing(cl)) if(is.numeric(cl)) cl = makeCluster(min(cl,detectCores() - 1))
  
  
  # Find optimal nu
  ## Compute w_k(1)
  if(missing(cl)){
    w_1 = 
      sapply(1:ndraws["posterior"],
             function(i){
               exp(-0.5 * tcrossprod(draws0$posterior[i,],draws0$posterior[i,] %*% I_m_P))
             })
    wtilde_1 = 
      sapply(1:ndraws["prior"],
             function(i){
               exp(-0.5 * tcrossprod(draws0$prior[i,],draws0$prior[i,] %*% I_m_P))
             })
  }else{
    clusterExport(cl,c("draws0","I_m_P"),envir = environment())
    w_1 = 
      parSapply(cl,
                1:ndraws,
                function(i){
                  exp(-0.5 * tcrossprod(draws0$posterior[i,],draws0$posterior[i,] %*% I_m_P))
                })
    wtilde_1 = 
      parSapply(cl,
                1:ndraws,
                function(i){
                  exp(-0.5 * tcrossprod(draws0$prior[i,],draws0$prior[i,] %*% I_m_P))
                })
  }
  
  ## Create function that creates w_k(\nu)
  if(missing(cl)){
    get_w_k = function(nu){
      w_1^nu
    }
    get_wtilde_k = function(nu){
      wtilde_1^nu
    }
  }else{
    clusterExport(cl,c("w_1","wtilde_1"),envir = environment())
    get_w_k = function(nu){
      parSapply(cl,w_1,function(x)x^nu)
    }
    get_wtilde_k = function(nu){
      parSapply(cl,wtilde_1,function(x)x^nu)
    }
  }
  
  if(missing(nu)){
    ## Find maximum \nu allowed by ESS
    find_max_nu = function(x){
      w_k = get_w_k(x)
      w_k_sum = sum(w_k)
      w_k = w_k / ifelse(w_k_sum == 0,1,w_k_sum)
      
      ESS = 1 / sum(w_k^2)
      
      (ESS - min_ESS)^2
    }
    if(missing(nu_max)){
      nu_max = nu_upper_bound = 1e3
      safety = 0
      if(verbose) cat("\n---Finding maximum nu to satisfy ESS constraints\n")
      while( ( abs(nu_max - nu_upper_bound) / nu_upper_bound < 1e-3) & safety < 100){
        nu_upper_bound = 2 * nu_upper_bound
        nu_max = 
          optimize(find_max_nu,
                   interval = c(0,nu_upper_bound))$min
        safety = safety + 1
      }
    }
    
    ## Get optimal nu according to bayes factor
    best_nu = function(nu){
      # Give log( BF(0,nu) )
      log(mean(get_wtilde_k(nu))) - 
        log(mean(get_w_k(nu)))
    }
    nu = optimize(best_nu,
                  interval = c(0,nu_max))$min
  }
  
  if(!missing(cl)) clusterExport(cl,"nu",envir=environment())
  
  # Get final IS weights
  w_k = get_w_k(nu)
  BF_nu_0 = 
    exp(
      log(mean(get_w_k(nu))) - 
        log(mean(get_wtilde_k(nu)))
    )
  w_k = w_k / sum(w_k)
  
  # Get ESS
  ESS = 1 / sum(w_k^2)
  
  # Get point estimates
  theta_hat = 
    apply(draws0$posterior,2,weighted.mean,w = w_k)
  
  # Get CI
  get_CI_from_IS = function(x,w){
    ord = order(x)
    w = w[ord]
    x = x[ord]
    w_cum = cumsum(w)
    
    c(lower = x[min(which(w_cum >= CI_level/2))],
      upper = x[min(which(w_cum >= 1 - CI_level/2))]
    )
  }
  
  # Put results together
  if(is.null(colnames(draws0$posterior))){
    if(!is.null(colnames(draws0$prior))){
      colnames(draws0$posterior) = 
        colnames(draws0$prior)
    }else{
      colnames(draws0$posterior) = 
        colnames(draws0$prior) = 
        paste("theta",1:p,sep="_")
    }
  }
  results = list()
  results$summary = 
    data.frame(variable = colnames(draws0$posterior),
               posterior_mean = theta_hat,
               lower = numeric(p),
               upper = numeric(p))
  for(j in 1:p){
    temp = 
      get_CI_from_IS(draws0$posterior[,j],w_k)
    results$summary$lower[j] = temp["lower"]
    results$summary$upper[j] = temp["upper"]
  }
  colnames(results$summary)[3:4] = 
    c(paste0(CI_level/2 * 100,"%"),
      paste0((1-CI_level/2) * 100,"%"))
  
  results$nu = nu
  results$ESS = ESS
  results$BF_favoring_nu = BF_nu_0
  results$proposal_draws = draws0$posterior
  results$is_weights = w_k
  results$base_prior_draws = draws0$prior
  
  class(results) = c("subset_fixed","subset_IS_fixed")
  
  rownames(results$summary) = NULL
  return(results)
}
