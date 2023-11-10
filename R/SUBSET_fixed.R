#' SUBspace Shrinkage via Exponential Tilting for fixed 
#' projection matrix
#' 
#' SUBSET_fixed() will compute the asymptotic normal 
#' approximation to the posterior using an exponential 
#' tilted prior that shrinks the  parameters towards a 
#' linear subspace.
#' 
#' SUBSET takes a base prior \eqn{\pi_0} and tilts it towards
#' a linear subspace via the multiplicative term 
#' \deqn{e^{-\frac{v}{2}\theta'(I-P)\theta}.}
#' The user must perform the initial Bayesian inference for 
#' this base prior, and compute the posterior mean and covariance. 
#' *It is imperative that the necessary conditions for an asymptotic 
#' normal approximation of the posterior hold for the base prior.* 
#' Then, for a fixed projection matrix \eqn{P}, SUBSET_fixed() will 
#' compute the asymptotic normal distribution corresponding to 
#' the tilted prior.  \code{v_sequence} is a user provided positive 
#' number or vector of positive numbers that control the shrinkage: 
#' higher values imply more shrinkage of the posterior while \eqn{v=0} 
#' corresponds to the base prior with no shrinkage.
#' 
#' 
#' @param Mean0 vector. The posterior mean corresponding
#' to the base prior (i.e., untilted prior).
#' @param Sigma0 matrix. The posterior covariance 
#' matrix corresponding to the base prior. NOTE: The
#' variable names are extracted from Sigma0, so use 
#' colnames(Sigma0) = ...
#' @param prior_draws matrix of prior draws (each row is a draw), 
#' corresponding to the base prior.
#' @param Proj_mat matrix.  Projection matrix for the 
#' desired linear subspace you wish to shrink towards.
#' @param CI_level numeric between (0,1). Level of the 
#' posterior credible interval.
#' @param nu shrinkage parameter.  If missing, the optimal nu via Bayes factors will be used
#' @param nu_max maximum value of nu to be considered.
#' @param cl optional object of class c("SOCKcluster", "cluster") generated from parallel::makeCluster.
#' 
#' @return Object of class "subset_fixed", "subset_asymp_fixed", "data.frame" with the 
#' following structure:
#' \itemize{
#'    \item variable - Variable names extracted from colnames(Sigma0)
#'    \item mean - Posterior mean
#'    \item sd - Posterior sd
#'    \item [(1-CI_level)/2]$% - lower CI bound
#'    \item [1 - (1-CI_level)/2]$% - upper CI bound
#'    \item v - SUBSET hyperparameter (shrinkage strength)
#' }
#' 
#' @examples 
#' # Simulate where truth lies far from subspace
#' set.seed(2023)
#' theta_true = c(1,1.5)
#' N = 100
#' ndraws = 1e5
#' # y ~ N(mu, I_2)
#' y = mvtnorm::rmvnorm(N,theta_true,diag(2))
#' # mu ~ N(0,2.5^2 * I_2)
#' prior_cov = 2.5^2 * diag(2)
#' # Draw from base prior
#' prior_draws = 
#'   mvtnorm::rmvnorm(ndraws, numeric(2), prior_cov)
#' # Get posterior parameters (under base prior)
#' Sigma0 = 
#'   chol2inv(chol( N * diag(2) + chol2inv(chol(prior_cov)) ))
#' Mean0 = 
#'   Sigma0 %*% diag(2) %*% (N * colMeans(y))
#' 
#' tilted_post_off_subspace = 
#'   SUBSET_fixed(Mean0 = Mean0,
#'                Sigma0 = Sigma0,
#'                prior_draws = prior_draws, 
#'                Proj_mat = Proj(c(1,1)))
#' tilted_post_off_subspace$nu # The shrinkage parameter selected by Bayes factor
#' tilted_post_off_subspace$BF_favoring_nu # The Bayes factor favoring nu = tilted_post$nu vs. nu = 0
#' tilted_post_off_subspace$summary
#' sum((Mean0 - theta_true)^2) # MSE for posterior mean under base prior
#' sum((tilted_post_off_subspace$summary$posterior_mean - theta_true)^2) # MSE for posterior mean under subset prior
#' 
#' # Now the truth lies near on the subspace
#' set.seed(2023)
#' theta_true = c(1.2,1.25)
#' N = 100
#' ndraws = 1e5
#' # y ~ N(mu, I_2)
#' y = mvtnorm::rmvnorm(N,theta_true,diag(2))
#' # mu ~ N(0,2.5^2 * I_2)
#' prior_cov = 2.5^2 * diag(2)
#' # Draw from base prior
#' prior_draws = 
#'   mvtnorm::rmvnorm(ndraws, numeric(2), prior_cov)
#' # Get posterior parameters (under base prior)
#' Sigma0 = 
#'   chol2inv(chol( N * diag(2) + chol2inv(chol(prior_cov)) ))
#' Mean0 = 
#'   Sigma0 %*% diag(2) %*% (N * colMeans(y))
#' 
#' tilted_post_near_subspace = 
#'   SUBSET_fixed(Mean0 = Mean0,
#'                Sigma0 = Sigma0,
#'                prior_draws = prior_draws, 
#'                Proj_mat = Proj(c(1,1)))
#' tilted_post_near_subspace$nu # The shrinkage parameter selected by Bayes factor
#' tilted_post_near_subspace$BF_favoring_nu # The Bayes factor favoring nu = tilted_post$nu vs. nu = 0
#' tilted_post_near_subspace$summary
#' sum((Mean0 - theta_true)^2) # MSE for posterior mean under base prior
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


SUBSET_fixed = function(Mean0,
                        Sigma0,
                        prior_draws,
                        Proj_mat = Proj(matrix(1,ncol(Sigma0),1)),
                        CI_level = 0.95,
                        nu,
                        nu_max,
                        cl){
  CI_level = 1 - CI_level
  p = ncol(Sigma0)
  Mean0 = drop(Mean0)
  
  Sigma0_inv = chol2inv(chol(Sigma0))
  
  
  if(!is.null(colnames(Sigma0))){
    results = 
      data.frame(variable = colnames(Sigma0),
                 posterior_mean = numeric(p),
                 sd = numeric(p),
                 lower = numeric(p),
                 upper = numeric(p),
                 v = numeric(p))
  }else{
    results = 
      data.frame(variable = paste0("theta",1:ncol(Sigma0)),
                 posterior_mean = numeric(p),
                 sd = numeric(p),
                 lower = numeric(p),
                 upper = numeric(p))
  }
  
  # Get optimal nu
  Sigma0_logdet = determinant(Sigma0)$modulus
  if(missing(nu)){
    
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
      BF_nu_0 = exp(opt$objective)
      
      Sigma_inv = Sigma0_inv + nu * I_m_P
      Sigma =  chol2inv(chol(Sigma_inv))
      Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
    }
  }else{
    Sigma_inv = Sigma0_inv + nu * I_m_P
    Sigma =  chol2inv(chol(Sigma_inv))
    Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
    
    BF_nu_0 =
      exp(
        -0.5 * Sigma0_logdet -
          0.5 * drop( crossprod(Mean0,Sigma0_inv %*% Mean0) ) -
          0.5 * determinant(Sigma_inv)$modulus +
          0.5 * drop( crossprod(Mean,Sigma_inv %*% Mean) ) - 
          log(mean(sapply(wtilde_1,function(x)x^nu)))
      )
  }
  
  results[,2] = Mean
  results[,3] = sqrt(diag(Sigma))
  results[,4] = qnorm(CI_level/2, Mean, sqrt(diag(Sigma)))
  results[,5] = qnorm(1 - CI_level/2, Mean, sqrt(diag(Sigma)))
  
  colnames(results)[4:5] = 
    c(paste0(CI_level/2 * 100,"%"),
      paste0((1-CI_level/2) * 100,"%"))
  
  results = list(summary = results)
  results$nu = nu
  results$BF_favoring_nu = as.numeric(BF_nu_0)
  
  
  class(results) = c("subset_fixed","subset_asymp_fixed")
  return(results)
}

