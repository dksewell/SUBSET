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
#' @param Proj_mat matrix.  Projection matrix for the 
#' desired linear subspace you wish to shrink towards.
#' @param v_sequence vector of positive values. See details.
#' @param CI_level numeric between (0,1). Level of the 
#' posterior credible interval.
#' @param verbose logical. Should any output be provided?  Note that
#' this is typically very fast, so usually not needed.
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
#' post_mean_from_pi0 = c(1,1.5)
#' post_cov_from_pi0 = diag(c(1,3))
#' P = Proj(rep(1,2))
#' SUBSET_fixed(Mean0 = post_mean_from_pi0, Sigma = post_cov_from_pi0, Proj_mat = P, v_sequence = c(0.5,100))
#' 
#' @export


SUBSET_fixed = function(Mean0,
                        Sigma0,
                        Proj_mat = Proj(matrix(1,ncol(Sigma0),1)),
                        v_sequence = seq(0.25,5,by = 0.25),
                        CI_level = 0.95,
                        verbose = FALSE){
  CI_level = 1 - CI_level
  p = ncol(Sigma0)
  v_len = length(v_sequence)
  Mean0 = drop(Mean0)
  
  Sigma0_inv = chol2inv(chol(Sigma0))
  
  results = list()
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
    
    if(!is.null(colnames(Sigma0))){
      results[[v]] = 
        data.frame(variable = colnames(Sigma0),
                   mean = numeric(p),
                   sd = numeric(p),
                   lower = numeric(p),
                   upper = numeric(p),
                   v = numeric(p))
    }else{
      results[[v]] = 
        data.frame(variable = paste0("v",1:ncol(Sigma0)),
                   mean = numeric(p),
                   sd = numeric(p),
                   lower = numeric(p),
                   upper = numeric(p),
                   v = numeric(p))
    }
    
    
    Sigma_inv = Sigma0_inv + v_sequence[v] * (diag(p) - Proj_mat)
    Sigma =  chol2inv(chol(Sigma_inv))
    Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
    
    results[[v]][,2] = Mean
    results[[v]][,3] = sqrt(diag(Sigma))
    results[[v]][,4] = qnorm(CI_level/2, Mean, sqrt(diag(Sigma)))
    results[[v]][,5] = qnorm(1 - CI_level/2, Mean, sqrt(diag(Sigma)))
    results[[v]][,6] = v_sequence[v]
  }
  
  results = do.call(rbind,results)
  if(!is.null(colnames(Sigma0))){
    results = 
      rbind(data.frame(variable = colnames(Sigma0),
                       mean = Mean0,
                       sd = sqrt(diag(Sigma0)),
                       lower = qnorm(CI_level/2,Mean0,sqrt(diag(Sigma0))),
                       upper = qnorm(1 - CI_level/2,Mean0,sqrt(diag(Sigma0))),
                       v = 0.0),
            results)
  }else{
    results = 
      rbind(data.frame(variable = paste0("v",1:ncol(Sigma0)),
                       mean = Mean0,
                       sd = sqrt(diag(Sigma0)),
                       lower = qnorm(CI_level/2,Mean0,sqrt(diag(Sigma0))),
                       upper = qnorm(1 - CI_level/2,Mean0,sqrt(diag(Sigma0))),
                       v = 0.0),
            results)
  }
  colnames(results)[4:5] = 
    c(paste0(CI_level/2 * 100,"%"),
      paste0((1-CI_level/2) * 100,"%"))
  class(results) = c("subset_fixed","subset_asymp_fixed","data.frame")
  rownames(results) = NULL
  return(results)
}

