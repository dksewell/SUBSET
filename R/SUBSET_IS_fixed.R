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
#' @examples 
#' post_mean_from_pi0 = c(1,1.5)
#' post_cov_from_pi0 = diag(c(1,3))
#' P = Proj(rep(1,2))
#' SUBSET_fixed(Mean0 = post_mean_from_pi0, Sigma = post_cov_from_pi0, Proj_mat = P, v_sequence = c(0.5,100))
#' 


SUBSET_IS_fixed = function(draws0,
                           Proj_mat = Proj(matrix(1,ncol(Sigma0),1)),
                           v_sequence = seq(0.25,5,by = 0.25),
                           CI_level = 0.95,
                           min_ESS = 1/2 * nrow(draws0),
                           n_u_values = 100,
                           verbose = FALSE){
  CI_level = 1 - CI_level
  p = ncol(draws0)
  v_len = length(v_sequence)
  I_m_P = diag(p) - Proj_mat
  
  TT = function(nu){
    sapply(1:ndraws,
           function(i){
             exp(-0.5 * nu * tcrossprod(draws0[i,],draws0[i,] %*% I_m_P))
           })
  }
  
  ESS = function(u,Tnu){
    w = Tnu^(u^2)
    w = w / sum(w)
    1 / sum(w^2)
  }
  
  draws_u = results = list()
  
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
    
    if(!is.null(colnames(draws0))){
      results[[v]] = 
        data.frame(variable = colnames(draws0),
                   mean = numeric(p),
                   lower = numeric(p),
                   upper = numeric(p),
                   v = v_sequence[v])
    }else{
      results[[v]] = 
        data.frame(variable = paste0("v",1:p),
                   mean = numeric(p),
                   lower = numeric(p),
                   upper = numeric(p),
                   v = v_sequence[v])
    }
    
    Tnu = TT(v_sequence[v])
    
    ESS_values = 
      sapply(seq(0,1,l = n_u_values+1)[-1],ESS,Tnu = Tnu)
    u = seq(0,1,l=n_u_values)[max(which(ESS_values > min_ESS))]
    results[[v]]$ESS = ESS_values[max(which(ESS_values > min_ESS))]
    results[[v]]$u = u
    
    draws_u[[v]] = 
      cbind(u * draws0 + (1 - u) * draws0 %*% Proj_mat,
            weight = Tnu / sum(Tnu))
    if(!is.null(colnames(draws0))) colnames(draws_u)[1:p] = colnames(draws0)
    
    results[[v]]$mean = 
      apply(draws_u[[v]][,1:p], 2, weighted.mean, w = draws_u[[v]][,"weight"])
    results[[v]]$mean = 
      apply(draws_u[[v]][,1:p], 2, weighted.mean, w = draws_u[[v]][,"weight"])
    
    for(j in 1:p){
      x = draws_u[[v]][,j]
      ord = order(x)
      w = draws_u[[v]][ord,"weight"]
      x = x[ord]
      rw <- cumsum(w)
      selection <- min(which(rw >= CI_level/2))
      if (rw[selection] == CI_level/2){
        results[[v]]$lower[j] = mean(x[selection:(selection + 1)])
      }else{
        results[[v]]$lower[j] = x[selection]
      }
      selection <- min(which(rw >= 1 - CI_level/2))
      if (rw[selection] == 1 - CI_level/2){
        results[[v]]$upper[j] = mean(x[selection:(selection + 1)])
      }else{
        results[[v]]$upper[j] = x[selection]
      }
    }
    
  }
  
  results = do.call(rbind,results)
  
  if(!is.null(colnames(draws0))){
    results = 
      rbind(data.frame(variable = colnames(draws0),
                       mean = colMeans(draws0),
                       lower = apply(draws0,2,quantile,prob = CI_level/2),
                       upper = apply(draws0,2,quantile,prob = 1 - CI_level/2),
                       v = 0.0,
                       ESS = nrow(draws0),
                       u = NA),
            results)
  }else{
    results = 
      rbind(data.frame(variable = paste0("v",1:ncol(draws0)),
                       mean = colMeans(draws0),
                       lower = apply(draws0,2,quantile,prob = CI_level/2),
                       upper = apply(draws0,2,quantile,prob = 1 - CI_level/2),
                       v = 0.0,
                       ESS = nrow(draws0),
                       u = NA),
            results)
  }
  colnames(results)[3:4] = 
    c(paste0(CI_level/2 * 100,"%"),
      paste0((1-CI_level/2) * 100,"%"))
  results = list(summary = results,
                 is_draws = draws_u)
  class(results) = c("subset_SI_fixed")
  rownames(results$summary) = NULL
  return(results)
}