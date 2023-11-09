#' SUBspace Shrinkage via Exponential Tilting for estimated 
#' projection matrix
#' 
#' @param object subset_gibbs object
#' @param level the level of credible intervals to be provided
#' @param v Optional. The desired value of v, i.e., the exponential tilting parameter. 
#' @export


summary.subset_gibbs = function(object,level = 0.95,v,...){
  
  alpha = 1 - level
  p = dim(object$samples)[2]
  
  if("subset_asymp_gibbs" %in% class(object)){
    summ = data.frame(Variable = rep(colnames(object$samples),dim(object$samples)[3]),
                      Estimate = 0.0,
                      Lower = 0.0,
                      Upper = 0.0,
                      v = rep(object$v_sequence,each = p))
    
    for(vv in 1:dim(object$samples)[3]){
      summ$Estimate[p * (vv - 1) + 1:p] = 
        colMeans(object$samples[,,vv])
      summ[p * (vv - 1) + 1:p,c("Lower","Upper")] = 
        t(apply(object$samples[,,vv],2,quantile, probs = c(alpha / 2, 1 - alpha/2)))
    }
    
    
    if(missing(v)){
      return(summ)
    }else{
      return(summ[which(summ$v == v),])
    }
  }
  
  
  
  if("subset_is_gibbs" %in% class(object)){
    summ = data.frame(Variable = colnames(object$samples),
                      Estimate = colMeans(object$samples),
                      Lower = rep(0.0,p),
                      Upper = rep(0.0,p),
                      v = rep(object$nu,p))
    summ[,c("Lower","Upper")] = 
      temp = 
      t(apply(object$samples,2,quantile, probs = c(alpha / 2, 1 - alpha/2)))
    colnames(summ)[3:4] = colnames(temp)
    rownames(summ) = NULL
    return(summ)
  }
}
