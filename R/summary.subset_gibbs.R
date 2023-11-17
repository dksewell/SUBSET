#' SUBspace Shrinkage via Exponential Tilting for estimated 
#' projection matrix
#' 
#' @param object subset_gibbs object
#' @param level the level of credible intervals to be provided
#' @param v Optional. The desired value of v, i.e., the exponential tilting parameter. 
#' @export summary.subset_gibbs
#' @export


summary.subset_gibbs = function(object,level = 0.95,v,...){
  
  alpha = 1 - level
  p = dim(object$samples)[2]
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
