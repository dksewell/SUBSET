#' SUBspace Shrinkage via Exponential Tilting for estimated 
#' projection matrix
#' 
#' @param object subset_gibbs object
#' @param level the level of credible intervals to be provided
#' @param v_index Optional Integer giving the index of 
#' the posterior samples corresponding to the desired 
#' value of v, i.e., the exponential tilting parameter. 


summary.subset_gibbs = function(object,level = 0.95,v_index,...){
  alpha = 1 - level
  p = dim(object$samples)[2]
  # summ = array(0.0, c(dim(object$samples)[2],4,dim(object$samples)[3]))
  summ = data.frame(Variable = rep(colnames(object$samples),dim(object$samples)[3]),
                    Estimate = 0.0,
                    Lower = 0.0,
                    Upper = 0.0,
                    v = rep(object$v_sequence,each = p))
  
  for(v in 1:dim(object$samples)[3]){
    summ$Estimate[p * (v - 1) + 1:p] = 
      colMeans(object$samples[,,v])
    summ[p * (v - 1) + 1:p,c("Lower","Upper")] = 
      t(apply(object$samples[,,v],2,quantile, probs = c(alpha / 2, 1 - alpha/2)))
  }
  
  if(missing(v_index)){
    return(summ)
  }else{
    return(summ[which(summ$v == object$v_sequence[v_index]),])
  }
}
