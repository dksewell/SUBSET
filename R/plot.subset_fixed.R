#' SUBspace Shrinkage via Exponential Tilting for estimated 
#' projection matrix
#' 
#' plot.subset_fixed will plot the transition of each 
#' parameter as v, the shrinkage hyperparameter, 
#' varies.
#' 
#' @param x Object of class "subset_gibbs", returned by 
#' the function SUBSET_gibbs().
#' @param include_CI logical.
#' @param ... Additional graphical parameters.
#' @export


plot.subset_fixed = function(x,
                             include_CI = TRUE,
                             ...){
  if("subset_IS_fixed" %in% class(x)){
    x = x$summary
    lb = unlist(x[,3])
    ub = unlist(x[,4])
  }else{
    lb = unlist(x[,4])
    ub = unlist(x[,5])
  }
  
  if(include_CI){
    yl = range(c(lb,ub))
  }else{
    yl = range(unlist(x$mean))
  }
  
  
  vars = unique(x$variable)
  p = length(vars)
  x = x[order(x$v),]
  
  plot(0,
       type = 'n',
       bty = 'l',
       ylab = '',
       xlab = 'v (SUBSET shrinkage strength)',
       xlim = range(x$v)*c(0.95,1.05),
       ylim = yl,
       xaxt = 'n')
  axis(1, at = pretty(sort(unique(x$v))))
  if(include_CI){
    for(j in 1:p){
      subset_index = which(x$variable == vars[j])
      polygon(x = c(x$v[subset_index],
                    rev(x$v[subset_index])),
              y = c(lb[subset_index],rev(ub[subset_index])),
              border = FALSE,
              col = adjustcolor(MetBrewer::met.brewer("Austria",p)[j],0.3))
    }
  }
  for(j in 1:p){
    subset_index = which(x$variable == vars[j])
    lines(x$mean[subset_index] ~ x$v[subset_index],
          lwd = 2,
          col = MetBrewer::met.brewer("Austria",p)[j])
  }
  
}
