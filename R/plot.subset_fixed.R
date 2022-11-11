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


plot.subset_fixed = function(x,
                             include_CI = TRUE,
                             ...){
  if(include_CI){
    yl = range(unlist(x[,c(2,4,5)]))
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
       xlab = '',
       xlim = range(x$v)*c(0.95,1.05),
       ylim = yl,
       xaxt = 'n')
  axis(1, at = pretty(sort(unique(x$v))))
  if(include_CI){
    for(j in 1:p){
      x_subset = x[which(x$variable == vars[j]),]
      polygon(x = c(x_subset$v,rev(x_subset$v)),
              y = c(x_subset$`2.5%`,rev(x_subset$`97.5%`)),
              border = FALSE,
              col = adjustcolor(MetBrewer::met.brewer("Austria",p)[j],0.3))
    }
  }
  for(j in 1:p){
    x_subset = x[which(x$variable == vars[j]),]
    lines(x_subset$mean ~ x_subset$v,
          lwd = 2,
          col = MetBrewer::met.brewer("Austria",p)[j])
  }
  
}