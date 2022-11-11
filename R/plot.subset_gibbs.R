#' SUBspace Shrinkage via Exponential Tilting for estimated 
#' projection matrix
#' 
#' plot.subset_gibbs will provide several plotting functions 
#' for an object returned by SUBSET_gibbs().  Obtain 
#' traceplots, cumulative mean plots, histograms of 
#' the model parameters, and histogram of the 
#' projection matrix parameter.
#' 
#' @param x Object of class "subset_gibbs", returned by 
#' the function SUBSET_gibbs().
#' @param type character denoting traceplots ("trace"), 
#' cumulative mean plots ("cummean"), histograms 
#' of the model parameters ("theta"), or a histogram 
#' of the projection matrix parameter ("phi"), or 
#' transition plot of the theta values over v ("transition").
#' @param v_index Integer giving the "slice" index of 
#' the posterior samples corresponding to the desired 
#' value of v, i.e., the exponential tilting parameter. 
#' Only used for type = "theta" or "phi".
#' @param include_prior Superimpose the marginal prior 
#' densities on the histograms if type = "theta" or "phi".
#' @param ... Additional graphical parameters.


plot.subset_gibbs = function(x,
                             type = c("trace","cummean","theta","phi","transition")[1],
                             v_index,
                             include_prior = TRUE,
                             ...){
  if(type == "trace"){
    for(j in 1:ncol(x$samples)){
      matplot(x$samples[,j,],
              type = 'l',
              bty = 'l',
              lty = 1,
              lwd = 2,
              col = gray(0.5,0.1),
              xlab = "Gibbs iteration",
              ylab = colnames(x$samples)[j],
              ...)
      
      if(j < ncol(x$samples)){
        cat("\nHit RETURN for next plot\n")
        readline()
      }
    }
  }
  
  if(type == "cummean"){
    for(j in 1:(ncol(x$samples) - 1)){
      cummean = function(z) cumsum(z) / c(1:NROW(z))
      
      cummeans = 
        apply(x$samples[,j,],2,cummean)
      matplot(cummeans,
              type = 'l',
              bty = 'l',
              lty = 1,
              lwd = 2,
              col = gray(seq(0.9,0.1,l = dim(x$samples)[3])),
              xlab = "Gibbs iteration",
              ylab = colnames(x$samples)[j],
              ...)
      
      if(j < (ncol(x$samples) - 1)){
        cat("\nHit RETURN for next plot\n")
        readline()
      }
    }
  }
  
  
  if(type == "theta"){
    if(missing(v_index)) stop("Must provide which value of v you want via the v_index argument.")
    for(j in 1:ncol(x$samples)){
      
      hist(x$samples[,j,v_index],
           freq = FALSE,
           yaxt = 'n',
           # bty = 'l',
           ylab = '',
           main = '',
           xlab = colnames(x$samples)[j],
           ...)
      if(include_prior){
        if( j < ncol(x$samples) ){
          lines(density(x$pi0_samples[,j],adjust = 2),
                lwd = 2,
                col = gray(0.3),
                ...)
        }else{
          lines(density(x$prior_phi$r(5e3),adjust = 2, from = 0),
                lwd = 2,
                col = gray(0.3),
                ...)
        }
      }
      
      if(j < ncol(x$samples)){
        cat("\nHit RETURN for next plot\n")
        readline()
      }
    }
  }
  
  if(type == "phi"){
    
    if(missing(v_index)){
      densities = 
        lapply(1:dim(x$samples)[3],
               function(v){
                 density(x$samples[,dim(x$samples)[2],v],
                         adjust = 2,
                         from = 0)
                 })
      
      xl = range(sapply(densities, function(z) range(z$x)))
      yl = range(sapply(densities, function(z) range(z$y)))
      
      plot(0,
           bty = 'l',
           type = 'n',
           xlab = expression(phi),
           ylab = '',
           yaxt = 'n',
           ylim = yl,
           xlim = xl,
           ...)
      for(v in 1:dim(x$samples)[3]){
        lines(densities[[v]],
              lwd = 2,
              col = gray(seq(0.8,0.1,l = dim(x$samples)[3])[v],0.5))
      }
      if(include_prior){
        lines(density(x$prior_phi$r(5e3),adjust = 2, from = 0),
              lwd = 2,
              col = gray(0.3),
              lty = 2)
        legend("topright",
               bty = 'n',
               pch = c(15,15,NA),
               pt.cex = 3,
               lty = c(NA,NA,2),
               lwd = c(NA,NA,3),
               col = gray(c(0.8,0.1,0.3)),
               legend = c(x$v_sequence[1],
                          x$v_sequence[length(x$v_sequence)],
                          "Prior"))
      }else{
        legend("topright",
               bty = 'n',
               pch = 15,
               pt.cex = 3,
               col = gray(c(0.8,0.1)),
               legend = c(x$v_sequence[1],
                          x$v_sequence[length(x$v_sequence)]))
      }
    }else{
      plot(density(x$samples[,dim(x$samples)[2],v_index],
                   adjust = 2,
                   from = 0),
           lwd = 3,
           bty = 'l',
           xlab = expression(phi),
           ylab = '',
           yaxt = 'n',
           main = "")
      if(include_prior){
        lines(density(x$prior_phi$r(5e3),adjust = 2, from = 0),
              lwd = 2,
              col = gray(0.3),
              lty = 2)
        legend("topright",
               bty = 'n',
               lwd = 3,
               lty = 2:1,
               pt.cex = 3,
               legend = c("Prior","Posterior"))
      }
    }
    
    
    
  }
  
  if(type == "transition"){
    summ = summary(x)
    
    yl = range(unlist(summ[,c(2:4)]))
    
    vars = unique(summ$Variable)
    p = length(vars)
    summ = summ[order(summ$v),]
    
    plot(0,
         type = 'n',
         bty = 'l',
         ylab = '',
         xlab = '',
         xlim = range(summ$v)*c(0.95,1.05),
         ylim = yl,
         xaxt = 'n')
    axis(1, at = pretty(sort(unique(x$v))))
    
    for(j in 1:p){
      x_subset = summ[which(summ$Variable == vars[j]),]
      polygon(x = c(x_subset$v,rev(x_subset$v)),
              y = c(x_subset$Lower,rev(x_subset$Upper)),
              border = FALSE,
              col = adjustcolor(MetBrewer::met.brewer("Austria",p)[j],0.3))
    }
    for(j in 1:p){
      x_subset = summ[which(summ$Variable == vars[j]),]
      lines(x_subset$Estimate ~ x_subset$v,
            lwd = 2,
            col = MetBrewer::met.brewer("Austria",p)[j])
    }
    
    
  }
  
}