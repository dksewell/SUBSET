#' SUBspace Shrinkage via Exponential Tilting for estimated 
#' projection matrix
#' 
#' plot.subset_gibbs will provide several plotting functions 
#' for an object returned by SUBSET_gibbs() or SUBSET_IS_gibbs().  
#' Obtain traceplots, cumulative mean plots, histograms of 
#' the model parameters, and histogram of the 
#' projection matrix parameter.
#' 
#' @param x Object of class "subset_gibbs", returned by 
#' the function SUBSET_gibbs().
#' @param type character denoting traceplots ("trace"), 
#' cumulative mean plots ("cummean"), or histograms 
#' of the model parameters ("theta").
#' @param burnin The number of initial samples to disregard before 
#' plotting histograms
#' @param include_prior Superimpose the marginal prior 
#' densities on the histograms if type = "theta".
#' @param plotting_args List of additional graphical parameters for primary plot.
#' @param legend_args List of additional graphical parameters for legend.
#' @param prior_args List of additional graphical parameters for superimposed prior plot.
#' @param prior_density_args List of arguments to be passed to density() for density estimation.
#' @export plot.subset_gibbs
#' @export


plot.subset_gibbs = function(x,
                             type = c("trace","cummean","theta")[1],
                             burnin,
                             include_prior = TRUE,
                             plotting_args,
                             legend_args,
                             prior_args,
                             prior_density_args){
  
  if(include_prior){
    if("subset_asymp_gibbs" %in% class(x)){
      prior_draws = x$prior_draws
    }else{
      prior_draws = x$draws0$prior
    }
  }
  
  if(type == "trace"){
    
    matplot_args = 
      list(type = 'l',
           bty = 'l',
           lty = 1,
           lwd = 2,
           xlab = "Gibbs iteration",
           col = gray(0.5,0.5))
    if(!missing(plotting_args)){
      for(j in names(plotting_args)) matplot_args[[j]] = plotting_args[[j]]
    }
    for(j in 1:ncol(x$samples)){
      matplot_args$x = 
        x$samples[,j]
      matplot_args$ylab = 
        colnames(x$samples)[j]
      
      do.call(matplot,matplot_args)
      
      if(j < ncol(x$samples)){
        cat("\nHit RETURN for next plot\n")
        readline()
      }
    }
    
  }
  
  if(type == "cummean"){
    
    cummean = function(z) cumsum(z) / c(1:NROW(z))
    cummeans = 
      apply(x$samples,2,cummean)
    cummeans = cummeans[,order(cummeans[nrow(cummeans),],decreasing=T)]
    margin = par()$mar
    layout(matrix(1:2,1,2),widths = c(3,1))
    par(mar = c(5.1,4.1,1,0.5))
    
    matplot_args = 
      list(x = cummeans,
           type = 'l',
           bty = 'l',
           lty = 1,
           lwd = 2,
           col = rainbow(ncol(cummeans),s=0.75, v = 0.75),
           xlab = "Gibbs iteration",
           ylab = "")
    if(!missing(plotting_args)){
      for(j in names(plotting_args)) matplot_args[[j]] = plotting_args[[j]]
    }
    
    do.call(matplot,matplot_args)
    
    par(mar = rep(0.1,4))
    plot(0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
    matplot_legend_args = 
      list(x = "topleft",
           bty = 'n',
           lwd = 4,
           col = rainbow(ncol(cummeans),s=0.75, v = 0.75),
           legend = colnames(cummeans))
    if(!missing(legend_args)){
      for(j in names(legend_args)) matplot_legend_args[[j]] = legend_args[[j]]
    }
    do.call(legend,matplot_legend_args)
    layout(matrix(1,1,1))
    par(mar = margin)
    
  }
  
  
  if(type == "theta"){
    
    if(missing(burnin)) burnin = 0
    
    for(j in 1:ncol(x$samples)){
      if( j < ncol(x$samples) ){
        
        hist_args = 
          list(x = x$samples[(burnin+1):nrow(x$samples),j],
               freq = FALSE,
               yaxt = 'n',
               ylab = '',
               main = '',
               xlab = colnames(x$samples)[j])
        
        if(!missing(plotting_args)){
          for(j in names(plotting_args)) hist_args[[j]] = plotting_args[[j]]
        }
        
        do.call(hist,hist_args)
        
        if(include_prior){
          d_args = 
            list(prior_draws[,j],
                 adjust = 2)
          if(!missing(prior_density_args)){
            for(j in names(prior_density_args)) d_args[[j]] = prior_density_args[[j]]
          }
          prior_density = 
            do.call(density,d_args)
          
          line_args = 
            list(x = prior_density,
                 lwd = 2,
                 col = gray(0.3))
          if(!missing(prior_args)){
            for(j in names(prior_args)) line_args[[j]] = prior_args[[j]]
          }
          do.call(lines,line_args)
          
        }
        
        cat("\nHit RETURN for next plot\n")
        readline()
      }else{
        
        if("subset_asymp_gibbs" %in% class(x)){
          
          hist_args = 
            list(x = x$samples[(burnin+1):nrow(x$samples),j],
                 freq = FALSE,
                 yaxt = 'n',
                 ylab = '',
                 main = '',
                 xlab = colnames(x$samples)[j])
          
          if(!missing(plotting_args)){
            for(j in names(plotting_args)) hist_args[[j]] = plotting_args[[j]]
          }
          
          do.call(hist,hist_args)
          
        }else{
          
          if(include_prior){
            margin = par()$mar
            layout(matrix(1:2,1,2),widths = c(3,1))
            par(mar = c(5.1,4.1,1,0.5))
            barplot(t(x$phi_pmf[,2:3]),
                    beside = TRUE,
                    col = gray(c(0.25,0.75)),
                    names.arg = prettyNum(x$phi_pmf[,1],digits = 3,format = "g"))
            par(mar = rep(0.1,4))
            plot(0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
            
            bp_legend_args = 
              list(x = "topleft",
                   bty = 'n',
                   pch = 15,
                   cex = 2,
                   pt.cex = 2,
                   col = gray(c(0.25,0.75)),
                   legend = c("Prior","Posterior"))
            if(!missing(legend_args)){
              for(j in names(legend_args)) bp_legend_args[[j]] = legend_args[[j]]
            }
            do.call(legend,bp_legend_args)
            layout(matrix(1,1,1))
            par(mar = margin)
          }else{
            barplot(x$phi_pmf[,3],
                    names.arg = prettyNum(x$phi_pmf[,1],digits = 3,format = "g"))
          }
          
        }
        
      }
    }
  }
  
}