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
#' cumulative mean plots ("cummean"), histograms 
#' of the model parameters ("theta"), or a histogram 
#' of the projection matrix parameter ("phi"), or 
#' transition plot of the theta values over v ("transition").
#' @param v The desired value of v, i.e., the exponential tilting parameter. 
#' Only used for type = "theta" or "phi".
#' @param include_prior Superimpose the marginal prior 
#' densities on the histograms if type = "theta" or "phi".
#' @param ... Additional graphical parameters.
#' @export



plot.subset_is_gibbs = function(x,
                                type = c("trace","cummean","theta","phi","transition")[1],
                                v,
                                include_prior = TRUE,
                                ...){
  
  if( (class(x) == "list") & (type="transition") )stop("Only type='transition' is compatible with a list of subset_is_gibbs objects.")
  
  
  if(type == "trace"){
    for(j in 1:ncol(x$samples)){
      matplot(x$samples[,j],
              type = 'l',
              bty = 'l',
              lty = 1,
              lwd = 2,
              col = gray(0.5,0.5),
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
    cummean = function(z) cumsum(z) / c(1:NROW(z))
    cummeans = 
      apply(x$samples,2,cummean)
    cummeans = cummeans[,order(cummeans[nrow(cummeans),],decreasing=T)]
    margin = par()$mar
    layout(matrix(1:2,1,2),widths = c(3,1))
    par(mar = c(5.1,4.1,1,0.5))
    matplot(cummeans,
            type = 'l',
            bty = 'l',
            lty = 1,
            lwd = 2,
            col = rainbow(ncol(cummeans),s=0.75, v = 0.75),
            xlab = "Gibbs iteration",
            ylab = "",
            ...)
    par(mar = rep(0.1,4))
    plot(0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
    legend("topleft",
           bty = 'n',
           lwd = 4,
           col = rainbow(ncol(cummeans),s=0.75, v = 0.75),
           legend = colnames(cummeans))
    layout(matrix(1,1,1))
    par(mar = margin)
  }
  
  
  if(type == "theta"){
    for(j in 1:ncol(x$samples)){
      if( j < ncol(x$samples) ){
        hist(x$samples[,j],
           freq = FALSE,
           yaxt = 'n',
           # bty = 'l',
           ylab = '',
           main = '',
           xlab = colnames(x$samples)[j],
           ...)
        if(include_prior){
          lines(density(x$draws0$prior[,j],adjust = 2),
                lwd = 2,
                col = gray(0.3),
                ...)
        }
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
          legend("topleft",
                 bty = 'n',
                 pch = 15,
                 cex = 2,
                 pt.cex = 2,
                 col = gray(c(0.25,0.75)),
                 legend = c("Prior","Posterior"))
          layout(matrix(1,1,1))
          par(mar = margin)
        }else{
          barplot(x$phi_pmf[,3],
                  names.arg = prettyNum(x$phi_pmf[,1],digits = 3,format = "g"))
        }
      }
      if(j < ncol(x$samples)){
        cat("\nHit RETURN for next plot\n")
        readline()
      }
    }
  }
  
  
  if(type == "transition"){
    
    summ = 
      do.call(rbind,lapply(x,summary))
    
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
    axis(1, at = pretty(sort(unique(summ$v))))
    
    for(j in 1:p){
      x_subset = summ[which(summ$Variable == vars[j]),]
      polygon(x = c(x_subset$v,rev(x_subset$v)),
              y = c(x_subset[,3],rev(x_subset[,4])),
              border = FALSE,
              col = adjustcolor(rainbow(p,s=0.75, v = 0.75)[j],0.3))
    }
    for(j in 1:p){
      x_subset = summ[which(summ$Variable == vars[j]),]
      lines(x_subset$Estimate ~ x_subset$v,
            lwd = 2,
            col = rainbow(p,s=0.75, v = 0.75)[j])
    }
    
  }
  
}