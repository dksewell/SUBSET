#' Create a projection matrix
#' 
#' @param x Either a matrix/numeric vector providing the column span of 
#' the projection matrix P, or else a list of such matrices/vecors.  If 
#' the latter, a block-diagonal projection matrix will be constructed.
#' 
#' @return A projection matrix.
#' 
#' @examples 
#' Proj(rep(1,3))
#' Proj(cbind(1,1:3))
#' Proj(list(diag(2),cbind(1,1:3),rep(1,4)))


Proj = function(x){
  proj_inner = function(y){
    if("matrix" %in% class(y)){
      return( tcrossprod(y %*% qr.solve(crossprod(y)),
                         y)
      )
    }else{
      return( tcrossprod(y, y) / drop(crossprod(y)) )
    }
  }
  if(isTRUE(class(x) ==  "list")){
    return(Matrix::as.matrix(Matrix::bdiag(lapply(x,proj_inner))))
  }else{
    return(proj_inner(x))
  }
}
