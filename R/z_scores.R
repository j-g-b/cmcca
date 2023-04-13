#'
#'
#' @export z_scores
z_scores <- function(y, ties.method="average") {
  u <- rank(y, na.last="keep", ties.method=ties.method)/(sum(!is.na(y))+1)
  z <- qnorm(u)
  names(z) <- names(y)
  m <- dim(y)
  if(length(m)==2){ 
    z<-matrix(z,nrow=m[1],ncol=m[2]) 
    dimnames(z)<-dimnames(y)
  }
  if(length(m)>=3){ 
    z<-array(z,dim=m) 
    dimnames(z)<-dimnames(y)
  }
  z
}