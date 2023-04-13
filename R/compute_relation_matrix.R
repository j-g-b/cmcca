#'
#'
#' @export compute_relation_matrix
compute_relation_matrix <- function(V, d, return_D = F){
  #
  D <- matrix(0, sum(d), sum(d))
  for(j in 1:length(d)){
    if(j == 1){
      t <- 1
    } else {
      t <- sum(d[1:(j-1)]) + 1
    }
    svd <- svd(V[t:sum(d[1:j]), t:sum(d[1:j])])
    D[t:sum(d[1:j]), t:sum(d[1:j])] <- svd$u%*%diag(sqrt(1/svd$d))%*%t(svd$u)
  }
  if(return_D){
    return(list(R = D%*%V%*%D, D = D))
  } else {
    return(D%*%V%*%D)
  }
}
