#'
#'
#' @export extract_association_matrices
extract_association_matrices <- function(R, d){
  m <- length(d)
  association_matrix_list <- list()
  for(j in 1:(m-1)){
    for(j_prime in (j+1):m){
      if(j == 1){
        k <- 1
        k_prime <- sum(d[1:(j_prime-1)]) + 1
      } else {
        k <- sum(d[1:(j-1)]) + 1
        k_prime <- sum(d[1:(j_prime-1)]) + 1
      }
      association_matrix_list[[paste0("(", j, ", ", j_prime, ")")]] <- R[k:(k+d[j]-1), k_prime:(k_prime+d[j_prime]-1)]
    }
  }
  return(association_matrix_list)
}
