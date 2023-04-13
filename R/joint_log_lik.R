#'
#'
#'
#' @export joint_log_lik
joint_log_lik <- function(X1, X2, Q1, Q2, Lambda){
  n <- nrow(X1)
  IL <- 1/((1/Lambda^2)-1)
  IL2 <- 1/(Lambda - (1/Lambda))
  -((n/2)*sum(log(1 - Lambda^2)) +
      (1/2)*sum(diag(t(X1)%*%X1)) +
      (1/2)*sum(diag(t(X2)%*%X2)) +
      (1/2)*sum(diag(diag(IL)%*%t(X1%*%Q1)%*%X1%*%Q1)) +
      (1/2)*sum(diag(diag(IL)%*%t(X2%*%Q2)%*%X2%*%Q2)) +
      (1/2)*sum(diag(diag(IL2)%*%t(X2%*%Q2)%*%X1%*%Q1)) +
      (1/2)*sum(diag(diag(IL2)%*%t(X1%*%Q1)%*%X2%*%Q2)) +
      n*(ncol(X1)+ncol(X2))*log(2*pi)/2)
}