#'
#'
#' @export sample_rbmf_stiefel_Q
sample_rbmf_stiefel_Q <- function(X1, X2, Q1, Q2, Lambda, W){
  #
  n <- nrow(X1)
  p <- ncol(X1)
  d <- length(Lambda)
  #
  conc <- n/3/p
  #
  A <- -(t(X1)%*%X1)/2
  B <- diag(1/(1/Lambda^2 - 1))
  C <- t(X1)%*%X2%*%Q2%*%diag(1/(1/Lambda - Lambda))
  #
  Q_new <- rstiefel::rmf.matrix(conc*Q1)
  log_r <- cmcca::Q_log_lik(Q_new, A, B, C) -
           cmcca::Q_log_lik(Q1, A, B, C)
  if(runif(1) < exp(log_r)){
    return(list(Q = Q_new, W = W))
  } else {
    return(list(Q = Q1, W = W))
  }
}
#'
#'
#' @export sample_rbmf_slice_Q
sample_rbmf_slice_Q <- function(X1, X2, Q1, Q2, Lambda, W){
  #
  n <- nrow(X1)
  p <- ncol(X1)
  d <- length(Lambda)
  #
  nu <- rnorm(p*d) %>% matrix(nrow = p)
  u <- runif(1)
  #
  A <- -(t(X1)%*%X1)/2
  B <- diag(1/(1/Lambda^2 - 1))
  C <- t(X1)%*%X2%*%Q2%*%diag(1/(1/Lambda - Lambda))
  #
  ly <- cmcca::Q_log_lik(Q1, A, B, C) + log(u)
  theta <- runif(1, 0, 2*pi)
  theta_min_max <- c(theta - 2*pi, theta)
  while(T){
    W_new <- W*cos(theta) + nu*sin(theta)
    svd_W_new <- svd(W_new)
    Q_new <- svd_W_new[["u"]]%*%t(svd_W_new[["v"]])
    l_new <- cmcca::Q_log_lik(Q_new, A, B, C)
    if(l_new > ly){
      return(list(Q = Q_new, W = W_new))
    } else {
      if(theta < 0){
        theta_min_max[1] <- theta
      } else {
        theta_min_max[2] <- theta
      }
      theta <- runif(1, theta_min_max[1], theta_min_max[2])
    }
  }
}
#'
#'
#' @export Q_log_lik
Q_log_lik <- function(Q, A, B, C){
  sum(diag(t(C)%*%Q + B%*%t(Q)%*%A%*%Q))
}
