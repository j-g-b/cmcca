#'
#'
#' @export sample_Lambda_cubic
sample_Lambda_cubic <- function(Lambda, X1, X2, Q1, Q2, plot = F){
  n <- nrow(X1)
  d <- ncol(Q1)
  #
  a <- diag(t(X1%*%Q1)%*%X1%*%Q1)
  b <- diag(t(X2%*%Q2)%*%X2%*%Q2)
  c <- diag(t(X2%*%Q2)%*%X1%*%Q1)
  #
  cubic_modes <- sapply(1:d, function(k){
    params <- c(1, -c[k]/n, (a[k] + b[k] - n)/n, -c[k]/n)
    real_root <- RConics::cubic(params) %>% magrittr::extract(1) %>% Re()
    return(real_root)
  })
  #
  rad <- sqrt(((a+b)/2 - abs(c)) / ((a+b)/2)) / sqrt(n)

  #
  for(k in 1:d){
    #
    Lambda_new <- Lambda
    if(k == 1){
      Lambda_new[k] <- truncdist::rtrunc(1, spec = "norm", a = Lambda[k+1], b = 1,
                                         mean = cubic_modes[k], sd = rad[k])
      #
      log_r <- sum(cmcca::Lambda_log_lik(Lambda_new, a, b, c, n)) -
               sum(cmcca::Lambda_log_lik(Lambda, a, b, c, n)) +
               truncdist::dtrunc(Lambda[k], spec = "norm", a = Lambda[k+1], b = 1,
                                 mean = cubic_modes[k], sd = rad[k], log = T) -
               truncdist::dtrunc(Lambda_new[k], spec = "norm", a = Lambda[k+1], b = 1,
                                 mean = cubic_modes[k], sd = rad[k], log = T)
    } else if(k == d){
      Lambda_new[k] <- truncdist::rtrunc(1, spec = "norm", a = 0, b = Lambda[k-1],
                                         mean = cubic_modes[k], sd = rad[k])
      #
      log_r <- sum(cmcca::Lambda_log_lik(Lambda_new, a, b, c, n)) -
               sum(cmcca::Lambda_log_lik(Lambda, a, b, c, n)) +
               truncdist::dtrunc(Lambda[k], spec = "norm", a = 0, b = Lambda[k-1],
                                 mean = cubic_modes[k], sd = rad[k], log = T) -
               truncdist::dtrunc(Lambda_new[k], spec = "norm", a = 0, b = Lambda[k-1],
                                 mean = cubic_modes[k], sd = rad[k], log = T)
    } else {
      Lambda_new[k] <- truncdist::rtrunc(1, spec = "norm", a = Lambda[k+1], b = Lambda[k-1],
                                         mean = cubic_modes[k], sd = rad[k])
      #
      log_r <- sum(cmcca::Lambda_log_lik(Lambda_new, a, b, c, n)) -
               sum(cmcca::Lambda_log_lik(Lambda, a, b, c, n)) +
               truncdist::dtrunc(Lambda[k], spec = "norm", a = Lambda[k+1], b = Lambda[k-1],
                                 mean = cubic_modes[k], sd = rad[k], log = T) -
               truncdist::dtrunc(Lambda_new[k], spec = "norm", a = Lambda[k+1], b = Lambda[k-1],
                                 mean = cubic_modes[k], sd = rad[k], log = T)
    }
    if(runif(1) < exp(log_r)){
      Lambda[k] <- Lambda_new[k]
    }
  }
  #
  if(plot){
    print(c)
    print((a+b)/2)
    l <- seq(-1, 1, length.out = 1000)
    plot(l, cmcca::Lambda_log_lik(l, a[1], b[1], c[1], n) %>% exp(), type = 'l')
    segments(x0 = cubic_modes[1] + 2*rad[1], x = cubic_modes[1] - 2*rad[1], y0 = 0, y = 0, col = 'red')
    Sys.sleep(0.2)
  }

  return(Lambda)
}
#'
#'
#' @export sample_Lambda
sample_Lambda <- function(Lambda, X1, X2, Q1, Q2){
  n <- nrow(X1)
  d <- ncol(Q1)
  conc <- 100
  #
  a <- diag(t(X1%*%Q1)%*%X1%*%Q1)
  b <- diag(t(X2%*%Q2)%*%X2%*%Q2)
  c <- diag(t(X2%*%Q2)%*%X1%*%Q1)
  #
  for(k in 1:length(Lambda)){
    if(k == 1){
      u_bound <- 1
      l_bound <- Lambda[k+1]
      lam_curr <- (Lambda[k] - l_bound) / (u_bound - l_bound)
      lam_prop <- rbeta(1, conc*lam_curr, conc*(1 - lam_curr))
      Lambda_prop <- lam_prop*(u_bound - l_bound) + l_bound
      if(Lambda_prop != 1){
        log_r <- cmcca::Lambda_log_lik(Lambda_prop, a[k], b[k], c[k], n) -
                cmcca::Lambda_log_lik(Lambda[k], a[k], b[k], c[k], n) +
                dbeta(lam_curr, conc*lam_prop, conc*(1 - lam_prop), log = T) -
                dbeta(lam_prop, conc*lam_curr, conc*(1 - lam_curr), log = T)
        if(runif(1) < exp(log_r)){
          Lambda[k] <- Lambda_prop
        }
      }
    } else if(k > 1 & k < d){
      u_bound <- Lambda[k-1]
      l_bound <- Lambda[k+1]
      lam_curr <- (Lambda[k] - l_bound) / (u_bound - l_bound)
      lam_prop <- rbeta(1, conc*lam_curr, conc*(1 - lam_curr))
      Lambda_prop <- lam_prop*(u_bound - l_bound) + l_bound
      log_r <- cmcca::Lambda_log_lik(Lambda_prop, a[k], b[k], c[k], n) -
        cmcca::Lambda_log_lik(Lambda[k], a[k], b[k], c[k], n) +
        dbeta(lam_curr, conc*lam_prop, conc*(1 - lam_prop), log = T) -
        dbeta(lam_prop, conc*lam_curr, conc*(1 - lam_curr), log = T)
      if(runif(1) < exp(log_r)){
        Lambda[k] <- Lambda_prop
      }
    } else {
      u_bound <- Lambda[k-1]
      l_bound <- 0
      lam_curr <- (Lambda[k] - l_bound) / (u_bound - l_bound)
      lam_prop <- rbeta(1, conc*lam_curr, conc*(1 - lam_curr))
      Lambda_prop <- lam_prop*(u_bound - l_bound) + l_bound
      if(Lambda_prop != 0){
        log_r <- cmcca::Lambda_log_lik(Lambda_prop, a[k], b[k], c[k], n) -
          cmcca::Lambda_log_lik(Lambda[k], a[k], b[k], c[k], n) +
          dbeta(lam_curr, conc*lam_prop, conc*(1 - lam_prop), log = T) -
          dbeta(lam_prop, conc*lam_curr, conc*(1 - lam_curr), log = T)
        if(runif(1) < exp(log_r)){
          Lambda[k] <- Lambda_prop
        }
      }
    }
  }
  return(Lambda)
}
#'
#'
#' @export sample_unrestricted_Lambda
sample_unrestricted_Lambda <- function(Lambda, X1, X2, Q1, Q2){
  #
  n <- nrow(X1)
  d <- ncol(Q1)
  conc <- 40
  #
  a <- diag(t(X1%*%Q1)%*%X1%*%Q1)
  b <- diag(t(X2%*%Q2)%*%X2%*%Q2)
  c <- diag(t(X2%*%Q2)%*%X1%*%Q1)
  #
  for(k in 1:length(Lambda)){
    lambda_prop <- truncdist::rtrunc(1, spec="norm", a=0, b=1, mean = Lambda[k], sd = 1/conc)
    log_r <- cmcca::Lambda_log_lik(lambda_prop, a[k], b[k], c[k], n) -
             cmcca::Lambda_log_lik(Lambda[k], a[k], b[k], c[k], n) +
             truncdist::dtrunc(Lambda[k], "norm", a=0, b=1, mean = lambda_prop, sd = 1/conc, log = T) -
             truncdist::dtrunc(lambda_prop, "norm", a=0, b=1, mean = Lambda[k], sd = 1/conc, log = T)
    if(runif(1) < exp(log_r)){
      Lambda[k] <- lambda_prop
    }
  }
  return(Lambda)
}
#'
#'
#' @export Lambda_log_lik
Lambda_log_lik <- function(Lambda, a, b, c, n){
 -(n/2)*log(1 + Lambda) - (n/2)*log(1 - Lambda) - ((a+b)/2 + c)/(2*(1+Lambda)) - ((a+b)/2 - c)/(2*(1-Lambda))
}
