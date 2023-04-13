#'
#'
#' @export slow_sample_X
slow_sample_X <- function(X1, X2, Q1, Q2, Lambda, Y){
  #
  n <- nrow(X1)
  p <- ncol(X1)
  #
  A <- Q1%*%diag(1/((1/Lambda^2)-1))%*%t(Q1)
  B <- t(Q1%*%diag(1/(Lambda - (1/Lambda)))%*%t(X2%*%Q2))
  #
  Xstar <- matrix(rnorm(n*p), nrow = n)
  assign <- transport::transport(transport::pp(Xstar),
                                  transport::pp(Y),
                                  p = 2)
  Xstar <- Xstar[order(assign[, 2]), ]
  #
  log_r <- (-0.5*sum(diag(A%*%t(Xstar)%*%Xstar)) - sum(diag(t(B)%*%Xstar))) -
           (-0.5*sum(diag(A%*%t(X1)%*%X1)) - sum(diag(t(B)%*%X1)))
  log_r2 <- cmcca::joint_log_lik(Xstar, X2, Q1, Q2, Lambda) -
            cmcca::joint_log_lik(X1, X2, Q1, Q2, Lambda) +
            sum(dnorm(X1, log = T)) -
            sum(dnorm(Xstar, log = T))
  stopifnot(abs(log_r - log_r2) < 1e-8)
  if(runif(1) < exp(log_r)){
    return(Xstar)
  } else {
    return(X1)
  }

}

#'
#'
#' @export cov_sample_X
cov_sample_X <- function(X1, X2, Q1, Q2, Lambda, Y){
  #
  n <- nrow(X1)
  p <- ncol(X1)
  ellipse_prob <- 0.9
  #
  A <- Q1%*%diag(1/((1/Lambda^2)-1))%*%t(Q1)
  B <- t(Q1%*%diag(1/(Lambda - (1/Lambda)))%*%t(X2%*%Q2))
  S <- X1%*%t(Y)
  S <- sweep(-S, 1, diag(S),"+")
  W <- rep(0, n)
  #
  for(i in 1:n){
    D <- diag(sapply(1:n, function(j){if(j==i){return(1)}else{return(1/(t(X1[i, ])%*%(Y[i, ] - Y[j, ]) - t(X1[j, ])%*%(Y[i, ] - Y[j, ]))^2)}}))
    SigmaInv <- t(rep(1,n)%*%t(Y[i,]) - Y)%*%D%*%(rep(1,n)%*%t(Y[i,]) - Y)
    Sigma <- MASS::ginv(SigmaInv)/qgamma((ellipse_prob), p/2, 1/2)
    L <- t(chol(Sigma))
    Xprop <- X1[i, ] + L%*%rnorm(p)
    Xstar <- X1
    Xstar[i, ] <- Xprop
    D1 <- diag(sapply(1:n, function(j){if(j==i){return(1)}else{return(1/(t(Xstar[i, ])%*%(Y[i, ] - Y[j, ]) - t(Xstar[j, ])%*%(Y[i, ] - Y[j, ]))^2)}}))
    SigmaInv1 <- t(rep(1,n)%*%t(Y[i,]) - Y)%*%D1%*%(rep(1,n)%*%t(Y[i,]) - Y)
    assign1 <- transport::transport(transport::pp(Xstar),
                                    transport::pp(Y),
                                    p = 2)
    if(!(all(assign1[, 2] == 1:nrow(Y)) == cmcca:::is_cm(Xstar, Y))){
      print(cmcca:::is_cm(Xstar, Y))
    }
    if(cmcca:::is_cm(Xstar, Y)){
      log_r <- (-0.5*t(Xstar[i, ])%*%A%*%Xstar[i, ] - t(B[i, ])%*%Xstar[i, ]) - (-0.5*t(X1[i, ])%*%A%*%X1[i, ] - t(B[i, ])%*%X1[i, ]) +
               sum(dnorm(Xstar[i, ], log = T)) - sum(dnorm(X1[i, ], log = T)) +
               (-0.5*qgamma((ellipse_prob), p/2, 1/2)*t(X1[i, ] - Xstar[i, ])%*%SigmaInv1%*%(X1[i, ] - Xstar[i, ]) + 0.5*log(det(qgamma((ellipse_prob), p/2, 1/2)*SigmaInv1))) - (-0.5*qgamma((ellipse_prob), p/2, 1/2)*t(Xstar[i, ] - X1[i, ])%*%SigmaInv%*%(Xstar[i, ] - X1[i, ]) + 0.5*log(det(qgamma((ellipse_prob), p/2, 1/2)*SigmaInv)))
      #print(log_r)
      if(runif(1) < exp(log_r)){
        X1[i, ] <- Xstar[i, ]
      }
    }
  }
  return(X1)
}

#'
#'
#' @export sample_X
sample_X <- function(X1, X2, Q1, Q2, Lambda, Y){
  #
  n <- nrow(X1)
  p <- ncol(X1)
  ellipse_prob <- 0.9
  #
  A <- Q1%*%diag(1/((1/Lambda^2)-1))%*%t(Q1)
  B <- t(Q1%*%diag(1/(Lambda - (1/Lambda)))%*%t(X2%*%Q2))
  S <- X1%*%t(Y)
  S <- sweep(-S, 1, diag(S),"+")
  W <- rep(0, n)
  #
  result <- cmcca:::eigen_sample_X(Y, X1, W, S, A, B)
  return(result$X)
}

#'
#'
#' @export is_cm
is_cm <- function(X, Y){
  XY<- X%*%t(Y)
  XYD<-sweep(-XY,1,diag(XY),"+")
  n<-nrow(X)
  w<-rep(0,n)
  count<-0
  obj <- min(XYD - outer(w,w,"-"))
  w_prev <- w
  while(obj < -1e-12)
  {
    for(i in 1:n){ if(w[i] < max( -XYD[,i] + w  )){w[i]<-max( -XYD[,i] + w  )} }
    obj_new <- min(XYD - outer(w,w,"-"))
    if(abs(obj_new - obj) < 1e-12){
      count<-count+1
    } else {
      count<-0
    }
    if(count > 10 & obj_new < -1e-12){
      return(F)
    }
    obj <- obj_new
    #plot(w, w_prev)
    #abline(0, 1)
    #Sys.sleep(0.1)
    w_prev <- w
  }
  return(T)
}
