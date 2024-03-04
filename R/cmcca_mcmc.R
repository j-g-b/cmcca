#' Function for semiparametric CCA model
#'
#' @export cmcca_mcmc
cmcca_mcmc <- function(Y1, Y2, iter = 1000, burn_in = 0, thin = 1, plot = F, d = NULL, maxn = 100){
  #
  require(magrittr)
  # Specify problem dimensions
  n <- nrow(Y1)
  if(is.null(d)){
    d <- min(ncol(Y1), ncol(Y2))
  }
  p1 <- ncol(Y1)
  p2 <- ncol(Y2)
  ss <- min(maxn, n)
  # Initialize parameter values
  X1 <- rnorm(n*p1) %>% matrix(nrow=n)
  X1 <- sweep(X1, 2, colMeans(X1))
  X2 <- rnorm(n*p2) %>% matrix(nrow=n)
  X2 <- sweep(X2, 2, colMeans(X2))
  assign1 <- transport::transport(transport::pp(X1),
                                  transport::pp(Y1),
                                  p = 2)
  assign2 <- transport::transport(transport::pp(X2),
                                  transport::pp(Y2),
                                  p = 2)
  X1 <- X1[order(assign1[, 2]), ]
  X2 <- X2[order(assign2[, 2]), ]
  # Initialize parameter values
  aX <- (t(cbind(X1, X2))%*%cbind(X1, X2) / n) %>%
    cmcca::extract_association_matrices(d=c(p1, p2)) %>%
    magrittr::extract2(1)
  Q1 <- svd(aX)[["u"]][, 1:d]
  Q2 <- svd(aX)[["v"]][, 1:d]
  W1 <- Q1
  W2 <- Q2
  mY <- colMeans(cbind(Y1, Y2))
  cY <- sweep(cbind(Y1, Y2), 2, mY)
  #
  Lambda <- svd(aX)[["d"]][1:d] %>%
    pmin(rep(0.99, d)) %>%
    pmax(rep(0.01, d)) %>%
    magrittr::add(rnorm(d, sd = 0.001)) %>%
    sort(decreasing = T)
  # Allocate parameter arrays
  Lambda_samps <- array(dim = c(d, floor((iter - burn_in)/thin)))
  Q1_samps <- array(dim = c(dim(Q1), floor((iter - burn_in)/thin)))
  Q2_samps <- array(dim = c(dim(Q2), floor((iter - burn_in)/thin)))
  #
  X1_mean <- matrix(0, nrow = nrow(X1), ncol = ncol(X1))
  X2_mean <- matrix(0, nrow = nrow(X2), ncol = ncol(X2))
  #
  X1_samps <- array(dim = c(dim(X1), floor((iter - burn_in)/thin)))
  X2_samps <- array(dim = c(dim(X2), floor((iter - burn_in)/thin)))
  # Run Markov Chain
  xlim <- c(min(c(range(Y1[, 1]), range(X1[,1]))), max(c(range(Y1[, 1]), range(X1[,1]))))
  ylim <- c(min(c(range(Y1[, 2]), range(X1[,2]))), max(c(range(Y1[, 2]), range(X1[,2]))))
  # Time it
  start_time <- Sys.time()
  for(s in 1:iter){
    if(s%%10 == 0){
      print(paste0(round(100*s / iter, 2), "% done..."))
    }
    #
    Lambda <- cmcca::sample_Lambda_cubic(Lambda, X1, X2, Q1, Q2)
    # Sample new parameter values
    #assign1 <- transport::transport(transport::pp(X1),
    #                                transport::pp(Y1),
    #                                p = 2)
    #if(!all(assign1[, 2] == 1:n)){
    #  View(assign1)
    #  stop()
    #}
    Sub <- sample(1:n, ss) - 1
    X1 <- cmcca::sample_X(X1, X2, Q1, Q2, Lambda, Y1, Sub)
    if(plot){
      #png(filename = paste0("cmmc_ani/p", s, ".png"), width = 8, height = 4, units = 'in', res = 300)
      #par(mfrow = c(1, 2))
      #plot(X1, xlim = xlim, ylim = ylim, col = NA, ann = F)
      #segments(X1[, 1], X1[, 2], Y1[, 1], Y1[, 2], lty = 3, col = 'blue')
      #points(Y1, col = 'black', pch=16)
      #points(X1, xlim = xlim, ylim = ylim, col = 'grey', ann = F)
      #plot(X1, xlim = xlim, ylim = ylim, col = 'grey', ann = F)
      #Sys.sleep(0.09)
      #dev.off()
      if(s==1){
        plot(X1, xlim = xlim, ylim = ylim, col = c('red', rep('black', n-1)), ann = F)
      } else {
        points(x = X1[1, 1], y = X1[1, 2], pch = 16, cex=0.2, xlim = xlim, ylim = ylim, col = 'red', ann = F)
        if(s%%10 == 0){
          Sys.sleep(0.09)
        }
      }
    }
    #
    Sub <- sample(1:n, ss) - 1
    X2 <- cmcca::sample_X(X2, X1, Q2, Q1, Lambda, Y2, Sub)
    #
    # Generalized Gibbs step
    IL <- 1/((1/Lambda^2)-1)
    IL2 <- 1/(Lambda - (1/Lambda))
    b <- sum(diag(t(X1)%*%X1)) +
         sum(diag(t(X2)%*%X2)) +
         sum(diag(diag(IL)%*%t(X1%*%Q1)%*%X1%*%Q1)) +
         sum(diag(diag(IL)%*%t(X2%*%Q2)%*%X2%*%Q2)) +
         sum(diag(diag(IL2)%*%t(X2%*%Q2)%*%X1%*%Q1)) +
         sum(diag(diag(IL2)%*%t(X1%*%Q1)%*%X2%*%Q2))
    u <- rgamma(1, (n*(p1+p2)) / 2, b / 2)
    g <- sqrt(u)
    X1 <- g*X1
    X2 <- g*X2
    #
    # Generalized Gibbs step
    Cmat <- rbind(cbind(diag(rep(1, p1)), Q1%*%diag(Lambda)%*%t(Q2)),
                  cbind(Q2%*%diag(Lambda)%*%t(Q1), diag(rep(1, p2))))
    Mean <- colMeans(cbind(X1, X2))
    g <- Mean + t(chol(Cmat/n))%*%rnorm(p1+p2)
    X1 <- sweep(X1, 2, g[1:p1], "-")
    X2 <- sweep(X2, 2, g[(p1+1):(p1+p2)], "-")
    #
    # Generalized Gibbs step
    Cmatinv <- solve(Cmat)
    iv <- sum(diag(Cmatinv%*%crossprod(cY)))
    mn <- sum(diag(Cmatinv%*%crossprod(cbind(X1, X2), cY)))
    print(mn/iv)
    if(mn < 0){
      g <- truncdist::rtrunc(1, "norm", b = 0, mean = mn/iv, sd = 1/sqrt(iv))
      X1 <- X1 - g*cY[, 1:p1]
      X2 <- X2 - g*cY[, (p1+1):(p1+p2)]
      print("G Gibbs")
    }
    #
    Q1 <- cmcca::sample_rbmf_slice_Q(X1, X2, Q1, Q2, Lambda, W1)
    W1 <- Q1[["W"]]
    Q1 <- Q1[["Q"]]
    #
    Q2 <- cmcca::sample_rbmf_slice_Q(X2, X1, Q2, Q1, Lambda, W2)
    W2 <- Q2[["W"]]
    Q2 <- Q2[["Q"]]
    #
    if(s > burn_in & (s-burn_in)%%thin == 0){
      # Store needed parameter values
      Lambda_samps[, (s-burn_in)/thin] <- Lambda
      Q1_samps[, , (s-burn_in)/thin] <- Q1
      Q2_samps[, , (s-burn_in)/thin] <- Q2
      X1_mean <- X1_mean + X1/floor((iter - burn_in)/thin)
      X2_mean <- X2_mean + X2/floor((iter - burn_in)/thin)
      X1_samps[, , (s-burn_in)/thin] <- X1
      X2_samps[, , (s-burn_in)/thin] <- X2
    }
  }
  end_time <- Sys.time()
  #
  return(list(Lambda = Lambda_samps, Q1 = Q1_samps, Q2 = Q2_samps,
              X1 = X1_mean, X2 = X2_mean, X1_samps = X1_samps, X2_samps = X2_samps,
              iter_time = (end_time - start_time)/iter))
}
