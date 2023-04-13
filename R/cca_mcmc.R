#' Template function for MCMC
#'
#' @export cca_mcmc
cca_mcmc <- function(X1, X2, iter = 1000, burn_in = 0, thin = 1, d = NULL){
  #
  require(magrittr)
  # Specify problem dimensions
  n <- nrow(X1)
  if(is.null(d)){
    d <- min(ncol(X1), ncol(X2))
  }
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  # Initialize parameter values
  aX <- (t(cbind(X1, X2))%*%cbind(X1, X2) / n) %>%
    sbcca::compute_relation_matrix(d=c(p1, p2)) %>%
    sbcca::extract_association_matrices(d=c(p1, p2)) %>%
    magrittr::extract2(1)
  Q1 <- svd(aX)[["u"]][, 1:d]
  Q2 <- svd(aX)[["v"]][, 1:d]
  W1 <- Q1
  W2 <- Q2
  #
  Lambda <- svd(aX)[["d"]][1:d] %>%
    pmin(rep(0.99, d)) %>%
    pmax(rep(0.01, d)) %>%
    magrittr::add(rnorm(d, sd = 0.001)) %>%
    sort(decreasing = T)
  # Allocate parameter arrays
  num_samps <- 0
  Lambda_samps <- array(dim = c(d, floor((iter - burn_in)/thin)))
  Q1_samps <- array(dim = c(dim(Q1), floor((iter - burn_in)/thin)))
  Q2_samps <- array(dim = c(dim(Q2), floor((iter - burn_in)/thin)))
  W1_samps <- array(dim = c(dim(W1), floor((iter - burn_in)/thin)))
  W2_samps <- array(dim = c(dim(W2), floor((iter - burn_in)/thin)))
  # Run Markov Chain
  for(s in 1:iter){
    if(s%%100 == 0){
      print(paste0(round(100*s / iter, 2), "% done..."))
    }
    #
    Q1 <- sbcca::sample_rbmf_stiefel_Q(X1, X2, Q1, Q2, Lambda, W1)
    W1 <- Q1[["W"]]
    Q1 <- Q1[["Q"]]
    #
    Q2 <- sbcca::sample_rbmf_stiefel_Q(X2, X1, Q2, Q1, Lambda, W2)
    W2 <- Q2[["W"]]
    Q2 <- Q2[["Q"]]
    #
    Lambda <- sbcca::sample_Lambda(Lambda, X1, X2, Q1, Q2)
    #
    if(s > burn_in & (s-burn_in)%%thin == 0){
      # Store needed parameter values
      Lambda_samps[, (s-burn_in)/thin] <- Lambda
      Q1_samps[, , (s-burn_in)/thin] <- Q1
      Q2_samps[, , (s-burn_in)/thin] <- Q2
      W1_samps[, , (s-burn_in)/thin] <- W1
      W2_samps[, , (s-burn_in)/thin] <- W2
      num_samps <- num_samps + 1
    }
  }
  return(list(Lambda = Lambda_samps, Q1 = Q1_samps, Q2 = Q2_samps,
              W1 = W1_samps, W2 = W2_samps))
}
