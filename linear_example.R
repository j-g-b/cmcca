#
library(magrittr)
library(tidyverse)
#
n <- 100
p1 <- 5
p2 <- 4
d <- min(p1, p2)
#
Lambda <- c(0.9, 0.7, 0.1, 0.05)#runif(d) %>% sort(decreasing = T)
#
S <- rnorm(p1*d) %>% matrix(nrow = p1)
Q1 <- svd(S)[["u"]]%*%t(svd(S)[["v"]])
S <- rnorm(p2*d) %>% matrix(nrow = p2)
Q2 <- svd(S)[["u"]]%*%t(svd(S)[["v"]])
#
X <- mvtnorm::rmvnorm(n, sigma = rbind(cbind(diag(rep(1, p1)), Q1%*%diag(Lambda)%*%t(Q2)), cbind(Q2%*%diag(Lambda)%*%t(Q1), diag(rep(1, p2)))))
X1 <- X[, 1:p1]
X2 <- X[, (p1+1):(p1+p2)]
#
W1 <- rWishart(1, p1+2, diag(rep(1, p1)))[,,1] %>% solve()
Y1 <- apply(X1, 1, function(x){
  W1%*%x
}) %>%
  t()
W2 <- rWishart(1, p2+2, diag(rep(1, p2)))[,,1] %>% solve()
Y2 <- apply(X2, 1, function(x){
  W2%*%x
}) %>%
  t()
#
sbmcmc_res <- cmcca::cmcca_mcmc(Y1, Y2, iter = 1000, burn_in = 100, thin = 2)
#
L_SBMCMC <- sapply(1:ncol(sbmcmc_res$Lambda), function(s){
  Q1 <- sbmcmc_res$Q1[,,s]
  Q2 <- sbmcmc_res$Q2[,,s]
  Lambda <- diag(sbmcmc_res$Lambda[,s])
  A <- Q1%*%Lambda%*%t(Q2)
  #return(c(A))
  svd(A)[["d"]]
})
#
X <- cbind(apply(Y1, 2, function(x){x - mean(x)}), apply(Y2, 2, function(x){x - mean(x)}))
A <- cmcca::compute_relation_matrix(t(X)%*%X, c(p1, p2)) %>%
  cmcca::extract_association_matrices(c(p1, p2)) %>%
  magrittr::extract2(1)
L <- svd(A)[['d']]
#
L_SBMCMC %>% t() %>% boxplot()
points(1:d, L, col = 'green')
points(1:d, Lambda, col = 'blue')

