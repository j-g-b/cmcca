#
set.seed(1)
#
library(magrittr)
library(tidyverse)
#
N <- c(50, 100, 250, 500, 1000)
n_sim <- 200
p1 <- 2
p2 <- 2
d <- min(p1, p2)
#
COP_MSE <- matrix(NA, length(N), n_sim)
CCA_MSE <- matrix(NA, length(N), n_sim)
PLUG_MSE <- matrix(NA, length(N), n_sim)
#
U <- matrix(rnorm(4, sd = 0.75), nrow = 2)
V <- matrix(rnorm(4, sd = 0.75), nrow = 2)
cm_trans <- function(x, U){
  t(U)%*%(qweibull(pnorm(U%*%x), 1, 1))
}
#
for(i in 1:length(N)){
  n <- N[i]
  for(sim in 1:n_sim){
    rm(.Random.seed, envir=globalenv())
    #
    l1 <- rbeta(1, 2, 1)
    l2 <- runif(1, 0, l1)
    Lambda <- c(l1, l2)
    #
    S <- rnorm(p1*d) %>% matrix(nrow = p1)
    Q1 <- svd(S)[["u"]]%*%t(svd(S)[["v"]])
    S <- rnorm(p2*d) %>% matrix(nrow = p2)
    Q2 <- svd(S)[["u"]]%*%t(svd(S)[["v"]])
    A <- Q1%*%diag(Lambda)%*%t(Q2)
    #
    Z <- sapply(1:n, function(i){
      rnorm(d, sd = sqrt(Lambda))
    }) %>% t()
    #
    Sigma1 <- diag(rep(1, p1)) - Q1%*%diag(Lambda)%*%t(Q1)
    Sigma2 <- diag(rep(1, p2)) - Q2%*%diag(Lambda)%*%t(Q2)
    X1 <- (Z%*%t(Q1)) + matrix(rnorm(n*p1), nrow=n)%*%chol(Sigma1)
    X2 <- (Z%*%t(Q2)) + matrix(rnorm(n*p2), nrow=n)%*%chol(Sigma2)
    #
    Y1 <- apply(X1, 1, function(x){
      cm_trans(x, U)
    }) %>%
      t()
    #
    Y2 <- apply(X2, 1, function(x){
      cm_trans(x, V)
    }) %>%
      t()
    #
    cop_res <- sbgcop::sbgcop.mcmc(cbind(Y1, Y2), plugin.marginal = rep(T, p1+p2), nsamp = 500)
    cop_samps <- cop_res$C.psamp %>%
      apply(3, function(C){
        C %>%
          cmcca::compute_relation_matrix(c(p1, p2)) %>%
          cmcca::extract_association_matrices(c(p1, p2)) %>%
          magrittr::extract2('(1, 2)') %>%
          c()
      })
    COP_MSE[i, sim] <- sum((rowMeans(cop_samps) - c(A))^2)
    #
    A_CCA <- cbind(Y1, Y2) %>%
      apply(2, function(y){
        y - mean(y)
      }) %>%
      magrittr::divide_by(sqrt(n)) %>%
      t() %>%
      magrittr::multiply_by_matrix(., t(.)) %>%
      cmcca::compute_relation_matrix(c(p1, p2)) %>%
      cmcca::extract_association_matrices(c(p1, p2)) %>%
      magrittr::extract2('(1, 2)') %>%
      c()
    CCA_MSE[i, sim] <- sum((A_CCA - c(A))^2)
    #
    X1 <- rnorm(n*p1) %>% matrix(nrow=n)
    X2 <- rnorm(n*p2) %>% matrix(nrow=n)
    assign1 <- transport::transport(transport::pp(X1),
                                    transport::pp(Y1),
                                    p = 2)
    assign2 <- transport::transport(transport::pp(X2),
                                    transport::pp(Y2),
                                    p = 2)
    X1 <- X1[order(assign1[, 2]), ]
    X2 <- X2[order(assign2[, 2]), ]
    A_CCA <- cbind(X1, X2) %>%
      apply(2, function(y){
        y - mean(y)
      }) %>%
      magrittr::divide_by(sqrt(n)) %>%
      t() %>%
      magrittr::multiply_by_matrix(., t(.)) %>%
      cmcca::compute_relation_matrix(c(p1, p2)) %>%
      cmcca::extract_association_matrices(c(p1, p2)) %>%
      magrittr::extract2('(1, 2)') %>%
      c()
    PLUG_MSE[i, sim] <- sum((A_CCA - c(A))^2)
  }
}
#
df <- data.frame(n = rep(N, 3),
                 Method = rep(c("CCA", "GCCCA", "CMCCA"), each = length(N)),
                 rbind(CCA_MSE, COP_MSE, PLUG_MSE)) %>%
  reshape2::melt(id.vars = c("n", "Method")) %>%
  dplyr::mutate(Method = factor(Method, levels = c('CCA', 'GCCCA', 'CMCCA')))

#
pdf("figures_tables/mse_cm.pdf", family="Times", height = 5, width = 7)
par(mar = c(4, 5, 4, 2))
boxplot(value ~ n + Method, data = df, names = c("50", "100", "250", "500", "1000",
                                                 "50", "100", "250", "500", "1000",
                                                 "50", "100", "250", "500", "1000"),
        at = c(1:5, 7:11, 13:17), xlab = "n", ylab = expression(L[W]),
        col = c("#1b9e77", "#1b9e77", "#1b9e77", "#1b9e77", "#1b9e77",
                "#d95f02", "#d95f02", "#d95f02", "#d95f02", "#d95f02",
                "#7570b3", "#7570b3", "#7570b3", "#7570b3", "#7570b3"), outline = F,
        ylim = c(0, 0.62), cex.lab = 1.2)
legend("topleft", fill = c("#1b9e77", "#d95f02", "#7570b3"), legend = c("CCA","GCCCA","CMCCA"), horiz = T)
dev.off()
