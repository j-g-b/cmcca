#
library(magrittr)
library(tidyverse)

#
DD <- readr::read_csv("data/DD.csv") %>%
      dplyr::filter(!grepl("2020|2021|2022", Date))
OLN <- readr::read_csv("data/OLN.csv") %>%
       dplyr::filter(!grepl("2020|2021|2022", Date))
BHP <- readr::read_csv("data/BHP.csv") %>%
       dplyr::filter(!grepl("2020|2021|2022", Date))
PPG <- readr::read_csv("data/PPG.csv") %>%
       dplyr::filter(!grepl("2020|2021|2022", Date))
LUMN <- readr::read_csv("data/LUMN.csv") %>%
        dplyr::filter(!grepl("2020|2021|2022", Date))
ATT <- readr::read_csv("data/T.csv") %>%
       dplyr::filter(!grepl("2020|2021|2022", Date))
VZ <- readr::read_csv("data/VZ.csv") %>%
      dplyr::filter(!grepl("2020|2021|2022", Date))
CMCSA <- readr::read_csv("data/CMCSA.csv") %>%
         dplyr::filter(!grepl("2020|2021|2022", Date))

#
MATS <- data.frame(DD = log(DD$`Adj Close`[2:nrow(DD)] / DD$`Adj Close`[1:(nrow(DD)-1)]),
                   OLN = log(OLN$`Adj Close`[2:nrow(OLN)] / OLN$`Adj Close`[1:(nrow(OLN)-1)]),
                   BHP = log(BHP$`Adj Close`[2:nrow(BHP)] / BHP$`Adj Close`[1:(nrow(BHP)-1)]),
                   PPG = log(PPG$`Adj Close`[2:nrow(PPG)] / PPG$`Adj Close`[1:(nrow(PPG)-1)]),
                   Date = DD$Date[2:nrow(DD)])
COMMS <- data.frame(LUMN = log(LUMN$`Adj Close`[2:nrow(LUMN)] / LUMN$`Adj Close`[1:(nrow(LUMN)-1)]),
                    ATT = log(ATT$`Adj Close`[2:nrow(ATT)] / ATT$`Adj Close`[1:(nrow(ATT)-1)]),
                    VZ = log(VZ$`Adj Close`[2:nrow(VZ)] / VZ$`Adj Close`[1:(nrow(VZ)-1)]),
                    CMCSA = log(CMCSA$`Adj Close`[2:nrow(CMCSA)] / CMCSA$`Adj Close`[1:(nrow(CMCSA)-1)]),
                    Date = LUMN$Date[2:nrow(LUMN)])

#
MVN::mvn(MATS[, 1:4] %>% apply(2, cmcca::z_scores), mvnTest = "energy")
MVN::mvn(COMMS[ 1:4] %>% apply(2, cmcca::z_scores), mvnTest = "energy")

#
MATS_REAGAN <- dplyr::filter(MATS, as.numeric(substr(Date, 1, 4)) < 1992)[, 1:4] %>% as.matrix()
COMMS_REAGAN <- dplyr::filter(COMMS, as.numeric(substr(Date, 1, 4)) < 1992)[, 1:4] %>% as.matrix()
reagan_cca <- cmcca::cmcca_mcmc(MATS_REAGAN, COMMS_REAGAN,
                                     burn_in = 1000, thin = 5, iter = 10000)

#
MATS_GLOBAL <- dplyr::filter(MATS, as.numeric(substr(Date, 1, 4)) >= 1992) %>%
               dplyr::filter(as.numeric(substr(Date, 1, 4)) < 2008) %>%
               magrittr::extract(1:4) %>%
               as.matrix()
COMMS_GLOBAL <- dplyr::filter(COMMS, as.numeric(substr(Date, 1, 4)) >= 1992) %>%
                dplyr::filter(as.numeric(substr(Date, 1, 4)) < 2008) %>%
                magrittr::extract(1:4) %>%
                as.matrix()
global_cca <- cmcca::cmcca_mcmc(MATS_GLOBAL, COMMS_GLOBAL,
                                     burn_in = 1000, thin = 5, iter = 10000)

#
MATS_RECESS <- dplyr::filter(MATS, as.numeric(substr(Date, 1, 4)) >= 2008) %>%
               magrittr::extract(1:4) %>%
               as.matrix()
COMMS_RECESS <- dplyr::filter(COMMS, as.numeric(substr(Date, 1, 4)) >= 2008) %>%
                magrittr::extract(1:4) %>%
                as.matrix()
recess_cca <- cmcca::cmcca_mcmc(MATS_RECESS, COMMS_RECESS,
                                burn_in = 1000, thin = 5, iter = 10000)

#
pdf("figures_tables/stock_canonical_corrs_posterior.pdf", family = "Times", width = 8, height = 6)
par(mfrow = c(3, 2), mar = c(2.5, 2.5, 2.5, 2.5))
Lambda1_range <- c(0, 1)
Lambda2_range <- c(0, 1)

reagan_cca$Lambda[1, ] %>%
  hist(breaks = 25, xlim = c(Lambda1_range[1], Lambda1_range[2]), freq = F,
       main = expression("1985-1991 " ~ lambda[1]), yaxt= "n", col = rgb(215,48,39, maxColorValue = 255))
reagan_first_hdi <- HDInterval::hdi(reagan_cca$Lambda[1,], 0.95)
abline(v = reagan_first_hdi[1], lt = 2)
abline(v = reagan_first_hdi[2], lt = 2)

reagan_cca$Lambda[2, ] %>%
  hist(breaks = 25, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("1985-1991 " ~ lambda[2]), yaxt = "n", col = rgb(215,48,39, maxColorValue = 255, alpha = 100))
reagan_second_hdi <- HDInterval::hdi(reagan_cca$Lambda[2,], 0.95)
abline(v = reagan_second_hdi[1], lt = 2)
abline(v = reagan_second_hdi[2], lt = 2)

global_cca$Lambda[1, ] %>%
  hist(breaks = 25, xlim = c(Lambda1_range[1], Lambda1_range[2]), freq = F,
       main = expression("1992-2007 " ~ lambda[1]), yaxt = "n", col = rgb(69,117,180, maxColorValue = 255))
global_first_hdi <- HDInterval::hdi(global_cca$Lambda[1,], 0.95)
abline(v = global_first_hdi[1], lt = 2)
abline(v = global_first_hdi[2], lt = 2)

global_cca$Lambda[2, ] %>%
  hist(breaks = 25, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("1992-2007 " ~ lambda[2]), yaxt = "n", col = rgb(69,117,180, maxColorValue = 255, alpha = 100))
global_second_hdi <- HDInterval::hdi(global_cca$Lambda[2,], 0.95)
abline(v = global_second_hdi[1], lt = 2)
abline(v = global_second_hdi[2], lt = 2)

recess_cca$Lambda[1, ] %>%
  hist(breaks = 25, xlim = c(Lambda1_range[1], Lambda1_range[2]), freq = F,
       main = expression("2008-2019 " ~ lambda[1]), yaxt= "n", col = rgb(254,224,144, maxColorValue = 255))
recess_first_hdi <- HDInterval::hdi(recess_cca$Lambda[1,], 0.95)
abline(v = recess_first_hdi[1], lt = 2)
abline(v = recess_first_hdi[2], lt = 2)

recess_cca$Lambda[2, ] %>%
  hist(breaks = 25, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("2008-2019 " ~ lambda[2]), yaxt = "n", col = rgb(254,224,144, maxColorValue = 255, alpha = 100))
recess_second_hdi <- HDInterval::hdi(recess_cca$Lambda[2,], 0.95)
abline(v = recess_second_hdi[1], lt = 2)
abline(v = recess_second_hdi[2], lt = 2)
dev.off()

#
pdf("figures_tables/stock_canonical_corrs_2d_posterior.pdf", family = "Times", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(4, 5, 4, 2), mgp = c(3, 1, 0))
Lambda1_range <- c(0, 1)
Lambda2_range <- c(0, 1)

plot(reagan_cca$Lambda[1, ], reagan_cca$Lambda[2, ], col = rgb(215,48,39, maxColorValue = 255, alpha = 100),
     xlim =  c(Lambda1_range[1], Lambda1_range[2]), ylim = c(Lambda2_range[1], Lambda2_range[2]),
     main = "1985-1991", xlab = expression(lambda[1]), ylab = expression(lambda[2]), cex.axis = 1.5, cex.lab = 2, cex.main = 2)
abline(0, 1)
z_reag <- MASS::kde2d(reagan_cca$Lambda[1, ],
                      reagan_cca$Lambda[2, ],
                      n=50)
contour(z_reag, drawlabels=FALSE, nlevels=10, col=rgb(215,48,39, maxColorValue = 255), add=TRUE)

plot(global_cca$Lambda[1, ], global_cca$Lambda[2, ], col = rgb(69,117,180, maxColorValue = 255, alpha = 100),
     xlim =  c(Lambda1_range[1], Lambda1_range[2]), ylim = c(Lambda2_range[1], Lambda2_range[2]),
     main = "1992-2007", xlab = expression(lambda[1]), ylab = expression(lambda[2]), cex.axis = 1.5, cex.lab = 2, cex.main = 2)
abline(0, 1)
z_global <- MASS::kde2d(global_cca$Lambda[1, ],
                      global_cca$Lambda[2, ],
                      n=50)
contour(z_global, drawlabels=FALSE, nlevels=10, col=rgb(69,117,180, maxColorValue = 255), add=TRUE)

plot(recess_cca$Lambda[1, ], recess_cca$Lambda[2, ], col = rgb(254,224,144, maxColorValue = 255, alpha = 100),
     xlim =  c(Lambda1_range[1], Lambda1_range[2]), ylim = c(Lambda2_range[1], Lambda2_range[2]),
     main = "2008-2019", xlab = expression(lambda[1]), ylab = expression(lambda[2]), cex.axis = 1.5, cex.lab = 2, cex.main = 2)
abline(0, 1)
z_recess <- MASS::kde2d(recess_cca$Lambda[1, ],
                        recess_cca$Lambda[2, ],
                        n=50)
contour(z_recess, drawlabels=FALSE, nlevels=10, col=rgb(254,224,144, maxColorValue = 255), add=TRUE)
dev.off()

#
reagan_first_cc <- sapply(1:ncol(reagan_cca$Lambda), function(s){
  reagan_cca$Q1[, 1, s]%*%t(reagan_cca$Q2[, 1, s])
}) %>%
  t()
global_first_cc <- sapply(1:ncol(global_cca$Lambda), function(s){
  global_cca$Q1[, 1, s]%*%t(global_cca$Q2[, 1, s])
}) %>%
  t()
recess_first_cc <- sapply(1:ncol(recess_cca$Lambda), function(s){
  recess_cca$Q1[, 1, s]%*%t(recess_cca$Q2[, 1, s])
}) %>%
  t()

reagan_first_cc_cr <- reagan_first_cc %>%
  apply(2, function(x){
    hdi <- HDInterval::hdi(x, 0.95)
    c(hdi[1], mean = mean(x), hdi[2])
  }) %>%
  t()
global_first_cc_cr <- global_first_cc %>%
  apply(2, function(x){
    hdi <- HDInterval::hdi(x, 0.95)
    c(hdi[1], mean = mean(x), hdi[2])
  }) %>%
  t()
recess_first_cc_cr <- recess_first_cc %>%
  apply(2, function(x){
    hdi <- HDInterval::hdi(x, 0.95)
    c(hdi[1], mean = mean(x), hdi[2])
  }) %>%
  t()

pdf("figures_tables/stock_group_similarities_posterior.pdf", family = "Times", width = 10, height = 6)
par(mfrow = c(1, 2), mar = c(6, 4, 4, 2))
first_cc_axes <- svd(rbind(reagan_first_cc, global_first_cc, recess_first_cc))
z_reag <- MASS::kde2d(first_cc_axes[["u"]][1:nrow(reagan_first_cc), 1]*first_cc_axes[["d"]][1],
                      first_cc_axes[["u"]][1:nrow(reagan_first_cc), 2]*first_cc_axes[["d"]][2],
                      n=50)
z_glob <- MASS::kde2d(first_cc_axes[["u"]][(nrow(reagan_first_cc)+1):(2*nrow(reagan_first_cc)), 1]*first_cc_axes[["d"]][1],
                      first_cc_axes[["u"]][(nrow(reagan_first_cc)+1):(2*nrow(reagan_first_cc)), 2]*first_cc_axes[["d"]][2],
                      n=50)
z_rec <- MASS::kde2d(first_cc_axes[["u"]][(2*nrow(reagan_first_cc)+1):(3*nrow(reagan_first_cc)), 1]*first_cc_axes[["d"]][1],
                      first_cc_axes[["u"]][(2*nrow(reagan_first_cc)+1):(3*nrow(reagan_first_cc)), 2]*first_cc_axes[["d"]][2],
                      n=50)
plot(first_cc_axes[["u"]][, 1]*first_cc_axes[["d"]][1], first_cc_axes[["u"]][, 2]*first_cc_axes[["d"]][2],
     col = c(rep(rgb(215,48,39, maxColorValue = 255, alpha = 50), nrow(reagan_first_cc)),
             rep(rgb(69,117,180, maxColorValue = 255, alpha = 50), nrow(global_first_cc)),
             rep(rgb(254,224,144, maxColorValue = 255, alpha = 50), nrow(recess_first_cc))),
     xlab = "PC1", ylab = "PC2")
contour(z_reag, drawlabels=FALSE, nlevels=10, col=rgb(215,48,39, maxColorValue = 255), add=TRUE)
contour(z_glob, drawlabels=FALSE, nlevels=10, col=rgb(69,117,180, maxColorValue = 255), add=TRUE)
contour(z_rec, drawlabels=FALSE, nlevels=10, col=rgb(254,224,144, maxColorValue = 255), add=TRUE)
legend(-0.3, -0.4, legend=c("1985-1991", "1992-2007", "2008-2019"),
       col=c(rgb(215,48,39, maxColorValue = 255),
             rgb(69,117,180, maxColorValue = 255),
             rgb(254,224,144, maxColorValue = 255)),
       lty=c(1,1,1), cex=0.8)

reag_glob_dot <- (reagan_first_cc*global_first_cc[sample(1:nrow(global_first_cc), nrow(global_first_cc)), ]) %>%
                 rowSums()
glob_rec_dot <- (global_first_cc*recess_first_cc[sample(1:nrow(recess_first_cc), nrow(recess_first_cc)), ]) %>%
                rowSums()
rec_reag_dot <- (recess_first_cc*reagan_first_cc[sample(1:nrow(reagan_first_cc), nrow(reagan_first_cc)), ]) %>%
                rowSums()
boxplot(cbind(reag_glob_dot, glob_rec_dot, rec_reag_dot) %>%
          magrittr::set_colnames(c("1985-1991\n vs\n 1992-2007", "1992-2007\n vs\n 2008-2019", "2008-2019\n vs\n 1985-1991")),
        las = 2, col = "white")
dev.off()

#
reagan_exemplar <- sapply(1:ncol(reagan_cca$Lambda), function(s){
  mean((t(reagan_cca$Q1[, 1, ])%*%reagan_cca$Q1[, 1, s])*(t(reagan_cca$Q2[, 1, ])%*%reagan_cca$Q2[, 1, s]))
})
reagan_exemplar <- cbind(reagan_cca$Q1[, 1, which(reagan_exemplar == max(reagan_exemplar))],
                         reagan_cca$Q2[, 1, which(reagan_exemplar == max(reagan_exemplar))])
reagan_first_sp_mat <- sapply(1:ncol(reagan_cca$Lambda), function(s){
  dot_pos <- sum(reagan_exemplar*cbind(reagan_cca$Q1[, 1, s], reagan_cca$Q2[, 1, s]))
  dot_neg <- sum(reagan_exemplar*cbind(-reagan_cca$Q1[, 1, s], -reagan_cca$Q2[, 1, s]))
  if(dot_pos > dot_neg){
    sign <- 1
  } else {
    sign <- -1
  }
  cor(MATS_REAGAN, sign*reagan_cca$X1_samps[,,s]%*%reagan_cca$Q1[, 1, s], method = "spearman")
}) %>%
  t()
reagan_first_sp_comm <- sapply(1:ncol(reagan_cca$Lambda), function(s){
  dot_pos <- sum(reagan_exemplar*cbind(reagan_cca$Q1[, 1, s], reagan_cca$Q2[, 1, s]))
  dot_neg <- sum(reagan_exemplar*cbind(-reagan_cca$Q1[, 1, s], -reagan_cca$Q2[, 1, s]))
  if(dot_pos > dot_neg){
    sign <- 1
  } else {
    sign <- -1
  }
  cor(COMMS_REAGAN, sign*reagan_cca$X2_samps[,,s]%*%reagan_cca$Q2[, 1, s], method = "spearman")
}) %>%
  t()

global_exemplar <- sapply(1:ncol(global_cca$Lambda), function(s){
  mean((t(global_cca$Q1[, 1, ])%*%global_cca$Q1[, 1, s])*(t(global_cca$Q2[, 1, ])%*%global_cca$Q2[, 1, s]))
})
global_exemplar <- cbind(global_cca$Q1[, 1, which(global_exemplar == max(global_exemplar))],
                         global_cca$Q2[, 1, which(global_exemplar == max(global_exemplar))])
global_first_sp_mat <- sapply(1:ncol(global_cca$Lambda), function(s){
  dot_pos <- sum(global_exemplar*cbind(global_cca$Q1[, 1, s], global_cca$Q2[, 1, s]))
  dot_neg <- sum(global_exemplar*cbind(-global_cca$Q1[, 1, s], -global_cca$Q2[, 1, s]))
  if(dot_pos > dot_neg){
    sign <- 1
  } else {
    sign <- -1
  }
  cor(MATS_GLOBAL, sign*global_cca$X1_samps[,,s]%*%global_cca$Q1[, 1, s], method = "spearman")
}) %>%
  t()
global_first_sp_comm <- sapply(1:ncol(global_cca$Lambda), function(s){
  dot_pos <- sum(global_exemplar*cbind(global_cca$Q1[, 1, s], global_cca$Q2[, 1, s]))
  dot_neg <- sum(global_exemplar*cbind(-global_cca$Q1[, 1, s], -global_cca$Q2[, 1, s]))
  if(dot_pos > dot_neg){
    sign <- 1
  } else {
    sign <- -1
  }
  cor(COMMS_GLOBAL, sign*global_cca$X2_samps[,,s]%*%global_cca$Q2[, 1, s], method = "spearman")
}) %>%
  t()

recess_exemplar <- sapply(1:ncol(recess_cca$Lambda), function(s){
  mean((t(recess_cca$Q1[, 1, ])%*%recess_cca$Q1[, 1, s])*(t(recess_cca$Q2[, 1, ])%*%recess_cca$Q2[, 1, s]))
})
recess_exemplar <- cbind(recess_cca$Q1[, 1, which(recess_exemplar == max(recess_exemplar))],
                         recess_cca$Q2[, 1, which(recess_exemplar == max(recess_exemplar))])
recess_first_sp_mat <- sapply(1:ncol(recess_cca$Lambda), function(s){
  dot_pos <- sum(recess_exemplar*cbind(recess_cca$Q1[, 1, s], recess_cca$Q2[, 1, s]))
  dot_neg <- sum(recess_exemplar*cbind(-recess_cca$Q1[, 1, s], -recess_cca$Q2[, 1, s]))
  if(dot_pos > dot_neg){
    sign <- 1
  } else {
    sign <- -1
  }
  cor(MATS_RECESS, sign*recess_cca$X1_samps[,,s]%*%recess_cca$Q1[, 1, s], method = "spearman")
}) %>%
  t()
recess_first_sp_comm <- sapply(1:ncol(recess_cca$Lambda), function(s){
  dot_pos <- sum(recess_exemplar*cbind(recess_cca$Q1[, 1, s], recess_cca$Q2[, 1, s]))
  dot_neg <- sum(recess_exemplar*cbind(-recess_cca$Q1[, 1, s], -recess_cca$Q2[, 1, s]))
  if(dot_pos > dot_neg){
    sign <- 1
  } else {
    sign <- -1
  }
  cor(COMMS_RECESS, sign*recess_cca$X2_samps[,,s]%*%recess_cca$Q2[, 1, s], method = "spearman")
}) %>%
  t()

mat_sp_cors <- rbind(paste0(round(apply(reagan_first_sp_mat, 2, mean), 3), ", [",
                            round(apply(reagan_first_sp_mat, 2, function(x){HDInterval::hdi(x, 0.95)[1]}), 3), ", ",
                            round(apply(reagan_first_sp_mat, 2, function(x){HDInterval::hdi(x, 0.95)[2]}), 3), "]"),
                     paste0(round(apply(global_first_sp_mat, 2, mean), 3), ", [",
                            round(apply(global_first_sp_mat, 2, function(x){HDInterval::hdi(x, 0.95)[1]}), 3), ", ",
                            round(apply(global_first_sp_mat, 2, function(x){HDInterval::hdi(x, 0.95)[2]}), 3), "]"),
                     paste0(round(apply(recess_first_sp_mat, 2, mean), 3), ", [",
                            round(apply(recess_first_sp_mat, 2, function(x){HDInterval::hdi(x, 0.95)[1]}), 3), ", ",
                            round(apply(recess_first_sp_mat, 2, function(x){HDInterval::hdi(x, 0.95)[2]}), 3), "]")) %>%
              magrittr::set_rownames(c("1985-1991", "1992-2007", "2008-2019")) %>%
              magrittr::set_colnames(c("DD", "OLN", "BHP", "PPG"))
xtable::xtable(mat_sp_cors) %>%
  xtable::autoformat() %>%
  print(file = "figures_tables/mat_sp_cors.txt")

comm_sp_cors <- rbind(paste0(round(apply(reagan_first_sp_comm, 2, mean), 3), ", [",
                             round(apply(reagan_first_sp_comm, 2, function(x){HDInterval::hdi(x, 0.95)[1]}), 3), ", ",
                             round(apply(reagan_first_sp_comm, 2, function(x){HDInterval::hdi(x, 0.95)[2]}), 3), "]"),
                     paste0(round(apply(global_first_sp_comm, 2, mean), 3), ", [",
                            round(apply(global_first_sp_comm, 2, function(x){HDInterval::hdi(x, 0.95)[1]}), 3), ", ",
                            round(apply(global_first_sp_comm, 2, function(x){HDInterval::hdi(x, 0.95)[2]}), 3), "]"),
                     paste0(round(apply(recess_first_sp_comm, 2, mean), 3), ", [",
                            round(apply(recess_first_sp_comm, 2, function(x){HDInterval::hdi(x, 0.95)[1]}), 3), ", ",
                            round(apply(recess_first_sp_comm, 2, function(x){HDInterval::hdi(x, 0.95)[2]}), 3), "]")) %>%
  magrittr::set_rownames(c("1985-1991", "1992-2007", "2008-2019")) %>%
  magrittr::set_colnames(c("LUMN", "ATT", "VZ", "CMCSA"))
xtable::xtable(comm_sp_cors) %>%
  xtable::autoformat() %>%
  print(file = "figures_tables/comm_sp_cors.txt")


