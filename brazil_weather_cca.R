#
#
library(tidyverse)
library(magrittr)

#
#
weather_data <- readr::read_delim("data/conventional_weather_stations_inmet_brazil_1961_2019.csv", delim = ";")
station_data <- readr::read_delim("data/weather_stations_codes.csv", delim = ";")

#
#
station_clusters <- ifelse(station_data$Latitude > -10 & station_data$Longitude < -46, 1,
                           ifelse(station_data$Longitude >= -48 & station_data$Latitude > -15, 2,
                                  ifelse(station_data$Longitude < -48 & station_data$Latitude > -24, 3,
                                         ifelse(station_data$Longitude >= -48 & station_data$Latitude <= -15, 4, 5))))
station_data[["cluster"]] <- station_clusters
station_data %<>% magrittr::set_colnames(c("Name", "Estacao", "Latitude", "Longitude",
                                           "Altitude", "Op Status", "Op Start", "cluster"))

#
#
weather_data %<>% dplyr::mutate(Day = substr(Data, 1, 2),
                                Month = substr(Data, 4, 5),
                                Year = substr(Data, 7, 10)) %>%
                  dplyr::left_join(dplyr::select(station_data, cluster, Estacao)) %>%
                  dplyr::group_by(Year, Month, cluster) %>%
                  dplyr::summarise(mean_temp = mean(`Temp Comp Media`, na.rm = T),
                                   mean_atm = mean(PressaoAtmEstacao, na.rm = T),
                                   mean_evap = mean(`Evaporacao Piche`, na.rm = T),
                                   mean_insol = mean(Insolacao, na.rm = T),
                                   mean_wind_vel = mean(VelocidadeVento, na.rm = T),
                                   mean_rel_humid = mean(`Umidade Relativa Media`, na.rm = T),
                                   mean_cloud = mean(Nebulosidade, na.rm = T),
                                   mean_precip = mean(Precipitacao, na.rm = T)) %>%
                  dplyr::ungroup()

#
#
detrended_weather_data <- weather_data %>%
                          dplyr::group_by(Year, cluster) %>%
                          dplyr::mutate(mean_temp = mean_temp - mean(mean_temp, na.rm = T),
                                        mean_atm = mean_atm - mean(mean_atm, na.rm = T),
                                        mean_evap = mean_evap - mean(mean_evap, na.rm = T),
                                        mean_insol = mean_insol - mean(mean_insol, na.rm = T),
                                        mean_wind_vel = mean_wind_vel - mean(mean_wind_vel, na.rm = T),
                                        mean_rel_humid = mean_rel_humid - mean(mean_rel_humid, na.rm = T),
                                        mean_cloud = mean_cloud - mean(mean_cloud, na.rm = T),
                                        mean_precip = mean_precip - mean(mean_precip, na.rm = T)) %>%
                          dplyr::ungroup()
#
cl_list <- list()
for(k in 1:5){
  #
  w <- dplyr::filter(detrended_weather_data, cluster == k)

  #
  dat_list <- list()
  y <- unique(w$Year)
  for(j in 2:length(y)){
    prev <- w[w$Year == y[j-1], ]
    curr <- w[w$Year == y[j], ]
    dat_mat <- log(as.matrix(curr[, 4:ncol(w)])/as.matrix(prev[ , 4:ncol(w)]))
    dat_list[[j-1]] <- dat_mat
  }
  cl_list[[k]] <- do.call(rbind, dat_list) %>%
                  cbind(rep(k, nrow(.)), .)
}
detrended_weather_data <- do.call(rbind, cl_list)
detrended_weather_data <- detrended_weather_data[apply(detrended_weather_data, 1, function(x){!any(is.na(x))}), ]

# Test Gaussian copula model
#
zscore_mvn_test <- matrix(NA, nrow = 5, ncol = 2)
for(k in 1:5){
  #
  cl_dat <- detrended_weather_data[detrended_weather_data[, 1] == k, 2:ncol(detrended_weather_data)]

  #
  zscore_mvn_test[k, ] <- MVN::mvn(cl_dat %>% apply(2, cmcca::z_scores))$multivariateNormality[2:3] %>% as.matrix() %>% c()
}

#
cl_dat <- detrended_weather_data[detrended_weather_data[, 1] == 1, 2:ncol(detrended_weather_data)]
pdf("figures_tables/brazil_heat_pairs.pdf", family = "Times", width = 6, height = 6)
cmcca:::cmcca_pairs(cl_dat[, 1:5] %>% magrittr::set_colnames(c("Temperature", "Atm. Pressure", "Evaporation",
                                                               "Insolation", "Wind Vel.")))
dev.off()
pdf("figures_tables/brazil_heat_pairs_zscore.pdf", family = "Times", width = 6, height = 6)
cmcca:::cmcca_pairs(cl_dat[, 1:5] %>% apply(2, cmcca::z_scores) %>% magrittr::set_colnames(c("Temperature", "Atm. Pressure", "Evaporation",
                                                                                             "Insolation", "Wind Vel.")))
dev.off()
pdf("figures_tables/brazil_h20_pairs.pdf", family = "Times", width = 6, height = 6)
cmcca:::cmcca_pairs(cl_dat[, 6:8] %>% magrittr::set_colnames(c("Rel. Humidity", "Cloudiness", "Precip.")))
dev.off()
pdf("figures_tables/brazil_h20_pairs_zscore.pdf", family = "Times", width = 6, height = 6)
cmcca:::cmcca_pairs(cl_dat[, 6:8] %>% apply(2, cmcca::z_scores) %>% magrittr::set_colnames(c("Rel. Humidity", "Cloudiness", "Precip.")))
dev.off()

# Test our model
#
mzscore_mvn_test <- matrix(NA, nrow = 5, ncol = 2)
for(k in 1:5){
  #
  cl_dat <- detrended_weather_data[detrended_weather_data[, 1] == k, 2:ncol(detrended_weather_data)]

  #
  Z1 <- matrix(rnorm(nrow(cl_dat)*5), nrow = nrow(cl_dat))
  Z2 <- matrix(rnorm(nrow(cl_dat)*3), nrow = nrow(cl_dat))

  #
  assign1 <- transport::transport(transport::pp(cl_dat[, 1:5]),
                                  transport::pp(Z1),
                                  p = 2)
  assign2 <- transport::transport(transport::pp(cl_dat[, 6:8]),
                                  transport::pp(Z2),
                                  p = 2)
  Z1 <- Z1[order(assign1[, 2]), ]
  Z2 <- Z2[order(assign2[, 2]), ]

  #
  mzscore_mvn_test[k, ] <- MVN::mvn(cbind(Z1, Z2))$multivariateNormality[2:3] %>% as.matrix() %>% c()
}
cbind(zscore_mvn_test, mzscore_mvn_test)  %>%
  t() %>%
  magrittr::set_colnames(paste0("Region ", 1:5)) %>%
  magrittr::set_rownames(c("H-Z statistic (normal scores)", "p-value (normal scores)",
                           "H-Z statistic (m.v. normal scores)", "p-value (m.v. normal scores)")) %>%
  xtable::xtable() %>%
  xtable::autoformat() %>%
  print(file = "figures_tables/brazil_mzscore_mvn_test.txt")

#
clust1_cca <- cmcca::cmcca_mcmc(detrended_weather_data[detrended_weather_data[, 1] == 1, 2:6],
                                detrended_weather_data[detrended_weather_data[, 1] == 1, 7:9],
                                iter = 10000, thin = 5, burn_in = 1000)
clust2_cca <- cmcca::cmcca_mcmc(detrended_weather_data[detrended_weather_data[, 1] == 2, 2:6],
                                detrended_weather_data[detrended_weather_data[, 1] == 2, 7:9],
                                iter = 10000, thin = 5, burn_in = 1000)
clust3_cca <- cmcca::cmcca_mcmc(detrended_weather_data[detrended_weather_data[, 1] == 3, 2:6],
                                detrended_weather_data[detrended_weather_data[, 1] == 3, 7:9],
                                iter = 10000, thin = 5, burn_in = 1000)
clust4_cca <- cmcca::cmcca_mcmc(detrended_weather_data[detrended_weather_data[, 1] == 4, 2:6],
                                detrended_weather_data[detrended_weather_data[, 1] == 4, 7:9],
                                iter = 10000, thin = 5, burn_in = 1000)
clust5_cca <- cmcca::cmcca_mcmc(detrended_weather_data[detrended_weather_data[, 1] == 5, 2:6],
                                detrended_weather_data[detrended_weather_data[, 1] == 5, 7:9],
                                iter = 10000, thin = 5, burn_in = 1000)

#
#
pdf("figures_tables/brazil_canonical_corrs_posterior.pdf", family = "Times", width = 8, height = 6)
par(mfrow = c(5, 3), mar = c(2.5, 2, 2.5, 2))
Lambda1_range <- c(0, 1)
Lambda2_range <- c(0, 1)

clust1_cca$Lambda[1, ] %>%
  hist(breaks = 20, xlim = c(Lambda1_range[1], Lambda1_range[2]), freq = F,
       main = expression("Region 1 " ~ lambda[1]), yaxt= "n", col = rgb(215,48,39, maxColorValue = 255))
clust1_first_hdi <- HDInterval::hdi(clust1_cca$Lambda[1,], 0.95)
abline(v = clust1_first_hdi[1], lt = 2)
abline(v = clust1_first_hdi[2], lt = 2)

clust1_cca$Lambda[2, ] %>%
  hist(breaks = 20, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("Region 1 " ~ lambda[2]), yaxt = "n", col = rgb(215,48,39, maxColorValue = 255, alpha = 180))
clust1_second_hdi <- HDInterval::hdi(clust1_cca$Lambda[2,], 0.95)
abline(v = clust1_second_hdi[1], lt = 2)
abline(v = clust1_second_hdi[2], lt = 2)

clust1_cca$Lambda[3, ] %>%
  hist(breaks = 20, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("Region 1 " ~ lambda[3]), yaxt = "n", col = rgb(215,48,39, maxColorValue = 255, alpha = 100))
clust1_third_hdi <- HDInterval::hdi(clust1_cca$Lambda[3,], 0.95)
abline(v = clust1_third_hdi[1], lt = 2)
abline(v = clust1_third_hdi[2], lt = 2)

clust2_cca$Lambda[1, ] %>%
  hist(breaks = 10, xlim = c(Lambda1_range[1], Lambda1_range[2]), freq = F,
       main = expression("Region 2 " ~ lambda[1]), yaxt= "n", col = rgb(153,112,171, maxColorValue = 255))
clust2_first_hdi <- HDInterval::hdi(clust2_cca$Lambda[1,], 0.95)
abline(v = clust2_first_hdi[1], lt = 2)
abline(v = clust2_first_hdi[2], lt = 2)

clust2_cca$Lambda[2, ] %>%
  hist(breaks = 20, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("Region 2 " ~ lambda[2]), yaxt = "n", col = rgb(153,112,171, maxColorValue = 255, alpha = 180))
clust2_second_hdi <- HDInterval::hdi(clust2_cca$Lambda[2,], 0.95)
abline(v = clust2_second_hdi[1], lt = 2)
abline(v = clust2_second_hdi[2], lt = 2)

clust2_cca$Lambda[3, ] %>%
  hist(breaks = 20, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("Region 2 " ~ lambda[3]), yaxt = "n", col = rgb(153,112,171, maxColorValue = 255, alpha = 100))
clust2_third_hdi <- HDInterval::hdi(clust2_cca$Lambda[3,], 0.95)
abline(v = clust2_third_hdi[1], lt = 2)
abline(v = clust2_third_hdi[2], lt = 2)

clust3_cca$Lambda[1, ] %>%
  hist(breaks = 20, xlim = c(Lambda1_range[1], Lambda1_range[2]), freq = F,
       main = expression("Region 3 " ~ lambda[1]), yaxt= "n", col = rgb(90,174,97, maxColorValue = 255))
clust3_first_hdi <- HDInterval::hdi(clust3_cca$Lambda[1,], 0.95)
abline(v = clust3_first_hdi[1], lt = 2)
abline(v = clust3_first_hdi[2], lt = 2)

clust3_cca$Lambda[2, ] %>%
  hist(breaks = 20, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("Region 3 " ~ lambda[2]), yaxt = "n", col = rgb(90,174,97, maxColorValue = 255, alpha = 180))
clust3_second_hdi <- HDInterval::hdi(clust3_cca$Lambda[2,], 0.95)
abline(v = clust3_second_hdi[1], lt = 2)
abline(v = clust3_second_hdi[2], lt = 2)

clust3_cca$Lambda[3, ] %>%
  hist(breaks = 20, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("Region 3 " ~ lambda[3]), yaxt = "n", col = rgb(90,174,97, maxColorValue = 255, alpha = 100))
clust3_third_hdi <- HDInterval::hdi(clust3_cca$Lambda[3,], 0.95)
abline(v = clust3_third_hdi[1], lt = 2)
abline(v = clust3_third_hdi[2], lt = 2)

clust4_cca$Lambda[1, ] %>%
  hist(breaks = 20, xlim = c(Lambda1_range[1], Lambda1_range[2]), freq = F,
       main = expression("Region 4 " ~ lambda[1]), yaxt= "n", col = rgb(254,224,144, maxColorValue = 255))
clust4_first_hdi <- HDInterval::hdi(clust4_cca$Lambda[1,], 0.95)
abline(v = clust4_first_hdi[1], lt = 2)
abline(v = clust4_first_hdi[2], lt = 2)

clust4_cca$Lambda[2, ] %>%
  hist(breaks = 20, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("Region 4 " ~ lambda[2]), yaxt = "n", col = rgb(254,224,144, maxColorValue = 255, alpha = 180))
clust4_second_hdi <- HDInterval::hdi(clust4_cca$Lambda[2,], 0.95)
abline(v = clust4_second_hdi[1], lt = 2)
abline(v = clust4_second_hdi[2], lt = 2)

clust4_cca$Lambda[3, ] %>%
  hist(breaks = 20, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("Region 4 " ~ lambda[3]), yaxt = "n", col = rgb(254,224,144, maxColorValue = 255, alpha = 100))
clust4_third_hdi <- HDInterval::hdi(clust4_cca$Lambda[3,], 0.95)
abline(v = clust4_third_hdi[1], lt = 2)
abline(v = clust4_third_hdi[2], lt = 2)

clust5_cca$Lambda[1, ] %>%
  hist(breaks = 20, xlim = c(Lambda1_range[1], Lambda1_range[2]), freq = F,
       main = expression("Region 5 " ~ lambda[1]), yaxt= "n", col = rgb(69,117,180, maxColorValue = 255))
clust5_first_hdi <- HDInterval::hdi(clust5_cca$Lambda[1,], 0.95)
abline(v = clust5_first_hdi[1], lt = 2)
abline(v = clust5_first_hdi[2], lt = 2)

clust5_cca$Lambda[2, ] %>%
  hist(breaks = 20, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("Region 5 " ~ lambda[2]), yaxt = "n", col = rgb(69,117,180, maxColorValue = 255, alpha = 180))
clust5_second_hdi <- HDInterval::hdi(clust5_cca$Lambda[2,], 0.95)
abline(v = clust5_second_hdi[1], lt = 2)
abline(v = clust5_second_hdi[2], lt = 2)

clust5_cca$Lambda[3, ] %>%
  hist(breaks = 20, xlim = c(Lambda2_range[1], Lambda2_range[2]), freq = F,
       main = expression("Region 5 " ~ lambda[3]), yaxt = "n", col = rgb(69,117,180, maxColorValue = 255, alpha = 100))
clust5_third_hdi <- HDInterval::hdi(clust5_cca$Lambda[3,], 0.95)
abline(v = clust5_third_hdi[1], lt = 2)
abline(v = clust5_third_hdi[2], lt = 2)
dev.off()

#
matrix(c(paste0(round(mean(clust1_cca$Lambda[1, ]), 3), " [",
              paste(round(HDInterval::hdi(clust1_cca$Lambda[1,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust1_cca$Lambda[2, ]), 3), " [",
                paste(round(HDInterval::hdi(clust1_cca$Lambda[2,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust1_cca$Lambda[3, ]), 3), " [",
                paste(round(HDInterval::hdi(clust1_cca$Lambda[3,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust2_cca$Lambda[1, ]), 3), " [",
                paste(round(HDInterval::hdi(clust2_cca$Lambda[1,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust2_cca$Lambda[2, ]), 3), " [",
                paste(round(HDInterval::hdi(clust2_cca$Lambda[2,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust2_cca$Lambda[3, ]), 3), " [",
                paste(round(HDInterval::hdi(clust2_cca$Lambda[3,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust3_cca$Lambda[1, ]), 3), " [",
                paste(round(HDInterval::hdi(clust3_cca$Lambda[1,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust3_cca$Lambda[2, ]), 3), " [",
                paste(round(HDInterval::hdi(clust3_cca$Lambda[2,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust3_cca$Lambda[3, ]), 3), " [",
                paste(round(HDInterval::hdi(clust3_cca$Lambda[3,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust4_cca$Lambda[1, ]), 3), " [",
                paste(round(HDInterval::hdi(clust4_cca$Lambda[1,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust4_cca$Lambda[2, ]), 3), " [",
                paste(round(HDInterval::hdi(clust4_cca$Lambda[2,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust4_cca$Lambda[3, ]), 3), " [",
                paste(round(HDInterval::hdi(clust4_cca$Lambda[3,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust5_cca$Lambda[1, ]), 3), " [",
                paste(round(HDInterval::hdi(clust5_cca$Lambda[1,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust5_cca$Lambda[2, ]), 3), " [",
                paste(round(HDInterval::hdi(clust5_cca$Lambda[2,], 0.95), 3), collapse=", "), "]"),
         paste0(round(mean(clust5_cca$Lambda[3, ]), 3), " [",
                paste(round(HDInterval::hdi(clust5_cca$Lambda[3,], 0.95), 3), collapse=", "), "]")),
       nrow = 3) %>%
  magrittr::set_rownames(c("Lambda 1", "Lambda 2", "Lambda 3")) %>%
  magrittr::set_colnames(paste0("Region ", 1:5)) %>%
  xtable::xtable() %>%
  xtable::autoformat() %>%
  print(file = "figures_tables/brazil_canonical_corrs_posterior.txt")

#
clust1_first_cc <- sapply(1:ncol(clust1_cca$Lambda), function(s){
  clust1_cca$Q1[, 1, s]%*%t(clust1_cca$Q2[, 1, s])
}) %>%
  t()
clust2_first_cc <- sapply(1:ncol(clust2_cca$Lambda), function(s){
  clust2_cca$Q1[, 1, s]%*%t(clust2_cca$Q2[, 1, s])
}) %>%
  t()
clust3_first_cc <- sapply(1:ncol(clust3_cca$Lambda), function(s){
  clust3_cca$Q1[, 1, s]%*%t(clust3_cca$Q2[, 1, s])
}) %>%
  t()
clust4_first_cc <- sapply(1:ncol(clust4_cca$Lambda), function(s){
  clust4_cca$Q1[, 1, s]%*%t(clust4_cca$Q2[, 1, s])
}) %>%
  t()
clust5_first_cc <- sapply(1:ncol(clust5_cca$Lambda), function(s){
  clust5_cca$Q1[, 1, s]%*%t(clust5_cca$Q2[, 1, s])
}) %>%
  t()

#
#
pdf("figures_tables/brazil_group_similarities_posterior.pdf", family = "Times", width = 6, height = 6)
par(mfrow = c(1, 1))
first_cc_axes <- svd(rbind(clust1_first_cc, clust2_first_cc, clust3_first_cc, clust4_first_cc, clust5_first_cc))
z_clust1 <- MASS::kde2d(-first_cc_axes[["u"]][1:nrow(clust1_first_cc), 1]*first_cc_axes[["d"]][1],
                        first_cc_axes[["u"]][1:nrow(clust1_first_cc), 2]*first_cc_axes[["d"]][2],
                        n=50)
z_clust2 <- MASS::kde2d(-first_cc_axes[["u"]][(nrow(clust1_first_cc)+1):(2*nrow(clust2_first_cc)), 1]*first_cc_axes[["d"]][1],
                        first_cc_axes[["u"]][(nrow(clust1_first_cc)+1):(2*nrow(clust2_first_cc)), 2]*first_cc_axes[["d"]][2],
                        n=50)
z_clust3 <- MASS::kde2d(-first_cc_axes[["u"]][(2*nrow(clust1_first_cc)+1):(3*nrow(clust2_first_cc)), 1]*first_cc_axes[["d"]][1],
                        first_cc_axes[["u"]][(2*nrow(clust1_first_cc)+1):(3*nrow(clust2_first_cc)), 2]*first_cc_axes[["d"]][2],
                        n=50)
z_clust4 <- MASS::kde2d(-first_cc_axes[["u"]][(3*nrow(clust1_first_cc)+1):(4*nrow(clust2_first_cc)), 1]*first_cc_axes[["d"]][1],
                        first_cc_axes[["u"]][(3*nrow(clust1_first_cc)+1):(4*nrow(clust2_first_cc)), 2]*first_cc_axes[["d"]][2],
                        n=50)
z_clust5 <- MASS::kde2d(-first_cc_axes[["u"]][(4*nrow(clust1_first_cc)+1):(5*nrow(clust2_first_cc)), 1]*first_cc_axes[["d"]][1],
                        first_cc_axes[["u"]][(4*nrow(clust1_first_cc)+1):(5*nrow(clust2_first_cc)), 2]*first_cc_axes[["d"]][2],
                        n=50)
plot(-first_cc_axes[["u"]][, 1]*first_cc_axes[["d"]][1], first_cc_axes[["u"]][, 2]*first_cc_axes[["d"]][2],
     col = c(rep(rgb(215,48,39, maxColorValue = 255, alpha = 50), nrow(clust1_first_cc)),
             rep(rgb(153,112,171, maxColorValue = 255, alpha = 50), nrow(clust2_first_cc)),
             rep(rgb(90,174,97, maxColorValue = 255, alpha = 50), nrow(clust3_first_cc)),
             rep(rgb(254,224,144, maxColorValue = 255, alpha = 50), nrow(clust4_first_cc)),
             rep(rgb(69,117,180, maxColorValue = 255, alpha = 50), nrow(clust5_first_cc))),
     xlab = "PC1", ylab = "PC2")
contour(z_clust1, drawlabels=FALSE, nlevels=10, col=rgb(215,48,39, maxColorValue = 255), add=TRUE)
contour(z_clust2, drawlabels=FALSE, nlevels=10, col=rgb(153,112,171, maxColorValue = 255), add=TRUE)
contour(z_clust3, drawlabels=FALSE, nlevels=10, col=rgb(90,174,97, maxColorValue = 255), add=TRUE)
contour(z_clust4, drawlabels=FALSE, nlevels=10, col=rgb(254,224,144, maxColorValue = 255), add=TRUE)
contour(z_clust5, drawlabels=FALSE, nlevels=10, col=rgb(69,117,180, maxColorValue = 255), add=TRUE)
legend(-0.5, 0.2, legend=c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5"),
       col=c(rgb(215,48,39, maxColorValue = 255),
             rgb(153,112,171, maxColorValue = 255),
             rgb(90,174,97, maxColorValue = 255),
             rgb(254,224,144, maxColorValue = 255),
             rgb(69,117,180, maxColorValue = 255)),
       lty=c(1,1,1,1,1), cex=0.8)
dev.off()

#
clust1_second_cc <- sapply(1:ncol(clust1_cca$Lambda), function(s){
  clust1_cca$Q1[, 2, s]%*%t(clust1_cca$Q2[, 2, s])
}) %>%
  t()
clust2_second_cc <- sapply(1:ncol(clust2_cca$Lambda), function(s){
  clust2_cca$Q1[, 2, s]%*%t(clust2_cca$Q2[, 2, s])
}) %>%
  t()
clust3_second_cc <- sapply(1:ncol(clust3_cca$Lambda), function(s){
  clust3_cca$Q1[, 2, s]%*%t(clust3_cca$Q2[, 2, s])
}) %>%
  t()
clust4_second_cc <- sapply(1:ncol(clust4_cca$Lambda), function(s){
  clust4_cca$Q1[, 2, s]%*%t(clust4_cca$Q2[, 2, s])
}) %>%
  t()
clust5_second_cc <- sapply(1:ncol(clust5_cca$Lambda), function(s){
  clust5_cca$Q1[, 2, s]%*%t(clust5_cca$Q2[, 2, s])
}) %>%
  t()

#
pdf("figures_tables/brazil_group_similarities_posterior2.pdf", family = "Times", width = 6, height = 6)
par(mfrow = c(1, 1))
second_cc_axes <- svd(rbind(clust1_second_cc, clust2_second_cc, clust3_second_cc, clust4_second_cc, clust5_second_cc))
z_clust1 <- MASS::kde2d(second_cc_axes[["u"]][1:nrow(clust1_second_cc), 1]*second_cc_axes[["d"]][1],
                        second_cc_axes[["u"]][1:nrow(clust1_second_cc), 2]*second_cc_axes[["d"]][2],
                        n=50)
z_clust2 <- MASS::kde2d(second_cc_axes[["u"]][(nrow(clust1_second_cc)+1):(2*nrow(clust2_second_cc)), 1]*second_cc_axes[["d"]][1],
                        second_cc_axes[["u"]][(nrow(clust1_second_cc)+1):(2*nrow(clust2_second_cc)), 2]*second_cc_axes[["d"]][2],
                        n=50)
z_clust3 <- MASS::kde2d(second_cc_axes[["u"]][(2*nrow(clust1_second_cc)+1):(3*nrow(clust2_second_cc)), 1]*second_cc_axes[["d"]][1],
                        second_cc_axes[["u"]][(2*nrow(clust1_second_cc)+1):(3*nrow(clust2_second_cc)), 2]*second_cc_axes[["d"]][2],
                        n=50)
z_clust4 <- MASS::kde2d(second_cc_axes[["u"]][(3*nrow(clust1_second_cc)+1):(4*nrow(clust2_second_cc)), 1]*second_cc_axes[["d"]][1],
                        second_cc_axes[["u"]][(3*nrow(clust1_second_cc)+1):(4*nrow(clust2_second_cc)), 2]*second_cc_axes[["d"]][2],
                        n=50)
z_clust5 <- MASS::kde2d(second_cc_axes[["u"]][(4*nrow(clust1_second_cc)+1):(5*nrow(clust2_second_cc)), 1]*second_cc_axes[["d"]][1],
                        second_cc_axes[["u"]][(4*nrow(clust1_second_cc)+1):(5*nrow(clust2_second_cc)), 2]*second_cc_axes[["d"]][2],
                        n=50)
plot(second_cc_axes[["u"]][, 1]*second_cc_axes[["d"]][1], second_cc_axes[["u"]][, 2]*second_cc_axes[["d"]][2],
     col = c(rep(rgb(215,48,39, maxColorValue = 255, alpha = 50), nrow(clust1_second_cc)),
             rep(rgb(153,112,171, maxColorValue = 255, alpha = 50), nrow(clust2_second_cc)),
             rep(rgb(90,174,97, maxColorValue = 255, alpha = 50), nrow(clust3_second_cc)),
             rep(rgb(254,224,144, maxColorValue = 255, alpha = 50), nrow(clust4_second_cc)),
             rep(rgb(69,117,180, maxColorValue = 255, alpha = 50), nrow(clust5_second_cc))),
     xlab = "PC1", ylab = "PC2")
contour(z_clust1, drawlabels=FALSE, nlevels=10, col=rgb(215,48,39, maxColorValue = 255), add=TRUE)
contour(z_clust2, drawlabels=FALSE, nlevels=10, col=rgb(153,112,171, maxColorValue = 255), add=TRUE)
contour(z_clust3, drawlabels=FALSE, nlevels=10, col=rgb(90,174,97, maxColorValue = 255), add=TRUE)
contour(z_clust4, drawlabels=FALSE, nlevels=10, col=rgb(254,224,144, maxColorValue = 255), add=TRUE)
contour(z_clust5, drawlabels=FALSE, nlevels=10, col=rgb(69,117,180, maxColorValue = 255), add=TRUE)
legend(-0.6, 0.7, legend=c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5"),
       col=c(rgb(215,48,39, maxColorValue = 255),
             rgb(153,112,171, maxColorValue = 255),
             rgb(90,174,97, maxColorValue = 255),
             rgb(254,224,144, maxColorValue = 255),
             rgb(69,117,180, maxColorValue = 255)),
       lty=c(1,1,1,1,1), cex=0.8)
dev.off()

pdf("figures_tables/brazil_weather_stations.pdf", height = 6, width = 6)
plot(station_data$Longitude, station_data$Latitude, col = c(rgb(215,48,39, maxColorValue = 255),
                                                            rgb(153,112,171, maxColorValue = 255),
                                                            rgb(90,174,97, maxColorValue = 255),
                                                            rgb(254,224,144, maxColorValue = 255),
                                                            rgb(69,117,180, maxColorValue = 255))[station_data$cluster],
     xlab = "Longitude", ylab = "Latitude")
legend(-70, -20, legend=c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5"),
       col=c(rgb(215,48,39, maxColorValue = 255),
             rgb(153,112,171, maxColorValue = 255),
             rgb(90,174,97, maxColorValue = 255),
             rgb(254,224,144, maxColorValue = 255),
             rgb(69,117,180, maxColorValue = 255)),
       lty=c(1,1,1,1,1), cex=0.8)
dev.off()

#
c12_dot <- (clust1_first_cc*clust2_first_cc[sample(1:nrow(clust2_first_cc), nrow(clust2_first_cc)), ]) %>%
           rowSums()
c13_dot <- (clust1_first_cc*clust3_first_cc[sample(1:nrow(clust3_first_cc), nrow(clust3_first_cc)), ]) %>%
           rowSums()
c14_dot <- (clust1_first_cc*clust4_first_cc[sample(1:nrow(clust4_first_cc), nrow(clust4_first_cc)), ]) %>%
           rowSums()
c15_dot <- (clust1_first_cc*clust5_first_cc[sample(1:nrow(clust5_first_cc), nrow(clust5_first_cc)), ]) %>%
           rowSums()
c23_dot <- (clust2_first_cc*clust3_first_cc[sample(1:nrow(clust3_first_cc), nrow(clust3_first_cc)), ]) %>%
           rowSums()
c24_dot <- (clust2_first_cc*clust4_first_cc[sample(1:nrow(clust4_first_cc), nrow(clust4_first_cc)), ]) %>%
           rowSums()
c25_dot <- (clust2_first_cc*clust5_first_cc[sample(1:nrow(clust5_first_cc), nrow(clust5_first_cc)), ]) %>%
           rowSums()
c34_dot <- (clust3_first_cc*clust4_first_cc[sample(1:nrow(clust4_first_cc), nrow(clust4_first_cc)), ]) %>%
           rowSums()
c35_dot <- (clust3_first_cc*clust5_first_cc[sample(1:nrow(clust5_first_cc), nrow(clust5_first_cc)), ]) %>%
           rowSums()
c45_dot <- (clust4_first_cc*clust5_first_cc[sample(1:nrow(clust5_first_cc), nrow(clust5_first_cc)), ]) %>%
           rowSums()

dot_mat <- rbind(c(paste0(round(mean(c12_dot), 3), ", [", round(HDInterval::hdi(c12_dot, 0.95)[1], 3), ", ", round(HDInterval::hdi(c12_dot, 0.95)[2], 3), "]"),
                   paste0(round(mean(c13_dot), 3), ", [", round(HDInterval::hdi(c13_dot, 0.95)[1], 3), ", ", round(HDInterval::hdi(c13_dot, 0.95)[2], 3), "]"),
                   paste0(round(mean(c14_dot), 3), ", [", round(HDInterval::hdi(c14_dot, 0.95)[1], 3), ", ", round(HDInterval::hdi(c14_dot, 0.95)[2], 3), "]"),
                   paste0(round(mean(c15_dot), 3), ", [", round(HDInterval::hdi(c15_dot, 0.95)[1], 3), ", ", round(HDInterval::hdi(c15_dot, 0.95)[2], 3), "]")),
                 c(".",
                   paste0(round(mean(c23_dot), 3), ", [", round(HDInterval::hdi(c23_dot, 0.95)[1], 3), ", ", round(HDInterval::hdi(c23_dot, 0.95)[2], 3), "]"),
                   paste0(round(mean(c24_dot), 3), ", [", round(HDInterval::hdi(c24_dot, 0.95)[1], 3), ", ", round(HDInterval::hdi(c24_dot, 0.95)[2], 3), "]"),
                   paste0(round(mean(c25_dot), 3), ", [", round(HDInterval::hdi(c25_dot, 0.95)[1], 3), ", ", round(HDInterval::hdi(c25_dot, 0.95)[2], 3), "]")),
                 c(".",
                   ".",
                   paste0(round(mean(c34_dot), 3), ", [", round(HDInterval::hdi(c34_dot, 0.95)[1], 3), ", ", round(HDInterval::hdi(c34_dot, 0.95)[2], 3), "]"),
                   paste0(round(mean(c35_dot), 3), ", [", round(HDInterval::hdi(c35_dot, 0.95)[1], 3), ", ", round(HDInterval::hdi(c35_dot, 0.95)[2], 3), "]")),
                 c(".",
                   ".",
                   ".",
                   paste0(round(mean(c45_dot), 3), ", [", round(HDInterval::hdi(c45_dot, 0.95)[1], 3), ", ", round(HDInterval::hdi(c45_dot, 0.95)[2], 3), "]"))) %>%
          magrittr::set_colnames(paste0("Region ", 2:5)) %>%
          magrittr::set_rownames(paste0("Region ", 1:4))
dot_mat %>% xtable::xtable() %>% xtable::autoformat() %>% print(file = "figures_tables/brazil_similarities.txt")

