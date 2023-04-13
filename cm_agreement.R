#
ns <- c(10, 50, 200, 500)
p <- 2
S <- 100
indep_results <- matrix(NA, length(ns), S)
cm_results <- matrix(NA, length(ns), S)
pert_results <- matrix(NA, length(ns), S)

for(i in 1:length(ns)){
  n <- ns[i]
  for(s in 1:S){
    # Independent rows
    X <- matrix(rnorm(n*p), n)
    Y <- matrix(rnorm(n*p), n)
    assignment <- transport::transport(transport::pp(X), transport::pp(Y), p = 2)
    is_cm <- all(assignment[, 2] == 1:n)
    indep_results[i, s] <- is_cm == cmcca::is_cm(X, Y)

    # Cyclically monotone transformation
    X <- matrix(rnorm(n*p), n)
    Sigma <- rWishart(1, p, diag(rep(1, p)))[,,1]
    Y <- X%*%Sigma
    assignment <- transport::transport(transport::pp(X), transport::pp(Y), p = 2)
    is_cm <- all(assignment[, 2] == 1:n)
    cm_results[i, s] <- is_cm == cmcca::is_cm(X, Y)

    # Small perturbation
    Y[1, ] <- Y[1, ] + rnorm(1, sd = 1/sqrt(n))
    assignment <- transport::transport(transport::pp(X), transport::pp(Y), p = 2)
    is_cm <- all(assignment[, 2] == 1:n)
    pert_results[i, s] <- is_cm == cmcca::is_cm(X, Y)
  }
  print(i)
}
