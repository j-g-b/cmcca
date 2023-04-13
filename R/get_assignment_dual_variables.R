#'
#'
#' @export get_assignment_dual_variables
get_assignment_dual_variables <- function(X, Y){
  require(transport)
  # Problem dims
  n <- nrow(Y)
  # Find optimal assignment and optimal dual variables
  assign <- transport::transport(transport::pp(X),
                                 transport::pp(Y),
                                 p = 2, fullreturn = T, method = "networkflow")
  dual_vars <- assign$dual[(n+1):(2*n),]
  c_dual_vars <- apply(X - Y, 1, function(x){sum(x^2)}) - dual_vars
  slack_matrix <- matrix(NA, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      slack_matrix[i, j] <- sum((X[i, ] - Y[j, ])^2) - c_dual_vars[i] - dual_vars[j]
    }
  }
  slack_matrix[slack_matrix < 0] <- 0
  slack_matrix[abs(slack_matrix) < 1e-10] <- 0
  return(list(dual_vars = dual_vars, c_dual_vars = c_dual_vars,
              slack_matrix = slack_matrix))
}
#'
#'
#' @export move_x_cm
move_x_cm <- function(Y, X, S, i, k, mean=0, sd=1, verbose=F){
  n <- nrow(Y)
  # Find constraints implied by row i
  row_const <- sapply(setdiff(1:n, i), function(j){
    disc <- (X[i, k] - Y[j, k])^2 - S[i, j]
    if(disc > 0){
      return(c(-(X[i, k] - Y[j, k]) - sqrt(disc), -(X[i, k] - Y[j, k]) + sqrt(disc)))
    } else {
      return(c(NA, NA))
    }
  })
  na_slots <- apply(row_const, 2, function(x){any(is.na(x))})
  # Find constraints implied by col i
  col_const <- sapply(setdiff(1:n, i), function(l){
    disc <- (X[i, k] - Y[i, k])^2 + S[l, i]
    return(c(-(X[i, k] - Y[i, k]) - sqrt(disc), -(X[i, k] - Y[i, k]) + sqrt(disc)))
  })
  if(!all(na_slots)){
    row_const <- row_const[, !na_slots, drop=F]
    row_const[abs(row_const) < 1e-10] <- 0
    check <- apply(row_const, 2, function(x){min(x) >= 0 | max(x) <= 0}) %>% all()
    if(!check){
      View(row_const)
      stop()
    }
    # Find interval
    l_bound_row <- max(row_const[2, row_const[2, ] <= 0])
    u_bound_row <- min(row_const[1, row_const[1, ] >= 0])
    l_bound_col <- max(col_const[1, ])
    u_bound_col <- min(col_const[2, ])
    #
    eps_interval <- c(max(l_bound_row, l_bound_col), min(u_bound_row, u_bound_col))
  } else {
    # Find interval
    l_bound_col <- max(col_const[1, ])
    u_bound_col <- min(col_const[2, ])
    #
    eps_interval <- c(l_bound_col, u_bound_col)
  }
  #
  if(eps_interval[1] == 0 & eps_interval[2] == 0){
    eps <- 0
  } else {
    eps <- runif(1, min = eps_interval[1], max = eps_interval[2])
    log_r <- dnorm(X[i,k] + eps, mean = mean, sd = 1, log = T) -
             dnorm(X[i,k], mean = mean, sd = 1, log = T)
    if(runif(1) < exp(log_r)){
      eps <- eps
    } else {
      eps <- 0
    }
  }
  if(verbose){
    print(paste0(eps_interval, collapse = ", "))
    print(eps)
    x_seq <- seq(X[i, k] + eps_interval[1], X[i, k] + eps_interval[2], length.out = 1000)
    y_seq <- dnorm(x_seq)
    plot(x_seq, y_seq, type = 'l')
    Sys.sleep(1)
  }
  #
  S[i, ] <- S[i, ] + eps^2 + 2*eps*(X[i, k] - Y[, k])
  S[, i] <- S[, i] - eps^2 - 2*eps*(X[i, k] - Y[i, k])
  X[i, k] <- X[i, k] + eps
  return(list(X = X, S = S, eps_interval = eps_interval))
}
