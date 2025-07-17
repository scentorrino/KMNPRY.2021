#' Dependent Wild Bootstrap
#'
#' This function performs dependent wild bootstrap on the input data using
#' the approach outlined in Gao et al. It utilizes a Bartlett kernel and
#' computes the necessary adjustments for the bootstrap procedure.
#'
#' @param y A numeric vector of observations for which the bootstrap is to be performed.
#' 
#' @return A numeric vector of bootstrapped values.
#' 
#' @import GCCfactor
#' @import stats
#' @export
DWB <- function(y) {
  
  T <- length(y)
  y_star <- rep(0, T)
  
  # Use the approach in Gao et al.
  # As you are using a Bartlett kernel, q = 1 and cq = 1
  # and int a(x)^2dx = 2/3
  # Truncation parameter
  QT <- ceiling(T^(2/9))
  nuq <- 1/3
  
  # Sequence of time points
  ss <- seq(1, T)
  afun <- matrix(0, T, T)
  
  for (tt in 1:T) {
    afun[, tt] <- sapply(ss, function(x) GCCfactor::Bartlett((x - ss[tt]) / T^nuq))
  }
  
  # Partial covariances of the error
  Delta1 <- 0
  
  for (k in 1:QT) {
    Delta1 <- Delta1 + 2 * (k / T) * mean(y[1:(T - k)] * y[(k + 1):T])
  }
  
  Delta2 <- (2 / 3) * (sum(apply((y %o% y) * afun, 1, mean)))^2
  lopt <- max(ceiling((T * Delta1^2 / Delta2)^(1 / nuq)), 10)
  
  Omega <- matrix(0, T, T)
  
  for (tt in 1:T) {
    Omega[, tt] <- sapply(ss, function(x) GCCfactor::Bartlett((x - ss[tt]) / lopt))
  }
  
  eig <- eigen(Omega)
  Omega_half <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  zeta <- as.vector(Omega_half %*% matrix(stats::rnorm(T, 0, 1), ncol = 1))
  
  y_star <- zeta * y
  
  return(y_star)
}
