#' Generate Counterfactual Scenarios Based on Temperature Changes
#'
#' This function generates counterfactual scenarios given a scenario of changes
#' in temperatures for all countries, based on specified trends and volatilities.
#'
#' @param m An integer representing the number of past observations used in the calculation.
#' @param beta A numeric value representing the base parameter for temperature changes.
#' @param sigma A numeric value representing the standard deviation of temperature changes.
#' @param N An integer representing the number of countries or observations.
#' @param psihat A numeric vector of psi coefficients.
#' @param delta.trend A scalar representing changes in temperature trends for each period.
#' @param H An integer specifying the number of time periods for the counterfactual scenarios. Default is 100.
#' 
#' @return A list containing:
#' \item{Delta}{A matrix of changes in temperatures for the counterfactual scenarios.}
#' \item{mu}{A matrix of mean expected temperature changes.}
#' \item{omega}{A matrix of standard deviations for temperature changes.}
#' \item{g0}{A numeric value representing the closed form for the expectation of the absolute change in temperatures.}
#' 
#' @export
cccounter <- function(m, beta, sigma, N, psihat, vec.new.trend, H = 100) {
  
  # Assume the shocks follow a normal distribution to find 
  # a closed form for the expectation of the absolute change in temperatures
  mu_Ti <- ((m + 1) / 2) * beta
  omega_Ti <- sqrt(sigma^2 * (1 + (1 / m)))
  
  g_0 <- mu_Ti * (pnorm((mu_Ti / omega_Ti)) - pnorm((-mu_Ti / omega_Ti))) + 
    2 * omega_Ti * dnorm((mu_Ti / omega_Ti))
  
  mu_Ti <- matrix(0, N, H)
  gvec <- matrix(0, N, H)
  Delta <- matrix(0, N, H)
  
  for (j in 1:H) {
    
    mu_Ti[, j] <- ((m + 1) / 2) * (beta + j * vec.new.trend)
    
    gvec[, j] <- mu_Ti[, j] * (pnorm((mu_Ti[, j] / omega_Ti)) - 
                                 pnorm((-mu_Ti[, j] / omega_Ti))) + 
      2 * omega_Ti * dnorm((mu_Ti[, j] / omega_Ti))
    
    for (k in 1:j) {
      ## MAIN EFFECTS
      Delta[, j] <- Delta[, j] + psihat[j - k + 1] * (gvec[, k] - g_0)
    }
  }
  
  return(list("Delta" = Delta, "mu" = mu_Ti, "omega" = omega_Ti, "g0" = g_0))
}
