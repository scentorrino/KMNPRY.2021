#' Generate Coefficients of the MA(∞) Representation of the ARDL Model
#'
#' This function computes the coefficients for the MA(∞) representation of the ARDL model
#' based on the input coefficients from a previous estimation.
#'
#' @param hpj.p.1 A list containing the coefficients from the ARDL model. 
#'        It should have a named component `coeff` with at least 9 coefficients.
#' @param n.years An integer specifying the number of years for which to generate coefficients.
#'        Default is 101.
#' 
#' @return A vector of coefficients for the MA(∞) representation for the specified number of years.
#' 
#' @export
psicoef <- function(hpj.p.1, n.years = 101) {
  
  ## Define weights for counterfactual scenarios
  phihat  <- hpj.p.1$coeff[1:4]
  betahat <- hpj.p.1$coeff[5:9]

  psihat  <- integer(n.years)
  psihat[1] <- betahat[1]
  psihat[2] <- betahat[2] + phihat[1] * psihat[1]
  psihat[3] <- betahat[3] + phihat[2] * psihat[1] + phihat[1] * psihat[2]
  psihat[4] <- betahat[4] + phihat[3] * psihat[1] + phihat[2] * psihat[2] + phihat[1] * psihat[3]
  psihat[5] <- betahat[5] + phihat[4] * psihat[1] + phihat[3] * psihat[2] + phihat[2] * psihat[3] + phihat[1] * psihat[4]

  for (j in 6:n.years) {
    psihat[j] <- phihat[4] * psihat[j - 4] + phihat[3] * psihat[j - 3] + phihat[2] * psihat[j - 2] + phihat[1] * psihat[j - 1]
  }
  
  return(psihat)
}