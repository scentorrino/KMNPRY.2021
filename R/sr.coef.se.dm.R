#' Calculate long-run model coefficients and their Standard Errors using the Delta Method
#'
#' This function calculates the long-run coefficients and their standard errors 
#' using the delta method based on a fitted model object. It takes a fitted model 
#' as input and computes the transformed coefficients as well as their standard errors.
#'
#' @param fitted.model A fitted model object that contains the coefficients and 
#' variance-covariance matrix (vcov) of the model.
#' @param p Order of the auto-regressive part. Default to 4.
#' @param q Lags of the exogenous variables. Default to p + 1.
#' @param numreg The number of regressors for which we wish to compute the long-run impacts and their .
#' standard errors. Default to 2.
#'
#' @return A list containing:
#' \item{lr.coef}{A vector of long-run coefficients.}
#' \item{lr.coef.se}{A vector of standard errors for the long-run coefficients.}
#'
#' @examples
#' # Assuming `model` is a fitted model object
#' result <- sr.coef.se.dm(model)
#'
#' @export

sr.coef.se.dm <- function(fitted.model,p = 4,q = p + 1,numreg = 4){
  
  ncoef <- length(fitted.model$coeff)
  phi11 <- 1 - sum(fitted.model$coeff[c(1:p)])
  sder.phi11 <- sqrt(array(1,c(1,p)) %*% fitted.model$vcov[c(1:p),c(1:p)] %*% array(1,c(p,1)))
  
  coeff.lr <- sum(fitted.model$coeff[(p+1):(p + q)])/phi11
  jacob.par <- c(fitted.model$coeff[c(1:p)]*sum(fitted.model$coeff[(p+1):(p + q)])/phi11^2,rep(1/phi11,q),rep(0,ncoef - p - q))
  for(kk in 2:numreg){
    coeff.lr <- c(coeff.lr,sum(fitted.model$coeff[(p + (kk-1)*q +1):(p + kk*q)])/phi11)
    jacob.par <- rbind(jacob.par,
                       c(fitted.model$coeff[c(1:p)]*sum(fitted.model$coeff[(p+(kk-1)*q + 1):(p + kk*q)])/phi11^2,rep(0,(kk-1)*q),rep(1/phi11,q),rep(0,ncoef - p - kk*q)))
  }

  coeff.lr <- c(coeff.lr,phi11)
  
  # Standard errors using delta method
  coeff.lr.se <- c(sqrt(diag(jacob.par %*% fitted.model$vcov %*% t(jacob.par))),sder.phi11)

  return(list("lr.coef" = coeff.lr,"lr.coef.se" = coeff.lr.se))
  
}