#' Calculate long-run model coefficients and their Standard Errors from Error Correction Model using the Delta Method
#'
#' This function calculates the long run coefficients of an ARDL and their standard errors 
#' using the delta method based on its error correction form. It takes a fitted model 
#' as input and computes the transformed coefficients as well as their standard errors.
#'
#' @param fitted.model A fitted model object that contains the coefficients and 
#' variance-covariance matrix (vcov) of the model.
#' @param numreg The number of regressors for which we wish to compute the long-run impacts and their .
#' standard errors. Default to 2.
#'
#' @return A list containing:
#' \item{lr.coef}{A vector of log-ratio coefficients.}
#' \item{lr.coef.se}{A vector of standard errors for the log-ratio coefficients.}
#'
#' @examples
#' # Assuming `model` is a fitted model object
#' result <- lr.coef.se.dm(model)
#'
#' @export

ec.coef.se.dm <- function(fitted.model,numreg = 4){
  
  coeff.lr <- c(-fitted.model$coeff[2:(numreg + 1)]/fitted.model$coeff[1],-fitted.model$coeff[1])

  jacob.par <- c(fitted.model$coeff[2]/fitted.model$coeff[1]^2,-1/fitted.model$coeff[1],rep(0,numreg - 1))
  for(kk in 2:numreg){
    jacob.par <- rbind(jacob.par,
                       c(fitted.model$coeff[2]/fitted.model$coeff[1]^2,rep(0,kk-1),-1/fitted.model$coeff[1],rep(0,numreg - kk)))
  }

  # Standard errors using delta method
  coeff.lr.se <- c(sqrt(diag(jacob.par %*% fitted.model$vcov[1:(numreg + 1),1:(numreg + 1)] %*% t(jacob.par))),
                        fitted.model$SE[1])
  
  return(list("lr.coef" = coeff.lr,"lr.coef.se" = coeff.lr.se))
}