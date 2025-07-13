#' Compute Long-Run Heterogeneous Coefficients
#'
#' This function computes long-run coefficients and their standard errors from 
#' heterogeneous panel data models. It processes coefficient matrices from different
#' panel estimations and applies bias correction methods to derive the long-run
#' multiplier and its variance-covariance matrix.
#'
#' @param het.coeff A list containing heterogeneous coefficient estimates with the following structure:
#'   \describe{
#'     \item{hpj_coeff}{Matrix of HPJ (Heterogeneous Panel J-estimator) coefficients}
#'     \item{all_coeff}{List containing coefficient matrices from different panel estimations:
#'       \describe{
#'         \item{full}{Full panel coefficient matrix}
#'         \item{panel.1}{First panel subset coefficient matrix}
#'         \item{panel.2}{Second panel subset coefficient matrix}
#'       }
#'     }
#'   }
#' @param country_iso A vector containing the iso3 codes of countries, 
#' @param p Integer. Number of autoregressive lags in the model. Default is 1.
#' @param q Integer. Number of explanatory variable lags in the model. Default is 4.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{lr_coeff}{Numeric. The bias-corrected long-run coefficient (theta hat)}
#'     \item{std_err}{Matrix. Standard errors of the long-run coefficient estimates}
#'   }
#'
#' @details 
#' The function performs the following steps:
#' \enumerate{
#'   \item Filters coefficient matrices to include only countries present in \code{country.vars$iso}
#'   \item Extracts autoregressive coefficients (phi) and explanatory variable coefficients (beta)
#'   \item Computes bias-corrected estimates using cross-sectional averages
#'   \item Calculates the long-run multiplier as: theta = beta / (1 - phi)
#'   \item Constructs variance-covariance matrix from the three panel estimations
#'   \item Applies delta method for standard error calculation of the long-run coefficient
#' }
#'
#' The variance calculation uses a complex combination of the variance-covariance matrices
#' from the full panel and two subpanels, accounting for cross-correlations between
#' different estimation methods.
#'
#' @note 
#' This function assumes a specific structure for the coefficient matrices where:
#' \itemize{
#'   \item Row 1 contains intercept terms
#'   \item Rows 2 to (p+1) contain autoregressive coefficients
#'   \item Rows (p+2) to (p+q+2) contain explanatory variable coefficients
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage with heterogeneous panel data
#' het_coeffs <- list(
#'   hpj_coeff = matrix(rnorm(50), nrow = 10, ncol = 5),
#'   all_coeff = list(
#'     full = matrix(rnorm(50), nrow = 10, ncol = 5),
#'     panel.1 = matrix(rnorm(50), nrow = 10, ncol = 5),
#'     panel.2 = matrix(rnorm(50), nrow = 10, ncol = 5)
#'   )
#' )
#' 
#' country_iso <- c("USA", "DEU", "JPN", "GBR", "FRA")
#' 
#' 
#' # Compute long-run coefficients
#' lr_results <- lr_het_coef(het_coeffs, country_iso, p = 1, q = 4)
#' 
#' # Access results
#' print(lr_results$lr_coeff)  # Long-run coefficient
#' print(lr_results$std_err)   # Standard errors
#' }
#'
#'
#' @export
#' 
lr_het_coef <- function(het.coeff, subgroup = NULL, ma = 30, p = 1, q = 4) {
  
  if(!is.null(subgroup)){
    het.coeff$hpj_coeff <- het.coeff$hpj_coeff[,colnames(het.coeff$hpj_coeff) %in% subgroup]
    het.coeff$all_coeff$full <- het.coeff$all_coeff$full[,colnames(het.coeff$all_coeff$full) %in% subgroup]
    
    het.coeff$all_coeff$panel.1 <- het.coeff$all_coeff$panel.1[,colnames(het.coeff$all_coeff$panel.1) %in% subgroup]
    het.coeff$all_coeff$panel.2 <- het.coeff$all_coeff$panel.2[,colnames(het.coeff$all_coeff$panel.2) %in% subgroup]
  }
  
  # Extract objects 
  phihatbc  <- sum(apply(matrix(het.coeff$hpj_coeff[2:(p + 1),],nrow = p, ncol = ncol(het.coeff$hpj_coeff)),1,mean))
  betahatbc <- ((ma+1)/2)*sum(apply(het.coeff$hpj_coeff[(p + 2):(p+q+2),],1,mean))
  
  thetahatbc <- betahatbc/(1 - phihatbc)
  
  # Now compute standard errors 
  vcov     <- cov(cbind(t(het.coeff$all_coeff$full[2:(p+q+2),]),
                    t(het.coeff$all_coeff$panel.1[2:(p+q+2),]),
                    t(het.coeff$all_coeff$panel.2[2:(p+q+2),])))
  
  totsize  <- dim(vcov)[1]/3
  vcov.all <- vcov[c(1:totsize),c(1:totsize)]
  vcov.1   <- vcov[c((totsize+1):(2*totsize)),c((totsize+1):(2*totsize))]
  vcov.2   <- vcov[c((2*totsize+1):(3*totsize)),c((2*totsize+1):(3*totsize))]
  
  cov.all.1<- vcov[c(1:totsize),c((totsize+1):(2*totsize))]
  cov.all.2<- vcov[c(1:totsize),c((2*totsize+1):(3*totsize))]
  cov.1.2  <- vcov[c((totsize+1):(2*totsize)),c((2*totsize+1):(3*totsize))]
  
  D <- matrix(c(thetahatbc/(1 - phihatbc),rep(1/(1-phihatbc),5)),ncol = 1)

  var.comp <- 4*vcov.all + 0.25*(vcov.1 + vcov.2 + cov.1.2 + t(cov.1.2)) - (cov.all.1 + t(cov.all.1) + cov.all.2 + t(cov.all.2))

  # Return coefficients
  tab.coef <- list(
                   "lr_coeff"= thetahatbc, 
                   "std_err" = sqrt(as.numeric(t(D) %*% var.comp %*% D)),
                   "vcov"    = var.comp,
                   "jacobian"= D 
                   )
  
  return(tab.coef)
}