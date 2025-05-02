#' Generate Bootstrap for the Psi Coefficients
#'
#' This function performs a dynamic wild bootstrap to generate estimates for the
#' psi coefficients based on the provided dataset and ARDL model.
#'
#' @param p.dataset A pdata.frame containing the dataset for the analysis.
#' @param formula.ardl A formula specifying the ARDL model.
#' @param thetaest A vector of estimated coefficients from the model.
#' @param residuals A vector of residuals from the model.
#' @param n.years An integer specifying the number of years for which to generate coefficients. Default is 101.
#' @param B An integer specifying the number of bootstrap repetitions. Default is 99.
#' 
#' @return A list containing:
#' \item{psihat}{A matrix of psi coefficients for each bootstrap sample.}
#' \item{thetahat}{A matrix of estimated coefficients for each bootstrap sample.}
#' 
#' @import plm
#' @export
psicoef.b <- function(p.dataset, formula.ardl, thetaest, residuals, n.years = 101, B = 99) {
  
  col.order <- gsub(" ", "", unlist(strsplit(as.character(formula.ardl)[3], "[+]")))
  orig.names <- names(p.dataset)
  orig.index <- names(index(p.dataset))
  
  psihat.b <- matrix(0, n.years, B)
  thetahat <- matrix(NA, B, length(col.order))
  
  indvar <- p.dataset[, names(p.dataset) %in% gsub(" ", "", unlist(strsplit(as.character(formula.ardl)[3], "[+]")))]
  col.index <- as.numeric(sapply(paste0("^", col.order, "$"), grep, colnames(indvar)))
  indvar <- indvar[, col.index]
  
  i.index <- as.numeric(index(p.dataset)[, 1])
  num.c <- max(i.index)
  num.T <- max(as.numeric(index(p.dataset)[, 2]))
  
  for (jj in 1:B) {
    
    ### DYNAMIC WILD BOOTSTRAP TO ACCOUNT FOR POTENTIAL CROSS-SECTIONAL CORRELATION
    tmpdb <- p.dataset
    tmpdb$y_star <- NA
    
    for (kk in 1:length(unique(i.index))) {
      
      tmp.resid <- residuals[i.index == kk]
      tmp.indvar <- indvar[i.index == kk, ]
      tmp.indvar[is.na(tmp.indvar)] <- 0
      
      if (sum(is.na(tmp.resid)) < num.T) {
        new.res <- DWB(tmp.resid[!is.na(tmp.resid)])
        
        ytmpboot <- tmpdb$y_star[i.index == kk]
        for (ii in 1:length(ytmpboot)) {
          ytmpboot[ii] <- ifelse(is.na(tmp.resid[ii]), NA, sum(tmp.indvar[ii, ] * thetaest) + new.res[ii])
          if (ii <= num.T - 1) {
            tmp.indvar$l1growth[ii + 1] <- ytmpboot[ii]
          }
          if (ii <= num.T - 2) {
            tmp.indvar$l2growth[ii + 2] <- ytmpboot[ii]
          }
          if (ii <= num.T - 3) {
            tmp.indvar$l3growth[ii + 3] <- ytmpboot[ii]
          }
          if (ii <= num.T - 4) {
            tmp.indvar$l4growth[ii + 4] <- ytmpboot[ii]
          }
        }
        
        tmpdb$y_star[i.index == kk] <- ytmpboot 
      }
    }
    
    # Remove the dependent variable and create its lags
    tmpdb$growth <- tmpdb$y_star
    tmpdb$l1growth <- plm::lag(tmpdb$growth)
    tmpdb$l2growth <- plm::lag(tmpdb$l1growth)
    tmpdb$l3growth <- plm::lag(tmpdb$l2growth)
    tmpdb$l4growth <- plm::lag(tmpdb$l3growth)
    
    tmpdb <- pdata.frame(tmpdb, index = orig.index)
    hpj.p.1.b <- hpj.fe(formula.ardl, tmpdb, effect = "individual", compute.std.err = FALSE)
    
    thetahat[jj, ] <- hpj.p.1.b$coeff
    psihat.b[, jj] <- psicoef(hpj.p.1.b, n.years = n.years)
  }
  
  return(list("psihat" = psihat.b, "thetahat" = thetahat))
}
