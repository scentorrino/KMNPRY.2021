#' hpj.fe
#'
#' A function to perform half-panel jacknife estimation of ARDL panel data models
#' using the plm package.
#'
#' @param formula A formula specifying the model to be estimated.
#' @param dataset A pdata frame that contains the data for the analysis.
#' @param effect A character string specifying the type of effect ("individual" or "time").
#' @param model A character string specifying the model type ("within", "random", "pooling").
#'        Default is "within".
#' @param compute.std.err A logical indicating whether to compute standard errors. Default is TRUE.
#'
#' @return A list containing coefficients, standard errors, fitted values, residuals, 
#'         variance-covariance matrix, and the data set used.
#' @export
hpj.fe <- function(formula, dataset, effect, model = "within", compute.std.err = TRUE) {
  
  # Check if the dataset is a pdata.frame
  if(class(dataset)[1] != "pdata.frame")
    stop(paste("Error: Data should be of class pdata.frame with individual and time indices.", "\n", sep = ""))
  
  # Extract and clean column order from the formula
  col.order <- gsub(" ", "", unlist(strsplit(as.character(formula)[3], "[+]")))
  col.order <- col.order[which(col.order %in% names(dataset))]
  
  # Define dependent and independent variables
  depvar <- dataset[, names(dataset) %in% as.character(formula)[2]]
  indvar <- dataset[, names(dataset) %in% col.order]
  
  ## Reorder independent variables to match the order in the formula
  col.index <- as.numeric(sapply(paste0("^", col.order, "$"), grep, colnames(indvar)))
  indvar <- indvar[, col.index]
  
  # Subset the dataset to remove NAs
  dataset <- subset(dataset, !is.na(depvar) & apply(!is.na(indvar), 1, prod) == 1)
  
  # Split the panel data based on the odd/even nature of T
  i.index <- index(dataset)[, 1]
  t.index <- index(dataset)[, 2]
  
  i.panel1 <- NULL
  i.panel2 <- NULL
  for(ii in unique(i.index)) {
    t.indexii <- length(t.index[i.index == ii])
    if(numbers::mod(t.indexii, 2) == 1) {
      i.panel1 <- c(i.panel1, 0, array(1, c((t.indexii - 1) / 2, 1)),
                    array(0, c((t.indexii - 1) / 2, 1)))
      i.panel2 <- c(i.panel2, 0, array(0, c((t.indexii - 1) / 2, 1)),
                    array(1, c((t.indexii - 1) / 2, 1)))
    } else {
      i.panel1 <- c(i.panel1, array(1, c(t.indexii / 2, 1)), 
                    array(0, c(t.indexii / 2, 1)))
      i.panel2 <- c(i.panel2, array(0, c(t.indexii / 2, 1)), 
                    array(1, c(t.indexii / 2, 1)))
    }
  }
  
  dataset1 <- plm::pdata.frame(subset(dataset, i.panel1 + i.panel2 == 1), index = names(index(dataset)))
  full.estimation <- plm::plm(formula, data = dataset1, effect = effect, model = model)
  
  panel.1 <- subset(dataset, i.panel1 == 1)
  panel.2 <- subset(dataset, i.panel2 == 1)
  
  ### Obtain the two subset estimators separately
  p1.estimation <- plm::plm(formula, data = panel.1, effect = effect, model = model)
  p2.estimation <- plm::plm(formula, data = panel.2, effect = effect, model = model)
  
  hpj.coeffs <- 2 * full.estimation$coefficients - 0.5 * (p1.estimation$coefficients + p2.estimation$coefficients)
  
  # Redefine the main elements to compute covariance matrix
  depvar <- dataset1[, names(dataset1) %in% as.character(formula)[2]]
  indvar <- dataset1[, names(dataset1) %in% gsub(" ", "", unlist(strsplit(as.character(formula)[3], "[+]")))]
  indvar <- indvar[, col.index]
  
  i.index <- as.numeric(index(dataset1)[, 1])
  t.index <- as.numeric(index(dataset1)[, 2])
  
  i.panel1new <- i.panel1[(i.panel1 + i.panel2) == 1]
  i.panel2new <- i.panel2[(i.panel1 + i.panel2) == 1]
  
  N <- length(unique(i.index))
  T <- length(unique(t.index))
  
  ## Compute overall means 
  if(!identical(model, "pooling")) {
    Xast <- matrix(NA, dim(indvar)[1], dim(indvar)[2])
    for(jj in 1:dim(indvar)[2])
      Xast[, jj] <- plm::Within(indvar[, jj], effect = effect, na.rm = TRUE)
    
    Yast <- plm::Within(depvar, effect = effect, na.rm = TRUE)
    
    fitted <- c(Xast %*% hpj.coeffs)
    resid <- c(Yast - Xast %*% hpj.coeffs)
  }
  
  if(isTRUE(compute.std.err)) {
    QM <- t(Xast) %*% Xast
    
    ## FIRST SUBSAMPLE
    Xast.1 <- matrix(NA, dim(indvar[i.panel1new == 1,])[1], dim(indvar)[2])
    for(jj in 1:dim(Xast)[2])
      Xast.1[, jj] <- plm::Within(indvar[i.panel1new == 1, jj], effect = effect, na.rm = TRUE)
    
    ### SECOND SUBSAMPLE
    Xast.2 <- matrix(NA, dim(indvar[i.panel2new == 1,])[1], dim(indvar)[2])
    for(jj in 1:dim(Xast)[2])
      Xast.2[, jj] <- plm::Within(indvar[i.panel2new == 1, jj], effect = effect, na.rm = TRUE)
    
    M.dast <- matrix(NA, dim(dataset1)[1], dim(indvar)[2])
    M.dast[i.panel1new == 1, ] <- Xast.1 + kronecker(array(1,c(dim(Xast.1)[1],1)),t(as.matrix(2*apply(indvar,2,mean)  - apply(indvar[i.panel1new == 1,],2,mean))))
    M.dast[i.panel2new == 1, ] <- Xast.2 + kronecker(array(1,c(dim(Xast.2)[1],1)),t(as.matrix(2*apply(indvar,2,mean)  - apply(indvar[i.panel2new == 1,],2,mean))))
    
    Dast <- 2 * Xast - M.dast
    
    vcovmat <- solve(QM) %*% (t(Dast) %*% diag(resid^2) %*% Dast) %*% solve(QM)
    std.err <- sqrt(diag(vcovmat))
  } else {
    vcovmat <- NULL
    std.err <- NULL
  }
  
  # Return results as a list
  tab.reg <- list("coeff" = hpj.coeffs, "SE" = std.err, "fitted" = fitted, "resid" = resid, "vcov" = vcovmat, "data" = dataset1)
  
  return(tab.reg)
}