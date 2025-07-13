#' hpj.het
#'
#' A function to perform half-panel jacknife estimation of ARDL heterogenous panel data models
#' using the plm package.
#'
#' @param formula A formula specifying the model to be estimated.
#' @param dataset A pdata frame that contains the data for the analysis.
#'
#' @return A list containing the half-panel jackknife heterogeneous estimates and the estimates 
#' from the full and each of the half panels.
#' 
#' @export
#' 
hpj.het <- function(formula, dataset) {
  
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
  
  full.estimation <- plm::pmg(formula, data = dataset1, model = "mg")
  
  panel.1 <- subset(dataset, i.panel1 == 1)
  panel.2 <- subset(dataset, i.panel2 == 1)
  
  ### Obtain the two subset estimators separately
  p1.estimation <- plm::pmg(formula, data = panel.1, model = "mg")
  p2.estimation <- plm::pmg(formula, data = panel.2, model = "mg")
  
  hpj.coeffs <- 2 * full.estimation$indcoef - 0.5 * (p1.estimation$indcoef + p2.estimation$indcoef)
  
  # Return coefficients
  tab.coef <- list(
                   "hpj_coeff" = hpj.coeffs, 
                   "all_coeff" = list("full" = full.estimation$indcoef, 
                                      "panel.1" = p1.estimation$indcoef, 
                                      "panel.2" = p2.estimation$indcoef)
                   )
  
  return(tab.coef)
}