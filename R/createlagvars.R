#' Create Lagged Variables
#'
#' This function creates lagged variables for various growth and precipitation metrics,
#' as well as computes counterfactual variables based on temperature differences.
#'
#' @param initdata A pdata.frame containing the data for which lagged variables will be created.
#' @param ma Number of years for full adaptation. Default to 30.
#' 
#' @return A pdata.frame with additional lagged variables and aggregate world GDP.
#' 
#' @import plm
#' @export
createlagvars <- function(initdata, ma=30) {
  
  # Check if the input data is a pdata.frame
  if (class(initdata)[1] != "pdata.frame") {
    stop("Data should be put in pdata.frame format")
  }
  
  # CREATE NEW LAGGED VARIABLES
  initdata$l2growth <- plm::lag(initdata$l1growth)
  initdata$l3growth <- plm::lag(initdata$l2growth)
  initdata$l4growth <- plm::lag(initdata$l3growth)
  
  initdata$l1dprecip_plus <- plm::lag(initdata$dprecip_plus)
  initdata$l2dprecip_plus <- plm::lag(initdata$l1dprecip_plus)
  initdata$l3dprecip_plus <- plm::lag(initdata$l2dprecip_plus)
  initdata$l4dprecip_plus <- plm::lag(initdata$l3dprecip_plus)
  
  initdata$l1dprecip_minus <- plm::lag(initdata$dprecip_minus)
  initdata$l2dprecip_minus <- plm::lag(initdata$l1dprecip_minus)
  initdata$l3dprecip_minus <- plm::lag(initdata$l2dprecip_minus)
  initdata$l4dprecip_minus <- plm::lag(initdata$l3dprecip_minus)
  
  # CREATE ABSOLUTE TEMPERATURE DIFFERENCES
  initdata$temp_abs <- abs(initdata$temp - initdata$temp_hn)
  initdata$dtempabs <- initdata$temp_abs - plm::lag(initdata$temp_abs)
  initdata$l1dtempabs <- plm::lag(initdata$dtempabs)
  initdata$l2dtempabs <- plm::lag(initdata$l1dtempabs)
  initdata$l3dtempabs <- plm::lag(initdata$l2dtempabs)
  initdata$l4dtempabs <- plm::lag(initdata$l3dtempabs)
  
  # # Obtain real per-capita World GDP growth
  # imf_datatools <- reticulate::import("imf_datatools")
  # 
  # # Metadata 
  # # co.names <- imf_datatools$ecos_sdmx_utilities$get_countries('WEO_WEO_PUBLISHED')
  # # var.names <- imf_datatools$ecos_sdmx_utilities$get_ecos_sdmx_metadata('WEO_WEO_PUBLISHED')
  # world.gdp.growth <- imf_datatools$get_ecos_sdmx_data('WEO_WEO_PUBLISHED', '001', 'NGDP_R_PPP_PC_PCH')
  # 
  # # Rename and adjust the world GDP growth data
  # names(world.gdp.growth) <- "l1wgrowth"
  # world.gdp.growth$year <- lubridate::year(rownames(world.gdp.growth)) - 1
  # world.gdp.growth$l1wgrowth <- world.gdp.growth$l1wgrowth / 100
  # 
  # numeric.year <- as.numeric(levels(initdata$year))
  # 
  # # Create a temporary data frame to hold missing GDP growth values
  # temp.dat <- data.frame("l1wgrowth" = rep(NA, length(min(numeric.year):(min(world.gdp.growth$year) - 1))),
  #                        "year" = min(numeric.year):(min(world.gdp.growth$year) - 1))
  # 
  # # Combine the temporary data with the world GDP growth data
  # world.gdp.growth <- rbind(temp.dat, world.gdp.growth[world.gdp.growth$year < 2015, ])
  # 
  # # Assign the GDP growth values to the initdata
  # initdata$l1wgrowth <- rep(world.gdp.growth$l1wgrowth, length(unique(index(initdata)[, 1])))
  
  return(initdata)
}