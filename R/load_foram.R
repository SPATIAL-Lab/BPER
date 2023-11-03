#' Checks and formats foraminiferal proxy data
#'
#' This function allows you to load, QC and format foram data to use in subsequent 'BPER' functions.
#' Function checks which data are present and whether measurement uncertainties are included.
#' If ages or uncertainties are not included the function stops. Function flags 'incomplete'
#' data object - i.e., missing Mg/Ca, or missing d18O, or missing d11B - but proceeds. Function also writes
#' "priors_in_foram.R" to parent directory which needs to be modified before 'priors_foram' function is used.
#' If no modifications to priors are needed, this can be turned off using 'load.priors' argument. 
#'
#' @param foram_data User inputs an 8 column R data object and assigns it to 'foram_data' in arguments. Object columns
#' should be organized as: age, d11B, d11B2s, d11Bspec, Mg/Ca, Mg/Ca2s, d18O, d18O2s. The order matters.
#' 
#' @param load.priors Logical. When TRUE, writes 'priors_in_foram.R' file to parent directory which can be modified 
#' by the user. Defaults to TRUE.
#'
#' @returns Returns list 'load_foram'. This list contains 'prox_in' data.frame which has been checked for data completeness
#' and compatibility with other package functions and 'obs_type' which indicates which measurements are included.
#'
#' @examples
#' load_foram(foram_data)
#' 
#' @export
load_foram <- function(foram_data = foram_data, load.priors = TRUE){
  # Check that the correct number of columns are input - i.e. follows the template
  if (ncol(foram_data) == 8){
    prox_in = data.frame(foram_data)
    prox_in = prox_in[,c(1:8)]
    names(prox_in) = c("age","d11B", "d11B2s", "d11Bspec", "MgCa", "MgCa2s", "d18O", "d18O2s")
  } else{
    stop("Number of columns in data object do not match template - 8 columns are required even if some are empty")
  }

  # Check that input data frame contains age info, 2 sigma uncertainties are included for each observation
  # and warn if d11B, Mg/Ca or d18O measurements are entirely missing
  if(length(prox_in$age) < 1){
    stop("Must include age data")
  }
  if(length(prox_in$d11B) < 1){
    warning("d11B data are missing or out of place")
  }
  if(length(prox_in$MgCa) < 1){
    warning("MgCa data are missing or out of place")
  }
  if(length(prox_in$d18O) < 1){
    warning("d18O data are missing or out of place")
  }
  if(length(prox_in$d11B2s) < 1 & length(prox_in$d11B) < 1){
    stop("Must include d11B2s (2 sigma uncertainty) if d11B data are being used")
  }
  if(length(prox_in$d11Bspec) < 1 & length(prox_in$d11B) < 1){
    warning("d11B vital effect calibration defaults to d11Bforam = d11Bborate if unspecified. Specify: 'Grub', 'Tsac', 'Ouni','custom', or 'borate'")
  }
  if(length(prox_in$d18O2s) < 1 & length(prox_in$d18O) < 1){
    stop("Must include d18O2s (2 sigma uncertainty) if d18O data are being used")
  }
  if(length(prox_in$MgCa2s) < 1 & length(prox_in$MgCa) < 1){
    stop("Must include MgCa2s (2 sigma uncertainty) if Mg/Ca data are being used")
  }

  # Assign obs_type based on which categories of data are input
  if(length(prox_in$d11B) > 0 & length(prox_in$MgCa) > 0 & length(prox_in$d18O) > 0){
    obs_type = 'BMO'
  } else if(length(prox_in$d11B) > 0 & length(prox_in$MgCa) > 0 & length(prox_in$d18O) < 1){
    obs_type = 'BM'
  } else if(length(prox_in$d11B) > 0 & length(prox_in$MgCa) < 1 & length(prox_in$d18O) > 0){
    obs_type = 'BO'
  } else if(length(prox_in$d11B) < 1 & length(prox_in$MgCa) > 0 & length(prox_in$d18O) > 0){
    obs_type = 'MO'
  } else if(length(prox_in$d11B) > 0 & length(prox_in$MgCa) < 1 & length(prox_in$d18O) < 1){
    obs_type = 'B'
  } else if(length(prox_in$d11B) < 1 & length(prox_in$MgCa) > 0 & length(prox_in$d18O) < 1){
    obs_type = 'M'
  } else{
    stop("Must include some type of the following data combinations (B = d11B, O = d18O, Mg = Mg/Ca): B+M+O, B+M, B+O, M+O, B, M")
  }
  
  # load the 'priors_in_foram.R' script for user to fill out for the 'priors_foram' function
  if(load.priors == TRUE){
  priors_in_foram <- system.file("extdata", "priors_in_foram.R", package = "BPER")
  save(priors_in_foram, "priors_in_foram.R")
  }
  
  load_foram = list(prox_in, obs_type)
  class(load_foram) = "load_foram"
  return(load_foram)
}
