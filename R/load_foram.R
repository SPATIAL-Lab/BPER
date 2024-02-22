#' Checks and formats foraminiferal proxy data
#'
#' This function allows you to load, QC and format foram data to use in subsequent 'BPER' functions.
#' Function checks which data are present and whether measurement uncertainties are included.
#' If ages or uncertainties are not included the function stops. Function flags 'incomplete'
#' data object - i.e., missing Mg/Ca, or missing d18O, or missing d11B - but proceeds. Dependencies: 'readxl'.
#'
#' @param foram_data User inputs an 8 column R data object and assigns it to 'foram_data' in arguments. Object columns
#' should be organized as: age, d11B, d11B2s, d11Bspec, Mg/Ca, Mg/Ca2s, d18O, d18O2s. The order matters. Current 
#' d11Bspec options are: 'Grub' (G. ruber), 'Tsac' (T. sacculifer), 'Ouni' (O. universa), 'custom' (specify values 
#' with 'priors_foram_adj' function), or 'borate'. If 'csv' or 'xlsx' is to be used instead of R object, provide a 
#' character string path to file from directory. 
#' 
#' @param sheet If '.xlsx' is being used provide a character string of the sheet name to read in. 
#'
#' @returns Returns list 'load_foram'. This list contains 'prox_in' data.frame which has been checked for data completeness
#' and compatibility with other package functions and 'obs_type' which indicates which measurements are included.
#'
#' @examples
#' load_foram(foram_data)
#' 
#' @export
load_foram <- function(foram_data, 
                       sheet = NULL){
  
  if(is.character(foram_data)){
    if(grepl(".csv", foram_data, fixed = TRUE)){
      foram_data = read.csv(foram_data)
    } else if(grepl(".xlsx", foram_data, fixed = TRUE)){
      foram_data = readxl::read_xlsx(path = file.path, sheet = sheet)
    } else{
      stop("Can only specify a file path to .csv or .xlsx file")
    }
  } else{
    foram_data = foram_data
  }
  
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

  
  load_foram = list("prox_in" = prox_in, "obs_type" = obs_type)
  class(load_foram) = "load_foram"
  return(load_foram)
}
