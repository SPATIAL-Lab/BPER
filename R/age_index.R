#' Age indices are assigned to the proxy data
#'
#' User specifies time stepping to be used by model, either 'regular' for regular interval time stepping
#' with step interval specified, i.e., 'step_int', OR 'every' for variable time steps at each age for which there is data.
#' Each row of data is assigned an index which will be used in the time series model.
#'
#' @param load_proxy User inputs the list that is returned by the 'load_foram', 'load_phyto', 'load_paleosol', or 'load_plants'
#' functions. Defauts to 'load_foram'.
#' 
#' @param age_units Indicate the age units used in the input data. Either 'kyr' or 'Myr'. Defaults to 'kyr'.
#' 
#' @param step_type Indicate either 'regular' for regular interval time stepping with step interval specified,
#' i.e., 'step_int', OR 'every' for variable time steps at each age for which there is data, and 'step_int' is not required, 
#' OR 'both' to have time steps at each age for which there is data in addition to regular time step intervals - 'step_int' must
#' be specified for 'both' step_type argument) OR 'custom' to provide your own vector of time step ages ('cust_steps' vector 
#' of ages must be included as argument for 'custom' step_type). Defaults to 'regular'.
#' 
#' @param step_int If 'step_type' is 'regular' or 'both', indicated the 'step_int' or regular step interval to be used. 
#' Defaults to '10'.
#'
#' @param cust_steps If 'custom' is specified for 'step_type' provide a vector of ages. There will be a time step for each 
#' age supplied. Defaults to null. 
#'
#' @returns Returns list 'age_index'. This list contains 1. 'prox_in_ai' data frame which has been age indexed by the function,
#' 2. 'ages_prox' which is a vector of the ages of each model time step to be used, 3. 'dt' which is a vector of the delta time, or
#' amount of time between each time step, and 4. 'obs_type' carried over from input list.
#'
#' @examples
#' age_index(load_proxy = load_foram, age_units = "kyr", step_type = "regular", step_int = 10)
#'
#' @export
age_index <- function(load_proxy = load_foram, 
                      age_units = "kyr", 
                      step_type = "regular", 
                      step_int = 10,
                      cust_steps){

  prox_in <- load_proxy[[1]]
  if(!is.null(load_proxy[[2]])){
    obs_type <- load_proxy[[2]]
  }

  if(age_units == "kyr"){
    prox_in$age = prox_in$age
  } else if(age_units == "Myr"){
    prox_in$age = prox_in$age*1e3
    step_int = step_int*1e3
    if(!is.null(cust_steps)){
      cust_steps = cust_steps*1e3
    }
  } else{
    stop("Must specify either 'kyr' or 'Myr' to indicate the input data 'age_units' in function arguments")
  }

  ages_prox = unique(prox_in$age)

  if(step_type == "every"){
    prox_in_ai = transform(prox_in, ai = as.numeric(factor(age*-1)))
    
  } else if(step_type == "regular"){
    if(step_int > 0 & step_int < ((max(ages_prox) - min(ages_prox))/2)){
      prox_in_ai = transform(prox_in, ai = as.numeric(1 + floor((max(prox_in$age) - prox_in$age) / step_int)))
      ages_prox = seq(from = max(ages_prox), to = min(ages_prox), by = -1*step_int) - 0.5*step_int
    } else{
      stop("Must specify positive value less than half the length of the total time interval for 'step_int' in function arguments if regular time step type is selected")
    }
    
  } else if(step_type =="both"){
    if(step_int > 0 & step_int < ((max(ages_prox) - min(ages_prox))/2)){
      ai = c(1:length(prox_in$age))
      ages_prox = unique(sort(append(prox_in$age,(seq(from = max(ages_prox), to = min(ages_prox), by = -1*step_int) - 0.5*step_int)), decreasing = TRUE))
      for (i in seq_along(prox_in$age)){
        ai[i] <- which(abs(ages_prox - prox_in$age[i]) == min(abs(ages_prox - prox_in$age[i])))
      }
      prox_in_ai = transform(prox_in, ai = sort(as.numeric(ai)))
    } else{
      stop("Must specify positive value less than half the length of the total time interval for 'step_int' in function arguments if 'both' time step type is selected")
    } 
    
  } else if (step_type =="custom"){
      if(is.null(cust_steps)){
        stop("Must include a vector of ages for 'cust_steps' argument if 'custom' is specified for 'step_type'.")
      }
    ai = c(1:length(prox_in$age))
    ages_prox = sort(cust_steps, decreasing = TRUE)
    for(i in seq_along(prox_in$age)){
      ai[i] <- which(abs(ages_prox - prox_in$age[i]) == min(abs(ages_prox - prox_in$age[i])))
    }
    prox_in_ai = transform(prox_in, ai = sort(as.numeric(ai)))
  } else{
    stop("Must specify time step type: 'regular', 'every', 'both', or 'custom'")
  }
  
  dt = abs(diff(ages_prox, lag=1))
  
  age_index <- list("prox_in_ai" = prox_in_ai, "ages_prox" = ages_prox, "dt" = dt, "obs_type" = obs_type, "step_type" = step_type)
  class(age_index) = "age_index"
  return(age_index)
}
