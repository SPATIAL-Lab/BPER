#' Formats foraminiferal observational data
#'
#' This function generates a list of observations, age indexes, and species indexes formatted for use in the model
#' inversion. It also species-indexes the d11B data which will allow for different 'vital effects' to be applied
#' to each sample data set if indicated in 'load_foram'.
#'
#' @param age_index User inputs the list that is returned by the 'age_index' function. Defaults to 'age_index'.

#' @returns Returns list 'clean_obs'. This list contains objects associated with observational data, age
#' indexing and species indexing which will be used in the model inversion.
#'
#' @examples
#' clean_foram(age_index = age_index)
#'
#' @export
clean_foram <- function(age_index = age_index){

  prox_in_ai <- age_index[[1]]
  ages_prox <- age_index[[2]]
  dt <- age_index[[3]]
  obs_type <- age_index[[4]]


  if(!inherits(dt, "numeric") | length(dt) < 2){
    stop("Must include numeric 'dt' vector, greater than length 2, in list object in arguments")
  }

  if(length(prox_in_ai$d11B) < 1){
    d11Bf.fs = 0
    warning("No d11B data are available to clean")
  } else{
    d11Bf.fs = 1
  }

  clean.d11B = prox_in_ai[complete.cases(prox_in_ai$d11B), ]
  clean.d11B = transform(clean.d11B, si=ifelse(clean.d11B$d11Bspec=="Grub",1,
                                               ifelse(clean.d11B$d11Bspec=="Tsac",2,
                                                      ifelse(clean.d11B$d11Bspec=="Ouni",3,
                                                             ifelse(clean.d11B$d11Bspec=="custom",4,5)))))
  ai.d11B = as.integer(c(clean.d11B$ai))
  si.d11B = as.integer(c(clean.d11B$si))
  d11Bf.data = c(clean.d11B$d11B)
  d11Bfu.data = c(clean.d11B$d11B2s/2)

  if(length(prox_in_ai$MgCa) < 1){
    mgcaf.fs = 0
    warning("No Mg/Ca data are available to clean")
  } else{
    mgcaf.fs = 1
  }

  clean.mgca = prox_in_ai[complete.cases(prox_in_ai$MgCa), ]
  ai.mgca = as.integer(c(clean.mgca$ai))
  mgcaf.data = c(clean.mgca$MgCa)
  mgcafu.data = c(clean.mgca$MgCa2s/2)

  if(length(prox_in_ai$d18O) < 1){
    d18Of.fs = 0
    warning("No d18O data are available to clean")
  } else{
    d18Of.fs = 1
  }

  clean.d18O = prox_in_ai[complete.cases(prox_in_ai$d18O), ]
  ai.d18O <- as.integer(c(clean.d18O$ai))
  d18Of.data = c(clean.d18O$d18O)
  d18Ofu.data = c(clean.d18O$d18O2s/2)

  ai.all = c(ai.d11B, ai.mgca, ai.d18O)
  ai.prox = unique(ai.all)
  ai.prox = sort(ai.prox, decreasing = FALSE)

  clean_obs = list("ai.d11B" = ai.d11B, "si.d11B" = si.d11B, "d11Bf.data" = d11Bf.data, "d11Bfu.data" = d11Bfu.data, "d11Bf.fs" = d11Bf.fs,
                   "ai.mgca" = ai.mgca, "mgcaf.data" = mgcaf.data, "mgcafu.data" = mgcafu.data, "mgcaf.fs" = mgcaf.fs,
                   "ai.d18O" = ai.d18O, "d18Of.data" = d18Of.data, "d18Ofu.data" = d18Ofu.data, "d18Of.fs" = d18Of.fs,
                   "ai.prox" = ai.prox, "dt" = dt, "ages.prox" = ages_prox, "obs_type" = obs_type)
  
  class(clean_obs) = "clean_obs"
  return(clean_obs)
}
