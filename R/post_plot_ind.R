#' Plots parameter posterior for a single time step
#'
#' This function generates a time series plot for a mcmc inversion output parameter of interest. This parameter must be
#' included in 'save.parms' argument in 'run_inversion' function. Two plot style options are included.
#'
#' @param inv_out User inputs the MCMC data object that is returned by the 'run_inversion' function. Defaults to 'inv_out'.
#'
#' @param parm Specify which parameter to plot. This parameter must be included in 'save.parms' argument in 'run_inversion'
#' function. Defaults to 'pco2'.
#'
#' @param tstep_age Indicate which time step age is desired if the parameter is time-dependent. This age much match one of
#' the input ages, compiled in 'ages_prox' object. No default.
#'
#' @param show.median Logical. Specify TRUE if you would like a line drawn through the median of the distribution. Defaults
#' to FALSE.
#'
#' @returns Returns a plot 'post_plot_out'.
#'
#' @examples
#' post_plot_ind(inv_out = inv_out, parm = "pco2", tstep_age, show.median = FALSE)
#'
#' @export
post_plot_ind <- function(inv_out = inv_out, parm = "pco2", tstep_age, show.median = FALSE){

  if(length(inv_out) != 3){
    stop("'inv_out' must be a list containing 3 elements from 'inv_out' function")
  }

  # load objects from 'inv_out' list
  jout <- inv_out[[1]]
  ages_prox <- inv_out[[2]]
  save.parms <- inv_out[[3]]

  parm_out <- jout$BUGSoutput$sims.list[[parm]]
  if(ncol(parm_out) > 1 & !(tstep_age %in% ages_prox)){
    stop("Parameter of interest is time-dependent, 'tstep_age' is required. Value must be in kyr and included in vector of time step ages, 'ages_prox'")
  }

  # set units for parameters
  units <- ""
  if(parm == "pco2"){
    units <- "(ppmv)"
  } else if(parm == "dic" | parm == "alk" | parm == "co3" | parm == "hco3"){
    parm_out <- parm_out*1e6
    parm_med <- parm_med*1e6
    units <- "(μmol/kg)"
  } else if(parm == "d11Bsw"){
    units = "(\u2030 SRM-951)"
  } else if(parm == "d18Osw"){
    units = "(\u2030 VSMOW)"
  } else if(parm == "tempC"){
    units = expression(paste("(", degree, "C", ")"))
  } else if(parm == "xca" | parm == "xso4" | parm == "xmg"){
    units = "(mmol/kg)"
  } else if(parm == "press"){
    units = "(bar)"
  }

  if(ncol(parm_out) > 1){
    match(tstep_age, ages_prox)
    post_plot_ind <- plot(density(parm_out[,match(tstep_age, ages_prox)]), main = paste(tstep_age, "ka"), xlab = paste(parm, units), col = "black", lwd=1.5)
    polygon(density(parm_out[,match(tstep_age, ages_prox)]), col = rgb(0,0,0, alpha = 0.4))
    if(show.median == TRUE){
      abline(v=median(parm_out[,match(tstep_age, ages_prox)]), col="deepskyblue4", lwd =1.5)
    }
  } else if(ncol(parm_out) == 1){
    post_plot_ind <- plot(density(parm_out), main ="", xlab = paste(parm, units), col = "black", lwd=1.5)
    polygon(density(parm_out), col = rgb(0,0,0, alpha = 0.4))
    if(show.median == TRUE){
      abline(v=median(parm_out), col="deepskyblue4", lwd =1.5)
    }
  }

  return(post_plot_ind)
}

