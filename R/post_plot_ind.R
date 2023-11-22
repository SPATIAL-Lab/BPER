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
#' @param show.prior Logical. Specify TRUE if you would like the prior distribution plotted along with the posterior. For 
#' time-dependent variables, the sampled distribution for time step 1 is shown. Defaults to TRUE.
#'
#' @param show.median Logical. Specify TRUE if you would like a line drawn through the median of the distribution. Defaults
#' to FALSE.
#' 
#' @param show.legend Logical. Specify TRUE if you would like a legend in the plot. Defaults to TRUE.
#' 
#' @param leg.pos Option to include a character string to adjust the legend position. Position options are: 'bottomright', 
#' 'bottom', 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right', 'center'. Defaults to "topleft".
#'
#' @returns Returns a plot 'post_plot_out'.
#'
#' @examples
#' post_plot_ind(inv_out = inv_out, parm = "pco2", tstep_age, show.prior = TRUE, show.median = FALSE, show.legend = TRUE)
#'
#' @export
post_plot_ind <- function(inv_out = inv_out, 
                          parm = "pH", 
                          tstep_age, 
                          show.prior = TRUE,
                          show.median = FALSE,
                          show.legend = TRUE,
                          leg.pos = "topleft"){

  if(length(inv_out) != 3){
    stop("'inv_out' must be a list containing 3 elements from 'inv_out' function")
  }

  # load objects from 'inv_out' list
  jout <- inv_out[[1]]
  ages_prox <- inv_out[[2]]
  save.parms <- inv_out[[3]]
  priors <- inv_out[[4]]

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
    units = "(\u00B0C)"
  } else if(parm == "xca" | parm == "xso4" | parm == "xmg"){
    units = "(mmol/kg)"
  } else if(parm == "press"){
    units = "(bar)"
  }

  parm.pri.m = as.character(paste0(parm, ".m"))
  parm.pri.m <- priors[[parm.pri.m]]
  parm.pri.sd = as.character(paste0(parm, ".sd"))
  parm.pri.sd <- priors[[parm.pri.sd]]
  pH.l <- priors$pH.l
  pH.u <- priors$pH.u
  
  if(!is.null(parm.pri.m)){
    xrange <- c((parm.pri.m - (2.5*parm.pri.sd)), (parm.pri.m + (2.5*parm.pri.sd)))
  } else if(parm =="pH"){
    xrange <- c((pH.l - 0.1), (pH.u + 0.1))
  } else{
    xrange = NULL
  }
  
  if(ncol(parm_out) > 1){
    post_plot_ind <- plot(density(parm_out[,match(tstep_age, ages_prox)]), xlim = xrange, main = paste(tstep_age, "ka"), xlab = paste(parm, units), col = "black", lwd=1.5)
    polygon(density(parm_out[,match(tstep_age, ages_prox)]), col = rgb(0,0,0, alpha = 0.4))
    if(isTRUE(show.prior) & !is.null(parm.pri.m) & parm != "pH"){
      polygon(density(rnorm(100000, mean = parm.pri.m, sd = parm.pri.sd)), col = rgb(1,0,0, alpha = 0.3))
    } else if(isTRUE(show.prior) & parm == "pH"){
      polygon(density(runif(100000, min = pH.l, max = pH.u)), col = rgb(1,0,0, alpha = 0.3))
    } else{
      warning("'show.prior' has been turned on and selected parameter is not defined by prior distribution")
    }
    if(isTRUE(show.median)){
      abline(v=median(parm_out[,match(tstep_age, ages_prox)]), col="deepskyblue4", lwd =1.5)
    }
    if(isTRUE(show.legend)){
      if(!is.null(parm.pri.m) | parm == "pH"){
        legend(x = leg.pos, fill = c("grey30", "red"), legend = c("posterior", "prior"))
      } else{
        legend(x = leg.pos, fill = c("grey30"), legend = c("posterior"))
      }
    } 
  } else if(ncol(parm_out) == 1){
    post_plot_ind <- plot(density(parm_out), xlim = xrange, main ="", xlab = paste(parm, units), col = "black", lwd=1.5)
    polygon(density(parm_out), col = rgb(0,0,0, alpha = 0.4))
    if(isTRUE(show.prior) & !is.null(parm.pri.m) & parm != "pH"){
      polygon(density(rnorm(100000, mean = parm.pri.m, sd = parm.pri.sd)))
    } else if(isTRUE(show.prior) & parm == "pH"){
      polygon(density(runif(100000, min = pH.l, max = pH.u)))
    } else{
      warning("'show.prior' has been turned on and selected parameter is not defined by prior distribution")
    }
    if(isTRUE(show.median)){
      abline(v=median(parm_out), col="deepskyblue4", lwd =1.5)
    }
    if(isTRUE(show.legend)){
      if(!is.null(parm.pri.m) | parm == "pH"){
        legend(x = leg.pos, fill = c("grey30", "red"), legend = c("posterior", "prior"))
      } else{
        legend(x = leg.pos, fill = c("grey30"), legend = c("posterior"))
      }
    }
  }

  return(post_plot_ind)
}

