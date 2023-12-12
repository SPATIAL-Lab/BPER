#' Plots time series parameter posteriors
#'
#' This function generates a time series plot for a mcmc inversion output parameter. This parameter must be
#' included in 'save.parms' argument in 'run_inversion' function. Two plot style options are included.
#'
#' @param inv_out User inputs the MCMC data object that is returned by the 'run_inversion' function. Defaults to 'inv_out'.
#'
#' @param parm Specify which parameter to plot. This parameter must be included in 'save.parms' argument in 'run_inversion'
#' function. Defaults to 'pco2'.
#'
#' @param type Indicate which plot style desired. Either 'CI' for 95 percent credible interval, or 'draws' for 500 draws of time
#' series posteriors. Defaults to 'CI'.
#' 
#' @param show.legend Logical. Specify TRUE if you would like a legend in the plot. Defaults to TRUE.
#' 
#' @param leg.pos Option to include a character string to adjust the legend position. Position options are: 'bottomright', 
#' 'bottom', 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right', 'center'. Defaults to "topleft".
#'
#' @returns Returns a plot called 'post_plot'.
#'
#' @examples
#' post_plot(inv_out = inv_out, parm = "pco2", type = "CI")
#'
#' @export
post_plot <- function(inv_out = inv_out, 
                      parm = "dic", 
                      type = "CI",
                      show.legend = TRUE,
                      leg.pos = "topleft"){

  if(length(inv_out) != 4){
    stop("'inv_out' must be a list containing 3 elements from 'run_inversion' function")
  }

  # load objects from 'inv_out' list
  jout <- inv_out[[1]]
  ages_prox <- inv_out[[2]]
  save.parms <- inv_out[[3]]
  parm_out <- jout$BUGSoutput$sims.list[[parm]]
  parm_med <- jout$BUGSoutput$median[[parm]]
  units <- ""

  # pull credible intervals for parameters
  step.vector <- seq(1, length(parm_med), by=1)
  parm_sum <- data.frame(jout$BUGSoutput$summary)
  parm.v <- paste(parm, "[", step.vector, "]", sep="")
  parm_sum <- parm_sum[c(parm.v),]
  parm_sum <- data.frame(ages_prox, parm_sum)
  parm_sum <- na.omit(parm_sum)

  # set units for parameters
  if(parm == "pco2"){
    units <- "(ppmv)"
  } else if(parm == "dic" | parm == "alk" | parm == "co3" | parm == "hco3"){
    parm_out <- parm_out*1e6
    parm_med <- parm_med*1e6
    parm_sum <- parm_sum*1e6
    units <- "(μmol/kg)"
  } else if(parm == "d11Bsw"){
    units = "(\u2030 SRM-951)"
  } else if(parm == "d18Osw"){
    units = "(\u2030 VSMOW)"
  } else if(parm == "tempC"){
    units = "(\u00B0C)" 
  } else if(parm == "xca" | parm == "xso4" | parm == "xmg"){
    units = "(mmol/kg)"
  }

  # stop if parm is not time-dependent or it is not in the save.parm list
  if(!(parm %in% save.parms) | length(parm_out) < 2){
    stop("Must input a time series parameter (i.e., 'parm' argument) that is in the 'save.parms' list used as argument for 'inv_out' function")
  }

  # generate time series plot (either random draws or confidence intervals) of parameter of interest
  if(type == "draws"){
    post_plot <- plot(ages_prox, parm_out[as.integer(runif(1,1,nrow(parm_out))),], type="l", xlab = "Age (kyr)", ylab = paste(parm, units), xlim = rev(range(ages_prox)), ylim = range(parm_out), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
    for (i in as.integer(runif(500,1,nrow(parm_out)))){
      lines(ages_prox, parm_out[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
    }
    lines(ages_prox, parm_med, col=rgb(red=1, green=0, blue=0), lwd=1.5)
    if(isTRUE(show.legend)){
      legend(x = leg.pos, fill = c("grey30", "red"), legend = c("sampled posterior draws", "median of posterior distributions"))
    }
  } else if(type =="CI"){
    post_plot <- plot(ages_prox, parm_med, type="l", xlab = "Age (kyr)", ylab = paste(parm, units), xlim = rev(range(ages_prox)), ylim = range(parm_out), col="black")
    polygon(c(ages_prox, rev(ages_prox)), c(parm_sum[,4], rev(parm_sum[,8])), col = "gray", lwd = 0.5)
    lines(ages_prox, parm_med, col="red", lwd=1.5)
    if(isTRUE(show.legend)){
      legend(x = leg.pos, fill = c("gray", "red"), legend = c("95% credible interval", "median of posterior distributions"))
    }
  } else{
    stop("Must specify type of plot to be drawn, either 'CI' for 50% and 95% confidence interval or 'draws' for 500 draws of time series posteriors")
  }
    

  return(post_plot)
}
