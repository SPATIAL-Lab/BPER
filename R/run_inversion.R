#' Runs the MCMC inversion using groomed data
#'
#' This function loads data to send to jags, writes the proxy system model text string to be used in the inversion
#' and runs the jags MCMC inversion. Dependencies: 'rjags' and 'R2jags'.
#'
#' @param priors User provides the list that is returned by 'priors_foram', 'priors_phyto', 'priors_paleosol', or
#' 'priors_plants' function. Defaults to 'priors_foram'.
#'
#' @param save.parms Provide a vector of character strings indicating each model parameter for which posteriors will be
#' recorded in the output. Defaults to c("sal", "tempC", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", "d18Osw.local",
#' "pco2", "dic", "pH", "press").
#'
#' @param n.iter Provide the number of iterations to perform in the MCMC inversion. See 'jags' user manual for details.
#' Defaults to '10000'.
#'
#' @param n.chains Provide the number of chains for the MCMC inversion. See 'jags' user manual for details. Defaults
#' to '3'.
#'
#' @param n.burnin Provide the number of iterations for inversion burn-in. See 'jags' user manual for details. Defaults
#' to '3000'.
#'
#' @param n.thin Provide the factor by which to thin the output data, i.e., '10' only saves every 10th iteration. See
#' 'jags' user manual for details. Defaults to '1'.
#'
#' @param parallel Logical. Specify TRUE if you want to enable parallel processing. Defaults to FALSE.
#'
#' @returns Returns list 'inv_out' with 'jout' MCMC data object, a vector of timestep ages 'ages_prox', and a list of 
#' the saved parameters 'save.parms'. 'jout' is the same output list that is generated when you call the 'jags' function.
#'
#' @examples
#' run_inversion(priors = priors_foram,
#' save.parms = c("sal", "tempC", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", "d18Osw.local","pco2", "dic", "pH", "press"),
#' n.iter = 10000, n.chains = 3, n.burnin = 3000, n.thin = 1, parallel = FALSE)
#'
#' @export
run_inversion <- function(priors = priors_foram,
                          save.parms = c("sal", "tempC", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", "d18Osw.local","pco2", "dic", "pH", "press"),
                          n.iter = 10000, 
                          n.chains = 3, 
                          n.burnin = 3000, 
                          n.thin = 1, 
                          parallel = FALSE){

  clean_pri <- priors_foram[[1]]
  clean_obs <- priors_foram[[2]]
  psm_type <- priors_foram[[3]]
  data <- append(clean_obs, clean_pri)

  # Determine forward PSM model file to use based on type of 2nd carb chem variable and prior
  # (time series or continuous time autoregressive)

  chunk1 <- chunks$chunk1
  chunk2_dic <- chunks$chunk2_dic 
  chunk2_alk <- chunks$chunk2_alk 
  chunk2_co3 <- chunks$chunk2_co3 
  chunk2_hco3 <- chunks$chunk2_hco3 
  chunk2_omegac <- chunks$chunk2_omegac
  chunk3_t1 <- chunks$chunk3_t1
  chunk3_ts <- chunks$chunk3_ts 

  if(psm_type == 'dic_t1'){
    model.string = paste(chunk1, chunk2_dic, chunk3_t1)
  } else if(psm_type == 'dic_ts'){
    model.string = paste(chunk1, chunk2_dic, chunk3_ts)
  } else if(psm_type == 'alk_t1'){
    model.string = paste(chunk1, chunk2_alk, chunk3_t1)
  } else if(psm_type == 'alk_ts'){
    model.string = paste(chunk1, chunk2_alk, chunk3_ts)
  } else if(psm_type == 'co3_t1'){
    model.string = paste(chunk1, chunk2_co3, chunk3_t1)
  } else if(psm_type == 'co3_ts'){
    model.string = paste(chunk1, chunk2_co3, chunk3_ts)
  } else if(psm_type == 'hco3_t1'){
    model.string = paste(chunk1, chunk2_hco3, chunk3_t1)
  } else if(psm_type == 'hco3_ts'){
    model.string = paste(chunk1, chunk2_hco3, chunk3_ts)
  } else if(psm_type == 'omegac_t1'){
    model.string = paste(chunk1, chunk2_omegac, chunk3_t1)
  } else if(psm_type == 'omegac_ts'){
    model.string = paste(chunk1, chunk2_omegac, chunk3_ts)
  } else{
    stop("Type of PSM not specified. 'psm_type' should be assigned first using 'priors_...' function")
  }
  
  # Save BUGS model string to temporary directory 
  txtPath <- tempfile(fileext = ".txt")
  writeLines(model.string, con = txtPath)
  
  # Call the jags function from R2jags and run the inversion; store output in 'inv_out' data object
  if(isTRUE(parallel)){
    jout = R2jags::jags.parallel(model.file = txtPath, parameters.to.save = save.parms,
                        data = data, inits = NULL, n.chains = n.chains, n.iter = n.iter,
                        n.burnin = n.burnin, n.thin = n.thin)
  } else{
    jout = R2jags::jags(model.file = txtPath, parameters.to.save = save.parms,
                       data = data, inits = NULL, n.chains = n.chains, n.iter = n.iter,
                       n.burnin = n.burnin, n.thin = n.thin)
  }

  ages_prox <- clean_obs[["ages.prox"]]
  inv_out = list("jout" = jout, "ages_prox" = ages_prox, "save.parms" = save.parms, "priors" = clean.pri)
  
  summarydf <- data.frame(jout$BUGSoutput$summary)
  
  if(any(summarydf$n.eff < 50)){
    warning("Some parameters have n.eff statistical parameter values less than 50, consider running more iterations")
  }
  
  if(any(summarydf$Rhat > 1.03)){
    warning("Some parameters have Rhat statistical parameter values greater than 1.03, consider running more iterations")
  }
  
  class(inv_out) = "inv_out"
  return(inv_out)
}

