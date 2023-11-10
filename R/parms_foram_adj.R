#' Revises prior distributions
#'
#' This function allows users to pull in the default system file which defines all user-accessible parameters, including priors,
#' and revise it as they see fit. This is a highly recommended step as prior distributions should reflect prior knowledge
#' for the specific study interval. 
#' 
#' @details The following is a list of the revisable parameters which define the prior distributions:
#' --- S, T, P ---
#' 'tempC.m' and 'tempC.sd' are the mean and stdev for temperature in degrees C 
#' 'sal.m' and 'sal.sd' are the mean and stdev for salinity 
#' 'press.m' and 'press.sd' are the mean and stdev for presure in bar
#' --- seawater composition ---
#' 'd11Bsw.m' and 'd11Bsw.sd' are the mean and stdev for seawater d11B in per mille, SRM-951 
#' 'd18Osw.m' and 'd18Osw.sd' are the mean and stdev for seawater d18O in per mille, VSMOW
#' 'xca.m' and 'xca.sd' are the mean and stdev for Ca concentration of seawater in mmol/kg 
#' 'xmg.m' and 'xmg.sd' are the mean and stdev for Mg concentration of seawater in mmol/kg 
#' 'xso4.m' and 'xso4.sd' are the mean and stdev for SO4 concentration of seawater in mmol/kg 
#' 'xca.lt' provides an option to prescribe linear change in Ca concentration of seawater as function of age, mmol/kg per kyr
#' 'xmg.lt' provides an option to prescribe linear change in Mg concentration of seawater as function of age, mmol/kg per kyr
#' 'xso4.lt' provides an option to prescribe linear change in SO4 concentration of seawater as function of age, mmol/kg per kyr
#' --- diagenesis d18O correction ---
#' 'seccal' and 'seccal.sd' are the mean and stdev for the percentage of secondary calcite
#' 'd18Oseccal' is the estiamted d18O, per mille VPDB, of secondary calcite
#' --- Mg/Ca temp calibration parameters ---
#' 'Hp.mean' and 'Hp.sd' are the mean and stdev for nonlinearity of the relationship b/w shell and Mg/Casw
#' 'Bmod.mean' and 'Bmod.sd' are the mean and stdev for modern, pre-corrected, pre-exponential constant in Mg/Ca-SST calibration
#' 'A.mean' and 'A.sd' are the mean and stdev for the exponential constant in Mg/Ca-SST calibration
#' 'pHpccorr' and 'pHpccorrsd' are the mean and stdev for the pH correction on Mg/Caf in % per tenth pH unit
#' --- d11B vital effect ---
#' 'm.custom', 'm.customsd', 'c.custom', and 'c.customsd' specify 'custom' vital effect slope and intercept 
#' 'Grub.coff', 'Tsac.coff', 'Ouni.coff', and 'borate.coff' are the 'c' intercept offsets for each modern species 'c' value; leave '0' for modern 
#' --- carboante chemistry ---
#' 'pH.u' and 'pH.l' are the upper and lower bounds on uniform distribution for time step 1 in total scale
#' 'carbchem2.m' and 'carbchem2.sd' are the mean and stdev for 2nd carbonate chemistry variable for time step 1*. 
#' Carbonate chemistry variables in the following units:
#' pH = 'total scale' equivalent 
#' DIC = μmol/kg
#' ALK = μmol/kg
#' CO3 = μmol/kg
#' HCO3 = μmol/kg
#' *Note that variable type, 'cc2ndparm.vt', should be specified in arguments in 'foram_priors' function. These values will be used if 
#' prior type, i.e. 'cc2ndparm.pt', is set to 't1' in 'foram_priors' argument. These values will not be used if prior type,
#' i.e. 'cc2ndparm.pt', is set to 'ts' in 'foram_priors' argument.
#'
#' @param parms2change Provide a vector of character strings specifying which parameters to revise. The revisable parameters are 
#' defined in 'Details'. 
#' 
#' @param change_values Provide a vector of values which correspond to 'priors2change' vector of parameter names. 
#'
#' @returns Returns list 'parms_foram_adj' which contains all of the values, revised or otherwise, for revisable parameters.
#'
#' @examples
#' parms_foram_adj(parms2change = parms2change, change_values = change_values)
#'
#' @export
parms_foram_adj <- function(parms2change = parms2change, change_values = change_values){
  
  if(length(parms2change) != length(change_values)){
    stop("Length of 'parms2change' vector much match that of 'change_values' vector")
  }
  
  parms_foram <- system.file("extdata", "parms_in_foram.R", package = "BPER")
  source(parms_foram)
  
  for(i in seq_along(parms2change)){
      assign(noquote(parms2change[i]), change_values[i])
  }
  
  parms_foram_adj <- list('tempC.m'=tempC.m, 'tempC.sd'=tempC.sd, 'sal.m'=sal.m, 'sal.sd'=sal.sd, 'press.m'=press.m, 'press.sd'=press.sd, 
                          'd11Bsw.m'=d11Bsw.m, 'd11Bsw.sd'=d11Bsw.sd, 'd18Osw.m'=d18Osw.m, 'd18Osw.sd'=d18Osw.sd, 'xca.m'=xca.m, 
                          'xca.sd'=xca.sd, 'xmg.m'=xmg.m, 'xmg.sd'=xmg.sd, 'xso4.m'=xso4.m, 'xso4.sd'=xso4.sd, 'xca.lt'=xca.lt,
                          'xmg.lt'=xmg.lt, 'xso4.lt'=xso4.lt, 'seccal'=seccal, 'seccal.sd'=seccal.sd, 'd18Oseccal'=d18Oseccal, 
                          'Hp.mean'=Hp.mean, 'Hp.sd'=Hp.sd, 'Bmod.mean'=Bmod.mean, 'Bmod.sd'=Bmod.sd, 'A.mean'=A.mean, 'A.sd'=A.sd, 
                          'pHpccorr'=pHpccorr, 'pHpccorrsd'=pHpccorrsd, 'm.custom'=m.custom, 'm.customsd'=m.customsd, 'c.custom'=c.custom, 
                          'c.customsd'=c.customsd, 'Grub.coff'=Grub.coff, 'Tsac.coff'=Tsac.coff, 'Ouni.coff'=Ouni.coff, 'borate.coff'=borate.coff, 
                          'pH.u'=pH.u, 'pH.l'=pH.l, 'carbchem2.m'=carbchem2.m, 'carbchem2.sd'=carbchem2.sd)
  
  class(parms_foram_adj) = "parms_foram_adj"
  return(parms_foram_adj)
}