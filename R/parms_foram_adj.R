#' Revises prior distributions and parameters
#'
#' This function allows users to pull in the default system file which defines all user-accessible parameters, including priors,
#' and revise it as they see fit. This is a highly recommended step as prior distributions should reflect prior knowledge
#' for the specific study interval. 
#' 
#' @details The following is a list of the revisable parameters (default values in parentheses) which define the prior distributions: \cr
#' \cr
#' --- S, T, P --- \cr
#' **'tempC.m'** (30) and **'tempC.sd'** (5) are the mean and stdev for temperature in degrees C  \cr
#' **'sal.m'** (35) and **'sal.sd'** (0.5) are the mean and stdev for salinity  \cr
#' **'press.m'** (6) and **'press.sd'** (1) are the mean and stdev for presure in bar \cr
#' \cr
#' --- seawater composition --- \cr
#' **'d11Bsw.m'** (38.45) and **'d11Bsw.sd'** (0.5) are the mean and stdev for seawater d11B in per mille, SRM-951  \cr
#' **'d18Osw.m'** (-1.2) and **'d18Osw.sd'** (0.2) are the mean and stdev for seawater d18O in per mille, VSMOW \cr
#' **'xca.m'** (21.41) and **'xca.sd'** (0.5) are the mean and stdev for Ca concentration of seawater in mmol/kg  \cr
#' **'xmg.m'** (68.51) and **'xmg.sd'** (0.5) are the mean and stdev for Mg concentration of seawater in mmol/kg  \cr
#' **'xso4.m'** (14) and **'xso4.sd'** (0.5) are the mean and stdev for SO4 concentration of seawater in mmol/kg  \cr
#' **'xca.lt'** (-1.9e-4) provides an option to prescribe linear change in Ca concentration of seawater as function of age, mmol/kg per kyr \cr
#' **'xmg.lt'** (0) provides an option to prescribe linear change in Mg concentration of seawater as function of age, mmol/kg per kyr \cr
#' **'xso4.lt'** (0) provides an option to prescribe linear change in SO4 concentration of seawater as function of age, mmol/kg per kyr \cr
#' \cr
#' --- diagenesis d18O correction --- \cr
#' **'seccal'** (0) and **'seccal.sd'** (0) are the mean and stdev for the percentage of secondary calcite \cr
#' **'d18Oseccal'** (0.85) is the estiamted d18O, per mille VPDB, of secondary calcite \cr
#' \cr
#' --- Mg/Ca temp calibration parameters --- \cr
#' **'Hp.mean'** (0.74) and **'Hp.sd'** (0.05) are the mean and stdev for nonlinearity of the relationship b/w shell and Mg/Casw \cr
#' **'Bmod.mean'** (0.38) and **'Bmod.sd'** (0.02) are the mean and stdev for modern, pre-corrected, pre-exponential constant in Mg/Ca-SST calibration \cr
#' **'A.mean'** (0.0757) and **'A.sd'** (0.0045) are the mean and stdev for the exponential constant in Mg/Ca-SST calibration \cr
#' **'pHpccorr'** (0) and **'pHpccorrsd'** (0) are the mean and stdev for the pH correction on Mg/Caf in percent per tenth pH unit \cr
#' \cr
#' --- d11B vital effect --- \cr
#' **'m.custom'** (0.9), **'m.customsd'** (0.2), **'c.custom'** (2), and **'c.customsd'** (2) specify 'custom' vital effect slope and intercept  \cr
#' **'Grub.coff'** (5.76), **'Tsac.coff'** (3.62), **'Ouni.coff'** (0), and **'borate.coff'** (0) are the 'c' intercept offsets for each modern species 'c' value; leave '0' for modern \cr
#' \cr
#' --- carbonate chemistry --- \cr
#' **'pH.u'** (7.75) and **'pH.l'** (7.45) are the upper and lower bounds on uniform distribution for time step 1 in total scale \cr
#' **'carbchem2.m'** (2000) and **'carbchem2.sd'** (150) are the mean and stdev for 2nd carbonate chemistry variable for time step 1*.  \cr
#' \cr
#' Carbonate chemistry variables in the following units: \cr
#' pH = 'total scale' equivalent  \cr
#' DIC = μmol/kg \cr
#' ALK = μmol/kg \cr
#' CO3 = μmol/kg \cr
#' HCO3 = μmol/kg \cr
#' \cr
#' *Note that variable type, 'cc2ndparm.vt', should be specified in arguments in 'foram_priors' function. These values will be used if
#' prior type, i.e. 'cc2ndparm.pt', is set to 't1' in 'foram_priors' argument. These values will not be used if prior type,
#' i.e. 'cc2ndparm.pt', is set to 'ts' in 'foram_priors' argument.
#' @md
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
parms_foram_adj <- function(parms2change = parms2change, 
                            change_values = change_values){
  
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
