#' Marine carbonate chemistry calculator
#'
#' This function computes any 6 of the master marine carbonate chemistry variables, plus omega for calcite, from pH data
#' plus an additional carbonate chem variable data set. The user must suply equal length vectors for each influencing variable
#' in the calculations: 'pH', 2nd carb chem variable 'ccparm2.vec', 'temp', 'sal', 'press', 'xca', 'xmg', and 'xso4'.
#' Units for each variable as follows: pH = total scale; ccparm2 = mmol/kg (dic, alk, co3, hco3) OR ppmv (pco2) OR unitless (omegac);
#' temp = deg C, sal = unitless; press = bar; elemental concentrations are in mmol/kg.
#'
#' @param ccparm2 User indicates which 2nd carbonate chemistry variable to use: 'dic', 'alk', 'pco2', 'co3', 'hco3', or 'omegac'.
#' Defaults to 'dic'.
#'
#' @param ccout Specify which parameter to compute: 'dic', 'alk', 'pco2', 'co3', 'hco3', or 'omegac'. Defaults to 'pco2'.
#' 
#' @param inv_out Use inversion output returned object from 'run_inversion' function to compute carbonate chemistry. Using 'inv_out' 
#' alone for this calculation requires that all the variables that go into the carbonate chemistry calculation were included in the 
#' 'save.parms' argument of the 'run_inversion' function. Otherwise, specify equal length vectors of values using individual variable 
#' arguments. Defaults to NULL.
#'
#' @param pH Provide a vector of pH data in total scale if 'inv_out' is not used. Defaults to NULL.
#'
#' @param ccparm2.vec Provide a vector of data of the type specified by 'ccparm2' if 'inv_out' is not used. Defaults to NULL.
#'
#' @param temp Provide a vector of temperature data in degrees C if 'inv_out' is not used. Defaults to NULL.
#'
#' @param sal Provide a vector of salinity data if 'inv_out' is not used. Defaults to NULL.
#'
#' @param press Provide a vector of pressure data in bar if 'inv_out' is not used. Defaults to NULL.
#'
#' @param xca Provide a vector of Ca ion concentration data in mmol/kg if 'inv_out' is not used. Defaults to NULL.
#'
#' @param xmg Provide a vector of Mg ion concentration data in mmol/kg if 'inv_out' is not used. Defaults to NULL.
#'
#' @param xso4 Provide a vector of SO4 ion concentration data in mmol/kg if 'inv_out' is not used. Defaults to NULL.
#'
#' @returns Returns vector of values for parameter of user interest specified by 'ccout' argument.
#'
#' @examples
#' carbchem(ccparm2 = "dic", ccout = "alk", inv_out = inv_out)
#'
#' @export
carbchem <- function(ccparm2 = "dic", 
                     ccout = "alk", 
                     inv_out,
                     pH, 
                     ccparm2.vec,
                     temp, 
                     sal, 
                     press,
                     xca, 
                     xmg, 
                     xso4){
  
  if(!(is.null(inv_out))){
    io_med <- inv_out$jout$BUGSoutput$median
    pH <- as.numeric(io_med$pH)
    ccparm2.vec <- as.numeric(io_med[[ccparm2]])
    temp <- as.numeric(io_med$tempC)
    sal <- as.numeric(io_med$sal)
    press <- rep_len(as.numeric(io_med$press), length.out = as.integer(length(sal)))
    xca <- as.numeric(io_med$xca)
    xmg <- as.numeric(io_med$xmg)
    xso4 <- as.numeric(io_med$xso4)
  }
  

  if(length(pH) != length(ccparm2.vec) | length(pH) != length(temp) | length(pH) != length(sal) | length(pH) != length(press) |
     length(pH) != length(xca) | length(pH) != length(xmg) | length(pH) != length(xso4)){
    stop("Supplied vectors must be equal length")
  }

  ccparm1 <- "pH"
  cc_in <- data.frame(pH, ccparm2.vec, temp, sal, press, xca, xmg, xso4)
  colnames(cc_in) <- c("pH","ccparm2.vec", "temp", "sal", "press", "xca", "xmg", "xso4")



  #################################################################################################
  # Initialize vectors in for loops
  tempC = c()
  temp = c()
  sal = c()
  BT = c()
  Ks1m_st = c()
  Ks2m_st = c()
  logKsspcm_st = c()
  Ksspcm_st = c()
  lnKsB_st = c()
  KsB_st = c()
  Ksw_st = c()
  K0 = c()
  delV1 = c()
  delV2 = c()
  delVspc = c()
  delVB = c()
  delVw = c()
  delk1 = c()
  delk2 = c()
  delkspc = c()
  delkB = c()
  delkw = c()
  press = c()
  Ks1m = c()
  Ks2m = c()
  Ksspcm = c()
  KsB = c()
  Ksw = c()
  xca = c()
  xmg = c()
  xso4 = c()
  Ks1 = c()
  Ks2 = c()
  Ksspc = c()
  pH = c()
  hyd = c()
  pco2 = c()
  fco2 = c()
  co2 = c()
  dic = c()
  co3 = c()
  hco3 = c()
  alk = c()
  omegac = c()


  #################################################################################################
  # Set modern concentrations for Mg, Ca, and SO4
  xcam <- 10.2821*1e-3 # modern [Ca] (mol kg^-1)
  xmgm <- 52.8171*1e-3 # modern [Mg] (mol kg^-1)
  xso4m <- 28.24*1e-3  # modern [SO4] (mol kg^-1)
  mgcaswm <- xmgm/xcam # modern Mg/Ca of seawater

  # Define ZT19 table 2 sensitivity parameters (si_j) for changing sw chemistry on carb chem
  s1_ca <- 5/1000
  s1_mg <- 17/1000
  s1_so4 <- 208/1000
  s2_ca <- 157/1000
  s2_mg <- 420/1000
  s2_so4 <- 176/1000
  sspc_ca <- 185/1000
  sspc_mg <- 518/1000
  sspc_so4 <- 106/1000

  R <- 83.131 # constant (cm^3 bar mol^-1 K^-1)


  #################################################################################################
  # compute equilibrium constants - follows Zeebe and Wolf-Gladrow (2001) unless otherwise noted
  for (i in seq_along(cc_in[,1])){
    # Calculate equil. constants using salinity and temp:
    tempC[i] <- cc_in$temp[i]
    temp[i] <- tempC[i]+273.15
    sal[i] <- cc_in$sal[i]
    BT[i] <- 4.3261e-4 * (sal[i]/35) # Lee et al. (2010)
    Ks1m_st[i] <-exp(2.83655-2307.1266/temp[i]-1.5529413*(log(temp[i]))-((0.20760841+4.0484/temp[i])*sqrt(sal[i]))+0.0846834*sal[i]-0.00654208*(sal[i]^1.5)+log(1-(0.001005*sal[i])))
    Ks2m_st[i] <- exp(-9.226508-3351.6106/temp[i]-0.2005743*(log(temp[i]))-((0.106901773+23.9722/temp[i])*sqrt(sal[i]))+0.1130822*sal[i]-0.00846934*(sal[i]^1.5)+log(1-(0.001005*sal[i])))
    logKsspcm_st[i] <- ((-171.9065-0.077993*temp[i]+2839.319/temp[i]+71.595*(log(temp[i])/log(10))+(-0.77712+0.0028426*temp[i]+178.34/temp[i])*(sal[i]^0.5)-0.07711*sal[i]+0.0041249*(sal[i]^1.5)))
    Ksspcm_st[i] <- 10^(logKsspcm_st[i])
    lnKsB_st[i] <- ((-8966.9-2890.53*sal[i]^0.5-77.942*sal[i]+1.728*sal[i]^1.5-0.0996*sal[i]^2)/temp[i])+148.0248+137.1942*sal[i]^0.5+1.62142*sal[i]-(24.4344+25.085*sal[i]^0.5+0.2474*sal[i])*(log(temp[i]))+(0.053105*sal[i]^0.5*temp[i])
    KsB_st[i] <- exp(lnKsB_st[i])
    Ksw_st[i] <- exp(148.96502-13847.26/temp[i]-23.6521*(log(temp[i]))+(118.67/temp[i]-5.977+1.0495*(log(temp[i])))*(sal[i]^0.5)-0.01615*sal[i])
    K0[i] <- exp(9345.17/temp[i]-60.2409+23.3585*(log(temp[i]/100))+sal[i]*(0.023517-0.00023656*temp[i]+0.0047036*((temp[i]/100)^2)))

    # Adjust equil. constants for the effect of pressure (Millero 1995):
    delV1[i] <- (-25.50)+0.1271*tempC[i]
    delV2[i] <- (-15.82)+(-0.0219*tempC[i])
    delVspc[i]<- (-48.76)+(0.5304*tempC[i])
    delVB[i] <- (-29.48)+0.1622*tempC[i]+(2.608/1000)*tempC[i]^2
    delVw[i] <- (-25.60)+0.2324*tempC[i]+(-3.6246/1000)*tempC[i]^2

    delk1[i] <- (-3.08/1000)+(0.0877/1000)*tempC[i]
    delk2[i] <- (1.13/1000)+(-0.1475/1000)*tempC[i]
    delkspc[i] <- (-11.76/1000)+(0.3692/1000)*tempC[i]
    delkB[i] <- -2.84/1000
    delkw[i] <- (-5.13/1000)+(0.0794/1000)*tempC[i]

    press[i] <- cc_in$press[i]
    Ks1m[i] <- (exp(-((delV1[i]/(R*temp[i]))*press[i])+((0.5*delk1[i])/(R*temp[i]))*press[i]^2))*Ks1m_st[i]
    Ks2m[i] <- (exp(-((delV2[i]/(R*temp[i]))*press[i])+((0.5*delk2[i])/(R*temp[i]))*press[i]^2))*Ks2m_st[i]
    Ksspcm[i] <- (exp(-((delVspc[i]/(R*temp[i]))*press[i])+((0.5*delkspc[i])/(R*temp[i]))*press[i]^2))*Ksspcm_st[i]
    KsB[i] <- (exp(-((delVB[i]/(R*temp[i]))*press[i])+((0.5*delkB[i])/(R*temp[i]))*press[i]^2))*KsB_st[i]
    Ksw[i] <- (exp(-((delVw[i]/(R*temp[i]))*press[i])+((0.5*delkw[i])/(R*temp[i]))*press[i]^2))*Ksw_st[i]

    # K*1, K*2, and K*spc are corrected for past seawater [Ca], [Mg], and [SO4] following ZT19
    xca[i] <- cc_in$xca[i]*1e-3
    xmg[i] <- cc_in$xmg[i]*1e-3
    xso4[i] <- cc_in$xso4[i]*1e-3
    Ks1[i] <- Ks1m[i]*(1+(s1_ca*(xca[i]/xcam-1)+s1_mg*(xmg[i]/xmgm-1)+s1_so4*(xso4[i]/xso4m-1)))
    Ks2[i] <- Ks2m[i]*(1+(s2_ca*(xca[i]/xcam-1)+s2_mg*(xmg[i]/xmgm-1)+s2_so4*(xso4[i]/xso4m-1)))
    Ksspc[i] <- Ksspcm[i]*(1+(sspc_ca*(xca[i]/xcam-1)+sspc_mg*(xmg[i]/xmgm-1)+sspc_so4*(xso4[i]/xso4m-1)))
  }


  #################################################################################################
  # compute carb chem
  if(ccparm1 == ccparm2){
    stop("Please input two different carb chem parameters for ccparm1 and ccparm2")
  }

  # (1)
  # pH and atmospheric CO2 are provided - Zeebe & Wolf Gladrow (2001) appendix B (1)
  if(ccparm1 == "pco2"| ccparm1 =="pH" & ccparm2 == "pco2"| ccparm2 =="pH"){
    for (i in seq_along(cc_in[,1])){

      pH[i] <- cc_in$pH[i]
      hyd[i] <- 10^(-(pH[i]))
      pco2[i] <- cc_in$ccparm2.vec[i] * 1e-6
      fco2[i] <- pco2[i] * 0.9968
      co2[i] <- fco2[i] * K0[i]

      # compute DIC
      dic[i] <- co2[i] * (1 + (Ks1[i] / hyd[i]) + ((Ks1[i] * Ks2[i]) / (hyd[i]^2)))
      # compute [CO3=]
      co3[i] <- dic[i] / (1 + (hyd[i] / Ks2[i]) + ((hyd[i]^2) / (Ks2[i] * Ks1[i])))
      # compute [HCO3-]
      hco3[i] <- dic[i] / (1 + (hyd[i] / Ks1[i]) + (Ks2[i] / hyd[i]))
      # compute TA
      alk[i] <- (co2[i] * ((Ks1[i] / hyd[i]) + (2 * ((Ks1[i] * Ks2[i]) / (hyd[i]^2))))) + ((BT[i]*KsB[i]) / (KsB[i] + hyd[i])) + (Ksw[i] / hyd[i]) - hyd[i]
      # compute omega calcite
      omegac[i] <- (cc_in$xca[i] * co3[i]) / Ksspc[i]
    }

    # (6)
    # pH and [HCO3-] are provided - Zeebe & Wolf Gladrow (2001) appendix B (6)
  } else if(ccparm1 == "pH"| ccparm1 =="hco3" & ccparm2 == "pH"| ccparm2 =="hco3"){

    for (i in seq_along(cc_in[,1])){

      pH[i] <- cc_in$pH[i]
      hyd[i] <- 10^(-(pH[i]))
      hco3[i] <- cc_in$ccparm2.vec[i] * 1e-3

      # compute DIC
      dic[i] <- hco3[i] * (1 + (hyd[i] / Ks1[i]) + ( Ks2[i] / hyd[i]))
      # compute CO2
      co2[i] <- dic[i] / (1 + (Ks1[i] / hyd[i]) + ((Ks1[i] * Ks2[i]) / (hyd[i]^2)))
      fco2[i] <- co2[i] / K0[i]
      pco2[i] <- (fco2[i] / 0.9968) * 1e6
      # compute [CO3=]
      co3[i] <- dic[i] / (1 + (hyd[i] / Ks2[i]) + ((hyd[i]^2) / (Ks1[i] * Ks2[i])))
      # compute TA
      alk[i] <- (co2[i] * ((Ks1[i] / hyd[i]) + (2 * ((Ks1[i] * Ks2[i]) / (hyd[i]^2))))) + ((BT[i]*KsB[i]) / (KsB[i] + hyd[i])) + (Ksw[i] / hyd[i]) - hyd[i]
      # compute omega calcite
      omegac[i] <- (cc_in$xca[i] * co3[i]) / Ksspc[i]
    }

    # (7)
    # pH and [CO3=] are provided - Zeebe & Wolf Gladrow (2001) appendix B (7)
  } else if(ccparm1 == "pH"| ccparm1 =="co3" & ccparm2 == "pH"| ccparm2 =="co3"){

    for (i in seq_along(cc_in[,1])){

      pH[i] <- cc_in$pH[i]
      hyd[i] <- 10^(-(pH[i]))
      co3[i] <- cc_in$ccparm2.vec[i] * 1e-3

      # compute DIC
      dic[i] <- co3[i] * (1 + (hyd[i] / Ks2[i]) + ((hyd[i]^2) / (Ks1[i] * Ks2[i]) ))
      # compute CO2
      co2[i] <- dic[i] / (1 + (Ks1[i] / hyd[i]) + ((Ks1[i] * Ks2[i]) / (hyd[i]^2)))
      fco2[i] <- co2[i] / K0[i]
      pco2[i] <- (fco2[i] / 0.9968) * 1e6
      # compute [HCO3-]
      hco3[i] <- dic[i] / (1 + (hyd[i] / Ks1[i]) + (Ks2[i] / hyd[i]))
      # compute TA
      alk[i] <- (co2[i] * ((Ks1[i] / hyd[i]) + (2 * ((Ks1[i] * Ks2[i]) / (hyd[i]^2))))) + ((BT[i]*KsB[i]) / (KsB[i] + hyd[i])) + (Ksw[i] / hyd[i]) - hyd[i]
      # compute omega calcite
      omegac[i] <- (cc_in$xca[i] * co3[i]) / Ksspc[i]
    }

    # (7 omega)
    # pH and omegac are provided - Zeebe & Wolf Gladrow (2001) appendix B (7)
  } else if(ccparm1 == "pH"| ccparm1 =="omegac" & ccparm2 == "pH"| ccparm2 =="omegac"){

    for (i in seq_along(cc_in[,1])){

      pH[i] <- cc_in$pH[i]
      hyd[i] <- 10^(-(pH[i]))
      omegac[i] <- cc_in$ccparm2.vec[i]

      #compute [CO3=]
      co3[i] <-(omegac[i] * Ksspc[i]) / cc_in$xca[i]
      # compute DIC
      dic[i] <- co3[i] * (1 + (hyd[i] / Ks2[i]) + ((hyd[i]^2) / (Ks1[i] * Ks2[i]) ))
      # compute CO2
      co2[i] <- dic[i] / (1 + (Ks1[i] / hyd[i]) + ((Ks1[i] * Ks2[i]) / (hyd[i]^2)))
      fco2[i] <- co2[i] / K0[i]
      pco2[i] <- (fco2[i] / 0.9968) * 1e6
      # compute [HCO3-]
      hco3[i] <- dic[i] / (1 + (hyd[i] / Ks1[i]) + (Ks2[i] / hyd[i]))
      # compute TA
      alk[i] <- (co2[i] * ((Ks1[i] / hyd[i]) + (2 * ((Ks1[i] * Ks2[i]) / (hyd[i]^2))))) + ((BT[i]*KsB[i]) / (KsB[i] + hyd[i])) + (Ksw[i] / hyd[i]) - hyd[i]
    }

    # (8)
    # pH and TA are provided - Zeebe & Wolf Gladrow (2001) appendix B (8)
  } else if(ccparm1 == "pH"| ccparm1 =="alk" & ccparm2 == "pH"| ccparm2 =="alk"){

    for (i in seq_along(cc_in[,1])){

      pH[i] <- cc_in$pH[i]
      hyd[i] <- 10^(-(pH[i]))
      alk[i] <- cc_in$ccparm2.vec[i] * 1e-3

      # compute CO2
      co2[i] <- (alk[i] - ((KsB[i] * BT[i]) / (KsB[i] + hyd[i])) - (Ksw[i] / hyd[i]) + hyd[i]) / ((Ks1[i] / hyd[i]) + (2*((Ks1[i] * Ks2[i])/(hyd[i]^2))))
      fco2[i] <- co2[i] / K0[i]
      pco2[i] <- (fco2[i] / 0.9968) * 1e6
      # compute DIC
      dic[i] <- co2[i] * (1 + (Ks1[i] / hyd[i]) + ((Ks1[i] * Ks2[i]) / (hyd[i]^2)))
      # compute [HCO3-]
      hco3[i] <- dic[i] / (1 + (hyd[i] / Ks1[i]) + (Ks2[i] / hyd[i]))
      # compute [CO3=]
      co3[i] <- dic[i] / (1 + (hyd[i] / Ks2[i]) + ((hyd[i]^2) / (Ks1[i] * Ks2[i])))
      # compute omega calcite
      omegac[i] <- (cc_in$xca[i] * co3[i]) / Ksspc[i]
    }

    # (9)
    # pH and DIC are provided  - Zeebe & Wolf Gladrow (2001) appendix B (9)
  } else if(ccparm1 == "pH"| ccparm1 =="dic" & ccparm2 == "pH"| ccparm2 =="dic"){

    for (i in seq_along(cc_in[,1])){

      pH[i] <- cc_in$pH[i]
      hyd[i] <- 10^(-(pH[i]))
      dic[i] <- cc_in$ccparm2.vec[i] * 1e-3

      # compute CO2
      co2[i] <- dic[i] / (1 + (Ks1[i] / hyd[i]) + ((Ks1[i] * Ks2[i]) / (hyd[i]^2)))
      fco2[i] <- co2[i] / K0[i]
      pco2[i] <- (fco2[i] / 0.9968) * 1e6
      # compute TA
      alk[i] <- (co2[i] * ((Ks1[i] / hyd[i]) + (2 * ((Ks1[i] * Ks2[i]) / (hyd[i]^2))))) + ((BT[i]*KsB[i]) / (KsB[i] + hyd[i])) + (Ksw[i] / hyd[i]) - hyd[i]
      # compute [HCO3-]
      hco3[i] <- dic[i] / (1 + (hyd[i] / Ks1[i]) + (Ks2[i] / hyd[i]))
      # compute [CO3=]
      co3[i] <- dic[i] / (1 + (hyd[i] / Ks2[i]) + ((hyd[i]^2) / (Ks1[i] * Ks2[i])))
      # compute omega calcite
      omegac[i] <- (cc_in$xca[i] * co3[i]) / Ksspc[i]
    }
  } else{
    stop("Parameter combination not yet built in, please select 'pH' as one of the ccparm inputs")
  }

  if(ccout == "pH"){
    return(pH)
  } else if (ccout == "dic"){
    return(dic)
  } else if (ccout == "alk"){
    return(alk)
  } else if (ccout == "pco2"){
    return(pco2)
  } else if (ccout == "co3"){
    return(co3)
  } else if (ccout == "hco3"){
    return(hco3)
  } else if (ccout == "omegac"){
    return(omegac)
  }
}

