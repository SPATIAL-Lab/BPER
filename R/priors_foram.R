#' Loads prior distributions for foram model
#'
#' This function generates a list of objects associated with the prior distributions to be used in the mcmc
#' inversion. It requires users to update a sheet defining all of the prior distributions, 'priors_in_foram.R'.
#' It determines d11B vital effects slope and intercept for each modern species by running an mcmc inversion of each
#' modern species' calibration data which is routinely updated. The function also generates the temp-, sal-, and press-
#' dependent equilibrium constant look-up arrays which are used by the proxy system model. These arrays are updated to
#' reflect sensitivities to major ion concentrations of Ca, Mg and SO4 in the proxy system model. Dependencies: 
#' 'rjags' and 'R2jags'.
#'
#' @param parms_foram_adj Specify the R object or character string file path to configuration script used to define parameters 
#' and priors. Defaults to system data file if no object or path is provided. 
#'
#' @param age_index User inputs the list that is returned by the 'age_index'. Defaults to 'age_index'.
#'
#' @param cc2ndparm.vt Second carbonate chemistry variable type to be used to compute CO2 in addition to pH.
#' Select from 'dic' 'alk' 'co3' hco3' or 'omegac'. Defaults to 'dic'.
#'
#' @param cc2ndparm.pt Second carbonate chemistry prior type select from 't1' for prior defined at time step 1
#' with autocorrelated time series model or 'ts' for time series input by user to define prior at each time step.
#'
#' @param cc2ndparmTS If 'ts' is selected for 'cc2ndparm.pt', user must provide time series R data object
#' consisting of 3 columns of data: age, mean value for 2nd carb chem variable, and sd for 2nd carb chem variable.
#'
#' @returns Returns list 'priors_foram'. This list contains a list 'clean_pri' which contains all of the prior
#' information used in the model inversion, and 'psm_type' which indicates which proxy model script to use.
#'
#' @examples
#' priors_foram(priors = "package/priors_in_foram.R", age_index = age_index, cc2ndparm.vt = 'dic',
#' cc2ndparm.pt = 'ts', cc2ndparmTS = cc2ndparmTS, save.output = FALSE)
#'
#' @export
priors_foram <- function(parms_foram_adj = NULL, 
                         age_index = age_index, 
                         cc2ndparm.vt = 'dic', 
                         cc2ndparm.pt = 'ts', 
                         cc2ndparmTS){

  ###########################################################################################################################
  ### Observational data ###
  
  prox_in_ai <- age_index[[1]]
  ages_prox <- age_index[[2]]
  dt <- age_index[[3]]
  obs_type <- age_index[[4]]
  step_type <- age_index[[5]]
  
  
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
  d11Bfu.data = c(as.numeric(clean.d11B$d11B2s)/2)
  
  if(length(prox_in_ai$MgCa) < 1){
    mgcaf.fs = 0
    warning("No Mg/Ca data are available to clean")
  } else{
    mgcaf.fs = 1
  }
  
  clean.mgca = prox_in_ai[complete.cases(prox_in_ai$MgCa), ]
  ai.mgca = as.integer(c(clean.mgca$ai))
  mgcaf.data = c(clean.mgca$MgCa)
  mgcafu.data = c(as.numeric(clean.mgca$MgCa2s)/2)
  
  if(length(prox_in_ai$d18O) < 1){
    d18Of.fs = 0
    warning("No d18O data are available to clean")
  } else{
    d18Of.fs = 1
  }
  
  clean.d18O = prox_in_ai[complete.cases(prox_in_ai$d18O), ]
  ai.d18O <- as.integer(c(clean.d18O$ai))
  d18Of.data = c(clean.d18O$d18O)
  d18Ofu.data = c(as.numeric(clean.d18O$d18O2s)/2)
  
  ai.all = c(ai.d11B, ai.mgca, ai.d18O)
  ai.prox = unique(ai.all)
  ai.prox = sort(ai.prox, decreasing = FALSE)
  
  # ai.d11B = match(ai.d11B, ai.prox)
  # ai.mgca = match(ai.mgca, ai.prox)
  # ai.d18O = match(ai.d18O, ai.prox)
  

  ai.prox = seq(from = 1, to = length(dt) + 1, by = 1)
  
  
  clean_obs = list("ai.d11B" = ai.d11B, "si.d11B" = si.d11B, "d11Bf.data" = d11Bf.data, "d11Bfu.data" = d11Bfu.data, "d11Bf.fs" = d11Bf.fs,
                   "ai.mgca" = ai.mgca, "mgcaf.data" = mgcaf.data, "mgcafu.data" = mgcafu.data, "mgcaf.fs" = mgcaf.fs,
                   "ai.d18O" = ai.d18O, "d18Of.data" = d18Of.data, "d18Ofu.data" = d18Ofu.data, "d18Of.fs" = d18Of.fs,
                   "ai.prox" = ai.prox, "dt" = dt, "ages.prox" = ages_prox)
  
  
  ###########################################################################################################################
  ###########################################################################################################################
  #### PRIORS ####
  
  # Load the prior values from user or default package data
  if(is.null(parms_foram_adj)){
    parms_foram <- system.file("extdata", "parms_in_foram.R", package = "BPER")
    source(parms_foram)
  } else if(is.character(parms_foram_adj)){
    parms_foram <- parms_foram_adj
    source(parms_foram)
  } else{
    parms_foram <- parms_foram_adj
    list2env(parms_foram,globalenv())
  }
  

  ###########################################################################################################################
  # Determine mean and sd for vital effect slope and intercept for modern species calibrations
  # using Bayesian inversion ('vital_inv.txt' in 'sysdata' folder) evaluated against modern calibration data ('sysdata' folder)

  model.string=
  "model{
  # Likelihood chunk
  for (i in 1:length(d11Bf)){
    d11Bf[i] ~ dnorm(d11Bforam[i], d11Bf.p[i])
    d11Bf.p[i] = 1/(d11Bf.sd[i])^2  
  }
  for (i in 1:length(d11Bb)){
    d11Bb[i] ~ dnorm(d11Bborate[i], d11Bb.p[i])
    d11Bb.p[i] = 1/(d11Bb.sd[i])^2  
  }
  
  # Model
  for (i in 1:length(d11Bf)){
  d11Bforam[i] = d11Bborate[i] * m.mod  + c.mod
  }
  
  # Priors
  for (i in 1:length(d11Bf)){
    d11Bborate[i] ~ dunif(12, 25)
  }
  m.mod ~ dnorm(m, 1/m.sd^2)  
  c.mod ~ dnorm(c, 1/c.sd^2) 

}
  "

  txtPath <- tempfile(fileext = ".txt")
  writeLines(model.string, con = txtPath)
  
  # G. ruber
  Grub_cal <- vitals$Grub_cal 
  lm_Grub_cal <- lm(Grub_cal$d11Bf ~ Grub_cal$d11Bborate, data = Grub_cal)
  sum_Grub_cal <- summary(lm_Grub_cal)
  m <- sum_Grub_cal$coefficients[2,1]
  m.sd <- 0.25
  c <- sum_Grub_cal$coefficients[1,1]
  c.sd <- 5
  data <- list("d11Bf" = Grub_cal$d11Bf, "d11Bf.sd" = Grub_cal$d11Bf.sd, "d11Bb" = Grub_cal$d11Bborate,
               "d11Bb.sd" = Grub_cal$d11Bborate.sd, "m" = m, "m.sd" = m.sd, "c" = c, "c.sd" = c.sd)
  vital_inv_out = R2jags::jags(model.file = txtPath, parameters.to.save = c("m.mod", "c.mod"),
                               data = data, inits = NULL, n.chains = 3, n.iter = 10000,
                               n.burnin = 3000, n.thin = 1)
  vital_sum <- data.frame(vital_inv_out$BUGSoutput$summary)
  m.Grub <- vital_sum["m.mod","mean"]
  m.Grubsd <- vital_sum["m.mod","sd"]
  c.Grub <- vital_sum["c.mod","mean"]
  c.Grubsd <- vital_sum["c.mod","sd"]

  # T. sacculifer
  Tsac_cal <- vitals$Tsac_cal
  lm_Tsac_cal <- lm(Tsac_cal$d11Bf ~ Tsac_cal$d11Bborate, data = Tsac_cal)
  sum_Tsac_cal <- summary(lm_Tsac_cal)
  m <- sum_Tsac_cal$coefficients[2,1]
  m.sd <- 0.25
  c <- sum_Tsac_cal$coefficients[1,1]
  c.sd <- 5
  data <- list("d11Bf" = Tsac_cal$d11Bf, "d11Bf.sd" = Tsac_cal$d11Bf.sd, "d11Bb" = Tsac_cal$d11Bborate,
               "d11Bb.sd" = Tsac_cal$d11Bborate.sd, "m" = m, "m.sd" = m.sd, "c" = c, "c.sd" = c.sd)
  vital_inv_out = R2jags::jags(model.file = txtPath, parameters.to.save = c("m.mod", "c.mod"),
                               data = data, inits = NULL, n.chains = 3, n.iter = 10000,
                               n.burnin = 3000, n.thin = 1)
  vital_sum <- data.frame(vital_inv_out$BUGSoutput$summary)
  m.Tsac <- vital_sum["m.mod","mean"]
  m.Tsacsd <- vital_sum["m.mod","sd"]
  c.Tsac <- vital_sum["c.mod","mean"]
  c.Tsacsd <- vital_sum["c.mod","sd"]

  # O universa
  Ouni_cal <- vitals$Ouni_cal
  lm_Ouni_cal <- lm(Ouni_cal$d11Bf ~ Ouni_cal$d11Bborate, data = Ouni_cal)
  sum_Ouni_cal <- summary(lm_Ouni_cal)
  m <- sum_Ouni_cal$coefficients[2,1]
  m.sd <- 0.25
  c <- sum_Ouni_cal$coefficients[1,1]
  c.sd <- 5
  data <- list("d11Bf" = Ouni_cal$d11Bf, "d11Bf.sd" = Ouni_cal$d11Bf.sd, "d11Bb" = Ouni_cal$d11Bborate,
               "d11Bb.sd" = Ouni_cal$d11Bborate.sd, "m" = m, "m.sd" = m.sd, "c" = c, "c.sd" = c.sd)
  vital_inv_out = R2jags::jags(model.file = txtPath, parameters.to.save = c("m.mod", "c.mod"),
                               data = data, inits = NULL, n.chains = 3, n.iter = 10000,
                               n.burnin = 3000, n.thin = 1)
  vital_sum <- data.frame(vital_inv_out$BUGSoutput$summary)
  m.Ouni <- vital_sum["m.mod","mean"]
  m.Ounisd <- vital_sum["m.mod","sd"]
  c.Ouni <- vital_sum["c.mod","mean"]
  c.Ounisd <- vital_sum["c.mod","sd"]

  # Borate "vitals" defaults to m=1, c=0 with very low uncertainty (1e-3) which can be adjusted here
  m.borsd = 1e-3       # prescribed uncertainty in "borate vitals" m=1
  c.borsd = 1e-3       # prescribed uncertainty in "borate vitals" c=0

  m.mean <- c(m.Grub, m.Tsac, m.Ouni, m.custom, 1)
  m.sd <- c(m.Grubsd, m.Tsacsd, m.Ounisd, m.customsd, m.borsd)

  c.mean <- c(c.Grub+Grub.coff, c.Tsac+Tsac.coff, c.Ouni+Ouni.coff, c.custom, 0 + borate.coff)
  c.sd <- c(c.Grubsd, c.Tsacsd, c.Ounisd, c.customsd, c.borsd)


  ###########################################################################################################################
  # Generate carb chem equil constant lookup arrays - follows Zeebe & Wolf-Gladrow 2001
  # Set upper and lower STP bounds for equil constant array
  tempC.lb = 0
  tempC.ub = 65
  sal.lb = 15
  sal.ub = 60
  press.lb = 0
  press.ub= 50

  # Step increments for sal (ppt) temp (degrees C) and press (bar)
  t.inc = 0.5 #0.25 min
  s.inc = 0.5 #0.25 min
  p.inc = 1 #1 min

  # Ranges of variables over which to evaluate
  tempC.vr = seq(tempC.lb, tempC.ub, by=t.inc)
  sal.vr = seq(sal.lb, sal.ub, by=s.inc)
  press.vr = seq(press.lb, press.ub, by=p.inc)

  # Initiate arrays
  temp.vr = c(1:length(tempC.vr))
  delV1 = c(1:length(tempC.vr))
  delV2 = c(1:length(tempC.vr))
  delVspc = c(1:length(tempC.vr))
  delVB = c(1:length(tempC.vr))
  delVw = c(1:length(tempC.vr))

  delk1 = c(1:length(tempC.vr))
  delk2 = c(1:length(tempC.vr))
  delkspc = c(1:length(tempC.vr))
  delkB = c(1:length(tempC.vr))
  delkw = c(1:length(tempC.vr))

  base2Darray = c(1:(length(tempC.vr)*length(sal.vr)))
  dim(base2Darray) = c((length(tempC.vr)), (length(sal.vr)))

  Ks1m_st = base2Darray
  Ks2m_st =  base2Darray
  logKsspcm_st = base2Darray
  Ksspcm_st = base2Darray
  lnKsB_st = base2Darray
  KsB_st = base2Darray
  Ksw_st = base2Darray
  K0a = base2Darray

  base3Darray = c(1:(length(tempC.vr)*length(sal.vr)*length(press.vr)))
  dim(base3Darray) = c((length(tempC.vr)), (length(sal.vr)), (length(press.vr)))

  Ks1a = base3Darray
  Ks2a = base3Darray
  Ksspca = base3Darray
  KsBa = base3Darray
  Kswa = base3Darray

  # Constant (cm^3 bar mol^-1 K^-1)
  R <- 83.131

  # Nested for loops to calculate 3D array for K*1a (i.e., array of potential K*1 values w/o any major ion effect)
  for (i in seq_along(tempC.vr)){
    for (j in seq_along(sal.vr)){
      for (k in seq_along(press.vr)){

        temp.vr[i] <- tempC.vr[i]+273.15

        Ks1m_st[i,j] <- exp(2.83655-2307.1266/temp.vr[i]-1.5529413*(log(temp.vr[i]))-((0.20760841+4.0484/temp.vr[i])*sqrt(sal.vr[j]))+0.0846834*sal.vr[j]-0.00654208*(sal.vr[j]^1.5)+log(1-(0.001005*sal.vr[j])))
        Ks2m_st[i,j] <- exp(-9.226508-3351.6106/temp.vr[i]-0.2005743*(log(temp.vr[i]))-((0.106901773+23.9722/temp.vr[i])*sqrt(sal.vr[j]))+0.1130822*sal.vr[j]-0.00846934*(sal.vr[j]^1.5)+log(1-(0.001005*sal.vr[j])))
        logKsspcm_st[i,j] <- ((-171.9065-0.077993*temp.vr[i]+2839.319/temp.vr[i]+71.595*(log(temp.vr[i])/log(10))+(-0.77712+0.0028426*temp.vr[i]+178.34/temp.vr[i])*(sal.vr[j]^0.5)-0.07711*sal.vr[j]+0.0041249*(sal.vr[j]^1.5)))
        Ksspcm_st[i,j] <- 10^(logKsspcm_st[i,j])
        lnKsB_st[i,j] <- ((-8966.9-2890.53*sal.vr[j]^0.5-77.942*sal.vr[j]+1.728*sal.vr[j]^1.5-0.0996*sal.vr[j]^2)/temp.vr[i])+148.0248+137.1942*sal.vr[j]^0.5+1.62142*sal.vr[j]-(24.4344+25.085*sal.vr[j]^0.5+0.2474*sal.vr[j])*(log(temp.vr[i]))+(0.053105*sal.vr[j]^0.5*temp.vr[i])
        KsB_st[i,j] <- exp(lnKsB_st[i,j])
        Ksw_st[i,j] <- exp(148.96502-13847.26/temp.vr[i]-23.6521*(log(temp.vr[i]))+(118.67/temp.vr[i]-5.977+1.0495*(log(temp.vr[i])))*(sal.vr[j]^0.5)-0.01615*sal.vr[j])
        K0a[i,j] <- exp(9345.17/temp.vr[i]-60.2409+23.3585*(log(temp.vr[i]/100))+sal.vr[j]*(0.023517-0.00023656*temp.vr[i]+0.0047036*((temp.vr[i]/100)^2)))

        delV1[i] <- (-25.50)+0.1271*tempC.vr[i]
        delV2[i] <- (-15.82)+(-0.0219*tempC.vr[i])
        delVspc[i]<- (-48.76)+(0.5304*tempC.vr[i])
        delVB[i] <- (-29.48)+0.1622*tempC.vr[i]+(2.608/1000)*tempC.vr[i]^2
        delVw[i] <- (-25.60)+0.2324*tempC.vr[i]+(-3.6246/1000)*tempC.vr[i]^2

        delk1[i] <- (-3.08/1000)+(0.0877/1000)*tempC.vr[i]
        delk2[i] <- (1.13/1000)+(-0.1475/1000)*tempC.vr[i]
        delkspc[i] <- (-11.76/1000)+(0.3692/1000)*tempC.vr[i]
        delkB[i] <- -2.84/1000
        delkw[i] <- (-5.13/1000)+(0.0794/1000)*tempC.vr[i]

        Ks1a[i,j,k] <- (exp(-((delV1[i]/(R*temp.vr[i]))*press.vr[k])+((0.5*delk1[i])/(R*temp.vr[i]))*press.vr[k]^2))*Ks1m_st[i,j]
        Ks2a[i,j,k] <- (exp(-((delV2[i]/(R*temp.vr[i]))*press.vr[k])+((0.5*delk2[i])/(R*temp.vr[i]))*press.vr[k]^2))*Ks2m_st[i,j]
        Ksspca[i,j,k] <- (exp(-((delVspc[i]/(R*temp.vr[i]))*press.vr[k])+((0.5*delkspc[i])/(R*temp.vr[i]))*press.vr[k]^2))*Ksspcm_st[i,j]
        KsBa[i,j,k] <- (exp(-((delVB[i]/(R*temp.vr[i]))*press.vr[k])+((0.5*delkB[i])/(R*temp.vr[i]))*press.vr[k]^2))*KsB_st[i,j]
        Kswa[i,j,k] <- (exp(-((delVw[i]/(R*temp.vr[i]))*press.vr[k])+((0.5*delkw[i])/(R*temp.vr[i]))*press.vr[k]^2))*Ksw_st[i,j]
      }
    }
  }


  ###########################################################################################################################
  # Check that prior type has been correctly specified
  if (cc2ndparm.pt != 'ts' & cc2ndparm.pt != 't1'){
    stop("Must specify 2nd carb chem variable prior type ('cc2ndparm.pt') as either 'ts' for time series or 't1' for initial value in function arguments")
  }

  # Check that variable type has been correctly specified
  if (cc2ndparm.vt != 'dic' & cc2ndparm.vt != 'alk' & cc2ndparm.vt != 'co3' & cc2ndparm.vt != 'hco3' & cc2ndparm.vt != 'omegac'){
    stop("Must specify 2nd carb chem variable type ('cc2ndparm.vt') as either 'dic', 'alk', 'co3' 'omegac' or 'hco3' in function arguments")
  }


  ###########################################################################################################################
  # Generate clean list of priors for t1 prior type
  if (cc2ndparm.pt == 't1'){
    clean_pri = list("tempC.m" = tempC.m, "tempC.sd" = tempC.sd, "sal.m" = sal.m, "sal.sd" = sal.sd, "press.m" = press.m, "press.sd" = press.sd,
                     "d11Bsw.m" = d11Bsw.m, "d11Bsw.sd" = d11Bsw.sd, "d18Osw.m" = d18Osw.m, "d18Osw.sd" = d18Osw.sd,
                     "xca.m" = xca.m, "xca.sd" = xca.sd, "xca.lt" = xca.lt, "xmg.m" = xmg.m, "xmg.sd" = xmg.sd, "xmg.lt" = xmg.lt,
                     "xso4.m" = xso4.m, "xso4.sd" = xso4.sd, "xso4.lt" = xso4.lt, "seccal" = seccal, "seccal.sd" = seccal.sd, "d18Oseccal" = d18Oseccal,
                     "Hp.mean" = Hp.mean, "Hp.sd" = Hp.sd, "Bmod.mean" = Bmod.mean, "Bmod.sd" = Bmod.sd, "A.mean" = A.mean, "A.sd" = A.sd,
                     "pHpccorr" = pHpccorr, "pHpccorrsd" = pHpccorrsd, "m.mean" = m.mean, "m.sd" = m.sd, "c.mean" = c.mean, "c.sd" = c.sd,
                     "pH.l" = pH.l, "pH.u" = pH.u, "carbchem2.m" = carbchem2.m, "carbchem2.sd" = carbchem2.sd,
                     "tempC.lb" = tempC.lb, "sal.lb" = sal.lb, "press.lb" = press.lb, "press.ub" = press.ub,
                     "t.inc" = t.inc, "s.inc" = s.inc, "p.inc" = p.inc,
                     "K0a" = K0a, "Ks1a" = Ks1a, "Ks2a" = Ks2a, "Ksspca" = Ksspca, "KsBa" = KsBa, "Kswa" = Kswa)
  }

  # Generate clean list of priors for time series prior type
  if (cc2ndparm.pt == 'ts'){
    if (length(cc2ndparmTS) < 1){
      stop("If 'ts' is selected, must provide time series R data object in arguments ('cc2ndparmTS') consisting of 3 columns of data: age, mean value for 2nd carb chem varaible, and sd for 2nd carb chem variable")
    } else{
      cc2ndparmTS = data.frame(cc2ndparmTS)
      cc2ndparmTS = cc2ndparmTS[,c(1:3)]
      names(cc2ndparmTS) = c("age.2cc","mean.2cc", "sd.2cc")

      # Linearly interpolate DIC (mean values and 1sd) for ages associated with each time step using input time series
      interp.2cc <- approx(cc2ndparmTS$age.2cc, cc2ndparmTS$mean.2cc, xout=ages_prox, method="linear")
      interp.2cc.sd <- approx(cc2ndparmTS$age.2cc, cc2ndparmTS$sd.2cc, xout=ages_prox, method="linear")
      carbchem2.m <- interp.2cc[["y"]]
      carbchem2.sd <- interp.2cc.sd[["y"]]

      clean_pri = list("tempC.m" = tempC.m, "tempC.sd" = tempC.sd, "sal.m" = sal.m, "sal.sd" = sal.sd, "press.m" = press.m, "press.sd" = press.sd,
                       "d11Bsw.m" = d11Bsw.m, "d11Bsw.sd" = d11Bsw.sd, "d18Osw.m" = d18Osw.m, "d18Osw.sd" = d18Osw.sd,
                       "xca.m" = xca.m, "xca.sd" = xca.sd, "xca.lt" = xca.lt, "xmg.m" = xmg.m, "xmg.sd" = xmg.sd, "xmg.lt" = xmg.lt,
                       "xso4.m" = xso4.m, "xso4.sd" = xso4.sd, "xso4.lt" = xso4.lt, "seccal" = seccal, "seccal.sd" = seccal.sd, "d18Oseccal" = d18Oseccal,
                       "Hp.mean" = Hp.mean, "Hp.sd" = Hp.sd, "Bmod.mean" = Bmod.mean, "Bmod.sd" = Bmod.sd, "A.mean" = A.mean, "A.sd" = A.sd,
                       "pHpccorr" = pHpccorr, "pHpccorrsd" = pHpccorrsd, "m.mean" = m.mean, "m.sd" = m.sd, "c.mean" = c.mean, "c.sd" = c.sd,
                       "pH.l" = pH.l, "pH.u" = pH.u, "carbchem2.m" = carbchem2.m, "carbchem2.sd" = carbchem2.sd,
                       "tempC.lb" = tempC.lb, "sal.lb" = sal.lb, "press.lb" = press.lb, "press.ub" = press.ub,
                       "t.inc" = t.inc, "s.inc" = s.inc, "p.inc" = p.inc,
                       "K0a" = K0a, "Ks1a" = Ks1a, "Ks2a" = Ks2a, "Ksspca" = Ksspca, "KsBa" = KsBa, "Kswa" = Kswa)
    }
  }

  if(clean_pri$pHpccorrsd <= 0 & clean_pri$pHpccorr != 0){
    stop("You have set 'pHpccorrsd' to zero. Only do this if you intend to turn off the pH correction on Mg/Caforam 
            entirely. If this is intended, you must also change 'pHpccorr' to zero")
  }
  
  if(clean_pri$pHpccorrsd <= 0){
    clean_pri$pHpccorrsd <- 1e-20
  }
  
  if(clean_pri$seccal.sd <= 0 & clean_pri$seccal != 0){
    stop("You have set 'seccal.sd' to zero. Only do this if you intend to turn off the diagenetic correction 
            entirely. If this is intended, you must also change 'seccal' to zero")
  }
  
  if(clean_pri$seccal.sd <= 0){
    clean_pri$seccal.sd <- 1e-20
  }
  
  if(clean_pri$tempC.sd <= 0 | clean_pri$sal.sd <= 0 | clean_pri$press.sd <= 0 | clean_pri$d11Bsw.sd <= 0 | clean_pri$d18Osw.sd <= 0 |
     clean_pri$xca.sd <= 0 | clean_pri$xmg.sd <= 0 | clean_pri$xso4.sd <= 0 |clean_pri$Hp.sd <= 0 | clean_pri$Bmod.sd <= 0 |
     clean_pri$A.sd <= 0){
    stop("Standard deviations of user-specified prior distributions must be positive and non-zero, unless turning off a correction")
  }
  
  psm_type <- paste(cc2ndparm.vt, cc2ndparm.pt, sep = "_")
  priors_foram <- list("clean_pri" = clean_pri, "clean_obs" = clean_obs, "psm_type" = psm_type)
  class(priors_foram) = "priors_foram"
  return(priors_foram)
}

