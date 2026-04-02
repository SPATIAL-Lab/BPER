model {
  
  # Likelihood block 
  ####################################################################################################

  # Benthic forams
  d13Cbf.analyt.sd <- 0.05
  d13Cbf.p <- 1/d13Cbf.analyt.sd^2
  for (i in 1:n.d13Cbf){
    d13Cbf.data[i] ~ dnorm(d13Cbf[ri.d13Cbf[i]], d13Cbf.p)
  }
  
  # Planktic forams
  d13Cpf.analyt.sd <- 0.05
  d13Cpf.p <- 1/d13Cpf.analyt.sd^2
  for (i in 1:n.d13Cpf){
    d13Cpf.data[i] ~ dnorm(d13Cpf[ri.d13Cpf[i]], d13Cpf.p)
  }
  
  # Brachiopods
  d13Cbrach.analyt.sd <- 0.05
  d13Cbrach.p <- 1/d13Cbrach.analyt.sd^2
  for (i in 1:n.d13Cbrach){
    d13Cbrach.data[i] ~ dnorm(d13Cbrach[ri.d13Cbrach[i]], d13Cbrach.p)
  }
  
  # Bivalves
  d13Cbivalve.analyt.sd <- 0.05
  d13Cbivalve.p <- 1/d13Cbivalve.analyt.sd^2
  for (i in 1:n.d13Cbivalve){
    d13Cbivalve.data[i] ~ dnorm(d13Cbivalve[ri.d13Cbivalve[i]], d13Cbivalve.p)
  }
  
  # Ammonite
  d13Camm.analyt.sd <- 0.05
  d13Camm.p <- 1/d13Camm.analyt.sd^2
  for (i in 1:n.d13Camm){
    d13Camm.data[i] ~ dnorm(d13Camm[ri.d13Camm[i]], d13Camm.p)
  }
  
  # Belemnite
  d13Cbel.analyt.sd <- 0.05
  d13Cbel.p <- 1/d13Cbel.analyt.sd^2
  for (i in 1:n.d13Cbel){
    d13Cbel.data[i] ~ dnorm(d13Cbel[ri.d13Cbel[i]], d13Cbel.p)
  }

  # micrite open ocean
  d13Cmicrite.analyt.sd <- 0.05
  d13Cmicrite.p <- 1/(d13Cmicrite.analyt.sd)^2
  for (i in 1:n.d13Cmicrite){
    d13Cmicrite.data[i] ~ dnorm(d13Cmicrite[ri.d13Cmicrite[i]], d13Cmicrite.p)
  }

  # Bulk carbonate open ocean 
  d13Cbulk.analyt.sd <- 0.05
  d13Cbulk.p <- 1/d13Cbulk.analyt.sd^2
  for (i in 1:n.d13Cbulk){
    d13Cbulk.data[i] ~ dnorm(d13Cbulk[ri.d13Cbulk[i]], d13Cbulk.p)
  }

  # Bulk carbonate semi restricted
  d13Cbulk_sr.analyt.sd <- 0.05
  d13Cbulk_sr.p <- 1/(d13Cbulk_sr.analyt.sd)^2
  for (i in 1:n.d13Cbulk_sr){
    d13Cbulk_sr.data[i] ~ dnorm(d13Cbulk_sr[ri.d13Cbulk_sr[i]], d13Cbulk_sr.p)
  }

  # Bulk carbonate marginal sea
  d13Cbulk_marg.analyt.sd <- 0.05
  d13Cbulk_marg.p <- 1/(d13Cbulk_marg.analyt.sd)^2
  for (i in 1:n.d13Cbulk_marg){
    d13Cbulk_marg.data[i] ~ dnorm(d13Cbulk_marg[ri.d13Cbulk_marg[i]], d13Cbulk_marg.p)
  }


  # Constants
  ####################################################################################################  
  # Equation 8 of Tipple et al. (2010)
  eps.dic_cc <- -1   
  
  # Bulk carbonate uncertainty terms from Romanek et al. (1992) 
  cc_co2_constant1 ~ dnorm(11.98, 1/0.13^2)T(11.72,12.24)
  cc_co2_coeff1 ~ dnorm(0.12, 1/0.01^2)T(0.1,0.14)
  
  # Benthic disequilibrium effects; totals to 2.8 (mean) from Tipple et al. (2010)
  asd ~ dnorm(1, 1/0.2^2)T(0,)  # air sea disequilibrium 
  bpump ~ dnorm(1.2, 1/0.4^2)T(0,) # biological pump
  remin ~ dnorm(0.6, 1/0.3^2)T(0,) # remineralization and oxidation 
  A <- asd + bpump + remin
  f_co3 ~ dnorm(0.12, 1/0.04^2)T(0.04,0.20)
  f_carbacid ~ dnorm(0.01, 1/0.005^2)T(0,0.02)
  
  # Non-secular bias uncertainty - hyperpriors set the across-site scatter (partial pooling)
  bf.nsb_mean ~ dnorm(bf.nsb.m, 1/bf.nsb.sd^2)         
  bf.nsb_tau ~ dgamma(1e3, 1e-3)
  pf.nsb_mean ~ dnorm(pf.nsb.m, 1/pf.nsb.sd^2)         
  pf.nsb_tau ~ dgamma(1e3, 1e-3)
  brach.nsb_mean ~ dnorm(brach.nsb.m, 1/brach.nsb.sd^2)         
  brach.nsb_tau ~ dgamma(1e3, 1e-3)
  bivalve.nsb_mean ~ dnorm(bivalve.nsb.m, 1/bivalve.nsb.sd^2)         
  bivalve.nsb_tau ~ dgamma(1e3, 1e-3)
  amm.nsb_mean ~ dnorm(amm.nsb.m, 1/amm.nsb.sd^2)         
  amm.nsb_tau ~ dgamma(1e3, 1e-3)
  bel.nsb_mean ~ dnorm(bel.nsb.m, 1/bel.nsb.sd^2)         
  bel.nsb_tau ~ dgamma(1e3, 1e-3)
  bulk.nsb_mean ~ dnorm(bulk.nsb.m, 1/bulk.nsb.sd^2)         
  bulk.nsb_tau ~ dgamma(1e3, 1e-3)
  micrite.nsb_mean ~ dnorm(micrite.nsb.m, 1/micrite.nsb.sd^2)         
  micrite.nsb_tau ~ dgamma(1e3, 1e-3)
  bulk_sr.nsb_mean ~ dnorm(bulk_sr.nsb.m, 1/bulk_sr.nsb.sd^2)         
  bulk_sr.nsb_tau ~ dgamma(1e3, 1e-3)
  bulk_marg.nsb_mean ~ dnorm(bulk_marg.nsb.m, 1/bulk_marg.nsb.sd^2)         
  bulk_marg.nsb_tau ~ dgamma(1e3, 1e-3)

  # Site-level NSBs
  for (i in 1:n.sites){
    bf.nsb_site[i] ~ dnorm(bf.nsb_mean, bf.nsb_tau)
    pf.nsb_site[i] ~ dnorm(pf.nsb_mean, pf.nsb_tau)
    brach.nsb_site[i] ~ dnorm(brach.nsb_mean, brach.nsb_tau)
    bivalve.nsb_site[i] ~ dnorm(bivalve.nsb_mean, bivalve.nsb_tau)
    amm.nsb_site[i] ~ dnorm(amm.nsb_mean, amm.nsb_tau)
    bel.nsb_site[i] ~ dnorm(bel.nsb_mean, bel.nsb_tau)
    bulk.nsb_site[i] ~ dnorm(bulk.nsb_mean, bulk.nsb_tau)
    micrite.nsb_site[i] ~ dnorm(micrite.nsb_mean, micrite.nsb_tau)
    bulk_sr.nsb_site[i] ~ dnorm(bulk_sr.nsb_mean, bulk_sr.nsb_tau)
    bulk_marg.nsb_site[i] ~ dnorm(bulk_marg.nsb_mean, bulk_marg.nsb_tau)
  }
  
  # Proxy system model 
  ####################################################################################################
  for (i in 1:length(ai.flat)){
    
      # Equations 4, 5, and 3 (respectively) of Tipple et al. (2010); bottom water carb chem fractionation
      eps.bicarb_co2_bot[i] <- -0.1141*tempC_bot[i] + 10.78                   
      eps.ci_co2_bot[i] <- 0.0049*tempC_bot[i] - 1.31
      eps.dic_co2_bot[i] <- (1-f_co3-f_carbacid)*(eps.bicarb_co2_bot[i]) + (f_co3*eps.ci_co2_bot[i])

      # Romanek et al. (1992) bulk calcite archives
      eps.cc_co2_surf[i] <- cc_co2_constant1 - cc_co2_coeff1*tempC[i]
      d13Cbulk_secular[i] <- d13CO2[ai.flat[i]] + eps.cc_co2_surf[i]

      # Calculate various component-derived carbonate archive d13C values
      d13Cpf[i] <- d13Cbulk_secular[i] + pf.nsb_site[si.flat[i]]
      d13Cbrach[i] <- d13Cbulk_secular[i] + brach.nsb_site[si.flat[i]]
      d13Cbivalve[i] <- d13Cbulk_secular[i] + bivalve.nsb_site[si.flat[i]]
      d13Camm[i] <- d13Cbulk_secular[i] + amm.nsb_site[si.flat[i]]
      d13Cbel[i] <- d13Cbulk_secular[i] + bel.nsb_site[si.flat[i]]
      
      # Benthic forams; equation 7 of Tipple et al. (2010)
      d13Cbf[i] <- ((d13CO2[ai.flat[i]]+1000)*((eps.dic_co2_bot[i]/1000)+1)) - eps.dic_cc 
      - A - 1000 + bf.nsb_site[si.flat[i]]
      
      # Calculate various bulk carbonate archive d13C values
      d13Cbulk[i] <- d13Cbulk_secular[i] + bulk.nsb_site[si.flat[i]]
      d13Cmicrite[i] <- d13Cbulk_secular[i] + micrite.nsb_site[si.flat[i]]
      d13Cbulk_sr[i] <- d13Cbulk_secular[i] + bulk_sr.nsb_site[si.flat[i]]
      d13Cbulk_marg[i] <- d13Cbulk_secular[i] + bulk_marg.nsb_site[si.flat[i]]
  }
  
  
  # Time evolution model  
  ####################################################################################################
  
  # Priors on process SDs
  GMST_sigma ~ dunif(0, 0.5)
  GMST_tau <- 1 / (GMST_sigma^2)
  
  BWT_sigma ~ dunif(0, 0.5)
  BWT_tau <- 1 / (BWT_sigma^2)
  
  d13CO2_sigma ~ dunif(0, 0.1)         
  d13CO2_tau <- 1 / (d13CO2_sigma^2)
  
  # Initial states 
  GMST[1] ~ dnorm(GMST.obs[1], 1/GMST.sd[1]^2)
  BWT[1] ~ dnorm(BWT.obs[1],  1/BWT.sd[1]^2)
  d13CO2[1] ~ dunif(d13CO2.l, d13CO2.u)
  
  # State evolution
  for (i in 2:n.steps){
    d13CO2[i] ~ dnorm(d13CO2[i-1], d13CO2_tau)
    GMST[i] ~ dnorm(GMST[i-1], GMST_tau)
    BWT[i] ~ dnorm(BWT[i-1], BWT_tau)
  }
  
  # Observation equations (these impose the evaluation against prescribed means)
  for (i in 2:n.steps){
    GMST.obs[i] ~ dnorm(GMST[i], 1/GMST.sd[i]^2)
    BWT.obs[i] ~ dnorm(BWT[i], 1/BWT.sd[i]^2)
  }
  
  # Spatial offsets
  for (i in 1:length(ai.flat)){
    toff[i]  ~ dnorm(toff.m[i], 1/toff.sd[i]^2)
    toff_bot[i] ~ dnorm(0, 1/toff_sd_uniform_bot^2)
    
    tempC[i] <- GMST[ai.flat[i]] + toff[i]
    tempC_bot[i] <- BWT[ai.flat[i]] + toff_bot[i]
  }
}
