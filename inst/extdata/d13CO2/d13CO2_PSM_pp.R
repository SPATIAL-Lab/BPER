model {

  # Likelihood block
  ####################################################################################################

  for (i in 1:n.d13Cbf){
    d13Cbf.data[i] ~ dnorm(d13Cbf[ri.d13Cbf[i]], 1/d13Cbf.sd[i]^2)
  }
  for (i in 1:n.d13Cpf){
    d13Cpf.data[i] ~ dnorm(d13Cpf[ri.d13Cpf[i]], 1/d13Cpf.sd[i]^2)
  }
  for (i in 1:n.d13Cbrach){
    d13Cbrach.data[i] ~ dnorm(d13Cbrach[ri.d13Cbrach[i]], 1/d13Cbrach.sd[i]^2)
  }
  for (i in 1:n.d13Cbivalve){
    d13Cbivalve.data[i] ~ dnorm(d13Cbivalve[ri.d13Cbivalve[i]], 1/d13Cbivalve.sd[i]^2)
  }
  for (i in 1:n.d13Camm){
    d13Camm.data[i] ~ dnorm(d13Camm[ri.d13Camm[i]], 1/d13Camm.sd[i]^2)
  }
  for (i in 1:n.d13Cbel){
    d13Cbel.data[i] ~ dnorm(d13Cbel[ri.d13Cbel[i]], 1/d13Cbel.sd[i]^2)
  }
  for (i in 1:n.d13Cmicrite){
    d13Cmicrite.data[i] ~ dnorm(d13Cmicrite[ri.d13Cmicrite[i]], 1/d13Cmicrite.sd[i]^2)
  }
  for (i in 1:n.d13Cbulk){
    d13Cbulk.data[i] ~ dnorm(d13Cbulk[ri.d13Cbulk[i]], 1/d13Cbulk.sd[i]^2)
  }
  for (i in 1:n.d13Cbulk_sr){
    d13Cbulk_sr.data[i] ~ dnorm(d13Cbulk_sr[ri.d13Cbulk_sr[i]], 1/d13Cbulk_sr.sd[i]^2)
  }
  for (i in 1:n.d13Cbulk_marg){
    d13Cbulk_marg.data[i] ~ dnorm(d13Cbulk_marg[ri.d13Cbulk_marg[i]], 1/d13Cbulk_marg.sd[i]^2)
  }

  # Constants
  ####################################################################################################

  eps.dic_cc <- -1

  cc_co2_constant1 ~ dnorm(11.98, 1/0.13^2)T(11.72,12.24)
  cc_co2_coeff1 ~ dnorm(0.12, 1/0.01^2)T(0.1,0.14)

  asd ~ dnorm(1, 1/0.2^2)T(0,)
  bpump ~ dnorm(1.2, 1/0.4^2)T(0,)
  remin ~ dnorm(0.6, 1/0.3^2)T(0,)
  Abf <- asd + bpump + remin
  Asurf <- asd
  f_co3 ~ dnorm(0.08, 1/0.04^2)T(0, 0.16)
  f_carbacid ~ dnorm(0.01, 1/0.005^2)T(0, 0.02)
  f_bicarb <- 1-f_co3-f_carbacid

  # Non-secular bias
  ####################################################################################################

  level_archive_mean[1] <- d13CO2_level_prior_mean
  level_archive_mean[2] <- d13CO2_level_prior_mean-Abf
  level_archive_mean[3] <- d13CO2_level_prior_mean+cc_co2_constant1-Asurf
  level_archive_mean[4] <- d13CO2_level_prior_mean+cc_co2_constant1-Asurf
  level_archive_mean[5] <- d13CO2_level_prior_mean+cc_co2_constant1-Asurf
  level_archive_mean[6] <- d13CO2_level_prior_mean+cc_co2_constant1-Asurf
  level_archive_mean[7] <- d13CO2_level_prior_mean+cc_co2_constant1-Asurf
  level_archive_mean[8] <- d13CO2_level_prior_mean+cc_co2_constant1-Asurf
  level_archive_mean[9] <- d13CO2_level_prior_mean+cc_co2_constant1-Asurf
  level_archive_mean[10] <- d13CO2_level_prior_mean+cc_co2_constant1-Asurf
  level_archive_mean[11] <- d13CO2_level_prior_mean+cc_co2_constant1-Asurf
  level_archive_block[1:11] ~ dmnorm(level_archive_mean[], level_archive_precision[,])

  d13CO2_level <- level_archive_block[1]
  bf.archive_level <- level_archive_block[2]
  pf.archive_level <- level_archive_block[3]
  brach.archive_level <- level_archive_block[4]
  bivalve.archive_level <- level_archive_block[5]
  amm.archive_level <- level_archive_block[6]
  bel.archive_level <- level_archive_block[7]
  micrite.archive_level <- level_archive_block[8]
  bulk.archive_level <- level_archive_block[9]
  bulk_sr.archive_level <- level_archive_block[10]
  bulk_marg.archive_level <- level_archive_block[11]

  bf.nsb_mean <- d13CO2_level-Abf-bf.archive_level
  bf.nsb_tau ~ dgamma(1e3, 1e-3)
  for (i in 1:n.sites.bf){
    bf.eta_std[i] ~ dnorm(0, 1)
    bf.eta_site[i] <- bf.nsb_mean + bf.eta_std[i]/sqrt(bf.nsb_tau)
    bf.nsb_site[i] <- bf.eta_site[i] + Abf
  }

  pf.nsb_mean <- d13CO2_level+cc_co2_constant1-Asurf-pf.archive_level
  pf.nsb_tau ~ dgamma(1e3, 1e-3)
  for (i in 1:n.sites.pf){
    pf.eta_std[i] ~ dnorm(0, 1)
    pf.eta_site[i] <- pf.nsb_mean + pf.eta_std[i]/sqrt(pf.nsb_tau)
    pf.nsb_site[i] <- pf.eta_site[i] + Asurf
  }

  brach.nsb_mean <- d13CO2_level+cc_co2_constant1-Asurf-brach.archive_level
  brach.nsb_tau ~ dgamma(1e3, 1e-3)
  for (i in 1:n.sites.brach){
    brach.eta_std[i] ~ dnorm(0, 1)
    brach.eta_site[i] <- brach.nsb_mean + brach.eta_std[i]/sqrt(brach.nsb_tau)
    brach.nsb_site[i] <- brach.eta_site[i] + Asurf
  }

  bivalve.nsb_mean <- d13CO2_level+cc_co2_constant1-Asurf-bivalve.archive_level
  bivalve.nsb_tau ~ dgamma(1e3, 1e-3)
  for (i in 1:n.sites.bivalve){
    bivalve.eta_std[i] ~ dnorm(0, 1)
    bivalve.eta_site[i] <- bivalve.nsb_mean + bivalve.eta_std[i]/sqrt(bivalve.nsb_tau)
    bivalve.nsb_site[i] <- bivalve.eta_site[i] + Asurf
  }

  amm.nsb_mean <- d13CO2_level+cc_co2_constant1-Asurf-amm.archive_level
  amm.nsb_tau ~ dgamma(1e3, 1e-3)
  for (i in 1:n.sites.amm){
    amm.eta_std[i] ~ dnorm(0, 1)
    amm.eta_site[i] <- amm.nsb_mean + amm.eta_std[i]/sqrt(amm.nsb_tau)
    amm.nsb_site[i] <- amm.eta_site[i] + Asurf
  }

  bel.nsb_mean <- d13CO2_level+cc_co2_constant1-Asurf-bel.archive_level
  bel.nsb_tau ~ dgamma(1e3, 1e-3)
  for (i in 1:n.sites.bel){
    bel.eta_std[i] ~ dnorm(0, 1)
    bel.eta_site[i] <- bel.nsb_mean + bel.eta_std[i]/sqrt(bel.nsb_tau)
    bel.nsb_site[i] <- bel.eta_site[i] + Asurf
  }

  micrite.nsb_mean <- d13CO2_level+cc_co2_constant1-Asurf-micrite.archive_level
  micrite.nsb_tau ~ dgamma(1e3, 1e-3)
  for (i in 1:n.sites.micrite){
    micrite.eta_std[i] ~ dnorm(0, 1)
    micrite.eta_site[i] <- micrite.nsb_mean + micrite.eta_std[i]/sqrt(micrite.nsb_tau)
    micrite.nsb_site[i] <- micrite.eta_site[i] + Asurf
  }

  bulk.nsb_mean <- d13CO2_level+cc_co2_constant1-Asurf-bulk.archive_level
  bulk.nsb_tau ~ dgamma(1e3, 1e-3)
  for (i in 1:n.sites.bulk){
    bulk.eta_std[i] ~ dnorm(0, 1)
    bulk.eta_site[i] <- bulk.nsb_mean + bulk.eta_std[i]/sqrt(bulk.nsb_tau)
    bulk.nsb_site[i] <- bulk.eta_site[i] + Asurf
  }

  bulk_sr.nsb_mean <- d13CO2_level+cc_co2_constant1-Asurf-bulk_sr.archive_level
  bulk_sr.nsb_tau ~ dgamma(1e3, 1e-3)
  for (i in 1:n.sites.bulk_sr){
    bulk_sr.eta_std[i] ~ dnorm(0, 1)
    bulk_sr.eta_site[i] <- bulk_sr.nsb_mean + bulk_sr.eta_std[i]/sqrt(bulk_sr.nsb_tau)
    bulk_sr.nsb_site[i] <- bulk_sr.eta_site[i] + Asurf
  }

  bulk_marg.nsb_mean <- d13CO2_level+cc_co2_constant1-Asurf-bulk_marg.archive_level
  bulk_marg.nsb_tau ~ dgamma(1e3, 1e-3)
  for (i in 1:n.sites.bulk_marg){
    bulk_marg.eta_std[i] ~ dnorm(0, 1)
    bulk_marg.eta_site[i] <- bulk_marg.nsb_mean + bulk_marg.eta_std[i]/sqrt(bulk_marg.nsb_tau)
    bulk_marg.nsb_site[i] <- bulk_marg.eta_site[i] + Asurf
  }

  # Proxy system model
  ####################################################################################################

  for (i in 1:n.flat.bf){
    eps.bicarb_co2_bot[i] <- -0.1141*tempC_bot[ri.flat.bf[i]] + 10.78
    eps.ci_co2_bot[i] <- 0.0049*tempC_bot[ri.flat.bf[i]] - 1.31
    eps.dic_co2_bot[i] <- f_bicarb*eps.bicarb_co2_bot[i] + f_co3*eps.ci_co2_bot[i]
    d13Cbf[i] <- ((d13CO2[ai.flat.bf[i]]+1000)*((eps.dic_co2_bot[i]/1000)+1)) - eps.dic_cc
    - 1000 - d13CO2_level + bf.archive_level - bf.eta_std[si.flat.bf[i]]/sqrt(bf.nsb_tau)
  }

  for (i in 1:n.flat.pf){
    d13Cpf[i] <- d13CO2_delta[ai.flat.pf[i]] + pf.archive_level
    - cc_co2_coeff1*tempC[ri.flat.pf[i]] - pf.eta_std[si.flat.pf[i]]/sqrt(pf.nsb_tau)
  }
  for (i in 1:n.flat.brach){
    d13Cbrach[i] <- d13CO2_delta[ai.flat.brach[i]] + brach.archive_level
    - cc_co2_coeff1*tempC[ri.flat.brach[i]] - brach.eta_std[si.flat.brach[i]]/sqrt(brach.nsb_tau)
  }
  for (i in 1:n.flat.bivalve){
    d13Cbivalve[i] <- d13CO2_delta[ai.flat.bivalve[i]] + bivalve.archive_level
    - cc_co2_coeff1*tempC[ri.flat.bivalve[i]] - bivalve.eta_std[si.flat.bivalve[i]]/sqrt(bivalve.nsb_tau)
  }
  for (i in 1:n.flat.amm){
    d13Camm[i] <- d13CO2_delta[ai.flat.amm[i]] + amm.archive_level
    - cc_co2_coeff1*tempC[ri.flat.amm[i]] - amm.eta_std[si.flat.amm[i]]/sqrt(amm.nsb_tau)
  }
  for (i in 1:n.flat.bel){
    d13Cbel[i] <- d13CO2_delta[ai.flat.bel[i]] + bel.archive_level
    - cc_co2_coeff1*tempC[ri.flat.bel[i]] - bel.eta_std[si.flat.bel[i]]/sqrt(bel.nsb_tau)
  }

  for (i in 1:n.flat.bulk){
    d13Cbulk[i] <- d13CO2_delta[ai.flat.bulk[i]] + bulk.archive_level
    - cc_co2_coeff1*tempC[ri.flat.bulk[i]] - bulk.eta_std[si.flat.bulk[i]]/sqrt(bulk.nsb_tau)
  }
  for (i in 1:n.flat.micrite){
    d13Cmicrite[i] <- d13CO2_delta[ai.flat.micrite[i]] + micrite.archive_level
    - cc_co2_coeff1*tempC[ri.flat.micrite[i]] - micrite.eta_std[si.flat.micrite[i]]/sqrt(micrite.nsb_tau)
  }
  for (i in 1:n.flat.bulk_sr){
    d13Cbulk_sr[i] <- d13CO2_delta[ai.flat.bulk_sr[i]] + bulk_sr.archive_level
    - cc_co2_coeff1*tempC[ri.flat.bulk_sr[i]] - bulk_sr.eta_std[si.flat.bulk_sr[i]]/sqrt(bulk_sr.nsb_tau)
  }
  for (i in 1:n.flat.bulk_marg){
    d13Cbulk_marg[i] <- d13CO2_delta[ai.flat.bulk_marg[i]] + bulk_marg.archive_level
    - cc_co2_coeff1*tempC[ri.flat.bulk_marg[i]] - bulk_marg.eta_std[si.flat.bulk_marg[i]]/sqrt(bulk_marg.nsb_tau)
  }

  # Time evolution model
  ####################################################################################################

  GMST_sigma ~ dunif(0, 0.5)
  GMST_tau <- 1/(GMST_sigma^2)
  BWT_sigma ~ dunif(0, 0.5)
  BWT_tau <- 1/(BWT_sigma^2)
  d13CO2_sigma ~ dunif(0, d13CO2_sigma_upper)
  d13CO2_tau <- 1/(d13CO2_sigma^2)

  GMST[1] ~ dnorm(GMST.obs[1], 1/GMST.sd[1]^2)
  BWT[1] ~ dnorm(BWT.obs[1], 1/BWT.sd[1]^2)
  d13CO2_delta[1] <- 0
  for (i in 2:n.steps){
    GMST[i] ~ dnorm(GMST[i-1], GMST_tau/dt.scale[i-1])
    BWT[i] ~ dnorm(BWT[i-1], BWT_tau/dt.scale[i-1])
    d13CO2_delta[i] ~ dnorm(d13CO2_delta[i-1], d13CO2_tau/dt.scale[i-1])
  }
  for (i in 1:n.steps){
    d13CO2[i] <- d13CO2_level + d13CO2_delta[i]
  }

  for (i in 2:n.steps){
    GMST.obs[i] ~ dnorm(GMST[i], 1/GMST.sd[i]^2)
    BWT.obs[i] ~ dnorm(BWT[i], 1/BWT.sd[i]^2)
  }

  for (i in 1:length(ai.flat)){
    toff[i] ~ dnorm(toff.m[i], 1/toff.sd[i]^2)
    toff_bot[i] ~ dnorm(0, 1/toff_sd_uniform_bot^2)
    tempC[i] <- GMST[ai.flat[i]] + toff[i]
    tempC_bot[i] <- BWT[ai.flat[i]] + toff_bot[i]
  }
}
