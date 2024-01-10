## ---- setup, include = FALSE--------------------------------------------------
knitr::opts_chunk$set(
  comment = "##",
  fig.width = 6, fig.height = 6
)

## ---- install, message=FALSE, warning=FALSE, results="hide"-------------------
library(devtools)
install_github("SPATIAL-Lab/BPER")
library(BPER)

## ---- loadData----------------------------------------------------------------
x <- data.frame(c(55956, 55956, 55954, 55932, 55926, 55924, 55901, 55899), 
                c(14.61, 15.46, NA, 14.49, 14.11, NA, 14.46, 13.94), 
                c(0.23, 0.22, NA, 0.34, 0.33, NA, 0.21, 0.29), 
                c("Tsac", "Grub", NA, "Grub", "Tsac", NA, "Grub", "Tsac"), 
                c(3.26, NA, 3.33, 5.04, NA, 4.83, NA, NA), 
                c(0.1, NA, 0.1, 0.15, NA, 0.14, NA, NA), 
                c(NA, NA, NA, -2.14, NA, -1.98, -2.19, -1.99), 
                c(NA, NA, NA, 0.16, NA, 0.16, 0.16, 0.16))

lf_out <- load_foram(foram_data = x)

## ---- ageIndex----------------------------------------------------------------
ai_out <- age_index(load_proxy = lf_out, age_units = 'kyr', step_type = 'regular', step_int = 10)

## ---- adjustParms-------------------------------------------------------------
parm_names <- c("seccal", "seccal.sd")
parm_values <- c(60, 5)

parms_foram_adj <- parms_foram_adj(parms2change = parm_names, change_values = parm_values)

## ---- loadPriors--------------------------------------------------------------

dic_pri <- data.frame(c(56200, 55934, 55928, 55919, 55900, 55883), 
                      c(2320, 2328, 2405, 2440, 2481, 2500), 
                      c(150, 150, 150, 150, 150, 150)) 

pf_out <- priors_foram(parms_foram_adj = parms_foram_adj, 
                       age_index = ai_out, 
                       cc2ndparm.vt = 'dic', 
                       cc2ndparm.pt = 'ts', 
                       cc2ndparmTS = dic_pri)

## ---- runInversion------------------------------------------------------------
inv_out <- run_inversion(priors = pf_out, 
                         save.parms = c("sal", "tempC", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", "d18Osw.local","pco2", "dic", "pH", "press"), 
                         n.iter = 100, 
                         n.chains = 3, 
                         n.burnin = 30, 
                         n.thin = 1)

## ---- viewSumm, results='hide'------------------------------------------------
inv_out$jout$BUGSoutput$summary 

## ---- tracePlot, fig.width = 3, fig.height = 3--------------------------------
trace_plot(inv_out = inv_out, parm = "pco2", n = 3)

## ---- tsPlots-----------------------------------------------------------------
post_plot(inv_out = inv_out, parm = "pco2", type = "CI")

post_plot(inv_out = inv_out, parm = "pco2", type = "draws", n.draws = 800)

## ---- densityPlots------------------------------------------------------------
post_plot_ind(inv_out = inv_out, parm = "pH", tstep_age = 55944, show.median = TRUE, show.legend = TRUE, show.prior = TRUE, leg.pos = "topright")

## ---- viewAges----------------------------------------------------------------
inv_out$ages_prox

## ---- carbChem----------------------------------------------------------------
carbchem(ccparm2 = "pco2", ccout = "alk", inv_out = inv_out)

