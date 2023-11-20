
# User-editable model parameters and prior distributions, with default values

# Default values are representative of Paleogene 'hothouse' climate and mid-latitude localities

########################################################################################
# SALINITY TEMPERATURE AND PRESSURE 
########################################################################################

# temperature (degrees C) mean value and 1 sd for time step 1 
tempC.m = 30
tempC.sd = 5

# salinity mean value and 1sd for time step 1
sal.m = 35
sal.sd = 0.5

# pressure (bar) mean value and 1sd for time step 1
press.m = 6
press.sd = 1


########################################################################################
# STABLE ISOTOPE COMPOSITION AND MAJOR ION CONCENTRATIONS
########################################################################################

# d11B of seawater (per mille SRM-951) mean value and 1sd for time step 1 (default from Anagnostou et al. 2016)
d11Bsw.m = 38.45
d11Bsw.sd = 0.5

# d18O of seawater (per mille SMOW) mean value and 1sd for time step 1
d18Osw.m = -1.2
d18Osw.sd = 0.2

# [Ca] of seawater (mmol/kg) mean value and 1sd for time step 1
xca.m = 21.41
xca.sd = 0.5
# Prescribe linear change in [Ca] of seawater (mmol/kg) as function of age (mmol/kg per kyr)
# negative = decrease over time, positive = increase over time; 0 = no prescribed change 
# Defaults to long-term linearly modeled decline in [Ca] over the past 120 Myr
xca.lt = -1.9e-4

# [Mg] of seawater (mmol/kg) mean value and 1sd for time step 1
xmg.m = 68.51
xmg.sd = 0.5
# Prescribe linear change in [Mg] of seawater (mmol/kg) as function of age (mmol/kg per kyr)
# negative = decrease over time, positive = increase over time; 0 = no prescribed change 
xmg.lt = 0

# [SO4] of seawater (mmol/kg) mean value and 1sd for time step 1
xso4.m = 14
xso4.sd = 0.5
# Prescribe linear change in [SO4] of seawater (mmol/kg) as function of age (mmol/kg per kyr)
# negative = decrease over time, positive = increase over time; 0 = no prescribed change 
xso4.lt = 0


########################################################################################
# DIAGENESIS EFFECT ON d18O 
########################################################################################

# Percent recrystallization (%) mean value and 1sd; i.e., % of foram d18O that is from secondary calcite
# Defaults to correction off 
seccal = 0
seccal.sd = 0
# d18O (per mille VPDB) of secondary calcite
d18Oseccal = 0.85


########################################################################################
# Mg/Ca - TEMPERATURE CALIBRATION COEFFECIENTS 
########################################################################################

# Hp mean value and 1 sd; nonlinearity of the relationship b/w shell and Mg/Casw (Default is Haynes et al. (2023), T. sacculifer) 
Hp.mean = 0.74
Hp.sd = 0.05

# B (modern) mean value and 1 sd; modern (pre-corrected) pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016)
Bmod.mean = 0.38
Bmod.sd = 0.02

# A mean value and 1 sd; exponential constant in Mg/Ca-SST calibration (Default is 'Eocene' value from Evans et al., 2016)
A.mean = 0.0757
A.sd = 0.0045

# pH correction on Mg/Caf (% per tenth pH unit) mean value and sd; defaults to correction off
pHpccorr = 0
pHpccorrsd = 0


########################################################################################
# D11B VITAL EFFECT 
########################################################################################

# specify custom vital effect slope and intercept 
m.custom = 0.9
m.customsd = 0.2
c.custom = 2
c.customsd = 2

# 'c' intercept offsets for each modern species 'c' value; leave '0' for modern 
Grub.coff = 5.76 # value from Harper et al., in review for PNAS
Tsac.coff = 3.62 # value from Harper et al., in review for PNAS
Ouni.coff = 0
borate.coff = 0


########################################################################################
# CARBONATE CHEMISTRY 
########################################################################################

# pH ('total scale equivalent') upper '.u' and lower '.l' bounds on uniform distribution for time step 1
pH.u = 7.75
pH.l = 7.45

# Provide 2nd carbonate chemistry variable mean value and 1sd for time step 1
# Variable type ('cc2ndparm.vt') to be specified in arguments in 'foram_priors' function
# These values will be used if prior type ('cc2ndparm.pt') is set to 't1' (argument in 'foram_priors' function)
# These values will not be used if prior type ('cc2ndparm.pt') is set to 'ts' (argument in 'foram_priors' function)
# Units for each variable are listed below.
carbchem2.m = 2000
carbchem2.sd = 150

# Carbonate chemistry variables in the following units:
# pH = 'total scale' equivalent 
# DIC = μmol/kg
# ALK = μmol/kg
# [CO3] = μmol/kg
# [HCO3] = μmol/kg


