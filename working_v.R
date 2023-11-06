

# First, load in an ix8 (8 column, i row(s)) data object, with columns ordered as:
# 

x <- read.csv("PETM_1209data.csv")

# Run the 'load_foram' function on the input data. 

load_foram <- load_foram(foram_data = x)

# Uake the returned list from 'load_foram' and use it for '' argument in 'age_index' function. Specify: the 'age_units'
# used in the data set ('kyr' or 'Myr'), the time step type ('step_type'): 'regular' for evenly spaced intervals or 'every' for a time step 
# for every unique age in the data set. If 'regular' is specified, specify the step interval ('step_int'). 

age_index <- age_index(load_proxy = load_foram, age_units = 'kyr', step_type = 'regular', step_int = 10)

# Use 'clean_foram' function on the return list from 'age_index' to compile and structure observational data for the 
# MCMC inversion. 

clean_obs <- clean_foram(age_index = age_index)

# Revise and load prior distributions for the MCMC inversion. Access the package data file 'priors_in_foram.R':

file.edit('priors_in_foram.R')

# This will open a '.R' script file in an editor window (should default to RStudio if using RStudio). Revise this script 
# as needed to define the prior distributions. More details are included as comments in the 'priors_in_foram.R' file. 
# Save this '.R' file in the working directory or along a file path which you will specify in the next function.

# Load in the priors from the 'priors_in_foram.R' file you just saved to the directory. Also include the 'aif' return object 
# from the 'age_index' function. Specify which 2nd carb chem varaible to use in pH calcs. Type ?priors_forams for list of options. 
# Also specify the type of prior 'cc2ndparm.pt' to use for the 2nd carb chem variable, either as an autocorrelated time series 
# model with prior distribution defined at time step 1 ('t1') or as a time series with the prior distribution defined at each time
# step ('ts'). If 'ts' is selected, the user must include time series data in argument 'cc2ndparmTS' consisting of 3 columns: 
# age, mean value for 2nd carb chem variable, and sd for 2nd carb chem variable. 

loscar_dic <- read.csv("LOSCAR_dic.csv")

priors_foram <- priors_foram(priors = 'priors_in_foram.R', age_index = aif, cc2ndparm.vt = 'dic', cc2ndparm.pt = 'ts', cc2ndparmTS = loscar_dic)

# Run the MCMC inversion by supplying in the 'clean_foram' and 'priors_foram' return objects. Specify a vector of character strings 
# indicating the parameters you would like saved in the MCMC output ('save.parms' argument). Also indicate the number of iterations
# ('n.iter'), the number of chains ('.n.chains'), the number of burn-in iterations ('n.burnin') and the amount of thinning (n.thin;
# i.e., saves every 'n'th iteration). 

run_inversion(clean_obs = clean_obs, priors_foram = priors_foram,
                          save.parms = c("sal", "tempC", "xca", "xmg", "xso4", "d11Bsw", "d18Osw", "d18Osw.local","pco2", "dic", "pH", "press"),
                          n.iter = 10000, n.chains = 3, n.burnin = 3000, n.thin = 1)

# A data set with on the order of hundreds Mg/Ca and d18O and tens of d11B data takes ~1 minute per 10000 iterations on my laptop.  

# A list is returned with MCMC data object 'jout' as the first object. It is recommended that you take a look at the statistical
# parameters Rhat and n.eff in inv_outjout$BUGSoutput$summary to make sure they are acceptable 

View(inv_out$jout$BUGSoutput$summary)

# It is also recommended that chain convergence is cofirmed using trace plots (i.e., built in 'trace_plot' function).







