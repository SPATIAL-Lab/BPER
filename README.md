# BPER
Bayesian Paleo-Environment Reconstructor

Dustin T. Harper & Gabriel J. Bowen 

Package developed for reconstructing the carbon cycle and climate from geochemical proxy data over geologic timescales. 

This package loads and grooms geochemical proxy data, writes a forward proxy system Bayesian heirarchical model, and runs 
a MCMC inversion to generate posterior distributions of user-specified proxy system parameters, and climate and carbon cycle 
variables. This version includes the foraminiferal proxy system, with plans to include phytoplankton, paleosol and terrestrial 
plant proxy systems. Marine equilibrium constant (Ks) are calculated following Zeebe & Tyrrell (2019) and marine carbonate 
chemistry calculations follow Zeebe & Wolf-Gladrow (2001). 