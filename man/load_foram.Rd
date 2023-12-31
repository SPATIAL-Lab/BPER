% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_foram.R
\name{load_foram}
\alias{load_foram}
\title{Checks and formats foraminiferal proxy data}
\usage{
load_foram(foram_data, sheet = NULL)
}
\arguments{
\item{foram_data}{User inputs an 8 column R data object and assigns it to 'foram_data' in arguments. Object columns
should be organized as: age, d11B, d11B2s, d11Bspec, Mg/Ca, Mg/Ca2s, d18O, d18O2s. The order matters. Current 
d11Bspec options are: 'Grub' (G. ruber), 'Tsac' (T. sacculifer), 'Ouni' (O. universa), 'custom' (specify values 
with 'priors_foram_adj' function), or 'borate'. If 'csv' or 'xlsx' is to be used instead of R object, provide a 
character string path to file from directory.}

\item{sheet}{If '.xlsx' is being used provide a character string of the sheet name to read in.}
}
\value{
Returns list 'load_foram'. This list contains 'prox_in' data.frame which has been checked for data completeness
and compatibility with other package functions and 'obs_type' which indicates which measurements are included.
}
\description{
This function allows you to load, QC and format foram data to use in subsequent 'BPER' functions.
Function checks which data are present and whether measurement uncertainties are included.
If ages or uncertainties are not included the function stops. Function flags 'incomplete'
data object - i.e., missing Mg/Ca, or missing d18O, or missing d11B - but proceeds. Function also writes
"priors_in_foram.R" to parent directory which needs to be modified before 'priors_foram' function is used.
If no modifications to priors are needed, this can be turned off using 'load.priors' argument. Dependencies: 'readxl'.
}
\examples{
load_foram(foram_data)

}
