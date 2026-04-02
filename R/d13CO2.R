#' Reconstruct Phanerozoic atmospheric d13CO2
#'
#' Runs the Phanerozoic atmospheric d13CO2 Bayesian hierarchical model adapted
#' from the `d13CO2` repository using packaged proxy, temperature, and
#' paleogeographic datasets distributed with `BPER`.
#'
#' @param age.min Minimum age in ka. Must be between 0 and 540000.
#' @param age.max Maximum age in ka. Must be between 0 and 540000.
#' @param step.int Time-step interval in ka. Must be between 1 and 200000.
#' @param GMST_model Character string. Either `"PhanDA"` or `"Scotese21"`.
#' @param temp_offset_model Character string. Either `"PhanDA"` or `"Li22"`.
#' @param plate_model Character string. One of `"PALEOMAP"`,
#'   `"TorsvikCocks2017"`, `"MERDITH2021"`, or `"CAO2024"`.
#' @param n.iter Number of MCMC iterations. Defaults to 1000.
#' @param n.chains Number of MCMC chains. Defaults to 3.
#' @param n.burnin Number of burn-in iterations. Defaults to 300.
#' @param n.thin Thinning interval. Defaults to 1.
#' @param parallel Logical. If `TRUE`, run `R2jags::jags.parallel()`. If
#'   `FALSE`, run `R2jags::jags()`. Defaults to `TRUE`.
#'
#' @returns Returns list `inv_out` with elements:
#' \itemize{
#'   \item `jout`: the fitted `R2jags` output object.
#'   \item `ages`: vector of model timestep ages in ka.
#'   \item `age.indices`: data frame linking model age indices to ages.
#'   \item `prox.in`: processed proxy input data used by the model.
#'   \item `flattened`: flattened age-site input grid passed to JAGS.
#'   \item `sites`: data frame mapping site names to site indices.
#'   \item `save.parms`: parameters saved from the JAGS model.
#'   \item `settings`: list of key function settings used in the run.
#' }
#'
#' @details
#' This function expects the following files to be installed with `BPER` under
#' `inst/extdata/d13CO2/`:
#' \itemize{
#'   \item `d13CO2_PSM_pp.R`
#'   \item `CenozoicBWT.csv`
#'   \item `PhanCompWithTemp_PALEOMAP.csv`
#'   \item `PhanCompWithTemp_TorsvikCocks2017.csv`
#'   \item `PhanCompWithTemp_MERDITH2021.csv`
#'   \item `PhanCompWithTemp_CAO2024.csv`
#' }
#'
#' After adding this file to `R/`, run `devtools::document()` (or
#' `roxygen2::roxygenise()`) to generate the `.Rd` file in `man/` and update
#' `NAMESPACE`.
#'
#' @examples
#' \dontrun{
#' out <- d13CO2()
#'
#' out <- d13CO2(age.min = 0,
#'               age.max = 66000,
#'               step.int = 1000,
#'               GMST_model = "PhanDA",
#'               temp_offset_model = "Li22",
#'               plate_model = "PALEOMAP")
#' }
#'
#' @export

d13CO2 <- function(age.min = 0,
                   age.max = 540000,
                   step.int = 1000,
                   GMST_model = "PhanDA",
                   temp_offset_model = "Li22",
                   plate_model = "PALEOMAP",
                   n.iter = 1000,
                   n.chains = 3,
                   n.burnin = 300,
                   n.thin = 1,
                   parallel = TRUE){

  ####################################################################################################
  ##### argument checks
  ####################################################################################################

  if(!is.numeric(age.min) || length(age.min) != 1 || is.na(age.min) || age.min < 0 || age.min > 540000){
    stop("'age.min' must be a single numeric value between 0 and 540000 ka")
  }

  if(!is.numeric(age.max) || length(age.max) != 1 || is.na(age.max) || age.max < 0 || age.max > 540000){
    stop("'age.max' must be a single numeric value between 0 and 540000 ka")
  }

  if(age.min >= age.max){
    stop("'age.min' must be less than 'age.max'")
  }

  if(!is.numeric(step.int) || length(step.int) != 1 || is.na(step.int) || step.int < 1 || step.int > 200000){
    stop("'step.int' must be a single numeric value between 1 and 200000 ka")
  }

  if(!GMST_model %in% c("PhanDA", "Scotese21")){
    stop("'GMST_model' must be either 'PhanDA' or 'Scotese21'")
  }

  if(!temp_offset_model %in% c("PhanDA", "Li22")){
    stop("'temp_offset_model' must be either 'PhanDA' or 'Li22'")
  }

  if(!plate_model %in% c("PALEOMAP", "TorsvikCocks2017", "MERDITH2021", "CAO2024")){
    stop("'plate_model' must be one of 'PALEOMAP', 'TorsvikCocks2017', 'MERDITH2021', or 'CAO2024'")
  }

  if(!is.numeric(n.iter) || length(n.iter) != 1 || is.na(n.iter) || n.iter <= 0){
    stop("'n.iter' must be a single positive numeric value")
  }

  if(!is.numeric(n.chains) || length(n.chains) != 1 || is.na(n.chains) || n.chains <= 0){
    stop("'n.chains' must be a single positive numeric value")
  }

  if(!is.numeric(n.burnin) || length(n.burnin) != 1 || is.na(n.burnin) || n.burnin < 0){
    stop("'n.burnin' must be a single non-negative numeric value")
  }

  if(!is.numeric(n.thin) || length(n.thin) != 1 || is.na(n.thin) || n.thin <= 0){
    stop("'n.thin' must be a single positive numeric value")
  }

  if(!is.logical(parallel) || length(parallel) != 1 || is.na(parallel)){
    stop("'parallel' must be either TRUE or FALSE")
  }


  ####################################################################################################
  ##### internal settings retained from d13CO2_Driver_pp.R
  ####################################################################################################

  n.spinup <- 10
  GMST_sd_Scotese21 <- 5
  toff_sd_uniform <- 2
  toff_sd_uniform_bot <- 1

  bf.nsb.m <- 0
  bf.nsb.sd <- 0.25

  pf.nsb.m <- 0
  pf.nsb.sd <- 0.25

  brach.nsb.m <- 0
  brach.nsb.sd <- 1

  bivalve.nsb.m <- 0
  bivalve.nsb.sd <- 1

  amm.nsb.m <- 0
  amm.nsb.sd <- 1

  bel.nsb.m <- 0
  bel.nsb.sd <- 1

  micrite.nsb.m <- 0
  micrite.nsb.sd <- 0.25

  bulk.nsb.m <- 0
  bulk.nsb.sd <- 0.5

  bulk_sr.nsb.m <- 0
  bulk_sr.nsb.sd <- 1

  bulk_marg.nsb.m <- 0
  bulk_marg.nsb.sd <- 0.75


  ####################################################################################################
  ##### locate packaged auxiliary files
  ####################################################################################################

  plate_file <- switch(plate_model,
                       "PALEOMAP" = "PhanCompWithTemp_PALEOMAP.csv",
                       "TorsvikCocks2017" = "PhanCompWithTemp_TorsvikCocks2017.csv",
                       "MERDITH2021" = "PhanCompWithTemp_MERDITH2021.csv",
                       "CAO2024" = "PhanCompWithTemp_CAO2024.csv")

  plate_path <- system.file("extdata", "d13CO2", plate_file, package = "BPER")
  bwt_path <- system.file("extdata", "d13CO2", "CenozoicBWT.csv", package = "BPER")
  model_path <- system.file("extdata", "d13CO2", "d13CO2_PSM_pp.R", package = "BPER")

  if(plate_path == ""){
    stop(paste0("Could not find packaged plate model file: ", plate_file,
                ". Place it in 'inst/extdata/d13CO2/' before building/installing BPER."))
  }

  if(bwt_path == ""){
    stop("Could not find packaged file 'CenozoicBWT.csv'. Place it in 'inst/extdata/d13CO2/' before building/installing BPER.")
  }

  if(model_path == ""){
    stop("Could not find packaged file 'd13CO2_PSM_pp.R'. Place it in 'inst/extdata/d13CO2/' before building/installing BPER.")
  }


  ####################################################################################################
  ##### load and prepare proxy and climate data, indexing vectors and matrices
  ####################################################################################################

  prox.in <- as.data.frame(utils::read.csv(file = plate_path))
  prox.in <- cbind(prox.in[,1:7], prox.in[,9:10], prox.in[,21:27], rep(x = toff_sd_uniform, times = nrow(prox.in)))
  names(prox.in) <- c("age", "d13C", "source", "site", "lat", "lon", "category",
                      "paleolon", "paleolat", "MAT", "GMST_Scotese21", "GMST_PhanDA", "GMST_PhanDA_hi",
                      "GMST_PhanDA_lo", "temp_offset", "temp_offset_PhanDA", "temp_offset_sd")

  prox.in$age <- prox.in$age*1e3
  prox.in <- prox.in[prox.in$age >= age.min & prox.in$age <= age.max,]

  if(nrow(prox.in) == 0){
    stop("No proxy rows remain after applying the selected age range")
  }

  prox.in <- transform(prox.in, ai = n.spinup + as.numeric(1 + floor((max(prox.in$age) - prox.in$age) / step.int)))
  prox.in <- prox.in[order(prox.in$age, decreasing = TRUE),]

  ages.short <- seq(from = max(prox.in$age), to = min(prox.in$age), by = -1*step.int) - 0.5*step.int
  ages <- seq(from = n.spinup*step.int + max(prox.in$age), to = min(prox.in$age), by = -1*step.int) - 0.5*step.int
  ai.all <- c(c(1:n.spinup), sort(unique(prox.in$ai), decreasing = FALSE))
  age.indices <- as.data.frame(cbind(ai.all, ages))
  names(age.indices) <- c("ai", "age")

  dt <- abs(diff(unique(ages), lag = 1))
  n.steps <- as.numeric(length(dt) + 1)
  age.max.spinup <- age.max + step.int*n.spinup

  PhanDA_sd <- ((prox.in$GMST_PhanDA_hi - prox.in$GMST_PhanDA) +
                  (prox.in$GMST_PhanDA - prox.in$GMST_PhanDA_lo)) / 2

  if(GMST_model == "PhanDA"){
    GMST.m <- stats::approx(prox.in$age, prox.in$GMST_PhanDA, xout = ages, rule = 2, ties = mean)$y
    GMST.sd <- stats::approx(prox.in$age, PhanDA_sd, xout = ages, rule = 2, ties = mean)$y
  } else if(GMST_model == "Scotese21"){
    GMST.m <- stats::approx(prox.in$age, prox.in$GMST_Scotese21, xout = ages, rule = 2, ties = mean)$y
    GMST.sd <- rep(x = GMST_sd_Scotese21, times = length(ages))
  }

  prox.in <- transform(prox.in, site.index = as.numeric(factor(site, ordered = is.ordered(site))))
  site.index <- c(prox.in$site.index)
  n.sites <- as.numeric(length(unique(site.index)))

  sites <- data.frame((sort(unique(prox.in$site), decreasing = FALSE)), seq(1:length(sort(unique(prox.in$site), decreasing = FALSE))))
  names(sites) <- c("site", "site.index")

  flattened <- unique(prox.in[c("ai", "site.index")])
  flattened <- flattened[order(flattened$ai, flattened$site.index), ]
  ai.flat <- flattened$ai

  flattened$ages <- age.indices$age[match(flattened$ai, age.indices$ai)]

  flattened$GMST_PhanDA_interp <- stats::approx(prox.in$age, prox.in$GMST_PhanDA, xout = flattened$ages, rule = 2, ties = mean)$y
  flattened$GMST_Scotese21_interp <- stats::approx(prox.in$age, prox.in$GMST_Scotese21, xout = flattened$ages, rule = 2, ties = mean)$y
  flattened$GMST_PhanDA_sd_interp <- stats::approx(prox.in$age, PhanDA_sd, xout = flattened$ages, rule = 2, ties = mean)$y
  flattened$temp_offset_interp <- stats::approx(prox.in$age, prox.in$temp_offset, xout = flattened$ages, rule = 2, ties = mean)$y
  flattened$temp_offset_PhanDA_interp <- stats::approx(prox.in$age, prox.in$temp_offset_PhanDA, xout = flattened$ages, rule = 2, ties = mean)$y
  flattened$temp_offset_sd_interp <- stats::approx(prox.in$age, prox.in$temp_offset_sd, xout = flattened$ages, rule = 2, ties = mean)$y
  flattened <- flattened[order(flattened$ai, flattened$site.index), ]
  flattened$row.index <- 1:nrow(flattened)
  rownames(flattened) <- NULL
  si.flat <- flattened$site.index


  ####################################################################################################
  ##### clean and prepare proxy data
  ####################################################################################################

  clean.d13C <- prox.in[complete.cases(prox.in$d13C), ]
  clean.d13Cbf <- clean.d13C[clean.d13C$category == "bf",]
  clean.d13Cpf <- clean.d13C[clean.d13C$category == "Planktonic foraminifera",]
  clean.d13Cbrach <- clean.d13C[clean.d13C$category == "Brachiopod calcite",]
  clean.d13Cbivalve <- clean.d13C[clean.d13C$category == "Bivalve",]
  clean.d13Camm <- clean.d13C[clean.d13C$category == "Ammonite",]
  clean.d13Cbel <- clean.d13C[clean.d13C$category == "Belemnite",]
  clean.d13Cmicrite <- clean.d13C[clean.d13C$category == "micrite open ocean",]
  clean.d13Cbulk <- clean.d13C[clean.d13C$category %in% c("bulk", "bulk open water", "bulk open ocean"), ]
  clean.d13Cbulk_sr <- clean.d13C[clean.d13C$category == "bulk semi restricted",]
  clean.d13Cbulk_marg <- clean.d13C[clean.d13C$category %in% c("bulk marginal sea", "bulk marginal sea restricting up section"), ]

  ai.d13Cbf <- sort(c(as.integer(clean.d13Cbf$ai)), decreasing = FALSE)
  si.d13Cbf <- clean.d13Cbf$site.index
  d13Cbf.data <- clean.d13Cbf$d13C
  n.d13Cbf <- length(d13Cbf.data)

  ai.d13Cpf <- sort(c(as.integer(clean.d13Cpf$ai)), decreasing = FALSE)
  si.d13Cpf <- clean.d13Cpf$site.index
  d13Cpf.data <- clean.d13Cpf$d13C
  n.d13Cpf <- length(d13Cpf.data)

  ai.d13Cbrach <- sort(c(as.integer(clean.d13Cbrach$ai)), decreasing = FALSE)
  si.d13Cbrach <- clean.d13Cbrach$site.index
  d13Cbrach.data <- clean.d13Cbrach$d13C
  n.d13Cbrach <- length(d13Cbrach.data)

  ai.d13Cbivalve <- sort(c(as.integer(clean.d13Cbivalve$ai)), decreasing = FALSE)
  si.d13Cbivalve <- clean.d13Cbivalve$site.index
  d13Cbivalve.data <- clean.d13Cbivalve$d13C
  n.d13Cbivalve <- length(d13Cbivalve.data)

  ai.d13Camm <- sort(c(as.integer(clean.d13Camm$ai)), decreasing = FALSE)
  si.d13Camm <- clean.d13Camm$site.index
  d13Camm.data <- clean.d13Camm$d13C
  n.d13Camm <- length(d13Camm.data)

  ai.d13Cbel <- sort(c(as.integer(clean.d13Cbel$ai)), decreasing = FALSE)
  si.d13Cbel <- clean.d13Cbel$site.index
  d13Cbel.data <- clean.d13Cbel$d13C
  n.d13Cbel <- length(d13Cbel.data)

  ai.d13Cmicrite <- sort(c(as.integer(clean.d13Cmicrite$ai)), decreasing = FALSE)
  si.d13Cmicrite <- clean.d13Cmicrite$site.index
  d13Cmicrite.data <- clean.d13Cmicrite$d13C
  n.d13Cmicrite <- length(d13Cmicrite.data)

  ai.d13Cbulk <- sort(c(as.integer(clean.d13Cbulk$ai)), decreasing = FALSE)
  si.d13Cbulk <- clean.d13Cbulk$site.index
  d13Cbulk.data <- clean.d13Cbulk$d13C
  n.d13Cbulk <- length(d13Cbulk.data)

  ai.d13Cbulk_sr <- sort(c(as.integer(clean.d13Cbulk_sr$ai)), decreasing = FALSE)
  si.d13Cbulk_sr <- clean.d13Cbulk_sr$site.index
  d13Cbulk_sr.data <- clean.d13Cbulk_sr$d13C
  n.d13Cbulk_sr <- length(d13Cbulk_sr.data)

  ai.d13Cbulk_marg <- sort(c(as.integer(clean.d13Cbulk_marg$ai)), decreasing = FALSE)
  si.d13Cbulk_marg <- clean.d13Cbulk_marg$site.index
  d13Cbulk_marg.data <- clean.d13Cbulk_marg$d13C
  n.d13Cbulk_marg <- length(d13Cbulk_marg.data)


  ####################################################################################################
  ##### index each row of data to flattened data frame combinations
  ####################################################################################################

  ri.d13Cbf <- flattened$row.index[
    ri.d13Cbf <- match(
      interaction(clean.d13Cbf$ai, clean.d13Cbf$site.index),
      interaction(flattened$ai, flattened$site.index))
  ]

  ri.d13Cpf <- flattened$row.index[
    ri.d13Cpf <- match(
      interaction(clean.d13Cpf$ai, clean.d13Cpf$site.index),
      interaction(flattened$ai, flattened$site.index))
  ]

  ri.d13Cbrach <- flattened$row.index[
    ri.d13Cbrach <- match(
      interaction(clean.d13Cbrach$ai, clean.d13Cbrach$site.index),
      interaction(flattened$ai, flattened$site.index))
  ]

  ri.d13Cbivalve <- flattened$row.index[
    ri.d13Cbivalve <- match(
      interaction(clean.d13Cbivalve$ai, clean.d13Cbivalve$site.index),
      interaction(flattened$ai, flattened$site.index))
  ]

  ri.d13Camm <- flattened$row.index[
    ri.d13Camm <- match(
      interaction(clean.d13Camm$ai, clean.d13Camm$site.index),
      interaction(flattened$ai, flattened$site.index))
  ]

  ri.d13Cbel <- flattened$row.index[
    ri.d13Cbel <- match(
      interaction(clean.d13Cbel$ai, clean.d13Cbel$site.index),
      interaction(flattened$ai, flattened$site.index))
  ]

  ri.d13Cmicrite <- flattened$row.index[
    ri.d13Cmicrite <- match(
      interaction(clean.d13Cmicrite$ai, clean.d13Cmicrite$site.index),
      interaction(flattened$ai, flattened$site.index))
  ]

  ri.d13Cbulk <- flattened$row.index[
    ri.d13Cbulk <- match(
      interaction(clean.d13Cbulk$ai, clean.d13Cbulk$site.index),
      interaction(flattened$ai, flattened$site.index))
  ]

  ri.d13Cbulk_sr <- flattened$row.index[
    ri.d13Cbulk_sr <- match(
      interaction(clean.d13Cbulk_sr$ai, clean.d13Cbulk_sr$site.index),
      interaction(flattened$ai, flattened$site.index))
  ]

  ri.d13Cbulk_marg <- flattened$row.index[
    ri.d13Cbulk_marg <- match(
      interaction(clean.d13Cbulk_marg$ai, clean.d13Cbulk_marg$site.index),
      interaction(flattened$ai, flattened$site.index))
  ]

  if(any(is.na(ri.d13Cbf)) || any(is.na(ri.d13Cbulk))){
    stop("Age-site indexing failed for one or more proxy observations")
  }

  stopifnot(all(ri.d13Cbf >= 1 & ri.d13Cbf <= nrow(flattened)))
  stopifnot(all(ri.d13Cpf >= 1 & ri.d13Cpf <= nrow(flattened)))
  stopifnot(all(ri.d13Cbrach >= 1 & ri.d13Cbrach <= nrow(flattened)))
  stopifnot(all(ri.d13Cbivalve >= 1 & ri.d13Cbivalve <= nrow(flattened)))
  stopifnot(all(ri.d13Camm >= 1 & ri.d13Camm <= nrow(flattened)))
  stopifnot(all(ri.d13Cbel >= 1 & ri.d13Cbel <= nrow(flattened)))
  stopifnot(all(ri.d13Cbulk >= 1 & ri.d13Cbulk <= nrow(flattened)))
  stopifnot(all(ri.d13Cmicrite >= 1 & ri.d13Cmicrite <= nrow(flattened)))
  stopifnot(all(ri.d13Cbulk_sr >= 1 & ri.d13Cbulk_sr <= nrow(flattened)))
  stopifnot(all(ri.d13Cbulk_marg >= 1 & ri.d13Cbulk_marg <= nrow(flattened)))
  stopifnot(length(unique(flattened$row.index)) == nrow(flattened))


  ####################################################################################################
  ##### BWT interpolated for time steps
  ####################################################################################################

  BWT.Cen <- as.data.frame(utils::read.csv(file = bwt_path))
  names(BWT.Cen) <- c("age", "BWT", "BWT_2sd")
  BWT.Cen$age <- BWT.Cen$age*1e3
  BWT.Cen <- BWT.Cen[order(BWT.Cen$age, decreasing = TRUE),]
  BWT.Cen_last <- cbind(age.max.spinup, BWT.Cen[1,2:3])
  names(BWT.Cen_last) <- c("age", "BWT", "BWT_2sd")
  BWT <- rbind(BWT.Cen, BWT.Cen_last)
  BWT <- BWT[order(BWT$age, decreasing = TRUE),]
  BWT.m <- stats::approx(BWT$age, BWT$BWT, xout = ages, method = "linear", ties = mean)
  BWT.sd <- stats::approx(BWT$age, BWT$BWT_2sd/2, xout = ages, method = "linear", ties = mean)
  BWT.m <- BWT.m[["y"]]
  BWT.sd <- BWT.sd[["y"]]


  ####################################################################################################
  ##### environmental prior and temperature offsets
  ####################################################################################################

  d13CO2.l <- -12
  d13CO2.u <- 0

  if(temp_offset_model == "Li22"){
    toff.m <- flattened$temp_offset_interp
    toff.sd <- flattened$temp_offset_sd_interp
  } else if(temp_offset_model == "PhanDA"){
    toff.m <- flattened$temp_offset_PhanDA_interp
    toff.sd <- flattened$temp_offset_sd_interp
  }

  stopifnot(!any(is.na(GMST.m)), !any(is.na(GMST.sd)))
  stopifnot(!any(is.na(BWT.m)), !any(is.na(BWT.sd)))


  ####################################################################################################
  ##### select objects to pass to jags
  ####################################################################################################

  data.pass <- list("n.steps" = n.steps,
                    "dt" = dt,
                    "n.sites" = n.sites,
                    "si.flat" = si.flat,
                    "ai.flat" = ai.flat,
                    "GMST.obs" = GMST.m,
                    "GMST.sd" = GMST.sd,
                    "BWT.obs" = BWT.m,
                    "BWT.sd" = BWT.sd,
                    "toff_sd_uniform_bot" = toff_sd_uniform_bot,
                    "toff.m" = toff.m,
                    "toff.sd" = toff.sd,
                    "d13CO2.l" = d13CO2.l,
                    "d13CO2.u" = d13CO2.u)

  data.pass.bf <- list("d13Cbf.data" = d13Cbf.data,
                       "ai.d13Cbf" = ai.d13Cbf,
                       "ri.d13Cbf" = ri.d13Cbf,
                       "n.d13Cbf" = n.d13Cbf,
                       "bf.nsb.m" = bf.nsb.m,
                       "bf.nsb.sd" = bf.nsb.sd)

  data.pass.pf <- list("d13Cpf.data" = d13Cpf.data,
                       "ai.d13Cpf" = ai.d13Cpf,
                       "ri.d13Cpf" = ri.d13Cpf,
                       "n.d13Cpf" = n.d13Cpf,
                       "pf.nsb.m" = pf.nsb.m,
                       "pf.nsb.sd" = pf.nsb.sd)

  data.pass.brach <- list("d13Cbrach.data" = d13Cbrach.data,
                          "ai.d13Cbrach" = ai.d13Cbrach,
                          "ri.d13Cbrach" = ri.d13Cbrach,
                          "n.d13Cbrach" = n.d13Cbrach,
                          "brach.nsb.m" = brach.nsb.m,
                          "brach.nsb.sd" = brach.nsb.sd)

  data.pass.bivalve <- list("d13Cbivalve.data" = d13Cbivalve.data,
                            "ai.d13Cbivalve" = ai.d13Cbivalve,
                            "ri.d13Cbivalve" = ri.d13Cbivalve,
                            "n.d13Cbivalve" = n.d13Cbivalve,
                            "bivalve.nsb.m" = bivalve.nsb.m,
                            "bivalve.nsb.sd" = bivalve.nsb.sd)

  data.pass.amm <- list("d13Camm.data" = d13Camm.data,
                        "ai.d13Camm" = ai.d13Camm,
                        "ri.d13Camm" = ri.d13Camm,
                        "n.d13Camm" = n.d13Camm,
                        "amm.nsb.m" = amm.nsb.m,
                        "amm.nsb.sd" = amm.nsb.sd)

  data.pass.bel <- list("d13Cbel.data" = d13Cbel.data,
                        "ai.d13Cbel" = ai.d13Cbel,
                        "ri.d13Cbel" = ri.d13Cbel,
                        "n.d13Cbel" = n.d13Cbel,
                        "bel.nsb.m" = bel.nsb.m,
                        "bel.nsb.sd" = bel.nsb.sd)

  data.pass.micrite <- list("d13Cmicrite.data" = d13Cmicrite.data,
                            "ai.d13Cmicrite" = ai.d13Cmicrite,
                            "ri.d13Cmicrite" = ri.d13Cmicrite,
                            "n.d13Cmicrite" = n.d13Cmicrite,
                            "micrite.nsb.m" = micrite.nsb.m,
                            "micrite.nsb.sd" = micrite.nsb.sd)

  data.pass.bulk <- list("d13Cbulk.data" = d13Cbulk.data,
                         "ai.d13Cbulk" = ai.d13Cbulk,
                         "ri.d13Cbulk" = ri.d13Cbulk,
                         "n.d13Cbulk" = n.d13Cbulk,
                         "bulk.nsb.m" = bulk.nsb.m,
                         "bulk.nsb.sd" = bulk.nsb.sd)

  data.pass.bulk_sr <- list("d13Cbulk_sr.data" = d13Cbulk_sr.data,
                            "ai.d13Cbulk_sr" = ai.d13Cbulk_sr,
                            "ri.d13Cbulk_sr" = ri.d13Cbulk_sr,
                            "n.d13Cbulk_sr" = n.d13Cbulk_sr,
                            "bulk_sr.nsb.m" = bulk_sr.nsb.m,
                            "bulk_sr.nsb.sd" = bulk_sr.nsb.sd)

  data.pass.bulk_marg <- list("d13Cbulk_marg.data" = d13Cbulk_marg.data,
                              "ai.d13Cbulk_marg" = ai.d13Cbulk_marg,
                              "ri.d13Cbulk_marg" = ri.d13Cbulk_marg,
                              "n.d13Cbulk_marg" = n.d13Cbulk_marg,
                              "bulk_marg.nsb.m" = bulk_marg.nsb.m,
                              "bulk_marg.nsb.sd" = bulk_marg.nsb.sd)

  data.pass <- c(data.pass, data.pass.bf, data.pass.pf, data.pass.brach, data.pass.bivalve, data.pass.amm,
                 data.pass.bel, data.pass.micrite, data.pass.bulk, data.pass.bulk_sr, data.pass.bulk_marg)


  ####################################################################################################
  ##### parameters to save as output
  ####################################################################################################

  parms <- c("d13CO2", "GMST", "BWT", "tempC", "tempC_bot", "toff", "toff_bot", "d13Cbf", "d13Cpf",
             "d13Cbrach", "d13Cbivalve", "d13Camm", "d13Cbel", "d13Cbulk", "d13Cbulk_sr", "d13Cbulk_marg",
             "d13Cmicrite", "bf.nsb_site", "pf.nsb_site", "brach.nsb_site", "bivalve.nsb_site", "amm.nsb_site",
             "bel.nsb_site", "bulk.nsb_site", "micrite.nsb_site", "bulk_sr.nsb_site", "bulk_marg.nsb_site")


  ####################################################################################################
  ##### run the inversion using jags
  ####################################################################################################

  if(isTRUE(parallel)){
    jout <- R2jags::jags.parallel(data = data.pass,
                                  model.file = model_path,
                                  parameters.to.save = parms,
                                  inits = NULL,
                                  n.chains = n.chains,
                                  n.iter = n.iter,
                                  n.burnin = n.burnin,
                                  n.thin = n.thin)
  } else{
    jout <- R2jags::jags(data = data.pass,
                         model.file = model_path,
                         parameters.to.save = parms,
                         inits = NULL,
                         n.chains = n.chains,
                         n.iter = n.iter,
                         n.burnin = n.burnin,
                         n.thin = n.thin)
  }


  ####################################################################################################
  ##### assemble output
  ####################################################################################################

  inv_out <- list("jout" = jout,
                  "ages" = ages,
                  "ages.short" = ages.short,
                  "age.indices" = age.indices,
                  "prox.in" = prox.in,
                  "flattened" = flattened,
                  "sites" = sites,
                  "save.parms" = parms,
                  "settings" = list("age.min" = age.min,
                                     "age.max" = age.max,
                                     "step.int" = step.int,
                                     "GMST_model" = GMST_model,
                                     "temp_offset_model" = temp_offset_model,
                                     "plate_model" = plate_model,
                                     "n.spinup" = n.spinup,
                                     "n.iter" = n.iter,
                                     "n.chains" = n.chains,
                                     "n.burnin" = n.burnin,
                                     "n.thin" = n.thin,
                                     "parallel" = parallel,
                                     "GMST_sd_Scotese21" = GMST_sd_Scotese21,
                                     "toff_sd_uniform" = toff_sd_uniform,
                                     "toff_sd_uniform_bot" = toff_sd_uniform_bot),
                  "data.pass" = data.pass)

  summarydf <- data.frame(jout$BUGSoutput$summary)

  if(any(summarydf$n.eff < 50)){
    warning("Some parameters have n.eff statistical parameter values less than 50; consider running more iterations")
  }

  if(any(summarydf$Rhat > 1.03, na.rm = TRUE)){
    warning("Some parameters have Rhat statistical parameter values greater than 1.03; consider running more iterations")
  }

  class(inv_out) <- c("d13CO2_out", "inv_out")
  return(inv_out)
}
