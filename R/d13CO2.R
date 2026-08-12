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
#' @param n.iter Number of MCMC iterations. Defaults to 3e5.
#' @param n.chains Number of MCMC chains. Defaults to 6.
#' @param n.burnin Number of burn-in iterations. Defaults to 1e5.
#' @param n.thin Thinning interval. Defaults to 500.
#' @param parallel Logical. If `TRUE`, run `R2jags::jags.parallel()`. If
#'   `FALSE`, run `R2jags::jags()`. Defaults to `TRUE`.
#' @param parallel Logical. If `TRUE`, you get a timeseries plot returned. If
#'   `FALSE`, you do not get a plot returned. Defaults to `TRUE`.
#'
#' @returns Returns time series plot (optional) and list `inv_out` with elements:
#' \itemize{
#'   \item `d13CO2_timeseries_out`: data frame of output ages and JAGS summary
#'   statistics for the `d13CO2` parameter.
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
                   n.iter = 3e5,
                   n.chains = 6,
                   n.burnin = 1e5,
                   n.thin = 500,
                   parallel = TRUE,
                   plot.out = TRUE){

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
  d13C.analyt.sd <- 0.1
  d13CO2_sigma_upper <- 0.1
  d13CO2_level_prior_mean <- -6

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

  nsb.priors <- data.frame(
    archive = c("bf", "pf", "brach", "bivalve", "amm", "bel", "micrite", "bulk", "bulk_sr", "bulk_marg"),
    mean = 0,
    sd = c(bf.nsb.sd, pf.nsb.sd, brach.nsb.sd, bivalve.nsb.sd, amm.nsb.sd,
           bel.nsb.sd, micrite.nsb.sd, bulk.nsb.sd, bulk_sr.nsb.sd, bulk_marg.nsb.sd),
    stringsAsFactors = FALSE
  )


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

  prox.raw <- as.data.frame(utils::read.csv(file = plate_path))
  paleogeo.cols <- paste0(c("plng_", "plat_"), plate_model)
  climate.cols <- c("MAT", "GMST_Li22", "GMST_PhanDA", "GMST_PhanDA_hi",
                    "GMST_PhanDA_lo", "temp_offset", "temp_offset_PhanDA")
  prox.in <- cbind(prox.raw[,1:7], prox.raw[,paleogeo.cols], prox.raw[,climate.cols],
                   rep(x = toff_sd_uniform, times = nrow(prox.raw)))
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
  ai.all <- seq_along(ages)
  age.indices <- as.data.frame(cbind(ai.all, ages))
  names(age.indices) <- c("ai", "age")

  dt.scale <- abs(diff(ages))/1000
  n.steps <- as.numeric(length(dt.scale) + 1)
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
  offset.column <- if(temp_offset_model == "Li22") "temp_offset" else "temp_offset_PhanDA"
  offset.data <- data.frame(ai = prox.in$ai, site.index = prox.in$site.index,
                            temp_offset = prox.in[[offset.column]],
                            temp_offset_sd = prox.in$temp_offset_sd)
  offset.data <- offset.data[is.finite(offset.data$temp_offset),]
  cell.offsets <- stats::aggregate(cbind(temp_offset, temp_offset_sd) ~ ai + site.index,
                                   data = offset.data, FUN = mean, na.rm = TRUE)
  cell.key <- paste(cell.offsets$ai, cell.offsets$site.index, sep = "_")
  flattened.key <- paste(flattened$ai, flattened$site.index, sep = "_")
  offset.rows <- match(flattened.key, cell.key)
  flattened$temp_offset_interp <- cell.offsets$temp_offset[offset.rows]
  flattened$temp_offset_sd_interp <- cell.offsets$temp_offset_sd[offset.rows]
  missing.offsets <- !is.finite(flattened$temp_offset_interp)
  if(any(missing.offsets)){
    offset.site <- stats::aggregate(cbind(temp_offset, temp_offset_sd) ~ site.index,
                                    data = offset.data, FUN = mean, na.rm = TRUE)
    offset.site.rows <- match(flattened$site.index[missing.offsets], offset.site$site.index)
    flattened$temp_offset_interp[missing.offsets] <- offset.site$temp_offset[offset.site.rows]
    flattened$temp_offset_sd_interp[missing.offsets] <- offset.site$temp_offset_sd[offset.site.rows]
  }
  missing.offsets <- !is.finite(flattened$temp_offset_interp)
  if(any(missing.offsets)){
    flattened$temp_offset_interp[missing.offsets] <- stats::approx(
      prox.in$age, prox.in[[offset.column]], xout = flattened$ages[missing.offsets],
      rule = 2, ties = mean)$y
    flattened$temp_offset_sd_interp[missing.offsets] <- toff_sd_uniform
  }
  flattened <- flattened[order(flattened$ai, flattened$site.index), ]
  flattened$row.index <- 1:nrow(flattened)
  rownames(flattened) <- NULL


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

  archive.data.raw <- list(bf = clean.d13Cbf, pf = clean.d13Cpf,
                           brach = clean.d13Cbrach, bivalve = clean.d13Cbivalve,
                           amm = clean.d13Camm, bel = clean.d13Cbel,
                           micrite = clean.d13Cmicrite, bulk = clean.d13Cbulk,
                           bulk_sr = clean.d13Cbulk_sr, bulk_marg = clean.d13Cbulk_marg)

  level.archive.covariance <- matrix(3^2, nrow = 11, ncol = 11)
  diag(level.archive.covariance) <- c(3^2, 3^2 + nsb.priors$sd^2)
  level_archive_precision <- solve(level.archive.covariance)

  aggregate.archive <- function(x){
    if(nrow(x) == 0){
      return(data.frame(ai = integer(0), site.index = integer(0), d13C = numeric(0),
                        d13C.sd = numeric(0)))
    }
    key <- interaction(x$ai, x$site.index, drop = TRUE)
    ans <- lapply(split(x, key), function(z){
      data.frame(ai = z$ai[1], site.index = z$site.index[1], d13C = mean(z$d13C),
                 d13C.sd = d13C.analyt.sd/sqrt(nrow(z)))
    })
    ans <- do.call(rbind, ans)
    ans[order(ans$ai, ans$site.index),]
  }
  archive.data <- lapply(archive.data.raw, aggregate.archive)

  make.archive.index <- function(x){
    if(nrow(x) == 0){
      return(list(n.sites = 0L, n.flat = 0L, ai.flat = integer(0), si.flat = integer(0),
                  ri.flat = integer(0), ri.obs = integer(0)))
    }
    archive.sites <- sort(unique(x$site.index))
    archive.flat <- unique(x[c("ai", "site.index")])
    archive.flat <- archive.flat[order(archive.flat$ai, archive.flat$site.index),]
    archive.key <- paste(archive.flat$ai, archive.flat$site.index, sep = "_")
    global.key <- paste(flattened$ai, flattened$site.index, sep = "_")
    observation.key <- paste(x$ai, x$site.index, sep = "_")
    list(n.sites = length(archive.sites), n.flat = nrow(archive.flat),
         ai.flat = archive.flat$ai,
         si.flat = match(archive.flat$site.index, archive.sites),
         ri.flat = match(archive.key, global.key),
         ri.obs = match(observation.key, archive.key))
  }
  archive.indexes <- lapply(archive.data, make.archive.index)


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
  BWT.m <- stats::approx(BWT$age, BWT$BWT, xout = ages, method = "linear", ties = mean, rule = 2)
  BWT.sd <- stats::approx(BWT$age, BWT$BWT_2sd/2, xout = ages, method = "linear", ties = mean, rule = 2)
  BWT.m <- BWT.m[["y"]]
  BWT.sd <- BWT.sd[["y"]]


  ####################################################################################################
  ##### environmental prior and temperature offsets
  ####################################################################################################

  toff.m <- flattened$temp_offset_interp
  toff.sd <- flattened$temp_offset_sd_interp

  stopifnot(!any(is.na(GMST.m)), !any(is.na(GMST.sd)))
  stopifnot(!any(is.na(BWT.m)), !any(is.na(BWT.sd)))


  ####################################################################################################
  ##### select objects to pass to jags
  ####################################################################################################

  data.pass <- list("n.steps" = n.steps,
                    "dt.scale" = dt.scale,
                    "ai.flat" = ai.flat,
                    "GMST.obs" = GMST.m,
                    "GMST.sd" = GMST.sd,
                    "BWT.obs" = BWT.m,
                    "BWT.sd" = BWT.sd,
                    "toff_sd_uniform_bot" = toff_sd_uniform_bot,
                    "toff.m" = toff.m,
                    "toff.sd" = toff.sd,
                    "d13CO2_sigma_upper" = d13CO2_sigma_upper,
                    "d13CO2_level_prior_mean" = d13CO2_level_prior_mean,
                    "level_archive_precision" = level_archive_precision)

  stems <- c(bf = "d13Cbf", pf = "d13Cpf", brach = "d13Cbrach",
             bivalve = "d13Cbivalve", amm = "d13Camm", bel = "d13Cbel",
             micrite = "d13Cmicrite", bulk = "d13Cbulk",
             bulk_sr = "d13Cbulk_sr", bulk_marg = "d13Cbulk_marg")
  for(archive in names(archive.data)){
    x <- archive.data[[archive]]
    idx <- archive.indexes[[archive]]
    stem <- stems[[archive]]
    data.pass[[paste0(stem, ".data")]] <- x$d13C
    data.pass[[paste0(stem, ".sd")]] <- x$d13C.sd
    data.pass[[paste0("ri.", stem)]] <- idx$ri.obs
    data.pass[[paste0("n.", stem)]] <- nrow(x)
    data.pass[[paste0("n.sites.", archive)]] <- idx$n.sites
    data.pass[[paste0("n.flat.", archive)]] <- idx$n.flat
    data.pass[[paste0("ai.flat.", archive)]] <- idx$ai.flat
    data.pass[[paste0("si.flat.", archive)]] <- idx$si.flat
    data.pass[[paste0("ri.flat.", archive)]] <- idx$ri.flat
  }


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

  summarydf <- as.data.frame(jout$BUGSoutput$summary, check.names = FALSE)
  
  d13CO2_rows <- grep("^d13CO2\\[", rownames(summarydf))
  
  if(length(d13CO2_rows) != length(ages)){
    stop("Could not match the full d13CO2 summary output to the model age vector")
  }
  
  # define spinup index correctly
  keep_idx <- seq.int(n.spinup + 1, length(ages))
  
  d13CO2_timeseries_out <- data.frame(
    age = ages[keep_idx],
    summarydf[d13CO2_rows[keep_idx],
              c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat", "n.eff"),
              drop = FALSE],
    row.names = NULL,
    check.names = FALSE
  )
  
  
  if(plot.out){
    
    ####################################################################################################
    ##### Plot d13CO2 reconstruction
    ####################################################################################################
    
    plot(d13CO2_timeseries_out$age,
         d13CO2_timeseries_out$`50%`,
         type = "n",
         xlim = c(max(d13CO2_timeseries_out$age), min(d13CO2_timeseries_out$age)),
         ylim = range(c(d13CO2_timeseries_out$`2.5%`,
                        d13CO2_timeseries_out$`97.5%`),
                      na.rm = TRUE),
         xlab = "Age (ka)",
         ylab = expression(delta^13*C[CO[2]]),
         bty = "l",
         las = 1)
    
    polygon(c(d13CO2_timeseries_out$age, rev(d13CO2_timeseries_out$age)),
            c(d13CO2_timeseries_out$`97.5%`, rev(d13CO2_timeseries_out$`2.5%`)),
            col = adjustcolor("grey70", alpha.f = 0.5),
            border = NA)
    
    polygon(c(d13CO2_timeseries_out$age, rev(d13CO2_timeseries_out$age)),
            c(d13CO2_timeseries_out$`75%`, rev(d13CO2_timeseries_out$`25%`)),
            col = adjustcolor("grey40", alpha.f = 0.6),
            border = NA)
    
    lines(d13CO2_timeseries_out$age,
          d13CO2_timeseries_out$`50%`,
          lwd = 2)
    
    box()
  }
  
  
  ####################################################################################################
  ##### Run inversion
  ####################################################################################################       

  inv_out <- list("d13CO2_timeseries_out" = d13CO2_timeseries_out,
                  "jout" = jout,
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

  if(any(summarydf$n.eff < 50)){
    warning("Some parameters have n.eff statistical parameter values less than 50; consider running more iterations")
  }

  if(any(summarydf$Rhat > 1.03, na.rm = TRUE)){
    warning("Some parameters have Rhat statistical parameter values greater than 1.03; consider running more iterations")
  }

  class(inv_out) <- c("d13CO2_out", "inv_out")
  return(inv_out)
}
