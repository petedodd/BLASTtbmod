##' Create parameter object for model
##'
##' Defaults for various things needed for model to run
##'
##' @title Create parameter object
##' @param start_year Earliest year of simulation
##' @param years N years in simulation
##' @param Dinit Matrix of initial prevalences
##' @param ari0 initial ARIs for initial state
##' @param hivoffset How many years ahead is Blantyre HIV incidence than MWI?
##' @param hivfac HIV incidence in Blantyre relative to MWI
##' @param hivdecline HIV decline rate from peak: overrides if >0
##' @param hiv_init_override if >0 provides initial HIV prevalence
##' @param ART_haz provides an ART initiation hazard
##' @param ART_init_override if >0 provides initial ART prevalence
##' @param hiv_checking if TRUE, prints HIV-related parameters and plots HIV incidence and ART initiation over time
##' @return list
##' @author Pete Dodd
##' @export
##' @import data.table
##' @import logitnorm
get.parms <- function(start_year,
                      years,
                      Dinit,
                      ari0,
                      hivoffset = 0,
                      hivfac = 2,
                      hivdecline = 0,
                      hiv_init_override = -1,
                      ART_haz = 0.5,
                      ART_init_override = -1,
                      hiv_checking = FALSE) {
  ########## Model dimensions & simulation parameters required for setup ############
  patch_dims <- 7 # number of patches = 3x3 grid     # Put in func
  age_dims <- 3 # Number of age groups             # put in func
  HIV_dims <- 3 # None, HIV, ART                   # put in func
  dt <- 1 / 12 # NOTE monthly timesteps hardcoded
  tt <- as.double(seq(0, 12 * years)) # Vector of timesteps
  sim_length <- length(tt) # Length of simulation
  offsetPops <- rep(1, patch_dims)

  ## birthrates
  birthrate_dat <- data.table::merge.data.table(
    BLASTtbmod::MWI$B[
      Year %in% start_year:(start_year + years),
      list(Year, Births)
    ],
    BLASTtbmod::MWI$N[
      Year %in% start_year:(start_year + years),
      list(N = sum(PopTotal)),
      by = Year
    ],
    by = "Year"
  )
  birthrate_dat[, birthrate := 1e3 * Births / N]
  births_int <- approx(birthrate_dat$birthrate, n = sim_length)$y

  ## initial age frax
  age.frax <- BLASTtbmod::MWI$N[
    Year == start_year,
    PopTotal
  ] # NOTE requires ages in order
  age.frax <- age.frax / sum(age.frax)

  ## Get correct years for HIV
  hiv_series <- data.table::data.table(
    hinc = BLASTtbmod::HIV_inc_1990_2021$Adult_incidence_1000_uninfected_est,
    year = BLASTtbmod::HIV_inc_1990_2021$Year
  )
  ## if hivdecline > 0, use a decline rate from the peak
  if (hivdecline > 0) {
    peakloc <- which.max(hiv_series$hinc)
    afterpeak <- peakloc:nrow(hiv_series)
    declinefac <- exp(-hivdecline * (afterpeak - peakloc))
    hiv_series$hinc[afterpeak] <- hiv_series$hinc[peakloc] * declinefac
  }
  offset <- hivoffset
  log_hinc <- Hmisc::approxExtrap(
    x = hiv_series$year,
    y = log(hiv_series$hinc),
    xout = (start_year + offset):(offset + start_year + years),
    method = "linear",
    rule = 2
  )$y
  hinc <- exp(log_hinc) * hivfac

  ## interpolate HIV
  HIV_int <- approx(
    hinc,
    n = sim_length
  )$y

  ## Time varying ART initiation - interpolate
  ## ART_int <- approx(BLASTtbmod::HIV_inc_1990_2021$ART_inc_est, n = sim_length)$y
  ART_int <- rep(ART_haz, sim_length)

  ## Risk-modifiers for TB based on HIV status
  TB_HIV_mod <- c(1, 1, 1) # infection by HIV
  Hirr <- c(1, 20, 4.4) # progression IRR by HIV

  ## Background death rates by HIV index
  ## set background death rate for ART same as no HIV
  tiz <- tz <- seq(
    from = start_year,
    to = (start_year + years), by = dt
  ) # Year timepoints
  tiz <- round(tiz)
  tiz <- as.integer(tiz) - start_year + 1 # which MU row to use

  ## create the death rate data structures
  ## NOTE: m_in still placeholder
  ## No longer interpolated - takes static value for corresponding year for each timestep
  mu_noHIV_int <- m_in_int <- matrix(data = 0, ncol = sim_length, nrow = age_dims)
  mu_ART_int <- mu_noHIV_int
  mu_HIV_int <- mu_noHIV_int
  MU <- BLASTtbmod::MWI$omega[
    Year %in% start_year:(start_year + years),
    list(Year, AgeGrp, mu = omegaT - r)
  ]
  MU <- as.matrix(
    dcast(MU, Year ~ AgeGrp, value.var = "mu")
  )[, -1] # mu per age group per year
  for (i in 1:ncol(mu_noHIV_int)) {
    for (j in 1:ncol(MU)) {
      val <- MU[tiz[i], j]
      mu_noHIV_int[j, i] <- ifelse(val > 0, val, 0) # net death/out
      m_in_int[j, i] <- ifelse(val < 0, -val, 0) # net in-migration
      mu_HIV_int[j, i] <- ifelse(j > 1, 0.1, 0)
      mu_ART_int[j, i] <- mu_noHIV_int[j, i]
    }
  }

  ## initial state for 'disease'
  if (missing(Dinit)) {
    Dinit <- matrix(0, nrow = patch_dims, ncol = age_dims)
    Dinit[, 2:age_dims] <- 1e-3
  }
  if (missing(ari0)) {
    ari0 <- qlnorm(0.5, log(2), 0.75) / 1e2
  }

  ## HIV initial state
  propinit_hiv <- array(0, c(
    patch_dims,
    age_dims,
    HIV_dims
  ))

  if (!start_year %in% hivp_mwi$Period) {
    stop("Need start year in HIV prevalence data range!")
  }

  if (ART_init_override >= 0) {
    artp <- ART_init_override
  } else {
    artp <- BLASTtbmod::hivp_mwi[
      Period == start_year & variable == "ARTpc", value
    ]
  }

  ## NOTE these HIV prevalences are in 2015
  ## if using different start years:
  ## adjust by corresponding HIV incidence relative to 2015
  if (hiv_init_override < 0) { # left as default
    adj_factor <- HIV_int[which(hiv_series$year == start_year)] /
      HIV_int[which(hiv_series$year == 2015)]
    hiv_init <- adj_factor * BLASTtbmod::blantyre$hivpre
  } else {
    ## use override but keep relative levels across patches
    hiv_init <- hiv_init_override *
      BLASTtbmod::blantyre$hivpre / mean(BLASTtbmod::blantyre$hivpre)
  }



  ## 15-49 & oldies
  for (j in 2:3) {
    propinit_hiv[, j, 2] <- hiv_init * (1 - artp)
    propinit_hiv[, j, 3] <- hiv_init * artp
  }
  ## complete HIV-
  for (j in 1:7) { # patch
    for (k in 1:3) { # zone
      propinit_hiv[j, k, 1] <- 1 - sum(propinit_hiv[j, k, 2:3])
    }
  }

  if (hiv_checking) {
    print(hiv_init)
    print(artp)
    print(range(HIV_int)) # NOTE per 1000 in stocm.R
    print(range(ART_int))
    plot(HIV_int, col = 2, ylim = c(0, 1.1 * max(HIV_int, ART_int)))
    points(ART_int, col = 3)
  }
  ## Set up list to pass to model
  parms <- list(
    dt = dt,
    tt = tt,
    Pi = pi,
    epsi = 0, # seasonality
    sim_length = sim_length,
    popinit = as.numeric(BLASTtbmod::blantyre$population), #2015
    patch_dims = patch_dims,
    HIV_dims = HIV_dims,
    age_dims = age_dims,
    age_rate = c(1 / 15, 1 / 35, 1e-6),
    agefracs = age.frax, # initial age fractions
    ageMids = c(15 / 2, (15 + 50) / 2, 60),
    initD = Dinit, # initial state
    propinit_hiv = propinit_hiv,
    births_int = births_int,
    HIV_int = HIV_int,
    ART_int = ART_int,
    TB_HIV_mod = TB_HIV_mod,
    Hirr = Hirr,
    mu_noHIV_int = mu_noHIV_int,
    mu_HIV_int = mu_HIV_int,
    mu_ART_int = mu_ART_int,
    m_in_int = m_in_int,
    IRR = rep(1, patch_dims),
    MM = (diag(patch_dims) * 3 + 0.5) / 3, # mixing matrix
    ## some HIV specifics
    HIV_dur_ratio = 6, # how much shorter TB in HIV+/ART-
    ART_det_OR = 2, # OR for detection in ART+
    ## TODO hyperparms separate data object
    progress_rate = 0.5, # Progression to clinical
    regress_rate = 0.05, # regression to subclinical
    beta = qlnorm(0.5, log(10), 0.75), # Effective contact rate - beta (bet)
    v = qbeta(0.5, 20.7, 77.9), # Partial infection protection - v (v)
    pLL = qlnorm(0.5, log(0.62), 0.068), # Stabilisation rate - pLL (arig)
    pDf = qlnorm(0.5, -2.837, 0.32), # Fast progression rate - pDf (pv)
    pDs = qlnorm(0.5, -6.89, 0.58), # Slow progression rate - pDs (eps)
    pDr = qlnorm(0.5, -3.95, 0.27), # Relapse rate - pDr (rel)
    cdr = logitnorm::qlogitnorm(0.5, 0, 0.3), # Case detection ratio - cdr
    cdr_SC = 0.05, # Case detection ratio for subclinical disease  (PLACEHOLDER)
    dur = qlnorm(0.5, 1.1, 0.2), # Duration untreated TB - dur
    tfr = qbeta(0.5, 2.71, 87.55), # CFR treated TB - tfr (txf)
    cfr = qbeta(0.5, 25.48, 33.78), # CFR untreated TB - cfr (cfrn)
    ari0 = ari0, # Initial condition parameter - ?? (ari0)
    ACFhaz0 = matrix(0.0, nrow = patch_dims, ncol = sim_length), # asymp ACF haz
    ACFhaz1 = matrix(0.0, nrow = patch_dims, ncol = sim_length) # symp ACF haz
  )
  return(parms)
}




##' Update parameters for given step using results of previous runs
##'
##' TODO
##'
##' @title Restart parameters at given step
##' @param parms Parameter list from previous runs, see 'get.parms()'
##' @param restart_step Step to restart at, in same units as 'tt' in parms
##' @param end_state Output from previous runs, in same format as output of 'stocm()' - should contain all timesteps up to and including restart_step
##' @return Updated parameter list with initial conditions and time-varying parameters updated to restart at given step
##' @author Pete Dodd
##' @export
##' @import data.table
restart_parms <- function(parms, restart_step, end_state) {
  cat("Creating new parameters to restart at given step...\n")
  ## gather snapshots for restart
  cat("Summarising end state...\n")
  denom <- extract.pops.multi(end_state, dim(end_state)[2], out_type = "N")
  numer <- extract.pops.multi(end_state, dim(end_state)[2], out_type = "D")
  hnumr <- denom[t == restart_step,
    .(N = sum(N)),
    by = .(patch, age, hiv, particle = chain_step)
  ]
  denom <- denom[t == restart_step,
    .(N = sum(N)),
    by = .(patch, age, particle = chain_step)
  ]
  numer <- numer[t == restart_step,
    .(D = sum(D)),
    by = .(patch, age, particle = chain_step)
  ]

  ## TB prevalence by patch and age
  both <- data.table::merge.data.table(
    denom, numer,
    by = c("patch", "age", "particle")
  )
  boths <- both[, .(prev = mean(D / N)), by = .(patch, age)]
  D00 <- data.table::dcast(boths, patch ~ age, value.var = "prev")
  D00 <- as.matrix(D00[, -1])

  ## HIV state
  hnumr[, tot := sum(N), by = .(patch, age, particle)]
  hnumr[, p := N / tot]
  hnumr <- hnumr[, .(p = mean(p)), by = .(patch, age, hiv)]
  H00 <- array(hnumr[order(hiv, age, patch)]$p,
    c(7, 3, 3),
    dimnames = list(
      patch = unique(hnumr$patch),
      age = unique(hnumr$age),
      hiv = unique(hnumr$hiv)
    )
  )

  cat("Updating parameter object...\n")
  ## rejig args:
  newparms <- parms
  newparms$initD <- D00
  newparms$propinit_hiv <- H00
  ## timed items:
  keep <- restart_step:length(parms$tt)
  newparms$tt <- newparms$tt[keep]
  newparms$tt <- newparms$tt - newparms$tt[1] # reset time to 0 at restart
  newparms$sim_length <- length(newparms$tt)
  newparms$births_int <- newparms$births_int[keep]
  newparms$HIV_int <- newparms$HIV_int[keep]
  newparms$ART_int <- newparms$ART_int[keep]
  newparms$mu_noHIV_int <- newparms$mu_noHIV_int[, keep]
  newparms$mu_HIV_int <- newparms$mu_HIV_int[, keep]
  newparms$mu_ART_int <- newparms$mu_ART_int[, keep]
  newparms$m_in_int <- newparms$m_in_int[, keep]
  newparms$ACFhaz0 <- newparms$ACFhaz0[, keep]
  newparms$ACFhaz1 <- newparms$ACFhaz1[, keep]

  return(newparms)
}
