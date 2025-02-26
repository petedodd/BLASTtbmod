##' Create parameter object for model
##'
##' TODO
##' @title Create parameter object
##' @param start_year Earliest year of simulation
##' @param years N years in simulation
##' @param Dinit Matrix of initial prevalences
##' @param ari0 initial ARIs for initial state
##' @return list
##' @author Pete Dodd
##' @export
##' @import data.table
##' @import logitnorm
get.parms <- function(start_year,
                      years,
                      Dinit,
                      ari0) {
  ########## Model dimensions & simulation parameters required for setup ############
  patch_dims <- 7 # number of patches = 3x3 grid     # Put in func
  age_dims <- 3 # Number of age groups             # put in func
  HIV_dims <- 3 # None, HIV, ART                   # put in func
  dt <- 1 / 12 # NOTE monthly timesteps hardcoded
  tt <- as.double(seq(0, 12 * years)) # Vector of timesteps
  sim_length <- length(tt) # Length of simulation
  offsetPops <- rep(1, patch_dims)

  ## birthrates
  birthrate_dat <- merge(
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
  age.frax <- BLASTtbmod::MWI$N[Year == start_year, PopTotal] # NOTE requires ages in order
  age.frax <- age.frax / sum(age.frax)

  ## Get correct years for HIV
  HIV_inc_1990_2021 <- BLASTtbmod::HIV_inc_1990_2021[which(HIV_inc_1990_2021$Year %in% start_year:(start_year + years)), ]
  HIV_int <- approx(BLASTtbmod::HIV_inc_1990_2021$Adult_incidence_1000_uninfected_est, n = sim_length)$y

  ## Time varying ART initiation - interpolate
  ART_int <- approx(BLASTtbmod::HIV_inc_1990_2021$ART_inc_est, n = sim_length)$y

  ## Risk-modifiers for TB based on HIV status
  TB_HIV_mod <- c(1, 1, 1) # infection by HIV
  Hirr <- c(1, 10, 2) # progression IRR by HIV

  ## Background death rates by HIV index
  ## set background death rate for ART same as no HIV
  tiz <- tz <- seq(from = start_year, to = (start_year + years), by = dt) # Year timepoints
  tiz <- round(tiz)
  tiz <- as.integer(tiz) - start_year + 1 # which MU row to use

  ## create the death rate data structures
  ## NOTE: m_in still placeholder
  ## No longer interpolated - takes static value for corresponding year for each timestep
  mu_noHIV_int <- m_in_int <- matrix(data = 0, ncol = sim_length, nrow = age_dims)
  MU <- BLASTtbmod::MWI$omega[
    Year %in% start_year:(start_year + years),
    list(Year, AgeGrp, mu = omegaT - r)
  ]
  MU <- as.matrix(dcast(MU, Year ~ AgeGrp, value.var = "mu"))[, -1] # mu per age group per year
  for (i in 1:ncol(mu_noHIV_int)) {
    for (j in 1:ncol(MU)) {
      val <- MU[tiz[i], j]
      mu_noHIV_int[j, i] <- ifelse(val > 0, val, 0) # net death/out
      m_in_int[j, i] <- ifelse(val < 0, -val, 0) # net in-migration
    }
  }

  ## TODO
  ## mu_HIV_int <- matrix( data=NA, ncol = sim_length, nrow = age_dims )
  ## for( i in 1:age_dims ){
  ##   mu_noHIV_int[i,] <- approx( runif( years, min=0.01, max=0.07 ), n=sim_length )$y
  ##   mu_HIV_int[i,] <- approx( runif( 30*age_dims, min=0.04, max=0.09 ), n=sim_length )$y
  ## }
  mu_ART_int <- mu_noHIV_int
  mu_HIV_int <- mu_noHIV_int

  ## initial state for 'disease'
  if (missing(Dinit)) {
    Dinit <- matrix(0, nrow = patch_dims, ncol = age_dims)
    Dinit[, 2:age_dims] <- 1e-3 # assuming prevalence ~ 0 for kids TODO check parms for infectiousness
  }
  if(missing(ari0)){ #TODO think about making vector?
    ari0 <- qlnorm(0.5, log(2), 0.75)
  }

  ## Set up list to pass to model
  parms <- list(
    dt = dt,
    tt = tt,
    Pi = pi,
    epsi = 0, # seasonality
    sim_length = sim_length,
    popinit = as.numeric(BLASTtbmod::blantyre$population),
    patch_dims = patch_dims,
    HIV_dims = HIV_dims,
    age_dims = age_dims,
    age_rate = c(1 / 15, 1 / 35, 1e-6),
    agefracs = age.frax, # initial age fractions
    ageMids = c(15 / 2, (15 + 50) / 2, 60),
    initD = Dinit, #initial state
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
    ACFhaz0 = matrix(0.0,nrow=patch_dims,ncol=sim_length), #asymp ACF haz
    ACFhaz1 = matrix(0.0,nrow=patch_dims,ncol=sim_length)  #symp ACF haz
  )
  return(parms)
}

