
#############################
## Parameters
#############################

## Simulation parameters
tt[] <- user()
sim_length <- user()
dt <- user()

##################################
## Set initial conditions
ari0 <- user()
popinit[] <- user()
ageMids[] <- user()
agefracs[] <- user()

## sampling initializations
## popinit_byage[,] <- rbinom(floor(popinit[i]),agefracs[j]/(sum(agefracs)+1e-10)) #division just in case
## popinit_byage[1:patch_dims, 1] <- rbinom(floor(popinit[i]), agefracs[1]/(sum(agefracs)+1e-10))
## popinit_byage[1:patch_dims, 2:age_dims] <- rbinom(floor(popinit[i]) - sum(popinit_byage[i, 1:(j - 1)]),agefracs[j]/(sum(agefracs[j:patch_dims])+1e-10)) # NO

popinit_byage[, ] <- floor(popinit[i] * agefracs[j] / (sum(agefracs) + 1e-10)) # division just in case
## popinit_byage[,] <- rbinom(floor(popinit[i]),agefracs[j]/(sum(agefracs)+1e-10)) #division just in case

## TB initial state
## Initial ratio of TBI to disease states
initD[, ] <- user() # NOTE now taken as input: note mix of D,SC,Tr,R below
initPrev[,] <- exp( -ari0*ageMids[j] ) * (1-5*initD[i,j]/2)             #non-LTBI=U
initLL[, ] <- (1.0 - exp( -ari0*ageMids[j] )) * exp(-2 * ari0) * (1-5*initD[i,j]/2)             #non-LTBI=U
initF[, ] <- (1.0 - exp( -ari0*ageMids[j] )) * (1.0 - exp(-2 * ari0)) * (1-5*initD[i,j]/2)             #non-LTBI=U

## ## fraction TBI initially R or LR: as ~2 year's of FOI/TBI
## initF[,] <- if(initPrev[i,j]>tol) (1-exp(-2*ari0)) / initPrev[i,j] else 0

## safety: should be 1
initDenom[, ] <- initPrev[i,j] + initF[i, j] + initLL[i, j] + 5 * initD[i, j] / 2

## print("initD: {initD[,,]}",when=sum(initD)>0)
## print("initPrev: {initPrev[,,]}",when=sum(initPrev)>0)
## print("initLL: {initLL[,,]}",when=sum(initLL)>0)
## print("initF: {initF[,,]}",when=sum(initF)>0)


## tbi_U[, ] <-  (initPrev[i, j]) / (initDenom[i, j]+tol)
## tbi_LR[,] <- (initF[i,j])/(initDenom[i,j]+tol)
## tbi_LL[,] <-  (initLL[i,j])/(initDenom[i,j]+tol)
## tbi_D[,] <-  (initD[i,j]/2)/(initDenom[i,j]+tol)
## tbi_SC[,] <-  (initD[i,j]/2)/(initDenom[i,j]+tol)
## tbi_Tr[,] <-  (initD[i,j]/2)/(initDenom[i,j]+tol)
## tbi_R[,] <- (initD[i,j])/(initDenom[i,j]+tol)

tbi_U[, ] <- if (initDenom[i, j] > tol) (initPrev[i, j]) / initDenom[i, j] else 0
tbi_LR[,] <- if(initDenom[i,j] > tol) (initF[i,j])/initDenom[i,j] else 0
tbi_LL[,] <- if(initDenom[i,j] > tol) (initLL[i,j])/initDenom[i,j] else 0
tbi_D[,] <- if(initDenom[i,j] > tol) (initD[i,j]/2)/initDenom[i,j] else 0
tbi_SC[,] <- if(initDenom[i,j] > tol) (initD[i,j]/2)/initDenom[i,j] else 0
tbi_Tr[,] <- if(initDenom[i,j] > tol) (initD[i,j]/2)/initDenom[i,j] else 0
tbi_R[,] <- if(initDenom[i,j] > tol) (initD[i,j])/initDenom[i,j] else 0


## NOTE total stochastic even given above: uncommented below & commented 3rd block below
## init_U[, ] <- rbinom(popinit_byage[i, j], tbi_U[i, j])
## init_LR[,] <- rbinom(popinit_byage[i,j],tbi_LR[i,j])
## init_LL[,] <- rbinom(popinit_byage[i,j],tbi_LL[i,j])
## init_D[,] <- rbinom(popinit_byage[i,j],tbi_D[i,j])
## init_SC[,] <- rbinom(popinit_byage[i,j],tbi_SC[i,j])
## init_Tr[,] <- rbinom(popinit_byage[i,j],tbi_Tr[i,j])
## init_R[,] <- rbinom(popinit_byage[i,j],tbi_R[i,j])


## init_U[, ] <- rbinom(popinit_byage[i, j], tbi_U[i, j])
## init_LR[, ] <- rbinom(popinit_byage[i, j] - init_U[i, j], tbi_LR[i, j] / (1 - tbi_U[i, j]))
## init_LL[, ] <- rbinom(popinit_byage[i, j] - init_U[i, j] - init_LR[i, j], tbi_LL[i, j] / (1 - tbi_U[i, j] - tbi_LR[i, j]))
## init_D[, ] <- rbinom(popinit_byage[i, j] - init_U[i, j] - init_LR[i, j] - init_LL[i, j], tbi_D[i, j] / (1 - tbi_U[i, j] - tbi_LR[i, j] - tbi_LL[i, j]))
## init_SC[, ] <- rbinom(popinit_byage[i, j] - init_U[i, j] - init_LR[i, j] - init_LL[i, j] - init_D[i, j], tbi_SC[i, j] / (1 - tbi_U[i, j] - tbi_LR[i, j] - tbi_LL[i, j]-tbi_D[i, j]))
## init_Tr[, ] <- rbinom(popinit_byage[i, j] - init_U[i, j] - init_LR[i, j] - init_LL[i, j] - init_D[i, j] - init_SC[i, j], tbi_Tr[i, j] / (1 - tbi_U[i, j] - tbi_LR[i, j] - tbi_LL[i, j] - tbi_D[i, j] - tbi_SC[i, j]))
## init_R[, ] <- rbinom(popinit_byage[i, j] - init_U[i, j] - init_LR[i, j] - init_LL[i, j] - init_D[i, j] - init_SC[i, j] - init_Tr[i, j], tbi_R[i, j] / (1 - tbi_U[i, j] - tbi_LR[i, j] - tbi_LL[i, j] - tbi_D[i, j] - tbi_SC[i, j]))

init_U[, ] <- round(popinit_byage[i, j] * tbi_U[i, j])
init_LR[,] <- round(popinit_byage[i,j] * tbi_LR[i,j])
init_LL[,] <- round(popinit_byage[i,j] * tbi_LL[i,j])
init_D[,] <- round(popinit_byage[i,j] * tbi_D[i,j])
init_SC[,] <- round(popinit_byage[i,j] * tbi_SC[i,j])
init_Tr[,] <- round(popinit_byage[i,j] * tbi_Tr[i,j])
init_R[,] <- round(popinit_byage[i,j] * tbi_R[i,j])


## TODO this will need changing to account for HIV split:
## less burnin than plannedw

initial(U[,,1]) <- init_U[i,j]
initial(U[,,2:3]) <- 0
initial(LR[,,1]) <- init_LR[i,j]
initial(LR[,,2:3]) <- 0
initial(LL[,,1]) <- init_LL[i,j]
initial(LL[,,2:3]) <- 0

# # For now: Split initial population 80:20 clinical:subclinical
initial(D[,,1]) <- init_D[i,j]
initial(D[,,2:3]) <- 0
initial(SC[,,1]) <- init_SC[i,j]
initial(SC[,,2:3]) <- 0

initial(Tr[,,1]) <- init_Tr[i,j]
initial(Tr[,,2:3]) <- 0
initial(R[,,1]) <- init_R[i,j]
initial(R[,,2:3]) <- 0


########################################
# Initial conditions for vars that were previously output,
# now update for compatibility with odin.dust

initial(N[,,1]) <- init_U[i,j] + init_LR[i,j] + init_LL[i,j] +
  init_D[i,j] + init_SC[i,j] + init_Tr[i,j] + init_R[i,j]
initial(N[,,2:3]) <- 0

# Not sure how to set initial vals for these as they depend on
# dynamics. 0 for now:
initial(TB_deaths[,,]) <- 0
initial(notifrate[,,]) <- 0
initial(notes[,,]) <- 0
initial(bg_deaths[,,]) <- 0
initial(incidence[,,]) <- 0
initial(tot_incidence) <- 0


## ########################
## test variables
initial(beta_test) <- beta

print("--- initial state done:{tot_incidence} ---")


###########################
## For population/compartment dynamics

v <- user()      # v - infection modifier (protection conferred from latent or previous infection?)
pDf <- user()    # pDf - progression rate early to prevalent infectious
pDs <- user()    # pDs - progression rate late to prevalent infectious
pDr <- user()    # pDr - Relapse rate from recovered to active disease
pLL <- user()    # pLL - progression rate early to late infection (stablisation)
tfr <- user()    # tfr - fatality rate on treatment
cfr <- user()    # cfr - case fatality rate (untreated)
dur <- user()    # dur - Duration of active disease
cdr <- user()    # cdr - case detection rate
cdr_SC <- user() # case detection rate for subclinical disease
beta <- user()   # beta - effective contact rate
t_dur <- user(0.5)               # t_dur - Duration of treatment (6 months)
TBd_rate <- cfr/dur              # TB death rate
slfcr_rate <- (1-cfr)/dur        # Selfcure rate
## m_in <- user(0.0)               # m_in - Migration rate in
age_rate[] <-user()
regress_rate <- user()
progress_rate <- user()

## === ACF
## Detection rates
dtct_rate <- cdr/(dur*(1-cdr))   # Detection rate
dtct_rate_SC <- cdr_SC / (dur * (1 - cdr_SC))
## see: https://mrc-ide.github.io/odin.dust/articles/porting.html
dtct_ratem[] <- if(as.integer(step) >= sim_length) dtct_rate + ACFhaz1[i,sim_length] else dtct_rate + ACFhaz1[i,step+1]  #symptomatic
dtct_rate_SCm[] <- if(as.integer(step) >= sim_length) dtct_rate_SC + ACFhaz0[i,sim_length] else dtct_rate_SC + ACFhaz0[i,step+1] #asymptomatic
dim(dtct_ratem) <- c(patch_dims)
dim(dtct_rate_SCm) <- c(patch_dims)
ACFhaz0[,] <- user() #asymptomatic
dim(ACFhaz0) <- c(patch_dims,sim_length)
ACFhaz1[,] <- user() #symptomatic
dim(ACFhaz1) <- c(patch_dims,sim_length)

########################
# Assigning interpolated inputs

# Assign birthrate
births_int[] <- user()
births_t <- if(as.integer(step) < length( births_int )) 
  births_int[step+1] else births_int[length(births_int)]
birth_rate_yr <- births_t/1000

# Assign HIV rate
# Time & age-dependent
HIV_int[] <- user()
HIV_t <- if(as.integer(step) < sim_length )
  HIV_int[step+1] else HIV_int[sim_length]
HIV_rate_yr[1] <- 0
HIV_rate_yr[2] <- HIV_t/1000
HIV_rate_yr[3] <- 0

# Assign ART initiation rate
# Time & age dependent
ART_int[] <- user()
ART_t <- if(as.integer(step) < sim_length )
  ART_int[step+1] else ART_int[sim_length]
ART_rate_yr[1] <- 0
ART_rate_yr[2:3] <- ART_t

## TB IRR placeholder
TB_HIV_mod[] <- user() #NOTE this is applied as 
Hirr[] <- user()       #NOTE IRR applied to progression

# Assign background mortality
mu_noHIV_int[,] <- user()
mu_noHIV_t[] <- if(as.integer(step) < sim_length )
  mu_noHIV_int[i,step+1] else mu_noHIV_int[i,sim_length]
mu_HIV_int[,] <- user()
mu_HIV_t[] <- if(as.integer(step) < sim_length )
  mu_HIV_int[i,step+1] else mu_HIV_int[i,sim_length]
mu_ART_int[,] <- user()
mu_ART_t[] <- if(as.integer(step) < sim_length )
  mu_ART_int[i,step+1] else mu_ART_int[i,sim_length]

## in migration
m_in_int[, ] <- user()
m_in_t[] <- if (as.integer(step) < sim_length) m_in_int[i, step + 1] else m_in_int[i, sim_length]


#################################
## Other

## Seasonality
beta_seas <- beta*( 1+epsi*( sin( 2*Pi*tt[step]/12 )))
Pi <- user()
epsi <- user()

## Mixing/infectivity
## Mij = contact rate from patch i to patch j
## foitemp[1:patch_dims, 1:patch_dims] <- beta_seas * MM[i, j] * InfPrev[j]
foitemp[1:patch_dims, 1:patch_dims] <- beta_seas * MM[i, j] * (sum(D[j, , ]) + relinf * sum(SC[j, , ])) / (sum(N[j, , ]) + 1e-15) #
## force of infection on patch i from patch j: contact rate x prev in patch j
foi[1:patch_dims] <- sum( foitemp[i,] ) #FOI on each patch i: sum over FOIs from each patch
MM[,] <- user()           # Mixing matrix
relinf <- user(1)         # Relative infectiousness of SC relative to CD
IRR[] <- user()

## flux variables
initial(cum_inf_flux[, ]) <- 0 # initialize cumulative fluxes
initial(cum_inf_ByPatch[ ]) <- 0 # initialize cumulative fluxes
initial(cum_note_flux[, ]) <- 0 # initialize cumulative fluxes
initial(Tijk[, , ]) <- 0
initial(Sij[, ]) <- 0
initial(PrevByPatch[]) <- 0


## patch trackers
initial(incidence_bypatch[]) <- 0
initial(prevalence_bypatch[]) <- 0


#########################
# Dimensions
#########################

# Patches = i dimension
# Age = j dimension
# HIV status = k dimension
patch_dims <- user()
age_dims <- user()
HIV_dims <- user()


## in migration event dimensions
## event counts
dim(U_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(LR_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(LL_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(SC_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(D_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(Tr_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(R_inmigr) <- c(patch_dims, age_dims, HIV_dims)
## prob counts
dim(pU_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(pLR_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(pLL_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(pSC_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(pD_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(pTr_inmigr) <- c(patch_dims, age_dims, HIV_dims)
dim(pR_inmigr) <- c(patch_dims, age_dims, HIV_dims)

## events
dim(NeventsU) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsU1) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsU2) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsU3) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsU4) <- c( patch_dims, age_dims, HIV_dims )

dim(NeventsLR) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsLR1) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsLR2) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsLR3) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsLR4) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsLR5) <- c( patch_dims, age_dims, HIV_dims )

dim(NeventsLL) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsLL1) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsLL2) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsLL3) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsLL4) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsLL5) <- c( patch_dims, age_dims, HIV_dims )

dim(NeventsSC) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsSC1) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsSC2) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsSC3) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsSC4) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsSC5) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsSC6) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsSC7) <- c( patch_dims, age_dims, HIV_dims )

dim(NeventsD) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsD1) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsD2) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsD3) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsD4) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsD5) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsD6) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsD7) <- c( patch_dims, age_dims, HIV_dims )

dim(NeventsT) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsT1) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsT2) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsT3) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsT4) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsT5) <- c( patch_dims, age_dims, HIV_dims )

dim(NeventsR) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsR1) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsR2) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsR3) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsR4) <- c( patch_dims, age_dims, HIV_dims )
dim(NeventsR5) <- c( patch_dims, age_dims, HIV_dims )

## Dims relating to U compartment (uninfected)
dim( U ) <- c( patch_dims, age_dims, HIV_dims )
dim( init_U ) <- c( patch_dims, age_dims )
dim( Udeaths) <- c( patch_dims, age_dims, HIV_dims )
dim( Uinfs) <- c( patch_dims, age_dims, HIV_dims )
## dim( m_in_U ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_in_U ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_out_U ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_in_U ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_out_U ) <- c( patch_dims, age_dims, HIV_dims )

# # Upgrade
dim( rate_U ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_Uinf ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_Uage ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_UHIV ) <- c( patch_dims, age_dims, HIV_dims )

# ## Dims relating to LR compartment (infected early)
dim( LR ) <- c( patch_dims, age_dims, HIV_dims )
dim( init_LR ) <- c( patch_dims, age_dims )
## dim( m_in_LR ) <- c( patch_dims, age_dims, HIV_dims )
dim( LRdeaths ) <- c( patch_dims, age_dims, HIV_dims )
dim( progFast ) <- c( patch_dims, age_dims, HIV_dims )
dim( rate_LR ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_progFast ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_stabilise ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_LRage ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_LRHIV ) <- c( patch_dims, age_dims, HIV_dims )
dim( stabilisations ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_in_LR ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_out_LR ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_in_LR ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_out_LR ) <- c( patch_dims, age_dims, HIV_dims )

## Dims relating to LL compartment (infected late)
dim( LL ) <- c( patch_dims, age_dims, HIV_dims )
dim( init_LL ) <- c( patch_dims, age_dims )
dim( initLL ) <- c( patch_dims, age_dims ) #this one is probs
dim( rate_LL ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_LLinfs ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_progSlow ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_LLage ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_LLHIV ) <- c( patch_dims, age_dims, HIV_dims )
dim( LLdeaths ) <- c( patch_dims, age_dims, HIV_dims )
dim( LLinfs ) <- c( patch_dims, age_dims, HIV_dims )
dim( progSlow ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_in_LL ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_out_LL ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_in_LL ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_out_LL ) <- c( patch_dims, age_dims, HIV_dims )

## Dims relating to D compartment (active disease)
dim( D ) <- c( patch_dims, age_dims, HIV_dims )
dim( init_D ) <- c( patch_dims, age_dims )
dim( Ddeaths ) <- c( patch_dims, age_dims, HIV_dims )
dim( TBdeaths ) <- c( patch_dims, age_dims, HIV_dims )
dim( detection ) <- c( patch_dims, age_dims, HIV_dims )
dim( selfcure ) <- c( patch_dims, age_dims, HIV_dims )
dim( regress ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_in_D ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_out_D ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_in_D ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_out_D ) <- c( patch_dims, age_dims, HIV_dims )
dim( rate_D ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_Dage ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_DHIV ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_detect ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_selfcr ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_regress ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_TBdeaths ) <- c( patch_dims, age_dims, HIV_dims )

## Dims relating to SC compartment (subclinical disease)
dim( SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( init_SC ) <- c( patch_dims, age_dims )
dim( SCdeaths ) <- c( patch_dims, age_dims, HIV_dims )
dim( TBdeaths_SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( detection_SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( selfcure_SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( progress ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_in_SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_out_SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( incidence ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_in_SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_out_SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( rate_SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_SCage ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_SCHIV ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_detect_SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_selfcr_SC ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_progress ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_TBdeaths_SC ) <- c( patch_dims, age_dims, HIV_dims )


## Dims relating to Tr compartment (treated)
dim( Tr ) <- c( patch_dims, age_dims, HIV_dims )
dim( init_Tr ) <- c( patch_dims, age_dims )
dim( Trdeaths ) <- c( patch_dims, age_dims, HIV_dims )
dim( Trsuccess ) <- c( patch_dims, age_dims, HIV_dims )
dim( Trxdeaths ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_in_Tr ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_out_Tr ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_in_Tr ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_out_Tr ) <- c( patch_dims, age_dims, HIV_dims )
dim( rate_Tr ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_Trage ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_TrHIV ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_Trsuccess ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_Trxdeaths ) <- c( patch_dims, age_dims, HIV_dims )


## Dims relating to R compartment (Recovered)
dim( R ) <- c( patch_dims, age_dims, HIV_dims )
dim( init_R ) <- c( patch_dims, age_dims )
dim( rate_R ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_Rinfs ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_relapse ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_Rage ) <- c( patch_dims, age_dims, HIV_dims )
dim( p_RHIV ) <- c( patch_dims, age_dims, HIV_dims )
dim( Rdeaths ) <- c( patch_dims, age_dims, HIV_dims )
dim( Rinfs ) <- c( patch_dims, age_dims, HIV_dims )
dim( relapse ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_in_R ) <- c( patch_dims, age_dims, HIV_dims )
dim( age_out_R ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_in_R ) <- c( patch_dims, age_dims, HIV_dims )
dim( HIV_out_R ) <- c( patch_dims, age_dims, HIV_dims )


## Additional dims
dim( age_rate ) <- age_dims
dim( births_int ) <- user()
dim( HIV_int ) <- user()
dim( HIV_rate_yr ) <- age_dims
dim( ART_int ) <- user()
dim( ART_rate_yr ) <- age_dims
dim(TB_HIV_mod) <- HIV_dims
dim(Hirr) <- HIV_dims
dim( mu_noHIV_int ) <- user()
dim( mu_noHIV_t ) <- c( age_dims )
dim( mu_HIV_int ) <- user()
dim( mu_HIV_t ) <- c( age_dims )
dim( mu_ART_int ) <- user()
dim(mu_ART_t) <- c(age_dims)
dim(m_in_int) <- user()
dim(m_in_t) <- c(age_dims)
dim( N ) <- c( patch_dims, age_dims, HIV_dims )
dim( births ) <- c( patch_dims, age_dims, HIV_dims )
dim( notifrate ) <- c( patch_dims, age_dims, HIV_dims )
dim( notes ) <- c( patch_dims, age_dims, HIV_dims )
dim( TB_deaths ) <- c( patch_dims, age_dims, HIV_dims )
dim( bg_deaths ) <- c( patch_dims, age_dims, HIV_dims )
dim( MM ) <- c( patch_dims, patch_dims )
dim( foitemp ) <- c( patch_dims, patch_dims )
dim( foi ) <- patch_dims
dim(popinit) <- patch_dims
dim(birthparm) <- patch_dims
dim( IRR ) <- patch_dims
dim( tt ) <- sim_length
dim( popinit_byage ) <- c( patch_dims, age_dims )
dim( initPrev ) <- c( patch_dims, age_dims )
dim( initF ) <- c( patch_dims, age_dims )
dim( initD ) <- c( patch_dims, age_dims )
# dim( ari0 ) <- patch_dims
dim( ageMids ) <- age_dims
dim( agefracs ) <- age_dims
dim(initDenom) <- c(patch_dims, age_dims)
dim(tbi_U) <- c(patch_dims, age_dims)
dim(tbi_LR) <- c( patch_dims, age_dims )
dim(tbi_LL) <- c( patch_dims, age_dims )
dim(tbi_D) <- c( patch_dims, age_dims )
dim(tbi_SC) <- c( patch_dims, age_dims )
dim(tbi_Tr) <- c( patch_dims, age_dims )
dim(tbi_R) <- c( patch_dims, age_dims )

## ## patch level aggregate indicators
## dim(pnotifrate) <- c(patch_dims)
## dim(pnotes) <- c(patch_dims)
dim(ppop) <- c(patch_dims)

## flux variables
dim(InfsByPatchPatch) <- c(patch_dims, patch_dims)
dim(NotesByPatchPatch) <- c(patch_dims, patch_dims)
dim(InfsByPatch) <- c(patch_dims)
dim(PrevByPatch) <- c(patch_dims)
dim(InfPrev) <- c(patch_dims)
dim(NotesByPatch) <- c(patch_dims)
dim(cum_inf_ByPatch) <- c(patch_dims)
dim(cum_inf_flux) <- c(patch_dims, patch_dims)
dim(cum_note_flux) <- c(patch_dims, patch_dims)
dim(A) <- 6
dim(zk) <- 6
dim(Tijk) <- c(patch_dims, patch_dims, 6)
dim(Sij) <- c(patch_dims, patch_dims)
dim(ellij) <- c(patch_dims, patch_dims)
dim(incidence_bypatch) <- patch_dims
dim(prevalence_bypatch) <- patch_dims


##########
## Additional outputs to record
##########

## ijk = patch, age, HIV
##########################
## Populations
##########################

# S (U) - susceptible (uninfected)
# E (LR) - early infection
# L (LL) - late infection
# I (D) - prevalent + infectious
# (SCD) - Subclinical disease - infectious
# T (Tr) - receiving treatment, assumed not infections
# R (R) - recovered at risk of relapse


#######################
## Notes on entry/exit events
######################

## BIRTHS
# Birth rate in patch i as % of total population in patch i
# FOR NOW assume all births in non-HIV index

## AGEING
# No ageing into first index; 
# (no ageing out of last index handled by age_rate array)

## HIV INFECTION
# No movement into first index
# Moving out of last index handled by HIV_rate array

## Overall exit rate (each compartment)
# This is total event rate (hazard)
# Exit events assigned according to relative probability 
# age & HIV now calculated first as is prob=0 for final dimension, so cannot be 'leftover' events


##############################
## Updates
##############################
tol <- 1e-10 # a safety: tolerance for denominator size

##########################################
####### 1: Uninfected (susceptible) ######
##########################################


###### UNINFECTED ENTRY EVENTS #######
### a) BIRTHS + (net?) MIGRATION
# Birth rate downloaded from https://population.un.org/wpp/

mxp <- 1000000 #patch max population
birthparm[1:patch_dims] <- if (birth_rate_yr * ppop[i] > 0 && ppop[i] < mxp) birth_rate_yr * ppop[i] * dt else 0
births[1:patch_dims, 1, 1] <- rpois(birthparm[i])
births[1:patch_dims, 2:age_dims, 2:HIV_dims] <- 0
## m_in_U[,,] <- rpois( m_in*U[i,j,k]*dt )

### b) AGEING
age_in_U[,2:age_dims,] <- age_out_U[i,j-1,k]
age_in_U[,1,]<-0

### c) HIV infection
HIV_in_U[,,2:HIV_dims] <- HIV_out_U[i,j,k-1]
HIV_in_U[,,1]<-0


###### UNINFECTED EXIT EVENTS #####
## Total rate & events
# print("U: {min(U[,,])}",when=min(U[,,])<0)


# a: Movement out of no-HIV compartments; age dependent; age 0-14 and 50+ Zero
rate_U[, , 1] <- m_in_t[j] + mu_noHIV_t[j] + foi[i] * TB_HIV_mod[1] + age_rate[j] + HIV_rate_yr[j]
rate_U[, , 2] <- m_in_t[j] + mu_HIV_t[j] + foi[i] * TB_HIV_mod[2] + age_rate[j] + ART_rate_yr[j]
rate_U[, , 3] <- m_in_t[j] + mu_ART_t[j] + foi[i] * TB_HIV_mod[3] + age_rate[j]
NeventsU[,,] <- if(U[i,j,k]>0) rbinom( U[i,j,k], 1-exp( -rate_U[i,j,k]*dt )) else 0

## NOTE denominator safety introduced

## a) AGEING
p_Uage[,,] <- if(rate_U[i,j,k] > tol) age_rate[j]/rate_U[i,j,k] else 0
age_out_U[,,] <- if( p_Uage[i,j,k]>0 && NeventsU[i, j, k]>0) rbinom( NeventsU[i,j,k], p_Uage[i,j,k] ) else 0#
NeventsU1[,,] <- NeventsU[i,j,k] - age_out_U[i,j,k]         #

## b) HIV INFECTION
## of remaining events (bar ageing), what proportion are HIV processes?
## NOTE safety
p_UHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + foi[i] * TB_HIV_mod[1] + HIV_rate_yr[j] > tol) HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + foi[i] * TB_HIV_mod[1] + HIV_rate_yr[j]) else 0
p_UHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + foi[i] * TB_HIV_mod[2] + ART_rate_yr[j] > tol) ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + foi[i] * TB_HIV_mod[2] + ART_rate_yr[j]) else 0
p_UHIV[,,3] <- 0
HIV_out_U[,,] <- if( p_UHIV[i,j,k]>0 && NeventsU1[i, j, k]>0) rbinom( NeventsU1[i,j,k], p_UHIV[i,j,k]) else 0#
NeventsU2[,,] <- NeventsU1[i,j,k] - HIV_out_U[i,j,k]         # update

## c) TB INFECTION
p_Uinf[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + foi[i] * TB_HIV_mod[1] > tol) foi[i] * TB_HIV_mod[1] / (m_in_t[j] + mu_noHIV_t[j] + foi[i] * TB_HIV_mod[1]) else 0
p_Uinf[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + foi[i] * TB_HIV_mod[2] > tol) foi[i] * TB_HIV_mod[2] / (m_in_t[j] + mu_HIV_t[j] + foi[i] * TB_HIV_mod[2]) else 0
p_Uinf[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + foi[i] * TB_HIV_mod[3] > tol) foi[i] * TB_HIV_mod[3] / (m_in_t[j] + mu_ART_t[j] + foi[i] * TB_HIV_mod[3]) else 0
Uinfs[,,] <- if(p_Uinf[i,j,k] >0 && NeventsU2[i, j, k]>0) rbinom( NeventsU2[i,j,k], p_Uinf[i,j,k]) else 0 #
NeventsU3[,,] <- NeventsU2[i,j,k] - Uinfs[i,j,k]         # update

## d) IN MIGRATION
pU_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / m_in_t[j] + mu_noHIV_t[j] else 0
pU_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pU_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
U_inmigr[, , ] <- if(pU_inmigr[i, j, k]>0 && NeventsU3[i, j, k] > 0) rbinom(NeventsU3[i, j, k],pU_inmigr[i, j, k]) else 0
NeventsU4[,,] <- NeventsU3[i,j,k] - U_inmigr[i,j,k]         # update

## e) DEATH
Udeaths[, , ] <- NeventsU4[i, j, k] #


## actual updates
update(U[,,]) <- U[i,j,k] + births[i,j,k] + age_in_U[i,j,k] + HIV_in_U[i,j,k] - ## m_in_U[i,j,k] -
  age_out_U[i, j, k] - HIV_out_U[i, j, k] - Udeaths[i, j, k] - Uinfs[i, j, k]+
  U_inmigr[i, j, k]


## ---------------------- extra calculations associated with infection fluxes
## number of infections in patch i this step
InfsByPatch[1:patch_dims] <- sum(Uinfs[i, , ]) + sum(LLinfs[i,,]) + sum(Rinfs[i,,])

## correct looped binom implementation of multinomial
InfsByPatchPatch[1:patch_dims, 1] <- if(InfsByPatch[i]>0) rbinom(InfsByPatch[i], foitemp[i, 1] / (foi[i] + tol)) else 0
InfsByPatchPatch[1:patch_dims, 2:patch_dims] <- if(InfsByPatch[i] - sum(InfsByPatchPatch[i, 1:(j - 1)]) > 0) rbinom(InfsByPatch[i] - sum(InfsByPatchPatch[i, 1:(j - 1)]), foitemp[i, j] / (sum(foitemp[i, j:patch_dims]) + tol)) else 0

## ## debug/test
## InfsByPatchPatch[1:patch_dims, 1:patch_dims] <- round(InfsByPatch[i] * foitemp[i, j] / (foi[i] + tol)) #debug


update(cum_inf_ByPatch[]) <- cum_inf_ByPatch[i] + InfsByPatch[i]
update(cum_inf_flux[, ]) <- cum_inf_flux[i, j] + InfsByPatchPatch[i, j]
InfPrev[1:patch_dims] <- (sum(D[i, , ]) + relinf * sum(SC[i, , ])) / (sum(N[i, , ]) + 1e-15)
update(PrevByPatch[1:patch_dims]) <- InfPrev[i]

NotesByPatch[] <- sum(notes[i,,]) #summing over age/HIV
NotesByPatchPatch[1:patch_dims, 1] <- if(NotesByPatch[i] > 0) rbinom(NotesByPatch[i], ellij[i, 1] / (sum(ellij[i,1:patch_dims]) + tol)) else 0
NotesByPatchPatch[1:patch_dims, 2:patch_dims] <- if(NotesByPatch[i] - sum(NotesByPatchPatch[i, 1:(j - 1)]) > 0) rbinom(NotesByPatch[i] - sum(NotesByPatchPatch[i, 1:(j - 1)]), ellij[i, j] / (sum(ellij[i, j:patch_dims]) + tol)) else 0

## ## test/debug
## NotesByPatchPatch[1:patch_dims, 1:patch_dims] <-  round(NotesByPatch[i] * ellij[i, j] / (sum(ellij[i, 1:patch_dims]) + tol)) #debug

update(cum_note_flux[, ]) <- cum_note_flux[i, j] + NotesByPatchPatch[i, j]


## Notes
## NotesByPatchPatch[1:patch_dims, 1] <- rbinom(NotesByPatch[i], foitemp[i, 1] / (foi[i] + tol))
## NotesByPatchPatch[1:patch_dims, 2:patch_dims] <- rbinom(NotesByPatch[i] - sum(NotesByPatchPatch[i, 1:(j - 1)]), foitemp[i, j] / (sum(foitemp[i, j:patch_dims]) + tol)) # OK
## update(cum_note_flux[, ]) <- cum_note_flux[i, j] + NotesByPatchPatch[i, j]
## update(cum_note_flux[, ]) <- cum_note_flux[i, j] + foitemp[i, j]
## update(cum_note_flux[, ]) <- cum_note_flux[i, j] + ellij[i, j] #OK


##############################
##### 2: Infected early ######
##############################

###### INFECTED EARLY: ENTRY EVENTS #####
rate_LR[, , 1] <- mu_noHIV_t[j] + pDf * IRR[i] * Hirr[1] + pLL + age_rate[j] + HIV_rate_yr[j]
rate_LR[, , 2] <- mu_HIV_t[j] + pDf * IRR[i] * Hirr[2] + pLL + age_rate[j] + ART_rate_yr[j]
rate_LR[, , 3] <- mu_ART_t[j] + pDf * IRR[i] * Hirr[3] + pLL + age_rate[j]
NeventsLR[, , ] <- if(LR[i, j, k]>0) rbinom(LR[i, j, k], 1 - exp(-rate_LR[i, j, k] * dt)) else 0

## b) AGEING
age_in_LR[, 2:age_dims, ] <- age_out_LR[i, j - 1, k]
age_in_LR[, 1, ] <- 0

## c) HIV
HIV_in_LR[, , 2:HIV_dims] <- HIV_out_LR[i, j, k - 1]
HIV_in_LR[, , 1] <- 0

## a) AGEING
p_LRage[,,] <- if(rate_LR[i,j,k] > tol) age_rate[j]/rate_LR[i,j,k] else 0
age_out_LR[,,] <- if(p_LRage[i,j,k] >0 && NeventsLR[i, j, k]>0) rbinom( NeventsLR[i,j,k], p_LRage[i,j,k] ) else 0
NeventsLR1[,,] <- NeventsLR[i,j,k] - age_out_LR[i,j,k]

## b) HIV INFECTION
p_LRHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + pDf * IRR[i] * Hirr[1] + pLL + HIV_rate_yr[j] > tol)
                    HIV_rate_yr[j] /
                      (m_in_t[j] + mu_noHIV_t[j] + pDf * IRR[i]* Hirr[1] + pLL + HIV_rate_yr[j]) else 0
p_LRHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + pDf * IRR[i] * Hirr[2] + pLL + ART_rate_yr[j] > tol)
                    ART_rate_yr[j] /
                      (m_in_t[j] + mu_HIV_t[j] + pDf * IRR[i] * Hirr[2] + pLL + ART_rate_yr[j]) else 0
p_LRHIV[,,3] <- 0
HIV_out_LR[,,] <- if(p_LRHIV[i,j,k] >0 && NeventsLR1[i, j, k]>0) rbinom( NeventsLR1[i,j,k], p_LRHIV[i,j,k]) else 0
NeventsLR2[,,] <- NeventsLR1[i,j,k] - HIV_out_LR[i,j,k]


## c) FAST PROGRESSION
p_progFast[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + pDf * IRR[i] * Hirr[1] + pLL > tol)
                       IRR[i] * pDf * Hirr[1] / (m_in_t[j] + mu_noHIV_t[j] + pDf * IRR[i] * Hirr[1] + pLL) else 0
p_progFast[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + pDf * IRR[i] * Hirr[2] + pLL > tol)
                       IRR[i] * pDf * Hirr[2] / (m_in_t[j] + mu_HIV_t[j] + pDf * IRR[i] * Hirr[2] + pLL) else 0
p_progFast[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + pDf * IRR[i] * Hirr[3] + pLL > tol)
                       IRR[i] * pDf * Hirr[3] / (m_in_t[j] + mu_ART_t[j] + pDf * IRR[i] * Hirr[3] + pLL) else 0
progFast[,,] <- if(p_progFast[i,j,k] >0 && NeventsLR2[i, j, k]>0) rbinom( NeventsLR2[i,j,k],p_progFast[i,j,k]) else 0
NeventsLR3[,,] <- NeventsLR2[i,j,k] - progFast[i,j,k]


## d) STABILISATION
## p_stabilise[,,] <- pLL/rate_LR[i,j,k]
p_stabilise[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + pLL > tol)
                        pLL / (m_in_t[j] + mu_noHIV_t[j] + pLL) else 0
p_stabilise[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + pLL > tol)
                        pLL / (m_in_t[j] + mu_HIV_t[j] + pLL) else 0
p_stabilise[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + pLL > tol)
                        pLL / (m_in_t[j] + mu_ART_t[j] + pLL) else 0
stabilisations[, , ] <- if(p_stabilise[i, j, k] >0 && NeventsLR3[i, j, k]>0) rbinom(NeventsLR3[i, j, k],p_stabilise[i, j, k]) else 0
NeventsLR4[,,] <- NeventsLR3[i,j,k] - stabilisations[i,j,k]


## e) IN MIGRATION
pLR_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pLR_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pLR_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
LR_inmigr[, , ] <- if(pLR_inmigr[i, j, k]>0 && NeventsLR4[i, j, k] > 0) rbinom(NeventsLR4[i, j, k],pLR_inmigr[i, j, k]) else 0
NeventsLR5[,,] <- NeventsLR4[i,j,k] - LR_inmigr[i,j,k]


## f) DEATH
LRdeaths[, , ] <- NeventsLR5[i, j, k]


## actual update
update(LR[, , ]) <- LR[i, j, k] + age_in_LR[i, j, k] + HIV_in_LR[i, j, k] + Uinfs[i, j, k] +
  LLinfs[i, j, k] + Rinfs[i, j, k] + LR_inmigr[i, j, k] -
  age_out_LR[i, j, k] - HIV_out_LR[i, j, k] - LRdeaths[i, j, k] -
  progFast[i, j, k] - stabilisations[i, j, k]


############################
##### 3: Infected late #####
############################

# print("LL: {min(LL[,,])}",when=min(LL[,,])<0)

###### INFECTED LATE: ENTRY EVENTS #####

## b) AGEING
age_in_LL[, 2:age_dims, ] <- age_out_LL[i, j - 1, k]
age_in_LL[, 1, ] <- 0


## c) HIV
HIV_in_LL[, , 2:HIV_dims] <- HIV_out_LL[i, j, k - 1]
HIV_in_LL[, , 1] <- 0


###### INFECTED LATE: EXIT EVENTS #####
rate_LL[, , 1] <- m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] +
  pDs * IRR[i] * Hirr[1] + age_rate[j] + HIV_rate_yr[j]
rate_LL[, , 2] <- m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] +
  pDs * IRR[i] * Hirr[2] + age_rate[j] + ART_rate_yr[j]
rate_LL[, , 3] <- m_in_t[j] + mu_ART_t[j] + v * foi[i] * TB_HIV_mod[3] + pDs * IRR[i] * Hirr[3] + age_rate[j]
NeventsLL[, , ] <- if(LL[i, j, k]>0) rbinom(LL[i, j, k], 1 - exp(-rate_LL[i, j, k] * dt)) else 0


## a) AGEING
p_LLage[, , ] <- if (rate_LL[i, j, k] > tol) age_rate[j] / rate_LL[i, j, k] else 0
age_out_LL[, , ] <- if(p_LLage[i, j, k] >0 && NeventsLL[i, j, k]>0) rbinom(NeventsLL[i, j, k],  p_LLage[i, j, k]) else 0
NeventsLL1[,,] <- NeventsLL[i,j,k] - age_out_LL[i,j,k] #update


## b) HIV INFECTION
p_LLHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] +
                      pDs * IRR[i] * Hirr[1] + HIV_rate_yr[j] > tol)
                    HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] +
                                      pDs * IRR[i] * Hirr[1] + HIV_rate_yr[j]) else 0
p_LLHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] +
                      pDs * IRR[i] * Hirr[2] + ART_rate_yr[j] > tol)
                    ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] +
                                      pDs * IRR[i] * Hirr[2] + ART_rate_yr[j]) else 0
p_LLHIV[, , 3] <- 0
HIV_out_LL[, , ] <- if(p_LLHIV[i, j, k] >0 && NeventsLL1[i, j, k]>0) rbinom(NeventsLL1[i, j, k], p_LLHIV[i, j, k]) else 0
NeventsLL2[,,] <- NeventsLL1[i,j,k] - HIV_out_LL[i,j,k] #update


## c) TB INFECTION
p_LLinfs[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] +
                       v * foi[i] * TB_HIV_mod[1] + pDs * IRR[i] * Hirr[1] > tol)
                     v * foi[i] * TB_HIV_mod[1] / (m_in_t[j] + mu_noHIV_t[j] +
                                                   v * foi[i] * TB_HIV_mod[1] + pDs * IRR[i] * Hirr[1]) else 0
p_LLinfs[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] +
                       v * foi[i] * TB_HIV_mod[2] + pDs * IRR[i] * Hirr[2] > tol)
                     v * foi[i] * TB_HIV_mod[2] / (m_in_t[j] + mu_HIV_t[j] +
                                                   v * foi[i] * TB_HIV_mod[2] + pDs * IRR[i] * Hirr[2]) else 0
p_LLinfs[, , 3] <- if (m_in_t[j] + mu_ART_t[j] +
                       v * foi[i] * TB_HIV_mod[3] + pDs * IRR[i] * Hirr[3] > tol)
                     v * foi[i] * TB_HIV_mod[3] / (m_in_t[j] + mu_ART_t[j] +
                                                   v * foi[i] * TB_HIV_mod[3] + pDs * IRR[i] * Hirr[3]) else 0
LLinfs[, , ] <- if(p_LLinfs[i, j, k] >0 && NeventsLL2[i, j, k]>0) rbinom(NeventsLL2[i, j, k],p_LLinfs[i, j, k]) else 0
NeventsLL3[,,] <- NeventsLL2[i,j,k] - LLinfs[i,j,k] #update


## d) SLOW PROGRESSION
p_progSlow[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + pDs * IRR[i] * Hirr[1] > tol)
                       IRR[i] * pDs * Hirr[1]/
                         (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] + pDs * IRR[i] * Hirr[1]) else 0
p_progSlow[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + pDs * IRR[i] * Hirr[2] > tol)
                       IRR[i] * pDs * Hirr[2] /
                         (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] + pDs * IRR[i] * Hirr[2]) else 0
p_progSlow[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + pDs * IRR[i] * Hirr[3] > tol)
                       IRR[i] * pDs * Hirr[3]/
                         (m_in_t[j] + mu_ART_t[j] + v * foi[i] * TB_HIV_mod[3] + pDs * IRR[i] * Hirr[3]) else 0
progSlow[, , ] <- if( p_progSlow[i, j, k] >0 && NeventsLL3[i, j, k]>0) rbinom(NeventsLL3[i, j, k],p_progSlow[i, j, k]) else 0
NeventsLL4[,,] <- NeventsLL3[i,j,k] - progSlow[i,j,k] #update


## e) IN MIGRATION
pLL_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pLL_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pLL_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
LL_inmigr[, , ] <- if(pLL_inmigr[i, j, k]>0 && NeventsLL4[i, j, k]>0) rbinom(NeventsLL4[i, j, k],pLL_inmigr[i, j, k]) else 0
NeventsLL5[,,] <- NeventsLL4[i,j,k] - LL_inmigr[i,j,k] #update


## f) DEATH
LLdeaths[, , ] <- NeventsLL5[i, j, k]

## actual update
update(LL[, , ]) <- LL[i, j, k] + age_in_LL[i, j, k] + HIV_in_LL[i, j, k] +
  stabilisations[i, j, k] + selfcure[i, j, k] + selfcure_SC[i, j, k] -
  age_out_LL[i, j, k] - HIV_out_LL[i, j, k] - LLdeaths[i, j, k] - LLinfs[i, j, k] - progSlow[i, j, k] + LL_inmigr[i, j, k]

##########################################
##### 4a: Subclinical active disease #####
##########################################

###### SC DISEASE: ENTRY EVENTS #####

## AGEING
age_in_SC[,2:age_dims,] <- age_out_SC[i,j-1,k]
age_in_SC[,1,]<-0

## HIV
HIV_in_SC[,,2:HIV_dims] <- HIV_out_SC[i,j,k-1]
HIV_in_SC[,,1]<-0

###### SC DISEASE: EXIT EVENTS #####
rate_SC[, , 1] <- m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate_SCm[i] + slfcr_rate +
  progress_rate + age_rate[j] + HIV_rate_yr[j]
rate_SC[, , 2] <- m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate_SCm[i] + slfcr_rate +
  progress_rate + age_rate[j] + ART_rate_yr[j]
rate_SC[, , 3] <- m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_rate_SCm[i] + slfcr_rate +
  progress_rate + age_rate[j]
NeventsSC[, , ] <- if(SC[i, j, k]>0) rbinom(SC[i, j, k], 1 - exp(-rate_SC[i, j, k] * dt)) else 0

## TBd_rate & slfcr_rate defined with active disease


## a) AGEING
p_SCage[, , ] <- if (rate_SC[i, j, k] > tol) age_rate[j] / rate_SC[i, j, k] else 0
age_out_SC[, , ] <- if(p_SCage[i, j, k] >0 && NeventsSC[i, j, k]>0) rbinom(NeventsSC[i, j, k], p_SCage[i, j, k]) else 0
NeventsSC1[,,] <- NeventsSC[i,j,k] - age_out_SC[i,j,k] #update


## b) HIV
p_SCHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate_SCm[i] +
                      slfcr_rate + progress_rate + HIV_rate_yr[j] > tol)
                    HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate_SCm[i] +
                                      slfcr_rate + progress_rate + HIV_rate_yr[j]) else 0
p_SCHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate_SCm[i] +
                      slfcr_rate + progress_rate + ART_rate_yr[j] > tol)
                    ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate_SCm[i] +
                                      slfcr_rate + progress_rate + ART_rate_yr[j]) else 0
p_SCHIV[, , 3] <- 0
HIV_out_SC[, , ] <- if(p_SCHIV[i, j, k] >0 && NeventsSC1[i, j, k]>0) rbinom(NeventsSC1[i, j, k], p_SCHIV[i, j, k]) else 0
NeventsSC2[,,] <- NeventsSC1[i,j,k] - HIV_out_SC[i,j,k] #update



## c) TB DETECTION
p_detect_SC[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate_SCm[i] +
                          slfcr_rate +progress_rate > tol)
                        dtct_rate_SCm[i] /(m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate_SCm[i] +
                                       slfcr_rate +progress_rate) else 0
p_detect_SC[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate_SCm[i] +
                          slfcr_rate +progress_rate > tol)
                        dtct_rate_SCm[i] /(m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate_SCm[i] +
                                       slfcr_rate +progress_rate) else 0
p_detect_SC[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_rate_SCm[i] +
                          slfcr_rate +progress_rate > tol)
                        dtct_rate_SCm[i] / (m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_rate_SCm[i] +
                                        slfcr_rate +progress_rate) else 0
detection_SC[, , ] <- if(p_detect_SC[i, j, k] >0 && NeventsSC2[i, j, k]>0) rbinom(NeventsSC2[i, j, k],p_detect_SC[i, j, k]) else 0
NeventsSC3[,,] <- NeventsSC2[i,j,k] - detection_SC[i,j,k] #update


## d) SELF-CURING
p_selfcr_SC[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + slfcr_rate + progress_rate > tol)
                        slfcr_rate /
                          (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + slfcr_rate + progress_rate) else 0
p_selfcr_SC[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + slfcr_rate + progress_rate > tol)
                        slfcr_rate /
                          (m_in_t[j] + mu_HIV_t[j] + TBd_rate + slfcr_rate + progress_rate) else 0
p_selfcr_SC[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + slfcr_rate + progress_rate > tol)
                        slfcr_rate /
                          (m_in_t[j] + mu_ART_t[j] + TBd_rate + slfcr_rate + progress_rate) else 0
selfcure_SC[, , ] <- if(p_selfcr_SC[i, j, k] >0 && NeventsSC3[i, j, k]>0) rbinom(NeventsSC3[i, j, k],p_selfcr_SC[i, j, k]) else 0
NeventsSC4[,,] <- NeventsSC3[i,j,k] - selfcure_SC[i,j,k] #update


## e) PROGRESSION to active disease
p_progress[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + progress_rate > tol)
                       progress_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + progress_rate) else 0
p_progress[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + progress_rate > tol)
                       progress_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + progress_rate) else 0
p_progress[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + progress_rate > tol)
                       progress_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate + progress_rate) else 0
progress[, , ] <- if(p_progress[i, j, k] >0 && NeventsSC4[i, j, k]>0) rbinom(NeventsSC4[i, j, k],p_progress[i, j, k]) else 0
NeventsSC5[,,] <- NeventsSC4[i,j,k] - progress[i,j,k] #update


## f) TB DEATH
p_TBdeaths_SC[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate > tol)
                          TBd_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate) else 0
p_TBdeaths_SC[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate > tol)
                          TBd_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate) else 0
p_TBdeaths_SC[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate > tol)
                          TBd_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate) else 0
TBdeaths_SC[, , ] <- if(p_TBdeaths_SC[i, j, k] >0 && NeventsSC5[i, j, k]>0) rbinom(NeventsSC5[i, j, k],p_TBdeaths_SC[i, j, k]) else 0
NeventsSC6[,,] <- NeventsSC5[i,j,k] - TBdeaths_SC[i,j,k] #update



## g) IN MIGRATION
pSC_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pSC_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pSC_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
SC_inmigr[, , ] <- if(pSC_inmigr[i, j, k] > 0 && NeventsSC6[i, j, k]> 0) rbinom(NeventsSC6[i, j, k],pSC_inmigr[i, j, k]) else 0
NeventsSC7[,,] <- NeventsSC6[i,j,k] - SC_inmigr[i,j,k] #update


## g) BACKGROUND DEATH
SCdeaths[, , ] <- NeventsSC7[i, j, k]


## actual update
## print("SC: {min(SC[,,])}",when=min(SC[,,])<0)
update(SC[, , ]) <- SC[i, j, k] + age_in_SC[i, j, k] + HIV_in_SC[i, j, k] +
  progFast[i, j, k] + progSlow[i, j, k] + relapse[i, j, k] +
  regress[i, j, k] + SC_inmigr[i, j, k] -
  age_out_SC[i, j, k] - HIV_out_SC[i, j, k] - SCdeaths[i, j, k] - TBdeaths_SC[i, j, k] -
  detection_SC[i, j, k] - selfcure_SC[i, j, k] - progress[i, j, k]

#################################################
##### 4: prevalent infectious/active disease ####
#################################################


# ###### ACTIVE DISEASE: ENTRY EVENTS #####

## AGEING
age_in_D[,2:age_dims,] <- age_out_D[i,j-1,k]
age_in_D[,1,]<-0

## HIV
HIV_in_D[,,2:HIV_dims] <- HIV_out_D[i,j,k-1]
HIV_in_D[,,1]<-0

# ###### ACTIVE DISEASE: EXIT EVENTS #####
rate_D[,,1] <- m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate +
  age_rate[j] + HIV_rate_yr[j]
rate_D[,,2] <- m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate +
  age_rate[j] + ART_rate_yr[j]
rate_D[, , 3] <- m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate + age_rate[j]
NeventsD[, , ] <- if(D[i, j, k]>0) rbinom(D[i, j, k], 1 - exp(-rate_D[i, j, k] * dt)) else 0


## a) AGEING
p_Dage[, , ] <- if (rate_D[i, j, k] > tol) age_rate[j] / rate_D[i, j, k] else 0
age_out_D[, , ] <- if(p_Dage[i, j, k] >0 && NeventsD[i, j, k]>0) rbinom(NeventsD[i, j, k], p_Dage[i, j, k]) else 0
NeventsD1[,,] <- NeventsD[i,j,k] - age_out_D[i,j,k] #update


## b) HIV
p_DHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate + HIV_rate_yr[j] > tol) HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate + HIV_rate_yr[j]) else 0
p_DHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate + ART_rate_yr[j] > tol) ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate + ART_rate_yr[j]) else 0
p_DHIV[, , 3] <- 0
HIV_out_D[, , ] <- if(p_DHIV[i, j, k] >0 && NeventsD1[i, j, k]>0) rbinom(NeventsD1[i, j, k], p_DHIV[i, j, k]) else 0
NeventsD2[,,] <- NeventsD1[i,j,k] - HIV_out_D[i,j,k] #update



# ## c) TB DETECTION
p_detect[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate > tol) dtct_ratem[i] / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate) else 0
p_detect[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate > tol) dtct_ratem[i] / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate) else 0
p_detect[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate > tol) dtct_ratem[i] / (m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_ratem[i] + slfcr_rate + regress_rate) else 0
detection[, , ] <- if(p_detect[i, j, k] >0 && NeventsD2[i, j, k]>0) rbinom(NeventsD2[i, j, k], p_detect[i, j, k]) else 0
NeventsD3[,,] <- NeventsD2[i,j,k] - detection[i,j,k] #update


# ## d) SELF-CURING
p_selfcr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + slfcr_rate + regress_rate > tol) slfcr_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + slfcr_rate + regress_rate) else 0
p_selfcr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + slfcr_rate + regress_rate > tol) slfcr_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + slfcr_rate + regress_rate) else 0
p_selfcr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + slfcr_rate + regress_rate > tol) slfcr_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate + slfcr_rate + regress_rate) else 0
selfcure[, , ] <- if(p_selfcr[i, j, k] >0 && NeventsD3[i, j, k]>0) rbinom(NeventsD3[i, j, k],p_selfcr[i, j, k]) else 0
NeventsD4[,,] <- NeventsD3[i,j,k] - selfcure[i,j,k] #update


# ## e) REGRESSION to Subclinical
p_regress[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + regress_rate > tol) regress_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + regress_rate) else 0
p_regress[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + regress_rate > tol) regress_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + regress_rate) else 0
p_regress[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + regress_rate > tol) regress_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate + regress_rate) else 0
regress[, , ] <- if(p_regress[i, j, k] >0 && NeventsD4[i, j, k]>0) rbinom(NeventsD4[i, j, k],p_regress[i, j, k]) else 0
NeventsD5[,,] <- NeventsD4[i,j,k] - regress[i,j,k] #update


# ## f) TB DEATH
p_TBdeaths[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate > tol) TBd_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate) else 0
p_TBdeaths[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate > tol) TBd_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate) else 0
p_TBdeaths[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate > tol) TBd_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate) else 0
TBdeaths[, , ] <- if(p_TBdeaths[i, j, k] >0 && NeventsD5[i, j, k]>0) rbinom(NeventsD5[i, j, k],p_TBdeaths[i, j, k]) else 0
NeventsD6[,,] <- NeventsD5[i,j,k] - TBdeaths[i,j,k] #update


## g) IN MIGRATION
pD_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pD_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pD_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
D_inmigr[, , ] <- if(pD_inmigr[i, j, k]>0 && NeventsD6[i, j, k]>0) rbinom(NeventsD6[i, j, k],pD_inmigr[i, j, k]) else 0
NeventsD7[,,] <- NeventsD6[i,j,k] - D_inmigr[i,j,k] #update


## h) BACKGROUND DEATH

Ddeaths[, , ] <- NeventsD7[i, j, k]

## actual update
## print("D: {min(D[,,])}",when=min(D[,,])<0)
update(D[, , ]) <- D[i, j, k] + age_in_D[i, j, k] + HIV_in_D[i, j, k] + progress[i, j, k] +
  D_inmigr[i, j, k] -
  age_out_D[i, j, k] - HIV_out_D[i, j, k] - Ddeaths[i, j, k] - TBdeaths[i, j, k] -
  detection[i, j, k] - selfcure[i, j, k] - regress[i, j, k]

###########################
##### 5: On treatment #####
###########################



###### ON TREATMENT: ENTRY EVENTS #####

## b) AGEING
age_in_Tr[,2:age_dims,] <- age_out_Tr[i,j-1,k]
age_in_Tr[,1,]<-0

## c) HIV
HIV_in_Tr[,,2:HIV_dims] <- HIV_out_Tr[i,j,k-1]
HIV_in_Tr[,,1]<-0

###### On treatment: Calculate exit events (3) #####
rate_Tr[, , 1] <- m_in_t[j] + mu_noHIV_t[j] + Trsucc_rate + Trxd_rate + age_rate[j] + HIV_rate_yr[j]
rate_Tr[, , 2] <- m_in_t[j] + mu_HIV_t[j] + Trsucc_rate + Trxd_rate + age_rate[j] + ART_rate_yr[j]
rate_Tr[, , 3] <- m_in_t[j] + mu_ART_t[j] + Trsucc_rate + Trxd_rate + age_rate[j]
Trsucc_rate <- (1 - tfr) / t_dur
Trxd_rate <- tfr / t_dur
NeventsT[, , ] <- if(Tr[i, j, k]>0) rbinom(Tr[i, j, k], 1 - exp(-rate_Tr[i, j, k] * dt)) else 0


## a) AGEING
p_Trage[, , ] <- if (rate_Tr[i, j, k] > tol) age_rate[j] / rate_Tr[i, j, k] else 0
age_out_Tr[, , ] <- if(p_Trage[i, j, k] >0 && NeventsT[i, j, k]>0) rbinom(NeventsT[i, j, k], p_Trage[i, j, k]) else 0
NeventsT1[,,] <- NeventsT[i,j,k] - age_out_Tr[i,j,k] #update


## b) HIV
p_TrHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + Trsucc_rate + Trxd_rate + HIV_rate_yr[j] > tol) HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + Trsucc_rate + Trxd_rate + HIV_rate_yr[j]) else 0
p_TrHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + Trsucc_rate + Trxd_rate + ART_rate_yr[j] > tol) ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + Trsucc_rate + Trxd_rate + ART_rate_yr[j]) else 0
p_TrHIV[, , 3] <- 0
HIV_out_Tr[, , ] <- if(p_TrHIV[i, j, k] >0 && NeventsT1[i, j, k]>0) rbinom(NeventsT1[i, j, k], p_TrHIV[i, j, k]) else 0
NeventsT2[,,] <- NeventsT1[i,j,k] - HIV_out_Tr[i,j,k] #update



## c) TREATMENT SUCCESS
p_Trsuccess[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + Trsucc_rate + Trxd_rate > tol) Trsucc_rate / (m_in_t[j] + mu_noHIV_t[j] + Trsucc_rate + Trxd_rate) else 0
p_Trsuccess[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + Trsucc_rate + Trxd_rate > tol) Trsucc_rate / (m_in_t[j] + mu_HIV_t[j] + Trsucc_rate + Trxd_rate) else 0
p_Trsuccess[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + Trsucc_rate + Trxd_rate > tol) Trsucc_rate / (m_in_t[j] + mu_ART_t[j] + Trsucc_rate + Trxd_rate) else 0
Trsuccess[, , ] <- if(p_Trsuccess[i, j, k] >0 && NeventsT2[i, j, k]>0) rbinom(NeventsT2[i, j, k],p_Trsuccess[i, j, k]) else 0
NeventsT3[,,] <- NeventsT2[i,j,k] - Trsuccess[i,j,k] #update



## d) TB DEATH
p_Trxdeaths[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + Trxd_rate > tol) Trxd_rate / (m_in_t[j] + mu_noHIV_t[j] + Trxd_rate) else 0
p_Trxdeaths[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + Trxd_rate > tol) Trxd_rate / (m_in_t[j] + mu_HIV_t[j] + Trxd_rate) else 0
p_Trxdeaths[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + Trxd_rate > tol) Trxd_rate / (m_in_t[j] + mu_ART_t[j] + Trxd_rate) else 0
Trxdeaths[, , ] <- if(p_Trxdeaths[i, j, k] >0 && NeventsT3[i, j, k]>0) rbinom(NeventsT3[i, j, k],p_Trxdeaths[i, j, k]) else 0
NeventsT4[,,] <- NeventsT3[i,j,k] - Trxdeaths[i,j,k] #update

## e) IN MIGRATION
pTr_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pTr_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pTr_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
Tr_inmigr[, , ] <- if(pTr_inmigr[i, j, k]>0 && NeventsT4[i, j, k]>0)rbinom(NeventsT4[i, j, k],pTr_inmigr[i, j, k]) else 0
NeventsT5[,,] <- NeventsT4[i,j,k] - Tr_inmigr[i,j,k] #update


## e) BACKGROUND DEATH
Trdeaths[, , ] <- NeventsT5[i, j, k]

## actual update
## print("Tr: {min(Tr[,,])}",when=min(Tr[,,])<0)
update(Tr[, , ]) <- Tr[i, j, k] + age_in_Tr[i, j, k] + HIV_in_Tr[i, j, k] +
  detection[i, j, k] + detection_SC[i, j, k] + Tr_inmigr[i, j, k] -
  age_out_Tr[i, j, k] - HIV_out_Tr[i, j, k] - Trdeaths[i, j, k] -
  Trsuccess[i, j, k] - Trxdeaths[i, j, k]

########################
##### 6: Recovered #####
########################


###### RECOVERED: ENTRY EVENTS #####

## b) AGEING
age_in_R[,2:age_dims,] <- age_out_R[i,j-1,k]
age_in_R[,1,]<-0

## c) HIV
HIV_in_R[,,2:HIV_dims] <- HIV_out_R[i,j,k-1]
HIV_in_R[,,1]<-0

###### RECOVERED: EXIT EVENTS #####
rate_R[, , 1] <- m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] + pDr * IRR[i] * Hirr[1] + age_rate[j] + HIV_rate_yr[j]
rate_R[, , 2] <- m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] + pDr * IRR[i] * Hirr[2] + age_rate[j] + ART_rate_yr[j]
rate_R[, , 3] <- m_in_t[j] + mu_ART_t[j] + v * foi[i] * TB_HIV_mod[3] + pDr * IRR[i] * Hirr[3] + age_rate[j]
NeventsR[, , ] <- if(R[i, j, k]>0) rbinom(R[i, j, k], 1 - exp(-rate_R[i, j, k] * dt)) else 0


## a) AGEING
p_Rage[, , ] <- if (rate_R[i, j, k] > tol) age_rate[j] / rate_R[i, j, k] else 0
age_out_R[, , ] <- if(p_Rage[i, j, k] >0 && NeventsR[i, j, k]>0) rbinom(NeventsR[i, j, k], p_Rage[i, j, k]) else 0
NeventsR1[,,] <- NeventsR[i,j,k] - age_out_R[i,j,k] #update


## b) HIV rate
p_RHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] + pDr * IRR[i] * Hirr[1] + HIV_rate_yr[j] > tol) HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] + pDr * IRR[i] * Hirr[1] + HIV_rate_yr[j]) else 0
p_RHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] + pDr * IRR[i] * Hirr[2] + ART_rate_yr[j] > tol) ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] + pDr * IRR[i] * Hirr[2] + ART_rate_yr[j]) else 0
p_RHIV[, , 3] <- 0
HIV_out_R[, , ] <- if(p_RHIV[i, j, k] >0 && NeventsR1[i, j, k]>0) rbinom(NeventsR1[i, j, k],p_RHIV[i, j, k]) else 0
NeventsR2[,,] <- NeventsR1[i,j,k] - HIV_out_R[i,j,k] #update


## c) REINFECTION
p_Rinfs[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] + pDr * IRR[i] * Hirr[1] > tol) v * foi[i] * TB_HIV_mod[1] / (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] + pDr * IRR[i]* Hirr[1]) else 0
p_Rinfs[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] + pDr * IRR[i] * Hirr[2] > tol) v * foi[i] * TB_HIV_mod[1] / (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] + pDr * IRR[i] * Hirr[2]) else 0
p_Rinfs[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + v * foi[i] * TB_HIV_mod[3] + pDr * IRR[i] * Hirr[3] > tol) v * foi[i] * TB_HIV_mod[1] / (m_in_t[j] + mu_ART_t[j] + v * foi[i] * TB_HIV_mod[3] + pDr * IRR[i] * Hirr[3]) else 0
Rinfs[, , ] <- if(p_Rinfs[i, j, k] >0 && NeventsR2[i, j, k]>0) rbinom(NeventsR2[i, j, k],p_Rinfs[i, j, k]) else 0
NeventsR3[,,] <- NeventsR2[i,j,k] - Rinfs[i,j,k] #update


## d) RELAPSE
p_relapse[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + pDr * IRR[i] * Hirr[1] > tol) IRR[i] * pDr * Hirr[1] / (m_in_t[j] + mu_noHIV_t[j] + pDr * IRR[i] * Hirr[1]) else 0
p_relapse[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + pDr * IRR[i] * Hirr[2] > tol) IRR[i] * pDr * Hirr[2] / (m_in_t[j] + mu_HIV_t[j] + pDr * IRR[i] * Hirr[2]) else 0
p_relapse[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + pDr * IRR[i] * Hirr[3] > tol) IRR[i] * pDr * Hirr[3] / (m_in_t[j] + mu_ART_t[j] + pDr * IRR[i] * Hirr[3]) else 0
relapse[, , ] <- if(p_relapse[i, j, k] >0 && NeventsR3[i, j, k]>0) rbinom(NeventsR3[i, j, k],p_relapse[i, j, k]) else 0
NeventsR4[,,] <- NeventsR3[i,j,k] - relapse[i,j,k] #update

## e) IN MIGRATION
pR_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pR_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pR_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
R_inmigr[, , ] <- if(pR_inmigr[i, j, k]>0 && NeventsR4[i, j, k]>0) rbinom(NeventsR4[i, j, k],pR_inmigr[i, j, k]) else 0
NeventsR5[,,] <- NeventsR4[i,j,k] - R_inmigr[i,j,k] #update

## e) BACKGROUND DEATH
Rdeaths[, , ] <- NeventsR5[i, j, k]

## actual update
## print("R: {min(R[,,])}",when=min(R[,,])<0)

update(R[, , ]) <- R[i, j, k] + age_in_R[i, j, k] + HIV_in_R[i, j, k] +
  Trsuccess[i, j, k] + R_inmigr[i, j, k] -
  age_out_R[i, j, k] - HIV_out_R[i, j, k] - Rdeaths[i, j, k] - Rinfs[i, j, k] - relapse[i, j, k]

###### Additional vars to record:
# NOTE: These will update based on values passed at the beginning of timestep
# (So are now effectively a step behind where they were prior to odin.dust which
# utilised output() rather than update())

update(notes[,,]) <- detection[i,j,k]+detection_SC[i,j,k]
update(notifrate[,,]) <- 1e5*(detection[i,j,k]+detection_SC[i,j,k])/( N[i,j,k] + 1e-10 )
update(TB_deaths[,,]) <- Trxdeaths[i,j,k] + TBdeaths[i,j,k] + TBdeaths_SC[i,j,k]
update(N[,,]) <- U[i,j,k] + LR[i,j,k] + LL[i,j,k] + D[i,j,k] + SC[i,j,k] + Tr[i,j,k] + R[i,j,k]
update(bg_deaths[,,]) <- Udeaths[i,j,k] + LRdeaths[i,j,k] + LLdeaths[i,j,k] +
  Ddeaths[i,j,k] + SCdeaths[i,j,k] + Trdeaths[i,j,k] + Rdeaths[i,j,k]
update(incidence[,,]) <- progFast[i,j,k] + progSlow[i,j,k] + relapse[i,j,k]
update(tot_incidence) <- sum(progFast[,,]) + sum(progSlow[,,]) + sum(relapse[,,])

update(incidence_bypatch[]) <- 1e5 * (sum(progFast[i,,]) + sum(progSlow[i,,]) + sum(relapse[i,,])) / (sum(N[i, , ]) + 1e-15)
update(prevalence_bypatch[]) <- 1e5 * (sum(D[i, , ]) + 1* sum(SC[i, , ])) / (sum(N[i, , ]) + 1e-15)

## ## patch-level total notifications,popn, and rates
## pnotes[] <- sum(detection[i, , ] + detection_SC[i, , ])
ppop[] <- sum(N[i, , ])
## pnotifrate[] <- 1e5 * pnotes[i] / (ppop[i] + 1e-10)

## #########################
## test variables
update(beta_test) <- beta_test #constant at initial value passed in



## =========== flux approximation calculations
update(Tijk[,,]) <- exp(-zk[k]) * Tijk[i,j,k] + InfsByPatchPatch[i, j]
update(Sij[, ]) <- exp(-zk[1]) * ( Tijk[i, j, 1] + Sij[i,j ] )

## \ell_{ij}(t) = \sum_{s=1}^t\Lambda_{ij}(t-s)p(s)
ellij[, ] <- abs(A0 * Sij[i, j] + A[1] * Tijk[i, j, 1] + A[2] * Tijk[i, j, 2] + A[3] * Tijk[i, j, 3] + A[4] * Tijk[i, j, 4] + A[5] * Tijk[i, j, 5] + A[6] * Tijk[i, j, 6])


## parameter mapping CHECK
## NOTES -> model
## omega  :  1/dur      #tb cessation
## delta  :  dtct_rate dtct_ratem[i] #detection rate
## mu     :  mu_noHIV_int[2:age_dims,1:sim_length] #mortality
## sigma  :  pLL #stabilisation
## alpha  :  pDf #fast progression
## epsilon:  pDs #slow progn
## gamma  :  progress_rate-regress_rate #symptom progn
## rho    :  pDr  #relapse rate
## tau    :  1/t_dur    #1/treatment dur
Pomega  <-  1/dur      #tb cessation
Pdelta  <-  dtct_rate  #detection rate TODO
Pmu     <-  mu_noHIV_int[2,1] #mortality: TODO possible to make mean?
Psigma  <-  pLL #stabilisation
Palpha  <-  pDf #fast progression
Pepsilon<-  pDs #slow progn
Pgamma  <-  progress_rate-regress_rate #symptom progn
Prho    <-  pDr  #relapse rate
Ptau    <-  1/t_dur    #1/treatment dur


## z coeffs
zk[1] <- Pomega + Pdelta + Pmu # om
zk[2] <- Psigma + Palpha + Pmu # sam
zk[3] <- Pepsilon + Pmu # em
zk[4] <- Pgamma + Pmu # gm
zk[5] <- Prho + Pmu # rr
zk[6] <- Ptau + Pmu # tm

## A coeffs
A0 <- Prho * Ptau * (ppJ /((Prho-Pomega-Pdelta)*(Ptau-Pomega-Pdelta)))
A[1] <- -Pgamma * ppC /((Pgamma-Psigma-Palpha)*(Pomega+Pdelta-Psigma-Palpha)) + Pgamma * ppB /((Pgamma-Pepsilon)*(Pomega+Pdelta-Pepsilon)) -  Pgamma * (ppB/(Pgamma-Pepsilon)-ppC/(Pgamma-Psigma-Palpha))  / (Pomega+Pdelta-Pgamma)- Prho * Ptau * (ppF*(1/(Pomega+Pdelta-Psigma-Palpha)-1/(Pomega+Pdelta-Prho))/((Prho-Psigma-Palpha)*(Ptau-Psigma-Palpha))+ppG*(1/(Pomega+Pdelta-Pepsilon)-1/(Pomega+Pdelta-Prho))/((Prho-Pepsilon)*(Ptau-Pepsilon)) + ppH*(1/(Pomega+Pdelta-Pgamma)-1/(Pomega+Pdelta-Prho))/((Prho-Pgamma)*(Ptau-Pgamma)) + ppJ*(-1/(Pomega+Pdelta-Prho))/((Prho-Pomega-Pdelta)*(Ptau-Pomega-Pdelta))+ ppK *(1/(Pomega+Pdelta-Ptau)-1/(Pomega+Pdelta-Prho))/(Prho-Ptau))
A[2] <- Pgamma * ppC /((Pgamma-Psigma-Palpha)*(Pomega+Pdelta-Psigma-Palpha)) + Prho * Ptau * (ppF/(Pomega+Pdelta-Psigma-Palpha))/((Prho-Psigma-Palpha)*(Ptau-Psigma-Palpha))
A[3] <- - Pgamma * ppB /((Pgamma-Pepsilon)*(Pomega+Pdelta-Pepsilon)) + Prho * Ptau *  ( ppG/(Pomega+Pdelta-Pepsilon))/((Prho-Pepsilon)*(Ptau-Pepsilon))
A[4] <- Pgamma * (ppB/(Pgamma-Pepsilon)-ppC/(Pgamma-Psigma-Palpha))/ (Pomega+Pdelta-Pgamma) + Prho * Ptau * (ppH/(Pomega+Pdelta-Pgamma) )/ ((Prho-Pgamma)*(Ptau-Pgamma))
A[5] <- -Prho * Ptau * ( ppF/((Prho-Psigma-Palpha)*(Ptau-Psigma-Palpha)) + ppG/ ((Prho-Pepsilon)*(Ptau-Pepsilon)) + ppH/((Prho-Pgamma)*(Ptau-Pgamma)) + ppJ/((Prho-Pomega-Pdelta)*(Ptau-Pomega-Pdelta)) + ppK/(Prho-Ptau) )/(Pomega+Pdelta-Prho)
A[6] <- Prho * Ptau * ((ppK / (Pomega + Pdelta - Ptau)) / (Prho - Ptau))

## interim variables
ppB <- Pepsilon * Psigma / (Pepsilon - Psigma - Palpha)
ppC <- ppB + Palpha
ppF <- Pdelta * Pgamma * ppC / ((Pgamma - Psigma - Palpha) * (Pomega + Pdelta - Psigma - Palpha))
ppG <- -Pdelta * Pgamma * ppB / ((Pgamma - Pepsilon) * (Pomega + Pdelta - Pepsilon))
ppH <- Pdelta * Pgamma * (ppB / (Pgamma - Pepsilon) - ppC / (Pgamma - Psigma - Palpha)) / (Pomega + Pdelta - Pgamma)
ppJ <- Pdelta * Pgamma * (ppC * (1 / (Pomega + Pdelta - Pgamma) - 1 / (Pomega + Pdelta - Psigma - Palpha)) / (Pgamma - Psigma - Palpha) + ppB * (1 / (Pomega + Pdelta - Pepsilon) - 1 / (Pomega + Pdelta - Pgamma)) / (Pgamma - Pepsilon))
ppK <- -ppF / (Ptau - Psigma - Palpha) - ppG / (Ptau - Pepsilon) - ppH / (Ptau - Pgamma) - ppJ / (Ptau - Pomega - Pdelta)

