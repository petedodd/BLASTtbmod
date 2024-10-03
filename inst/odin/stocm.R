
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
popinit_byage[,] <- rbinom(floor(popinit[i]),agefracs[j]/(sum(agefracs)+1e-10)) #division just in case


## TB initial state
initPrev[,] <- 1-exp( -ari0*ageMids[j] )             #LTBI


## fraction TBI initially R or LR: as ~2 year's of FOI/TBI
initF[,] <- if(initPrev[i,j]>tol) (1-exp(-2*ari0)) / initPrev[i,j] else 0
## Initial ratio of TBI to disease states
initD[,] <- dur*( initF[i,j]*(pDf-pDs) + pDs )
initLL[,] <- 1.0 #basically rest given denom
initDenom[,] <- initF[i,j] + initLL[i,j] + initD[i,j]
tbi_LR[,] <- if(initDenom[i,j] > tol) (initF[i,j])/initDenom[i,j] else 0
tbi_LL[,] <- if(initDenom[i,j] > tol) (initLL[i,j])/initDenom[i,j] else 0
tbi_D[,] <- if(initDenom[i,j] > tol) (initD[i,j]/4)/initDenom[i,j] else 0
tbi_SC[,] <- if(initDenom[i,j] > tol) (initD[i,j]/4)/initDenom[i,j] else 0
tbi_Tr[,] <- if(initDenom[i,j] > tol) (initD[i,j]/4)/initDenom[i,j] else 0
tbi_R[,] <- if(initDenom[i,j] > tol) (initD[i,j]/4)/initDenom[i,j] else 0


## NOTE total stochastic even given above
init_U[,] <- rbinom(popinit_byage[i,j],(1-initPrev[i,j]))
init_LR[,] <- rbinom(popinit_byage[i,j],initPrev[i,j]*tbi_LR[i,j])
init_LL[,] <- rbinom(popinit_byage[i,j],initPrev[i,j]*tbi_LL[i,j])
init_D[,] <- rbinom(popinit_byage[i,j],initPrev[i,j]*tbi_D[i,j])
init_SC[,] <- rbinom(popinit_byage[i,j],initPrev[i,j]*tbi_SC[i,j])
init_Tr[,] <- rbinom(popinit_byage[i,j],initPrev[i,j]*tbi_Tr[i,j])
init_R[,] <- rbinom(popinit_byage[i,j],initPrev[i,j]*tbi_R[i,j])



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
dtct_rate <- cdr/(dur*(1-cdr))   # Detection rate
slfcr_rate <- (1-cfr)/dur        # Selfcure rate
## m_in <- user(0.0)               # m_in - Migration rate in
age_rate[] <-user()
regress_rate <- user()
progress_rate <- user()  


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
# HIV_rate_yr[1:3] <- 0
# HIV_rate_yr[4:10] <- HIV_t/1000
# HIV_rate_yr[11:18] <- 0
HIV_rate_yr[1] <- 0
HIV_rate_yr[2] <- HIV_t/1000
HIV_rate_yr[3] <- 0

# Assign ART initiation rate
# Time & age dependent
ART_int[] <- user()
ART_t <- if(as.integer(step) < sim_length )
  ART_int[step+1] else ART_int[sim_length]
# ART_rate_yr[1:3] <- 0
# ART_rate_yr[4:18] <- ART_t
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
foitemp[, ] <- beta_seas * MM[i, j] * (sum(D[j, , ]) + relinf * sum(SC[j, , ])) / (sum(N[j, , ]) + 1e-15) # force of infection on patch i from patch j: contact rate x prev in patch j
foi[] <- sum( foitemp[i,] ) #FOI on each patch i: sum over FOIs from each patch
MM[,] <- user()           # Mixing matrix 
relinf <- user(1)         # Relative infectiousness of SC relative to CD
IRR[] <- user()

initial(cum_inf_flux[, ]) <- 0 #initialize cumulative fluxes


## TODO TB_HIV_mod acts only on infection: need HIV specific progression
## TODO IRR acts only on patch
## both the above need to be 1 really


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


## Dims relating to U compartment (uninfected)
dim( U ) <- c( patch_dims, age_dims, HIV_dims )
dim( init_U ) <- c( patch_dims, age_dims )
dim( Udeaths) <- c( patch_dims, age_dims, HIV_dims )
dim( Uinfs) <- c( patch_dims, age_dims, HIV_dims )
## dim( m_in_U ) <- c( patch_dims, age_dims, HIV_dims )
dim( Nevents_U ) <- c( patch_dims, age_dims, HIV_dims )
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
dim( Nevents_LR ) <- c( patch_dims, age_dims, HIV_dims )
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
dim( Nevents_LL ) <- c( patch_dims, age_dims, HIV_dims )
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
dim( Nevents_D ) <- c( patch_dims, age_dims, HIV_dims )
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
dim( Nevents_SC ) <- c( patch_dims, age_dims, HIV_dims )
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
dim( Nevents_Tr ) <- c( patch_dims, age_dims, HIV_dims )
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
dim( Nevents_R ) <- c( patch_dims, age_dims, HIV_dims )
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
dim( popinit ) <- patch_dims
dim( IRR ) <- patch_dims
dim( tt ) <- sim_length
dim( popinit_byage ) <- c( patch_dims, age_dims )
dim( initPrev ) <- c( patch_dims, age_dims )
dim( initF ) <- c( patch_dims, age_dims )
dim( initD ) <- c( patch_dims, age_dims )
# dim( ari0 ) <- patch_dims
dim( ageMids ) <- age_dims
dim( agefracs ) <- age_dims


dim(initDenom) <- c( patch_dims, age_dims )
dim(tbi_LR) <- c( patch_dims, age_dims )
dim(tbi_LL) <- c( patch_dims, age_dims )
dim(tbi_D) <- c( patch_dims, age_dims )
dim(tbi_SC) <- c( patch_dims, age_dims )
dim(tbi_Tr) <- c( patch_dims, age_dims )
dim(tbi_R) <- c( patch_dims, age_dims )

## ## patch level aggregate indicators
## dim(pnotifrate) <- c(patch_dims)
## dim(pnotes) <- c(patch_dims)
## dim(ppop) <- c(patch_dims)

## flux variables
dim(InfsByPatchPatch) <- c(patch_dims, patch_dims)
dim(InfsByPatch) <- c(patch_dims)
dim(cum_inf_flux) <- c(patch_dims, patch_dims)


##########
## Additional outputs to record
##########

# All suppressed for now
# Can use this file later to 'mirror' any vars needed as output

# output(N[,,]) <- TRUE
# output(TB_deaths) <- TRUE
# output(notifrate) <- TRUE
# output(progress) <- TRUE
# output(regress) <- TRUE
# output( bg_deaths ) <- TRUE

################
# Temporary outputs for diagnostics
################
# output( age_in_U ) <- TRUE
# output( age_out_U ) <- TRUE
# output( births ) <- TRUE
# output( age_out_LR ) <- TRUE
# output( Nevents_LR ) <- TRUE
# output( LRdeaths ) <- TRUE
# output( progFast ) <- TRUE
# output( stabilisations ) <- TRUE
# output( beta_seas ) <- TRUE
# output( p_detect) <- TRUE
# output( HIV_int ) <- TRUE
# output( HIV_rate_yr ) <- TRUE
# output( ART_rate_yr ) <- TRUE
# output( HIV_out_U[,,] ) <- TRUE
# output( HIV_in_U[,,] ) <- TRUE
# output( births[] ) <- TRUE

## output(pnotifrate) <- TRUE
## output(pnotes) <- TRUE
## output(ppop) <- TRUE

##########################
## Populations
##########################

# S (U) - susceptible (uninfected)
# E (LR) - early infection
# L (LL) - late infection
# I (D) - prevalent + infectious
# *new* (SCD) - Subclinical disease - infectious?
# T (Tr) - receiving treatment, assumed not infections
# R (R) - recovered at risk of relapse


## importance/priority groups: 1/-3/

## 1/
## TODO HIV/ART dependent mortality
## TODO make IRRs time-dependent

## 2/
## TODO sublinical D and clinical D state split

## 3/
## TODO between patch migration? perhaps shelve until discussion


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

update(U[,,]) <- U[i,j,k] + births[i,j,k] + age_in_U[i,j,k] + HIV_in_U[i,j,k] - ## m_in_U[i,j,k] -
  age_out_U[i, j, k] - HIV_out_U[i, j, k] - Udeaths[i, j, k] - Uinfs[i, j, k]+
  U_inmigr[i, j, k]


###### UNINFECTED ENTRY EVENTS #######
### a) BIRTHS + (net?) MIGRATION
# Birth rate downloaded from https://population.un.org/wpp/

births[ 1:patch_dims, 1, 1 ]  <- rpois( birth_rate_yr*sum(N[i,,])*dt )
births[ 1:patch_dims, 2:age_dims, 2:HIV_dims ] <- 0
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
Nevents_U[,,] <- rbinom( U[i,j,k], 1-exp( -rate_U[i,j,k]*dt ))

## NOTE denominator safety introduced

## a) AGEING
p_Uage[,,] <- if(rate_U[i,j,k] > tol) age_rate[j]/rate_U[i,j,k] else 0
age_out_U[,,] <- rbinom( Nevents_U[i,j,k], p_Uage[i,j,k] )

## b) HIV INFECTION
## of remaining events (bar ageing), what proportion are HIV processes?
## NOTE safety
p_UHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + foi[i] * TB_HIV_mod[1] + HIV_rate_yr[j] > tol) HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + foi[i] * TB_HIV_mod[1] + HIV_rate_yr[j]) else 0
p_UHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + foi[i] * TB_HIV_mod[2] + ART_rate_yr[j] > tol) ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + foi[i] * TB_HIV_mod[2] + ART_rate_yr[j]) else 0
p_UHIV[,,3] <- 0
HIV_out_U[,,] <- rbinom( Nevents_U[i,j,k]-age_out_U[i,j,k], p_UHIV[i,j,k])

## c) TB INFECTION
p_Uinf[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + foi[i] * TB_HIV_mod[1] > tol) foi[i] * TB_HIV_mod[1] / (m_in_t[j] + mu_noHIV_t[j] + foi[i] * TB_HIV_mod[1]) else 0
p_Uinf[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + foi[i] * TB_HIV_mod[2] > tol) foi[i] * TB_HIV_mod[2] / (m_in_t[j] + mu_HIV_t[j] + foi[i] * TB_HIV_mod[2]) else 0
p_Uinf[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + foi[i] * TB_HIV_mod[3] > tol) foi[i] * TB_HIV_mod[3] / (m_in_t[j] + mu_ART_t[j] + foi[i] * TB_HIV_mod[3]) else 0
Uinfs[,,] <- rbinom( Nevents_U[i,j,k]-age_out_U[i,j,k]-HIV_out_U[i,j,k], p_Uinf[i,j,k])

## d) IN MIGRATION
pU_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / m_in_t[j] + mu_noHIV_t[j] else 0
pU_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pU_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
## print("mint: {min(m_in_t[])}")
## print("pU: {min(pU_inmigr[,,])}")
## print("pU: {max(pU_inmigr[,,])}", when = min(pU_inmigr[, , ]) > 0)
U_inmigr[, , ] <- if(pU_inmigr[i, j, k]>0 && Nevents_U[i, j, k] - age_out_U[i, j, k] - HIV_out_U[i, j, k] - Uinfs[i, j, k] > 0) rbinom(
  Nevents_U[i, j, k] - age_out_U[i, j, k] - HIV_out_U[i, j, k] - Uinfs[i, j, k],
  pU_inmigr[i, j, k]
) else 0
## print("min U: {min(U_inmigr[,,])}")
## print("max U: {max(U_inmigr[,,])}")


## e) DEATH
Udeaths[, , ] <- Nevents_U[i, j, k] - age_out_U[i, j, k] - HIV_out_U[i, j, k] - Uinfs[i, j, k]-
  U_inmigr[i, j, k]

## ---------------------- extra calculations associated with infection fluxes
InfsByPatch[] <- sum(Uinfs[i, , ]) # number of infections in patch i this step
InfsByPatchPatch[, ] <- rbinom(InfsByPatch[i], foitemp[i, j] / (foi[i] + tol))
update(cum_inf_flux[, ]) <- cum_inf_flux[i, j] + InfsByPatchPatch[i, j]

##############################
##### 2: Infected early ######
##############################

# print("LR: {min(LR[,,])}",when=min(LR[,,])<0)
update(LR[, , ]) <- LR[i, j, k] + age_in_LR[i, j, k] + HIV_in_LR[i, j, k] + Uinfs[i, j, k] +
  LLinfs[i, j, k] + Rinfs[i, j, k] + LR_inmigr[i, j, k] -
  age_out_LR[i, j, k] - HIV_out_LR[i, j, k] - LRdeaths[i, j, k] -
  progFast[i, j, k] - stabilisations[i, j, k]

###### INFECTED EARLY: ENTRY EVENTS #####
rate_LR[, , 1] <- mu_noHIV_t[j] + pDf * IRR[i] * Hirr[1] + pLL + age_rate[j] + HIV_rate_yr[j]
rate_LR[, , 2] <- mu_HIV_t[j] + pDf * IRR[i] * Hirr[2] + pLL + age_rate[j] + ART_rate_yr[j]
rate_LR[, , 3] <- mu_ART_t[j] + pDf * IRR[i] * Hirr[3] + pLL + age_rate[j]
Nevents_LR[, , ] <- rbinom(LR[i, j, k], 1 - exp(-rate_LR[i, j, k] * dt))

## b) AGEING
age_in_LR[, 2:age_dims, ] <- age_out_LR[i, j - 1, k]
age_in_LR[, 1, ] <- 0

## c) HIV
HIV_in_LR[, , 2:HIV_dims] <- HIV_out_LR[i, j, k - 1]
HIV_in_LR[, , 1] <- 0

## a) AGEING
p_LRage[,,] <- if(rate_LR[i,j,k] > tol) age_rate[j]/rate_LR[i,j,k] else 0
age_out_LR[,,] <- rbinom( Nevents_LR[i,j,k], p_LRage[i,j,k] )

## b) HIV INFECTION
p_LRHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + pDf * IRR[i] * Hirr[1] + pLL + HIV_rate_yr[j] > tol)
                    HIV_rate_yr[j] /
                      (m_in_t[j] + mu_noHIV_t[j] + pDf * IRR[i]* Hirr[1] + pLL + HIV_rate_yr[j]) else 0
p_LRHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + pDf * IRR[i] * Hirr[2] + pLL + ART_rate_yr[j] > tol)
                    ART_rate_yr[j] /
                      (m_in_t[j] + mu_HIV_t[j] + pDf * IRR[i] * Hirr[2] + pLL + ART_rate_yr[j]) else 0
p_LRHIV[,,3] <- 0
HIV_out_LR[,,] <- rbinom( Nevents_LR[i,j,k]-age_out_LR[i,j,k], p_LRHIV[i,j,k])

## c) FAST PROGRESSION
p_progFast[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + pDf * IRR[i] * Hirr[1] + pLL > tol)
                       IRR[i] * pDf * Hirr[1] / (m_in_t[j] + mu_noHIV_t[j] + pDf * IRR[i] * Hirr[1] + pLL) else 0
p_progFast[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + pDf * IRR[i] * Hirr[2] + pLL > tol)
                       IRR[i] * pDf * Hirr[2] / (m_in_t[j] + mu_HIV_t[j] + pDf * IRR[i] * Hirr[2] + pLL) else 0
p_progFast[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + pDf * IRR[i] * Hirr[3] + pLL > tol)
                       IRR[i] * pDf * Hirr[3] / (m_in_t[j] + mu_ART_t[j] + pDf * IRR[i] * Hirr[3] + pLL) else 0
progFast[,,] <- rbinom( Nevents_LR[i,j,k]-age_out_LR[i,j,k]-HIV_out_LR[i,j,k],p_progFast[i,j,k])


## d) STABILISATION
## p_stabilise[,,] <- pLL/rate_LR[i,j,k]
p_stabilise[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + pLL > tol)
                        pLL / (m_in_t[j] + mu_noHIV_t[j] + pLL) else 0
p_stabilise[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + pLL > tol)
                        pLL / (m_in_t[j] + mu_HIV_t[j] + pLL) else 0
p_stabilise[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + pLL > tol)
                        pLL / (m_in_t[j] + mu_ART_t[j] + pLL) else 0
stabilisations[, , ] <- rbinom(
  Nevents_LR[i, j, k] - age_out_LR[i, j, k] - HIV_out_LR[i, j, k] - progFast[i, j, k],
  p_stabilise[i, j, k]
)


## e) IN MIGRATION
pLR_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pLR_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pLR_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
LR_inmigr[, , ] <- if(pLR_inmigr[i, j, k]>0 && Nevents_LR[i, j, k] - age_out_LR[i, j, k] - HIV_out_LR[i, j, k] -  progFast[i, j, k] - stabilisations[i, j, k] > 0) rbinom(
  Nevents_LR[i, j, k] - age_out_LR[i, j, k] - HIV_out_LR[i, j, k] -
  progFast[i, j, k] - stabilisations[i, j, k],
  pLR_inmigr[i, j, k]
) else 0

## f) DEATH
LRdeaths[, , ] <- Nevents_LR[i, j, k] - age_out_LR[i, j, k] - HIV_out_LR[i, j, k] - progFast[i, j, k] - stabilisations[i, j, k] - LR_inmigr[i, j, k]


############################
##### 3: Infected late #####
############################

# print("LL: {min(LL[,,])}",when=min(LL[,,])<0)

update(LL[, , ]) <- LL[i, j, k] + age_in_LL[i, j, k] + HIV_in_LL[i, j, k] +
  stabilisations[i, j, k] + selfcure[i, j, k] + selfcure_SC[i, j, k] -
  age_out_LL[i, j, k] - HIV_out_LL[i, j, k] - LLdeaths[i, j, k] - LLinfs[i, j, k] - progSlow[i, j, k] + LL_inmigr[i, j, k]

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
Nevents_LL[, , ] <- rbinom(LL[i, j, k], 1 - exp(-rate_LL[i, j, k] * dt))


## a) AGEING
p_LLage[, , ] <- if (rate_LL[i, j, k] > tol) age_rate[j] / rate_LL[i, j, k] else 0
age_out_LL[, , ] <- rbinom(Nevents_LL[i, j, k], p_LLage[i, j, k])


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
HIV_out_LL[, , ] <- rbinom(Nevents_LL[i, j, k] - age_out_LL[i, j, k], p_LLHIV[i, j, k])


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
LLinfs[, , ] <- rbinom(
  Nevents_LL[i, j, k] - age_out_LL[i, j, k] - HIV_out_LL[i, j, k],
  p_LLinfs[i, j, k]
)


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
progSlow[, , ] <- rbinom(
  Nevents_LL[i, j, k] - age_out_LL[i, j, k] - HIV_out_LL[i, j, k] - LLinfs[i, j, k],
  p_progSlow[i, j, k]
)


## e) IN MIGRATION
pLL_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pLL_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pLL_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
LL_inmigr[, , ] <- if(pLL_inmigr[i, j, k]>0 && Nevents_LL[i, j, k] - age_out_LL[i, j, k] - HIV_out_LL[i, j, k] - LLinfs[i, j, k] - progSlow[i, i, k]>0) rbinom(
  Nevents_LL[i, j, k] - age_out_LL[i, j, k] - HIV_out_LL[i, j, k] - LLinfs[i, j, k] - progSlow[i, i, k],
  pLL_inmigr[i, j, k]
) else 0
## print("min LL: {min(LL_inmigr[,,])}")
## print("max LL: {max(LL_inmigr[,,])}")

## f) DEATH
LLdeaths[, , ] <- Nevents_LL[i, j, k] - age_out_LL[i, j, k] - HIV_out_LL[i, j, k] -
  LLinfs[i, j, k] - progSlow[i, j, k] - LL_inmigr[i, j, k]



##########################################
##### 4a: Subclinical active disease #####
##########################################

# print("SC: {min(SC[,,])}",when=min(SC[,,])<0)
update(SC[, , ]) <- SC[i, j, k] + age_in_SC[i, j, k] + HIV_in_SC[i, j, k] +
  progFast[i, j, k] + progSlow[i, j, k] + relapse[i, j, k] +
  regress[i, j, k] + SC_inmigr[i, j, k] -
  age_out_SC[i, j, k] - HIV_out_SC[i, j, k] - SCdeaths[i, j, k] - TBdeaths_SC[i, j, k] -
  detection_SC[i, j, k] - selfcure_SC[i, j, k] - progress[i, j, k]

###### SC DISEASE: ENTRY EVENTS #####

## AGEING
age_in_SC[,2:age_dims,] <- age_out_SC[i,j-1,k]
age_in_SC[,1,]<-0

## HIV
HIV_in_SC[,,2:HIV_dims] <- HIV_out_SC[i,j,k-1]
HIV_in_SC[,,1]<-0

###### SC DISEASE: EXIT EVENTS #####
rate_SC[, , 1] <- m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate_SC + slfcr_rate +
  progress_rate + age_rate[j] + HIV_rate_yr[j]
rate_SC[, , 2] <- m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate_SC + slfcr_rate +
  progress_rate + age_rate[j] + ART_rate_yr[j]
rate_SC[, , 3] <- m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_rate_SC + slfcr_rate +
  progress_rate + age_rate[j]
Nevents_SC[, , ] <- rbinom(SC[i, j, k], 1 - exp(-rate_SC[i, j, k] * dt))
dtct_rate_SC <- cdr_SC / (dur * (1 - cdr_SC))

## TBd_rate & slfcr_rate defined with active disease


## a) AGEING
p_SCage[, , ] <- if (rate_SC[i, j, k] > tol) age_rate[j] / rate_SC[i, j, k] else 0
age_out_SC[, , ] <- rbinom(Nevents_SC[i, j, k], p_SCage[i, j, k])


## b) HIV
p_SCHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate_SC +
                      slfcr_rate + progress_rate + HIV_rate_yr[j] > tol)
                    HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate_SC +
                                      slfcr_rate + progress_rate + HIV_rate_yr[j]) else 0
p_SCHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate_SC +
                      slfcr_rate + progress_rate + ART_rate_yr[j] > tol)
                    ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate_SC +
                                      slfcr_rate + progress_rate + ART_rate_yr[j]) else 0
p_SCHIV[, , 3] <- 0
HIV_out_SC[, , ] <- rbinom(Nevents_SC[i, j, k] - age_out_SC[i, j, k], p_SCHIV[i, j, k])



## c) TB DETECTION
p_detect_SC[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate_SC +
                          slfcr_rate +progress_rate > tol)
                        dtct_rate_SC /(m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate_SC +
                                       slfcr_rate +progress_rate) else 0
p_detect_SC[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate_SC +
                          slfcr_rate +progress_rate > tol)
                        dtct_rate_SC /(m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate_SC +
                                       slfcr_rate +progress_rate) else 0
p_detect_SC[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_rate_SC +
                          slfcr_rate +progress_rate > tol)
                        dtct_rate_SC / (m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_rate_SC +
                                        slfcr_rate +progress_rate) else 0
detection_SC[, , ] <- rbinom(
  Nevents_SC[i, j, k] - age_out_SC[i, j, k] - HIV_out_SC[i, j, k],
  p_detect_SC[i, j, k]
)


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
selfcure_SC[, , ] <- rbinom(
  Nevents_SC[i, j, k] - age_out_SC[i, j, k] - HIV_out_SC[i, j, k] - detection_SC[i, j, k],
  p_selfcr_SC[i, j, k]
)


## e) PROGRESSION to active disease
p_progress[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + progress_rate > tol)
                       progress_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + progress_rate) else 0
p_progress[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + progress_rate > tol)
                       progress_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + progress_rate) else 0
p_progress[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + progress_rate > tol)
                       progress_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate + progress_rate) else 0
progress[, , ] <- rbinom(
  Nevents_SC[i, j, k] - age_out_SC[i, j, k] - HIV_out_SC[i, j, k] -
    detection_SC[i, j, k] - selfcure_SC[i, j, k],
  p_progress[i, j, k]
)


## f) TB DEATH
p_TBdeaths_SC[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate > tol)
                          TBd_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate) else 0
p_TBdeaths_SC[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate > tol)
                          TBd_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate) else 0
p_TBdeaths_SC[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate > tol)
                          TBd_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate) else 0
TBdeaths_SC[, , ] <- rbinom(
  Nevents_SC[i, j, k] - age_out_SC[i, j, k] - HIV_out_SC[i, j, k] -
    detection_SC[i, j, k] - selfcure_SC[i, j, k] - progress[i, j, k],
  p_TBdeaths_SC[i, j, k]
)



## g) IN MIGRATION
pSC_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pSC_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pSC_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
SC_inmigr[, , ] <- if(pSC_inmigr[i, j, k] > 0 && Nevents_SC[i, j, k] - age_out_SC[i, j, k] - HIV_out_SC[i, j, k] -detection_SC[i, j, k] - selfcure_SC[i, j, k] - progress[i, j, k] - TBdeaths_SC[i, j, k] > 0) rbinom(
  Nevents_SC[i, j, k] - age_out_SC[i, j, k] - HIV_out_SC[i, j, k] -
  detection_SC[i, j, k] - selfcure_SC[i, j, k] - progress[i, j, k] - TBdeaths_SC[i, j, k],
  pSC_inmigr[i, j, k]
) else 0


## g) BACKGROUND DEATH
SCdeaths[, , ] <- Nevents_SC[i, j, k] - age_out_SC[i, j, k] -
  HIV_out_SC[i, j, k] - detection_SC[i, j, k] -
  selfcure_SC[i, j, k] - progress[i, j, k] - TBdeaths_SC[i, j, k] - SC_inmigr[i, j, k]



#################################################
##### 4: prevalent infectious/active disease ####
#################################################
# print("D: {min(D[,,])}",when=min(D[,,])<0)
update(D[, , ]) <- D[i, j, k] + age_in_D[i, j, k] + HIV_in_D[i, j, k] + progress[i, j, k] +
  D_inmigr[i, j, k] -
  age_out_D[i, j, k] - HIV_out_D[i, j, k] - Ddeaths[i, j, k] - TBdeaths[i, j, k] -
  detection[i, j, k] - selfcure[i, j, k] - regress[i, j, k]

# ###### ACTIVE DISEASE: ENTRY EVENTS #####

## AGEING
age_in_D[,2:age_dims,] <- age_out_D[i,j-1,k]
age_in_D[,1,]<-0

## HIV
HIV_in_D[,,2:HIV_dims] <- HIV_out_D[i,j,k-1]
HIV_in_D[,,1]<-0

# ###### ACTIVE DISEASE: EXIT EVENTS #####
rate_D[,,1] <- m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate +
  age_rate[j] + HIV_rate_yr[j]
rate_D[,,2] <- m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate +
  age_rate[j] + ART_rate_yr[j]
rate_D[, , 3] <- m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate + age_rate[j]
Nevents_D[, , ] <- rbinom(D[i, j, k], 1 - exp(-rate_D[i, j, k] * dt))


## a) AGEING
p_Dage[, , ] <- if (rate_D[i, j, k] > tol) age_rate[j] / rate_D[i, j, k] else 0
age_out_D[, , ] <- rbinom(Nevents_D[i, j, k], p_Dage[i, j, k])


## b) HIV
p_DHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate + HIV_rate_yr[j] > tol) HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate + HIV_rate_yr[j]) else 0
p_DHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate + ART_rate_yr[j] > tol) ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate + ART_rate_yr[j]) else 0
p_DHIV[, , 3] <- 0
HIV_out_D[, , ] <- rbinom(Nevents_D[i, j, k] - age_out_D[i, j, k], p_DHIV[i, j, k])



# ## c) TB DETECTION
p_detect[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate > tol) dtct_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate) else 0
p_detect[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate > tol) dtct_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate) else 0
p_detect[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate > tol) dtct_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate + dtct_rate + slfcr_rate + regress_rate) else 0
detection[, , ] <- rbinom(Nevents_D[i, j, k] - age_out_D[i, j, k] - HIV_out_D[i, j, k], p_detect[i, j, k])


# ## d) SELF-CURING
p_selfcr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + slfcr_rate + regress_rate > tol) slfcr_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + slfcr_rate + regress_rate) else 0
p_selfcr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + slfcr_rate + regress_rate > tol) slfcr_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + slfcr_rate + regress_rate) else 0
p_selfcr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + slfcr_rate + regress_rate > tol) slfcr_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate + slfcr_rate + regress_rate) else 0
selfcure[, , ] <- rbinom(
  Nevents_D[i, j, k] - age_out_D[i, j, k] - HIV_out_D[i, j, k] - detection[i, j, k],
  p_selfcr[i, j, k]
)


# ## e) REGRESSION to Subclinical
p_regress[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + regress_rate > tol) regress_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate + regress_rate) else 0
p_regress[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate + regress_rate > tol) regress_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate + regress_rate) else 0
p_regress[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate + regress_rate > tol) regress_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate + regress_rate) else 0
regress[, , ] <- rbinom(
  Nevents_D[i, j, k] - age_out_D[i, j, k] - HIV_out_D[i, j, k] -
    detection[i, j, k] - selfcure[i, j, k],
  p_regress[i, j, k]
)



# ## f) TB DEATH
p_TBdeaths[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + TBd_rate > tol) TBd_rate / (m_in_t[j] + mu_noHIV_t[j] + TBd_rate) else 0
p_TBdeaths[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + TBd_rate > tol) TBd_rate / (m_in_t[j] + mu_HIV_t[j] + TBd_rate) else 0
p_TBdeaths[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + TBd_rate > tol) TBd_rate / (m_in_t[j] + mu_ART_t[j] + TBd_rate) else 0
TBdeaths[, , ] <- rbinom(
  Nevents_D[i, j, k] - age_out_D[i, j, k] - HIV_out_D[i, j, k] - detection[i, j, k] -
    selfcure[i, j, k] - regress[i, j, k],
  p_TBdeaths[i, j, k]
)


## g) IN MIGRATION
pD_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pD_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pD_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
D_inmigr[, , ] <- if(pD_inmigr[i, j, k]>0 && Nevents_D[i, j, k] - age_out_D[i, j, k] - HIV_out_D[i, j, k] - detection[i, j, k] -selfcure[i, j, k] - regress[i, j, k] - TBdeaths[i, j, k]>0) rbinom(
  Nevents_D[i, j, k] - age_out_D[i, j, k] - HIV_out_D[i, j, k] - detection[i, j, k] -
  selfcure[i, j, k] - regress[i, j, k] - TBdeaths[i, j, k],
  pD_inmigr[i, j, k]
) else 0


# ## h) BACKGROUND DEATH
Ddeaths[, , ] <- Nevents_D[i, j, k] - age_out_D[i, j, k] - HIV_out_D[i, j, k] - detection[i, j, k] -
  selfcure[i, j, k] - regress[i, j, k] - TBdeaths[i, j, k] - D_inmigr[i, j, k]


###########################
##### 5: On treatment #####
###########################
# print("Tr: {min(Tr[,,])}",when=min(Tr[,,])<0)
update(Tr[, , ]) <- Tr[i, j, k] + age_in_Tr[i, j, k] + HIV_in_Tr[i, j, k] +
  detection[i, j, k] + detection_SC[i, j, k] + Tr_inmigr[i, j, k] -
  age_out_Tr[i, j, k] - HIV_out_Tr[i, j, k] - Trdeaths[i, j, k] -
  Trsuccess[i, j, k] - Trxdeaths[i, j, k]

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
Nevents_Tr[, , ] <- rbinom(Tr[i, j, k], 1 - exp(-rate_Tr[i, j, k] * dt))


## a) AGEING
p_Trage[, , ] <- if (rate_Tr[i, j, k] > tol) age_rate[j] / rate_Tr[i, j, k] else 0
age_out_Tr[, , ] <- rbinom(Nevents_Tr[i, j, k], p_Trage[i, j, k])


## b) HIV
p_TrHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + Trsucc_rate + Trxd_rate + HIV_rate_yr[j] > tol) HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + Trsucc_rate + Trxd_rate + HIV_rate_yr[j]) else 0
p_TrHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + Trsucc_rate + Trxd_rate + ART_rate_yr[j] > tol) ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + Trsucc_rate + Trxd_rate + ART_rate_yr[j]) else 0
p_TrHIV[, , 3] <- 0
HIV_out_Tr[, , ] <- rbinom(Nevents_Tr[i, j, k] - age_out_Tr[i, j, k], p_TrHIV[i, j, k])



## c) TREATMENT SUCCESS
p_Trsuccess[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + Trsucc_rate + Trxd_rate > tol) Trsucc_rate / (m_in_t[j] + mu_noHIV_t[j] + Trsucc_rate + Trxd_rate) else 0
p_Trsuccess[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + Trsucc_rate + Trxd_rate > tol) Trsucc_rate / (m_in_t[j] + mu_HIV_t[j] + Trsucc_rate + Trxd_rate) else 0
p_Trsuccess[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + Trsucc_rate + Trxd_rate > tol) Trsucc_rate / (m_in_t[j] + mu_ART_t[j] + Trsucc_rate + Trxd_rate) else 0
Trsuccess[, , ] <- rbinom(
  Nevents_Tr[i, j, k] - age_out_Tr[i, j, k] - HIV_out_Tr[i, j, k],
  p_Trsuccess[i, j, k]
)



## d) TB DEATH
p_Trxdeaths[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + Trxd_rate > tol) Trxd_rate / (m_in_t[j] + mu_noHIV_t[j] + Trxd_rate) else 0
p_Trxdeaths[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + Trxd_rate > tol) Trxd_rate / (m_in_t[j] + mu_HIV_t[j] + Trxd_rate) else 0
p_Trxdeaths[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + Trxd_rate > tol) Trxd_rate / (m_in_t[j] + mu_ART_t[j] + Trxd_rate) else 0
Trxdeaths[, , ] <- rbinom(
  Nevents_Tr[i, j, k] - age_out_Tr[i, j, k] - HIV_out_Tr[i, j, k] - Trsuccess[i, j, k],
  p_Trxdeaths[i, j, k]
)

## e) IN MIGRATION
pTr_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pTr_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pTr_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
Tr_inmigr[, , ] <- if(pTr_inmigr[i, j, k]>0 && Nevents_Tr[i, j, k] - age_out_Tr[i, j, k] - HIV_out_Tr[i, j, k] - Trsuccess[i, j, k] - Trxdeaths[i, j, k]>0)rbinom(
  Nevents_Tr[i, j, k] - age_out_Tr[i, j, k] - HIV_out_Tr[i, j, k] -
  Trsuccess[i, j, k] - Trxdeaths[i, j, k],
  pTr_inmigr[i, j, k]
) else 0


## e) BACKGROUND DEATH
Trdeaths[, , ] <- Nevents_Tr[i, j, k] - age_out_Tr[i, j, k] - HIV_out_Tr[i, j, k] -
  Trsuccess[i, j, k] - Trxdeaths[i, j, k] - Tr_inmigr[i, j, k]

########################
##### 6: Recovered #####
########################
# print("R: {min(R[,,])}",when=min(R[,,])<0)

update(R[, , ]) <- R[i, j, k] + age_in_R[i, j, k] + HIV_in_R[i, j, k] +
  Trsuccess[i, j, k] + R_inmigr[i, j, k] -
  age_out_R[i, j, k] - HIV_out_R[i, j, k] - Rdeaths[i, j, k] - Rinfs[i, j, k] - relapse[i, j, k]

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
Nevents_R[, , ] <- rbinom(R[i, j, k], 1 - exp(-rate_R[i, j, k] * dt))


## a) AGEING
p_Rage[, , ] <- if (rate_R[i, j, k] > tol) age_rate[j] / rate_R[i, j, k] else 0
age_out_R[, , ] <- rbinom(Nevents_R[i, j, k], p_Rage[i, j, k])


## b) HIV rate
p_RHIV[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] + pDr * IRR[i] * Hirr[1] + HIV_rate_yr[j] > tol) HIV_rate_yr[j] / (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] + pDr * IRR[i] * Hirr[1] + HIV_rate_yr[j]) else 0
p_RHIV[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] + pDr * IRR[i] * Hirr[2] + ART_rate_yr[j] > tol) ART_rate_yr[j] / (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] + pDr * IRR[i] * Hirr[2] + ART_rate_yr[j]) else 0
p_RHIV[, , 3] <- 0
HIV_out_R[, , ] <- rbinom(
  Nevents_R[i, j, k] - age_out_R[i, j, k],
  p_RHIV[i, j, k]
)


## c) REINFECTION
p_Rinfs[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] + pDr * IRR[i] * Hirr[1] > tol) v * foi[i] * TB_HIV_mod[1] / (m_in_t[j] + mu_noHIV_t[j] + v * foi[i] * TB_HIV_mod[1] + pDr * IRR[i]* Hirr[1]) else 0
p_Rinfs[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] + pDr * IRR[i] * Hirr[2] > tol) v * foi[i] * TB_HIV_mod[1] / (m_in_t[j] + mu_HIV_t[j] + v * foi[i] * TB_HIV_mod[2] + pDr * IRR[i] * Hirr[2]) else 0
p_Rinfs[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + v * foi[i] * TB_HIV_mod[3] + pDr * IRR[i] * Hirr[3] > tol) v * foi[i] * TB_HIV_mod[1] / (m_in_t[j] + mu_ART_t[j] + v * foi[i] * TB_HIV_mod[3] + pDr * IRR[i] * Hirr[3]) else 0
Rinfs[, , ] <- rbinom(
  Nevents_R[i, j, k] - age_out_R[i, j, k] - HIV_out_R[i, j, k],
  p_Rinfs[i, j, k]
)


## d) RELAPSE
p_relapse[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] + pDr * IRR[i] * Hirr[1] > tol) IRR[i] * pDr * Hirr[1] / (m_in_t[j] + mu_noHIV_t[j] + pDr * IRR[i] * Hirr[1]) else 0
p_relapse[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] + pDr * IRR[i] * Hirr[2] > tol) IRR[i] * pDr * Hirr[2] / (m_in_t[j] + mu_HIV_t[j] + pDr * IRR[i] * Hirr[2]) else 0
p_relapse[, , 3] <- if (m_in_t[j] + mu_ART_t[j] + pDr * IRR[i] * Hirr[3] > tol) IRR[i] * pDr * Hirr[3] / (m_in_t[j] + mu_ART_t[j] + pDr * IRR[i] * Hirr[3]) else 0
relapse[, , ] <- rbinom(
  Nevents_R[i, j, k] - age_out_R[i, j, k] - HIV_out_R[i, j, k] - Rinfs[i, j, k],
  p_relapse[i, j, k]
)

## e) IN MIGRATION
pR_inmigr[, , 1] <- if (m_in_t[j] + mu_noHIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_noHIV_t[j]) else 0
pR_inmigr[, , 2] <- if (m_in_t[j] + mu_HIV_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_HIV_t[j]) else 0
pR_inmigr[, , 3] <- if (m_in_t[j] + mu_ART_t[j] > tol) m_in_t[j] / (m_in_t[j] + mu_ART_t[j]) else 0
R_inmigr[, , ] <- if(pR_inmigr[i, j, k]>0 && Nevents_R[i, j, k] - age_out_R[i, j, k] - HIV_out_R[i, j, k] - Rinfs[i, j, k] - relapse[i, j, k]>0) rbinom(
  Nevents_R[i, j, k] - age_out_R[i, j, k] - HIV_out_R[i, j, k] - Rinfs[i, j, k] - relapse[i, j, k],
  pR_inmigr[i, j, k]
) else 0

## e) BACKGROUND DEATH
Rdeaths[, , ] <- Nevents_R[i, j, k] - age_out_R[i, j, k] - HIV_out_R[i, j, k] -
  Rinfs[i, j, k] - relapse[i, j, k] - R_inmigr[i, j, k]


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

## ## patch-level total notifications,popn, and rates
## pnotes[] <- sum(detection[i, , ] + detection_SC[i, , ])
## ppop[] <- sum(N[i, , ])
## pnotifrate[] <- 1e5 * pnotes[i] / (ppop[i] + 1e-10)

## #########################
## test variables
update(beta_test) <- beta_test #constant at initial value passed in
