##' This is the main function for performing particle filter inference
##'
##' TODO
##' @title Run PMCMC inference on particle filter
##' @param particle.filter an mcstate particle filter object
##' @param parms base parameters list
##' @param prior.list priors for inference parameters
##' @param initial.proposal.matrix starting point for PMCMC proposal matric
##' @param n.steps number of PMCMC iterations
##' @param n.burnin number of interations to discard
##' @param n.chains number of PMCMC chains to run
##' @param n.threads number of threads to use (not tested)
##' @param n.epochs number of times to re-run PMCMC updating the proposal matrix based on previous run
##' @param save_restart where to save the restart
##' @param transform if not supplied will default to transform = function(theta) c(parms, as.list(theta))
##' @param mcmc_pars if not supplied, will create from prior.list & initial.proposal.matrix (& optionally transform)
##' @param returnall whether to return a list of things or just processed chains
##' @return mcstate processed chains object (unless returnall=TRUE, in which case a list)
##' @author Pete Dodd
##' @import mcstate
##' @export
run.pmcmc <- function(particle.filter,
                      parms,
                      prior.list=NULL,
                      initial.proposal.matrix=NULL,
                      n.steps, n.burnin, n.chains,
                      n.threads = 4,
                      n.epochs = 1,
                      save_restart=NULL,
                      transform=NULL,
                      mcmc_pars=NULL,
                      returnall=FALSE) {
  ## if not using mcmc_pars directly, and missing: create a default transform
  if(is.null(transform) & is.null(mcmc_pars)) transform <- function(theta) c(parms, as.list(theta))
  if(is.null(mcmc_pars) & !(is.null(prior.list) & is.null(initial.proposal.matrix) & is.null(transform)) )
    warning("Since mcmc_pars is supplied, will be ignoring prior.list, initial.proposal.matrix, and transform")
  if(is.null(mcmc_pars) & is.null(prior.list) & is.null(initial.proposal.matrix) & is.null(transform))
    stop("Need either mcmc_pars, or both prior.list & initial.proposal.matrix (& optionally transform)!")
  if(!is.null(initial.proposal.matrix)){
    proposal.matrix <- initial.proposal.matrix
  }
  ## start epoch loop
  for(epoch in 1:n.epochs){
    cat("------ starting epoch ", epoch, " / ", n.epochs, " ------\n")
    ## create PF parameters
    if( !is.null(mcmc_pars) & epoch==1 ){ #extract this triplet from provided mcmc_pars for potential future reference
      prior.list <- mcmc_pars$.__enclos_env__$private$parameters
      proposal.matrix <- mcmc_pars$.__enclos_env__$private$proposal_kernel
      transform <- mcmc_pars$.__enclos_env__$private$transform
    } else { #create mcmc_pars (or overwrite)
      mcmc_pars <- mcstate::pmcmc_parameters$new(
                                               prior.list,
                                               proposal.matrix,
                                               transform = transform
                                             )
    }
    ## PF control object
    control <- mcstate::pmcmc_control(
                          n_steps = n.steps,
                          save_state = TRUE,
                          save_trajectories = TRUE,
                          save_restart = save_restart,
                          restart_match = TRUE,
                          n_chains = n.chains,
                          n_threads_total = n.threads,
                          progress = TRUE
                        )
    ## run filter
    pmcmc_run <- mcstate::pmcmc(mcmc_pars,
                                particle.filter,
                                control = control
                                )
    ## Remove burnin and thin if desired:
    processed_chains <- mcstate::pmcmc_thin(pmcmc_run, burnin = n.burnin)
    ## make new proposal matrix
    proposal.matrix <- cov(pmcmc_run$pars)
  }
  if(returnall){
    return(list(pmcmc_run=pmcmc_run,
                processed_chains=processed_chains,
                proposal_matrix=proposal.matrix))
  }
  processed_chains
}
