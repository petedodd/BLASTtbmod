##' This is the main function for performing particle filter inference
##'
##' TODO
##' @title Run PMCMC inference on particle filter
##' @param particle.filter TODO
##' @param parms TODO
##' @param prior.list TODO
##' @param initial.proposal.matrix TODO
##' @param n.steps TODO
##' @param n.burnin TODO
##' @param n.chains TODO
##' @param n.threads TODO
##' @param n.epochs TODO
##' @param save_restart where to save the restart
##' @param returnall whether to return a list of things or just processed chains
##' @return mcstate processed chains object (unless returnall=TRUE, in which case a list)
##' @author Pete Dodd
##' @import mcstate
##' @export
run.pmcmc <- function(particle.filter,
                      parms,
                      prior.list,
                      initial.proposal.matrix,
                      n.steps, n.burnin, n.chains,
                      n.threads = 4,
                      n.epochs = 1,
                      save_restart=NULL,
                      returnall=FALSE) {
  proposal.matrix <- initial.proposal.matrix
  for(epoch in 1:n.epochs){
    cat("------ starting epoch ", epoch, " / ", n.epochs, " ------\n")
    ## create PF parameters
    mcmc_pars <- mcstate::pmcmc_parameters$new(
                                             prior.list,
                                             proposal.matrix,
                                             transform = function(theta) c(parms, as.list(theta))
                                           )
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
