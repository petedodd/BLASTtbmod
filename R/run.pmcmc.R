##' This is the main function for performing particle filter inference
##'
##' TODO
##' @title Run PMCMC inference on particle filter
##' @param particle.filter 
##' @param parms 
##' @param prior.list 
##' @param initial.proposal.matrix 
##' @param n.steps 
##' @param n.burnin 
##' @param n.chains 
##' @param n.threads 
##' @param n.epochs 
##' @return mcstate processed chains object
##' @author Pete Dodd
##' @import mcstate
##' @export
run.pmcmc <- function(particle.filter,
                      parms,
                      prior.list,
                      initial.proposal.matrix,
                      n.steps, n.burnin, n.chains,
                      n.threads = 4, n.epochs = 1) {
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
  processed_chains
}
