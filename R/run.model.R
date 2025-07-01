##' This will run the stochastic version of the model without reference to data
##'
##' TODO
##' @title Run dust version of model
##' @param parm 
##' @param times 
##' @param n.particles 
##' @param convert 
##' @return array
##' @author Pete Dodd
##' @import dust
##' @export
run.model <- function(parm,
                      times,
                      n.particles = 1,
                      convert = FALSE) {
  ## create dust model
  dmod <- BLASTtbmod:::stocm$new(
    pars = parm,
    time = min(times) + 1, # TODO
    n_particles = n.particles
  )
  ## model loop over time
  y_looped <- array(NA, dim = c(dmod$info()$len, n.particles, max(times)))
  for (t in 1:max(times)) {
    y_looped[, , t] <- dmod$run(t)
  }
  if (convert == TRUE) { # convert to df
    y_looped <- as.data.frame(t(y_looped[, 1, ])) # TODO cope with multiple parms?
    ## Make it readable
    if (ncol(y_looped) == length(BLASTtbmod::get_cols)) {
      colnames(y_looped) <- BLASTtbmod::get_cols
    }
  }
  y_looped
}

##' For simulating counterfactuals
##'
##' TODO
##' @title run counterfactuals
##' @param parm base case parameters
##' @param cfparm counterfactual parameters
##' @param times times to run over
##' @param n.particles number of particles
##' @param pars_multi TODO
##' @param convert to dataframe?
##' @return array or data.frame
##' @import dust
##' @export
##' @author Pete Dodd
run.CF <- function(parm,
                   cfparm,
                   times,
                   n.particles = 1,
                   pars_multi=FALSE,
                   convert = FALSE) {
  ## create basecase dust model
  bcmod <- BLASTtbmod:::stocm$new(
                                pars = parm,
                                time = min(times) + 1,
                                n_particles = 1, seed = 1,
                                deterministic = TRUE,
                                pars_multi = pars_multi
                              )
  ## create counterfactual dust model
  cfmod <- BLASTtbmod:::stocm$new(
                                pars = cfparm,
                                time = min(times) + 1,
                                n_particles = 1, seed = 1,
                                deterministic = TRUE,
                                pars_multi = pars_multi
                              )
  ## simulate
  ## bcres <- bcmod$simulate(times)
  ## cfres <- cfmod$simulate(times)
  ## model loop over time
  bcres <- cfres <- array(NA, dim = c(bcmod$info()$len, n.particles, max(times)))
  for (t in 1:max(times)) {
    bcres[, , t] <- bcmod$run(t)
    cfres[, , t] <- cfmod$run(t)
  }
  if (convert == TRUE) { # convert to df
    bcres <- as.data.frame(t(bcres[, 1, ])) # TODO cope with multiple parms?
    cfres <- as.data.frame(t(cfres[, 1, ])) # TODO cope with multiple parms?
    ## Make it readable
    if (ncol(bcres) == length(BLASTtbmod::get_cols)) {
      colnames(bcres) <- colnames(cfres) <- BLASTtbmod::get_cols
    }
  }
  list(bcres,cfres)
}
